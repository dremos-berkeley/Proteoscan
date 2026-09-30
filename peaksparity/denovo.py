"""De novo sequencing with Casanovo, for ion-trap (low-resolution MS2) HCD data.

PEAKS's de novo result feeds its ``denovoOnly`` export and the Tag Length / CAA columns of the
DB tables. Our open substitute is Casanovo (MIT/Apache) fine-tuned on ion-trap HCD PSMs; see
docs/DENOVO.md for the training recipe and the measured accuracy.

Main entry points

* :func:`to_casanovo` / :func:`from_casanovo`: PEAKS <-> Casanovo (ProForma) peptide notation.
* :func:`training_psms` + :func:`write_annotated_mgf`: build Casanovo training MGFs
  (``SEQ=`` lines) from a PEAKS-format ``db.psms.csv`` and the run's mzML.
* :func:`sequence_mzml`: run ``casanovo sequence`` and return one row per spectrum
  (scan, ProForma sequence, score, per-residue scores).
* :func:`evaluate`: full-sequence accuracy, Casanovo AA precision/recall, accuracy vs score.
* :func:`denovo_only_table`: rows in the PEAKS ``<Sample>.db.denovoOnly.csv`` column layout.

Scores: Casanovo's peptide score is the product of the per-residue softmax probabilities.
PEAKS-style ALC is computed here as mean(per-residue score) * 100 and local confidence as
per-residue score * 100 (rounded to integers).
"""
from __future__ import annotations

import logging
import os
import re
import subprocess
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .seqnorm import canonical, format_delta, parse, strip_mods

log = logging.getLogger(__name__)

PROTON = 1.007276466
H2O = 18.010565
C13 = 1.003355

#: Casanovo 5.2 default residue vocabulary (config.yaml ``residues``); identical to the
#: tokenizer stored in casanovo_orbitrap_v5-2-0.ckpt (which also has plain ``C``).
CASANOVO_RESIDUES = {
    "G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414,
    "T": 101.047670, "C": 103.009185, "C[Carbamidomethyl]": 160.030649, "L": 113.084064,
    "I": 113.084064, "N": 114.042927, "D": 115.026943, "Q": 128.058578, "K": 128.094963,
    "E": 129.042593, "M": 131.040485, "H": 137.058912, "F": 147.068414, "R": 156.101111,
    "Y": 163.063329, "W": 186.079313,
    "M[Oxidation]": 147.035400, "N[Deamidated]": 115.026943, "Q[Deamidated]": 129.042594,
    "[Acetyl]-": 42.010565, "[Carbamyl]-": 43.005814, "[Ammonia-loss]-": -17.026549,
    "[+25.980265]-": 25.980265,
}

#: (residue, PEAKS 2-decimal delta) -> Casanovo token
_RESIDUE_MOD_TOKENS = {
    ("C", "+57.02"): "C[Carbamidomethyl]",
    ("M", "+15.99"): "M[Oxidation]",
    ("N", "+0.98"): "N[Deamidated]",
    ("Q", "+0.98"): "Q[Deamidated]",
}
_NTERM_TOKENS = {
    "+42.01": "[Acetyl]-",
    "+43.01": "[Carbamyl]-",
    "-17.03": "[Ammonia-loss]-",
    "+25.98": "[+25.980265]-",
}
_TOKEN_RE = re.compile(r"\[[^\]]+\]-|[A-Z](?:\[[^\]]+\])?")


# --------------------------------------------------------------------------- notation
def to_casanovo(peptide: str, fixed_cam: bool = False) -> str | None:
    """PEAKS/ProForma/any seqnorm-parsable peptide -> Casanovo ProForma tokens.

    ``PEPC(+57.02)M(+15.99)N(+0.98)K`` -> ``PEPC[Carbamidomethyl]M[Oxidation]N[Deamidated]K``.
    ``fixed_cam=True`` adds carbamidomethyl to every bare C (for tools that omit fixed mods;
    PEAKS writes it explicitly). Returns None if a modification is outside the vocabulary.
    """
    stripped, mods = parse(peptide)
    md = {p: format_delta(d) for p, d in mods}
    n = len(stripped)
    out = []
    if 0 in md:
        tok = _NTERM_TOKENS.get(md.pop(0))
        if tok is None:
            return None
        out.append(tok)
    if n + 1 in md:
        return None
    for k, aa in enumerate(stripped, start=1):
        if k in md:
            tok = _RESIDUE_MOD_TOKENS.get((aa, md[k]))
            if tok is None:
                return None
            out.append(tok)
        elif aa == "C" and fixed_cam:
            out.append("C[Carbamidomethyl]")
        else:
            if aa not in CASANOVO_RESIDUES:
                return None
            out.append(aa)
    return "".join(out)


def from_casanovo(seq: str) -> str:
    """Casanovo ProForma -> PEAKS-style canonical string (``M[Oxidation]`` -> ``M(+15.99)``)."""
    return canonical(seq)


def tokens(seq: str) -> list[str]:
    """Split a Casanovo ProForma sequence into vocabulary tokens (N-terminal mod separate)."""
    return _TOKEN_RE.findall(seq)


def peptide_mass(seq: str) -> float:
    """Neutral monoisotopic mass of a Casanovo ProForma peptide."""
    return sum(CASANOVO_RESIDUES[t] for t in tokens(seq)) + H2O


def precursor_ppm(seq: str, mz: float, z: int, isotopes: Iterable[int] = (0, 1)) -> float:
    """Smallest |ppm| between the peptide and the observed precursor over isotope offsets."""
    theo = peptide_mass(seq)
    obs = (mz - PROTON) * z
    return min((abs(obs - i * C13 - theo) / theo * 1e6 for i in isotopes))


# --------------------------------------------------------------------------- training data
def training_psms(psms_csv: str | Path, min_score: float = 30.0,
                  unique_peptide_charge: bool = True) -> pd.DataFrame:
    """Training labels from a PEAKS-format ``db.psms.csv``.

    Keeps PSMs with -10LgP >= ``min_score`` whose modifications are in Casanovo's vocabulary,
    collapses repeated (scan, peptide, z) rows (chimeric feature duplicates), and, if
    ``unique_peptide_charge``, keeps only the best-scoring PSM per (peptide, z).
    Returns columns scan, seq (Casanovo), peptide (PEAKS), charge, mz, feature_rt (min), score.
    """
    d = pd.read_csv(psms_csv, usecols=["Peptide", "-10LgP", "m/z", "z", "RT", "Scan"])
    d = d[d["-10LgP"] >= min_score].copy()
    d["seq"] = d["Peptide"].map(to_casanovo)
    n_bad = int(d["seq"].isna().sum())
    if n_bad:
        log.info("%s: %d PSMs with out-of-vocabulary modifications dropped", psms_csv, n_bad)
    d = d.dropna(subset=["seq"])
    d = d.sort_values(["-10LgP", "Scan"], ascending=[False, True], kind="stable")
    # I/L variant rows of one PSM (PEAKS reports both) are one training example
    d["_il"] = d["Peptide"].str.replace("I", "L")
    d = d.drop_duplicates(["Scan", "_il", "z"])
    if unique_peptide_charge:
        d = d.drop_duplicates(["_il", "z"])
    out = pd.DataFrame({
        "scan": d["Scan"].astype(int).to_numpy(), "seq": d["seq"].to_numpy(),
        "peptide": d["Peptide"].to_numpy(), "charge": d["z"].astype(int).to_numpy(),
        "mz": d["m/z"].astype(float).to_numpy(), "feature_rt": d["RT"].astype(float).to_numpy(),
        "score": d["-10LgP"].astype(float).to_numpy(),
    })
    return out.sort_values("scan", kind="stable").reset_index(drop=True)


def _ms2_index(mzml: str | Path):
    import pyopenms as oms
    od = oms.OnDiscMSExperiment()
    if not od.openFile(str(mzml)):
        raise FileNotFoundError(mzml)
    meta = od.getMetaData()
    index = {}
    for i in range(meta.getNrSpectra()):
        s = meta[i]
        if s.getMSLevel() != 2:
            continue
        m = re.search(r"scan=(\d+)", s.getNativeID())
        if m:
            index.setdefault(int(m.group(1)), i)
    return od, meta, index


def mzml_precursors(mzml: str | Path) -> pd.DataFrame:
    """scan, mz, charge, rt (s) of every MS2 spectrum, from the mzML precursor records."""
    od, meta, index = _ms2_index(mzml)
    rows = []
    for scan, i in index.items():
        s = meta[i]
        pcs = s.getPrecursors()
        mz = pcs[0].getMZ() if pcs else np.nan
        z = pcs[0].getCharge() if pcs else 0
        rows.append((scan, mz, z, s.getRT()))
    return pd.DataFrame(rows, columns=["scan", "mz", "charge", "rt"]).sort_values("scan").reset_index(drop=True)


def write_annotated_mgf(entries: pd.DataFrame, mzml: str | Path, out_mgf: str | Path,
                        title_prefix: str | None = None, return_rows: bool = False):
    """Write an MGF (optionally annotated) of the scans in ``entries``.

    ``entries`` columns: scan, and optionally charge, mz, rt, seq. Missing precursor values are
    taken from the mzML. A ``seq`` column produces ``SEQ=`` lines (Casanovo training format).
    Entries sharing a scan (chimeric precursors) are written as separate spectra.
    Returns the number of spectra written, or with ``return_rows`` the positions (0-based,
    in ``entries``) of the written rows, in MGF order (MGF ``index=k`` <-> ``rows[k]``).
    """
    od, meta, index = _ms2_index(mzml)
    prefix = title_prefix or Path(mzml).stem
    n = 0
    written = []
    with open(out_mgf, "w") as fh:
        for pos, row in enumerate(entries.itertuples(index=False)):
            i = index.get(int(row.scan))
            if i is None:
                continue
            spec = od.getSpectrum(i)
            mz_arr, it_arr = spec.get_peaks()
            if len(mz_arr) == 0:
                continue
            pcs = meta[i].getPrecursors()
            pmz = getattr(row, "mz", None)
            if pmz is None or not np.isfinite(pmz):
                pmz = pcs[0].getMZ()
            z = getattr(row, "charge", None)
            if z is None or not np.isfinite(z) or int(z) <= 0:
                z = pcs[0].getCharge()
            rt = getattr(row, "rt", None)
            if rt is None or not np.isfinite(rt):
                rt = meta[i].getRT()
            lines = [
                "BEGIN IONS",
                f"TITLE={prefix}.{int(row.scan)}.{int(row.scan)}.{int(z)}",
                f"PEPMASS={pmz:.6f}",
                f"CHARGE={int(z)}+",
                f"SCANS={int(row.scan)}",
                f"RTINSECONDS={rt:.3f}",
            ]
            seq = getattr(row, "seq", None)
            if isinstance(seq, str) and seq:
                lines.append(f"SEQ={seq}")
            lines += [f"{m:.5f} {v:.2f}" for m, v in zip(mz_arr, it_arr)]
            lines.append("END IONS")
            fh.write("\n".join(lines) + "\n")
            n += 1
            written.append(pos)
    return written if return_rows else n


# --------------------------------------------------------------------------- inference
def read_mztab(path: str | Path, mgf: str | Path | None = None) -> pd.DataFrame:
    """Casanovo mzTab -> DataFrame(scan, sequence, peptide, score, per_residue_scores, charge,
    exp_mz, calc_mz, rank). ``scan`` is the raw scan: parsed from the mzML native id, or from the
    MGF ``SCANS=`` column Casanovo writes (``mgf`` may be given to map ``index=N`` otherwise)."""
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("PSH"):
                header = line.rstrip("\n").split("\t")
            elif line.startswith("PSM"):
                rows.append(line.rstrip("\n").split("\t"))
    if header is None:
        return pd.DataFrame(columns=["scan", "sequence", "peptide", "score", "per_residue_scores",
                                     "charge", "exp_mz", "calc_mz", "rank"])
    d = pd.DataFrame(rows, columns=header)
    ref = d["spectra_ref"].str.split(":", n=1).str[1]
    scan = ref.str.extract(r"scan=(\d+)")[0]
    if "opt_global_cv_MS:1003057_scan_number" in d.columns:
        scan = scan.fillna(d["opt_global_cv_MS:1003057_scan_number"].str.extract(r"(\d+)$")[0])
    if scan.isna().any() and mgf is not None:
        idx2scan = {}
        k = 0
        with open(mgf) as fh:
            cur = None
            for line in fh:
                if line.startswith("SCANS="):
                    cur = line.strip().split("=", 1)[1]
                elif line.startswith("END IONS"):
                    idx2scan[f"index={k}"] = cur
                    k += 1
        scan = scan.fillna(ref.map(idx2scan))
    seq = d["opt_global_cv_MS:1003169_proforma_peptidoform_sequence"]
    out = pd.DataFrame({
        "scan": pd.to_numeric(scan, errors="coerce").astype("Int64"),
        "spectra_ref": ref,
        "sequence": seq,
        "score": pd.to_numeric(d["search_engine_score[1]"], errors="coerce"),
        "per_residue_scores": d["opt_global_aa_scores"].map(
            lambda s: np.array([float(x) for x in s.split(",")]) if s and s != "null" else np.array([])),
        "charge": pd.to_numeric(d["charge"], errors="coerce").astype("Int64"),
        "exp_mz": pd.to_numeric(d["exp_mass_to_charge"], errors="coerce"),
        "calc_mz": pd.to_numeric(d["calc_mass_to_charge"], errors="coerce"),
    })
    out["peptide"] = out["sequence"].map(_safe_canonical)
    out["rank"] = out.groupby("spectra_ref", sort=False).cumcount() + 1
    return out


def _safe_canonical(s):
    try:
        return from_casanovo(s)
    except Exception:
        return None


def pick_precursor_consistent(df: pd.DataFrame, tol_ppm: float = 15.0,
                              isotopes: Iterable[int] = (0, 1)) -> pd.DataFrame:
    """From several candidates per spectrum (``top_match`` > 1), keep the highest-scoring one
    whose mass matches the precursor within ``tol_ppm`` (isotope offsets allowed); spectra with
    no consistent candidate keep their top candidate. Adds ``precursor_ok`` and ``ppm``."""
    d = df.copy()
    iso = tuple(isotopes)
    d["ppm"] = [precursor_ppm(s, m, int(z), iso) if isinstance(s, str) else np.inf
                for s, m, z in zip(d["sequence"], d["exp_mz"], d["charge"])]
    d["precursor_ok"] = d["ppm"] <= tol_ppm
    d = d.sort_values(["spectra_ref", "precursor_ok", "score"], ascending=[True, False, False],
                      kind="stable")
    return d.drop_duplicates("spectra_ref").reset_index(drop=True)


def write_config(path: str | Path, **overrides) -> Path:
    """Write a Casanovo YAML config: the package default with ``overrides`` applied."""
    import yaml
    import casanovo
    base = Path(casanovo.__file__).with_name("config.yaml")
    cfg = yaml.safe_load(base.read_text())
    cfg.update(overrides)
    Path(path).write_text(yaml.safe_dump(cfg, sort_keys=False))
    return Path(path)


def sequence_mzml(mzml: str | Path, model_ckpt: str | Path, out_mztab: str | Path,
                  scans: Iterable[int] | None = None, precursors: pd.DataFrame | None = None,
                  n_beams: int = 1, top_match: int = 1, precursor_tol_ppm: float = 15.0,
                  isotopes: Iterable[int] = (0, 1), config: str | Path | None = None,
                  batch_size: int = 1024, gpu: str | None = None,
                  casanovo_bin: str = "casanovo", reuse: bool = True) -> pd.DataFrame:
    """De novo sequence an mzML with Casanovo.

    Parameters
    ----------
    mzml : input mzML (or MGF) path.
    model_ckpt : Casanovo checkpoint (e.g. the fine-tuned ion-trap model).
    out_mztab : output mzTab path (Casanovo writes ``<stem>.mztab`` next to it).
    scans : optional subset of raw scan numbers; a temporary MGF with those scans is written.
    precursors : optional DataFrame(scan, charge, mz[, rt (s)], ...) overriding the mzML
        precursor (e.g. corrected monoisotopic m/z from feature detection; one spectrum per row,
        so a chimeric scan can be sequenced once per precursor). Implies an MGF. Its other
        columns (e.g. feature_id, area, feature_rt) are carried over to the returned rows.
    n_beams, top_match : beam search width and candidates returned per spectrum. With
        top_match > 1 the highest-scoring candidate within ``precursor_tol_ppm`` of the precursor
        is chosen (:func:`pick_precursor_consistent`). Casanovo 5.2 itself applies no
        precursor filter in de novo mode.
    gpu : value for CUDA_VISIBLE_DEVICES (None = inherit).

    Returns DataFrame(scan, sequence (ProForma), score, per_residue_scores, peptide (PEAKS
    notation), charge, exp_mz, calc_mz, precursor_ok, ppm).
    """
    out_mztab = Path(out_mztab)
    out_mztab.parent.mkdir(parents=True, exist_ok=True)
    stem = out_mztab.name[:-6] if out_mztab.name.endswith(".mztab") else out_mztab.name
    result = out_mztab.parent / f"{stem}.mztab"
    peak_path = Path(mzml)
    mgf = None
    if precursors is not None or scans is not None:
        if precursors is not None:
            ent = precursors.copy()
            if scans is not None:
                ent = ent[ent["scan"].isin(set(int(s) for s in scans))]
        else:
            ent = pd.DataFrame({"scan": sorted(set(int(s) for s in scans))})
        ent = ent.reset_index(drop=True)
        mgf = out_mztab.parent / f"{stem}.input.mgf"
        rows_tsv = out_mztab.parent / f"{stem}.input.rows.tsv"
        if not (reuse and mgf.exists() and rows_tsv.exists()):
            written = write_annotated_mgf(ent.drop(columns=["seq"], errors="ignore"), mzml, mgf,
                                          return_rows=True)
            ent.iloc[written].assign(spectra_ref=[f"index={k}" for k in range(len(written))]).to_csv(
                rows_tsv, sep="\t", index=False)
        peak_path = mgf
    if not (reuse and result.exists()):
        cfg_path = out_mztab.parent / f"{stem}.config.yaml"
        if config is None:
            write_config(cfg_path, n_beams=n_beams, top_match=top_match,
                         predict_batch_size=batch_size)
        else:
            import yaml
            cfg = yaml.safe_load(Path(config).read_text())
            cfg.update(n_beams=n_beams, top_match=top_match, predict_batch_size=batch_size)
            cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        env = dict(os.environ)
        if gpu is not None:
            env["CUDA_VISIBLE_DEVICES"] = str(gpu)
        cmd = [casanovo_bin, "sequence", "-m", str(model_ckpt), "-c", str(cfg_path),
               "-d", str(out_mztab.parent), "-o", stem, "-f", "-v", "warning", str(peak_path)]
        log.info("running %s", " ".join(cmd))
        subprocess.run(cmd, check=True, env=env)
    df = read_mztab(result, mgf=mgf)
    extra = []
    if mgf is not None and precursors is not None:
        # carry the precursor table's own columns (feature_id, rt, area, ...) to each prediction
        rows = pd.read_csv(out_mztab.parent / f"{stem}.input.rows.tsv", sep="\t")
        extra = [c for c in rows.columns if c not in ("scan", "charge", "mz", "seq", "spectra_ref")]
        df = df.merge(rows[["spectra_ref"] + extra], on="spectra_ref", how="left")
    if top_match > 1:
        df = pick_precursor_consistent(df, precursor_tol_ppm, isotopes)
    else:
        iso = tuple(isotopes)
        df["ppm"] = [precursor_ppm(s, m, int(z), iso) for s, m, z in
                     zip(df["sequence"], df["exp_mz"], df["charge"])]
        df["precursor_ok"] = df["ppm"] <= precursor_tol_ppm
    cols = ["scan", "sequence", "score", "per_residue_scores", "peptide", "charge", "exp_mz",
            "calc_mz", "precursor_ok", "ppm", "spectra_ref"] + extra
    return df[cols].reset_index(drop=True)


# --------------------------------------------------------------------------- evaluation
def _eval_tokens(seq: str) -> list[str]:
    return [("L" if t == "I" else t) for t in tokens(seq)]


def evaluate(pred: pd.DataFrame, truth: pd.DataFrame, key: str | list[str] = "scan") -> dict:
    """Compare de novo predictions with reference identifications.

    ``truth``: columns ``key`` and ``seq`` (Casanovo notation) or ``peptide`` (PEAKS notation).
    ``pred``: output of :func:`sequence_mzml` / :func:`read_mztab` (one row per ``key``).
    Every truth row counts; spectra without a prediction count as wrong (as in Casanovo).

    Returns dict with n, n_predicted, exact (exact modified sequence, I != L), exact_il
    (I = L), casanovo_pep (Casanovo's mass-based peptide match), aa_precision, aa_recall
    (Casanovo definitions), and ``curve``: DataFrame(coverage, score, precision_exact_il,
    precision_casanovo) sorted by descending score.
    """
    from casanovo.denovo.evaluate import aa_match_batch, aa_match_metrics
    t = truth.copy()
    if "seq" not in t.columns:
        t["seq"] = t["peptide"].map(to_casanovo)
    keys = [key] if isinstance(key, str) else list(key)
    m = t[keys + ["seq"]].merge(pred[keys + ["sequence", "score"]], on=keys, how="left")
    has = m["sequence"].notna()
    true_c = m["seq"].map(canonical)
    pred_c = m["sequence"].map(lambda s: _safe_canonical(s) if isinstance(s, str) else None)
    exact = (pred_c == true_c) & has
    exact_il = (pred_c.fillna("").str.replace("I", "L") == true_c.str.replace("I", "L")) & has
    t_tok = [_eval_tokens(s) for s in m["seq"]]
    p_tok = [(_eval_tokens(s) if isinstance(s, str) else None) for s in m["sequence"]]
    aa_dict = dict(CASANOVO_RESIDUES)
    matches, n_true, n_pred = aa_match_batch(t_tok, p_tok, aa_dict)
    aa_prec, aa_rec, pep_prec = aa_match_metrics(matches, n_true, n_pred)
    pep_match = np.array([mm[1] for mm in matches])
    order = np.argsort(-m["score"].fillna(-np.inf).to_numpy(), kind="stable")
    n = len(m)
    ex_sorted = exact_il.to_numpy()[order]
    pm_sorted = pep_match[order]
    k = np.arange(1, n + 1)
    curve = pd.DataFrame({
        "coverage": k / n,
        "score": m["score"].to_numpy()[order],
        "precision_exact_il": np.cumsum(ex_sorted) / k,
        "precision_casanovo": np.cumsum(pm_sorted) / k,
    })
    return {
        "n": n, "n_predicted": int(has.sum()),
        "exact": float(exact.mean()), "exact_il": float(exact_il.mean()),
        "casanovo_pep": float(pep_match.mean()),
        "aa_precision": float(aa_prec), "aa_recall": float(aa_rec),
        "curve": curve,
        "per_spectrum": pd.DataFrame({**{k_: m[k_] for k_ in keys}, "true": true_c,
                                      "pred": pred_c, "score": m["score"],
                                      "exact": exact, "exact_il": exact_il,
                                      "casanovo_pep": pep_match}),
    }


def curve_points(curve: pd.DataFrame, coverages=(0.1, 0.25, 0.5, 0.75, 1.0)) -> pd.DataFrame:
    """Precision at selected coverage levels from :func:`evaluate`'s curve."""
    rows = []
    n = len(curve)
    for c in coverages:
        i = max(0, min(n - 1, int(np.ceil(c * n)) - 1))
        r = curve.iloc[i]
        rows.append({"coverage": c, "score_threshold": r["score"],
                     "precision_exact_il": r["precision_exact_il"],
                     "precision_casanovo": r["precision_casanovo"]})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- PEAKS table
def _ptm_string(peptide: str) -> str:
    """PEAKS denovoOnly PTM column: one name per modified site, sorted, ``; ``-joined."""
    from .config import SearchConfig
    from .export_builder import ptm_names
    names = ptm_names(peptide, SearchConfig())[0]
    return "; ".join(sorted(n for n in names.split(";") if n)) if names else ""


class _ProteinIndex:
    """Exact substring lookup of peptides in a FASTA (I and L equivalent)."""

    K = 5

    def __init__(self, fasta: str | Path):
        from pyteomics import fasta as pfasta
        self.names, self.seqs = [], []
        for header, seq in pfasta.read(str(fasta)):
            acc = header.split()[0]
            parts = acc.split("|")
            if len(parts) >= 3 and parts[0] in ("sp", "tr"):
                acc = "|".join(parts[1:3])
            if acc.startswith("#DECOY#") or acc.startswith("rev_"):
                continue
            self.names.append(acc)
            self.seqs.append(seq.replace("I", "L"))
        self.kmers: dict[str, set[int]] = {}
        for j, s in enumerate(self.seqs):
            for i in range(len(s) - self.K + 1):
                self.kmers.setdefault(s[i:i + self.K], set()).add(j)

    def find(self, stripped: str) -> list[str]:
        p = stripped.replace("I", "L")
        cand = self.kmers.get(p[:self.K], set()) if len(p) >= self.K else range(len(self.seqs))
        return [self.names[j] for j in sorted(cand) if p in self.seqs[j]]


def denovo_only_table(denovo_df: pd.DataFrame, accepted_psm_scans: Iterable,
                      source_file: str, precursors: pd.DataFrame | None = None,
                      min_alc: float = 50.0, mode: str = "HCD",
                      precursor_tol_ppm: float | None = 15.0,
                      fasta: str | Path | None = None) -> pd.DataFrame:
    """Rows in the PEAKS ``<Sample>.db.denovoOnly.csv`` layout.

    ``denovo_df``: output of :func:`sequence_mzml` (scan, sequence, per_residue_scores, charge,
    exp_mz; optional feature_id, rt_min / rt (s), area, mz columns are used when present).
    ``accepted_psm_scans``: scans (ints) or (scan, feature_id) pairs that have an accepted DB
    PSM; those spectra are excluded (PEAKS: "de novo only" = no DB hit).
    ``precursors``: optional DataFrame(scan[, feature_id], mz, rt_min (PEAKS RT, minutes) or
    rt (seconds), area) with the precursor/feature attributes to report, joined on
    (scan, feature_id) when both sides have feature_id, else on scan; otherwise the values
    carried by ``denovo_df`` are used (Area empty, Feature Id 0 when absent).
    ALC = mean(per-residue score) * 100; local confidence = per-residue score * 100 (ints);
    rows with ALC < ``min_alc`` are dropped; sorted by ALC descending, as PEAKS does.
    ``precursor_tol_ppm``: drop predictions whose mass is not within this tolerance of the
    reported precursor (PEAKS's denovoOnly rows all satisfy |ppm| <= the precursor tolerance,
    15 ppm in the reference; Casanovo itself does not enforce it). None keeps all.
    ``fasta``: optional database to fill Accession (exact I=L substring matches, ``;``-joined).
    """
    from .peaks_export import canonical_columns
    d = denovo_df.copy()
    d = d[d["sequence"].notna() & (d["sequence"] != "")]
    acc = list(accepted_psm_scans)
    pair_mode = bool(acc) and isinstance(acc[0], (tuple, list))
    if pair_mode and "feature_id" in d.columns:
        excl = {(int(a), int(b)) for a, b in acc}
        mask = [(int(s), int(f)) not in excl for s, f in zip(d["scan"], d["feature_id"])]
    else:
        excl = {int(a[0]) if isinstance(a, (tuple, list)) else int(a) for a in acc}
        mask = [int(s) not in excl for s in d["scan"]]
    d = d[np.asarray(mask, dtype=bool)].copy()
    if precursors is not None:
        p = precursors.rename(columns={"z": "charge"}).copy()
        on = ["scan", "feature_id"] if ("feature_id" in p.columns and "feature_id" in d.columns) else ["scan"]
        if on == ["scan"]:
            p = p.drop_duplicates("scan")
        cols = ["mz", "rt", "rt_min", "area", "feature_id"]
        d = d.drop(columns=[c for c in cols if c in d.columns and c not in on]).merge(
            p[[c for c in p.columns if c in on + cols]], on=on, how="left")
    if "mz" not in d.columns:
        d["mz"] = d["exp_mz"]
    d["mz"] = d["mz"].fillna(d["exp_mz"])
    if "rt_min" not in d.columns:
        d["rt_min"] = d["rt"].astype(float) / 60.0 if "rt" in d.columns else np.nan
    for c, v in (("area", np.nan), ("feature_id", 0)):
        if c not in d.columns:
            d[c] = v
    lc = d["per_residue_scores"].map(lambda a: np.clip(np.asarray(a, dtype=float), 0, 1) * 100)
    d["alc"] = lc.map(lambda a: float(np.mean(a)) if len(a) else 0.0)
    d = d[d["alc"] >= min_alc].copy()
    lc = lc.loc[d.index]
    peptide = d["sequence"].map(from_casanovo)
    stripped = peptide.map(strip_mods)
    mass = d["sequence"].map(peptide_mass)
    z = d["charge"].astype(int)
    obs = (d["mz"].astype(float) - PROTON) * z
    ppm = (obs - mass) / mass * 1e6
    if precursor_tol_ppm is not None:
        keep = (ppm.abs() <= precursor_tol_ppm).to_numpy()
        d, lc, peptide, stripped, mass, z, ppm = (x[keep] for x in (d, lc, peptide, stripped, mass, z, ppm))
    accession = ""
    if fasta is not None:
        idx = _ProteinIndex(fasta)
        accession = stripped.map(lambda s: ";".join(idx.find(s)))
    out = pd.DataFrame({
        "Source File": source_file,
        "Scan": d["scan"].astype(int).to_numpy(),
        "Peptide": peptide.to_numpy(),
        "Tag length": stripped.str.len().to_numpy(),
        "ALC (%)": d["alc"].to_numpy(),
        "Length": stripped.str.len().to_numpy(),
        "m/z": d["mz"].astype(float).to_numpy(),
        "z": z.to_numpy(),
        "RT": d["rt_min"].astype(float).to_numpy(),
        "Area": d["area"].astype(float).to_numpy(),
        "Mass": mass.to_numpy(),
        "ppm": ppm.to_numpy(),
        "PTM": peptide.map(_ptm_string).to_numpy(),
        "local confidence (%)": lc.map(lambda a: " ".join(str(int(round(x))) for x in a)).to_numpy(),
        "mode": mode,
        "tag(>=0%)": peptide.to_numpy(),
        "Accession": accession if isinstance(accession, str) else accession.to_numpy(),
        "Feature Id": d["feature_id"].fillna(0).astype(int).to_numpy(),
    })
    out = out.sort_values(["ALC (%)", "Scan"], ascending=[False, True], kind="stable")
    return out[canonical_columns("denovo_only")].reset_index(drop=True)


# --------------------------------------------------------------------------- training helpers
def finetune_checkpoint(src: str | Path, dst: str | Path) -> Path:
    """Copy a Casanovo checkpoint so ``casanovo train --model`` can fine-tune from it.

    The released casanovo_orbitrap_v5-2-0.ckpt stores the DB-search hyperparameters
    ``precursor_mass_tol`` and ``isotope_error_range``; Casanovo 5.2.1's Spec2Pep passes unknown
    hyperparameters to ``torch.optim.Adam`` and training fails with
    ``TypeError: Adam.__init__() got an unexpected keyword argument 'precursor_mass_tol'``.
    This drops those keys (and the optimizer state, unused without ``--load_all_states``);
    the weights are unchanged.
    """
    import torch
    ck = torch.load(str(src), map_location="cpu", weights_only=False)
    for k in ("precursor_mass_tol", "isotope_error_range"):
        ck["hyper_parameters"].pop(k, None)
    ck["optimizer_states"] = []
    ck["lr_schedulers"] = []
    torch.save(ck, str(dst))
    return Path(dst)
