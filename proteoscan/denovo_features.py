"""De novo agreement features for rescoring, imitating PEAKS DB's CAA score and Tag Length.

PEAKS DB scores every candidate partly by its agreement with the de novo sequence of the
spectrum (docs/PEAKS_ALGORITHMS.md sections 2.1 and 6): residues whose de novo local confidence is
below 30 % are collapsed into mass segments, the remaining tag is aligned to the candidate so
that the mass of each aligned block is equal on both sides, and the CAA score is the number of
common amino acids. The DB exports carry ``Tag Length`` and ``CAA (%)``. This module computes
the same kind of quantities from our Casanovo calls (docs/DENOVO.md, docs/DENOVO_FEATURES.md).

Main entry points

* :func:`run_denovo_queries`: Casanovo on one spectrum per (scan, precursor) query (the
  pipeline's expanded mzML), sharded over several processes/GPUs; returns every beam candidate
  with the precursor-consistent pick flagged.
* :func:`denovo_features`: per-candidate features against that table (vectorised; numba).
* :func:`add_denovo_features`: convenience wrapper that joins them onto a candidate table.

Feature definitions (candidate c, de novo call d of the same query, confidence threshold 0.3,
mass tolerance 0.02 Da, I = L, modifications included as mass and identity):

``caa``               confident residues of d that align to the same residue of c: the residue's
                      start and end cumulative masses both coincide with residue boundaries of c,
                      counted from the N terminus (prefix masses) or from the C terminus (suffix
                      masses), so a mass offset inside the peptide (wrong precursor, isotope
                      error, unknown modification) costs only the residues around it.
``caa_pct``           100 * caa / len(c)  (PEAKS ``CAA (%)`` scale).
``tag_len``           longest run of consecutive aligned confident residues of d.
``max_caa_any_beam``  max caa over all beam candidates of the query.
``denovo_conf_mean``  mean de novo confidence of the aligned residues (0 if none).
``denovo_exact``      1 if the chosen de novo call equals c residue by residue (I = L, mods as
                      rounded mass deltas, N-terminal mods folded into the first residue).

Queries without a de novo call (charge > 4, fewer than ``min_peaks`` peaks, ...) get NaN.
"""
from __future__ import annotations

import logging
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .seqnorm import parse

log = logging.getLogger(__name__)

DENOVO_FEATURE_COLUMNS = ["caa", "caa_pct", "tag_len", "max_caa_any_beam", "denovo_conf_mean",
                          "denovo_exact"]

MONO = {
    "G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414,
    "T": 101.047670, "C": 103.009185, "L": 113.084064, "I": 113.084064, "N": 114.042927,
    "D": 115.026943, "Q": 128.058578, "K": 128.094963, "E": 129.042593, "M": 131.040485,
    "H": 137.058912, "F": 147.068414, "R": 156.101111, "Y": 163.063329, "W": 186.079313,
    "U": 150.953636, "O": 237.147727,
}


# --------------------------------------------------------------------------- encoding
class _Encoder:
    """Peptide -> per-residue (identity code, mass). Identity = residue letter (I -> L) plus the
    modification delta rounded to 0.01 Da; an N-terminal modification is folded into residue 1
    (both its mass and its identity), a C-terminal one into the last residue."""

    def __init__(self):
        self.mod_ids: dict[int, int] = {0: 0}
        self.cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def _code(self, aa: str, delta: float) -> int:
        key = int(round(delta * 100))
        mid = self.mod_ids.setdefault(key, len(self.mod_ids))
        letter = ord("L" if aa == "I" else aa) - 65
        return letter + 32 * mid

    def encode(self, seq: str):
        hit = self.cache.get(seq)
        if hit is not None:
            return hit
        stripped, mods = parse(seq)
        n = len(stripped)
        md = {}
        for p, dlt in mods:
            p = min(max(p, 1), n)  # fold N-/C-terminal mods into the terminal residue
            md[p] = md.get(p, 0.0) + dlt
        codes = np.empty(n, dtype=np.int32)
        masses = np.empty(n, dtype=np.float64)
        for k, aa in enumerate(stripped, start=1):
            dlt = md.get(k, 0.0)
            codes[k - 1] = self._code(aa, dlt)
            masses[k - 1] = MONO.get(aa, 0.0) + dlt
        self.cache[seq] = (codes, masses)
        return codes, masses

    def encode_many(self, seqs: Sequence[str], max_len: int | None = None):
        """Padded arrays (codes -1 beyond the length), lengths."""
        uniq, inv = np.unique(np.asarray(seqs, dtype=object).astype(str), return_inverse=True)
        enc = [self.encode(s) for s in uniq]
        lens_u = np.array([len(c) for c, _ in enc], dtype=np.int32)
        L = int(lens_u.max()) if len(lens_u) else 1
        if max_len is not None:
            L = min(L, max_len)
        codes_u = np.full((len(uniq), L), -1, dtype=np.int32)
        mass_u = np.zeros((len(uniq), L), dtype=np.float64)
        for i, (c, m) in enumerate(enc):
            k = min(len(c), L)
            codes_u[i, :k] = c[:k]
            mass_u[i, :k] = m[:k]
        lens_u = np.minimum(lens_u, L)
        return codes_u[inv], mass_u[inv], lens_u[inv]


def _residue_confidences(sequence: str, scores, n_res: int) -> np.ndarray:
    """Per-residue confidences aligned to residues. Casanovo emits an N-terminal modification
    as its own token with its own score; that score is merged (min) into residue 1."""
    a = np.asarray(scores, dtype=np.float64)
    if len(a) == n_res:
        return a
    if len(a) == n_res + 1:
        out = a[1:].copy()
        out[0] = min(out[0], a[0])
        return out
    out = np.zeros(n_res)
    k = min(len(a), n_res)
    out[:k] = a[:k]
    return out


# --------------------------------------------------------------------------- kernel
try:
    import numba as _nb
    _njit = _nb.njit(cache=True, parallel=True, fastmath=False)
    _prange = _nb.prange
except Exception:  # pragma: no cover - numba is in the env; slow pure-python fallback
    _njit = lambda f: f  # noqa: E731
    _prange = range


@_njit
def _kernel(c_codes, c_mass, c_len, beam_start, beam_n, pick_off,
            d_codes, d_mass, d_conf, d_len, tol, thr,
            o_caa, o_tag, o_conf, o_max, o_exact):
    n = c_codes.shape[0]
    Lc = c_codes.shape[1]
    Ld = d_codes.shape[1]
    for r in _prange(n):
        b0 = beam_start[r]
        if b0 < 0:
            o_caa[r] = np.nan; o_tag[r] = np.nan; o_conf[r] = np.nan
            o_max[r] = np.nan; o_exact[r] = np.nan
            continue
        lc = c_len[r]
        pc = np.zeros(Lc + 1)
        for i in range(lc):
            pc[i + 1] = pc[i] + c_mass[r, i]
        tot_c = pc[lc]
        best_any = 0
        mN = np.empty(Ld, dtype=np.int64)
        mC = np.empty(Ld, dtype=np.int64)
        sufC = np.empty(Ld + 1, dtype=np.int64)
        for b in range(beam_n[r]):
            q = b0 + b
            ld = d_len[q]
            pd_ = np.zeros(Ld + 1)
            for j in range(ld):
                pd_[j + 1] = pd_[j] + d_mass[q, j]
            tot_d = pd_[ld]
            # candidate residue aligned to each confident de novo residue, anchored at the
            # N terminus (equal prefix mass) and at the C terminus (equal suffix mass)
            i_n = 0
            i_c = 0
            for j in range(ld):
                mN[j] = -1
                mC[j] = -1
                if d_conf[q, j] < thr:
                    continue
                code = d_codes[q, j]
                start = pd_[j]
                while i_n < lc and pc[i_n] < start - tol:
                    i_n += 1
                if i_n < lc and abs(pc[i_n] - start) <= tol and c_codes[r, i_n] == code:
                    mN[j] = i_n
                suf = tot_d - pd_[j + 1]
                while i_c < lc and (tot_c - pc[i_c + 1]) > suf + tol:
                    i_c += 1
                if i_c < lc and abs((tot_c - pc[i_c + 1]) - suf) <= tol and c_codes[r, i_c] == code:
                    mC[j] = i_c
            # one monotone alignment: N-anchored matches before a split point, C-anchored
            # matches after it, each candidate residue used once
            sufC[ld] = 0
            for j in range(ld - 1, -1, -1):
                sufC[j] = sufC[j + 1] + (1 if mC[j] >= 0 else 0)
            nN = 0; maxN = -1; p = 0
            best = -1; best_s = 0; best_p = 0
            for s_ in range(ld + 1):
                if s_ > 0 and mN[s_ - 1] >= 0:
                    nN += 1
                    maxN = mN[s_ - 1]
                while p < ld and mC[p] <= maxN:
                    p += 1
                k = s_ if s_ > p else p
                tot = nN + sufC[k]
                if tot > best:
                    best = tot; best_s = s_; best_p = k
            caa = 0; run = 0; tag = 0; csum = 0.0
            for j in range(ld):
                ok = (j < best_s and mN[j] >= 0) or (j >= best_p and mC[j] >= 0)
                if ok:
                    caa += 1; run += 1; csum += d_conf[q, j]
                    if run > tag:
                        tag = run
                else:
                    run = 0
            if caa > best_any:
                best_any = caa
            if b == pick_off[r]:
                o_caa[r] = caa; o_tag[r] = tag
                o_conf[r] = csum / caa if caa > 0 else 0.0
                ex = 1.0 if ld == lc else 0.0
                if ex == 1.0:
                    for j in range(ld):
                        if d_codes[q, j] != c_codes[r, j]:
                            ex = 0.0
                            break
                o_exact[r] = ex
        o_max[r] = best_any


# --------------------------------------------------------------------------- features
def denovo_features(cands: pd.DataFrame, denovo_all: pd.DataFrame, key: str = "query_id",
                    tol: float = 0.02, conf_threshold: float = 0.3,
                    threads: int | None = 16) -> pd.DataFrame:
    """Per-candidate de novo agreement features.

    cands : columns ``key`` and ``peptide_sage`` (any seqnorm-parsable notation).
    denovo_all : one row per de novo beam candidate: ``key``, ``sequence`` (Casanovo ProForma or
        any parsable notation), ``per_residue_scores`` (array per row), and ``picked`` (bool,
        the precursor-consistent choice; if absent, the first row of each key is used).
    Returns a DataFrame indexed like ``cands`` with :data:`DENOVO_FEATURE_COLUMNS`.
    """
    t0 = time.time()
    if threads is not None:
        try:
            import numba
            numba.set_num_threads(max(1, min(threads, numba.config.NUMBA_NUM_THREADS)))
        except Exception:
            pass
    d = denovo_all[denovo_all["sequence"].notna() & (denovo_all["sequence"].astype(str) != "")]
    d = d.copy()
    if "picked" not in d.columns:
        d["picked"] = ~d.duplicated(key)
    d = d.sort_values([key, "picked"], ascending=[True, False], kind="stable").reset_index(drop=True)
    enc = _Encoder()
    d_codes, d_mass, d_len = enc.encode_many(d["sequence"].tolist())
    Ld = d_codes.shape[1]
    d_conf = np.zeros(d_codes.shape, dtype=np.float64)
    for i, (s, sc, n) in enumerate(zip(d["sequence"].to_numpy(), d["per_residue_scores"].to_numpy(), d_len)):
        c = _residue_confidences(s, sc if sc is not None else [], int(n))
        d_conf[i, :n] = c[:n]
    keys = d[key].to_numpy()
    uk, first, counts = np.unique(keys, return_index=True, return_counts=True)
    # the picked beam is first within its key after the sort above
    ck = cands[key].to_numpy()
    pos = np.searchsorted(uk, ck)
    pos_c = np.clip(pos, 0, max(len(uk) - 1, 0))
    found = (len(uk) > 0) & (pos < len(uk)) & (uk[pos_c] == ck) if len(uk) else np.zeros(len(ck), bool)
    beam_start = np.where(found, first[pos_c], -1).astype(np.int64)
    beam_n = np.where(found, counts[pos_c], 0).astype(np.int64)
    pick_off = np.zeros(len(ck), dtype=np.int64)
    c_codes, c_mass, c_len = enc.encode_many(cands["peptide_sage"].tolist())
    n = len(cands)
    out = {c: np.empty(n) for c in ("caa", "tag", "conf", "max", "exact")}
    t1 = time.time()
    _kernel(c_codes, c_mass, c_len.astype(np.int64), beam_start, beam_n, pick_off,
            d_codes, d_mass, d_conf, d_len.astype(np.int64), float(tol), float(conf_threshold),
            out["caa"], out["tag"], out["conf"], out["max"], out["exact"])
    log.info("de novo features: %d candidates, %d de novo rows (%d queries); encode %.1fs, kernel %.1fs",
             n, len(d), len(uk), t1 - t0, time.time() - t1)
    res = pd.DataFrame({
        "caa": out["caa"],
        "caa_pct": 100.0 * out["caa"] / np.maximum(c_len, 1),
        "tag_len": out["tag"],
        "max_caa_any_beam": out["max"],
        "denovo_conf_mean": out["conf"],
        "denovo_exact": out["exact"],
    }, index=cands.index)
    return res


def add_denovo_features(cands: pd.DataFrame, denovo_all: pd.DataFrame, key: str = "query_id",
                        **kw) -> pd.DataFrame:
    """Copy of ``cands`` with :data:`DENOVO_FEATURE_COLUMNS` appended (existing ones replaced)."""
    f = denovo_features(cands, denovo_all, key=key, **kw)
    out = cands.drop(columns=[c for c in DENOVO_FEATURE_COLUMNS if c in cands.columns])
    for c in f.columns:
        out[c] = f[c].to_numpy()
    return out


# --------------------------------------------------------------------------- running Casanovo
def _run_shard(args):
    cmd, env, log_path = args
    with open(log_path, "w") as fh:
        t = time.time()
        subprocess.run(cmd, check=True, env=env, stdout=fh, stderr=subprocess.STDOUT)
        return time.time() - t


def run_denovo_queries(mzml: str | Path, queries: pd.DataFrame, model_ckpt: str | Path,
                       out_dir: str | Path, stem: str, n_shards: int = 4,
                       gpus: Sequence[str] = ("2", "3"), threads_per_shard: int = 4,
                       n_beams: int = 5, top_match: int = 5, max_charge: int = 4,
                       max_peptide_len: int = 50, batch_size: int = 512,
                       precursor_tol_ppm: float = 15.0, isotopes: Iterable[int] = (0, 1),
                       precision: str = "32-true", casanovo_bin: str | None = None, reuse: bool = True) -> pd.DataFrame:
    """De novo sequence one spectrum per query with Casanovo, sharded over processes and GPUs.

    mzml : spectra addressable by ``scan=<id>`` native ids. For the pipeline's expanded mzML the
        scan number of each spectrum is the ``query_id`` (one spectrum per (scan, precursor)).
    queries : DataFrame(id, mz, charge[, scan]) -- ``id`` is looked up in ``mzml``; mz/charge
        override the precursor written to the MGF. Other columns are ignored.
    Shards are written as MGFs (SCANS = id), run as independent ``casanovo sequence`` processes
    (shard k on ``gpus[k % len(gpus)]``, ``threads_per_shard`` CPU threads each) and merged.
    ``casanovo_bin`` None runs :mod:`proteoscan.casanovo_fastio` (same predictions, mzTab
    written without Casanovo's slow per-PSM ProForma parsing). ``precision="16-mixed"`` halves
    GPU time; on 2,000 Jurkat spectra 95 % of picked calls are identical to fp32 (98 % of the
    top-scoring half) and the precursor-consistent fraction is unchanged.

    Returns every beam candidate: id, rank, sequence (ProForma), peptide (PEAKS notation), score,
    per_residue_scores, charge, exp_mz, ppm (best of isotopes), precursor_ok, picked (the
    highest-scoring candidate within ``precursor_tol_ppm``, else the top one), shard.
    """
    from . import denovo as dn
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    q = queries[(queries["charge"] >= 1) & (queries["charge"] <= max_charge)].copy()
    q = q.sort_values("id").reset_index(drop=True)
    shard_of = np.arange(len(q)) % n_shards
    jobs = []
    t0 = time.time()
    mgfs = []
    for k in range(n_shards):
        mgf = out_dir / f"{stem}.shard{k}.mgf"
        rows = out_dir / f"{stem}.shard{k}.rows.tsv"
        mgfs.append((mgf, rows, q[shard_of == k]))
    todo = [(m, r, e) for m, r, e in mgfs if not (reuse and m.exists() and r.exists())]
    if todo:
        # MGF writing is I/O + pyopenms bound: one worker per shard (processes)
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=min(len(todo), threads_per_shard * len(gpus))) as ex:
            futs = [ex.submit(_write_shard, str(mzml), str(m), str(r), e[["id", "mz", "charge"]])
                    for m, r, e in todo]
            for f in futs:
                f.result()
    log.info("MGFs written in %.0fs", time.time() - t0)
    for k, (mgf, rows, _) in enumerate(mgfs):
        shard_stem = f"{stem}.shard{k}"
        result = out_dir / f"{shard_stem}.mztab"
        if reuse and result.exists():
            continue
        cfg = out_dir / f"{shard_stem}.config.yaml"
        dn.write_config(cfg, n_beams=n_beams, top_match=top_match, predict_batch_size=batch_size,
                        max_peptide_len=max_peptide_len, max_charge=max_charge,
                        precision=precision)
        env = dict(os.environ)
        env["CUDA_VISIBLE_DEVICES"] = str(gpus[k % len(gpus)])
        for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[v] = str(threads_per_shard)
        if casanovo_bin is None:  # the casanovo CLI with fast mzTab writing (casanovo_fastio)
            launcher = [sys.executable, "-m", "proteoscan.casanovo_fastio"]
            env["PYTHONPATH"] = os.pathsep.join(
                [str(Path(__file__).resolve().parents[1])] + [x for x in [env.get("PYTHONPATH")] if x])
        else:
            launcher = [casanovo_bin]
        cmd = launcher + ["sequence", "-m", str(model_ckpt), "-c", str(cfg), "-d", str(out_dir),
               "-o", shard_stem, "-f", "-v", "info", str(mgf)]
        jobs.append((cmd, env, str(out_dir / f"{shard_stem}.stdout.log")))
    if jobs:
        t1 = time.time()
        with ThreadPoolExecutor(max_workers=len(jobs)) as ex:
            times = list(ex.map(_run_shard, jobs))
        log.info("Casanovo: %d shards in %.0fs (per shard %s)", len(jobs), time.time() - t1,
                 [round(t) for t in times])
    parts = []
    for k, (mgf, rows, _) in enumerate(mgfs):
        df = dn.read_mztab(out_dir / f"{stem}.shard{k}.mztab", mgf=mgf)
        r = pd.read_csv(rows, sep="\t")
        df = df.drop(columns=["scan"]).merge(r[["spectra_ref", "id"]], on="spectra_ref", how="left")
        df["shard"] = k
        parts.append(df)
    d = pd.concat(parts, ignore_index=True)
    iso = tuple(isotopes)
    d["ppm"] = [dn.precursor_ppm(s, m, int(z), iso) if isinstance(s, str) and s else np.inf
                for s, m, z in zip(d["sequence"], d["exp_mz"], d["charge"])]
    d["precursor_ok"] = d["ppm"] <= precursor_tol_ppm
    order = d.sort_values(["id", "precursor_ok", "score"], ascending=[True, False, False], kind="stable")
    d["picked"] = False
    d.loc[order.drop_duplicates("id").index, "picked"] = True
    cols = ["id", "rank", "sequence", "peptide", "score", "per_residue_scores", "charge", "exp_mz",
            "ppm", "precursor_ok", "picked", "shard"]
    return d[cols].sort_values(["id", "rank"]).reset_index(drop=True)


def _write_shard(mzml: str, mgf: str, rows_tsv: str, entries: pd.DataFrame):
    from . import denovo as dn
    ent = entries.rename(columns={"id": "scan"}).reset_index(drop=True)
    written = dn.write_annotated_mgf(ent, mzml, mgf + ".tmp", return_rows=True)
    ent.iloc[written].rename(columns={"scan": "id"}).assign(
        spectra_ref=[f"index={k}" for k in range(len(written))]).to_csv(rows_tsv, sep="\t", index=False)
    os.replace(mgf + ".tmp", mgf)
    return len(written)
