"""PEAKS-style HTML export pages for a PEAKS-format DB-search export (`summary`).

Usage::

    python -m proteoscan.summary --export <export dir> [--mzml run.mzML ...] [--fasta db.fasta]
        [--psms-scored psms_scored.parquet] [--features ms1_features.parquet ...]
        [--precursors precursors.parquet ...] [--config search.yml] [--mgf run.mgf] --out <dir>

Writes ``db.summary.html``, ``db.proteins.html``, ``db.peptides.html``, ``export.css`` and an
``img/`` folder of PNG figures into ``--out``.  The pages follow the structure of PEAKS
Studio's HTML export (sections "1. Notes" to "4. Other Information", Figures 1-6 and
Tables 1-6 with PEAKS's captions, one block per protein with a coverage map, the peptide
list) so they read the same way.  Stylesheet and figures are our own (matplotlib).

Inputs other than the export are optional; every figure or number that needs a missing
input is omitted with a note in "1. Notes".  When the export sits inside a pipeline run
directory (``<run>/export``), ``psms_scored.parquet``, ``<run>/<stem>/precursors.parquet``,
``<run>/<stem>/ms1_features.parquet`` and ``<run>/search_config.yml`` are found
automatically, and a ``*.mgf`` in the export directory (PEAKS writes one) supplies the
chimeric-scan count.  See docs/UI.md.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import glob
import html
import io
import json
import math
import os
import re
import shutil
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from . import __version__
from . import seqnorm
from .fdr import fdr_curve, threshold_summary
from .peaks_export import read_export
from .proteins import Enzyme, group_proteins, protein_score, read_protein_db

__all__ = ["build_summary", "annotate_fragments", "compute_statistics", "main"]

PROTON = 1.007276466
H2O = 18.010565
MONO = {
    "G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414,
    "T": 101.047679, "C": 103.009185, "L": 113.084064, "I": 113.084064, "N": 114.042927,
    "D": 115.026943, "Q": 128.058578, "K": 128.094963, "E": 129.042593, "M": 131.040485,
    "H": 137.058912, "F": 147.068414, "R": 156.101111, "Y": 163.063329, "W": 186.079313,
    "U": 150.953636, "O": 237.147727,
}

# validated categorical slots (dataviz reference palette, light mode, first three all-pairs safe)
C_BLUE, C_ORANGE, C_AQUA = "#2a78d6", "#eb6834", "#1baf7a"
C_VIOLET, C_MAGENTA, C_YELLOW, C_GREEN, C_RED = "#6250d6", "#e87ba4", "#eda100", "#008300", "#e34948"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
B_COLOR, Y_COLOR = C_BLUE, C_RED
MOD_COLORS = [C_ORANGE, C_AQUA, C_VIOLET, C_MAGENTA, C_YELLOW, C_GREEN]
NA = "N/A"


def _log(msg):
    print(f"[summary] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------


def _isna(v):
    if v is None or v is pd.NA:
        return True
    if isinstance(v, str):
        return v == ""
    try:
        return bool(np.isnan(v))
    except (TypeError, ValueError):
        return False


def _trim(v, n):
    """Fixed n decimals with trailing zeros removed (PEAKS HTML style: 44.8, 955.545)."""
    if _isna(v):
        return ""
    s = f"{float(v):.{n}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def _fixed(v, n):
    return "" if _isna(v) else f"{float(v):.{n}f}"


def _area(v):
    """PEAKS HTML area: 4.97e6, 2e5; missing -> 0."""
    if _isna(v) or float(v) == 0:
        return "0"
    m, e = f"{float(v):.2e}".split("e")
    m = m.rstrip("0").rstrip(".")
    return f"{m}e{int(e)}"


def _dec2(v):
    """Java DecimalFormat("#.##"): .98, 15.99 (as PEAKS prints the PTM mass)."""
    s = _trim(v, 2)
    return s[1:] if s.startswith("0.") else ("-" + s[2:] if s.startswith("-0.") else s)


def _int(v):
    return "" if _isna(v) else str(int(v))


def _esc(s):
    return html.escape("" if s is None else str(s), quote=True)


def _pct(x, nd=1):
    return NA if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{100.0 * x:.{nd}f}%"


def _norm_source(s):
    b = os.path.basename(str(s))
    for ext in (".raw", ".RAW", ".mzML", ".mzml", ".mgf", ".d"):
        if b.endswith(ext):
            return b[: -len(ext)]
    return b


def _truthy(v):
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    return str(v).strip().lower() in ("true", "1", "y", "yes")


# --------------------------------------------------------------------------
# input discovery
# --------------------------------------------------------------------------


def resolve_inputs(export_dir, run_dir=None, psms_scored=None, features=None, precursors=None,
                   config=None, mgf=None):
    """Fill in optional inputs from a pipeline run directory / the export directory."""
    export_dir = os.path.abspath(export_dir)
    if run_dir is None:
        parent = os.path.dirname(export_dir)
        if any(os.path.exists(os.path.join(parent, f)) for f in ("psms_scored.parquet", "run_stats.json")):
            run_dir = parent
    r = {"run_dir": run_dir, "psms_scored": psms_scored, "features": list(features or []),
         "precursors": list(precursors or []), "config": config, "mgf": mgf, "run_stats": None}
    if run_dir:
        if r["psms_scored"] is None and os.path.exists(os.path.join(run_dir, "psms_scored.parquet")):
            r["psms_scored"] = os.path.join(run_dir, "psms_scored.parquet")
        if not r["features"]:
            r["features"] = sorted(glob.glob(os.path.join(glob.escape(run_dir), "*", "ms1_features.parquet")))
        if not r["precursors"]:
            r["precursors"] = sorted(glob.glob(os.path.join(glob.escape(run_dir), "*", "precursors.parquet")))
        if r["config"] is None and os.path.exists(os.path.join(run_dir, "search_config.yml")):
            r["config"] = os.path.join(run_dir, "search_config.yml")
        if os.path.exists(os.path.join(run_dir, "run_stats.json")):
            r["run_stats"] = os.path.join(run_dir, "run_stats.json")
    if r["mgf"] is None:
        hits = sorted(glob.glob(os.path.join(glob.escape(export_dir), "*.mgf")))
        r["mgf"] = hits if hits else None
    elif isinstance(r["mgf"], str):
        r["mgf"] = [r["mgf"]]
    return r


# --------------------------------------------------------------------------
# search parameters (Table 6) and filters (Table 2)
# --------------------------------------------------------------------------


def _count_fasta(path):
    n = 0
    try:
        with open(path, "rb") as fh:
            for line in fh:
                if line.startswith(b">") and not line.startswith(b">#DECOY#"):
                    n += 1
    except OSError:
        return None
    return n


def params_from_config(cfg, fasta=None):
    """(Table 6 lines, filter dict, enzyme params) from a :class:`SearchConfig`."""
    db = cfg.database or fasta or ""
    n_entries = _count_fasta(db) if db and os.path.exists(db) else None
    fixed = [f"{m.name} ({m.delta:+.2f})" for m in cfg.fixed_modifications]
    var = [f"{m.name} ({m.residues}) ({m.delta:+.2f})" for m in cfg.variable_modifications]
    lines = [
        ("Search Engine Name", f"proteoscan {__version__}"),
        ("Parent Mass Error Tolerance", f"{cfg.parent_mass_error_tolerance_ppm:.1f} ppm"),
        ("Fragment Mass Error Tolerance", f"{cfg.fragment_mass_error_tolerance_da:g} Da"),
        ("Precursor Mass Search Type", cfg.precursor_mass_search_type),
        ("Enzyme", cfg.enzyme),
        ("Max Missed Cleavages", str(cfg.max_missed_cleavages)),
        ("Digest Mode", cfg.digest_mode),
        ("Peptide Length Range", f"{cfg.peptide_length_min} - {cfg.peptide_length_max}"),
        ("Fixed Modifications", fixed),
        ("Variable Modifications", var),
        ("Max Variable PTM Per Peptide", str(cfg.max_variable_ptm_per_peptide)),
        ("Database", os.path.splitext(os.path.basename(db))[0] if db else NA),
        ("Taxon", "all species"),
        ("Searched Entries", str(n_entries) if n_entries is not None else NA),
        ("Deep Learning Boost", "Yes" if cfg.deep_learning_boost else "No"),
        ("FDR Estimation", "Enabled" if cfg.fdr_estimation else "Disabled"),
        ("Confident Amino Acid Threshold", f"{cfg.confident_amino_acid_threshold_pct:.2f}%"),
        ("Find Unspecified PTMs with PEAKS PTM", "No"),
        ("Find Mutations and Sequence Variants with SPIDER", "No"),
    ]
    filters = {"peptide_fdr": cfg.peptide_fdr, "protein_min_10lgp": cfg.protein_score_min,
               "protein_min_unique_peptides": cfg.protein_unique_peptides_min,
               "denovo_min_alc_pct": cfg.denovo_alc_min_pct,
               "fragment_tolerance": cfg.fragment_mass_error_tolerance_da, "fragment_tolerance_unit": "Da",
               "precursor_tolerance": cfg.parent_mass_error_tolerance_ppm,
               "fixed_deltas": {f"{m.delta:+.2f}" for m in cfg.fixed_modifications}}
    enz = {"enzyme": cfg.enzyme, "digest_mode": cfg.digest_mode, "max_missed_cleavages": cfg.max_missed_cleavages}
    return lines, filters, enz


def params_from_export(params):
    """Same triple from ``parse_summary_html`` output (an export that has a db.summary.html)."""
    raw6 = (params.get("raw") or {}).get("search_parameters") or {}
    lines = [(k, v) for k, v in raw6.items()]
    fixed = {f"{m['delta']:+.2f}" for m in params.get("fixed_modifications") or [] if m.get("delta") is not None}
    filters = {"peptide_fdr": 0.01,
               "peptide_min_10lgp": params.get("peptide_min_10lgp"),
               "protein_min_10lgp": params.get("protein_min_10lgp", 15),
               "protein_min_unique_peptides": params.get("protein_min_unique_peptides", 1),
               "denovo_min_alc_pct": params.get("denovo_min_alc_pct", 50),
               "fragment_tolerance": params.get("fragment_tolerance", 0.5),
               "fragment_tolerance_unit": params.get("fragment_tolerance_unit", "Da"),
               "precursor_tolerance": params.get("precursor_tolerance"),
               "fixed_deltas": fixed}
    return lines, filters, params


# --------------------------------------------------------------------------
# raw data counts (mzML, precursor list, MGF, feature table)
# --------------------------------------------------------------------------

_SCAN_RE = re.compile(r"\bscan=(\d+)\b")


class MzmlSpectra:
    """MS level counts and random access to spectra by raw scan number (pyopenms)."""

    def __init__(self, path):
        import pyopenms as oms
        self.path = path
        self.od = oms.OnDiscMSExperiment()
        if not self.od.openFile(path):
            raise IOError(f"pyopenms could not open {path} (needs an indexed mzML)")
        meta = self.od.getMetaData()
        self.levels = Counter()
        self.index = {}
        self.rt = {}
        for i in range(meta.getNrSpectra()):
            s = meta.getSpectrum(i)
            lvl = s.getMSLevel()
            self.levels[lvl] += 1
            m = _SCAN_RE.search(s.getNativeID())
            sc = int(m.group(1)) if m else i + 1
            self.index[sc] = i
            self.rt[sc] = s.getRT() / 60.0

    def peaks(self, scan):
        i = self.index.get(int(scan))
        if i is None:
            return None
        s = self.od.getSpectrum(i)
        mz, it = s.get_peaks()
        mz = np.asarray(mz, dtype=np.float64)
        it = np.asarray(it, dtype=np.float64)
        keep = it > 0
        return mz[keep], it[keep]


def mgf_scan_counts(path):
    """(#distinct RAWSCANS, #scans with more than one MGF entry) of a PEAKS MGF."""
    c = Counter()
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b"RAWSCANS="):
                c[line[9:].strip()] += 1
    return len(c), sum(1 for v in c.values() if v > 1)


def precursor_scan_counts(paths):
    n_scans = n_chim = 0
    for p in paths:
        e = pd.read_parquet(p, columns=["scan"])
        vc = e["scan"].value_counts()
        n_scans += len(vc)
        n_chim += int((vc > 1).sum())
    return n_scans, n_chim


# --------------------------------------------------------------------------
# FDR from scored targets + decoys
# --------------------------------------------------------------------------


def _decoy_protein_groups(best, protein_min, min_unique):
    """Protein groups passing the protein filter, built separately for target and decoy
    sequences with the PEAKS rules (score = sum_k s_k/k over sequences, greedy grouping by
    sequence count, >= min_unique group-unique sequences).  Returns (#target, #decoy)."""
    out = {}
    for decoy in (False, True):
        sub = best[best["is_decoy"] == decoy]
        fam_score = dict(zip(sub["seqkey"], sub["minus10lgp"]))
        prot_fams = defaultdict(set)
        for fam, prots in zip(sub["seqkey"], sub["proteins"].astype(str)):
            for p in prots.split(";"):
                p = p.strip()
                if not p or p.startswith("#DECOY#") != decoy:
                    continue
                prot_fams[p].add(fam)
        if not prot_fams:
            out[decoy] = 0
            continue
        order = sorted(prot_fams, key=lambda p: (-len(prot_fams[p]), p))
        dom = {p: frozenset(f) for p, f in prot_fams.items()}
        group_of, groups, _top = group_proteins(prot_fams, dom, order)
        fam_groups = defaultdict(set)
        for p, g in group_of.items():
            for f in prot_fams[p]:
                fam_groups[f].add(g)
        n = 0
        for g, members in enumerate(groups):
            lead = members[0]
            sc = protein_score(fam_score[f] for f in prot_fams[lead])
            uniq = sum(1 for f in prot_fams[lead] if len(fam_groups[f]) == 1)
            if sc >= protein_min and uniq >= min_unique:
                n += 1
        out[decoy] = n
    return out[False], out[True]


def scored_statistics(path, peptide_fdr=0.01, protein_min=15.0, min_unique=1, psm_floor=20.0):
    """Threshold, FDR curves and Table 3 FDRs from ``psms_scored.parquet`` (best candidate per
    query, targets and decoys, column ``minus10lgp``).  Mirrors ``pipeline.peaks_filter``."""
    want = ["source_file", "scan", "charge", "peptide_sage", "stripped", "proteins", "is_decoy",
            "minus10lgp", "precursor_ppm"]
    import pyarrow.parquet as pq
    have = set(pq.ParquetFile(path).schema.names)
    t = pd.read_parquet(path, columns=[c for c in want if c in have])
    t = t[t["minus10lgp"].notna()].copy()
    t["is_decoy"] = t["is_decoy"].astype(bool)
    if "stripped" not in t:
        t["stripped"] = t["peptide_sage"].map(seqnorm.strip_mods)
    t["seqkey"] = t["stripped"].str.replace("I", "L")
    best = t.sort_values("minus10lgp", ascending=False, kind="mergesort").drop_duplicates("seqkey")
    summ = threshold_summary(best["minus10lgp"], best["is_decoy"], peptide_fdr)
    thr = summ["threshold"]
    ok = set(best.loc[best["minus10lgp"] >= thr, "seqkey"])
    floor = min(thr, psm_floor)
    acc = t[t["seqkey"].isin(ok) & (t["minus10lgp"] >= floor)]
    keys = [c for c in ("source_file", "scan", "peptide_sage", "charge") if c in t]
    acc = acc.sort_values("minus10lgp", ascending=False).drop_duplicates(keys)
    n_t, n_d = int((~acc["is_decoy"]).sum()), int(acc["is_decoy"].sum())
    # PSM unit for the curve and plots: one row per (file, scan, peptide, charge), as exported
    tq = t.sort_values("minus10lgp", ascending=False, kind="mergesort").drop_duplicates(keys)
    passing = best[best["minus10lgp"] >= thr]
    g_t, g_d = _decoy_protein_groups(passing, protein_min, min_unique) if math.isfinite(thr) else (0, 0)
    return {
        "threshold": thr if math.isfinite(thr) else None,
        "n_target_seq": summ["n_target"], "n_decoy_seq": summ["n_decoy"],
        "fdr_seq": summ["fdr"] if math.isfinite(thr) else None,
        "n_target_psm": n_t, "n_decoy_psm": n_d, "fdr_psm": n_d / n_t if n_t else None,
        "n_target_groups": g_t, "n_decoy_groups": g_d, "fdr_group": g_d / g_t if g_t else None,
        "curve_seq": fdr_curve(best["minus10lgp"], best["is_decoy"]),
        "curve_psm": fdr_curve(tq["minus10lgp"], tq["is_decoy"]),
        "scores": tq["minus10lgp"].to_numpy(), "decoy": tq["is_decoy"].to_numpy(),
        "ppm": tq["precursor_ppm"].to_numpy() if "precursor_ppm" in tq else None,
        "n_queries": len(t),
    }


# --------------------------------------------------------------------------
# Table 1 / Table 5 statistics
# --------------------------------------------------------------------------


def _sample_names(exp):
    names = []
    for df in (exp.peptides, exp.proteins):
        if df is None:
            continue
        for c in df.columns:
            if c.startswith("Area ") and c[5:] not in names:
                names.append(c[5:])
        if names:
            break
    return names or [exp.sample_name or "Sample 1"]


def identification_counts(exp):
    """Identified #PSMs/#Scans/#Features, #Peptides, #Sequences and #Proteins from the CSVs
    (definitions: docs/PEAKS_EXPORT_FORMAT.md)."""
    c = {}
    ps = exp.psms
    if ps is not None and len(ps):
        src = [_norm_source(x) for x in ps["Source File"].tolist()] if "Source File" in ps else [""] * len(ps)
        scans = ps["Scan"].tolist()
        c["identified_scans"] = len(set(zip(src, scans)))
        if "Feature Id" in ps:
            fid = ps["Feature Id"].fillna(0).astype(int).tolist()
            c["psms"] = len(set(zip(src, scans, fid)))
            c["identified_features"] = len({(s, f) for s, f in zip(src, fid) if f > 0})
    elif ps is not None:
        c.update(identified_scans=0, psms=0, identified_features=0)
    peps = exp.peptides["Peptide"].tolist() if exp.peptides is not None else (
        sorted(set(ps["Peptide"])) if ps is not None else [])
    c["peptides"] = len({seqnorm.il_collapse(p) for p in peps})
    seqs = set()
    for p in peps:
        try:
            seqs.add(seqnorm.il_collapse(seqnorm.strip_mods(p)))
        except ValueError:
            seqs.add(p)
    c["sequences"] = len(seqs)
    pr = exp.proteins
    if pr is not None:
        c["proteins_all"] = len(pr)
        c["protein_groups"] = int(pr["Protein Group"].nunique()) if "Protein Group" in pr else None
        c["proteins_top"] = int(pr["Top"].map(_truthy).sum()) if "Top" in pr else None
    return c


def missed_cleavage_table(exp, enzyme, samples):
    rows = {}
    pep = exp.peptides
    if pep is None:
        return rows
    for s in samples:
        df = pep
        col = f"#Feature {s}"
        if len(samples) > 1 and col in pep:
            df = pep[pep[col].fillna(0).astype(float) > 0]
        c = Counter()
        for p in df["Peptide"].tolist():
            try:
                st = seqnorm.strip_mods(p)
            except ValueError:
                continue
            n = 0
            if enzyme.cleave_after:
                for j in range(1, len(st)):
                    if st[j - 1] in enzyme.cleave_after and st[j] not in enzyme.not_before:
                        n += 1
            c["4+" if n >= 4 else str(n)] += 1
        rows[s] = {k: int(c.get(k, 0)) for k in ("0", "1", "2", "3", "4+")}
    return rows


def compute_statistics(exp, inputs, mzml_readers, notes):
    """Table 1 values (Total row) and the sources they came from."""
    st = {k: None for k in ("ms1_scans", "ms2_scans", "chimera_scans", "features")}
    src = {}
    if mzml_readers:
        st["ms1_scans"] = sum(r.levels.get(1, 0) for r in mzml_readers.values())
        st["ms2_scans"] = sum(r.levels.get(2, 0) for r in mzml_readers.values())
        src["ms1_scans"] = src["ms2_scans"] = "mzML"
    if inputs.get("precursors"):
        n_sc, n_ch = precursor_scan_counts(inputs["precursors"])
        st["chimera_scans"] = n_ch
        src["chimera_scans"] = "precursors.parquet"
        if st["ms2_scans"] is None:
            st["ms2_scans"], src["ms2_scans"] = n_sc, "precursors.parquet"
    elif inputs.get("mgf"):
        n_sc = n_ch = 0
        for m in inputs["mgf"]:
            _log(f"counting MGF entries in {m}")
            a, b = mgf_scan_counts(m)
            n_sc += a
            n_ch += b
        st["chimera_scans"], src["chimera_scans"] = n_ch, "MGF"
        if st["ms2_scans"] is None:
            st["ms2_scans"], src["ms2_scans"] = n_sc, "MGF"
    if inputs.get("features"):
        st["features"] = int(sum(len(pd.read_parquet(f, columns=["feature_id"])) for f in inputs["features"]))
        src["features"] = "ms1_features.parquet"
    st.update(identification_counts(exp))
    missing = [k for k in ("ms1_scans", "ms2_scans", "chimera_scans", "features") if st[k] is None]
    if missing:
        label = {"ms1_scans": "#MS1 (needs --mzml)", "ms2_scans": "#MS/MS (needs --mzml or a precursor list)",
                 "chimera_scans": "#Chimera (needs precursors.parquet or the PEAKS MGF)",
                 "features": "#Features (needs ms1_features.parquet)"}
        notes.append("Table 1: " + ", ".join(label[k] for k in missing) + " not available; shown as -.")
    return st, src


# --------------------------------------------------------------------------
# figures (matplotlib, Agg)
# --------------------------------------------------------------------------


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": "#b9b8b2",
        "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
        "xtick.major.size": 3, "ytick.major.size": 0, "legend.frameon": False,
        "svg.fonttype": "none",
    })
    return plt


def _save(fig, path, scale=2.0, colors=96):
    """Save a PNG (palette-quantised to keep the img folder small)."""
    plt = _mpl()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100 * scale, facecolor="white")
    plt.close(fig)
    try:
        from PIL import Image
        buf.seek(0)
        im = Image.open(buf).convert("RGB")
        im = im.quantize(colors=colors, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE)
        im.save(path, format="PNG", optimize=True)
    except Exception:  # PIL missing or failed: keep matplotlib's PNG
        with open(path, "wb") as fh:
            fh.write(buf.getvalue())


def _ygrid(ax):
    ax.grid(axis="y", color=GRID, linewidth=0.7, linestyle=(0, (4, 3)))
    ax.set_axisbelow(True)


def _legend_below(ax, ncol):
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=ncol, fontsize=7.5, handlelength=1,
              handleheight=1, columnspacing=1.2)


def _thousands(ax, axis="y"):
    from matplotlib.ticker import FuncFormatter
    f = FuncFormatter(lambda v, _p: f"{v:,.0f}")
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_formatter(f)


def fig_fdr_curve(path, sc, alpha):
    """(left) FDR vs number of target peptide sequences kept, (right) the same for PSMs;
    dashed line = the current -10lgP threshold."""
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 2.5))
    ymax = max(0.06, 6 * alpha)
    thr = sc["threshold"]
    for ax, curve, label, color in ((axes[0], sc["curve_seq"], "Number of peptide sequences", C_BLUE),
                                    (axes[1], sc["curve_psm"], "Number of PSMs", C_BLUE)):
        c = curve[curve["fdr"] <= ymax * 1.05]
        if len(c):
            ax.plot(c["n_target"], 100 * c["fdr"], color=color, lw=2)
        if thr is not None and len(curve):
            at = curve[curve["threshold"] >= thr]
            if len(at):
                n_at, f_at = int(at["n_target"].iloc[-1]), float(at["fdr"].iloc[-1])
                ax.axvline(n_at, color=MUTED, lw=1.2, ls=(0, (5, 3)))
                kind = "Sequences" if "sequence" in label else "PSMs"
                ax.text(0.03, 0.95, f"{kind} : {n_at:,}\nFDR : {_pct(f_at)}", transform=ax.transAxes,
                        va="top", ha="left", color=INK, fontsize=8)
        ax.set_ylim(0, 100 * ymax)
        ax.set_xlim(0, None)
        ax.set_xlabel(label)
        ax.set_ylabel("FDR (%)")
        _ygrid(ax)
        _thousands(ax, "x")
    if thr is not None:
        fig.suptitle(f"current threshold: -10lgP \u2265 {thr:.2f}", fontsize=8, color=INK2, x=0.99, ha="right", y=0.99)
    fig.tight_layout()
    _save(fig, path)


def fig_score_hist(path, scores, decoy, thr):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.0, 3.2))
    scores = np.asarray(scores, float)
    hi = max(5.0, math.ceil(np.nanmax(scores)) + 1) if len(scores) else 50.0
    bins = np.arange(0, hi + 1, 1.0)
    if decoy is not None and np.any(decoy):
        ht, _ = np.histogram(scores[~decoy], bins)
        hd, _ = np.histogram(scores[decoy], bins)
        ax.bar(bins[:-1], hd, width=0.9, align="edge", color=C_ORANGE, label="Decoy", linewidth=0)
        ax.bar(bins[:-1], ht, width=0.9, align="edge", bottom=hd, color=C_BLUE, label="Target", linewidth=0)
        _legend_below(ax, 2)
    else:
        ht, _ = np.histogram(scores, bins)
        ax.bar(bins[:-1], ht, width=0.9, align="edge", color=C_BLUE, linewidth=0)
    if thr is not None:
        ax.axvline(thr, color=MUTED, lw=1.2, ls=(0, (5, 3)))
    ax.set_xlabel("PSM score (-10lgP)")
    ax.set_ylabel("Number of PSMs")
    ax.set_xlim(0, hi)
    _ygrid(ax)
    _thousands(ax)
    fig.tight_layout()
    _save(fig, path)


def fig_score_ppm(path, scores, ppm, decoy, thr, max_points=60000):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.0, 3.2 if decoy is not None else 2.9))
    s, p = np.asarray(scores, float), np.asarray(ppm, float)
    d = np.zeros(len(s), bool) if decoy is None else np.asarray(decoy, bool)
    rng = np.random.default_rng(0)
    for mask, color, label in ((~d, C_BLUE, "Target"), (d, C_ORANGE, "Decoy")):
        idx = np.nonzero(mask & np.isfinite(s) & np.isfinite(p))[0]
        if len(idx) == 0:
            continue
        if len(idx) > max_points:
            idx = rng.choice(idx, max_points, replace=False)
        ax.scatter(s[idx], p[idx], s=4, marker="s", color=color, alpha=0.25, linewidths=0,
                   rasterized=True, label=label)
    if thr is not None:
        ax.axvline(thr, color=MUTED, lw=1.2, ls=(0, (5, 3)))
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("PSM score (-10lgP)")
    ax.set_ylabel("ppm")
    if np.any(d):
        leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.5, markerscale=2.5)
        for h in leg.legend_handles:
            h.set_alpha(1)
    _ygrid(ax)
    fig.tight_layout()
    _save(fig, path)


def fig_overlaid_hist(path, series, bins, xlabel, ylabel="Number of features", log_x=False):
    """Overlaid histograms of nested sets (largest first, drawn behind)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.0, 3.2))
    for label, vals, color in series:
        v = np.asarray(vals, float)
        v = v[np.isfinite(v)]
        h, _ = np.histogram(v, bins)
        w = np.diff(bins)
        ax.bar(bins[:-1], h, width=w * 0.9, align="edge", color=color, label=f"{label} ({len(v):,})", linewidth=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if log_x:
        lo, hi = int(math.floor(bins[0])), int(math.ceil(bins[-1]))
        ax.set_xticks(range(lo, hi + 1))
        ax.set_xticklabels([f"1E{k}" for k in range(lo, hi + 1)])
    _ygrid(ax)
    _thousands(ax)
    _legend_below(ax, len(series))
    fig.tight_layout()
    _save(fig, path)


def fig_ppm_hist(path, series, lim):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.0, 3.2 if len(series) > 1 else 2.9))
    bins = np.arange(-lim, lim + 0.5, 0.5)
    for label, vals, color in series:
        v = np.asarray(vals, float)
        v = v[np.isfinite(v)]
        h, _ = np.histogram(np.clip(v, -lim, lim - 1e-9), bins)
        ax.bar(bins[:-1], h, width=0.45, align="edge", color=color, alpha=0.85 if len(series) > 1 else 1,
               label=label, linewidth=0)
    ax.set_xlabel("mass error (ppm)")
    ax.set_ylabel("Number of PSMs")
    _ygrid(ax)
    _thousands(ax)
    if len(series) > 1:
        _legend_below(ax, 2)
    fig.tight_layout()
    _save(fig, path)


def fig_ppm_mz(path, mz, ppm, lim):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.0, 2.9))
    ax.scatter(mz, ppm, s=5, marker="s", color=C_BLUE, alpha=0.2, linewidths=0, rasterized=True)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("m/z")
    ax.set_ylabel("mass error (ppm)")
    _ygrid(ax)
    _thousands(ax, "x")
    fig.tight_layout()
    _save(fig, path)


def fig_venn(path, sets, names):
    """Venn diagram for 2 or 3 sets (region counts; circles not area-proportional)."""
    plt = _mpl()
    from matplotlib.patches import Circle
    fig, ax = plt.subplots(figsize=(2.7, 2.5))
    ax.set_aspect("equal")
    ax.axis("off")
    k = len(sets)
    centers = [(-0.5, 0.0), (0.5, 0.0)] if k == 2 else [(-0.5, 0.3), (0.5, 0.3), (0.0, -0.55)]
    for (x, y), c in zip(centers, (C_BLUE, C_ORANGE, C_AQUA)):
        ax.add_patch(Circle((x, y), 0.95, facecolor=c, alpha=0.25, edgecolor=c, lw=1.5))
    allx = set().union(*sets)
    regions = Counter()
    for x in allx:
        regions[tuple(x in s for s in sets)] += 1
    if k == 2:
        pos = {(True, False): (-0.85, 0), (False, True): (0.85, 0), (True, True): (0, 0)}
        lab = [(-0.6, 1.1), (0.6, 1.1)]
    else:
        pos = {(True, False, False): (-0.85, 0.55), (False, True, False): (0.85, 0.55),
               (False, False, True): (0, -1.0), (True, True, False): (0, 0.7),
               (True, False, True): (-0.55, -0.4), (False, True, True): (0.55, -0.4), (True, True, True): (0, 0.05)}
        lab = [(-1.0, 1.35), (1.0, 1.35), (0, -1.65)]
    for key, (x, y) in pos.items():
        ax.text(x, y, f"{regions.get(key, 0):,}", ha="center", va="center", fontsize=8, color=INK)
    for (x, y), n in zip(lab, names):
        ax.text(x, y, n, ha="center", va="center", fontsize=8, color=INK2)
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.8 if k == 3 else -1.1, 1.5)
    fig.tight_layout()
    _save(fig, path)


# --------------------------------------------------------------------------
# annotated spectra
# --------------------------------------------------------------------------


def _residue_masses(stripped, mods):
    r = [MONO[a] for a in stripped]
    n = len(r)
    for pos, d in mods:
        if pos <= 1:
            r[0] += d
        elif pos > n:
            r[-1] += d
        else:
            r[pos - 1] += d
    return r


def annotate_fragments(stripped, mods, z, mz, inten, tol=0.5, unit="Da"):
    """Match b/y ions (1+, and 2+ when the precursor charge is >= 3) to peaks.

    ``mods``: [(position, exact delta)] with 0 = N-term and n+1 = C-term.  Each theoretical
    ion takes the most intense peak within ``tol``; a peak keeps only its first assignment
    (1+ before 2+).  Returns a list of dicts (ion, kind, index, charge, theo, obs, inten, err).
    """
    mz = np.asarray(mz, float)
    inten = np.asarray(inten, float)
    order = np.argsort(mz)
    mz, inten = mz[order], inten[order]
    r = _residue_masses(stripped, mods)
    n = len(r)
    cum = np.cumsum(r)
    total = cum[-1]
    charges = [1, 2] if (z or 2) >= 3 else [1]
    out, used = [], set()
    for c in charges:
        for i in range(1, n):
            for kind, neutral in (("b", cum[i - 1]), ("y", total - cum[n - i - 1] + H2O)):
                theo = (neutral + c * PROTON) / c
                w = tol if unit == "Da" else theo * tol * 1e-6
                lo, hi = np.searchsorted(mz, theo - w), np.searchsorted(mz, theo + w, side="right")
                if hi <= lo:
                    continue
                j = lo + int(np.argmax(inten[lo:hi]))
                if j in used:
                    continue
                used.add(j)
                err = mz[j] - theo if unit == "Da" else (mz[j] - theo) / theo * 1e6
                out.append({"ion": f"{kind}{i}" + (f"[{c}+]" if c > 1 else ""), "kind": kind, "index": i,
                            "charge": c, "theo": theo, "obs": float(mz[j]), "inten": float(inten[j]), "err": float(err)})
    return out


def render_spectrum(job):
    """Annotated MS2 spectrum PNG: sequence ladder, peaks with b/y labels, fragment errors."""
    try:
        plt = _mpl()
        from matplotlib.collections import LineCollection
        stripped, mods, z = job["stripped"], job["mods"], job["z"]
        mz, it = np.asarray(job["mz"], float), np.asarray(job["inten"], float)
        tol, unit = job["tol"], job["unit"]
        rel = 100.0 * it / it.max() if len(it) and it.max() > 0 else it
        ann = annotate_fragments(stripped, mods, z, mz, rel, tol, unit)
        fig = plt.figure(figsize=(12.0, 3.3))
        gs = fig.add_gridspec(3, 1, height_ratios=[0.55, 3.0, 0.9], hspace=0.08,
                              left=0.055, right=0.985, top=0.9, bottom=0.12)
        axl, ax, axe = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
        fig.text(0.055, 0.955, job["title"], ha="left", va="center", fontsize=9, color=INK)
        # ladder
        n = len(stripped)
        axl.axis("off")
        axl.set_xlim(-0.8, max(n, 25) - 0.2)
        axl.set_ylim(0, 1)
        bset = {a["index"] for a in ann if a["kind"] == "b"}
        yset = {a["index"] for a in ann if a["kind"] == "y"}
        for i, aa in enumerate(stripped):
            modded = any((p == i + 1) or (p == 0 and i == 0) or (p == n + 1 and i == n - 1) for p, _d in mods)
            axl.text(i, 0.5, aa, ha="center", va="center", fontsize=10, family="DejaVu Sans Mono",
                     weight="bold", color=C_ORANGE if modded else INK)
        for i in range(1, n):
            x = i - 0.5
            if i in bset:
                axl.plot([x, x, x - 0.35], [0.5, 0.05, 0.05], color=B_COLOR, lw=1.3, solid_capstyle="butt")
            if (n - i) in yset:
                axl.plot([x, x, x + 0.35], [0.5, 0.95, 0.95], color=Y_COLOR, lw=1.3, solid_capstyle="butt")
        # peaks
        xmax = max(float(mz.max()) if len(mz) else 1000.0, job.get("pmz") or 0) * 1.03
        segs = [((m, 0), (m, h)) for m, h in zip(mz, rel)]
        ax.add_collection(LineCollection(segs, colors="#a3a29c", linewidths=0.6))
        for kind, color in (("b", B_COLOR), ("y", Y_COLOR)):
            a = [x for x in ann if x["kind"] == kind]
            if a:
                ax.add_collection(LineCollection([((x["obs"], 0), (x["obs"], x["inten"])) for x in a],
                                                 colors=color, linewidths=1.2))
        # labels: most intense matched first, skip if crowded
        taken = []
        for a in sorted(ann, key=lambda x: -x["inten"]):
            if a["inten"] < 1.5 or len(taken) >= 45:
                continue
            x = a["obs"]
            if any(abs(x - t) < xmax * 0.012 and abs(a["inten"] - h) < 12 for t, h in taken):
                continue
            taken.append((x, a["inten"]))
            ax.text(x, a["inten"] + 2, a["ion"], ha="center", va="bottom", fontsize=7.5,
                    color=B_COLOR if a["kind"] == "b" else Y_COLOR)
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, 118)
        ax.set_yticks([0, 50, 100])
        ax.set_ylabel("Intensity (%)")
        ax.tick_params(axis="x", labelbottom=False)
        _ygrid(ax)
        # fragment errors
        for kind, color in (("b", B_COLOR), ("y", Y_COLOR)):
            a = [x for x in ann if x["kind"] == kind]
            if a:
                axe.scatter([x["theo"] for x in a], [x["err"] for x in a], s=9, marker="s", color=color, linewidths=0)
        axe.axhline(0, color=INK2, lw=0.7)
        axe.set_xlim(0, xmax)
        axe.set_ylim(-tol * 1.15, tol * 1.15)
        axe.set_yticks([-tol, 0, tol])
        axe.set_yticklabels([f"{-tol:g}", "0", f"{tol:g}"])
        axe.set_ylabel(unit, rotation=0, ha="right", va="center")
        axe.set_xlabel("m/z", labelpad=1)
        nb, ny = len(bset), len(yset)
        fig.text(0.985, 0.955, f"b {nb}/{n - 1}   y {ny}/{n - 1}   tol {tol:g} {unit}", ha="right",
                 va="center", fontsize=8, color=INK2)
        _save(fig, job["path"], scale=1.0, colors=48)
        return job["path"], None
    except Exception as e:  # never fail the page for one image
        return job["path"], f"{type(e).__name__}: {e}"


# --------------------------------------------------------------------------
# coverage maps
# --------------------------------------------------------------------------

_FS = 9.0                       # sequence font size (pt)
_CW = 0.6021 * _FS * 100 / 72   # DejaVu Sans Mono advance (px at 100 dpi)
_PER_ROW = 100
_X0 = 58
_W = 1150


def _xpos(i):
    """x (px) of residue i (0-based, within the row) incl. a space every 10 residues."""
    return _X0 + (i + i // 10) * _CW


def render_coverage(job):
    """Sequence coverage map PNG: covered residues shaded, peptides as bars (blue unique,
    grey shared), variable modifications as coloured squares, legend on the right."""
    try:
        plt = _mpl()
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Rectangle
        seq = job["seq"]
        L = len(seq)
        peps = job["peps"]          # [(start0, end0_excl, unique, [(abs0, modname)])]
        colors = job["mod_colors"]  # {name: color}
        cov = np.zeros(L, bool)
        modpos = {}
        for s, e, _u, mods in peps:
            cov[max(0, s):min(L, e)] = True
            for p, name in mods:
                if 0 <= p < L:
                    modpos.setdefault(p, name)
        nrows = (L + _PER_ROW - 1) // _PER_ROW
        row_cov = [bool(cov[r * _PER_ROW:(r + 1) * _PER_ROW].any()) for r in range(nrows)]
        # collapse long uncovered stretches of long proteins
        items, r = [], 0
        while r < nrows:
            if L > 1500 and not row_cov[r]:
                k = r
                while k < nrows and not row_cov[k]:
                    k += 1
                if k - r >= 3:
                    items.append(("gap", r, k))
                    r = k
                    continue
            items.append(("row", r, r + 1))
            r += 1
        # lanes per row
        layout, y = [], 30.0
        for kind, r0, r1 in items:
            if kind == "gap":
                layout.append((kind, r0, r1, y, []))
                y += 20
                continue
            a, b = r0 * _PER_ROW, min(L, r1 * _PER_ROW)
            ivs = sorted((max(s, a), min(e, b), u, mods) for s, e, u, mods in peps if s < b and e > a)
            lanes_end, lanes = [], []
            for s, e, u, mods in ivs:
                for li, end in enumerate(lanes_end):
                    if s >= end:
                        lanes_end[li] = e
                        lanes.append((li, s, e, u, mods))
                        break
                else:
                    lanes_end.append(e)
                    lanes.append((len(lanes_end) - 1, s, e, u, mods))
            layout.append((kind, r0, r1, y, lanes))
            y += 26 + 7 * len(lanes_end) + 8
        H = int(y + 6)
        fig = plt.figure(figsize=(_W / 100, H / 100))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, _W)
        ax.set_ylim(H, 0)
        ax.axis("off")
        ax.text(6, 14, job["header"][:150], ha="left", va="center", fontsize=8.5, color=INK)
        shade, bars_u, bars_s, marks = [], [], [], defaultdict(list)
        for kind, r0, r1, y0, lanes in layout:
            if kind == "gap":
                ax.text(_X0, y0 + 8, f"… residues {r0 * _PER_ROW + 1}–{min(L, r1 * _PER_ROW)}: "
                        "no identified peptides …", ha="left", va="center", fontsize=8,
                        style="italic", color=MUTED)
                continue
            a, b = r0 * _PER_ROW, min(L, r1 * _PER_ROW)
            base = y0 + 14
            ax.text(_X0 - 10, base, str(a + 1), ha="right", va="baseline", fontsize=8.5, color=MUTED)
            # covered shading (contiguous runs)
            i = a
            while i < b:
                if cov[i]:
                    j = i
                    while j < b and cov[j]:
                        j += 1
                    x1 = _xpos(i - a) - 1
                    x2 = _xpos(j - 1 - a) + _CW + 1
                    shade.append(Rectangle((x1, base - 11), x2 - x1, 15))
                    i = j
                else:
                    i += 1
            # text runs (bold when covered); spaces inserted at block boundaries
            i = a
            while i < b:
                c = cov[i]
                j = i
                while j < b and cov[j] == c:
                    j += 1
                txt = "".join((" " if (k - a) % 10 == 0 and k != i else "") + seq[k] for k in range(i, j))
                ax.text(_xpos(i - a), base, txt, ha="left", va="baseline", fontsize=_FS,
                        family="DejaVu Sans Mono", weight="bold" if c else "normal", color=INK if c else INK2)
                i = j
            for p in range(a, b):
                if p in modpos:
                    marks[modpos[p]].append(Rectangle((_xpos(p - a) + 0.5, base - 19), _CW - 1, 6))
            for li, s, e, u, mods in lanes:
                yb = y0 + 22 + 7 * li
                x1, x2 = _xpos(s - a), _xpos(e - 1 - a) + _CW
                (bars_u if u else bars_s).append(Rectangle((x1, yb), x2 - x1, 4))
                for p, name in mods:
                    if s <= p < e:
                        marks[name].append(Rectangle((_xpos(p - a) + 0.5, yb - 1), _CW - 1, 6))
        ax.add_collection(PatchCollection(shade, facecolor="#dcdcd8", edgecolor="none", zorder=0))
        ax.add_collection(PatchCollection(bars_u, facecolor=C_BLUE, edgecolor="none"))
        ax.add_collection(PatchCollection(bars_s, facecolor="#a3a29c", edgecolor="none"))
        for name, rects in marks.items():
            ax.add_collection(PatchCollection(rects, facecolor=colors.get(name, C_VIOLET), edgecolor="white",
                                              linewidth=0.5, zorder=3))
        # legend
        lx, ly = _X0 + 112 * _CW, 34
        for label, color in [("unique peptide", C_BLUE), ("shared peptide", "#a3a29c")] + \
                [(n, colors.get(n, C_VIOLET)) for n in sorted({m for v in modpos.values() for m in [v]})]:
            ax.add_patch(Rectangle((lx, ly - 4), 10, 8, facecolor=color, edgecolor="none"))
            ax.text(lx + 15, ly, label, ha="left", va="center", fontsize=7.5, color=INK2)
            ly += 15
        _save(fig, job["path"], scale=1.0, colors=64)
        return job["path"], None
    except Exception as e:
        return job["path"], f"{type(e).__name__}: {e}"


def _run_jobs(fn, jobs, processes):
    if not jobs:
        return []
    if processes <= 1 or len(jobs) < 8:
        return [fn(j) for j in jobs]
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures.process import BrokenProcessPool
    # fork where available: the workers only run matplotlib on the job data, and fork does not
    # re-import __main__ (spawn fails when the caller runs from stdin or a notebook)
    method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    chunk = max(1, min(32, len(jobs) // (processes * 4) or 1))
    try:
        with ProcessPoolExecutor(max_workers=processes, mp_context=mp.get_context(method)) as ex:
            return list(ex.map(fn, jobs, chunksize=chunk))
    except BrokenProcessPool as e:
        _log(f"worker pool failed ({e}); drawing in this process")
        return [fn(j) for j in jobs]


# --------------------------------------------------------------------------
# HTML
# --------------------------------------------------------------------------

CSS = """/* proteoscan export pages */
:root { --ink:#0b0b0b; --ink2:#52514e; --muted:#8a8984; --line:#cfcec8; --head:#eef3fb;
        --bar:#dfe9f8; --accent:#1f5fae; --row:#f6f6f4; --bg:#ffffff; }
body { margin:0; background:var(--bg); color:var(--ink);
       font:13px/1.45 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }
.container { max-width:860px; margin:0 auto; padding:16px; }
.wide { padding:12px 16px; }
h1 { font-size:20px; margin:28px 0 6px; padding-bottom:4px; border-bottom:2px solid var(--line); }
.section { padding:6px 0 22px; }
.caption { margin:10px 0 8px; color:var(--ink); }
.caption b { margin-right:4px; }
.figrow { display:flex; flex-wrap:wrap; gap:8px; }
.figrow > div { flex:1 1 260px; text-align:center; }
.sub { font-weight:bold; text-align:center; margin-bottom:2px; }
img { max-width:100%; height:auto; }
.na { padding:18px 8px; color:var(--ink2); font-weight:bold; }
.note { color:var(--ink2); }
ul.notes { margin:6px 0; padding-left:20px; }
table.stats, table.grid { border-collapse:collapse; margin:6px 0; }
table.stats td, table.stats th { padding:3px 14px 3px 4px; text-align:left; white-space:nowrap;
       border-bottom:1px solid #ecebe7; }
table.grid td, table.grid th { border:1px solid var(--line); padding:3px 7px; text-align:center; }
table.grid th { background:var(--head); font-weight:600; }
.foot { font-size:11px; color:var(--ink2); }
.align-left { text-align:left; line-height:1.6; }
.bar { display:flex; justify-content:space-between; align-items:baseline; margin-top:22px;
       padding:5px 10px; background:var(--bar); color:#0d3566; font-weight:bold; font-size:15px; }
.bar a { color:#0d3566; font-weight:normal; font-size:12px; text-decoration:none; }
.bar .tit a { font-weight:bold; font-size:15px; }
.lst { font-weight:bold; color:#0d3566; margin:10px 0 4px; }
.lst a { color:#0d3566; }
.meta { color:var(--ink2); margin:4px 0; }
table.mytable { border-collapse:collapse; font-size:11.5px; margin:4px 0; }
table.mytable th { background:var(--head); position:sticky; top:0; }
table.mytable td, table.mytable th { border:1px solid var(--line); padding:2px 5px; white-space:nowrap; }
table.mytable td.wt, table.mytable td.cl:last-child { white-space:normal; }
table.mytable tbody tr:nth-child(even) { background:var(--row); }
table.mytable tbody tr:hover { background:#e8eef8; }
.cl { text-align:left; } .cc { text-align:center; } .cr { text-align:right; }
.wt { text-align:left; word-break:break-all; max-width:340px; }
.mono { font-family:Menlo, Consolas, "DejaVu Sans Mono", monospace; }
.dim { color:var(--muted); }
a { color:var(--accent); text-decoration:none; }
a:hover { text-decoration:underline; }
.generator { margin-top:30px; font-size:11px; color:var(--muted); border-top:1px solid var(--line); padding-top:6px; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) { --ink:#f0efec; --ink2:#c3c2b7; --muted:#9b9a92; --line:#45443f;
     --head:#24303f; --bar:#22324a; --row:#1f1f1e; --bg:#1a1a19; --accent:#7fb0ef; }
  :root:not([data-theme="light"]) .bar, :root:not([data-theme="light"]) .bar a,
  :root:not([data-theme="light"]) .lst, :root:not([data-theme="light"]) .lst a { color:#cfe0f7; }
  :root:not([data-theme="light"]) img { background:#fff; }
}
:root[data-theme="dark"] { --ink:#f0efec; --ink2:#c3c2b7; --muted:#9b9a92; --line:#45443f;
     --head:#24303f; --bar:#22324a; --row:#1f1f1e; --bg:#1a1a19; --accent:#7fb0ef; }
@media print { .bar a { display:none; } table.mytable th { position:static; } }
"""


def _page(title, body, wide=False):
    return ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            f"<meta name=\"generator\" content=\"proteoscan {__version__}\">\n"
            f"<title>{_esc(title)}</title>\n<link rel=\"stylesheet\" href=\"export.css\">\n</head>\n<body>\n"
            f"<div class=\"{'wide' if wide else 'container'}\">\n{body}\n"
            f"<div class=\"generator\">Prepared with proteoscan {__version__} (open source, MIT); "
            f"layout follows the PEAKS Studio HTML export. {_dt.datetime.now():%Y-%m-%d %H:%M}</div>\n"
            "</div>\n</body>\n</html>\n")


def _img(src, alt, width=None, lazy=False):
    w = f' width="{int(width)}"' if width else ""
    lz = ' loading="lazy"' if lazy else ""
    return f'<img src="{_esc(src)}" alt="{_esc(alt)}"{w}{lz}>'


def _v(x):
    return "-" if x is None else (f"{x}" if isinstance(x, (int, np.integer)) else _esc(x))


def html_table1(stats, samples):
    th = 'style="text-align:center"'
    h1 = ("<tr><th rowspan=\"2\"></th><th colspan=\"3\">#Scans</th><th rowspan=\"2\">#Features</th>"
          "<th colspan=\"3\">Identified</th><th rowspan=\"2\">#Peptides</th><th rowspan=\"2\">#Sequences</th>"
          "<th colspan=\"3\">#Proteins<sup>*</sup></th></tr>")
    h2 = ("<tr><th>MS1</th><th>MS/MS</th><th>#Chimera</th><th>#PSMs</th><th>#Scans</th>"
          "<th>#Features<sup>**</sup></th><th>Groups</th><th>All</th><th>Top</th></tr>")
    keys = ["ms1_scans", "ms2_scans", "chimera_scans", "features", "psms", "identified_scans",
            "identified_features", "peptides", "sequences", "protein_groups", "proteins_all", "proteins_top"]
    rows = [("Total", stats)]
    if len(samples) == 1:
        rows.append((samples[0], stats))
    body = "".join("<tr><td>" + _esc(n) + "</td>" + "".join(f"<td>{_v(s.get(k))}</td>" for k in keys) + "</tr>"
                   for n, s in rows)
    del th
    return (f"<table class=\"grid\">{h1}{h2}{body}</table>"
            "<div class=\"foot\"><sup>* proteins with significant peptides are used in counts.</sup><br>"
            "<sup>** features are identified by DB search only.</sup></div>")


def _stats_table(rows):
    return "<table class=\"stats\">" + "".join(
        f"<tr><td>{_esc(a)}</td><td>{b}</td></tr>" for a, b in rows) + "</table>"


def _fig_section(num, caption, inner, note=None):
    n = f'<div class="note">{note}</div>' if note else ""
    return (f'<div class="section"><div class="caption"><b>Figure {num}.</b> {caption}</div>{n}{inner}</div>')


# --------------------------------------------------------------------------
# main builder
# --------------------------------------------------------------------------


def _mod_name_map(exp, filters):
    """{delta string: PTM name} from the PTM columns of the export."""
    names = {}
    for df in (exp.psms, exp.peptides):
        if df is None or "PTM" not in df or "Peptide" not in df:
            continue
        for pep, ptm in zip(df["Peptide"].tolist()[:30000], df["PTM"].tolist()[:30000]):
            if not isinstance(ptm, str) or not ptm:
                continue
            nm = [x.strip() for x in re.split(r";", ptm) if x.strip()]
            try:
                _s, mods = seqnorm.parse(pep)
            except (ValueError, TypeError):
                continue
            if len(mods) == len(nm):
                for (_p, d), n in zip(mods, nm):
                    names.setdefault(seqnorm.format_delta(d), n)
    return names


def _exact(d):
    key = seqnorm.format_delta(d)
    for v in seqnorm.MOD_MASSES.values():
        if seqnorm.format_delta(v) == key:
            return float(v)
    return float(d)


_FIGURE_NAMES = {"fdr_curve", "score_histogram", "score_vs_ppm", "feature_mz", "feature_rt", "feature_abundance",
                 "precursor_error_histogram", "precursor_error_vs_mz", "venn_all_proteins", "venn_top_proteins",
                 "venn_peptides"}


def build_summary(export_dir, out_dir, mzml=(), fasta=None, psms_scored=None, features=None,
                  precursors=None, config=None, mgf=None, run_dir=None, max_spectra=2000,
                  processes=8, title=None, images=True):
    """Write db.summary.html, db.proteins.html, db.peptides.html, export.css and img/ into
    ``out_dir``.  ``config``: a :class:`SearchConfig`, a YAML path, or None.  Returns a dict
    with sizes, counts and the computed Table 1/3/5 values."""
    t0 = time.time()
    export_dir = os.path.abspath(export_dir)
    out_dir = os.path.abspath(out_dir)
    img_dir = os.path.join(out_dir, "img")
    os.makedirs(img_dir, exist_ok=True)
    for old in glob.glob(os.path.join(glob.escape(img_dir), "*.png")):  # our own figures from an earlier run
        b = os.path.basename(old)
        if re.match(r"(spec|cov)_\d{5}\.png$", b) or b[:-4] in _FIGURE_NAMES:
            os.remove(old)
    if isinstance(mzml, str):
        mzml = [mzml]
    mzml = [os.path.abspath(p) for p in (mzml or [])]
    inputs = resolve_inputs(export_dir, run_dir, psms_scored, features, precursors, config, mgf)
    notes = []
    _log(f"reading export {export_dir}")
    exp = read_export(export_dir, tables=("psms", "peptides", "proteins", "protein_peptides", "features",
                                          "ptm_profile"), verify=False)
    samples = _sample_names(exp)

    # ---- parameters
    cfg = inputs["config"]
    if isinstance(cfg, str):
        from .config import SearchConfig
        cfg = SearchConfig.from_yaml(cfg)
    if cfg is not None:
        lines6, filters, enz_params = params_from_config(cfg, fasta)
        param_src = "the pipeline configuration"
    elif exp.params:
        lines6, filters, enz_params = params_from_export(exp.params)
        param_src = "the export's own db.summary.html"
    else:
        lines6, filters, enz_params = [], {"peptide_fdr": 0.01, "protein_min_10lgp": 15,
                                           "protein_min_unique_peptides": 1, "denovo_min_alc_pct": 50,
                                           "fragment_tolerance": 0.5, "fragment_tolerance_unit": "Da",
                                           "fixed_deltas": set()}, {}
        param_src = None
        notes.append("No search configuration (--config) and no db.summary.html in the export: Table 6 is "
                     "empty and default filters (protein -10lgP >= 15, 1 unique peptide) are assumed.")
    enzyme = Enzyme.from_params(enz_params)

    # ---- scored decoys
    sc = None
    if inputs["psms_scored"] and os.path.exists(inputs["psms_scored"]):
        _log(f"reading {inputs['psms_scored']}")
        sc = scored_statistics(inputs["psms_scored"], filters.get("peptide_fdr", 0.01),
                               float(filters.get("protein_min_10lgp") or 15),
                               int(filters.get("protein_min_unique_peptides") or 1))
    else:
        notes.append("No scored target/decoy PSM table (psms_scored.parquet): Figure 1 is omitted, Figure 2 "
                     "shows the exported (target) PSMs only, and Table 3 FDRs are not recomputed.")

    # ---- mzML
    readers = {}
    for p in mzml:
        _log(f"indexing {p}")
        readers[_norm_source(p)] = MzmlSpectra(p)
    if not mzml:
        notes.append("No mzML given (--mzml): #MS1 is not counted and no annotated spectra are drawn.")

    stats, stat_src = compute_statistics(exp, inputs, readers, notes)

    # ---- filters, Table 3
    if sc is not None and sc["threshold"] is not None:
        pep_thr = sc["threshold"]
    elif filters.get("peptide_min_10lgp") is not None:
        pep_thr = float(filters["peptide_min_10lgp"])
    elif exp.peptides is not None and len(exp.peptides):
        pep_thr = float(exp.peptides["-10LgP"].min())
        notes.append("Peptide -10lgP threshold taken as the lowest peptide score in db.peptides.csv.")
    else:
        pep_thr = None
    denovo_n = None
    dn = glob.glob(os.path.join(glob.escape(export_dir), "*.db.denovoOnly.csv"))
    if dn:
        with open(dn[0], "rb") as fh:
            denovo_n = max(0, sum(1 for _ in fh) - 1)
    reported = exp.params if (exp.params and cfg is None) else {}

    def _rep(key):
        v = reported.get(key)
        return f"{NA} (export reports {v}%)" if v is not None else NA
    if sc is not None:
        t3 = {"psm": _pct(sc["fdr_psm"]), "seq": _pct(sc["fdr_seq"]), "group": _pct(sc["fdr_group"])}
    else:
        t3 = {"psm": _rep("fdr_psm_pct"), "seq": _rep("fdr_peptide_pct"), "group": _rep("fdr_protein_group_pct")}
    t3["denovo"] = str(denovo_n) if denovo_n is not None else f"{NA} (no de novo results in the export)"
    mc = missed_cleavage_table(exp, enzyme, samples)

    # ---- figures
    figs = {}

    def fpath(name):
        figs[name] = f"img/{name}.png"
        return os.path.join(img_dir, f"{name}.png")
    alpha = float(filters.get("peptide_fdr") or 0.01)
    if sc is not None and len(sc["curve_seq"]):
        fig_fdr_curve(fpath("fdr_curve"), sc, alpha)
    if sc is not None:
        fig_score_hist(fpath("score_histogram"), sc["scores"], sc["decoy"], pep_thr)
        if sc["ppm"] is not None:
            fig_score_ppm(fpath("score_vs_ppm"), sc["scores"], sc["ppm"], sc["decoy"], pep_thr)
    elif exp.psms is not None and len(exp.psms):
        fig_score_hist(fpath("score_histogram"), exp.psms["-10LgP"].to_numpy(float), None, pep_thr)
        fig_score_ppm(fpath("score_vs_ppm"), exp.psms["-10LgP"].to_numpy(float), exp.psms["ppm"].to_numpy(float),
                      None, pep_thr)

    # Figure 4/5 feature distributions
    feat_series = {"mz": [], "rt": [], "area": []}
    allf = None
    if inputs["features"]:
        allf = pd.concat([pd.read_parquet(f, columns=["mz", "rt", "area"]) for f in inputs["features"]],
                         ignore_index=True)
        feat_series["mz"].append(("all", allf["mz"], C_BLUE))
        feat_series["rt"].append(("all", allf["rt"], C_BLUE))
        feat_series["area"].append(("all", np.log10(allf["area"].clip(lower=1e-9)), C_BLUE))
    ft = exp.features
    if ft is not None and len(ft):
        isdb = ft["DB Peptide"].astype(str).str.len() > 0 if "DB Peptide" in ft else pd.Series(True, index=ft.index)
        has_dn = "Denovo Peptide" in ft and (ft["Denovo Peptide"].astype(str).str.len() > 0).any()
        if has_dn:
            feat_series["mz"].append(("de novo+db", ft["m/z"], C_AQUA))
            feat_series["rt"].append(("de novo+db", ft["RT"], C_AQUA))
            feat_series["area"].append(("de novo+db", np.log10(ft["Area"].astype(float).clip(lower=1e-9)), C_AQUA))
        feat_series["mz"].append(("db", ft.loc[isdb, "m/z"], C_ORANGE))
        feat_series["rt"].append(("db", ft.loc[isdb, "RT"], C_ORANGE))
        feat_series["area"].append(("db", np.log10(ft.loc[isdb, "Area"].astype(float).clip(lower=1e-9)), C_ORANGE))
    if feat_series["mz"]:
        mzv = np.concatenate([np.asarray(v, float) for _l, v, _c in feat_series["mz"]])
        lo, hi = math.floor(np.nanmin(mzv) / 100) * 100, math.ceil(np.nanmax(mzv) / 100) * 100
        fig_overlaid_hist(fpath("feature_mz"), feat_series["mz"], np.arange(lo, hi + 20, 20), "Feature m/z")
        rtv = np.concatenate([np.asarray(v, float) for _l, v, _c in feat_series["rt"]])
        fig_overlaid_hist(fpath("feature_rt"), feat_series["rt"],
                          np.linspace(0, math.ceil(np.nanmax(rtv) / 10) * 10, 61), "Feature RT (min)")
        av = np.concatenate([np.asarray(v, float) for _l, v, _c in feat_series["area"]])
        av = av[np.isfinite(av) & (av > -8)]
        a_lo, a_hi = np.percentile(av, [0.05, 99.95])
        fig_overlaid_hist(fpath("feature_abundance"), feat_series["area"],
                          np.arange(math.floor(a_lo), math.ceil(a_hi) + 0.1, 0.1),
                          "Feature abundance (area)", log_x=True)
        if allf is None:
            notes.append("Figures 4 and 5 show the exported (identified) features only; the full feature list "
                         "(ms1_features.parquet) was not given.")
    else:
        notes.append("No feature table: Figures 4 and 5 are omitted.")

    # Figure 6 precursor error
    ppm_lim = None
    if exp.psms is not None and len(exp.psms):
        ps = exp.psms
        after = ps["ppm"].to_numpy(float)
        mass = ps["Mass"].to_numpy(float)
        before = ((ps["m/z"].to_numpy(float) - PROTON) * ps["z"].to_numpy(float) - mass) / mass * 1e6
        tolp = float(filters.get("precursor_tolerance") or 15.0)
        ppm_lim = math.ceil(max(tolp, np.nanpercentile(np.abs(after), 99.5) if len(after) else tolp) / 5) * 5 + 5
        diff = np.nanmedian(before - after) if len(after) else 0.0
        if abs(diff) > 1.0:
            ser = [("Before s/w calibration", before, C_ORANGE), ("After s/w calibration",
                                                                            after, C_BLUE)]
        else:
            ser = [("PSM precursor error", after, C_BLUE)]
        fig_ppm_hist(fpath("precursor_error_histogram"), ser, ppm_lim)
        fig_ppm_mz(fpath("precursor_error_vs_mz"), ps["m/z"].to_numpy(float), after, ppm_lim)
        fig6_note = None if abs(diff) > 1.0 else "The ppm column is not recalibrated (it agrees with the error computed from m/z)."
    else:
        fig6_note = None

    # Figure 3 Venn
    venn = None
    if len(samples) > 1 and 2 <= len(samples) <= 3:
        venn = {}
        pr, pp = exp.proteins, exp.peptides
        def present(df, s, key):
            for col in (f"#Spec {s}", f"#Feature {s}", f"Coverage(%) {s}", f"Area {s}"):
                if df is not None and col in df:
                    return set(df.loc[df[col].fillna(0).astype(float) > 0, key])
            return set()
        if pr is not None:
            fig_venn(fpath("venn_all_proteins"), [present(pr, s, "Accession") for s in samples], samples)
            top = pr[pr["Top"].map(_truthy)]
            fig_venn(fpath("venn_top_proteins"), [present(top, s, "Accession") for s in samples], samples)
        if pp is not None:
            fig_venn(fpath("venn_peptides"), [present(pp, s, "Peptide") for s in samples], samples)

    # ---- spectra and coverage images
    pep_df = exp.peptides
    spec_img = {}   # peptide -> img path
    spec_info = {"requested": 0, "drawn": 0, "not_found": 0, "failed": 0, "beyond_cap": 0}
    tol = float(filters.get("fragment_tolerance") or 0.5)
    tunit = "ppm" if str(filters.get("fragment_tolerance_unit", "Da")).lower() == "ppm" else "Da"
    fixed_deltas = set(filters.get("fixed_deltas") or ())
    names = _mod_name_map(exp, filters)
    pp = exp.protein_peptides
    prot_best = {}
    if pp is not None and len(pp):
        g = pp.assign(_u=pp["Unique"].map(_truthy)).sort_values(["_u", "-10LgP"], ascending=[False, False],
                                                               kind="mergesort")
        prot_best = dict(zip(*[g.drop_duplicates("Accession")[c] for c in ("Accession", "Peptide")]))
    if pep_df is not None and len(pep_df):
        pep_df = pep_df.sort_values("-10LgP", ascending=False, kind="mergesort").reset_index(drop=True)
    if images and readers and pep_df is not None and len(pep_df):
        top_peps = pep_df["Peptide"].tolist()[:max_spectra]
        spec_info["beyond_cap"] = max(0, len(pep_df) - max_spectra)
        want = list(dict.fromkeys(top_peps + [p for p in prot_best.values()]))
        spec_info["extra_for_proteins"] = len(want) - len(top_peps)
        rowmap = {p: i for i, p in enumerate(pep_df["Peptide"].tolist())}
        jobs = []
        single = len(readers) == 1
        for k, pep in enumerate(want):
            i = rowmap.get(pep)
            if i is None:
                continue
            r = pep_df.iloc[i]
            src = _norm_source(r.get("Source File", ""))
            rd = readers.get(src) or (next(iter(readers.values())) if single else None)
            spec_info["requested"] += 1
            if rd is None:
                spec_info["not_found"] += 1
                continue
            pk = rd.peaks(int(r["Scan"]))
            if pk is None or len(pk[0]) == 0:
                spec_info["not_found"] += 1
                continue
            try:
                stripped, mods = seqnorm.parse(pep)
            except ValueError:
                continue
            if any(a not in MONO for a in stripped):
                continue
            mods = [(p, _exact(d)) for p, d in mods]
            fname = f"spec_{k + 1:05d}.png"
            stitle = (f"{pep}    Scan {r.get('Source File', '')}:{int(r['Scan'])}   m/z={_trim(r.get('m/z'), 4)}"
                     f"   z={_int(r.get('z'))}   RT={_fixed(r.get('RT'), 2)}   ppm={_trim(r.get('ppm'), 1)}"
                     f"   -10lgP={_fixed(r.get('-10LgP'), 2)}")
            jobs.append({"path": os.path.join(img_dir, fname), "stripped": stripped, "mods": mods,
                         "z": int(r["z"]) if not _isna(r.get("z")) else 2, "mz": pk[0], "inten": pk[1],
                         "tol": tol, "unit": tunit, "title": stitle, "pmz": float(r.get("m/z") or 0),
                         "peptide": pep})
        _log(f"drawing {len(jobs)} spectra with {processes} processes")
        res = _run_jobs(render_spectrum, [{k: v for k, v in j.items() if k != "peptide"} for j in jobs], processes)
        for j, (path, err) in zip(jobs, res):
            if err:
                spec_info["failed"] += 1
                if spec_info["failed"] <= 3:
                    _log(f"spectrum {j['peptide']}: {err}")
            else:
                spec_img[j["peptide"]] = "img/" + os.path.basename(path)
        spec_info["drawn"] = len(spec_img)

    cov_img = {}
    sequences = {}
    fasta_used = fasta or (cfg.database if cfg is not None and cfg.database else None)
    if not fasta_used and os.path.exists(os.path.join(export_dir, "db.proteins.fasta")):
        fasta_used = os.path.join(export_dir, "db.proteins.fasta")
    pr = exp.proteins
    if pr is not None and len(pr):
        pr = pr.assign(_top=pr["Top"].map(_truthy) if "Top" in pr else True)
        pr = pr.sort_values(["Protein Group", "_top", "-10LgP"], ascending=[True, False, False],
                            kind="mergesort").reset_index(drop=True)
    if images and fasta_used and os.path.exists(fasta_used) and pr is not None and len(pr):
        _log(f"reading FASTA {fasta_used}")
        db = read_protein_db(fasta_used)
        for acc in pr["Accession"].tolist():
            i = db.index.get(acc)
            if i is not None:
                sequences[acc] = db.sequence[i]
        miss = len(pr) - len(sequences)
        if miss:
            notes.append(f"{miss} protein accession(s) are not in {os.path.basename(fasta_used)}: no coverage map.")
        mod_order = [n for _d, n in sorted(names.items()) if _d not in fixed_deltas]
        mod_colors = {n: MOD_COLORS[i % len(MOD_COLORS)] for i, n in enumerate(dict.fromkeys(mod_order))}
        by_acc = defaultdict(list)
        if pp is not None:
            for row in pp[["Accession", "Peptide", "Start", "Unique"]].itertuples(index=False):
                by_acc[row[0]].append(row)
        jobs = []
        for k, r in enumerate(pr.itertuples(index=False)):
            acc = getattr(r, "Accession")
            seq = sequences.get(acc)
            if not seq:
                continue
            peps = []
            for row in by_acc.get(acc, []):
                pep, st, uq = row[1], row[2], _truthy(row[3])
                try:
                    stripped, mods = seqnorm.parse(pep)
                except ValueError:
                    continue
                s0 = int(st) - 1 if not _isna(st) else seq.find(stripped)
                if s0 < 0:
                    continue
                n = len(stripped)
                ml = []
                for p, d in mods:
                    ds = seqnorm.format_delta(d)
                    if ds in fixed_deltas:
                        continue
                    ml.append((s0 + min(max(p, 1), n) - 1, names.get(ds, ds)))
                peps.append((s0, s0 + n, uq, ml))
            desc = pr.loc[k, "Description"] if "Description" in pr else ""
            jobs.append({"path": os.path.join(img_dir, f"cov_{k + 1:05d}.png"), "seq": seq, "peps": peps,
                         "header": f">{acc} {desc}", "mod_colors": mod_colors, "acc": acc})
        _log(f"drawing {len(jobs)} coverage maps")
        res = _run_jobs(render_coverage, [{k: v for k, v in j.items() if k != "acc"} for j in jobs], processes)
        for j, (path, err) in zip(jobs, res):
            if err:
                _log(f"coverage {j['acc']}: {err}")
            else:
                cov_img[j["acc"]] = "img/" + os.path.basename(path)
    elif images and pr is not None and len(pr):
        notes.append("No FASTA (--fasta): protein coverage maps are not drawn.")

    if readers and spec_info["beyond_cap"]:
        notes.append(f"Annotated spectra are drawn for the {max_spectra:,} highest-scoring peptides "
                     f"(plus the best unique PSM of each protein); {spec_info['beyond_cap']:,} lower-scoring "
                     "peptides are listed in db.peptides.html without a spectrum image.")
    if spec_info["not_found"]:
        notes.append(f"{spec_info['not_found']} best-PSM scans were not found in the mzML file(s).")

    # ---- write pages
    with open(os.path.join(out_dir, "export.css"), "w", encoding="utf-8") as fh:
        fh.write(CSS)
    summary_html = _summary_page(exp, samples, stats, stat_src, filters, pep_thr, t3, mc, figs, notes, lines6,
                                 param_src, sc, inputs, mzml, fasta_used, export_dir, fig6_note, venn,
                                 title or _default_title(export_dir))
    with open(os.path.join(out_dir, "db.summary.html"), "w", encoding="utf-8") as fh:
        fh.write(summary_html)
    with open(os.path.join(out_dir, "db.proteins.html"), "w", encoding="utf-8") as fh:
        fh.write(_proteins_page(pr, pp, pep_df, samples, cov_img, spec_img, prot_best))
    with open(os.path.join(out_dir, "db.peptides.html"), "w", encoding="utf-8") as fh:
        fh.write(_peptides_page(pep_df, samples, spec_img, spec_info, max_spectra))

    sizes = {f: os.path.getsize(os.path.join(out_dir, f)) for f in
             ("db.summary.html", "db.proteins.html", "db.peptides.html", "export.css")}
    img_files = os.listdir(img_dir)
    sizes["img"] = sum(os.path.getsize(os.path.join(img_dir, f)) for f in img_files)
    info = {"out": out_dir, "bytes": sizes, "n_images": len(img_files), "spectra": spec_info,
            "coverage_maps": len(cov_img), "table1": {k: (int(v) if isinstance(v, (int, np.integer)) else v)
                                                     for k, v in stats.items()},
            "table1_sources": stat_src, "table3": t3, "table5": mc,
            "peptide_threshold": pep_thr, "notes": notes, "elapsed_s": round(time.time() - t0, 1)}
    _log(f"wrote {out_dir} ({sum(sizes.values()) / 1e6:.1f} MB, {len(img_files)} images, {info['elapsed_s']} s)")
    return info


def _default_title(export_dir):
    b = os.path.basename(export_dir.rstrip(os.sep))
    if b.lower().startswith("export"):  # pipeline layout <run>/export
        b = os.path.basename(os.path.dirname(export_dir.rstrip(os.sep))) + "/" + b
    return b


def _summary_page(exp, samples, stats, stat_src, filters, pep_thr, t3, mc, figs, notes, lines6, param_src,
                  sc, inputs, mzml, fasta_used, export_dir, fig6_note, venn, title):
    out = []
    # 1. Notes
    out.append("<h1>1. Notes</h1>")
    src = [f"Export: {_esc(export_dir)}"]
    if mzml:
        src.append("mzML: " + ", ".join(_esc(p) for p in mzml))
    if inputs.get("psms_scored"):
        src.append(f"Scored targets/decoys: {_esc(inputs['psms_scored'])}")
    if inputs.get("precursors"):
        src.append("Precursor list: " + ", ".join(_esc(p) for p in inputs["precursors"]))
    elif inputs.get("mgf"):
        src.append("Precursor list (PEAKS MGF): " + ", ".join(_esc(p) for p in inputs["mgf"]))
    if inputs.get("features"):
        src.append("Feature table: " + ", ".join(_esc(p) for p in inputs["features"]))
    if fasta_used:
        src.append(f"FASTA: {_esc(fasta_used)}")
    if param_src:
        src.append(f"Search parameters from {param_src}.")
    out.append("<ul class=\"notes\">" + "".join(f"<li>{s}</li>" for s in src + [_esc(n) for n in notes]) + "</ul>")
    # 2. Result statistics
    out.append("<h1>2. Result Statistics</h1>")
    if "fdr_curve" in figs:
        out.append(_fig_section(1, "False discovery rate (FDR) curve. X axis is the number of peptides being kept. "
                                   "Y axis is the corresponding FDR.", _img(figs["fdr_curve"], "FDR curve", 800),
                                "Left: unique peptide sequences (I = L), the unit of the peptide FDR; right: PSMs "
                                "(best candidate per query). The dashed line marks the current -10lgP threshold."))
    else:
        out.append(_fig_section(1, "False discovery rate (FDR) curve. X axis is the number of peptides being kept. "
                                   "Y axis is the corresponding FDR.",
                                '<div class="na">Not available: needs scored target and decoy PSMs '
                                '(psms_scored.parquet).</div>'))
    f2 = []
    for key, sub, alt in (("score_histogram", "(a)", "score histogram"), ("score_vs_ppm", "(b)", "score vs ppm")):
        if key in figs:
            f2.append(f'<div><div class="sub">{sub}</div>{_img(figs[key], alt, 400)}</div>')
    out.append(_fig_section(2, "Score distribution. <b>(a)</b> Distribution of PEAKS peptide score; <b>(b)</b> "
                               "Scatterplot of PEAKS peptide score versus precursor mass error.",
                            f'<div class="figrow">{"".join(f2)}</div>' if f2 else '<div class="na">No PSMs.</div>',
                            None if sc is not None else "Exported (target) PSMs only; no decoys available."))
    out.append('<div class="section"><div class="caption"><b>Table 1.</b> Statistics of data.</div>'
               + html_table1(stats, samples) + "</div>")
    v = []
    for key, sub in (("venn_all_proteins", "(a)"), ("venn_top_proteins", "(b)"), ("venn_peptides", "(c)")):
        if venn is not None and key in figs:
            v.append(f'<div><div class="sub">{sub}</div>{_img(figs[key], key, 270)}</div>')
        elif len(samples) > 3:
            v.append(f'<div><div class="sub">{sub}</div><div class="na">Drawn for up to 3 samples</div></div>')
        else:
            v.append(f'<div><div class="sub">{sub}</div><div class="na">Not applicable to only one sample</div></div>')
    out.append(_fig_section(3, "Sample overlap for Proteins and Peptides (up to 8 samples). <b>(a)</b> All Proteins; "
                               "<b>(b)</b> Top Proteins; <b>(c)</b> Peptides;", f'<div class="figrow">{"".join(v)}</div>'))
    f4 = [f'<div>{_img(figs[k], k, 400)}</div>' for k in ("feature_mz", "feature_rt") if k in figs]
    out.append(_fig_section(4, "Distribution of peptide feature detection. <b>(a)</b> Feature m/z distribution; "
                               "<b>(b)</b> Feature RT distribution.",
                            f'<div class="figrow">{"".join(f4)}</div>' if f4 else '<div class="na">Not available.</div>',
                            "all = every detected MS1 feature; de novo+db = exported features with a DB or de novo "
                            "peptide; db = features identified by DB search." if f4 else None))
    out.append(_fig_section(5, "Distribution of identified peptide features. <b>(a)</b> Feature abundance distribution;",
                            _img(figs["feature_abundance"], "feature abundance", 400) if "feature_abundance" in figs
                            else '<div class="na">Not available.</div>'))
    t2 = [("Peptide -10lgP", f"&ge; {pep_thr:.2f}" if pep_thr is not None else NA),
          ("Protein -10lgP", f"&ge; {float(filters.get('protein_min_10lgp') or 15):g}"),
          ("Proteins unique peptides", f"&ge; {int(filters.get('protein_min_unique_peptides') or 1)}"),
          ("De novo score(%)", f"&ge; {float(filters.get('denovo_min_alc_pct') or 50):g}%")]
    t3r = [("FDR (Peptide-Spectrum Matches)", _esc(t3["psm"])), ("FDR (Peptide Sequences)", _esc(t3["seq"])),
           ("FDR (Protein Group)", _esc(t3["group"])), ("De Novo Only Spectra", _esc(t3["denovo"]))]
    t3note = ""
    if sc is not None:
        t3note = (f'<div class="foot">From {sc["n_queries"]:,} scored queries: PSMs {sc["n_decoy_psm"]}/'
                  f'{sc["n_target_psm"]:,} decoy/target; sequences {sc["n_decoy_seq"]}/{sc["n_target_seq"]:,}; '
                  f'protein groups {sc["n_decoy_groups"]}/{sc["n_target_groups"]:,}. PSMs: every PSM of an accepted sequence '
                  "down to -10lgP 20 (as exported); protein groups: decoy groups built with the same grouping rules "
                  "from decoy sequences above the threshold.</div>")
    out.append('<div class="section"><div class="figrow" style="text-align:left">'
               '<div style="text-align:left"><div class="caption"><b>Table 2.</b> Result filtration parameters.</div>'
               + _stats_table(t2) + '</div><div style="text-align:left"><div class="caption"><b>Table 3.</b> '
               "Statistics of filtered result.</div>" + _stats_table(t3r) + t3note + "</div></div></div>")
    ptm = exp.ptm_profile
    if ptm is not None and len(ptm):
        rows = "".join(
            f"<tr><td>{_esc(r['Name'])}</td><td style=\"text-align:right\">{_dec2(r['Mass (Da)'])}</td>"
            f"<td style=\"text-align:center\">{_esc(r['Position'])}</td><td style=\"text-align:right\">"
            f"{_int(r['# peptides'])}</td><td style=\"text-align:right\">{_int(r['# proteins'])}</td>"
            f"<td style=\"text-align:right\">{_int(r.get('# confident modification sites'))} "
            f"({_int(r.get('# modification sites'))})</td></tr>" for _i, r in ptm.iterrows())
        t4 = ("<table class=\"stats\"><tr><th>Name</th><th>&#8710;Mass</th><th>Position</th><th># peptides</th>"
              "<th># proteins</th><th># modification<br>sites</th></tr>" + rows + "</table>"
              '<div class="foot">Modification sites: confident (all).</div>')
    else:
        t4 = '<div class="na">No db.ptmProfile.csv in the export.</div>'
    out.append(f'<div class="section"><div class="caption"><b>Table 4.</b> PTM profile.</div>{t4}</div>')
    # 3. Experiment control
    out.append("<h1>3. Experiment Control</h1>")
    f6 = []
    for key, sub in (("precursor_error_histogram", "(a)"), ("precursor_error_vs_mz", "(b)")):
        if key in figs:
            f6.append(f'<div><div class="sub">{sub}</div>{_img(figs[key], key, 400)}</div>')
    out.append(_fig_section(6, "Precursor mass error of peptide-spectrum matches (PSM) in filtered result. <b>(a)</b> "
                               "Distribution of precursor mass error in ppm; <b>(b)</b> Scatterplot of precursor m/z "
                               "versus precursor mass error in ppm.",
                            f'<div class="figrow">{"".join(f6)}</div>' if f6 else '<div class="na">No PSMs.</div>',
                            fig6_note))
    hdr = "<tr><th>Missed Cleavages</th>" + "".join(f"<th>{k}</th>" for k in ("0", "1", "2", "3", "4+")) + "</tr>"
    body = "".join(f"<tr><td>{_esc(s)}</td>" + "".join(f"<td>{v}</td>" for v in row.values()) + "</tr>"
                   for s, row in mc.items())
    out.append('<div class="section"><div class="caption"><b>Table 5.</b> Number of identified peptides in each '
               f'sample by the number of missed cleavages.</div><table class="grid">{hdr}{body}</table></div>')
    # 4. Other information
    out.append("<h1>4. Other Information</h1>")
    items = []
    for label, val in lines6:
        if isinstance(val, list):
            items.append(f"{_esc(label)}: <br>" + "".join(f"&nbsp; {_esc(x)}<br>" for x in val))
        else:
            items.append(f"{_esc(label)}: {_esc(val)}<br>")
    out.append('<div class="section"><div class="caption"><b>Table 6.</b> Search parameters.</div>'
               f'<div class="align-left">{"".join(items) or "Not available."}</div></div>')
    return _page("Protein ID Summary", f"<div style=\"font-size:18px;font-weight:bold\">Protein ID Summary"
                 f" <span class=\"dim\" style=\"font-weight:normal\">{_esc(title)}</span></div>" + "".join(out))


def _ptm_list(v):
    if not isinstance(v, str) or not v:
        return ""
    return ", ".join(x.strip() for x in re.split(r";", v) if x.strip())


def _proteins_page(pr, pp, pep_df, samples, cov_img, spec_img, prot_best):
    s = samples[0]
    out = ['<div class="bar"><div class="tit"><a id="proteinTableLink">Protein List</a></div></div>']
    if pr is None or not len(pr):
        out.append('<div class="na">No db.proteins.csv in the export.</div>')
        return _page("proteins", "".join(out), wide=True)
    cols = ["Protein Group", "Protein ID", "Accession", "Gene", "-10lgP", "Coverage (%)", f"Coverage (%) {s}",
            f"Area {s}", "#Peptides", "#Unique", f"#Spec {s}", "PTM", "Avg. Mass", "Description"]
    rows = []
    for k, r in pr.iterrows():
        aid = f"p{k + 1}"
        rows.append(
            f'<tr class="{"" if r["_top"] else "dim"}"><td class="cc">{_int(r["Protein Group"])}</td><td class="cc">{k + 1}</td>'
            f'<td class="cl"><a href="#{aid}">{_esc(r["Accession"])}</a></td><td class="cc">{_esc(r.get("Gene", ""))}</td>'
            f'<td class="cc">{_trim(r["-10LgP"], 2)}</td><td class="cc">{_int(round(float(r["Coverage(%)"])))}</td>'
            f'<td class="cc">{_int(round(float(r.get(f"Coverage(%) {s}", r["Coverage(%)"]))))}</td>'
            f'<td class="cc">{_area(r.get(f"Area {s}"))}</td><td class="cc">{_int(r["#Peptides"])}</td>'
            f'<td class="cc">{_int(r["#Unique"])}</td><td class="cc">{_int(r.get(f"#Spec {s}"))}</td>'
            f'<td class="cc">{_esc(_ptm_list(r.get("PTM", "")))}</td><td class="cr">{_trim(r.get("Average Mass"), 4)}</td>'
            f'<td class="cl">{_esc(r.get("Description", ""))}</td></tr>')
    out.append('<table class="mytable"><thead><tr>' + "".join(f'<th class="cc">{_esc(c)}</th>' for c in cols)
               + "</tr></thead><tbody>" + "\n".join(rows) + "</tbody>"
               f'<tfoot><tr><th colspan="{len(cols)}" class="cl">total {len(pr)} proteins '
               "(non-top proteins of a group are dimmed)</th></tr></tfoot></table>")
    by_acc = defaultdict(list)
    if pp is not None:
        for _i, r in pp.iterrows():
            by_acc[r["Accession"]].append(r)
    pep_rows = {p: r for p, r in zip(pep_df["Peptide"], pep_df.to_dict("records"))} if pep_df is not None else {}
    pcols = ["Peptide", "Uniq", "-10lgP", "Mass", "Length", "ppm", "m/z", "z", "RT", "Scan", f"Area {s}",
             "#Feature", "Start", "End", "PTM", "AScore", "Found By"]
    for k, r in pr.iterrows():
        aid = f"p{k + 1}"
        acc = r["Accession"]
        out.append(f'<div class="bar"><div class="tit">{_esc(acc)}</div><a id="{aid}" href="#proteinTableLink">back to list</a></div>')
        out.append(f'<div class="meta">{_esc(r.get("Description", ""))}</div>')
        out.append(f'<div class="meta">Group {_int(r["Protein Group"])}{" (top)" if r["_top"] else ""} &middot; '
                   f'-10lgP {_trim(r["-10LgP"], 2)} &middot; coverage {_fixed(r["Coverage(%)"], 1)}% &middot; '
                   f'#peptides {_int(r["#Peptides"])} &middot; #unique {_int(r["#Unique"])}</div>')
        out.append(f'<div class="lst">| <a href="#cov{aid}">Protein Coverage</a> | <a href="#pep{aid}">Supporting '
                   f'Peptides</a> | <a href="#psm{aid}">Best Unique PSM</a> |</div>')
        out.append(f'<div class="lst"><a id="cov{aid}">Protein Coverage:</a></div>')
        out.append(_img(cov_img[acc], f"coverage {acc}", lazy=True) if acc in cov_img else '<div class="dim">No coverage map (sequence not available).</div>')
        out.append(f'<div class="lst"><a id="pep{aid}">Supporting Peptides:</a></div>')
        prs = sorted(by_acc.get(acc, []), key=lambda x: (float(x["Start"]) if not _isna(x.get("Start")) else 1e9,
                                                         -float(x["-10LgP"])))
        trs = []
        for x in prs:
            trs.append(
                f'<tr><td class="cl mono">{_esc(x.get("Extended Peptide") or x["Peptide"])}</td>'
                f'<td class="cc">{"Y" if _truthy(x.get("Unique")) else "N"}</td><td class="cc">{_trim(x["-10LgP"], 2)}</td>'
                f'<td class="cc">{_trim(x.get("Mass"), 4)}</td><td class="cc">{_int(x.get("Length"))}</td>'
                f'<td class="cc">{_trim(x.get("ppm"), 2)}</td><td class="cc">{_trim(x.get("m/z"), 4)}</td>'
                f'<td class="cc">{_int(x.get("z"))}</td><td class="cc">{_fixed(x.get("RT"), 2)}</td>'
                f'<td class="cc">{_int(x.get("Scan"))}</td><td class="cc">{_area(x.get(f"Area {s}"))}</td>'
                f'<td class="cc">{_int(x.get("#Feature"))}</td><td class="cc">{_int(x.get("Start"))}</td>'
                f'<td class="cc">{_int(x.get("End"))}</td><td class="cc">{_esc(_ptm_list(x.get("PTM", "")))}</td>'
                f'<td class="wt">{_esc(x.get("AScore", ""))}</td><td class="cc">{_esc(x.get("Found By", ""))}</td></tr>')
        out.append('<table class="mytable"><thead><tr>' + "".join(f'<th class="cc">{_esc(c)}</th>' for c in pcols)
                   + "</tr></thead><tbody>" + "".join(trs) + "</tbody>"
                   f'<tfoot><tr><th colspan="{len(pcols)}" class="cl">total {len(prs)} peptides</th></tr></tfoot></table>')
        best = prot_best.get(acc)
        pr_ = pep_rows.get(best)
        if pr_ is not None:
            line = (f"Best Unique PSM ({_esc(best)}, Scan {_esc(pr_.get('Source File', ''))}:{_int(pr_.get('Scan'))}, "
                    f"m/z={_trim(pr_.get('m/z'), 4)}, z={_int(pr_.get('z'))}, RT={_fixed(pr_.get('RT'), 2)}, "
                    f"ppm={_trim(pr_.get('ppm'), 2)}):")
            out.append(f'<div class="lst"><a id="psm{aid}">{line}</a></div>')
            out.append(_img(spec_img[best], f"spectrum {best}", lazy=True) if best in spec_img
                       else '<div class="dim">No spectrum image (no --mzml, or scan not found).</div>')
        else:
            out.append(f'<div class="lst"><a id="psm{aid}">Best Unique PSM:</a></div><div class="dim">None.</div>')
    return _page("proteins", "\n".join(out), wide=True)


def _peptides_page(pep_df, samples, spec_img, spec_info, max_spectra):
    s = samples[0]
    out = ['<div class="bar"><div class="tit"><a id="peptideTableLink">Peptide List</a></div></div>']
    if pep_df is None or not len(pep_df):
        out.append('<div class="na">No db.peptides.csv in the export.</div>')
        return _page("peptide list", "".join(out), wide=True)
    n_img = sum(1 for p in pep_df["Peptide"] if p in spec_img)
    msg = f"{n_img:,} of {len(pep_df):,} peptides have an annotated spectrum of their best PSM (click the peptide)."
    if spec_info.get("beyond_cap"):
        msg += (f" Images are drawn for the {max_spectra:,} highest-scoring peptides only (plus proteins' best "
                f"unique PSMs) to keep the img folder small; the remaining peptides have no image.")
    if not spec_img:
        msg = "No spectrum images (no --mzml given)."
    out.append(f'<div class="meta">{_esc(msg)}</div>')
    cols = ["Peptide", "-10lgP", "CAA (%)", "Mass", "Length", "ppm", "m/z", "z", "RT", "delta RT", "MS2 Correlation",
            "Scan", f"Area {s}", "Accession", "PTM", "AScore", "Found By"]
    rows, blocks = [], []
    for k, r in enumerate(pep_df.to_dict("records")):
        pep = r["Peptide"]
        aid = f"s{k + 1}"
        cell = f'<a href="#{aid}">{_esc(pep)}</a>' if pep in spec_img else _esc(pep)
        acc = str(r.get("Accession") or "")
        rows.append(
            f'<tr><td class="cl mono">{cell}</td><td class="cc">{_trim(r["-10LgP"], 2)}</td>'
            f'<td class="cc">{_fixed(r.get("CAA (%)"), 1)}</td><td class="cc">{_trim(r.get("Mass"), 4)}</td>'
            f'<td class="cc">{_int(r.get("Length"))}</td><td class="cc">{_trim(r.get("ppm"), 2)}</td>'
            f'<td class="cc">{_trim(r.get("m/z"), 4)}</td><td class="cc">{_int(r.get("z"))}</td>'
            f'<td class="cc">{_fixed(r.get("RT"), 2)}</td><td class="cc">{_trim(r.get("Delta RT"), 2)}</td>'
            f'<td class="cc">{_trim(r.get("MS2 Correlation"), 2)}</td><td class="cc">{_int(r.get("Scan"))}</td>'
            f'<td class="cc">{_area(r.get(f"Area {s}"))}</td><td class="wt">{_esc(acc)}</td>'
            f'<td class="cc">{_esc(_ptm_list(r.get("PTM", "")))}</td><td class="wt">{_esc(r.get("AScore", ""))}</td>'
            f'<td class="cc">{_esc(r.get("Found By", ""))}</td></tr>')
        if pep in spec_img:
            blocks.append(
                f'<div class="bar"><div class="tit mono">{_esc(pep)}</div><a id="{aid}" href="#peptideTableLink">'
                f'back to list</a></div><div class="lst">Best PSM (Scan {_esc(r.get("Source File", ""))}:'
                f'{_int(r.get("Scan"))}, m/z={_trim(r.get("m/z"), 4)}, z={_int(r.get("z"))}, RT={_fixed(r.get("RT"), 2)}, '
                f'ppm={_trim(r.get("ppm"), 2)}, -10lgP={_trim(r["-10LgP"], 2)}):</div>'
                + _img(spec_img[pep], f"spectrum {pep}", lazy=True))
    out.append('<table class="mytable"><thead><tr>' + "".join(f'<th class="cc">{_esc(c)}</th>' for c in cols)
               + "</tr></thead><tbody>" + "\n".join(rows) + "</tbody>"
               f'<tfoot><tr><th colspan="{len(cols)}" class="cl">total {len(pep_df)} peptides</th></tr></tfoot></table>')
    out.extend(blocks)
    return _page("peptide list", "\n".join(out), wide=True)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m proteoscan.summary",
                                 description="Write PEAKS-style HTML export pages (db.summary.html, db.proteins.html, "
                                             "db.peptides.html, img/, export.css) for a PEAKS-format export.")
    ap.add_argument("--export", required=True, help="export directory (PEAKS DB-search CSV layout)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--mzml", nargs="*", default=[], help="indexed mzML file(s): MS1/MS2 counts and spectra")
    ap.add_argument("--fasta", help="protein FASTA for coverage maps (default: config database, "
                                    "else <export>/db.proteins.fasta)")
    ap.add_argument("--psms-scored", help="psms_scored.parquet (targets + decoys) for Figure 1, 2 and Table 3")
    ap.add_argument("--features", nargs="*", default=None, help="ms1_features.parquet file(s) (all detected features)")
    ap.add_argument("--precursors", nargs="*", default=None, help="precursors.parquet file(s) (#Chimera)")
    ap.add_argument("--config", help="pipeline search config YAML (Table 2/6)")
    ap.add_argument("--mgf", nargs="*", default=None, help="PEAKS MGF (#Chimera when no precursor list)")
    ap.add_argument("--run-dir", help="pipeline run directory (default: the export's parent, if it is one)")
    ap.add_argument("--max-spectra", type=int, default=2000, help="spectrum images for the top-N peptides")
    ap.add_argument("--processes", type=int, default=8, help="worker processes for the images (<= 8)")
    ap.add_argument("--no-images", action="store_true", help="skip coverage maps and spectra")
    ap.add_argument("--title", help="page title (default: export directory name)")
    a = ap.parse_args(argv)
    info = build_summary(a.export, a.out, mzml=a.mzml, fasta=a.fasta, psms_scored=a.psms_scored,
                         features=a.features, precursors=a.precursors, config=a.config,
                         mgf=a.mgf if a.mgf else None, run_dir=a.run_dir, max_spectra=a.max_spectra,
                         processes=max(1, min(8, a.processes)), title=a.title, images=not a.no_images)
    print(json.dumps(info, default=str, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
