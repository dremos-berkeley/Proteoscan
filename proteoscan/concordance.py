"""Concordance between a reference PEAKS export and a test export in PEAKS format.

Usage::

    python -m proteoscan.concordance --ref <peaks_export_dir> --test <our_export_dir> --out <dir>
        [--features] [--feature-ppm 10] [--feature-rt 0.5] [--ignore-source-file]
        [--fixed-mod C:+57.02] [--merge-nterm] [--min-10lgp FLOAT|ref]

Levels and keys
---------------
psm                   (source file, raw scan, canonical modified sequence, z); rows
                      are de-duplicated on this key (PEAKS repeats a PSM once
                      per matched feature)
scan                  (source file, raw scan) identified at all (diagnostic)
peptide               canonical modified sequence (from db.peptides.csv, or
                      from PSMs when that file is absent)
peptide_sequence      stripped sequence (diagnostic)
peptide_IL_collapsed  canonical modified sequence with I->L (diagnostic; this
                      is what PEAKS counts as "#Peptides" in db.summary.html)
protein_group         groups matched by any shared accession (connected
                      components of the ref-group / test-group / accession graph)
protein_top_accession accessions flagged Top=true
protein_accession     all accessions in db.proteins.csv (diagnostic)
feature (optional)    same source file and z, m/z within 10 ppm, RT within 0.5 min,
                      one-to-one greedy matching by normalised distance

For each level: n_ref, n_test, shared, ref_only, test_only,
recall = shared/n_ref, precision = shared/n_test, jaccard.  Exact matching
everywhere: I != L and modification sites matter.

The headline is PSM recall/precision with both sides taken as exported, i.e.
at the reference's own filter (PEAKS exports every PSM of each peptide that
passed the peptide-level FDR threshold, including PSMs scoring below it).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd

from . import seqnorm
from .peaks_export import PeaksExport, read_export

__all__ = ["compare", "write_report", "main", "set_metrics", "norm_source", "norm_accession"]

_SRC_EXT = re.compile(r"(\.(raw|mzml|mzxml|mgf|d|wiff|ms2|gz|bz2|zip))+$", re.I)


def norm_source(name, ignore=False):
    """Source file key: basename without raw/mzML/mgf/... extension(s)."""
    if ignore:
        return ""
    base = os.path.basename(str(name).replace("\\", "/"))
    return _SRC_EXT.sub("", base)


_ACC_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|(.+)$")


def norm_accession(acc):
    """``sp|P02768|ALBU_HUMAN`` -> ``P02768|ALBU_HUMAN`` (PEAKS style); others unchanged."""
    a = str(acc).strip()
    m = _ACC_RE.match(a)
    return f"{m.group(1)}|{m.group(2)}" if m else a


def _ratio(a, b):
    return float(a) / b if b else None


def set_metrics(ref, test):
    ref, test = set(ref), set(test)
    shared = len(ref & test)
    union = len(ref | test)
    return {
        "n_ref": len(ref),
        "n_test": len(test),
        "shared": shared,
        "ref_only": len(ref - test),
        "test_only": len(test - ref),
        "recall": _ratio(shared, len(ref)),
        "precision": _ratio(shared, len(test)),
        "jaccard": _ratio(shared, union),
    }


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


class _Canon:
    """Memoised canonicaliser that records unparseable strings instead of failing."""

    def __init__(self, fixed=None, merge_nterm=False):
        self.fixed = fixed or None
        self.merge_nterm = merge_nterm
        self.cache = {}
        self.failed = Counter()

    def __call__(self, s):
        s = "" if s is None or (isinstance(s, float) and s != s) else str(s)
        r = self.cache.get(s)
        if r is None:
            try:
                r = seqnorm.canonical(s, fixed=self.fixed, nterm_to_first_residue=self.merge_nterm)
            except ValueError:
                self.failed[s] += 1
                r = "UNPARSED:" + s
            self.cache[s] = r
        return r


def _strip(modseq):
    if modseq.startswith("UNPARSED:"):
        return modseq
    return seqnorm.strip_mods(modseq)


def _require(df, cols, what):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{what}: missing required column(s) {missing}")


def _scan_int(series):
    if pd.api.types.is_integer_dtype(series.dtype):
        return series.astype("int64")
    s = series.astype(str).str.replace(r"^.*:", "", regex=True)  # "F1:1234" style
    return pd.to_numeric(s, errors="coerce").fillna(-1).astype("int64")


def _psm_frame(df, canon, ignore_source, min_score, what):
    cols = ["source", "scan", "peptide", "modseq", "stripped", "z", "score"]
    if df is None:
        return pd.DataFrame({c: [] for c in cols})
    _require(df, ["Peptide", "Scan", "z", "Source File"], what)
    out = pd.DataFrame({
        "source": [norm_source(x, ignore_source) for x in df["Source File"].tolist()],
        "scan": _scan_int(df["Scan"]).to_numpy(),
        "peptide": df["Peptide"].astype(str).to_numpy(),
        "z": pd.to_numeric(df["z"], errors="coerce").fillna(0).astype("int64").to_numpy(),
        "score": (pd.to_numeric(df["-10LgP"], errors="coerce").to_numpy()
                  if "-10LgP" in df.columns else np.full(len(df), np.nan)),
    })
    out["modseq"] = [canon(p) for p in out["peptide"].tolist()]
    out["stripped"] = [_strip(m) for m in out["modseq"].tolist()]
    if min_score is not None:
        out = out[out["score"] >= min_score]
    return out[cols].reset_index(drop=True)


PSM_CATEGORIES = [
    "same_scan_same_peptide_other_charge",
    "same_scan_IL_variant",
    "same_scan_same_sequence_other_mods",
    "same_scan_IL_variant_other_mods",
    "same_scan_other_peptide",
    "other_scan_same_peptide",
    "other_scan_same_sequence",
    "not_found",
]
SAME_SCAN_DIFFERENT_PEPTIDE = PSM_CATEGORIES[1:5]


def _psm_index(fr):
    by_scan = {}
    for s, sc, m, z, st in zip(fr["source"], fr["scan"], fr["modseq"], fr["z"], fr["stripped"]):
        by_scan.setdefault((s, int(sc)), []).append((m, int(z), st))
    return {
        "by_scan": by_scan,
        "modseqs": set(fr["modseq"]),
        "stripped": set(fr["stripped"]),
    }


def _classify_psm(src, scan, modseq, z, stripped, other):
    entries = other["by_scan"].get((src, int(scan)))
    if entries:
        if any(m == modseq for m, _z, _s in entries):
            return "same_scan_same_peptide_other_charge"
        il = seqnorm.il_collapse(modseq)
        if any(seqnorm.il_collapse(m) == il for m, _z, _s in entries):
            return "same_scan_IL_variant"
        if any(s == stripped for _m, _z, s in entries):
            return "same_scan_same_sequence_other_mods"
        ils = seqnorm.il_collapse(stripped)
        if any(seqnorm.il_collapse(s) == ils for _m, _z, s in entries):
            return "same_scan_IL_variant_other_mods"
        return "same_scan_other_peptide"
    if modseq in other["modseqs"]:
        return "other_scan_same_peptide"
    if stripped in other["stripped"]:
        return "other_scan_same_sequence"
    return "not_found"


def _other_peptides_at_scan(src, scan, other):
    entries = other["by_scan"].get((src, int(scan))) or []
    return ";".join(sorted({f"{m}/{z}" for m, z, _s in entries}))


def _psm_only_table(fr, keys_only, other):
    key = list(zip(fr["source"], fr["scan"], fr["modseq"], fr["z"]))
    mask = np.array([k in keys_only for k in key], dtype=bool)
    sub = fr[mask].copy()
    # one row per key; keep the best-scoring row
    sub = sub.sort_values("score", ascending=False, kind="stable").drop_duplicates(
        ["source", "scan", "modseq", "z"])
    sub["category"] = [
        _classify_psm(s, sc, m, z, st, other)
        for s, sc, m, z, st in zip(sub["source"], sub["scan"], sub["modseq"], sub["z"], sub["stripped"])
    ]
    sub["other_side_at_scan"] = [
        _other_peptides_at_scan(s, sc, other) for s, sc in zip(sub["source"], sub["scan"])
    ]
    return sub.sort_values(["source", "scan", "modseq", "z"]).reset_index(drop=True)


def _peptide_list(exp, canon):
    """Canonical peptides of an export: db.peptides.csv, else PSM peptides."""
    if exp.peptides is not None:
        _require(exp.peptides, ["Peptide"], "peptides")
        df = exp.peptides
        src = "peptides"
    elif exp.psms is not None:
        df = exp.psms
        src = "psms"
    else:
        return pd.DataFrame({"peptide": [], "modseq": [], "score": []}), "none"
    out = pd.DataFrame({
        "peptide": df["Peptide"].astype(str).to_numpy(),
        "score": (pd.to_numeric(df["-10LgP"], errors="coerce").to_numpy()
                  if "-10LgP" in df.columns else np.full(len(df), np.nan)),
    })
    out["modseq"] = [canon(p) for p in out["peptide"].tolist()]
    out = out.sort_values("score", ascending=False, kind="stable").drop_duplicates("modseq")
    return out.reset_index(drop=True), src


def _classify_peptide(modseq, other_modseqs_il, other_stripped, other_stripped_il):
    st = _strip(modseq)
    if seqnorm.il_collapse(modseq) in other_modseqs_il:
        return "IL_variant"
    if st in other_stripped:
        return "same_sequence_other_mods"
    if seqnorm.il_collapse(st) in other_stripped_il:
        return "IL_variant_other_mods"
    return "not_found"


def _peptide_only_table(pl, only, other_modseqs):
    oil = {seqnorm.il_collapse(m) for m in other_modseqs}
    ost = {_strip(m) for m in other_modseqs}
    ostil = {seqnorm.il_collapse(s) for s in ost}
    sub = pl[pl["modseq"].isin(only)].copy()
    sub["stripped"] = [_strip(m) for m in sub["modseq"]]
    sub["category"] = [_classify_peptide(m, oil, ost, ostil) for m in sub["modseq"]]
    return sub[["modseq", "peptide", "stripped", "score", "category"]].sort_values("modseq").reset_index(drop=True)


# ---- proteins ------------------------------------------------------------


def _protein_groups(df):
    """{group id: {"accessions": set, "top": set}} from db.proteins.csv."""
    _require(df, ["Accession"], "proteins")
    gid = df["Protein Group"].tolist() if "Protein Group" in df.columns else list(range(len(df)))
    top = df["Top"].tolist() if "Top" in df.columns else [True] * len(df)
    groups = {}
    for g, a, t in zip(gid, df["Accession"].tolist(), top):
        acc = norm_accession(a)
        e = groups.setdefault(g, {"accessions": set(), "top": set()})
        e["accessions"].add(acc)
        if t is True or str(t).lower() == "true":
            e["top"].add(acc)
    return groups


def _protein_group_compare(rg, tg):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for side, groups in (("R", rg), ("T", tg)):
        for g, e in groups.items():
            node = (side, g)
            find(node)
            for a in e["accessions"]:
                union(node, ("A", a))
    comps = {}
    for side, groups in (("R", rg), ("T", tg)):
        for g in groups:
            comps.setdefault(find((side, g)), []).append((side, g))
    mixed = [c for c in comps.values() if {s for s, _ in c} == {"R", "T"}]
    ref_matched = {g for c in mixed for s, g in c if s == "R"}
    test_matched = {g for c in mixed for s, g in c if s == "T"}
    ref_only = [g for g in rg if g not in ref_matched]
    test_only = [g for g in tg if g not in test_matched]
    n_comp = len(mixed) + len(ref_only) + len(test_only)
    test_sets = {frozenset(e["accessions"]) for e in tg.values()}
    metrics = {
        "n_ref": len(rg),
        "n_test": len(tg),
        "shared": len(mixed),
        "ref_only": len(ref_only),
        "test_only": len(test_only),
        "recall": _ratio(len(ref_matched), len(rg)),
        "precision": _ratio(len(test_matched), len(tg)),
        "jaccard": _ratio(len(mixed), n_comp),
        "ref_groups_matched": len(ref_matched),
        "test_groups_matched": len(test_matched),
        "one_to_one_components": sum(1 for c in mixed if len(c) == 2),
        "ref_groups_identical_accession_set": sum(
            1 for e in rg.values() if frozenset(e["accessions"]) in test_sets),
        "note": "shared = connected components containing both ref and test groups; "
                "recall/precision = fraction of ref/test groups sharing any accession",
    }

    def table(groups, ids):
        return pd.DataFrame({
            "protein_group": list(ids),
            "top_accessions": [";".join(sorted(groups[g]["top"])) for g in ids],
            "accessions": [";".join(sorted(groups[g]["accessions"])) for g in ids],
            "n_accessions": [len(groups[g]["accessions"]) for g in ids],
        })

    return metrics, table(rg, ref_only), table(tg, test_only)


# ---- features ------------------------------------------------------------


def _feature_frame(df, ignore_source):
    _require(df, ["m/z", "z", "RT"], "features")
    src = df["Source File"].tolist() if "Source File" in df.columns else [""] * len(df)
    return pd.DataFrame({
        "source": [norm_source(s, ignore_source) for s in src],
        "feature_id": df["Feature Id"].to_numpy() if "Feature Id" in df.columns else np.arange(len(df)),
        "mz": pd.to_numeric(df["m/z"], errors="coerce").to_numpy(),
        "z": pd.to_numeric(df["z"], errors="coerce").fillna(0).astype("int64").to_numpy(),
        "rt": pd.to_numeric(df["RT"], errors="coerce").to_numpy(),
        "area": (pd.to_numeric(df["Area"], errors="coerce").to_numpy()
                 if "Area" in df.columns else np.full(len(df), np.nan)),
    })


def _match_features(rf, tf, ppm_tol, rt_tol):
    """Greedy one-to-one matching; returns array of (ref_idx, test_idx, dppm, drt)."""
    cands = []
    for (src, z), rsub in rf.groupby(["source", "z"], sort=False):
        tsub = tf[(tf["source"] == src) & (tf["z"] == z)]
        if tsub.empty:
            continue
        order = np.argsort(tsub["mz"].to_numpy(), kind="stable")
        tmz = tsub["mz"].to_numpy()[order]
        trt = tsub["rt"].to_numpy()[order]
        tidx = tsub.index.to_numpy()[order]
        rmz = rsub["mz"].to_numpy()
        rrt = rsub["rt"].to_numpy()
        lo = np.searchsorted(tmz, rmz * (1 - ppm_tol * 1e-6), side="left")
        hi = np.searchsorted(tmz, rmz * (1 + ppm_tol * 1e-6), side="right")
        for i, ridx in enumerate(rsub.index.to_numpy()):
            if lo[i] >= hi[i]:
                continue
            sl = slice(lo[i], hi[i])
            dppm = (tmz[sl] - rmz[i]) / rmz[i] * 1e6
            drt = trt[sl] - rrt[i]
            ok = (np.abs(dppm) <= ppm_tol) & (np.abs(drt) <= rt_tol)
            for tj, dp, dr in zip(tidx[sl][ok], dppm[ok], drt[ok]):
                cands.append(((dp / ppm_tol) ** 2 + (dr / rt_tol) ** 2, ridx, tj, dp, dr))
    cands.sort(key=lambda c: (c[0], c[1], c[2]))
    used_r, used_t, pairs = set(), set(), []
    for _d, r, t, dp, dr in cands:
        if r in used_r or t in used_t:
            continue
        used_r.add(r)
        used_t.add(t)
        pairs.append((r, t, dp, dr))
    return pairs


def _feature_compare(ref_df, test_df, ignore_source, ppm_tol, rt_tol):
    rf = _feature_frame(ref_df, ignore_source)
    tf = _feature_frame(test_df, ignore_source)
    pairs = _match_features(rf, tf, ppm_tol, rt_tol)
    shared = len(pairs)
    mr = {p[0] for p in pairs}
    mt = {p[1] for p in pairs}
    n_r, n_t = len(rf), len(tf)
    m = {
        "n_ref": n_r,
        "n_test": n_t,
        "shared": shared,
        "ref_only": n_r - shared,
        "test_only": n_t - shared,
        "recall": _ratio(shared, n_r),
        "precision": _ratio(shared, n_t),
        "jaccard": _ratio(shared, n_r + n_t - shared),
        "ppm_tolerance": ppm_tol,
        "rt_tolerance_min": rt_tol,
        "median_ppm_diff": float(np.median([p[2] for p in pairs])) if pairs else None,
        "median_rt_diff_min": float(np.median([p[3] for p in pairs])) if pairs else None,
    }
    return m, rf[~rf.index.isin(mr)].reset_index(drop=True), tf[~tf.index.isin(mt)].reset_index(drop=True)


# --------------------------------------------------------------------------
# main comparison
# --------------------------------------------------------------------------

_REF_PARAM_KEYS = ["peptide_min_10lgp", "protein_min_10lgp", "protein_min_unique_peptides",
                   "denovo_min_alc_pct", "fdr_psm_pct", "fdr_peptide_pct", "fdr_protein_group_pct",
                   "enzyme", "digest_mode", "max_missed_cleavages", "precursor_tolerance",
                   "precursor_tolerance_unit", "fragment_tolerance", "fragment_tolerance_unit",
                   "fixed_modifications", "variable_modifications", "database", "searched_entries",
                   "deep_learning_boost"]


def _load(x, tables):
    if isinstance(x, PeaksExport):
        return x
    return read_export(x, tables=tables, verify=False)


def compare(ref, test, features=False, ignore_source=False, fixed=None, merge_nterm=False,
            min_10lgp=None, feature_ppm=10.0, feature_rt=0.5):
    """Compare two exports (directories or :class:`PeaksExport`).

    ``min_10lgp``: None (as exported), a float, or "ref" (the reference's
    peptide -10lgP threshold), applied to PSMs of both sides.
    Returns ``{"metrics": dict, "tables": {name: DataFrame}}``.
    """
    tables = ["psms", "peptides", "proteins"] + (["features"] if features else [])
    R = _load(ref, tables)
    T = _load(test, tables)
    warnings = []
    canon = _Canon(fixed=fixed, merge_nterm=merge_nterm)

    thr = min_10lgp
    if isinstance(thr, str):
        if thr != "ref":
            thr = float(thr)
        else:
            thr = R.params.get("peptide_min_10lgp")
            if thr is None:
                warnings.append("--min-10lgp ref requested but reference has no db.summary.html threshold")

    levels = {}
    out_tables = {}

    # ---- PSMs
    if R.psms is None:
        raise ValueError(f"reference has no db.psms.csv ({R.source_dir})")
    if T.psms is None:
        warnings.append("test export has no db.psms.csv")
    rp = _psm_frame(R.psms, canon, ignore_source, thr, "ref psms")
    tp = _psm_frame(T.psms, canon, ignore_source, thr, "test psms")
    rkeys = set(zip(rp["source"], rp["scan"], rp["modseq"], rp["z"]))
    tkeys = set(zip(tp["source"], tp["scan"], tp["modseq"], tp["z"]))
    psm = set_metrics(rkeys, tkeys)
    psm["n_ref_rows"] = len(rp)
    psm["n_test_rows"] = len(tp)
    levels["psm"] = psm
    rindex, tindex = _psm_index(rp), _psm_index(tp)
    ref_only_psm = _psm_only_table(rp, rkeys - tkeys, tindex)
    test_only_psm = _psm_only_table(tp, tkeys - rkeys, rindex)
    out_tables["psm_ref_only"] = ref_only_psm
    out_tables["psm_test_only"] = test_only_psm
    rb = Counter(ref_only_psm["category"])
    tb = Counter(test_only_psm["category"])
    ref_breakdown = {c: int(rb.get(c, 0)) for c in PSM_CATEGORIES}
    test_breakdown = {c: int(tb.get(c, 0)) for c in PSM_CATEGORIES}

    # ---- scans
    rs = set(zip(rp["source"], rp["scan"]))
    ts = set(zip(tp["source"], tp["scan"]))
    levels["scan"] = set_metrics(rs, ts)
    agree = agree_il = 0
    for k in rs & ts:
        a = {m for m, _z, _s in rindex["by_scan"][k]}
        b = {m for m, _z, _s in tindex["by_scan"][k]}
        if a & b:
            agree += 1
        elif {seqnorm.il_collapse(x) for x in a} & {seqnorm.il_collapse(x) for x in b}:
            agree_il += 1
    n_shared_scans = len(rs & ts)
    scan_agreement = {
        "shared_scans": n_shared_scans,
        "same_peptide": agree,
        "IL_variant_only": agree_il,
        "different_peptide": n_shared_scans - agree - agree_il,
    }

    # ---- peptides
    rpl, rsrc = _peptide_list(R, canon)
    tpl, tsrc = _peptide_list(T, canon)
    rpep, tpep = set(rpl["modseq"]), set(tpl["modseq"])
    levels["peptide"] = set_metrics(rpep, tpep)
    levels["peptide"]["source_ref"] = rsrc
    levels["peptide"]["source_test"] = tsrc
    levels["peptide_sequence"] = set_metrics({_strip(m) for m in rpep}, {_strip(m) for m in tpep})
    levels["peptide_IL_collapsed"] = set_metrics({seqnorm.il_collapse(m) for m in rpep},
                                                 {seqnorm.il_collapse(m) for m in tpep})
    out_tables["peptide_ref_only"] = _peptide_only_table(rpl, rpep - tpep, tpep)
    out_tables["peptide_test_only"] = _peptide_only_table(tpl, tpep - rpep, rpep)

    # ---- proteins
    if R.proteins is not None and T.proteins is not None:
        rg, tg = _protein_groups(R.proteins), _protein_groups(T.proteins)
        m, ro, to = _protein_group_compare(rg, tg)
        levels["protein_group"] = m
        out_tables["protein_group_ref_only"] = ro
        out_tables["protein_group_test_only"] = to
        rtop = set().union(*(e["top"] for e in rg.values())) if rg else set()
        ttop = set().union(*(e["top"] for e in tg.values())) if tg else set()
        levels["protein_top_accession"] = set_metrics(rtop, ttop)
        out_tables["protein_top_ref_only"] = pd.DataFrame({"accession": sorted(rtop - ttop)})
        out_tables["protein_top_test_only"] = pd.DataFrame({"accession": sorted(ttop - rtop)})
        rall = set().union(*(e["accessions"] for e in rg.values())) if rg else set()
        tall = set().union(*(e["accessions"] for e in tg.values())) if tg else set()
        levels["protein_accession"] = set_metrics(rall, tall)
    else:
        warnings.append("protein levels skipped: db.proteins.csv missing on "
                        + ("ref" if R.proteins is None else "test") + " side")

    # ---- features
    if features:
        if R.features is not None and T.features is not None:
            m, ro, to = _feature_compare(R.features, T.features, ignore_source, feature_ppm, feature_rt)
            levels["feature"] = m
            out_tables["feature_ref_only"] = ro
            out_tables["feature_test_only"] = to
        else:
            warnings.append("feature level skipped: features.csv missing on "
                            + ("ref" if R.features is None else "test") + " side")

    if canon.failed:
        warnings.append(f"{len(canon.failed)} peptide strings could not be parsed "
                        f"(e.g. {list(canon.failed)[:3]}); kept as UNPARSED:<string>")

    metrics = {
        "ref": R.source_dir,
        "test": T.source_dir,
        "options": {
            "features": features, "ignore_source_file": ignore_source, "fixed_mods": fixed,
            "merge_nterm": merge_nterm, "min_10lgp": min_10lgp, "min_10lgp_applied": thr,
            "feature_ppm": feature_ppm, "feature_rt_min": feature_rt,
        },
        "ref_params": {k: R.params.get(k) for k in _REF_PARAM_KEYS if k in R.params},
        "ref_statistics": R.params.get("statistics", {}),
        "headline": {
            "psm_recall": psm["recall"],
            "psm_precision": psm["precision"],
            "psm_jaccard": psm["jaccard"],
            "psm_shared": psm["shared"],
            "psm_n_ref": psm["n_ref"],
            "psm_n_test": psm["n_test"],
        },
        "levels": levels,
        "psm_same_scan_different_peptide": {
            "ref_only": sum(ref_breakdown[c] for c in SAME_SCAN_DIFFERENT_PEPTIDE),
            "test_only": sum(test_breakdown[c] for c in SAME_SCAN_DIFFERENT_PEPTIDE),
        },
        "psm_ref_only_breakdown": ref_breakdown,
        "psm_test_only_breakdown": test_breakdown,
        "scan_agreement": scan_agreement,
        "warnings": warnings,
    }
    return {"metrics": metrics, "tables": out_tables}


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------


def _pct(x):
    return "n/a" if x is None else f"{100 * x:.2f}%"


_LEVEL_ORDER = ["psm", "scan", "peptide", "peptide_sequence", "peptide_IL_collapsed",
                "protein_group", "protein_top_accession", "protein_accession", "feature"]


def _markdown(metrics):
    L = []
    h = metrics["headline"]
    L.append("# Concordance report\n")
    L.append(f"- Reference: `{metrics['ref']}`")
    L.append(f"- Test: `{metrics['test']}`")
    opts = {k: v for k, v in metrics["options"].items() if v not in (None, False)}
    L.append(f"- Options: `{json.dumps(opts)}`\n")
    L.append("## Headline (PSM level, reference's own filter)\n")
    L.append(f"**PSM recall {_pct(h['psm_recall'])}**, **precision {_pct(h['psm_precision'])}** "
             f"({h['psm_shared']} shared of {h['psm_n_ref']} ref / {h['psm_n_test']} test PSMs; "
             f"Jaccard {_pct(h['psm_jaccard'])}).\n")
    L.append("## Levels\n")
    L.append("| level | ref | test | shared | ref-only | test-only | recall | precision | Jaccard |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name in _LEVEL_ORDER:
        m = metrics["levels"].get(name)
        if not m:
            continue
        L.append(f"| {name} | {m['n_ref']} | {m['n_test']} | {m['shared']} | {m['ref_only']} | "
                 f"{m['test_only']} | {_pct(m['recall'])} | {_pct(m['precision'])} | {_pct(m['jaccard'])} |")
    L.append("")
    L.append("PSM key: source file (extension stripped) + raw scan + canonical modified sequence + charge. "
             "Peptide key: canonical modified sequence (I != L, sites matter). Protein groups match on any "
             "shared accession (shared = connected components).\n")
    L.append("## PSM discrepancies\n")
    ssd = metrics["psm_same_scan_different_peptide"]
    L.append(f"Same scan, different peptide: {ssd['ref_only']} ref-only PSMs / {ssd['test_only']} test-only PSMs.\n")
    L.append("| category | ref-only | test-only |")
    L.append("|---|---:|---:|")
    for c in PSM_CATEGORIES:
        L.append(f"| {c} | {metrics['psm_ref_only_breakdown'][c]} | {metrics['psm_test_only_breakdown'][c]} |")
    L.append("")
    sa = metrics["scan_agreement"]
    L.append(f"Scans identified on both sides: {sa['shared_scans']} (same peptide {sa['same_peptide']}, "
             f"I/L variant only {sa['IL_variant_only']}, different peptide {sa['different_peptide']}).\n")
    rp = metrics.get("ref_params") or {}
    if rp:
        L.append("## Reference filter (db.summary.html)\n")
        for k in ("peptide_min_10lgp", "protein_min_10lgp", "fdr_psm_pct", "fdr_peptide_pct",
                  "fdr_protein_group_pct", "enzyme", "digest_mode", "database"):
            if k in rp:
                L.append(f"- {k}: {rp[k]}")
        L.append("")
    if metrics["warnings"]:
        L.append("## Warnings\n")
        L.extend(f"- {w}" for w in metrics["warnings"])
        L.append("")
    return "\n".join(L)


def write_report(result, out_dir):
    """Write metrics.json, summary.md and one CSV per ref-only/test-only table."""
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "metrics.json"), "w") as fh:
        json.dump(result["metrics"], fh, indent=2, default=str)
    for name, df in result["tables"].items():
        df.to_csv(os.path.join(out_dir, f"{name}.csv"), index=False)
    with open(os.path.join(out_dir, "summary.md"), "w") as fh:
        fh.write(_markdown(result["metrics"]))
    return out_dir


def _parse_fixed(items):
    if not items:
        return None
    out = {}
    for it in items:
        aa, _, val = it.partition(":")
        if not aa or not val:
            raise argparse.ArgumentTypeError(f"--fixed-mod expects AA:DELTA or AA:Name, got {it!r}")
        out[aa] = val
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m proteoscan.concordance", description=__doc__.split("\n")[0])
    ap.add_argument("--ref", required=True, help="reference PEAKS export directory")
    ap.add_argument("--test", required=True, help="test export directory (PEAKS format)")
    ap.add_argument("--out", required=True, help="output directory for the report")
    ap.add_argument("--features", action="store_true", help="also compare features.csv")
    ap.add_argument("--feature-ppm", type=float, default=10.0)
    ap.add_argument("--feature-rt", type=float, default=0.5, help="RT tolerance in minutes")
    ap.add_argument("--ignore-source-file", action="store_true",
                    help="drop the source file from PSM/feature keys (single-run comparisons with renamed files)")
    ap.add_argument("--fixed-mod", action="append", metavar="AA:DELTA",
                    help="add a fixed mod to unmodified residues on both sides, e.g. C:+57.02")
    ap.add_argument("--merge-nterm", action="store_true",
                    help="fold N-terminal mods into the first residue before comparing")
    ap.add_argument("--min-10lgp", default=None,
                    help="only PSMs with -10LgP >= value on both sides ('ref' = reference peptide threshold)")
    a = ap.parse_args(argv)
    res = compare(a.ref, a.test, features=a.features, ignore_source=a.ignore_source_file,
                  fixed=_parse_fixed(a.fixed_mod), merge_nterm=a.merge_nterm, min_10lgp=a.min_10lgp,
                  feature_ppm=a.feature_ppm, feature_rt=a.feature_rt)
    write_report(res, a.out)
    h = res["metrics"]["headline"]
    print(f"PSM recall {_pct(h['psm_recall'])}  precision {_pct(h['psm_precision'])}  "
          f"({h['psm_shared']}/{h['psm_n_ref']} ref, {h['psm_n_test']} test) -> {a.out}")
    for w in res["metrics"]["warnings"]:
        print("warning:", w, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
