#!/usr/bin/env python
"""Compare detected MS1 features (proteoscan.features parquet) with PEAKS's features.csv.

Usage::

    python scripts/eval_features.py --features <run>/<stem>.features.parquet \
        --ref "<PEAKS export>/Sample 4.db.features.csv" --mzml <run.mzML> \
        [--ms2 <run>/<stem>.ms2_features.parquet] [--mgf "<PEAKS export>/<run>.mgf"] \
        [--out <dir>] [--label NAME] [--ppm 10] [--rt 0.3]

Matching: same charge, monoisotopic m/z within --ppm, and RT ranges overlap or apexes
within --rt minutes ("recall_strict": overlap AND apexes within --rt).  Each PEAKS feature takes its best candidate (smallest
|ppm|/ppm_tol + |dRT|/rt_tol); several PEAKS features may share one detected feature.

Reported
  recall                PEAKS features with a match / all PEAKS features
  agreement             Spearman of log area, median signed/abs ppm, median apex dRT,
                        median log10 area ratio (scale), begin/end differences
  miss reasons          first applicable of: charge_differs, isotope_off_ours_lower /
                        isotope_off_ours_higher (same z, our mono 1-2 isotopes below/above), in_envelope (PEAKS m/z is a non-mono peak of a
                        detected feature of any charge), rt_mismatch (same z and m/z,
                        apex within 2 min), mz_10_25ppm, no_feature_signal (no feature,
                        but mono in >=3 and M+1 in >=2 MS1 scans of the PEAKS range),
                        no_signal (no feature and hardly any MS1 centroids there)
  precision proxy       detected features that were fragmented (>=1 MS2 scan whose
                        isolation window holds an isotope of the feature while it elutes),
                        and how many of them are not in features.csv / not in the PEAKS
                        MGF precursor list
  MS2 queries (--mgf)   PEAKS writes one MGF entry per (MS2 scan, co-isolated feature);
                        entries with PEPMASS = feature m/z and RTINSECONDS = feature apex.
                        We compare those (scan, m/z, z) queries with ours.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan.features import (C13_DELTA, _lookup, _ms1_arrays, assign_ms2_to_features,  # noqa: E402
                                  load_ms1, read_ms2_metadata)


def read_ref(path) -> pd.DataFrame:
    r = pd.read_csv(path)
    return pd.DataFrame({
        "ref_id": r["Feature Id"].to_numpy(), "mz": r["m/z"].to_numpy(float),
        "charge": r["z"].to_numpy(int), "rt": r["RT"].to_numpy(float),
        "rt_begin": r["RT Begin"].to_numpy(float), "rt_end": r["RT End"].to_numpy(float),
        "area": r["Area"].to_numpy(float),
        "identified": (r["DB Peptide"].notna() | r["Denovo Peptide"].notna()).to_numpy(),
    })


def read_mgf_headers(path, cache: Path | None = None) -> pd.DataFrame:
    """(query, raw_scan, mz, z, rt) from a PEAKS MGF export (header lines only)."""
    if cache is not None and cache.exists():
        return pd.read_parquet(cache)
    rows, cur = [], {}
    with open(path) as fh:
        for line in fh:
            c = line[:1]
            if c.isdigit() or c == "\n":
                continue
            line = line.strip()
            if line == "BEGIN IONS":
                cur = {}
            elif line == "END IONS":
                rows.append(cur)
            elif "=" in line:
                k, v = line.split("=", 1)
                cur[k] = v
    d = pd.DataFrame(rows)
    out = pd.DataFrame({
        "query": d["SCANS"].astype(int), "raw_scan": d["RAWSCANS"].astype(int),
        "mz": d["PEPMASS"].str.split().str[0].astype(float),
        "z": d["CHARGE"].str.rstrip("+").astype(int),
        "rt": d["RTINSECONDS"].astype(float) / 60.0,
    })
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        out.to_parquet(cache, index=False)
    return out


def pairs_within(a_mz: np.ndarray, b_mz: np.ndarray, ppm: float):
    """All (i, j) with |a_i - b_j| <= ppm * a_i; returns index arrays and signed ppm."""
    o = np.argsort(b_mz)
    bs = b_mz[o]
    lo = np.searchsorted(bs, a_mz * (1 - ppm * 1e-6), side="left")
    hi = np.searchsorted(bs, a_mz * (1 + ppm * 1e-6), side="right")
    n = hi - lo
    i = np.repeat(np.arange(len(a_mz)), n)
    st = np.repeat(np.cumsum(n) - n, n)
    j = o[np.repeat(lo, n) + (np.arange(n.sum()) - st)]
    return i, j, (b_mz[j] - a_mz[i]) / a_mz[i] * 1e6


def rt_ok(ref, det, i, j, rt_tol):
    overlap = (det["rt_begin"].to_numpy()[j] <= ref["rt_end"].to_numpy()[i]) & \
              (det["rt_end"].to_numpy()[j] >= ref["rt_begin"].to_numpy()[i])
    d = np.abs(det["rt"].to_numpy()[j] - ref["rt"].to_numpy()[i])
    return overlap | (d <= rt_tol), d


def match(ref: pd.DataFrame, det: pd.DataFrame, ppm: float, rt_tol: float,
          strict: bool = False) -> pd.DataFrame:
    """Best detected feature per reference feature.  strict: apex within rt_tol AND overlap."""
    i, j, dppm = pairs_within(ref["mz"].to_numpy(), det["mz"].to_numpy(), ppm)
    ok, drt = rt_ok(ref, det, i, j, rt_tol)
    if strict:
        overlap = (det["rt_begin"].to_numpy()[j] <= ref["rt_end"].to_numpy()[i]) & \
                  (det["rt_end"].to_numpy()[j] >= ref["rt_begin"].to_numpy()[i])
        ok = overlap & (drt <= rt_tol)
    ok &= ref["charge"].to_numpy()[i] == det["charge"].to_numpy()[j]
    m = pd.DataFrame({"i": i[ok], "j": j[ok], "ppm": dppm[ok], "drt": drt[ok]})
    m["score"] = m["ppm"].abs() / ppm + m["drt"] / rt_tol
    m = m.sort_values("score").drop_duplicates("i")
    return m


def miss_reasons(ref, det, missed_idx, ppm, rt_tol, ms1) -> pd.Series:
    """Classify unmatched PEAKS features (first applicable reason)."""
    r = ref.iloc[missed_idx].reset_index(drop=True)
    reason = np.full(len(r), "", dtype=object)
    rz = r["charge"].to_numpy()
    dz = det["charge"].to_numpy()

    def with_rt(i, j):
        ok, _ = rt_ok(r, det, i, j, rt_tol)
        return i[ok], j[ok]

    # charge differs: same m/z (10 ppm), RT ok, different charge
    i, j, _ = pairs_within(r["mz"].to_numpy(), det["mz"].to_numpy(), ppm)
    i, j = with_rt(i, j)
    cd = np.zeros(len(r), bool)
    cd[i[dz[j] != rz[i]]] = True
    # isotope off: same charge, mono shifted by +-1, +-2 isotopes
    iso = np.zeros(len(r), bool)
    iso_k = np.zeros(len(r), int)
    for k in (-2, -1, 1, 2):
        i, j, _ = pairs_within(r["mz"].to_numpy() + k * C13_DELTA / rz, det["mz"].to_numpy(), ppm)
        i, j = with_rt(i, j)
        sel = i[dz[j] == rz[i]]
        iso_k[sel[~iso[sel]]] = k
        iso[sel] = True
    # in envelope: PEAKS m/z equals a non-mono isotope peak of a detected feature (any z)
    env = np.zeros(len(r), bool)
    for k in (1, 2, 3):
        i, j, _ = pairs_within(r["mz"].to_numpy(), det["mz"].to_numpy() + k * C13_DELTA / dz, ppm)
        i, j = with_rt(i, j)
        env[i] = True
    # RT mismatch: same z, m/z ok, apex within 2 min but no overlap / apex > rt_tol
    i, j, _ = pairs_within(r["mz"].to_numpy(), det["mz"].to_numpy(), ppm)
    same = (dz[j] == rz[i]) & (np.abs(det["rt"].to_numpy()[j] - r["rt"].to_numpy()[i]) <= 2.0)
    rtm = np.zeros(len(r), bool)
    rtm[i[same]] = True
    # m/z 10-25 ppm
    i, j, _ = pairs_within(r["mz"].to_numpy(), det["mz"].to_numpy(), 25.0)
    i, j = with_rt(i, j)
    mzo = np.zeros(len(r), bool)
    mzo[i[dz[j] == rz[i]]] = True
    # MS1 signal at the PEAKS m/z within the PEAKS RT range
    rt = ms1["rt"]
    s0 = np.searchsorted(rt, r["rt_begin"].to_numpy(), side="left")
    s1 = np.maximum(np.searchsorted(rt, r["rt_end"].to_numpy(), side="right") - 1, s0)
    n = s1 - s0 + 1
    fid = np.repeat(np.arange(len(r)), n)
    st = np.repeat(np.cumsum(n) - n, n)
    scan = np.minimum(np.repeat(s0, n) + (np.arange(n.sum()) - st), len(rt) - 1)
    inten, _ = _lookup(ms1, scan, r["mz"].to_numpy()[fid], ppm)
    inten1, _ = _lookup(ms1, scan, (r["mz"].to_numpy() + C13_DELTA / rz)[fid], ppm)
    nsc = np.bincount(fid, weights=(inten > 0), minlength=len(r))
    nsc1 = np.bincount(fid, weights=(inten1 > 0), minlength=len(r))
    sig = (nsc >= 3) & (nsc1 >= 2)
    for name, mask in [("charge_differs", cd), ("isotope_off_ours_lower", iso & (iso_k < 0)),
                       ("isotope_off_ours_higher", iso & (iso_k > 0)), ("in_envelope", env),
                       ("rt_mismatch", rtm), ("mz_10_25ppm", mzo),
                       ("no_feature_signal", sig), ("no_signal", np.ones(len(r), bool))]:
        reason = np.where((reason == "") & mask, name, reason)
    return pd.Series(reason, index=missed_idx), pd.Series(nsc, index=missed_idx)


def evaluate(det: pd.DataFrame, ref: pd.DataFrame, ms1, ms2=None, mgf=None, ms2meta=None,
             ppm=10.0, rt_tol=0.3) -> tuple[dict, pd.DataFrame]:
    res: dict = {"n_detected": len(det), "n_ref": len(ref)}
    m = match(ref, det, ppm, rt_tol)
    res["recall"] = len(m) / len(ref)
    ms = match(ref, det, ppm, rt_tol, strict=True)
    res["recall_strict"] = len(ms) / len(ref)
    res["n_matched_ref"] = len(m)
    res["n_matched_det"] = int(m["j"].nunique())
    ident = ref["identified"].to_numpy()
    res["recall_identified"] = float(np.isin(np.where(ident)[0], m["i"]).mean()) if ident.any() else None
    for z in (1, 2, 3, 4):
        sel = np.where(ref["charge"].to_numpy() == z)[0]
        res[f"recall_z{z}"] = float(np.isin(sel, m["i"]).mean()) if len(sel) else None
    sel = np.where(ref["charge"].to_numpy() >= 5)[0]
    res["recall_z5+"] = float(np.isin(sel, m["i"]).mean()) if len(sel) else None
    ra = ref["area"].to_numpy()[m["i"]]
    da = det["area"].to_numpy()[m["j"]]
    good = (ra > 0) & (da > 0)
    res["spearman_log_area"] = float(spearmanr(np.log(ra[good]), np.log(da[good]))[0])
    res["median_log10_area_ratio"] = float(np.median(np.log10(da[good] / ra[good])))
    res["iqr_log10_area_ratio"] = float(np.subtract(*np.percentile(np.log10(da[good] / ra[good]), [75, 25])))
    res["median_ppm"] = float(np.median(m["ppm"]))
    res["median_abs_ppm"] = float(np.median(np.abs(m["ppm"])))
    rrt = ref["rt"].to_numpy()[m["i"]]
    drt_s = det["rt"].to_numpy()[m["j"]] - rrt
    res["median_rt_apex_diff_min"] = float(np.median(drt_s))
    res["median_abs_rt_apex_diff_min"] = float(np.median(np.abs(drt_s)))
    res["median_rt_begin_diff_min"] = float(np.median(det["rt_begin"].to_numpy()[m["j"]] - ref["rt_begin"].to_numpy()[m["i"]]))
    res["median_rt_end_diff_min"] = float(np.median(det["rt_end"].to_numpy()[m["j"]] - ref["rt_end"].to_numpy()[m["i"]]))
    rw = (ref["rt_end"] - ref["rt_begin"]).to_numpy()[m["i"]]
    dw = (det["rt_end"] - det["rt_begin"]).to_numpy()[m["j"]]
    res["median_width_ratio"] = float(np.median((dw + 1e-3) / (rw + 1e-3)))
    # misses
    missed = np.setdiff1d(np.arange(len(ref)), m["i"].to_numpy())
    reasons, nsc = miss_reasons(ref, det, missed, ppm, rt_tol, ms1)
    res["miss_reasons"] = reasons.value_counts().to_dict()
    # per-ref table
    tab = ref.copy()
    tab["det_feature_id"] = -1
    tab.loc[m["i"].to_numpy(), "det_feature_id"] = det["feature_id"].to_numpy()[m["j"]]
    tab["ppm"] = np.nan
    tab.loc[m["i"].to_numpy(), "ppm"] = m["ppm"].to_numpy()
    tab["det_rt"] = np.nan
    tab.loc[m["i"].to_numpy(), "det_rt"] = det["rt"].to_numpy()[m["j"]]
    tab["det_area"] = np.nan
    tab.loc[m["i"].to_numpy(), "det_area"] = da
    tab["miss_reason"] = ""
    tab.loc[missed, "miss_reason"] = reasons.to_numpy()
    tab["ms1_scans_at_ref_mz"] = np.nan
    tab.loc[missed, "ms1_scans_at_ref_mz"] = nsc.to_numpy()

    # precision proxy: fragmented detected features
    if ms2 is not None:
        frag = ms2.loc[ms2["feature_id"] >= 0, "feature_id"].unique()
        is_frag = det["feature_id"].isin(frag).to_numpy()
        matched_det = np.zeros(len(det), bool)
        matched_det[m["j"].to_numpy()] = True
        res["n_det_fragmented"] = int(is_frag.sum())
        res["n_det_fragmented_not_in_ref"] = int((is_frag & ~matched_det).sum())
        res["frac_fragmented_in_ref"] = float((is_frag & matched_det).sum() / max(is_frag.sum(), 1))
        res["n_det_top1_fragmented"] = int(ms2.loc[ms2["rank"] == 1, "feature_id"].nunique())
    if mgf is not None:
        res.update(mgf_compare(det, mgf, ms2, ms2meta, ppm, rt_tol))
    return res, tab


def mgf_compare(det, mgf, ms2, ms2meta, ppm, rt_tol) -> dict:
    """Compare against PEAKS's per-(scan, feature) MGF queries."""
    out = {}
    q = mgf.copy()
    # a feature apex is an MS1 scan time, so an entry whose RT equals its MS2 scan's own
    # RT is a PEAKS fallback query (no feature; m/z may still be mono-corrected)
    if ms2meta is not None:
        mm = ms2meta.set_index("scan")
        srt = mm["rt"].reindex(q["raw_scan"]).to_numpy()
        q["fallback"] = np.abs(q["rt"] - srt) < 2e-4
    else:
        q["fallback"] = False
    out["mgf_queries"] = len(q)
    out["mgf_scans"] = int(q["raw_scan"].nunique())
    out["mgf_fallback_queries"] = int(q["fallback"].sum())
    # unique PEAKS features with MS2 (m/z, z, apex rt)
    fq = q[~q["fallback"]].drop_duplicates(["mz", "z", "rt"]).reset_index(drop=True)
    out["mgf_unique_features"] = len(fq)
    i, j, dppm = pairs_within(fq["mz"].to_numpy(), det["mz"].to_numpy(), ppm)
    ok = fq["z"].to_numpy()[i] == det["charge"].to_numpy()[j]
    rr = fq["rt"].to_numpy()[i]
    ok &= ((det["rt_begin"].to_numpy()[j] <= rr) & (det["rt_end"].to_numpy()[j] >= rr)) | \
          (np.abs(det["rt"].to_numpy()[j] - rr) <= rt_tol)
    hit = np.zeros(len(fq), bool)
    hit[i[ok]] = True
    out["mgf_feature_recall"] = float(hit.mean())
    matched_det = np.zeros(len(det), bool)
    matched_det[j[ok]] = True
    if ms2 is not None:
        frag = det["feature_id"].isin(ms2.loc[ms2["feature_id"] >= 0, "feature_id"].unique()).to_numpy()
        out["frac_fragmented_in_mgf"] = float((frag & matched_det).sum() / max(frag.sum(), 1))
        out["n_det_fragmented_not_in_mgf"] = int((frag & ~matched_det).sum())
        # per-query comparison (scan, m/z, z)
        a = ms2[ms2["feature_id"] >= 0]
        out["our_queries"] = len(a)
        out["our_scans_with_feature"] = int(a["scan"].nunique())
        qq = q[~q["fallback"]].reset_index(drop=True)
        key_q = qq["raw_scan"].to_numpy()
        i, j, _ = pairs_within(qq["mz"].to_numpy(), a["feature_mz"].to_numpy(), ppm)
        ok = (key_q[i] == a["scan"].to_numpy()[j]) & (qq["z"].to_numpy()[i] == a["feature_charge"].to_numpy()[j])
        out["query_recall"] = float(np.isin(np.arange(len(qq)), i[ok]).mean())
        out["query_precision"] = float(np.isin(np.arange(len(a)), j[ok]).mean())
        # top-ranked feature per scan vs PEAKS's first query of the scan
        first = qq.sort_values("query").drop_duplicates("raw_scan")
        top = a[a["rank"] == 1].set_index("scan")
        tmz = top["feature_mz"].reindex(first["raw_scan"]).to_numpy()
        tz = top["feature_charge"].reindex(first["raw_scan"]).to_numpy()
        agree = (np.abs(tmz - first["mz"].to_numpy()) <= first["mz"].to_numpy() * ppm * 1e-6) & \
                (tz == first["z"].to_numpy())
        out["first_query_top1_agreement"] = float(np.nanmean(agree))
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--features", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--mzml", required=True)
    ap.add_argument("--ms2", help="ms2_features.parquet (default: next to --features)")
    ap.add_argument("--mgf", help="PEAKS MGF export of the same run")
    ap.add_argument("--cache", default="/drive_a/drew/proteoscan-work/analysis/features/cache")
    ap.add_argument("--out")
    ap.add_argument("--label", default="")
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--rt", type=float, default=0.3)
    a = ap.parse_args(argv)
    cache = Path(a.cache)
    cache.mkdir(parents=True, exist_ok=True)
    stem = Path(a.mzml).stem
    det = pd.read_parquet(a.features)
    ref = read_ref(a.ref)
    ms1 = _ms1_arrays(load_ms1(a.mzml))
    meta_p = cache / f"{stem}.ms2meta.parquet"
    if meta_p.exists():
        ms2meta = pd.read_parquet(meta_p)
    else:
        ms2meta = read_ms2_metadata(a.mzml)
        ms2meta.to_parquet(meta_p, index=False)
    ms2_path = a.ms2 or str(a.features).replace(".features.parquet", ".ms2_features.parquet")
    ms2 = pd.read_parquet(ms2_path) if Path(ms2_path).exists() else \
        assign_ms2_to_features(det, ms2meta, rt_tolerance=0.02, min_window_fraction=0.2)
    mgf = read_mgf_headers(a.mgf, cache / f"{stem}.mgf_headers.parquet") if a.mgf else None
    res, tab = evaluate(det, ref, ms1, ms2, mgf, ms2meta, a.ppm, a.rt)
    res = {"label": a.label, "features": str(a.features), **res}
    print(json.dumps(res, indent=1, default=float))
    if a.out:
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "eval.json").write_text(json.dumps(res, indent=1, default=float))
        tab.to_parquet(out / "ref_matches.parquet", index=False)


if __name__ == "__main__":
    main()
