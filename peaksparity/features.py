"""Stage 2: MS1 feature detection (isotope envelopes, charge, monoisotopic m/z, RT range, area).

Public API
----------
detect_features(mzml_path, out_featurexml=None, **params) -> pandas.DataFrame
    Columns: feature_id, mz, charge, rt (apex, min), rt_begin, rt_end, area, quality,
    n_isotopes, ms1_scans, plus diagnostics (intensity = the OpenMS feature intensity,
    rt_model = the algorithm's own apex, area_mono, apex_intensity).
assign_ms2_to_features(features_df, mzml_path, isolation_window_da=None) -> pandas.DataFrame
    One row per (MS2 scan, co-isolated feature), sorted by feature area within each scan.
load_config(path) -> dict of keyword arguments for detect_features.

Algorithms (``algorithm=``)
    centroided  pyopenms FeatureFinderAlgorithmPicked (FeatureFinderCentroided).  Run in
                parallel over overlapping m/z slabs (``slabs``) because it is slow on
                large Orbitrap runs.
    metabo      MassTraceDetection -> ElutionPeakDetection -> FeatureFindingMetabo with
                the "peptides" isotope model (fast, trace based).
    multiplex   FeatureFinderMultiplexAlgorithm in label-free mode.

Whatever algorithm is used, every feature is re-quantified on the raw MS1 centroids
(``_reextract``): the isotope traces mono + k*1.003355/z (k < n_isotopes, at most
``area_isotopes``) are extracted at ``xic_ppm`` inside the feature's RT range and the
area is the trapezoid integral of the summed XIC with RT in minutes x ``area_scale``.
This makes ``area`` comparable between algorithms and with PEAKS "Area".

CLI::

    python -m peaksparity.features <mzML> --config configs/features_lumos_120k.yml \
        --out <dir> [--threads 16] [--isolation-window 1.6]

writes <dir>/<stem>.features.parquet, <stem>.ms2_features.parquet and <stem>.featureXML.
"""
from __future__ import annotations

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import multiprocessing as mp

import numpy as np
import pandas as pd
import yaml

# pyopenms uses OpenMP with all cores unless told otherwise; the server is shared.
os.environ.setdefault("OMP_NUM_THREADS", "16")

C13_DELTA = 1.0033548378
PROTON = 1.00727646688

__all__ = ["detect_features", "assign_ms2_to_features", "load_config", "load_ms1",
           "read_ms2_metadata", "DEFAULT_OPENMS", "C13_DELTA"]

# Parameters tuned for a Lumos 120k Orbitrap MS1 (centroided by ThermoRawFileParser).
DEFAULT_OPENMS: dict[str, dict[str, Any]] = {
    "centroided": {
        "intensity:bins": 10,
        "mass_trace:mz_tolerance": 0.006,
        "mass_trace:min_spectra": 4,
        "mass_trace:max_missing": 2,
        "mass_trace:slope_bound": 0.1,
        "isotopic_pattern:charge_low": 1,
        "isotopic_pattern:charge_high": 6,
        "isotopic_pattern:mz_tolerance": 0.008,
        "seed:min_score": 0.7,
        "feature:min_score": 0.6,
        "feature:min_isotope_fit": 0.7,
        "feature:min_trace_score": 0.4,
        "feature:min_rt_span": 0.2,
        "feature:max_rt_span": 4.0,
        "feature:max_intersection": 0.35,
    },
    "metabo": {
        # MassTraceDetection
        "mtd:mass_error_ppm": 8.0,
        "mtd:noise_threshold_int": 0.0,
        "mtd:chrom_peak_snr": 2.0,
        "mtd:min_sample_rate": 0.3,
        "mtd:min_trace_length": 3.0,
        "mtd:max_trace_length": 600.0,
        "mtd:trace_termination_criterion": "sample_rate",
        "mtd:trace_termination_outliers": 3,
        # ElutionPeakDetection
        "epd:enabled": "true",
        "epd:width_filtering": "fixed",
        "epd:min_fwhm": 1.0,
        "epd:max_fwhm": 90.0,
        "epd:chrom_fwhm": 8.0,
        "epd:chrom_peak_snr": 2.0,
        "epd:masstrace_snr_filtering": "false",
        # FeatureFindingMetabo
        "ffm:charge_lower_bound": 1,
        "ffm:charge_upper_bound": 6,
        "ffm:isotope_filtering_model": "peptides",
        "ffm:mz_scoring_13C": "true",
        "ffm:use_smoothed_intensities": "true",
        "ffm:report_convex_hulls": "true",
        "ffm:remove_single_traces": "true",
        "ffm:local_rt_range": 10.0,
        "ffm:local_mz_range": 6.5,
        "ffm:chrom_fwhm": 8.0,
        "ffm:enable_RT_filtering": "true",
        "ffm:mz_scoring_by_elements": "false",
        "ffm:report_chromatograms": "false",
    },
    "multiplex": {
        "algorithm:labels": "",
        "algorithm:charge": "1:6",
        "algorithm:isotopes_per_peptide": "3:6",
        "algorithm:rt_typical": 20.0,
        "algorithm:rt_band": 0.0,
        "algorithm:rt_min": 3.0,
        "algorithm:mz_tolerance": 8.0,
        "algorithm:mz_unit": "ppm",
        "algorithm:intensity_cutoff": 1000.0,
        "algorithm:peptide_similarity": 0.5,
        "algorithm:averagine_similarity": 0.4,
        "algorithm:averagine_similarity_scaling": 0.95,
        "algorithm:missed_cleavages": 0,
        "algorithm:spectrum_type": "centroid",
        "algorithm:averagine_type": "peptide",
        "algorithm:knock_out": "false",
    },
}


# --------------------------------------------------------------------------------------
# I/O helpers
# --------------------------------------------------------------------------------------

def load_config(path: str | os.PathLike) -> dict[str, Any]:
    """YAML config -> keyword arguments for :func:`detect_features`."""
    raw = yaml.safe_load(Path(path).read_text()) or {}
    return dict(raw)


def load_ms1(mzml_path: str | os.PathLike):
    """Load only the MS1 spectra of an mzML into a pyopenms MSExperiment (RT in seconds)."""
    import pyopenms as oms

    f = oms.MzMLFile()
    opt = f.getOptions()
    opt.setMSLevels([1])
    f.setOptions(opt)
    exp = oms.MSExperiment()
    f.load(str(mzml_path), exp)
    exp.sortSpectra(True)
    exp.updateRanges()
    return exp


def _ms1_arrays(exp) -> dict[str, np.ndarray]:
    """Flatten MS1 centroids: global key = scan_index * 4096 + m/z (sorted)."""
    rts, mzs, ints, offs = [], [], [], [0]
    for s in exp:
        mz, it = s.get_peaks()
        o = np.argsort(mz, kind="stable")
        mzs.append(mz[o].astype(np.float64))
        ints.append(it[o].astype(np.float64))
        rts.append(s.getRT() / 60.0)
        offs.append(offs[-1] + len(mz))
    scan = np.repeat(np.arange(len(rts)), np.diff(offs))
    mz = np.concatenate(mzs) if mzs else np.zeros(0)
    return {
        "rt": np.asarray(rts),
        "mz": mz,
        "int": np.concatenate(ints) if ints else np.zeros(0),
        "key": scan * 4096.0 + mz,
        "offs": np.asarray(offs),
    }


def read_ms2_metadata(mzml_path: str | os.PathLike) -> pd.DataFrame:
    """MS2 scan table: scan, rt (min), precursor m/z & charge, isolation window, MS1 parent."""
    import pyopenms as oms

    f = oms.MzMLFile()
    opt = f.getOptions()
    opt.setMSLevels([1, 2])
    opt.setFillData(False)
    f.setOptions(opt)
    exp = oms.MSExperiment()
    f.load(str(mzml_path), exp)
    rows = []
    last_ms1 = -1
    for s in exp:
        nid = s.getNativeID()
        nid = nid.decode() if isinstance(nid, bytes) else nid
        scan = int(nid.rsplit("scan=", 1)[-1]) if "scan=" in nid else -1
        if s.getMSLevel() == 1:
            last_ms1 = scan
            continue
        pcs = s.getPrecursors()
        if not pcs:
            continue
        p = pcs[0]
        rows.append((scan, s.getRT() / 60.0, p.getMZ(), p.getCharge(),
                     p.getIsolationWindowLowerOffset(), p.getIsolationWindowUpperOffset(),
                     p.getIntensity(), last_ms1))
    return pd.DataFrame(rows, columns=["scan", "rt", "precursor_mz", "precursor_charge",
                                       "iso_lower", "iso_upper", "precursor_intensity",
                                       "ms1_scan"])


# --------------------------------------------------------------------------------------
# Algorithms
# --------------------------------------------------------------------------------------

def _param(defaults, overrides: dict[str, Any]):
    p = defaults
    for k, v in overrides.items():
        kb = k.encode() if not isinstance(k, bytes) else k
        if not p.exists(kb):
            raise KeyError(f"unknown OpenMS parameter {k!r}")
        cur = p.getValue(kb)
        if isinstance(cur, float) and isinstance(v, int) and not isinstance(v, bool):
            v = float(v)
        if isinstance(v, bool):
            v = "true" if v else "false"
        p.setValue(kb, v)
    return p


def _ffc_slab(args):
    """Run FeatureFinderCentroided on one m/z slab (worker process)."""
    mzml_path, lo, hi, core_lo, core_hi, openms, omp = args
    os.environ["OMP_NUM_THREADS"] = str(omp)
    import pyopenms as oms

    exp = load_ms1(mzml_path)
    sub = oms.MSExperiment()
    for s in exp:
        mz, it = s.get_peaks()
        m = (mz >= lo) & (mz < hi)
        ns = oms.MSSpectrum(s)
        ns.set_peaks((mz[m], it[m]))
        sub.addSpectrum(ns)
    del exp
    sub.updateRanges()
    ff = oms.FeatureFinderAlgorithmPicked()
    p = _param(ff.getDefaults(), openms)
    out = oms.FeatureMap()
    ff.run(sub, out, p, oms.FeatureMap())
    recs = _featuremap_records(out)
    recs = [r for r in recs if core_lo <= r["mz"] < core_hi]
    return recs


def _featuremap_records(fm) -> list[dict[str, Any]]:
    recs = []
    for f in fm:
        hulls = f.getConvexHulls()
        n_iso = len(hulls)
        if n_iso:
            bbs = [h.getBoundingBox() for h in hulls]
            rt0 = min(b.minPosition()[0] for b in bbs)
            rt1 = max(b.maxPosition()[0] for b in bbs)
            mono = bbs[0]
            mrt0, mrt1 = mono.minPosition()[0], mono.maxPosition()[0]
        else:
            bb = f.getConvexHull().getBoundingBox()
            rt0, rt1 = bb.minPosition()[0], bb.maxPosition()[0]
            mrt0, mrt1 = rt0, rt1
        recs.append({
            "mz": f.getMZ(),
            "charge": f.getCharge(),
            "rt_model": f.getRT() / 60.0,
            "hull_rt_begin": rt0 / 60.0,
            "hull_rt_end": rt1 / 60.0,
            "mono_rt_begin": mrt0 / 60.0,
            "mono_rt_end": mrt1 / 60.0,
            "intensity": f.getIntensity(),
            "quality": f.getOverallQuality(),
            "n_isotopes": n_iso,
        })
    return recs


def _run_centroided(mzml_path, openms, threads, slabs, slab_overlap):
    """FFC in parallel over ``slabs`` overlapping m/z ranges (features kept by mono m/z)."""
    exp = load_ms1(mzml_path)
    lo = max(0.0, exp.getMinMZ() - 1)
    hi = exp.getMaxMZ() + 1
    del exp
    edges = np.linspace(lo, hi, slabs + 1)
    jobs = []
    for i in range(slabs):
        a, b = edges[i], edges[i + 1]
        # mono inside the core; its isotopes (and a -1 isotope check) inside the slab
        jobs.append((str(mzml_path), max(0.0, a - 1.5), b + slab_overlap, a, b, openms,
                     max(1, threads // slabs)))
    recs: list[dict[str, Any]] = []
    workers = max(1, min(slabs, threads))
    if workers == 1:
        for j in jobs:
            recs += _ffc_slab(j)
    else:
        # spawn (not fork): each worker must initialise OpenMP with its own thread count
        old = os.environ.get("OMP_NUM_THREADS")
        os.environ["OMP_NUM_THREADS"] = str(max(1, threads // workers))
        try:
            with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
                for r in ex.map(_ffc_slab, jobs):
                    recs += r
        finally:
            if old is not None:
                os.environ["OMP_NUM_THREADS"] = old
    return recs


def _run_metabo(exp, openms):
    import pyopenms as oms

    groups = {"mtd": {}, "epd": {}, "ffm": {}}
    for k, v in openms.items():
        g, _, name = k.partition(":")
        groups[g][name] = v
    mtd = oms.MassTraceDetection()
    mtd.setParameters(_param(mtd.getDefaults(), groups["mtd"]))
    traces = []
    t = time.time()
    mtd.run(exp, traces, 0)
    print(f"[features] metabo: {len(traces)} mass traces ({time.time() - t:.0f}s)", flush=True)
    epd_on = str(groups["epd"].pop("enabled", "true")).lower() == "true"
    if epd_on:
        epd = oms.ElutionPeakDetection()
        epd.setParameters(_param(epd.getDefaults(), groups["epd"]))
        split = []
        epd.detectPeaks(traces, split)
        if str(groups["epd"].get("width_filtering", "fixed")) != "off":
            filt = []
            epd.filterByPeakWidth(split, filt)
            split = filt
        traces = split
        print(f"[features] metabo: {len(traces)} elution peaks ({time.time() - t:.0f}s)", flush=True)
    ffm = oms.FeatureFindingMetabo()
    ffm.setParameters(_param(ffm.getDefaults(), groups["ffm"]))
    fm = oms.FeatureMap()
    chroms = []
    ffm.run(traces, fm, chroms)
    print(f"[features] metabo: {fm.size()} features ({time.time() - t:.0f}s)", flush=True)
    return fm


def _run_multiplex(exp, openms):
    import pyopenms as oms

    alg = oms.FeatureFinderMultiplexAlgorithm()
    alg.setParameters(_param(alg.getDefaults(), openms))
    alg.run(exp, False)
    return alg.getFeatureMap()


# --------------------------------------------------------------------------------------
# Re-quantification on raw MS1 centroids
# --------------------------------------------------------------------------------------

def _lookup(ms1, scan_idx, target_mz, ppm):
    """Intensity (and m/z) of the nearest centroid within ppm of target in each scan."""
    q = scan_idx * 4096.0 + target_mz
    pos = np.searchsorted(ms1["key"], q)
    n = len(ms1["key"])
    best_int = np.zeros(len(q))
    best_mz = np.full(len(q), np.nan)
    best_d = np.full(len(q), np.inf)
    for cand in (pos - 1, pos):
        ok = (cand >= 0) & (cand < n)
        c = np.where(ok, cand, 0)
        d = np.abs(ms1["key"][c] - q)
        ok &= d <= target_mz * ppm * 1e-6
        better = ok & (d < best_d)
        best_d = np.where(better, d, best_d)
        best_int = np.where(better, ms1["int"][c], best_int)
        best_mz = np.where(better, ms1["mz"][c], best_mz)
    return best_int, best_mz


def _reextract(df: pd.DataFrame, ms1, ppm: float, max_iso: int, area_scale: float,
               apex_from: str, boundary_frac: float) -> pd.DataFrame:
    """Re-extract isotope XICs inside each feature's RT range; area, apex, scan count."""
    rt = ms1["rt"]
    s0 = np.searchsorted(rt, df["rt_begin"].to_numpy() - 1e-9, side="left")
    s1 = np.searchsorted(rt, df["rt_end"].to_numpy() + 1e-9, side="right") - 1
    s1 = np.maximum(s1, s0)
    s1 = np.minimum(s1, len(rt) - 1)
    s0 = np.minimum(s0, s1)
    n = (s1 - s0 + 1).astype(np.int64)
    fid = np.repeat(np.arange(len(df)), n)
    start = np.repeat(np.cumsum(n) - n, n)
    scan = np.repeat(s0, n) + (np.arange(n.sum()) - start)
    mono = df["mz"].to_numpy()[fid]
    z = np.maximum(df["charge"].to_numpy()[fid], 1)
    niso = np.clip(df["n_isotopes"].to_numpy(), 1, max_iso)[fid]
    xic = np.zeros(len(fid))
    xic0 = None
    for k in range(max_iso):
        inten, _ = _lookup(ms1, scan, mono + k * C13_DELTA / z, ppm)
        inten = np.where(k < niso, inten, 0.0)
        xic += inten
        if k == 0:
            xic0 = inten
    t = rt[scan]
    # optional boundary trimming at a fraction of the apex (per feature)
    starts = np.cumsum(n) - n
    apex_val = np.maximum.reduceat(xic, starts) if len(xic) else np.zeros(0)
    keep = np.ones(len(xic), bool)
    if boundary_frac > 0 and len(xic):
        keep = xic >= boundary_frac * apex_val[fid]
    # trapezoid per feature over consecutive kept scans
    same = np.r_[fid[1:] == fid[:-1], False]
    dt = np.r_[t[1:] - t[:-1], 0.0]
    both = same & keep & np.r_[keep[1:], False]
    seg = np.where(both, 0.5 * (xic + np.r_[xic[1:], 0.0]) * dt, 0.0)
    seg0 = np.where(both, 0.5 * (xic0 + np.r_[xic0[1:], 0.0]) * dt, 0.0)
    area = np.add.reduceat(seg, starts) * area_scale
    area0 = np.add.reduceat(seg0, starts) * area_scale
    ms1_scans = np.add.reduceat(((xic0 > 0) & keep).astype(int), starts)
    # apex = scan of max summed XIC
    order = np.lexsort((-xic, fid))
    first = np.r_[True, fid[order][1:] != fid[order][:-1]]
    apex_rt = t[order][first]
    # kept range
    tk = np.where(keep & (xic > 0), t, np.nan)
    with np.errstate(all="ignore"):
        rb = pd.Series(tk).groupby(fid).min().reindex(range(len(df))).to_numpy()
        re_ = pd.Series(tk).groupby(fid).max().reindex(range(len(df))).to_numpy()
    out = df.copy()
    out["area"] = area
    out["area_mono"] = area0
    out["apex_intensity"] = apex_val
    out["ms1_scans"] = ms1_scans
    out["rt_xic_apex"] = apex_rt
    if boundary_frac > 0:
        out["rt_begin"] = np.where(np.isnan(rb), out["rt_begin"], rb)
        out["rt_end"] = np.where(np.isnan(re_), out["rt_end"], re_)
    out["rt"] = out["rt_xic_apex"] if apex_from == "xic" else out["rt_model"]
    return out


# --------------------------------------------------------------------------------------
# Main entry points
# --------------------------------------------------------------------------------------

def detect_features(mzml_path: str | os.PathLike, out_featurexml: str | os.PathLike | None = None,
                    algorithm: str = "centroided", openms: dict[str, Any] | None = None,
                    threads: int = 16, slabs: int = 16, slab_overlap: float = 6.0,
                    rt_range: str = "hull", xic_ppm: float = 10.0, area_isotopes: int = 3,
                    area_scale: float = 1.0, apex_from: str = "model",
                    boundary_frac: float = 0.0, min_charge: int = 1, max_charge: int = 9,
                    min_ms1_scans: int = 1, verbose: bool = True, **params) -> pd.DataFrame:
    """Detect MS1 features in an mzML and return one row per feature.

    ``openms`` (and any extra keyword ``section__name=value``, mapped to
    ``section:name``) override :data:`DEFAULT_OPENMS` for the chosen algorithm.
    ``rt_range`` = hull (union of isotope hulls) | mono (monoisotopic trace hull).
    ``apex_from`` = model (the algorithm's RT) | xic (scan of max re-extracted intensity).
    ``boundary_frac`` > 0 trims rt_begin/rt_end to scans whose summed XIC is at least this
    fraction of the apex.  RT columns are in minutes.
    """
    t0 = time.time()
    algorithm = algorithm.lower()
    if algorithm not in DEFAULT_OPENMS:
        raise ValueError(f"unknown algorithm {algorithm!r}")
    op = dict(DEFAULT_OPENMS[algorithm])
    op.update(openms or {})
    op.update({k.replace("__", ":"): v for k, v in params.items()})
    os.environ["OMP_NUM_THREADS"] = str(threads)  # effective only if pyopenms not yet loaded

    fm = None
    if algorithm == "centroided":
        recs = _run_centroided(mzml_path, op, threads, slabs, slab_overlap)
        exp = load_ms1(mzml_path)
    else:
        exp = load_ms1(mzml_path)
        fm = _run_metabo(exp, op) if algorithm == "metabo" else _run_multiplex(exp, op)
        recs = _featuremap_records(fm)
    t_detect = time.time() - t0
    df = pd.DataFrame(recs, columns=["mz", "charge", "rt_model", "hull_rt_begin", "hull_rt_end",
                                     "mono_rt_begin", "mono_rt_end", "intensity", "quality",
                                     "n_isotopes"])
    df = df[(df["charge"] >= min_charge) & (df["charge"] <= max_charge)].copy()
    if rt_range == "mono":
        df["rt_begin"], df["rt_end"] = df["mono_rt_begin"], df["mono_rt_end"]
    else:
        df["rt_begin"], df["rt_end"] = df["hull_rt_begin"], df["hull_rt_end"]
    df = df.sort_values(["rt_model", "mz"]).reset_index(drop=True)
    ms1 = _ms1_arrays(exp)
    df = _reextract(df, ms1, xic_ppm, area_isotopes, area_scale, apex_from, boundary_frac)
    df = df[df["ms1_scans"] >= min_ms1_scans].reset_index(drop=True)
    df.insert(0, "feature_id", np.arange(1, len(df) + 1))
    cols = ["feature_id", "mz", "charge", "rt", "rt_begin", "rt_end", "area", "quality",
            "n_isotopes", "ms1_scans"]
    df = df[cols + [c for c in df.columns if c not in cols]]
    df.attrs["runtime_s"] = time.time() - t0
    df.attrs["detect_s"] = t_detect
    if verbose:
        print(f"[features] {algorithm}: {len(df)} features in {time.time() - t0:.0f}s "
              f"(detection {t_detect:.0f}s)", flush=True)
    if out_featurexml is not None:
        _write_featurexml(df, out_featurexml)
    return df


def _write_featurexml(df: pd.DataFrame, path) -> None:
    """Store the final feature table as featureXML (RT in seconds, id = feature_id)."""
    import pyopenms as oms

    fm = oms.FeatureMap()
    for r in df.itertuples(index=False):
        f = oms.Feature()
        f.setMZ(float(r.mz))
        f.setRT(float(r.rt) * 60.0)
        f.setCharge(int(r.charge))
        f.setIntensity(float(r.area))
        f.setOverallQuality(float(r.quality))
        f.setUniqueId(int(r.feature_id))
        h = oms.ConvexHull2D()
        h.setHullPoints(np.array([[r.rt_begin * 60.0, r.mz],
                                  [r.rt_end * 60.0, r.mz],
                                  [r.rt_end * 60.0, r.mz + (max(r.n_isotopes, 1) - 1) * C13_DELTA / r.charge],
                                  [r.rt_begin * 60.0, r.mz + (max(r.n_isotopes, 1) - 1) * C13_DELTA / r.charge]]))
        f.setConvexHulls([h])
        f.setMetaValue("n_isotopes", int(r.n_isotopes))
        f.setMetaValue("ms1_scans", int(r.ms1_scans))
        fm.push_back(f)
    fm.setUniqueIds()
    oms.FeatureXMLFile().store(str(path), fm)


def assign_ms2_to_features(features_df: pd.DataFrame, mzml_path: str | os.PathLike | pd.DataFrame,
                           isolation_window_da: float | None = None, n_isotopes: int = 4,
                           rt_tolerance: float = 0.0) -> pd.DataFrame:
    """List, for every MS2 scan, the features co-isolated with it.

    A feature is co-isolated when rt_begin - tol <= scan RT <= rt_end + tol and any of its
    isotope peaks mono + k*1.003355/z (k = 0..n_isotopes-1) lies inside the isolation
    window [target - w/2, target + w/2].  ``isolation_window_da`` is the full width; None
    uses the offsets written in the mzML.  ``mzml_path`` may also be a DataFrame from
    :func:`read_ms2_metadata`.  Rows are sorted by scan, then feature area (descending);
    ``rank`` is 1 for the most intense feature.  MS2 scans with no feature get one row
    with feature_id = -1.
    """
    ms2 = mzml_path if isinstance(mzml_path, pd.DataFrame) else read_ms2_metadata(mzml_path)
    if isolation_window_da is None:
        lo_off = ms2["iso_lower"].to_numpy()
        hi_off = ms2["iso_upper"].to_numpy()
    else:
        lo_off = hi_off = np.full(len(ms2), isolation_window_da / 2.0)
    target = ms2["precursor_mz"].to_numpy()
    lo, hi = target - lo_off, target + hi_off
    srt = ms2["rt"].to_numpy()

    f = features_df.reset_index(drop=True)
    fz = f["charge"].to_numpy().clip(min=1)
    k = np.arange(n_isotopes)
    peak_mz = (f["mz"].to_numpy()[:, None] + k[None, :] * C13_DELTA / fz[:, None]).ravel()
    peak_f = np.repeat(np.arange(len(f)), n_isotopes)
    o = np.argsort(peak_mz)
    peak_mz, peak_f = peak_mz[o], peak_f[o]
    a = np.searchsorted(peak_mz, lo, side="left")
    b = np.searchsorted(peak_mz, hi, side="right")
    cnt = b - a
    scan_i = np.repeat(np.arange(len(ms2)), cnt)
    st = np.repeat(np.cumsum(cnt) - cnt, cnt)
    idx = np.repeat(a, cnt) + (np.arange(cnt.sum()) - st)
    fi = peak_f[idx]
    rb = f["rt_begin"].to_numpy()[fi] - rt_tolerance
    re_ = f["rt_end"].to_numpy()[fi] + rt_tolerance
    ok = (srt[scan_i] >= rb) & (srt[scan_i] <= re_)
    pairs = pd.DataFrame({"i": scan_i[ok], "f": fi[ok]})
    pairs = pairs.groupby(["i", "f"]).size().rename("n_isotopes_in_window").reset_index()
    out = pd.DataFrame({
        "scan": ms2["scan"].to_numpy()[pairs["i"]],
        "ms2_rt": srt[pairs["i"]],
        "isolation_target": target[pairs["i"]],
        "precursor_charge": ms2["precursor_charge"].to_numpy()[pairs["i"]],
        "feature_id": f["feature_id"].to_numpy()[pairs["f"]],
        "feature_mz": f["mz"].to_numpy()[pairs["f"]],
        "feature_charge": f["charge"].to_numpy()[pairs["f"]],
        "feature_rt": f["rt"].to_numpy()[pairs["f"]],
        "feature_area": f["area"].to_numpy()[pairs["f"]],
        "n_isotopes_in_window": pairs["n_isotopes_in_window"].to_numpy(),
    })
    missing = np.setdiff1d(np.arange(len(ms2)), pairs["i"].unique())
    if len(missing):
        empty = pd.DataFrame({
            "scan": ms2["scan"].to_numpy()[missing], "ms2_rt": srt[missing],
            "isolation_target": target[missing],
            "precursor_charge": ms2["precursor_charge"].to_numpy()[missing],
            "feature_id": -1, "feature_mz": np.nan, "feature_charge": 0,
            "feature_rt": np.nan, "feature_area": np.nan, "n_isotopes_in_window": 0})
        out = pd.concat([out, empty], ignore_index=True)
    out = out.sort_values(["scan", "feature_area"], ascending=[True, False],
                          na_position="last").reset_index(drop=True)
    out["rank"] = out.groupby("scan").cumcount() + 1
    out.loc[out["feature_id"] < 0, "rank"] = 0
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mzml")
    ap.add_argument("--config", help="YAML config (configs/features_lumos_120k.yml)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--isolation-window", type=float, default=None,
                    help="full isolation width in Th (default: from the mzML)")
    a = ap.parse_args(argv)
    kw = load_config(a.config) if a.config else {}
    if a.threads:
        kw["threads"] = a.threads
    iso = kw.pop("isolation_window_da", None) if a.isolation_window is None else a.isolation_window
    kw.pop("isolation_window_da", None)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(a.mzml).name
    for ext in (".gz", ".mzML", ".mzml"):
        stem = stem[: -len(ext)] if stem.endswith(ext) else stem
    df = detect_features(a.mzml, out / f"{stem}.featureXML", **kw)
    df.to_parquet(out / f"{stem}.features.parquet", index=False)
    ms2 = assign_ms2_to_features(df, a.mzml, iso)
    ms2.to_parquet(out / f"{stem}.ms2_features.parquet", index=False)
    print(f"[features] wrote {out}/{stem}.features.parquet ({len(df)}) and "
          f".ms2_features.parquet ({len(ms2)} rows, "
          f"{ms2.loc[ms2.feature_id >= 0, 'scan'].nunique()} scans with a feature)")


if __name__ == "__main__":
    main()
