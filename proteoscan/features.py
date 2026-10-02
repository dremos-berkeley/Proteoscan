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
    traces      proteoscan.feature_traces: numba mass-trace linking, valley splitting and
                averagine isotope grouping (no pyopenms algorithm; fast).

Whatever algorithm is used, every feature is re-quantified on the raw MS1 centroids
(``_reextract``): the isotope traces mono + k*1.003355/z (k < n_isotopes, at most
``area_isotopes``) are extracted at ``xic_ppm`` inside the feature's RT range and the
area is the trapezoid integral of the summed XIC with RT in minutes x ``area_scale``.
This makes ``area`` comparable between algorithms and with PEAKS "Area".

CLI::

    python -m proteoscan.features <mzML> --config configs/features_lumos_120k.yml \
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
    "traces": {},  # see proteoscan.feature_traces.DEFAULTS
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
    """Run FeatureFinderCentroided on one m/z slab (worker process).

    args = (rts_seconds, [(mz, int) per spectrum], core_lo, core_hi, openms params, omp)
    """
    rts, peaks, core_lo, core_hi, openms, omp = args
    os.environ["OMP_NUM_THREADS"] = str(omp)
    import pyopenms as oms

    sub = oms.MSExperiment()
    for rt, (mz, it) in zip(rts, peaks):
        sp = oms.MSSpectrum()
        sp.setRT(float(rt))
        sp.setMSLevel(1)
        sp.set_peaks((mz, it))
        sub.addSpectrum(sp)
    sub.updateRanges()
    ff = oms.FeatureFinderAlgorithmPicked()
    p = _param(ff.getDefaults(), openms)
    out = oms.FeatureMap()
    ff.run(sub, out, p, oms.FeatureMap())
    recs = _featuremap_records(out)
    return [r for r in recs if core_lo <= r["mz"] < core_hi]


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


def _run_centroided(exp, openms, threads, slabs, slab_overlap):
    """FFC in parallel over ``slabs`` overlapping m/z slabs of equal peak count.

    A feature is kept by the slab whose core holds its monoisotopic m/z; each slab extends
    1.5 Th below and ``slab_overlap`` Th above its core so isotope envelopes (and FFC's
    check for a -1 isotope) are complete.
    """
    rts = np.array([s.getRT() for s in exp])
    spectra = [s.get_peaks() for s in exp]
    allmz = np.concatenate([m for m, _ in spectra])
    edges = np.quantile(allmz, np.linspace(0, 1, slabs + 1))
    edges[0] -= 0.01
    edges[-1] += 0.01
    del allmz
    jobs = []
    workers = max(1, min(slabs, threads))
    omp = max(1, threads // workers)
    for i in range(slabs):
        a, b = edges[i], edges[i + 1]
        lo, hi = a - 1.5, b + slab_overlap
        pk = []
        for mz, it in spectra:
            j0, j1 = np.searchsorted(mz, [lo, hi])
            pk.append((mz[j0:j1].astype(np.float64), it[j0:j1].astype(np.float32)))
        jobs.append((rts, pk, a, b, openms, omp))
    del spectra
    recs: list[dict[str, Any]] = []
    if workers == 1:
        for j in jobs:
            recs += _ffc_slab(j)
        return recs
    # spawn (not fork): each worker must initialise OpenMP with its own thread count
    old = os.environ.get("OMP_NUM_THREADS")
    os.environ["OMP_NUM_THREADS"] = str(omp)
    try:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
            # largest slabs first for load balance
            order = sorted(range(slabs), key=lambda i: -sum(len(m) for m, _ in jobs[i][1]))
            futs = [ex.submit(_ffc_slab, jobs[i]) for i in order]
            for f in futs:
                recs += f.result()
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


def _walk_bounds_py(x, offs, apex_idx, frac, max_gap, valley):
    n_f = len(offs) - 1
    lo = np.empty(n_f, np.int64)
    hi = np.empty(n_f, np.int64)
    for f in range(n_f):
        a, b = offs[f], offs[f + 1]
        ap = apex_idx[f]
        top = x[ap]
        thr = frac * top
        for direction in (-1, 1):
            last = ap
            gap = 0
            run_min = top
            i = ap
            while True:
                i += direction
                if i < a or i >= b:
                    break
                v = x[i]
                if v >= thr and v > 0:
                    # valley split: signal rises again after dropping below valley*top
                    if valley > 0 and run_min < valley * top and v > 2.0 * run_min and v > thr * 2:
                        break
                    last = i
                    gap = 0
                    if v < run_min:
                        run_min = v
                else:
                    gap += 1
                    if v < run_min:
                        run_min = v
                    if gap > max_gap:
                        break
            if direction < 0:
                lo[f] = last
            else:
                hi[f] = last
    return lo, hi


try:  # numba is optional; the pure-python walk is ~50x slower
    from numba import njit as _njit

    _walk_bounds = _njit(cache=True)(_walk_bounds_py)
except Exception:  # pragma: no cover
    _walk_bounds = _walk_bounds_py


def _reextract(df: pd.DataFrame, ms1, xic_ppm: float = 5.0, area_isotopes: int = 2,
               area_scale: float = 1.0, apex_from: str = "xic", boundary: str = "walk",
               boundary_isotopes: int = 2, boundary_frac: float = 0.1, boundary_gap: int = 1,
               boundary_window: float = 1.0, boundary_valley: float = 0.0) -> pd.DataFrame:
    """Re-extract isotope XICs around each feature and re-derive RT range, apex and area.

    x(t) = summed XIC of the first ``boundary_isotopes`` isotope peaks (``xic_ppm``).
    The apex is the max of x inside the algorithm's RT range.  With boundary="walk",
    rt_begin/rt_end are found by walking outward from the apex while x >= boundary_frac *
    apex, tolerating up to ``boundary_gap`` consecutive scans below it, never beyond
    ``boundary_window`` minutes outside the algorithm's range (and optionally splitting at
    valleys deeper than ``boundary_valley`` * apex).  boundary="hull" keeps the algorithm's
    range.  area = trapezoid (RT in minutes) of the summed XIC of the first
    ``area_isotopes`` isotopes over [rt_begin, rt_end], times ``area_scale``.
    """
    rt = ms1["rt"]
    nsc = len(rt)
    rb0 = df["rt_begin"].to_numpy()
    re0 = df["rt_end"].to_numpy()
    w = boundary_window if boundary == "walk" else 0.0
    s0 = np.searchsorted(rt, rb0 - w - 1e-9, side="left")
    s1 = np.searchsorted(rt, re0 + w + 1e-9, side="right") - 1
    s1 = np.clip(s1, 0, nsc - 1)
    s0 = np.minimum(s0, s1)
    n = (s1 - s0 + 1).astype(np.int64)
    offs = np.r_[0, np.cumsum(n)]
    fid = np.repeat(np.arange(len(df)), n)
    scan = np.repeat(s0, n) + (np.arange(offs[-1]) - offs[:-1][fid])
    t = rt[scan]
    mono = df["mz"].to_numpy()[fid]
    z = np.maximum(df["charge"].to_numpy()[fid], 1)
    K = max(area_isotopes, boundary_isotopes)
    X = np.zeros((K, len(fid)))
    for k in range(K):
        X[k], _ = _lookup(ms1, scan, mono + k * C13_DELTA / z, xic_ppm)
    xb = X[:boundary_isotopes].sum(0)
    xa = X[:area_isotopes].sum(0)
    # apex inside the original range
    inside = (t >= rb0[fid] - 1e-9) & (t <= re0[fid] + 1e-9)
    score = np.where(inside, xb, -1.0)
    order = np.lexsort((-score, fid))
    first = np.r_[True, fid[order][1:] != fid[order][:-1]]
    apex_idx = order[first]
    if boundary == "walk":
        lo, hi = _walk_bounds(xb, offs, apex_idx, float(boundary_frac), int(boundary_gap),
                              float(boundary_valley))
    else:
        idx_in = np.where(inside)[0]
        lo = np.full(len(df), -1)
        hi = np.full(len(df), -1)
        g = pd.Series(idx_in).groupby(fid[idx_in])
        lo[g.min().index] = g.min().to_numpy()
        hi[g.max().index] = g.max().to_numpy()
        lo = np.where(lo < 0, apex_idx, lo)
        hi = np.where(hi < 0, apex_idx, hi)
    keep = (np.arange(len(fid)) >= lo[fid]) & (np.arange(len(fid)) <= hi[fid])
    nxt_same = np.r_[fid[1:] == fid[:-1], False]
    both = keep & np.r_[keep[1:], False] & nxt_same
    dt = np.r_[t[1:] - t[:-1], 0.0]
    seg = np.where(both, 0.5 * (xa + np.r_[xa[1:], 0.0]) * dt, 0.0)
    seg0 = np.where(both, 0.5 * (X[0] + np.r_[X[0][1:], 0.0]) * dt, 0.0)
    out = df.copy()
    out["area"] = np.bincount(fid, weights=seg, minlength=len(df)) * area_scale
    out["area_mono"] = np.bincount(fid, weights=seg0, minlength=len(df)) * area_scale
    out["apex_intensity"] = xb[apex_idx]
    out["rt_xic_apex"] = t[apex_idx]
    out["rt_begin"] = t[lo]
    out["rt_end"] = t[hi]
    out["ms1_scans"] = (hi - lo + 1).astype(int)
    out["ms1_scans_mono"] = np.bincount(fid, weights=(keep & (X[0] > 0)), minlength=len(df)).astype(int)
    out["rt"] = out["rt_xic_apex"] if apex_from == "xic" else out["rt_model"].clip(out["rt_begin"], out["rt_end"])
    return out


# --------------------------------------------------------------------------------------
# Main entry points
# --------------------------------------------------------------------------------------

_RAW_COLS = ["mz", "charge", "rt_model", "hull_rt_begin", "hull_rt_end", "mono_rt_begin",
             "mono_rt_end", "intensity", "quality", "n_isotopes"]


def detect_features(mzml_path: str | os.PathLike, out_featurexml: str | os.PathLike | None = None,
                    algorithm: str = "centroided", openms: dict[str, Any] | None = None,
                    threads: int = 16, slabs: int = 32, slab_overlap: float = 6.0,
                    rt_range: str = "hull", xic_ppm: float = 5.0, area_isotopes: int = 2,
                    area_scale: float = 1.0, apex_from: str = "xic", boundary: str = "walk",
                    boundary_isotopes: int = 2, boundary_frac: float = 0.1,
                    boundary_gap: int = 1, boundary_window: float = 1.0,
                    boundary_valley: float = 0.0, min_charge: int = 1, max_charge: int = 9,
                    min_ms1_scans: int = 1, min_quality: float = 0.0,
                    min_apex_intensity: float = 0.0, min_isotopes: int = 1,
                    high_charge: int = 99, high_charge_min_isotopes: int = 1,
                    dedupe_ppm: float = 5.0,
                    raw_cache: str | os.PathLike | None = None, verbose: bool = True,
                    ms1_arrays: dict | None = None, **params) -> pd.DataFrame:
    """Detect MS1 features in an mzML and return one row per feature.

    ``openms`` (and any extra keyword ``section__name=value``, mapped to
    ``section:name``) override :data:`DEFAULT_OPENMS` for the chosen algorithm.
    The algorithm's features are then re-quantified on the MS1 centroids by
    :func:`_reextract` (see there for ``xic_ppm``, ``area_*``, ``boundary*``, ``apex_from``).
    ``rt_range`` = hull | mono: which algorithm hull seeds the apex search.
    Filters: ``min_ms1_scans`` (scans in the final RT range), ``min_quality`` (algorithm
    quality; isotope-pattern fit for ``traces``), ``min_apex_intensity``, ``min_isotopes``,
    and ``high_charge_min_isotopes`` for charges >= ``high_charge``.
    Features that end up with the same charge, m/z within ``dedupe_ppm`` and identical RT
    range are merged (largest area kept).  ``raw_cache``: parquet of the algorithm's raw
    features; read instead of re-running detection when it exists, written otherwise.
    RT columns are in minutes.
    """
    t0 = time.time()
    algorithm = algorithm.lower()
    if algorithm not in DEFAULT_OPENMS:
        raise ValueError(f"unknown algorithm {algorithm!r}")
    op = dict(DEFAULT_OPENMS[algorithm])
    op.update(openms or {})
    op.update({k.replace("__", ":"): v for k, v in params.items()})
    os.environ["OMP_NUM_THREADS"] = str(threads)  # effective only if pyopenms not yet loaded

    ms1 = ms1_arrays  # optional pre-loaded _ms1_arrays(load_ms1(mzml_path))
    if raw_cache is not None and Path(raw_cache).exists():
        df = pd.read_parquet(raw_cache)
    elif algorithm == "traces":
        from .feature_traces import run_traces

        if ms1 is None:
            ms1 = _ms1_arrays(load_ms1(mzml_path))
        recs = run_traces(ms1, op, verbose=verbose)
    else:
        exp = load_ms1(mzml_path)
        if ms1 is None:
            ms1 = _ms1_arrays(exp)  # before the run: FeatureFinderMultiplex empties the map
        if algorithm == "centroided":
            recs = _run_centroided(exp, op, threads, slabs, slab_overlap)
        else:
            fm = _run_metabo(exp, op) if algorithm == "metabo" else _run_multiplex(exp, op)
            recs = _featuremap_records(fm)
        del exp
    if raw_cache is None or not Path(raw_cache).exists():
        df = pd.DataFrame(recs, columns=_RAW_COLS)
        if raw_cache is not None:
            Path(raw_cache).parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(raw_cache, index=False)
    t_detect = time.time() - t0
    n_raw = len(df)
    df = df[(df["charge"] >= min_charge) & (df["charge"] <= max_charge)].copy()
    if rt_range == "mono":
        df["rt_begin"], df["rt_end"] = df["mono_rt_begin"], df["mono_rt_end"]
    else:
        df["rt_begin"], df["rt_end"] = df["hull_rt_begin"], df["hull_rt_end"]
    df = df.sort_values(["rt_model", "mz"]).reset_index(drop=True)
    if ms1 is None:
        ms1 = _ms1_arrays(load_ms1(mzml_path))
    df = _reextract(df, ms1, xic_ppm=xic_ppm, area_isotopes=area_isotopes, area_scale=area_scale,
                    apex_from=apex_from, boundary=boundary, boundary_isotopes=boundary_isotopes,
                    boundary_frac=boundary_frac, boundary_gap=boundary_gap,
                    boundary_window=boundary_window, boundary_valley=boundary_valley)
    keep = (df["ms1_scans"] >= min_ms1_scans) & (df["quality"] >= min_quality) & \
        (df["apex_intensity"] >= min_apex_intensity) & (df["n_isotopes"] >= min_isotopes) & \
        ((df["charge"] < high_charge) | (df["n_isotopes"] >= high_charge_min_isotopes))
    df = df[keep]
    if dedupe_ppm > 0 and len(df):
        df = df.sort_values(["charge", "rt_begin", "rt_end", "mz"]).reset_index(drop=True)
        same = (df["charge"].diff() == 0) & (df["rt_begin"].diff() == 0) & \
               (df["rt_end"].diff() == 0) & (df["mz"].diff().abs() <= df["mz"] * dedupe_ppm * 1e-6)
        grp = (~same).cumsum()
        df = df.loc[df.groupby(grp)["area"].idxmax()]
    df = df.sort_values(["rt", "mz"]).reset_index(drop=True)
    df.insert(0, "feature_id", np.arange(1, len(df) + 1))
    cols = ["feature_id", "mz", "charge", "rt", "rt_begin", "rt_end", "area", "quality",
            "n_isotopes", "ms1_scans"]
    df = df[cols + [c for c in df.columns if c not in cols]]
    df.attrs["runtime_s"] = time.time() - t0
    df.attrs["detect_s"] = t_detect
    if verbose:
        print(f"[features] {algorithm}: {n_raw} raw -> {len(df)} features in "
              f"{time.time() - t0:.0f}s (detection {t_detect:.0f}s)", flush=True)
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
                                  [r.rt_begin * 60.0, r.mz + (max(r.n_isotopes, 1) - 1) * C13_DELTA / r.charge]],
                                 dtype=np.float32))
        f.setConvexHulls([h])
        f.setMetaValue("n_isotopes", int(r.n_isotopes))
        f.setMetaValue("ms1_scans", int(r.ms1_scans))
        fm.push_back(f)
    fm.setUniqueIds()
    oms.FeatureXMLFile().store(str(path), fm)


def _averagine_fractions(mz: np.ndarray, z: np.ndarray, kmax: int) -> np.ndarray:
    """Poisson (averagine) relative isotope abundances, shape (n, kmax), rows sum to <= 1."""
    lam = np.clip((mz - PROTON) * z, 0, None) * 0.000594
    k = np.arange(kmax)
    logp = -lam[:, None] + k[None, :] * np.log(np.maximum(lam, 1e-12))[:, None] \
        - np.cumsum(np.r_[0.0, np.log(np.arange(1, kmax))])[None, :]
    return np.exp(logp)


def assign_ms2_to_features(features_df: pd.DataFrame, mzml_path: str | os.PathLike | pd.DataFrame,
                           isolation_window_da: float | None = None, n_isotopes: int = 4,
                           rt_tolerance: float = 0.0, min_window_fraction: float = 0.0) -> pd.DataFrame:
    """List, for every MS2 scan, the features co-isolated with it.

    A feature is co-isolated when rt_begin - tol <= scan RT <= rt_end + tol and any of its
    isotope peaks mono + k*1.003355/z (k = 0..n_isotopes-1) lies inside the isolation
    window [target - w/2, target + w/2].  ``isolation_window_da`` is the full width; None
    uses the offsets written in the mzML.  ``mzml_path`` may also be a DataFrame from
    :func:`read_ms2_metadata`.  ``window_fraction`` is the averagine share of the
    feature's envelope inside the window; ``min_window_fraction`` > 0 drops pairs below it
    (PEAKS keeps almost no pair under 0.2; with ``rt_tolerance`` = 0.02 min this rule
    reproduces ~95% of PEAKS's (scan, feature) pairs, see the features REPORT).  Rows are
    sorted by scan, then feature area (descending); ``rank`` is 1 for the most intense
    feature.  MS2 scans with no feature get one row with feature_id = -1.
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
    fmz = f["mz"].to_numpy()
    kmax = max(n_isotopes, 8)
    frac_k = _averagine_fractions(fmz, fz, kmax)
    k = np.arange(n_isotopes)
    peak_mz = (fmz[:, None] + k[None, :] * C13_DELTA / fz[:, None]).ravel()
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
    # averagine share of the envelope inside the window (all kmax isotopes)
    pi, pf = pairs["i"].to_numpy(), pairs["f"].to_numpy()
    kk = np.arange(kmax)
    pm = fmz[pf][:, None] + kk[None, :] * C13_DELTA / fz[pf][:, None]
    inside = (pm >= lo[pi][:, None]) & (pm <= hi[pi][:, None])
    wfrac = (frac_k[pf] * inside).sum(1)
    keep = wfrac >= min_window_fraction
    pairs = pairs[keep]
    pi, pf, wfrac = pi[keep], pf[keep], wfrac[keep]
    out = pd.DataFrame({
        "scan": ms2["scan"].to_numpy()[pi],
        "ms2_rt": srt[pi],
        "isolation_target": target[pi],
        "precursor_charge": ms2["precursor_charge"].to_numpy()[pi],
        "feature_id": f["feature_id"].to_numpy()[pf],
        "feature_mz": fmz[pf],
        "feature_charge": f["charge"].to_numpy()[pf],
        "feature_rt": f["rt"].to_numpy()[pf],
        "feature_area": f["area"].to_numpy()[pf],
        "n_isotopes_in_window": pairs["n_isotopes_in_window"].to_numpy(),
        "window_fraction": wfrac,
    })
    missing = np.setdiff1d(np.arange(len(ms2)), np.unique(pi))
    if len(missing):
        empty = pd.DataFrame({
            "scan": ms2["scan"].to_numpy()[missing], "ms2_rt": srt[missing],
            "isolation_target": target[missing],
            "precursor_charge": ms2["precursor_charge"].to_numpy()[missing],
            "feature_id": -1, "feature_mz": np.nan, "feature_charge": 0,
            "feature_rt": np.nan, "feature_area": np.nan, "n_isotopes_in_window": 0,
            "window_fraction": 0.0})
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
    ms2_kw = {k[4:]: kw.pop(k) for k in list(kw) if k.startswith("ms2_")}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(a.mzml).name
    for ext in (".gz", ".mzML", ".mzml"):
        stem = stem[: -len(ext)] if stem.endswith(ext) else stem
    df = detect_features(a.mzml, out / f"{stem}.featureXML", **kw)
    df.to_parquet(out / f"{stem}.features.parquet", index=False)
    ms2 = assign_ms2_to_features(df, a.mzml, iso, **ms2_kw)
    ms2.to_parquet(out / f"{stem}.ms2_features.parquet", index=False)
    print(f"[features] wrote {out}/{stem}.features.parquet ({len(df)}) and "
          f".ms2_features.parquet ({len(ms2)} rows, "
          f"{ms2.loc[ms2.feature_id >= 0, 'scan'].nunique()} scans with a feature)")


if __name__ == "__main__":
    main()
