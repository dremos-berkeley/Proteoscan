"""Stage 2b: PEAKS-style precursor list ("data refinement").

For every MS2 scan PEAKS emits one query per candidate precursor (docs/PEAKS_PRECURSOR_RULES.md):
  * feature-derived entries: every MS1 feature (from features.detect_features) whose elution
    range covers the scan RT (+/-0.02 min) and whose monoisotopic m/z lies inside the isolation
    window (recall 0.90 of PEAKS's links; precision ~0.5). m/z, charge and RT are the feature's.
  * a featureless entry: the instrument charge with an averagine-corrected monoisotopic m/z from
    the master MS1 (monoiso.correct_precursor), unless a feature entry of the same charge lies
    within 12 ppm.
The list is materialised as an "expanded" mzML in which each MS2 spectrum is repeated once per
precursor, with native id "... scan=N precursor=K", so any engine can search it.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms as oms

from .monoiso import C13, correct_precursor

SCAN_RE = re.compile(r"scan=(\d+)")
TARGET_RE = re.compile(r"ms2 ([0-9.]+)@")


def _scan(nid: str) -> int:
    return int(SCAN_RE.search(nid).group(1))


def load_experiment(mzml: str | Path) -> oms.MSExperiment:
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(mzml), exp)
    return exp


def ms2_metadata(exp: oms.MSExperiment) -> pd.DataFrame:
    """One row per MS2: scan, rt_min, target_mz, half_width, charge, trailer_mono, master_ms1 (index)."""
    rows = []
    ms1_indices: list[int] = []
    ms1_by_scan: dict[int, int] = {}
    for i in range(exp.getNrSpectra()):
        s = exp[i]
        if s.getMSLevel() == 1:
            ms1_indices.append(i)
            ms1_by_scan[_scan(s.getNativeID())] = i
            continue
        if s.getMSLevel() != 2 or not s.getPrecursors():
            continue
        p = s.getPrecursors()[0]
        nid = s.getNativeID()
        filt = s.getMetaValue("filter string") if s.metaValueExists("filter string") else None
        if filt is None and s.getAcquisitionInfo().size():
            a = s.getAcquisitionInfo()[0]
            filt = a.getMetaValue("filter string") if a.metaValueExists("filter string") else None
        m = TARGET_RE.search(str(filt)) if filt else None
        target = float(m.group(1)) if m else p.getMZ()
        ref = p.getMetaValue("spectrum_ref") if p.metaValueExists("spectrum_ref") else None
        master = None
        if ref:
            rs = SCAN_RE.search(str(ref))
            if rs:
                master = ms1_by_scan.get(int(rs.group(1)))
        if master is None:
            master = ms1_indices[-2] if len(ms1_indices) >= 2 else (ms1_indices[-1] if ms1_indices else None)
        trailer = None
        for key in ("[Thermo Trailer Extra]Monoisotopic M/Z:",):
            if s.metaValueExists(key):
                trailer = float(s.getMetaValue(key))
        rows.append(dict(index=i, scan=_scan(nid), rt_min=s.getRT() / 60.0, precursor_mz=p.getMZ(),
                         target_mz=target, half_width=max(p.getIsolationWindowLowerOffset(), p.getIsolationWindowUpperOffset()) or 0.65,
                         charge=p.getCharge(), trailer_mono=trailer, master_ms1=master))
    return pd.DataFrame(rows)


def link_features(ms2: pd.DataFrame, features: pd.DataFrame, rt_slack: float = 0.02,
                  require_mono_in_window: bool = True, max_isotope: int = 4) -> pd.DataFrame:
    """Return rows (scan, feature_id, mz, charge, rt, area) for every feature linked to an MS2."""
    f = features.sort_values("mz").reset_index(drop=True)
    fmz = f["mz"].to_numpy()
    out = []
    for r in ms2.itertuples(index=False):
        lo, hi = r.target_mz - r.half_width, r.target_mz + r.half_width
        # candidate features: mono m/z in [lo - max_isotope*C13, hi]
        i0 = np.searchsorted(fmz, lo - max_isotope * C13)
        i1 = np.searchsorted(fmz, hi, side="right")
        for j in range(i0, i1):
            fr = f.iloc[j]
            if not (fr.rt_begin - rt_slack <= r.rt_min <= fr.rt_end + rt_slack):
                continue
            z = int(fr.charge)
            if require_mono_in_window:
                ok = lo <= fr.mz <= hi
            else:
                ok = any(lo <= fr.mz + k * C13 / z <= hi for k in range(0, max_isotope + 1))
            if not ok:
                continue
            out.append(dict(scan=int(r.scan), feature_id=int(fr.feature_id), mz=float(fr.mz), charge=z,
                            rt=float(fr.rt), area=float(fr.area)))
    return pd.DataFrame(out, columns=["scan", "feature_id", "mz", "charge", "rt", "area"])


def build_precursor_list(exp: oms.MSExperiment, features: pd.DataFrame | None, ms2: pd.DataFrame | None = None,
                         dedup_ppm: float = 12.0, correct_featureless: bool = True) -> pd.DataFrame:
    if ms2 is None:
        ms2 = ms2_metadata(exp)
    linked = link_features(ms2, features) if features is not None and len(features) else pd.DataFrame(
        columns=["scan", "feature_id", "mz", "charge", "rt", "area"])
    by_scan = {k: g for k, g in linked.groupby("scan")} if len(linked) else {}
    ms1_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    entries = []
    for r in ms2.itertuples(index=False):
        k = 0
        feats = by_scan.get(int(r.scan))
        if feats is not None:
            for fr in feats.itertuples(index=False):
                entries.append(dict(scan=int(r.scan), precursor_idx=k, mz=fr.mz, charge=fr.charge, rt=fr.rt,
                                    feature_id=fr.feature_id, area=fr.area, source="feature", shift=0))
                k += 1
        # featureless entry
        z = int(r.charge) if r.charge else 2
        mz, shift = float(r.precursor_mz), 0
        if correct_featureless and r.master_ms1 is not None and not pd.isna(r.master_ms1):
            mi = int(r.master_ms1)
            if mi not in ms1_cache:
                m, it = exp[mi].get_peaks()
                ms1_cache[mi] = (np.asarray(m, dtype=float), np.asarray(it, dtype=float))
                if len(ms1_cache) > 8:
                    ms1_cache.pop(next(iter(ms1_cache)))
            m, it = ms1_cache[mi]
            c = correct_precursor(m, it, float(r.target_mz), z)
            mz, shift = c.mz, c.shift
        dup = False
        if feats is not None:
            same = feats[feats["charge"] == z]
            if len(same) and (np.abs(same["mz"].to_numpy() - mz) / mz * 1e6 <= dedup_ppm).any():
                dup = True
        if not dup:
            entries.append(dict(scan=int(r.scan), precursor_idx=k, mz=mz, charge=z, rt=float(r.rt_min),
                                feature_id=0, area=np.nan, source="featureless", shift=shift))
    return pd.DataFrame(entries)


def write_expanded_mzml(exp: oms.MSExperiment, entries: pd.DataFrame, out_path: str | Path,
                        remove_precursor_region: bool = True, region_da: float = 0.5) -> Path:
    """Write an mzML with every MS2 repeated once per precursor entry (MS1 spectra are kept)."""
    by_scan = {k: g.sort_values("precursor_idx") for k, g in entries.groupby("scan")}
    out = oms.MSExperiment()
    for i in range(exp.getNrSpectra()):
        s = exp[i]
        if s.getMSLevel() != 2:
            out.addSpectrum(s)
            continue
        scan = _scan(s.getNativeID())
        g = by_scan.get(scan)
        if g is None:
            continue
        for e in g.itertuples(index=False):
            c = oms.MSSpectrum(s)
            p = c.getPrecursors()[0]
            p.setMZ(float(e.mz))
            p.setCharge(int(e.charge))
            c.setPrecursors([p])
            c.setNativeID(f"{s.getNativeID()} precursor={int(e.precursor_idx)}")
            if remove_precursor_region:
                mzs, ints = c.get_peaks()
                keep = np.abs(mzs - float(e.mz)) > region_da + 0.7
                c.set_peaks((mzs[keep], ints[keep]))
            out.addSpectrum(c)
    out.updateRanges()
    oms.MzMLFile().store(str(out_path), out)
    return Path(out_path)
