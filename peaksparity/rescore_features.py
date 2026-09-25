"""Stage 5a: per-PSM features for rescoring, imitating PEAKS's Deep Learning Boost inputs.

Given the canonical PSM candidate table (from sage_search.read_sage / merge_candidates) and the
mzML, compute:
  ms2_corr        Pearson correlation between observed and MS2PIP-predicted b/y intensities
                  (PEAKS: "MS2 Correlation")
  pred_rt, delta_rt, abs_delta_rt   DeepLC prediction calibrated per run (PEAKS: "Delta RT")
  sequence / search features: length, n_mods, missed cleavages, semi-enzymatic flag, charge,
                  precursor ppm, isotope error, hyperscore, delta_next, delta_best, matched
                  peaks, longest b/y series, matched intensity fraction, scored candidates.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
from psm_utils import PSM, PSMList

log = logging.getLogger(__name__)

SAGE_FEATURES = [
    "hyperscore", "delta_next", "delta_best", "matched_peaks", "longest_b", "longest_y",
    "longest_y_pct", "matched_intensity_pct", "scored_candidates", "poisson", "precursor_ppm",
    "fragment_ppm", "isotope_error", "missed_cleavages", "semi_enzymatic", "peptide_len", "charge",
]


def to_psmlist(df: pd.DataFrame, run_name: str | None = None) -> PSMList:
    psms = []
    for r in df.itertuples(index=False):
        psms.append(PSM(
            peptidoform=f"{r.peptide_sage}/{int(r.charge)}",
            spectrum_id=str(int(r.scan)),
            run=run_name or str(r.source_file),
            is_decoy=bool(r.is_decoy),
            score=float(r.hyperscore),
            precursor_mz=float(r.exp_mz) if hasattr(r, "exp_mz") else None,
            retention_time=float(r.rt) if hasattr(r, "rt") else None,
            rank=int(r.rank) if hasattr(r, "rank") else None,
        ))
    return PSMList(psm_list=psms)


def ms2pip_correlations(df: pd.DataFrame, mzml: str | Path, model: str = "HCD2021",
                        ms2_tolerance: float = 0.5, processes: int = 8) -> np.ndarray:
    """Return an array aligned with df of Pearson correlations (NaN when unavailable)."""
    from ms2pip.core import correlate
    psms = to_psmlist(df)
    results = correlate(psms, spectrum_file=str(mzml), spectrum_id_pattern=r"scan=(\d+)",
                        compute_correlations=True, model=model, ms2_tolerance=ms2_tolerance,
                        ms2_tolerance_mode="Da", processes=processes)
    out = np.full(len(df), np.nan)
    for i, res in enumerate(results):
        if res is not None and res.correlation is not None:
            out[i] = res.correlation
    return out


def deeplc_predictions(df: pd.DataFrame, calibrate_on: pd.Series | None = None) -> np.ndarray:
    """Predicted RT (minutes) per PSM, calibrated on the confident target PSMs of the same run.
    DeepLC's auto-calibration picks the best target PSMs by score when no reference is given."""
    import deeplc
    psms = to_psmlist(df)
    if calibrate_on is not None and calibrate_on.sum() >= 50:
        ref = PSMList(psm_list=[p for p, keep in zip(psms, calibrate_on.to_numpy()) if keep])
        preds = deeplc.predict_and_calibrate(psms, psm_list_reference=ref)
    else:
        preds = deeplc.predict_and_calibrate(psms)
    return np.asarray(preds, dtype=float)


def sequence_features(df: pd.DataFrame) -> pd.DataFrame:
    f = pd.DataFrame(index=df.index)
    f["n_mods"] = df["peptide_sage"].str.count(r"\[")
    f["n_ox"] = df["peptide_sage"].str.count(r"M\[\+15\.99")
    f["n_deam"] = df["peptide_sage"].str.count(r"[NQ]\[\+0\.98")
    f["has_cys"] = df["stripped"].str.contains("C").astype(int)
    f["len_per_mc"] = df["peptide_len"] / (df["missed_cleavages"] + 1)
    f["len_per_mod"] = df["peptide_len"] / (f["n_mods"] + 1)
    f["abs_ppm"] = df["precursor_ppm"].abs()
    f["log_scored_candidates"] = np.log1p(df["scored_candidates"].clip(lower=0))
    return f


def compute_features(df: pd.DataFrame, mzml: str | Path, ms2pip_model: str = "HCD2021",
                     ms2_tolerance: float = 0.5, processes: int = 8, do_rt: bool = True,
                     do_ms2: bool = True) -> pd.DataFrame:
    """Append feature columns to a copy of df (one run)."""
    out = df.copy()
    feats = sequence_features(out)
    for c in feats.columns:
        out[c] = feats[c]
    if do_ms2:
        log.info("ms2pip correlations (%s, tol %.2f Da) for %d PSMs", ms2pip_model, ms2_tolerance, len(out))
        out["ms2_corr"] = ms2pip_correlations(out, mzml, ms2pip_model, ms2_tolerance, processes)
        out["ms2_corr"] = out["ms2_corr"].fillna(0.0)
    if do_rt:
        log.info("DeepLC predictions for %d PSMs", len(out))
        top = (out["rank"] == 1) & (~out["is_decoy"]) if "rank" in out else ~out["is_decoy"]
        conf = top & (out["hyperscore"] >= out.loc[top, "hyperscore"].quantile(0.7))
        out["pred_rt"] = deeplc_predictions(out, calibrate_on=conf)
        out["delta_rt"] = out["rt"] - out["pred_rt"]
        out["abs_delta_rt"] = out["delta_rt"].abs()
    return out


FEATURE_COLUMNS = SAGE_FEATURES + ["n_mods", "n_ox", "n_deam", "has_cys", "len_per_mc",
                                   "len_per_mod", "abs_ppm", "log_scored_candidates",
                                   "ms2_corr", "delta_rt", "abs_delta_rt"]
