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


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def _longest_run(mask: np.ndarray) -> int:
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


class SpectrumStore:
    """Random access to MS2 peak lists of an mzML by raw scan number (pyopenms on-disc)."""

    def __init__(self, mzml: str | Path):
        import pyopenms as oms
        self._od = oms.OnDiscMSExperiment()
        if not self._od.openFile(str(mzml)):
            raise FileNotFoundError(mzml)
        meta = self._od.getMetaData()
        self._index = {}
        for i in range(meta.getNrSpectra()):
            m = re.search(r"scan=(\d+)", meta[i].getNativeID())
            if m and meta[i].getMSLevel() == 2:
                self._index.setdefault(int(m.group(1)), i)

    def peaks(self, scan: int):
        i = self._index.get(int(scan))
        if i is None:
            return None
        mz, it = self._od.getSpectrum(i).get_peaks()
        return np.asarray(mz, dtype=float), np.asarray(it, dtype=float)


def observed_log2(theo_mz: np.ndarray, mz: np.ndarray, it: np.ndarray, tol: float, tic: float) -> np.ndarray:
    """MS2PIP-style observed targets: log2(max intensity within +/-tol of each ion / TIC + 0.001)."""
    lo = np.searchsorted(mz, theo_mz - tol, side="left")
    hi = np.searchsorted(mz, theo_mz + tol, side="right")
    out = np.full(len(theo_mz), 0.001)
    for k in range(len(theo_mz)):
        if hi[k] > lo[k]:
            out[k] += it[lo[k]:hi[k]].max() / tic
    return np.log2(out)


def ms2pip_correlations(df: pd.DataFrame, mzml: str | Path, model: str = "HCD2021",
                        ms2_tolerance: float = 0.5, processes: int = 8) -> pd.DataFrame:
    """Predicted-vs-observed fragment features aligned with df.
    Predictions and theoretical m/z come from MS2PIP's predict_batch (deterministic and
    correctly attached to each PSM); the observed intensities are extracted here, because in
    ms2pip 4.2 `correlate` attaches observed/predicted vectors inconsistently when several
    PSMs share a spectrum.
    Columns: ms2_corr (Pearson on log2 targets over all b/y ions), ms2_cos (cosine on linear
    intensities), ms2_frac_pred_found, longest_b_obs, longest_y_obs, n_ions_obs, residue_coverage."""
    from ms2pip.core import predict_batch
    psms = to_psmlist(df)
    preds = predict_batch(psms, model=model, processes=processes)
    store = SpectrumStore(mzml)
    n = len(df)
    corr = np.zeros(n); cos = np.zeros(n); frac = np.zeros(n)
    longest_b = np.zeros(n); longest_y = np.zeros(n); n_obs = np.zeros(n); cov_both = np.zeros(n)
    scans = df["scan"].to_numpy()
    order = np.argsort(scans, kind="stable")
    cache_scan, cache = None, None
    by_index = {r.psm_index: r for r in preds}
    for i in order:
        r = by_index.get(int(i))
        if r is None or r.predicted_intensity is None or r.theoretical_mz is None:
            continue
        sc = int(scans[i])
        if sc != cache_scan:
            cache_scan, cache = sc, store.peaks(sc)
        if cache is None or len(cache[0]) == 0:
            continue
        mz, it = cache
        tic = float(it.sum()) or 1.0
        ions = [k for k in ("b", "y") if k in r.predicted_intensity]
        pred = np.concatenate([np.asarray(r.predicted_intensity[k], dtype=float) for k in ions])
        theo = np.concatenate([np.asarray(r.theoretical_mz[k], dtype=float) for k in ions])
        obs = observed_log2(theo, mz, it, ms2_tolerance, tic)
        corr[i] = _pearson(obs, pred)
        o_lin = np.clip(np.power(2.0, obs) - 0.001, 0, None); p_lin = np.clip(np.power(2.0, pred) - 0.001, 0, None)
        denom = np.linalg.norm(o_lin) * np.linalg.norm(p_lin)
        cos[i] = float(np.dot(o_lin, p_lin) / denom) if denom > 0 else 0.0
        strong = p_lin >= 0.05 * (p_lin.max() if p_lin.max() > 0 else 1)
        frac[i] = float((o_lin[strong] > 0).mean()) if strong.any() else 0.0
        nb = len(r.predicted_intensity.get("b", []))
        ob = obs[:nb] > -9.9; oy = obs[nb:] > -9.9
        longest_b[i] = _longest_run(ob); longest_y[i] = _longest_run(oy)
        n_obs[i] = int(ob.sum() + oy.sum())
        if nb and len(oy):
            nres = nb + 1
            b_bound = np.concatenate([[True], ob, [True]])
            y_bound = np.concatenate([[True], oy[::-1], [True]])
            covered = [(b_bound[k - 1] or y_bound[k - 1]) and (b_bound[k] or y_bound[k]) for k in range(1, nres + 1)]
            cov_both[i] = float(np.mean(covered))
    return pd.DataFrame({"ms2_corr": corr, "ms2_cos": cos, "ms2_frac_pred_found": frac,
                         "longest_b_obs": longest_b, "longest_y_obs": longest_y, "n_ions_obs": n_obs,
                         "residue_coverage": cov_both}, index=df.index)


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
        ms2 = ms2pip_correlations(out, mzml, ms2pip_model, ms2_tolerance, processes)
        for c in ms2.columns:
            out[c] = ms2[c].fillna(0.0).to_numpy()
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
                                   "ms2_corr", "ms2_cos", "ms2_frac_pred_found", "longest_b_obs", "longest_y_obs",
                                   "n_ions_obs", "residue_coverage", "delta_rt", "abs_delta_rt"]
