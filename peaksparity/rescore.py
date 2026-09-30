"""Stage 5b: rescoring.

Two modes:
  tda       semi-supervised target-decoy rescoring (mokapot; Percolator-style), the honest
            default when no PEAKS reference exists.
  imitate   supervised model trained to reproduce PEAKS's accepted-PSM set: label = 1 if PEAKS
            reported the same (scan, modified sequence, charge) in its export, else 0. Trained
            on reference runs and applied to new runs; never trained and evaluated on the same
            run. Decoys are kept so that the resulting score still supports decoy-based FDR.

Both produce `score` (higher is better) and `minus10lgp`, a -10*log10(p) style transform where
p is the empirical fraction of decoy PSMs scoring at least as well (PEAKS reports a p-value
based significance; we imitate its scale so that thresholds such as "protein -10lgP >= 15"
mean roughly the same thing).
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .rescore_features import FEATURE_COLUMNS

log = logging.getLogger(__name__)


def _features(df: pd.DataFrame, cols=None) -> tuple[pd.DataFrame, list[str]]:
    cols = [c for c in (cols or FEATURE_COLUMNS) if c in df.columns]
    X = df[cols].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return X, cols


def rescore_tda(df: pd.DataFrame, feature_cols=None, seed: int = 0) -> pd.DataFrame:
    """mokapot semi-supervised rescoring; adds score, mokapot_q (PSM level)."""
    import mokapot
    X, cols = _features(df, feature_cols)
    work = df.copy()
    work["specid"] = np.arange(len(work))
    work["target"] = (~work["is_decoy"]).astype(int)
    work["scan_key"] = work["source_file"].astype(str) + ":" + work["scan"].astype(str)
    for c in cols:
        work[c] = X[c].to_numpy()
    ds = mokapot.LinearPsmDataset(
        work[["specid", "target", "scan_key", "peptide_sage", "proteins"] + cols],
        target_column="target", spectrum_columns="scan_key", peptide_column="peptide_sage",
        protein_column="proteins", feature_columns=cols, copy_data=False)
    results, models = mokapot.brew(ds, test_fdr=0.01, max_workers=4, rng=seed)
    psms = results.confidence_estimates["psms"] if hasattr(results, "confidence_estimates") else results.psms
    # mokapot returns only targets in its confidence tables; recover scores for every row from the model
    scores = np.zeros(len(work))
    for m in models:
        scores += m.predict(ds)
    scores /= max(len(models), 1)
    work["score"] = scores
    work["minus10lgp"] = minus10lgp(work["score"], work["is_decoy"])
    return work.drop(columns=["specid", "target", "scan_key"])


def train_imitation(df: pd.DataFrame, label: pd.Series, feature_cols=None, seed: int = 0):
    """LightGBM classifier: label 1 = PEAKS reported this PSM. Decoys get label 0."""
    import lightgbm as lgb
    X, cols = _features(df, feature_cols)
    model = lgb.LGBMClassifier(n_estimators=600, learning_rate=0.03, num_leaves=63,
                               subsample=0.8, colsample_bytree=0.8, min_child_samples=50,
                               random_state=seed, verbose=-1)
    model.fit(X, label.astype(int))
    model.feature_cols_ = cols
    return model


def apply_imitation(model, df: pd.DataFrame) -> pd.DataFrame:
    X, _ = _features(df, model.feature_cols_)
    missing = [c for c in model.feature_cols_ if c not in X.columns]
    if missing:  # e.g. a model trained with de novo features applied to a run without them
        log.warning("imitation model: %d feature(s) absent from this run, set to 0: %s", len(missing), missing)
    X = X.reindex(columns=model.feature_cols_, fill_value=0.0)
    work = df.copy()
    p = model.predict_proba(X)[:, 1]
    work["score"] = np.log(p + 1e-9) - np.log(1 - p + 1e-9)
    work["minus10lgp"] = minus10lgp(work["score"], work["is_decoy"])
    return work


def minus10lgp(score: pd.Series, is_decoy: pd.Series) -> np.ndarray:
    """-10*log10(p), p = fraction of decoy scores >= score (empirical null with +1 smoothing)."""
    s = np.asarray(score, dtype=float)
    d = np.sort(np.asarray(score[is_decoy.to_numpy()], dtype=float))
    n = len(d)
    if n == 0:
        return np.full(len(s), np.nan)
    ge = n - np.searchsorted(d, s, side="left")  # decoys with score >= s
    p = (ge + 1.0) / (n + 1.0)
    return -10.0 * np.log10(p)


def best_per_spectrum(df: pd.DataFrame, score_col: str = "score") -> pd.DataFrame:
    """PEAKS keeps one peptide per spectrum (per precursor instance)."""
    keys = ["source_file", "scan"]
    if "precursor_idx" in df.columns:
        keys.append("precursor_idx")
    idx = df.groupby(keys)[score_col].idxmax()
    return df.loc[idx].reset_index(drop=True)
