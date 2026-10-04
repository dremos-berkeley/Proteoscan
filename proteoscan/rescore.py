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


def _cv_scores(ds, models, test_fdr: float, seed: int) -> np.ndarray:
    """Held-out scores exactly as mokapot.brew assigns them: brew splits the PSMs into folds with a
    generator seeded by `seed` (the split is its first use of that generator) and scores each fold with
    the model trained on the other folds, calibrated per fold (Granholm et al. 2012). Re-creating the
    split here relies on mokapot 0.10 internals (`_split`, `_calibrate_scores`, `_data`); the caller
    checks the result against brew's own scores."""
    import copy
    ds.rng = np.random.default_rng(seed)
    folds = ds._split(len(models))
    scores = np.full(len(ds.data), np.nan)
    sub = copy.copy(ds)
    for idx, m in zip(folds, models):
        idx = list(idx)
        sub._data = ds.data.loc[idx, :]
        scores[idx] = sub._calibrate_scores(m.predict(sub), test_fdr)
    return scores


def rescore_tda(df: pd.DataFrame, feature_cols=None, seed: int = 0) -> pd.DataFrame:
    """mokapot semi-supervised rescoring; adds `score` (cross-validated: every PSM is scored by the fold
    model that did not train on it, as mokapot reports) and `mokapot_q` (mokapot's PSM-level q-value for
    the PSMs that won their spectrum in mokapot's competition; NaN for the others). The returned frame
    has a fresh RangeIndex."""
    import mokapot
    work = df.reset_index(drop=True)  # mokapot indexes folds positionally via .loc
    X, cols = _features(work, feature_cols)
    work = work.copy()
    work["specid"] = np.arange(len(work))
    work["target"] = (~work["is_decoy"]).astype(int)
    work["scan_key"] = work["source_file"].astype(str) + ":" + work["scan"].astype(str)
    for c in cols:
        work[c] = X[c].to_numpy()
    log.info("rescore_tda: %d PSMs, %d features: %s", len(work), len(cols), ", ".join(cols))
    ds = mokapot.LinearPsmDataset(
        work[["specid", "target", "scan_key", "peptide_sage", "proteins"] + cols],
        target_column="target", spectrum_columns="scan_key", peptide_column="peptide_sage",
        protein_column="proteins", feature_columns=cols, copy_data=False)
    scores = None
    mokapot_q = np.full(len(work), np.nan)
    for test_fdr in (0.01, 0.05, 0.10):
        try:
            results, models = mokapot.brew(ds, test_fdr=test_fdr, max_workers=4, rng=seed)
        except RuntimeError as e:  # "No PSMs found below the eval_fdr" on small/weak inputs
            log.warning("mokapot failed at test_fdr %.2f (%s); relaxing", test_fdr, e)
            continue
        win = pd.concat([results.confidence_estimates["psms"], results.decoy_confidence_estimates["psms"]])
        mokapot_q[win["specid"].to_numpy(int)] = win["mokapot q-value"].to_numpy(float)
        brew_score = pd.Series(win["mokapot score"].to_numpy(float), index=win["specid"].to_numpy(int))
        try:
            cv = _cv_scores(ds, models, test_fdr, seed)
            ok = np.isfinite(cv).all() and np.allclose(cv[brew_score.index.to_numpy()], brew_score.to_numpy(),
                                                       rtol=1e-6, atol=1e-9)
        except Exception as e:  # private mokapot API changed
            log.warning("could not reproduce mokapot's cross-validated scores (%s)", e)
            ok = False
        if ok:
            scores = cv
            log.info("rescore_tda: cross-validated scores reproduce mokapot's (%d spectrum winners checked)",
                     len(brew_score))
        else:
            # brew fell back to its best single feature / a reset model, or the private API moved:
            # keep the previous behaviour (mean of the fold models, in-sample) and say so
            log.warning("rescore_tda: cross-validated scores do NOT match mokapot's; falling back to the "
                        "mean of the fold models (in-sample scores)")
            scores = np.zeros(len(work))
            for m in models:
                scores += m.predict(ds)
            scores /= max(len(models), 1)
        break
    if scores is None:
        # last resort: a target/decoy logistic regression on standardised features
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        Xs = StandardScaler().fit_transform(X.to_numpy())
        clf = LogisticRegression(max_iter=500).fit(Xs, work["target"].to_numpy())
        scores = clf.decision_function(Xs)
        log.warning("mokapot unusable; used logistic regression target/decoy rescoring")
    work["score"] = scores
    work["mokapot_q"] = mokapot_q
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


TAIL_MIN_DECOYS = 20       # decoys used for the tail fit: max(TAIL_MIN_DECOYS, TAIL_FRACTION * n)
TAIL_FRACTION = 0.01
TAIL_MIN_SLOPE = 1e-6      # keeps the extrapolation strictly increasing if the fitted slope is <= 0


def minus10lgp(score: pd.Series, is_decoy: pd.Series) -> np.ndarray:
    """-10*log10(p), p = fraction of decoy scores >= score (empirical null with +1 smoothing).

    At or below the best decoy this is the empirical rank transform. Above it the empirical p would be
    stuck at 1/(n+1), so every such target would tie at the ceiling 10*log10(n+1); instead the decoy
    tail is extrapolated exponentially: a line of -10*log10(p_emp) against score is fitted to the top
    k = max(20, 1% of n) decoys (least squares; slope clipped to a small positive minimum) and continued
    from the ceiling at the best decoy, so the transform is continuous from the right there and
    monotone. With fewer than 2 decoys (or identical top-decoy scores) the old flat ceiling is kept;
    with no decoys the result is NaN.
    """
    s = np.asarray(score, dtype=float)
    d = np.sort(np.asarray(score[np.asarray(is_decoy, dtype=bool)], dtype=float))
    d = d[~np.isnan(d)]
    n = len(d)
    if n == 0:
        return np.full(len(s), np.nan)
    ge = n - np.searchsorted(d, s, side="left")  # decoys with score >= s
    p = (ge + 1.0) / (n + 1.0)
    out = -10.0 * np.log10(p)
    k = min(n, max(TAIL_MIN_DECOYS, int(np.ceil(TAIL_FRACTION * n))))
    top = d[n - k:]
    if k >= 2 and np.ptp(top) > 0:
        y = -10.0 * np.log10((n - np.searchsorted(d, top, side="left") + 1.0) / (n + 1.0))
        slope = max(float(np.polyfit(top, y, 1)[0]), TAIL_MIN_SLOPE)
        above = s > d[-1]
        ceiling = 10.0 * np.log10(n + 1.0)
        out[above] = ceiling + slope * (s[above] - d[-1])
    return out


def best_per_spectrum(df: pd.DataFrame, score_col: str = "score") -> pd.DataFrame:
    """PEAKS keeps one peptide per spectrum (per precursor instance)."""
    keys = ["source_file", "scan"]
    if "precursor_idx" in df.columns:
        keys.append("precursor_idx")
    idx = df.groupby(keys)[score_col].idxmax()
    return df.loc[idx].reset_index(drop=True)
