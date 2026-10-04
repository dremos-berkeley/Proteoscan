"""Tests for proteoscan.rescore: tail-extrapolated -10lgP and cross-validated mokapot scores."""
import numpy as np
import pandas as pd
import pytest

from proteoscan.rescore import minus10lgp


def _old_minus10lgp(s, d):
    s = np.asarray(s, float); dd = np.sort(s[np.asarray(d, bool)]); n = len(dd)
    ge = n - np.searchsorted(dd, s, side="left")
    return -10 * np.log10((ge + 1.0) / (n + 1.0))


def _table(n_t=3000, n_d=2000, seed=0):
    rng = np.random.default_rng(seed)
    s = np.r_[rng.normal(2.0, 1.5, n_t), rng.normal(0.0, 1.0, n_d)]
    d = np.r_[np.zeros(n_t, bool), np.ones(n_d, bool)]
    return pd.Series(s), pd.Series(d)


def test_unchanged_at_or_below_best_decoy():
    s, d = _table()
    new, old = minus10lgp(s, d), _old_minus10lgp(s, d)
    below = s.to_numpy() <= s[d].max()
    assert below.sum() > 0 and (~below).sum() > 0
    np.testing.assert_allclose(new[below], old[below])


def test_continuous_and_monotone_above():
    s, d = _table()
    dmax, n = s[d].max(), int(d.sum())
    grid = np.linspace(s.min(), s.max() + 5, 4000)
    gs = pd.Series(np.r_[s.to_numpy(), grid]); gd = pd.Series(np.r_[d.to_numpy(), np.zeros(len(grid), bool)])
    v = minus10lgp(gs, gd)[len(s):]
    assert np.all(np.diff(v) >= -1e-12)  # monotone non-decreasing
    above = grid > dmax
    assert np.all(np.diff(v[above]) > 0)  # strictly increasing above the best decoy (no ceiling ties)
    ceiling = 10 * np.log10(n + 1)
    eps = 1e-9
    just_above = minus10lgp(pd.Series(np.r_[s.to_numpy(), dmax + eps]), pd.Series(np.r_[d.to_numpy(), False]))[-1]
    assert just_above == pytest.approx(ceiling, abs=1e-6)  # joins the old ceiling (right-continuous)
    assert v[above].min() >= ceiling


def test_no_decoys_nan():
    v = minus10lgp(pd.Series([1.0, 2.0]), pd.Series([False, False]))
    assert np.isnan(v).all()


def test_too_few_decoys_flat():
    s = pd.Series([5.0, 6.0, 7.0, 0.0]); d = pd.Series([False, False, False, True])
    v = minus10lgp(s, d)
    assert v[0] == v[1] == v[2] == pytest.approx(10 * np.log10(2))


def _synthetic_psms(n_spec=3000, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_spec):
        good = rng.random() < 0.5
        for rank in (1, 2):
            decoy = bool(rng.random() < 0.5) if (rank == 2 or not good) else False
            mu = 4.0 if (good and rank == 1 and not decoy) else 0.0
            rows.append(dict(source_file="f.mzML", scan=i, peptide_sage=f"PEP{i}K{rank}", stripped=f"PEP{i}K",
                             proteins=("#DECOY#" if decoy else "") + f"sp|P{i}|X", is_decoy=decoy,
                             f1=rng.normal(mu, 1.0), f2=rng.normal(mu * 0.6, 1.0), f3=rng.normal(0, 1.0)))
    return pd.DataFrame(rows)


def test_rescore_tda_cv_scores_match_brew(caplog):
    mokapot = pytest.importorskip("mokapot")
    from proteoscan.rescore import _cv_scores, rescore_tda
    df = _synthetic_psms()
    df.index = df.index + 1000  # non-positional index: rescore_tda must reset it
    with caplog.at_level("INFO", logger="proteoscan.rescore"):
        out = rescore_tda(df, feature_cols=["f1", "f2", "f3"], seed=0)
    assert "reproduce mokapot" in caplog.text and "do NOT match" not in caplog.text
    assert list(out.index) == list(range(len(df)))
    assert out["score"].notna().all() and out["mokapot_q"].notna().any()
    # scores differ from the in-sample mean of the fold models
    w = df.reset_index(drop=True).copy()
    w["specid"] = np.arange(len(w)); w["target"] = (~w["is_decoy"]).astype(int)
    w["scan_key"] = w["source_file"] + ":" + w["scan"].astype(str)
    ds = mokapot.LinearPsmDataset(w[["specid", "target", "scan_key", "peptide_sage", "proteins", "f1", "f2", "f3"]],
                                  target_column="target", spectrum_columns="scan_key", peptide_column="peptide_sage",
                                  protein_column="proteins", feature_columns=["f1", "f2", "f3"], copy_data=False)
    res, models = mokapot.brew(ds, test_fdr=0.01, max_workers=4, rng=0)
    cv = _cv_scores(ds, models, 0.01, 0)
    est = res.confidence_estimates["psms"]
    np.testing.assert_allclose(cv[est["specid"].to_numpy(int)], est["mokapot score"].to_numpy(float), rtol=1e-6)
    # rescore_tda's own self-check passed (log above); its q-values exist for exactly one PSM per spectrum
    # (mokapot's spectrum winners, targets and decoys). Note: two brew calls with the same seed are not
    # bit-identical when max_workers > 1 (the fold models share one generator across threads), so the
    # scores of `out` are not compared with this second brew call.
    assert out["mokapot_q"].notna().sum() == df["scan"].nunique()
    assert out.groupby("scan")["mokapot_q"].count().eq(1).all()
