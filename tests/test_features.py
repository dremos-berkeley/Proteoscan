"""Synthetic tests for peaksparity.features / feature_traces (no reference data needed)."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("numba")

from peaksparity.feature_traces import run_traces  # noqa: E402
from peaksparity.features import C13_DELTA, _reextract, assign_ms2_to_features  # noqa: E402


def _synthetic_ms1(peptides, n_scans=200, dt_min=0.01, noise=50, seed=0):
    """peptides: list of (mono_mz, z, apex_scan, sigma_scans, height)."""
    rng = np.random.default_rng(seed)
    rts = np.arange(n_scans) * dt_min
    mzs, ints, offs = [], [], [0]
    for s in range(n_scans):
        pk_mz, pk_int = [], []
        for mono, z, apex, sig, h in peptides:
            mass = (mono - 1.00728) * z
            lam = mass * 0.000594
            env = np.exp(-lam) * lam ** np.arange(5) / np.array([1, 1, 2, 6, 24])
            prof = h * np.exp(-0.5 * ((s - apex) / sig) ** 2)
            for k in range(5):
                v = prof * env[k] / env.max()
                if v > 5 * noise:
                    pk_mz.append(mono + k * C13_DELTA / z + rng.normal(0, mono * 1e-6))
                    pk_int.append(v)
        for _ in range(20):  # random noise centroids
            pk_mz.append(rng.uniform(400, 1000))
            pk_int.append(rng.uniform(noise, 5 * noise))
        o = np.argsort(pk_mz)
        mzs.append(np.asarray(pk_mz)[o])
        ints.append(np.asarray(pk_int)[o])
        offs.append(offs[-1] + len(o))
    scan = np.repeat(np.arange(n_scans), np.diff(offs))
    mz = np.concatenate(mzs)
    return {"rt": rts, "mz": mz, "int": np.concatenate(ints), "key": scan * 4096.0 + mz,
            "offs": np.asarray(offs)}


def test_traces_finds_charge_and_mono():
    peps = [(500.2500, 2, 60, 4, 1e6), (650.3300, 3, 120, 5, 5e5), (720.4000, 1, 90, 4, 8e5)]
    ms1 = _synthetic_ms1(peps)
    recs = pd.DataFrame(run_traces(ms1, {"min_corr": 0.5}, verbose=False))
    for mono, z, apex, _, _ in peps:
        hit = recs[(recs.charge == z) & (np.abs(recs.mz - mono) < mono * 5e-6)]
        assert len(hit) >= 1, (mono, z, recs)
        assert abs(hit.iloc[0].rt_model - apex * 0.01) < 0.03


def test_reextract_area_and_bounds():
    ms1 = _synthetic_ms1([(500.25, 2, 60, 4, 1e6)])
    df = pd.DataFrame({"mz": [500.25], "charge": [2], "rt_model": [0.60], "rt_begin": [0.55],
                       "rt_end": [0.65], "n_isotopes": [3], "quality": [1.0]})
    out = _reextract(df, ms1, xic_ppm=10, area_isotopes=2, boundary="walk", boundary_frac=0.05,
                     boundary_gap=1, boundary_window=0.5)
    r = out.iloc[0]
    assert r.rt_begin < 0.55 and r.rt_end > 0.65          # walk widened the range
    assert abs(r.rt - 0.60) < 0.011
    # Gaussian integral of mono + M+1 (in minutes): h * sigma * sqrt(2 pi) * (1 + ratio)
    lam = (500.25 - 1.00728) * 2 * 0.000594
    expected = 1e6 * 0.04 * np.sqrt(2 * np.pi) * (1 + lam)
    assert r.area == pytest.approx(expected, rel=0.1)


def test_assign_ms2_window_and_rt():
    feats = pd.DataFrame({"feature_id": [1, 2, 3], "mz": [500.0, 500.6, 800.0], "charge": [2, 2, 2],
                          "rt": [1.0, 1.0, 1.0], "rt_begin": [0.9, 0.9, 0.9],
                          "rt_end": [1.1, 1.1, 1.1], "area": [10.0, 100.0, 5.0]})
    ms2 = pd.DataFrame({"scan": [10, 11, 12], "rt": [1.0, 1.0, 2.0],
                        "precursor_mz": [500.3, 800.0, 500.3], "precursor_charge": [2, 2, 2],
                        "iso_lower": [0.65] * 3, "iso_upper": [0.65] * 3,
                        "precursor_intensity": [0.0] * 3, "ms1_scan": [9, 9, 9]})
    a = assign_ms2_to_features(feats, ms2)
    s10 = a[a.scan == 10]
    assert list(s10.feature_id) == [2, 1]                 # sorted by area, both co-isolated
    assert list(a[a.scan == 11].feature_id) == [3]
    assert list(a[a.scan == 12].feature_id) == [-1]       # outside every RT range
    b = assign_ms2_to_features(feats, ms2, isolation_window_da=0.4)   # [500.1, 500.5]
    assert list(b[b.scan == 10].feature_id) == [-1]       # 500.0 / 500.502 / 500.6 all outside
    c = assign_ms2_to_features(feats, ms2, rt_tolerance=1.0)
    assert set(c[c.scan == 12].feature_id) == {1, 2}
