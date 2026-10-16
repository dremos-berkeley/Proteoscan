"""Monoisotopic precursor correction for MS2 scans without a detected MS1 feature
("featureless" entries in PEAKS terms; see docs/PEAKS_PRECURSOR_RULES.md, rule 3).

PEAKS takes the master-MS1 centroid nearest the isolation target and moves the precursor
DOWN by k * 1.00335/z (k = 0..-4) when the isotope pattern says the target is not the
monoisotopic peak; it never moves up.  The rule below was fitted to PEAKS's exported
PEPMASS values (never to its binary) and reproduces PEAKS's k for ~91% of S4 and ~83% of
S1 featureless entries; see scripts/eval_monoiso.py and
/drive_a/drew/proteoscan-work/analysis/monoiso_eval/REPORT.md.

1. Anchor: the centroid nearest the target within ``tol_ppm``.
2. Isotope ladder: for positions j = -max_shift-1 .. +4, the centroid nearest
   anchor + j*C13/z within ``tol_da`` (absolute; PEAKS cuts off sharply at 0.01 Da).
3. Cluster: the contiguous run of ladder peaks that contains the anchor, j from -max_shift
   up to +3.  With fewer than ``min_peaks`` peaks there is no evidence and k = 0.
4. Each hypothesis "the monoisotopic peak is at k" (k = 0 .. -max_shift) is scored by the
   Pearson correlation between observed intensities and the averagine envelope of that
   mono mass, over the window [first cluster peak - 1, last cluster peak + 1]; positions
   without a cluster peak count as 0.  A hypothesis is dropped when the averagine abundance
   of the target isotope is below ``min_ratio`` x the envelope maximum (the instrument would
   not have picked such a small peak).  A mono below the first cluster peak (unobserved,
   reported arithmetically) is allowed when the mono mass is >= ``arith_min_mass``.
5. From k = 0 downwards, a deeper shift replaces the current best only when its score is
   higher by more than ``margin``.
6. m/z: the observed centroid of the chosen isotope if one lies within ``tol_da``
   ("centroid"), else anchor + k*C13/z ("arithmetic").
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

C13 = 1.003355
PROTON = 1.007276466
# averagine: C4.9384 H7.7583 N1.3577 O1.4773 S0.0417 per 111.1254 Da
_AVG = {"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}
_ISO = {  # (mass of element unit, [abundances of M, M+1, M+2])
    "C": [0.9893, 0.0107, 0.0],
    "H": [0.999885, 0.000115, 0.0],
    "N": [0.99632, 0.00368, 0.0],
    "O": [0.99757, 0.00038, 0.00205],
    "S": [0.9493, 0.0076, 0.0429],
}
_AVG_MASS = 111.1254

# cluster / window geometry, in isotope positions relative to the anchor
_CLUSTER_HI = 3      # highest position that can join the cluster
_WIN_BELOW = 1       # scoring window starts this many positions below the first cluster peak
_WIN_ABOVE = 1       # ... and ends this many above the last one
# averagine lookup table
_ENV_STEP = 10.0     # Da
_ENV_MAX = 20000.0
_ENV_PEAKS = 16


def averagine_envelope(neutral_mass: float, n_peaks: int = 6) -> np.ndarray:
    """Relative isotope intensities (M, M+1, ...) for an averagine peptide of this mass."""
    scale = max(neutral_mass, 50.0) / _AVG_MASS
    dist = np.zeros(n_peaks)
    dist[0] = 1.0
    for el, count in _AVG.items():
        n = int(round(count * scale))
        p = np.array(_ISO[el])
        # n-fold convolution of the element distribution by repeated squaring
        el_dist = np.zeros(n_peaks); el_dist[0] = 1.0
        base = np.zeros(n_peaks); base[:3] = p[:3]
        k = n
        while k:
            if k & 1:
                el_dist = np.convolve(el_dist, base)[:n_peaks]
            base = np.convolve(base, base)[:n_peaks]
            k >>= 1
        dist = np.convolve(dist, el_dist)[:n_peaks]
    return dist / dist.max()


@lru_cache(maxsize=1)
def _envelope_table() -> np.ndarray:
    """Max-normalised averagine envelopes on a 10 Da grid (0 .. 20 kDa); built once (~0.3 s)."""
    grid = np.arange(0.0, _ENV_MAX + _ENV_STEP, _ENV_STEP)
    return np.stack([averagine_envelope(m, _ENV_PEAKS) for m in grid])


def _envelope(neutral_mass: float) -> np.ndarray:
    t = _envelope_table()
    return t[min(max(int(round(neutral_mass / _ENV_STEP)), 0), len(t) - 1)]


@dataclass
class Corrected:
    mz: float
    charge: int
    shift: int          # k <= 0; number of isotopes stepped down
    source: str         # centroid | arithmetic | target
    score: float


def _nearest(mzs: np.ndarray, ints: np.ndarray, mz: float, tol_da: float):
    i = np.searchsorted(mzs, mz)
    best, bestd = -1, tol_da
    for j in (i - 1, i):
        if 0 <= j < len(mzs):
            d = abs(mzs[j] - mz)
            if d <= bestd:
                best, bestd = j, d
    return best


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else -1.0


def correct_precursor(ms1_mz: np.ndarray, ms1_int: np.ndarray, target_mz: float, charge: int,
                      max_shift: int = 4, tol_ppm: float = 10.0, min_ratio: float = 0.1,
                      tol_da: float = 0.01, margin: float = 0.03, min_peaks: int = 2,
                      arith_min_mass: float = 0.0) -> Corrected:
    """Choose the monoisotopic m/z for an isolation at target_mz with the given charge.

    ms1_mz/ms1_int: centroids of the master MS1, m/z ascending.
    max_shift       deepest isotope step considered (k >= -max_shift)
    tol_ppm         anchor search, target -> nearest centroid
    min_ratio       minimum averagine abundance of the target isotope (relative to the
                    envelope maximum) for a hypothesis to be considered
    tol_da          isotope-ladder matching tolerance (absolute m/z)
    margin          Pearson-score margin by which a deeper shift must win
    min_peaks       minimum cluster size (anchor included) before any shift is considered
    arith_min_mass  minimum neutral mono mass at which an unobserved mono may be chosen
    """
    z = max(int(charge), 1)
    ms1_mz = np.asarray(ms1_mz)
    n = len(ms1_mz)
    i0 = _nearest(ms1_mz, ms1_int, target_mz, target_mz * tol_ppm * 1e-6)
    if i0 < 0:
        return Corrected(float(target_mz), z, 0, "target", 0.0)
    anchor = float(ms1_mz[i0])
    step = C13 / z

    # isotope ladder around the anchor
    lo_j = -max_shift - _WIN_BELOW
    js = np.arange(lo_j, _CLUSTER_HI + _WIN_ABOVE + 1)
    a0 = -lo_j                                   # index of the anchor in js
    x = anchor + js * step
    ii = np.searchsorted(ms1_mz, x)
    left = np.clip(ii - 1, 0, n - 1)
    right = np.clip(ii, 0, n - 1)
    pick = np.where(np.abs(ms1_mz[right] - x) < np.abs(ms1_mz[left] - x), right, left)
    ok = np.abs(ms1_mz[pick] - x) <= tol_da
    ok[a0], pick[a0] = True, i0
    inten = np.where(ok, np.asarray(ms1_int, dtype=float)[pick], 0.0)

    # contiguous cluster containing the anchor
    first = 0
    while first - 1 >= -max_shift and ok[a0 + first - 1]:
        first -= 1
    last = 0
    while last + 1 <= _CLUSTER_HI and ok[a0 + last + 1]:
        last += 1
    if last - first + 1 < min_peaks:
        return Corrected(anchor, z, 0, "centroid", 0.0)

    # scoring window [first - 1, last + 1]; only cluster peaks are observed
    w0, w1 = a0 + first - _WIN_BELOW, a0 + last + _WIN_ABOVE
    obs = np.zeros(w1 - w0 + 1)
    obs[_WIN_BELOW:_WIN_BELOW + last - first + 1] = inten[a0 + first:a0 + last + 1]
    wj = js[w0:w1 + 1]
    neutral0 = (anchor - PROTON) * z
    best_k, best_s = None, -np.inf
    for k in range(0, -max_shift - 1, -1):
        mono_mass = neutral0 + k * C13
        if k < first and mono_mass < arith_min_mass:
            continue
        env = _envelope(mono_mass)
        if env[-k] < min_ratio:
            continue
        idx = wj - k
        mdl = np.where(idx >= 0, env[np.clip(idx, 0, len(env) - 1)], 0.0)
        s = _pearson(obs, mdl)
        if best_k is None or s > best_s + margin:
            best_k, best_s = k, s
    if not best_k:
        return Corrected(anchor, z, 0, "centroid", float(best_s) if best_k == 0 else 0.0)
    gi = a0 + best_k
    if ok[gi]:
        return Corrected(float(ms1_mz[pick[gi]]), z, best_k, "centroid", best_s)
    return Corrected(anchor + best_k * step, z, best_k, "arithmetic", best_s)
