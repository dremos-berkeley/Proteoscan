"""Monoisotopic precursor correction for MS2 scans without a detected MS1 feature
("featureless" entries in PEAKS terms; see docs/PEAKS_PRECURSOR_RULES.md, rule 3).

PEAKS takes the master-MS1 centroid nearest the isolation target and steps the precursor
DOWN by 1.00335/z while the isotope pattern supports it (never up). We do the same with an
averagine test: for each candidate shift k in {0,-1,-2,-3,-4} the implied monoisotopic peak
is at target - k*C13/z; we score how well the observed intensities at mono, mono+1, mono+2, ...
fit the averagine isotope envelope for that neutral mass, and take the best-supported shift.
"""
from __future__ import annotations

from dataclasses import dataclass

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


def averagine_envelope(neutral_mass: float, n_peaks: int = 6) -> np.ndarray:
    """Relative isotope intensities (M, M+1, ...) for an averagine peptide of this mass."""
    scale = max(neutral_mass, 50.0) / _AVG_MASS
    dist = np.zeros(n_peaks)
    dist[0] = 1.0
    for el, count in _AVG.items():
        n = int(round(count * scale))
        p = np.array(_ISO[el])
        # convolve n copies of the element distribution (binomial-ish via repeated FFT-free convolution)
        el_dist = np.zeros(n_peaks); el_dist[0] = 1.0
        # fast power by repeated squaring
        base = np.zeros(n_peaks); base[:3] = p[:3]
        k = n
        while k:
            if k & 1:
                el_dist = np.convolve(el_dist, base)[:n_peaks]
            base = np.convolve(base, base)[:n_peaks]
            k >>= 1
        dist = np.convolve(dist, el_dist)[:n_peaks]
    return dist / dist.max()


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


def correct_precursor(ms1_mz: np.ndarray, ms1_int: np.ndarray, target_mz: float, charge: int,
                      max_shift: int = 4, tol_ppm: float = 10.0, min_ratio: float = 0.3) -> Corrected:
    """Choose the monoisotopic m/z for an isolation at target_mz with the given charge."""
    z = max(int(charge), 1)
    tol = max(target_mz * tol_ppm * 1e-6, 0.005)
    i0 = _nearest(ms1_mz, ms1_int, target_mz, tol)
    if i0 < 0:
        return Corrected(target_mz, z, 0, "target", 0.0)
    anchor_mz = float(ms1_mz[i0])
    # observed intensities on the isotope grid around the anchor: index j = -max_shift .. +3
    grid = np.arange(-max_shift, 4)
    obs = np.zeros(len(grid))
    present = np.zeros(len(grid), dtype=bool)
    for gi, j in enumerate(grid):
        idx = _nearest(ms1_mz, ms1_int, anchor_mz + j * C13 / z, tol)
        if idx >= 0:
            obs[gi] = ms1_int[idx]
            present[gi] = True
    best = Corrected(anchor_mz, z, 0, "centroid", -np.inf)
    for k in range(0, -max_shift - 1, -1):
        gi_mono = int(np.where(grid == k)[0][0])
        if k < 0 and not present[gi_mono]:
            # allow one missing mono peak only for heavy precursors (arithmetic fallback)
            if not (k == best.shift - 1 and (anchor_mz * z) > 2000):
                continue
        # every isotope between mono and anchor must be present
        if k < 0 and not present[gi_mono:int(np.where(grid == 0)[0][0])].all():
            continue
        mono_mz = anchor_mz + k * C13 / z
        neutral = mono_mz * z - z * PROTON
        env = averagine_envelope(neutral, n_peaks=len(grid) - gi_mono)
        o = obs[gi_mono:]
        if o.max() <= 0:
            continue
        env = env / env.max()
        o = o / o.max()
        # cosine similarity of envelopes over the peaks up to +3 beyond the anchor
        score = float(np.dot(env, o) / (np.linalg.norm(env) * np.linalg.norm(o) + 1e-12))
        # require the candidate mono peak to be a plausible fraction of the anchor intensity
        anchor_gi = int(np.where(grid == 0)[0][0])
        ratio = obs[gi_mono] / (obs[anchor_gi] + 1e-9) if present[gi_mono] else 0.0
        expected = env[0] / env[-k] if -k < len(env) and env[-k] > 0 else 0.0
        if k < 0 and ratio < min_ratio * expected:
            continue
        if score > best.score + (0.0 if k == 0 else 0.02):  # prefer no shift unless clearly better
            src = "centroid" if present[gi_mono] else "arithmetic"
            best = Corrected(float(ms1_mz[_nearest(ms1_mz, ms1_int, mono_mz, tol)]) if present[gi_mono] else mono_mz,
                             z, k, src, score)
    return best
