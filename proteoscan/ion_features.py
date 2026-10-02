"""Extra fragment-ion evidence per candidate: neutral-loss and doubly charged series.

PEAKS's DL model scores b, b-H2O, b-NH3, b2+, y, y-H2O, y-NH3 and y2+ (its model config);
MS2PIP predicts only b/y 1+. Here we count matched ions and their intensity share for the other
six series from the observed spectrum directly (no prediction), within the fragment tolerance.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .rescore_features import SpectrumStore

MONO = {"G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414, "T": 101.047679,
        "C": 103.009185, "L": 113.084064, "I": 113.084064, "N": 114.042927, "D": 115.026943, "Q": 128.058578,
        "K": 128.094963, "E": 129.042593, "M": 131.040485, "H": 137.058912, "F": 147.068414, "R": 156.101111,
        "Y": 163.063329, "W": 186.079313, "U": 150.953636, "O": 237.147727, "X": 0.0}
H2O, NH3, PROTON = 18.010565, 17.026549, 1.007276466
_TOK = re.compile(r"([A-Z])(\[[^\]]*\])?")
_NTERM = re.compile(r"^\[([^\]]*)\]-")
_MODVAL = re.compile(r"[+-]?\d+\.?\d*")


def residue_masses(proforma: str) -> np.ndarray:
    """Residue masses including modification deltas (mass-only or a few named mods)."""
    named = {"Carbamidomethyl": 57.021464, "Oxidation": 15.994915, "Deamidated": 0.984016, "Acetyl": 42.010565,
             "Propionamide": 71.037114, "Gln->pyro-Glu": -17.026549, "Glu->pyro-Glu": -18.010565, "Carbamyl": 43.005814}
    def delta(txt):
        t = txt.strip("[]")
        if t in named:
            return named[t]
        m = _MODVAL.search(t)
        return float(m.group(0)) if m else 0.0
    nterm = 0.0
    m = _NTERM.match(proforma)
    if m:
        nterm = delta(m.group(1)); proforma = proforma[m.end():]
    masses = []
    for aa, mod in _TOK.findall(proforma):
        masses.append(MONO.get(aa, 0.0) + (delta(mod) if mod else 0.0))
    if masses:
        masses[0] += nterm
    return np.asarray(masses)


def series_mz(res: np.ndarray, charge: int) -> dict[str, np.ndarray]:
    b = np.cumsum(res)[:-1] + PROTON
    y = np.cumsum(res[::-1])[:-1] + H2O + PROTON
    out = {"b-H2O": b - H2O, "b-NH3": b - NH3, "y-H2O": y - H2O, "y-NH3": y - NH3}
    if charge >= 2:
        out["b2+"] = (b + PROTON) / 2
        out["y2+"] = (y + PROTON) / 2
    return out


def ion_features(cands: pd.DataFrame, mzml: str, tol: float = 0.5) -> pd.DataFrame:
    """Per candidate: n_loss_matched, loss_intensity_frac, n_2plus_matched, twoplus_intensity_frac,
    frac_explained (all series incl. b/y) — aligned with cands.index."""
    store = SpectrumStore(mzml)
    n = len(cands)
    nl = np.zeros(n); li = np.zeros(n); n2 = np.zeros(n); i2 = np.zeros(n); expl = np.zeros(n)
    scans = cands["scan"].to_numpy(); order = np.argsort(scans, kind="stable")
    cache_scan, cache = None, None
    peps = cands["peptide_sage"].to_numpy(); zs = cands["charge"].to_numpy()
    for i in order:
        sc = int(scans[i])
        if sc != cache_scan:
            cache_scan, cache = sc, store.peaks(sc)
        if cache is None or len(cache[0]) == 0:
            continue
        mz, it = cache; tic = it.sum() or 1.0
        res = residue_masses(peps[i])
        if len(res) < 2:
            continue
        ser = series_mz(res, int(zs[i]))
        b = np.cumsum(res)[:-1] + PROTON; y = np.cumsum(res[::-1])[:-1] + H2O + PROTON
        matched_idx = set()
        def match(theo):
            lo = np.searchsorted(mz, theo - tol, "left"); hi = np.searchsorted(mz, theo + tol, "right")
            cnt = 0; inten = 0.0
            for a, bb in zip(lo, hi):
                if bb > a:
                    k = a + int(np.argmax(it[a:bb])); matched_idx.add(k); cnt += 1; inten += it[k]
            return cnt, inten
        c1 = i1 = c2 = ii2 = 0.0
        for name, theo in ser.items():
            c, inten = match(theo)
            if name.endswith("2+"):
                c2 += c; ii2 += inten
            else:
                c1 += c; i1 += inten
        match(b); match(y)
        nl[i] = c1; li[i] = i1 / tic; n2[i] = c2; i2[i] = ii2 / tic
        expl[i] = sum(it[k] for k in matched_idx) / tic
    return pd.DataFrame({"n_loss_matched": nl, "loss_intensity_frac": li, "n_2plus_matched": n2,
                         "twoplus_intensity_frac": i2, "frac_explained": expl}, index=cands.index)


ION_FEATURE_COLUMNS = ["n_loss_matched", "loss_intensity_frac", "n_2plus_matched", "twoplus_intensity_frac", "frac_explained"]
