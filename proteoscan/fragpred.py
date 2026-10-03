"""Fragment-intensity prediction with AlphaPeptDeep (peptdeep) MS2 models, fine-tuned for ion-trap
spectra (docs/IONTRAP_MODEL.md), plus the shared fragment-annotation helpers used to build the
training data and the rescoring features.

Layout convention (alphabase): a peptide of n residues has n-1 fragment rows; row i is the cleavage
after residue i, so b = b_{i+1} and y = y_{n-1-i}. Columns are charged fragment types
("b_z1", "y_z2", "b_H2O_z1", ...).

Checkpoint = <name>.pth (peptdeep state dict) + <name>.json:
  {"charged_frag_types": [...], "nce": 28.0, "instrument": "Lumos", ...}
`PeptdeepPredictor("path/to/model.pth")` reads both.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

PROTON = 1.007276466
H2O = 18.010565
NH3 = 17.026549
MONO = {"G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414, "T": 101.047679,
        "C": 103.009185, "L": 113.084064, "I": 113.084064, "N": 114.042927, "D": 115.026943, "Q": 128.058578,
        "K": 128.094963, "E": 129.042593, "M": 131.040485, "H": 137.058912, "F": 147.068414, "R": 156.101111,
        "Y": 163.063329, "W": 186.079313, "U": 150.953636, "O": 237.147727}

# mass delta -> (alphabase name per residue, delta); residue "^" = peptide N-terminus
_MODS = [
    (57.021464, {"C": "Carbamidomethyl@C", "^": "Carbamidomethyl@Any_N-term"}),
    (15.994915, {"M": "Oxidation@M"}),
    (0.984016, {"N": "Deamidated@N", "Q": "Deamidated@Q"}),
    (42.010565, {"^": "Acetyl@Any_N-term"}),
    (71.037114, {"C": "Propionamide@C", "^": "Propionamide@Any_N-term"}),
    (-17.026549, {"^Q": "Gln->pyro-Glu@Q^Any_N-term"}),
    (-18.010565, {"^E": "Glu->pyro-Glu@E^Any_N-term"}),
    (43.005814, {"^": "Carbamyl@Any_N-term", "K": "Carbamyl@K"}),
]
_NAMED = {"Carbamidomethyl": 57.021464, "Oxidation": 15.994915, "Deamidated": 0.984016, "Acetyl": 42.010565,
          "Propionamide": 71.037114, "Gln->pyro-Glu": -17.026549, "Glu->pyro-Glu": -18.010565, "Carbamyl": 43.005814}
_TOKEN = re.compile(r"([A-Z])((?:\[[^\]]*\]|\([^)]*\))*)")
_NTERM = re.compile(r"^((?:\[[^\]]*\]|\([^)]*\))+)-?")
_NUM = re.compile(r"[+-]?\d+(?:\.\d+)?")

HCD_TYPES = ["b_z1", "b_z2", "y_z1", "y_z2"]
LOSS_TYPES = ["b_H2O_z1", "b_NH3_z1", "y_H2O_z1", "y_NH3_z1"]
DEFAULT_FRAG_TYPES = HCD_TYPES + LOSS_TYPES


def _delta(txt: str) -> float:
    total = 0.0
    for t in re.findall(r"\[([^\]]*)\]|\(([^)]*)\)", txt):
        t = t[0] or t[1]
        if t in _NAMED:
            total += _NAMED[t]
        else:
            m = _NUM.search(t)
            total += float(m.group(0)) if m else 0.0
    return total


def parse_peptide(modseq: str):
    """'PEPC[+57.0215]M(+15.99)K' (mass deltas in [] or (), or Unimod names) ->
    (sequence, residue masses incl. mods, alphabase mods str, mod_sites str, ok).
    ok is False when a modification has no alphabase name (masses are still exact)."""
    s = modseq.strip()
    nterm = 0.0
    m = _NTERM.match(s)
    if m:
        nterm = _snap(_delta(m.group(1))); s = s[m.end():]
    seq, masses, names, sites, ok = [], [], [], [], True
    for i, (aa, mods) in enumerate(_TOKEN.findall(s)):
        seq.append(aa)
        d = _delta(mods) if mods else 0.0
        if mods:
            d = _snap(d)
        masses.append(MONO.get(aa, 0.0) + d)
        if mods:
            nm = _mod_name(aa, d, False)
            if nm is None:
                ok = False
            else:
                names.append(nm); sites.append(str(i + 1))
    if nterm:
        nm = _mod_name(seq[0] if seq else "", nterm, True)
        if nm is None:
            ok = False
        else:
            names.append(nm); sites.append("0")
        if masses:
            masses[0] += nterm
    return "".join(seq), np.asarray(masses, dtype=float), ";".join(names), ";".join(sites), ok


def _snap(d: float) -> float:
    """Replace a rounded delta (e.g. PEAKS's +15.99) by the exact mass of the known modification."""
    for delta, _ in _MODS:
        if abs(d - delta) < 0.01:
            return delta
    return d


def _mod_name(aa: str, d: float, nterm: bool):
    for delta, names in _MODS:
        if abs(d - delta) < 0.01:
            if nterm:
                return names.get("^" + aa) or names.get("^")
            return names.get(aa)
    return None


def fragment_mz(res: np.ndarray, frag_types: list[str]) -> np.ndarray:
    """(n-1, len(frag_types)) theoretical m/z in alphabase row layout."""
    b = np.cumsum(res)[:-1]
    y = np.cumsum(res[::-1])[:-1][::-1] + H2O
    out = np.empty((len(b), len(frag_types)))
    for k, ft in enumerate(frag_types):
        base, z = ft.rsplit("_z", 1)
        z = int(z)
        parts = base.split("_")
        m = (b if parts[0] == "b" else y).copy()
        if len(parts) > 1:
            m -= {"H2O": H2O, "NH3": NH3}[parts[1]]
        out[:, k] = (m + z * PROTON) / z
    return out


def annotate(theo: np.ndarray, mz: np.ndarray, it: np.ndarray, tol: float) -> np.ndarray:
    """Max observed intensity within +/- tol of each theoretical m/z (0 if none); any shape."""
    flat = theo.ravel()
    lo = np.searchsorted(mz, flat - tol, side="left")
    hi = np.searchsorted(mz, flat + tol, side="right")
    out = np.zeros(len(flat))
    for k in np.nonzero(hi > lo)[0]:
        out[k] = it[lo[k]:hi[k]].max()
    return out.reshape(theo.shape)


def valid_mask(theo: np.ndarray, frag_types: list[str], charge: int, lo_mz: float, hi_mz: float) -> np.ndarray:
    """True where an ion is observable: fragment charge <= precursor charge (z2 needs charge >= 2;
    PEAKS/peptdeep convention) and m/z inside the scan window."""
    zs = np.array([int(ft.rsplit("_z", 1)[1]) for ft in frag_types])
    ok = (zs[None, :] <= max(int(charge), 1)) & (theo >= lo_mz) & (theo <= hi_mz)
    return ok


class PeptdeepPredictor:
    """Deterministic, batched MS2 intensity prediction from a (fine-tuned) peptdeep checkpoint.
    ckpt = "released" uses the released generic HCD model."""

    def __init__(self, ckpt: str = "released", device: str | None = None, nce: float | None = None,
                 instrument: str | None = None):
        import torch
        from peptdeep.model.ms2 import pDeepModel
        torch.manual_seed(0)
        self.meta = {"charged_frag_types": ["b_z1", "b_z2", "y_z1", "y_z2"], "nce": 30.0, "instrument": "Lumos"}
        if ckpt != "released":
            mp = Path(ckpt).with_suffix(".json")
            if mp.exists():
                self.meta.update(json.loads(mp.read_text()))
        if nce is not None:
            self.meta["nce"] = float(nce)
        if instrument is not None:
            self.meta["instrument"] = instrument
        dev = device or ("gpu" if torch.cuda.is_available() else "cpu")
        self.frag_types = list(self.meta["charged_frag_types"])
        if ckpt == "released":
            from peptdeep.pretrained_models import ModelManager
            mm = ModelManager(device=dev)
            mm.load_installed_models()
            self.model = mm.ms2_model
            self.model.charged_frag_types = self.frag_types
        else:
            self.model = pDeepModel(charged_frag_types=self.frag_types, device=dev)
            self.model.load(str(ckpt))
        self.model.model.eval()
        try:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
        except Exception:
            pass
        self.ckpt = ckpt

    def predict(self, sequences, mods, mod_sites, charges, nce=None, instrument=None,
                batch_size: int = 512) -> list[np.ndarray]:
        """List of (n-1, len(frag_types)) arrays (max-normalised, linear) in input order."""
        df = pd.DataFrame({"sequence": list(sequences), "mods": list(mods), "mod_sites": list(mod_sites),
                           "charge": np.asarray(charges, dtype=int)})
        df["nAA"] = df["sequence"].str.len()
        df["nce"] = float(nce if nce is not None else self.meta["nce"]) if np.isscalar(nce) or nce is None else np.asarray(nce, float)
        df["instrument"] = instrument or self.meta["instrument"]
        df["orig"] = np.arange(len(df))
        df = df.sort_values(["nAA", "orig"], kind="stable").reset_index(drop=True)
        pred = self.model.predict(df, batch_size=batch_size, verbose=False)
        vals = pred[self.frag_types].to_numpy(dtype=np.float32)
        st = df["frag_start_idx"].to_numpy(); sp = df["frag_stop_idx"].to_numpy()
        out = [None] * len(df)
        for k, o in enumerate(df["orig"].to_numpy()):
            out[o] = vals[st[k]:sp[k]]
        return out
