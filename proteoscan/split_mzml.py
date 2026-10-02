"""Split an mzML by MS2 scan type (analyzer and activation), e.g. ProteomeTools runs that
interleave CID-IT, HCD-IT and HCD-OT scans of the same precursor. MS1 scans are kept in every
output. Output names: <stem>.<type>.mzML with type in {hcd_it, cid_it, hcd_ft, ...}."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pyopenms as oms

_FILTER = re.compile(r"^(FTMS|ITMS).*?@(hcd|cid|etd|ethcd|etcid)([0-9.]+)", re.I)


def scan_type(spec: oms.MSSpectrum) -> str | None:
    filt = None
    if spec.metaValueExists("filter string"):
        filt = str(spec.getMetaValue("filter string"))
    elif spec.getAcquisitionInfo().size():
        a = spec.getAcquisitionInfo()[0]
        if a.metaValueExists("filter string"):
            filt = str(a.getMetaValue("filter string"))
    if filt:
        m = _FILTER.search(filt)
        if m:
            analyzer = "ft" if m.group(1).upper() == "FTMS" else "it"
            return f"{m.group(2).lower()}_{analyzer}_{m.group(3).rstrip('0').rstrip('.')}"
    # fallback: instrument settings
    try:
        act = spec.getPrecursors()[0].getActivationMethods()
        inst = spec.getInstrumentSettings()
    except Exception:
        return None
    return None


def split(mzml: str | Path, out_dir: str | Path, by_energy: bool = False) -> dict[str, Path]:
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(mzml), exp)
    groups: dict[str, oms.MSExperiment] = {}
    ms1 = []
    for i in range(exp.getNrSpectra()):
        s = exp[i]
        if s.getMSLevel() == 1:
            ms1.append(s); continue
        t = scan_type(s)
        if t is None:
            continue
        if not by_energy:
            t = "_".join(t.split("_")[:2])
        groups.setdefault(t, oms.MSExperiment())
        groups[t].addSpectrum(s)
    out = {}
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    for t, g in groups.items():
        full = oms.MSExperiment()
        for s in ms1:
            full.addSpectrum(s)
        for i in range(g.getNrSpectra()):
            full.addSpectrum(g[i])
        full.sortSpectra(True)
        p = out_dir / f"{Path(mzml).stem}.{t}.mzML"
        oms.MzMLFile().store(str(p), full)
        out[t] = p
        print(f"{t}: {g.getNrSpectra()} MS2 -> {p.name}")
    return out


if __name__ == "__main__":
    split(sys.argv[1], sys.argv[2], by_energy="--by-energy" in sys.argv)
