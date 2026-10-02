"""Stage 1: Thermo .raw -> centroided mzML via ThermoRawFileParser (mono on macOS)."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def find_trfp() -> list[str]:
    for name in ("ThermoRawFileParser", "thermorawfileparser"):
        exe = shutil.which(name)
        if exe:
            return [exe]
    raise FileNotFoundError("ThermoRawFileParser not on PATH (conda install -c bioconda thermorawfileparser)")


def raw_to_mzml(raw: str | Path, out_dir: str | Path, overwrite: bool = False) -> Path:
    raw = Path(raw)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / (raw.stem + ".mzML")
    if out.exists() and out.stat().st_size > 0 and not overwrite:
        return out
    cmd = find_trfp() + ["-i", str(raw), "-o", str(out_dir), "-f", "2"]
    subprocess.run(cmd, check=True)
    if not out.exists():
        raise RuntimeError(f"conversion produced no file: {out}")
    return out
