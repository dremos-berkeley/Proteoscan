"""Stage 4: database search with Sage, driven by a SearchConfig.

Produces a canonical PSM table (pandas) with the columns every later stage relies on:
  source_file, scan, charge, exp_mz, calc_mass, exp_mass, rt, peptide_sage (Sage notation),
  stripped, proteins (';'-joined), is_decoy, rank, isotope_error, precursor_ppm, plus all
  Sage score columns (hyperscore, delta_next, delta_best, matched_peaks, longest_y, ...).
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from .config import SearchConfig

ENZYMES = {
    "trypsin": {"cleave_at": "KR", "restrict": "P", "c_terminal": True},
    "trypsin/p": {"cleave_at": "KR", "restrict": None, "c_terminal": True},
    "lys-c": {"cleave_at": "K", "restrict": "P", "c_terminal": True},
    "chymotrypsin": {"cleave_at": "FWYL", "restrict": "P", "c_terminal": True},
    "glu-c": {"cleave_at": "DE", "restrict": None, "c_terminal": True},
    "asp-n": {"cleave_at": "D", "restrict": None, "c_terminal": False},
}

SCAN_RE = re.compile(r"scan=(\d+)")


def _mods_for_sage(mods) -> dict:
    out: dict[str, list[float] | float] = {}
    for m in mods:
        res = m.residues
        if res.lower() in ("n-term", "nterm"):
            key = "^" if m.position != "protein_nterm" else "["
            out.setdefault(key, []).append(m.delta)
            continue
        if res.lower() in ("c-term", "cterm"):
            out.setdefault("$", []).append(m.delta)
            continue
        for r in res:
            out.setdefault(r, []).append(m.delta)
    return out


def build_sage_config(cfg: SearchConfig, mzmls: list[str | Path], out_dir: str | Path,
                      predict_rt: bool = True, chimera: bool = False, report_psms: int = 1,
                      decoy_tag: str = "rev_") -> dict:
    enz = ENZYMES[cfg.enzyme.lower()]
    static = {k: v[0] for k, v in _mods_for_sage(cfg.fixed_modifications).items()}
    variable = _mods_for_sage(cfg.variable_modifications)
    tol = cfg.parent_mass_error_tolerance_ppm
    ftol = cfg.fragment_mass_error_tolerance_da
    return {
        "database": {
            "bucket_size": 8192 if cfg.semi_enzymatic else 32768,
            "enzyme": {
                "missed_cleavages": cfg.max_missed_cleavages,
                "min_len": cfg.peptide_length_min,
                "max_len": cfg.peptide_length_max,
                "cleave_at": enz["cleave_at"],
                "restrict": enz["restrict"],
                "c_terminal": enz["c_terminal"],
                "semi_enzymatic": cfg.semi_enzymatic,
            },
            "fragment_min_mz": 100.0,
            "fragment_max_mz": 2000.0,
            "peptide_min_mass": 400.0,
            "peptide_max_mass": 6000.0,
            "ion_kinds": ["b", "y"],
            "min_ion_index": 1,
            "static_mods": static,
            "variable_mods": variable,
            "max_variable_mods": cfg.max_variable_ptm_per_peptide,
            "decoy_tag": decoy_tag,
            "generate_decoys": True,
            "fasta": str(cfg.database),
        },
        "precursor_tol": {"ppm": [-tol, tol]},
        "fragment_tol": {"da": [-ftol, ftol]},
        "precursor_charge": list(cfg.precursor_charges),
        "isotope_errors": list(cfg.isotope_errors),
        "deisotope": False,
        "chimera": chimera,
        "wide_window": False,
        "predict_rt": predict_rt,
        "min_peaks": 10,
        "max_peaks": 300,
        "min_matched_peaks": 3,
        "max_fragment_charge": 1 if cfg.low_res_fragments else 2,
        "report_psms": report_psms,
        "output_directory": str(out_dir),
        "mzml_paths": [str(p) for p in mzmls],
    }


def run_sage(cfg: SearchConfig, mzmls: list[str | Path], out_dir: str | Path, **kw) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sage_json = build_sage_config(cfg, mzmls, out_dir, **kw)
    (out_dir / "sage.json").write_text(json.dumps(sage_json, indent=2))
    exe = shutil.which("sage")
    if not exe:
        raise FileNotFoundError("sage not on PATH (conda install -c bioconda sage-proteomics)")
    cmd = [exe, str(out_dir / "sage.json"), "--batch-size", "4"]
    env = {"RAYON_NUM_THREADS": str(cfg.threads)}
    import os
    subprocess.run(cmd, check=True, env={**os.environ, **env})
    return out_dir / "results.sage.tsv"


def read_sage(tsv: str | Path, decoy_tag: str = "rev_") -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t")
    df["scan"] = df["scannr"].astype(str).str.extract(SCAN_RE)[0].astype(int)
    df["source_file"] = df["filename"].map(lambda f: Path(f).name)
    df["is_decoy"] = df["label"] == -1
    df["stripped"] = df["peptide"].str.replace(r"\[[^\]]*\]|-", "", regex=True)
    df = df.rename(columns={"peptide": "peptide_sage", "expmass": "exp_mass", "calcmass": "calc_mass"})
    df["exp_mz"] = (df["exp_mass"] + df["charge"] * 1.007276466) / df["charge"]
    return df
