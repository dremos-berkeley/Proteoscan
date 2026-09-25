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
    df["precursor_idx"] = df["scannr"].astype(str).str.extract(r"precursor=(\d+)")[0].fillna(0).astype(int)
    df["source_file"] = df["filename"].map(lambda f: Path(f).name)
    df["is_decoy"] = df["label"] == -1
    df["stripped"] = df["peptide"].str.replace(r"\[[^\]]*\]|-", "", regex=True)
    df = df.rename(columns={"peptide": "peptide_sage", "expmass": "exp_mass", "calcmass": "calc_mass"})
    df["exp_mz"] = (df["exp_mass"] + df["charge"] * 1.007276466) / df["charge"]
    return df


# ----------------------------------------------------------------------------------------------
# Memory-bounded chunked search.
#
# Sage keeps the whole fragment index in RAM. A semi-tryptic human search with 3 variable mods
# needs several hundred GB, far beyond a laptop. We split the FASTA round-robin into chunks,
# search each chunk with report_psms candidates per spectrum, and merge candidates by hyperscore.
# Hyperscore does not depend on the database, so merging is exact for the candidate list; the
# rank-dependent columns (delta_next, delta_best, rank) and Sage's own discriminant/q-values are
# recomputed or dropped here and replaced by our rescoring stage.
# ----------------------------------------------------------------------------------------------
import dataclasses as _dc

from .decoys import split_fasta

SAGE_RANK_COLS = ["sage_discriminant_score", "posterior_error", "spectrum_q", "peptide_q", "protein_q"]


def run_sage_chunked(cfg: SearchConfig, mzmls: list[str | Path], out_dir: str | Path,
                     n_chunks: int = 16, report_psms: int = 5, **kw) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks = split_fasta(cfg.database, out_dir / "fasta_chunks", n_chunks)
    merged_parts = []
    for i, chunk in enumerate(chunks):
        cdir = out_dir / f"chunk_{i:03d}"
        tsv = cdir / "results.sage.tsv"
        if not tsv.exists():
            ccfg = _dc.replace(cfg, database=str(chunk))
            run_sage(ccfg, mzmls, cdir, report_psms=report_psms, **kw)
        part = pd.read_csv(tsv, sep="\t")
        part["chunk"] = i
        merged_parts.append(part)
    merged = pd.concat(merged_parts, ignore_index=True)
    merged = merge_candidates(merged, report_psms)
    out = out_dir / "results.sage.tsv"
    merged.to_csv(out, sep="\t", index=False)
    return out


def merge_candidates(df: pd.DataFrame, keep: int) -> pd.DataFrame:
    """Keep the top `keep` candidates per (filename, scannr) by hyperscore and recompute
    rank/delta columns. Identical peptides found in two chunks (cannot happen with a round-robin
    split of distinct proteins, but shared sequences across proteins can) are collapsed, merging
    their protein lists."""
    df = df.copy()
    key = ["filename", "scannr", "peptide", "charge"]  # scannr carries the precursor index
    agg = {c: "first" for c in df.columns if c not in key}
    agg["proteins"] = lambda s: ";".join(sorted(set(";".join(s).split(";"))))
    agg["num_proteins"] = "sum"
    df = df.groupby(key, as_index=False, sort=False).agg(agg)
    df["label"] = df["proteins"].map(lambda p: -1 if all(x.startswith("rev_") or x.startswith("#DECOY#") for x in p.split(";")) else 1)
    df = df.sort_values(["filename", "scannr", "hyperscore"], ascending=[True, True, False])
    df["rank"] = df.groupby(["filename", "scannr"]).cumcount() + 1
    df = df[df["rank"] <= keep].copy()
    best = df.groupby(["filename", "scannr"])["hyperscore"].transform("max")
    nxt = df.groupby(["filename", "scannr"])["hyperscore"].shift(-1)
    df["delta_best"] = best - df["hyperscore"]
    df["delta_next"] = (df["hyperscore"] - nxt).fillna(df["hyperscore"])
    df = df.drop(columns=[c for c in SAGE_RANK_COLS if c in df.columns])
    return df.reset_index(drop=True)
