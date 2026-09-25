"""Stage 8: build PEAKS-format tables (psms, peptides) from the rescored PSM table.

Protein tables come from `proteins.infer_proteins`; features / denovoOnly from their stages.
Column semantics follow docs/PEAKS_EXPORT_FORMAT.md.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .config import SearchConfig
from .peaks_export import PeaksExport, canonical_columns
from .seqnorm import canonical, mods_of, strip_mods

PROTON = 1.007276466
C13 = 1.003355


def accession_peaks(protein_id: str) -> str:
    """'sp|P62987|RL40_HUMAN' -> 'P62987|RL40_HUMAN'; decoy tags are kept in front."""
    p = protein_id
    tag = ""
    for t in ("#DECOY#", "rev_"):
        if p.startswith(t):
            tag, p = "#DECOY#", p[len(t):]
    parts = p.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr"):
        p = "|".join(parts[1:3])
    return tag + p


def ptm_names(modseq: str, cfg: SearchConfig) -> tuple[str, str]:
    """PEAKS 'PTM' and 'AScore' strings for a canonical modified sequence."""
    lookup: dict[str, tuple[str, bool]] = {}
    for m in cfg.fixed_modifications:
        lookup[f"{m.delta:+.2f}"] = (m.name, True)
    for m in cfg.variable_modifications:
        lookup[f"{m.delta:+.2f}"] = (f"{m.name} ({m.residues})", False)
    names, ascores = [], []
    for pos, site, delta in mods_of(modseq):
        name, fixed = lookup.get(delta, (f"Unknown ({delta})", False))
        names.append(name)
        ascores.append(f"{site}{pos}:{name}:1000.00")
    return ";".join(names), ";".join(ascores)


def build_psms(df: pd.DataFrame, cfg: SearchConfig, sample_name: str = "Sample 1",
               raw_name_map: dict[str, str] | None = None) -> pd.DataFrame:
    """df: one row per accepted PSM with columns from sage_search/rescore stages:
    source_file, scan, charge, exp_mass, calc_mass, rt, peptide_sage, proteins, isotope_error,
    minus10lgp, and optional ms2_corr, delta_rt, feature_id, feature_area, raw_mz, raw_rt,
    tag_length, caa_pct, positional_confidence, tag."""
    n = len(df)
    modseq = df["peptide_sage"].map(canonical)
    stripped = modseq.map(strip_mods)
    z = df["charge"].astype(int)
    mono_mass = df["exp_mass"] - df["isotope_error"].fillna(0) * C13
    mz = (mono_mass + z * PROTON) / z
    raw_mz = df["raw_mz"] if "raw_mz" in df else (df["exp_mass"] + z * PROTON) / z
    ppm = (mono_mass - df["calc_mass"]) / df["calc_mass"] * 1e6
    ptm = modseq.map(lambda s: ptm_names(s, cfg))
    src = df["source_file"].astype(str)
    if raw_name_map:
        src = src.map(lambda s: raw_name_map.get(s, s))
    else:
        src = src.str.replace(r"\.mzML$", ".raw", regex=True)
    out = pd.DataFrame({
        "Peptide": modseq,
        "-10LgP": df["minus10lgp"].round(2),
        "Mass": df["calc_mass"].round(4),
        "Tag Length": df["tag_length"] if "tag_length" in df else 0,
        "CAA (%)": df["caa_pct"] if "caa_pct" in df else 0.0,
        "Length": stripped.str.len(),
        "Delta RT": df["delta_rt"].round(4) if "delta_rt" in df else np.nan,
        "MS2 Correlation": df["ms2_corr"].round(2) if "ms2_corr" in df else np.nan,
        "ppm": ppm.round(1),
        "m/z": mz.round(5),
        "raw m/z": pd.Series(raw_mz).round(5),
        "z": z,
        "RT": df["rt"].round(4),
        "raw RT": (df["raw_rt"] if "raw_rt" in df else df["rt"]).round(4),
        "Area": df["feature_area"] if "feature_area" in df else np.nan,
        "Scan": df["scan"].astype(int),
        "Source File": src,
        "Feature Id": (df["feature_id"] if "feature_id" in df else 0),
        "Accession": df["proteins"].map(lambda p: ";".join(accession_peaks(x) for x in str(p).split(";"))),
        "PTM": [p[0] for p in ptm],
        "AScore": [p[1] for p in ptm],
        "Positional Confidence": df["positional_confidence"] if "positional_confidence" in df else stripped.map(lambda s: " ".join(["0"] * len(s))),
        "tag(>=0.0%)": df["tag"] if "tag" in df else stripped,
    })
    out["Feature Id"] = out["Feature Id"].fillna(0).astype(int)
    out = out[canonical_columns("psms", (sample_name,))]
    return out.sort_values(["Accession", "Peptide", "Scan"]).reset_index(drop=True)


def build_peptides(psms: pd.DataFrame, sample_name: str = "Sample 1") -> pd.DataFrame:
    """One row per modified peptide: its best PSM, summed Area, feature counts."""
    p = psms.copy()
    p["_area"] = pd.to_numeric(p["Area"], errors="coerce")
    p["_feat"] = (p["Feature Id"].fillna(0).astype(int) > 0).astype(int)
    order = p.sort_values("-10LgP", ascending=False)
    best = order.drop_duplicates("Peptide").set_index("Peptide").drop(columns=["_area", "_feat"])
    agg = p.groupby("Peptide").agg(_area=("_area", "sum"), _nfeat=("_feat", "sum"), _anyarea=("_area", lambda s: s.notna().any()))
    best = best.join(agg)
    out = pd.DataFrame({
        "Peptide": best.index,
        "-10LgP": best["-10LgP"].to_numpy(),
        "Mass": best["Mass"].to_numpy(),
        "Tag Length": best["Tag Length"].to_numpy(),
        "CAA (%)": best["CAA (%)"].to_numpy(),
        "Length": best["Length"].to_numpy(),
        "Delta RT": best["Delta RT"].to_numpy(),
        "MS2 Correlation": best["MS2 Correlation"].to_numpy(),
        "ppm": best["ppm"].to_numpy(),
        "m/z": best["m/z"].to_numpy(),
        "z": best["z"].to_numpy(),
        "RT": best["RT"].to_numpy(),
        f"Area {sample_name}": np.where(best["_anyarea"], best["_area"], np.nan),
        "Scan": best["Scan"].to_numpy(),
        "Source File": best["Source File"].to_numpy(),
        "#Feature": best["_nfeat"].to_numpy(),
        f"#Feature {sample_name}": best["_nfeat"].to_numpy(),
        "Accession": best["Accession"].str.replace(";", ":").to_numpy(),
        "PTM": best["PTM"].to_numpy(),
        "AScore": best["AScore"].to_numpy(),
        "Positional Confidence": best["Positional Confidence"].to_numpy(),
        "tag(>=0.0%)": best["tag(>=0.0%)"].to_numpy(),
        "Found By": "DB Search",
    })
    out = out[canonical_columns("peptides", (sample_name,))]
    return out.sort_values(["Accession", "Peptide"]).reset_index(drop=True)


def assemble_export(psms: pd.DataFrame, peptides: pd.DataFrame, proteins=None, protein_peptides=None,
                    features=None, denovo_only=None, ptm_profile=None, modified_peptides=None,
                    sample_name: str = "Sample 1", params: dict | None = None) -> PeaksExport:
    return PeaksExport(psms=psms, peptides=peptides, proteins=proteins, protein_peptides=protein_peptides,
                       features=features, denovo_only=denovo_only, ptm_profile=ptm_profile,
                       modified_peptides=modified_peptides, params=params or {}, sample_name=sample_name)
