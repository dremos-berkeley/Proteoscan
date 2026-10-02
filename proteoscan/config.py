"""Search configuration, expressed with PEAKS Studio's parameter names.

A config is a YAML file such as configs/peaks_default.yml. Field names deliberately mirror
the "Search parameters" table PEAKS writes into db.summary.html so a PEAKS run can be
copied verbatim.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

# Unimod monoisotopic deltas for the modifications PEAKS names in its exports.
MOD_DELTAS: dict[str, float] = {
    "Carbamidomethylation": 57.021464,
    "Oxidation": 15.994915,
    "Deamidation": 0.984016,
    "Acetylation": 42.010565,
    "Propionamide": 71.037114,
    "Pyro-glu from Q": -17.026549,
    "Pyro-glu from E": -18.010565,
    "Carbamylation": 43.005814,
    "Phosphorylation": 79.966331,
    "Methylation": 14.015650,
}


@dataclass
class Modification:
    name: str            # PEAKS name, e.g. "Oxidation"
    residues: str        # residues, e.g. "M" or "NQ"; "N-term" / "C-term" for termini
    delta: float | None = None  # monoisotopic mass delta; looked up from MOD_DELTAS if None
    position: str = "anywhere"  # anywhere | peptide_nterm | protein_nterm | peptide_cterm

    def __post_init__(self) -> None:
        if self.delta is None:
            if self.name not in MOD_DELTAS:
                raise ValueError(f"unknown modification {self.name!r}; give delta explicitly")
            self.delta = MOD_DELTAS[self.name]


@dataclass
class SearchConfig:
    # --- PEAKS "Search parameters" table -------------------------------------------------
    parent_mass_error_tolerance_ppm: float = 15.0
    fragment_mass_error_tolerance_da: float = 0.5
    precursor_mass_search_type: str = "monoisotopic"
    enzyme: str = "Trypsin"
    max_missed_cleavages: int = 2
    digest_mode: str = "Semi-Specific"       # Specific | Semi-Specific | Unspecific
    peptide_length_min: int = 6
    peptide_length_max: int = 45
    fixed_modifications: list[Modification] = field(
        default_factory=lambda: [Modification("Carbamidomethylation", "C")])
    variable_modifications: list[Modification] = field(
        default_factory=lambda: [Modification("Deamidation", "NQ"), Modification("Oxidation", "M")])
    max_variable_ptm_per_peptide: int = 3
    database: str = ""                        # FASTA path (targets only; decoys generated)
    deep_learning_boost: bool = True
    fdr_estimation: bool = True
    confident_amino_acid_threshold_pct: float = 2.0
    # --- PEAKS "Result filtration parameters" ---------------------------------------------
    peptide_fdr: float = 0.01                 # PEAKS picks the -10lgP threshold at this peptide FDR
    protein_score_min: float = 15.0           # protein -10lgP >= 15
    protein_unique_peptides_min: int = 1
    denovo_alc_min_pct: float = 50.0
    # --- engine details not exposed by PEAKS --------------------------------------------
    precursor_charges: tuple[int, int] = (1, 4)
    isotope_errors: tuple[int, int] = (-1, 3)
    threads: int = 16
    fragment_bin_offset: float = 0.4          # only used by Comet-style low-res scoring

    @classmethod
    def from_yaml(cls, path: str | Path) -> "SearchConfig":
        raw: dict[str, Any] = yaml.safe_load(Path(path).read_text()) or {}
        return cls.from_dict(raw)

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> "SearchConfig":
        kw = dict(raw)
        for key in ("fixed_modifications", "variable_modifications"):
            if key in kw:
                kw[key] = [m if isinstance(m, Modification) else Modification(**m) for m in kw[key]]
        for key in ("precursor_charges", "isotope_errors"):
            if key in kw and kw[key] is not None:
                kw[key] = tuple(kw[key])
        unknown = set(kw) - {f.name for f in dataclasses.fields(cls)}
        if unknown:
            raise ValueError(f"unknown config keys: {sorted(unknown)}")
        return cls(**kw)

    def to_dict(self) -> dict[str, Any]:
        d = dataclasses.asdict(self)
        d["precursor_charges"] = list(self.precursor_charges)
        d["isotope_errors"] = list(self.isotope_errors)
        return d

    def to_yaml(self, path: str | Path) -> None:
        Path(path).write_text(yaml.safe_dump(self.to_dict(), sort_keys=False))

    @property
    def semi_enzymatic(self) -> bool:
        return self.digest_mode.lower().startswith("semi")

    @property
    def low_res_fragments(self) -> bool:
        return self.fragment_mass_error_tolerance_da >= 0.1
