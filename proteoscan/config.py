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


POSITIONS = {"anywhere", "peptide_nterm", "protein_nterm", "peptide_cterm"}
# (lowest, highest) isotope error -> Comet `isotope_error` option
ISOTOPE_ERROR_OPTIONS = {(0, 0): 0, (0, 1): 1, (0, 2): 2, (0, 3): 3, (-1, 3): 4, (-1, 1): 5, (-3, 3): 6}


@dataclass
class Modification:
    name: str            # PEAKS name, e.g. "Oxidation"
    residues: str        # residues, e.g. "M" or "NQ"; "N-term" / "C-term" for termini
    delta: float | None = None  # monoisotopic mass delta; looked up from MOD_DELTAS if None
    position: str = "anywhere"  # anywhere | peptide_nterm | protein_nterm | peptide_cterm

    def __post_init__(self) -> None:
        if self.position not in POSITIONS:
            raise ValueError(f"modification {self.name!r}: position {self.position!r} not in {sorted(POSITIONS)}")
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
    # PEAKS treats removal of the initiator Met as an enzymatic N-terminus (docs/PROTEIN_INFERENCE.md), so Comet
    # also searches protein sequences without it by default
    clip_nterm_methionine: bool = True
    precursor_mass_offset_ppm: float = 0.0    # measured precursor m/z error; searched m/z = m/z / (1 + ppm * 1e-6)
    # Precursor mass calibration when no offset is given. "auto" (default): a quick first search measures the run's
    # precursor error (median and spread of confident PSMs); the real search applies the offset and narrows the
    # precursor tolerance to max(5 ppm, 10 robust SD), never wider than parent_mass_error_tolerance_ppm (as
    # MaxQuant and MSFragger narrow the tolerance after calibration; it stops a 13C isotope peak from being reported
    # as a deamidated peptide, 19 mDa away). "offset": apply the offset, keep the configured tolerance.
    # "off": search the measured m/z with the configured tolerance.
    precursor_recalibration: str = "auto"
    # which MS1 features become precursor entries of an MS2 scan: "all" = every feature whose monoisotopic m/z is in
    # the isolation window (default; the rescorer is told which entry is the instrument's own precursor),
    # "primary" = only the instrument's precursor and its isotopes (fewer second identifications per scan)
    feature_linking: str = "all"
    # validation only: proteins with this accession prefix are entrapment twins of the target of the same name;
    # the second pass searches the twin of every accepted protein and the twin inherits its protein support
    entrapment_prefix: str = ""
    # m/z Comet searches for each precursor entry: "entry" = the entry's own m/z (feature-based or monoisotope-corrected,
    # offset included; the intended behaviour, from 2026-11-09); "trailer" = the instrument's trailer monoisotopic m/z
    # for every entry of a scan (what every search before 2026-11-09 did; offset applied to it)
    comet_precursor_mz_source: str = "entry"

    def __post_init__(self) -> None:
        if self.comet_precursor_mz_source not in ("entry", "trailer"):
            raise ValueError(f"comet_precursor_mz_source {self.comet_precursor_mz_source!r} not in ('entry', 'trailer')")
        if isinstance(self.precursor_recalibration, bool):   # YAML reads a bare `off` / `on` as a boolean
            self.precursor_recalibration = "auto" if self.precursor_recalibration else "off"
        if self.precursor_recalibration not in ("auto", "offset", "off"):
            raise ValueError(f"precursor_recalibration {self.precursor_recalibration!r} not in ('auto', 'offset', 'off')")
        if self.feature_linking not in ("all", "primary"):
            raise ValueError(f"feature_linking {self.feature_linking!r} not in ('all', 'primary')")
        if tuple(self.isotope_errors) not in ISOTOPE_ERROR_OPTIONS:
            raise ValueError(f"isotope_errors {tuple(self.isotope_errors)} not supported by Comet; choose one of "
                             f"{sorted(ISOTOPE_ERROR_OPTIONS)}")
        # A variable mod on a residue that also carries a fixed mod is stacked by Comet on top of the
        # static delta (e.g. fixed CAM + variable CAM on C searches C+114), so refuse the combination.
        def sites(m: Modification) -> set[str]:
            return {m.residues} if m.residues in ("N-term", "C-term") else set(m.residues)
        fixed = {s: m.name for m in self.fixed_modifications for s in sites(m)}
        for m in self.variable_modifications:
            clash = sorted(sites(m) & set(fixed))
            if clash:
                raise ValueError(f"variable modification {m.name!r} on {clash} overlaps fixed modification(s) "
                                 f"{sorted({fixed[c] for c in clash})}; Comet would stack them. Make the residue "
                                 "either fixed or variable.")

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
