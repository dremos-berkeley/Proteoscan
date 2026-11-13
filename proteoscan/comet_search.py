"""Stage 4 (alternative engine): Comet. Low memory (about 3 GB for a semi-tryptic human search
with 3 variable modifications), so it is the default engine for laptops. Produces the same
canonical PSM table as sage_search.read_sage so the later stages are engine-agnostic:
  source_file, scan, precursor_idx, charge, exp_mass, calc_mass, exp_mz, rt, peptide_sage
  (ProForma-like "PEPTM[+15.9949]IDE"), stripped, proteins, is_decoy, rank, isotope_error,
  precursor_ppm, hyperscore (=xcorr), delta_next (=delta_cn*xcorr), plus comet_* columns.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .config import ISOTOPE_ERROR_OPTIONS, SearchConfig
from .monoiso import C13, PROTON

ENZYME_NUMBERS = {"trypsin": 1, "trypsin/p": 2, "lys-c": 3, "lys-n": 4, "arg-c": 5, "asp-n": 6,
                  "cnbr": 7, "glu-c": 8, "pepsina": 9, "chymotrypsin": 10, "none": 0}
SCAN_RE = re.compile(r"scan=(\d+)")


def _default_params() -> str:
    exe = shutil.which("comet")
    if not exe:
        raise FileNotFoundError("comet not on PATH (conda install -c bioconda comet-ms)")
    with tempfile.TemporaryDirectory(prefix="comet_params_") as tmp:  # portable (no hard-coded /tmp)
        subprocess.run([exe, "-p"], cwd=tmp, check=True, capture_output=True)
        return (Path(tmp) / "comet.params.new").read_text()


AA_NAMES = {"G": "glycine", "A": "alanine", "S": "serine", "P": "proline", "V": "valine", "T": "threonine",
            "C": "cysteine", "L": "leucine", "I": "isoleucine", "N": "asparagine", "D": "aspartic_acid",
            "Q": "glutamine", "K": "lysine", "E": "glutamic_acid", "M": "methionine", "H": "histidine",
            "F": "phenylalanine", "U": "selenocysteine", "R": "arginine", "Y": "tyrosine", "W": "tryptophan",
            "O": "pyrrolysine", "B": "user_amino_acid", "J": "user_amino_acid", "X": "user_amino_acid",
            "Z": "user_amino_acid"}
TERMINI = {"N-term": "add_Nterm_peptide", "C-term": "add_Cterm_peptide"}
# Comet's unused variable-mod slot, as written by `comet -p`
EMPTY_VARIABLE_MOD = "0.0 X 0 3 -1 0 0 0.0"
# Modification.position -> Comet (term_distance, n/c-term flag); flag 0/1/2/3 = protein N/C, peptide N/C
COMET_TERM = {"anywhere": (-1, 0), "protein_nterm": (0, 0), "peptide_nterm": (0, 2), "peptide_cterm": (0, 3)}
TERMINUS_DEFAULT_POSITION = {"N-term": "peptide_nterm", "C-term": "peptide_cterm"}
MAX_VARIABLE_MOD_LINES = 15


def fixed_mod_deltas(cfg: SearchConfig) -> dict[str, float]:
    """Effective fixed (static) modifications: residue letter, "N-term" or "C-term" -> mass delta.

    This is what build_comet_params writes as add_* lines and what read_comet must make explicit
    in ProForma (Comet's output omits static mods). Raises on a residue Comet has no add_ line for.
    """
    out: dict[str, float] = {}
    for m in cfg.fixed_modifications:
        keys = [m.residues] if m.residues in TERMINI else list(m.residues)
        for k in keys:
            if k not in AA_NAMES and k not in TERMINI:
                raise ValueError(f"fixed modification {m.name!r}: residue {k!r} has no Comet static-mod parameter")
            if k in out and abs(out[k] - m.delta) > 1e-9:
                raise ValueError(f"two fixed modifications on {k!r} ({out[k]} and {m.delta})")
            out[k] = float(m.delta)
    return out


def _static_param_key(key: str) -> str:
    return TERMINI[key] if key in TERMINI else f"add_{key}_{AA_NAMES[key]}"


def build_comet_params(cfg: SearchConfig, num_output_lines: int = 5, decoy_prefix: str = "DECOY_",
                       generate_decoys: bool = True) -> str:
    text = _default_params()
    # `comet -p` ships non-zero defaults (add_C_cysteine = 57.021464, variable_mod01 = 15.9949 M):
    # zero every static mod and every variable slot so that only the configured mods are searched
    text = re.sub(r"^(add_\w+) = \S+", r"\1 = 0.0", text, flags=re.M)
    text = re.sub(r"^(variable_mod\d\d) = .*$", rf"\1 = {EMPTY_VARIABLE_MOD}", text, flags=re.M)
    lo, hi = cfg.isotope_errors
    iso = ISOTOPE_ERROR_OPTIONS[(lo, hi)]  # SearchConfig refuses unsupported pairs
    settings = {
        "database_name": str(cfg.database),
        "decoy_search": 1 if generate_decoys else 0,
        "decoy_prefix": decoy_prefix,
        "num_threads": cfg.threads,
        "peptide_mass_tolerance_upper": cfg.parent_mass_error_tolerance_ppm,
        "peptide_mass_tolerance_lower": -cfg.parent_mass_error_tolerance_ppm,
        "peptide_mass_units": 2,
        "precursor_tolerance_type": 1,
        "isotope_error": iso,
        "search_enzyme_number": ENZYME_NUMBERS[cfg.enzyme.lower()],
        "num_enzyme_termini": 1 if cfg.semi_enzymatic else 2,
        "allowed_missed_cleavage": cfg.max_missed_cleavages,
        "max_variable_mods_in_peptide": cfg.max_variable_ptm_per_peptide,
        "fragment_bin_tol": 1.0005 if cfg.low_res_fragments else 0.02,
        "fragment_bin_offset": 0.4 if cfg.low_res_fragments else 0.0,
        "theoretical_fragment_ions": 1 if cfg.low_res_fragments else 0,
        "use_A_ions": 0, "use_B_ions": 1, "use_C_ions": 0, "use_X_ions": 0, "use_Y_ions": 1, "use_Z_ions": 0, "use_NL_ions": 0,
        "output_txtfile": 1, "output_pepxmlfile": 0,
        "num_output_lines": num_output_lines,
        "digest_mass_range": "400.0 6000.0",
        "peptide_length_range": f"{cfg.peptide_length_min} {cfg.peptide_length_max}",
        "precursor_charge": f"{cfg.precursor_charges[0]} {cfg.precursor_charges[1]}",
        "max_precursor_charge": cfg.precursor_charges[1],
        "max_fragment_charge": 1 if cfg.low_res_fragments else 2,
        "minimum_peaks": 10,
        "activation_method": "ALL",
    }
    # variable mods: one line per (delta, position); termini as Comet 'n'/'c' with term_distance and flag
    residues_by: dict[tuple[float, str], str] = {}
    for m in cfg.variable_modifications:
        res = {"N-term": "n", "C-term": "c"}.get(m.residues, m.residues)
        pos = m.position
        if pos == "anywhere" and m.residues in TERMINUS_DEFAULT_POSITION:
            pos = TERMINUS_DEFAULT_POSITION[m.residues]
        residues_by[(m.delta, pos)] = residues_by.get((m.delta, pos), "") + res
    if len(residues_by) > MAX_VARIABLE_MOD_LINES:
        raise ValueError(f"{len(residues_by)} variable-mod lines; Comet supports {MAX_VARIABLE_MOD_LINES}")
    for i, ((delta, pos), res) in enumerate(sorted(residues_by.items()), start=1):
        dist, flag = COMET_TERM[pos]
        settings[f"variable_mod{i:02d}"] = f"{delta:.6f} {res} 0 {cfg.max_variable_ptm_per_peptide} {dist} {flag} 0 0.0"
    settings["clip_nterm_methionine"] = int(cfg.clip_nterm_methionine)
    # fixed (static) mods: every add_* line was zeroed above, so only configured mods are applied
    for key, delta in fixed_mod_deltas(cfg).items():
        settings[_static_param_key(key)] = delta
    missing = []
    for k, v in settings.items():
        text, n = re.subn(rf"^{re.escape(k)} = .*$", f"{k} = {v}", text, flags=re.M)
        if n != 1:
            missing.append(k)
    if missing:
        raise RuntimeError(f"comet.params keys not found for this Comet version: {missing}")
    return text


def run_comet(cfg: SearchConfig, mzmls: list[str | Path], out_dir: str | Path, num_output_lines: int = 5) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    params = out_dir / "comet.params"
    params.write_text(build_comet_params(cfg, num_output_lines))
    outs = []
    for m in mzmls:
        base = Path(m).stem
        out_txt = out_dir / f"{base}.txt"
        if not out_txt.exists():
            subprocess.run([shutil.which("comet"), f"-P{params}", f"-N{out_dir / base}", str(m)], check=True)
        outs.append(out_txt)
    return outs


def _proforma(modified_peptide: str) -> str:
    """'K.PEPTM[15.9949]IDE.R' -> 'PEPTM[+15.9949]IDE'; 'n[42.0106]PEP' -> '[+42.0106]-PEP'."""
    s = modified_peptide
    if len(s) > 4 and s[1] == "." and s[-2] == ".":
        s = s[2:-2]
    s = re.sub(r"\[(-?\d)", lambda m: "[" + ("+" if not m.group(1).startswith("-") else "") + m.group(1), s)
    s = re.sub(r"^n\[([^\]]+)\]", r"[\1]-", s)
    s = re.sub(r"c\[([^\]]+)\]$", r"-[\1]", s)
    return s


_NTERM_RE = re.compile(r"^\[([+-]?[\d.]+)\]-")
_CTERM_RE = re.compile(r"-\[([+-]?[\d.]+)\]$")


def _add_terminal(pep: str, which: str, delta: float) -> str:
    """Add a fixed terminal delta to a ProForma-like peptide; an existing (variable) mod on the same
    terminus is summed with it so the written mass stays correct."""
    rx = _NTERM_RE if which == "N-term" else _CTERM_RE
    m = rx.search(pep)
    total = delta + (float(m.group(1)) if m else 0.0)
    tag = f"[{total:+.4f}]"
    if which == "N-term":
        return (tag + "-" + pep[m.end():]) if m else (tag + "-" + pep)
    return (pep[:m.start()] + "-" + tag) if m else (pep + "-" + tag)


def read_comet(txt: str | Path, source_file: str | None = None, decoy_prefix: str = "DECOY_",
               fixed_mods: dict[str, float] | None = None, clipped_met_enzymatic: bool = False) -> pd.DataFrame:
    df = pd.read_csv(txt, sep="\t", skiprows=1)
    fixed_mods = fixed_mods or {}
    src = source_file or Path(txt).stem
    out = pd.DataFrame({
        "source_file": src,
        "scan": df["scan"].astype(int),
        "precursor_idx": 0,
        "charge": df["charge"].astype(int),
        "exp_mass": df["exp_neutral_mass"].astype(float),
        "calc_mass": df["calc_neutral_mass"].astype(float),
        "rt": df["retention_time_sec"].astype(float) / 60.0,
        "peptide_sage": df["modified_peptide"].map(_proforma),
        "stripped": df["plain_peptide"],
        "proteins": df["protein"].astype(str).str.replace(",", ";"),
        "rank": df["num"].astype(int),
        "comet_evalue": df["e-value"].astype(float),
        "hyperscore": df["xcorr"].astype(float),
        "comet_deltacn": df["delta_cn"].astype(float),
        "comet_sp": df["sp_score"].astype(float),
        "comet_sp_rank": df["sp_rank"].astype(float) if "sp_rank" in df else np.nan,
        "matched_peaks": df["ions_matched"].astype(int),
        "comet_ions_total": df["ions_total"].astype(int),
        "prev_aa": df["prev_aa"], "next_aa": df["next_aa"],
    })
    # fixed mods are implicit in Comet output; make them explicit for ProForma / canonical forms.
    # Keys as returned by fixed_mod_deltas: residue letters, "N-term", "C-term".
    for res, delta in fixed_mods.items():
        if res in TERMINI:
            out["peptide_sage"] = out["peptide_sage"].map(lambda p, r=res, d=delta: _add_terminal(p, r, d))
        else:
            out["peptide_sage"] = out["peptide_sage"].str.replace(res, f"{res}[+{delta:.4f}]", regex=False)
    out["is_decoy"] = out["proteins"].map(lambda p: all(x.startswith(decoy_prefix) for x in p.split(";")))
    out["exp_mz"] = (out["exp_mass"] + out["charge"] * PROTON) / out["charge"]
    diff = out["exp_mass"] - out["calc_mass"]
    out["isotope_error"] = np.rint(diff / C13).astype(int)
    out["precursor_ppm"] = (diff - out["isotope_error"] * C13) / out["calc_mass"] * 1e6
    out["peptide_len"] = out["stripped"].str.len()
    out["missed_cleavages"] = out["stripped"].str[:-1].str.count(r"[KR](?!P)")
    nterm_ok = out["prev_aa"].isin(["K", "R", "-"]) & ~out["stripped"].str.startswith("P")
    if clipped_met_enzymatic:  # fully specific search with clip_nterm_methionine: prev_aa M = initiator Met removed
        nterm_ok |= out["prev_aa"] == "M"
    cterm_ok = out["stripped"].str[-1].isin(["K", "R"]) | (out["next_aa"] == "-")
    out["semi_enzymatic"] = (~(nterm_ok & cterm_ok)).astype(int)
    out["matched_intensity_pct"] = np.nan
    out["longest_b"] = np.nan; out["longest_y"] = np.nan; out["longest_y_pct"] = np.nan
    out["scored_candidates"] = np.nan; out["poisson"] = np.nan; out["fragment_ppm"] = np.nan
    out["delta_next"] = out["comet_deltacn"] * out["hyperscore"]
    best = out.groupby(["source_file", "scan"])["hyperscore"].transform("max")
    out["delta_best"] = best - out["hyperscore"]
    return out


def canonicalize_clipped_nterm(cands: pd.DataFrame, fasta: str | Path,
                               decoy_prefixes: tuple[str, ...] = ("DECOY_", "#DECOY#")) -> pd.DataFrame:
    """Give every peptide that starts at residue 2 of a protein beginning with Met the flank '-' (protein
    N-terminus after initiator-Met removal) and recompute `semi_enzymatic` for it.

    With clip_nterm_methionine Comet reports such a peptide either as the clipped protein N-terminus (prev_aa
    '-') or as a cleavage after Met-1 (prev_aa 'M'): mostly 'M' for targets, about half '-' for decoys, and not
    the same way in two runs of one search. Decoy peptides are the target peptide reversed with the C-terminal
    residue kept, so a decoy is looked up through its target peptide. Modifies and returns `cands`."""
    from .decoys import read_fasta
    rows = np.flatnonzero((cands["prev_aa"] == "M").to_numpy())
    if not len(rows):
        return cands

    def split(protein_id):
        for t in decoy_prefixes:
            if protein_id.startswith(t):
                return protein_id[len(t):], True
        return protein_id, False

    stripped = cands["stripped"].to_numpy()
    proteins = cands["proteins"].to_numpy()
    need = {split(x)[0] for i in rows for x in str(proteins[i]).split(";")}
    start = {}
    for h, seq in read_fasta(fasta):
        pid = h.split()[0]
        if pid in need and seq[:1] == "M":
            start[pid] = seq[1:]
    hit = []
    for i in rows:
        pep = stripped[i]
        for x in str(proteins[i]).split(";"):
            pid, decoy = split(x)
            rest = start.get(pid)
            if rest is not None and rest.startswith(pep[:-1][::-1] + pep[-1] if decoy else pep):
                hit.append(i)
                break
    if hit:
        hit = np.asarray(hit)
        col = cands.columns.get_loc
        cands.iloc[hit, col("prev_aa")] = "-"
        if "semi_enzymatic" in cands.columns:
            sub = cands.iloc[hit]
            nterm_ok = ~sub["stripped"].str.startswith("P")
            cterm_ok = sub["stripped"].str[-1].isin(["K", "R"]) | (sub["next_aa"] == "-")
            cands.iloc[hit, col("semi_enzymatic")] = (~(nterm_ok & cterm_ok)).astype(int).to_numpy()
    return cands

