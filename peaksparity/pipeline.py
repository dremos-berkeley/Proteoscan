"""End-to-end driver: mzML (or Thermo raw) -> PEAKS-format export directory.

    python -m peaksparity.pipeline --config configs/peaks_facility_human.yml \
        --input run.mzML --out work/run1 [--sample-name "Sample 1"] [--ref <PEAKS export dir>]
        [--mode tda|imitate --model model.pkl] [--features features.parquet]

Stages (docs/ARCHITECTURE.md):
  1 convert      raw -> mzML (skipped for mzML input)
  2 precursors   MS1 features + PEAKS-style precursor list -> expanded mzML
  3 search       Comet (low memory) on the expanded mzML
  4 features     fragment-intensity / retention-time / ion-series evidence per candidate
  5 rescore      target-decoy (mokapot) or PEAKS-imitation model
  6 filter       one peptide per (scan, precursor); threshold at the peptide-sequence FDR
  7 export       PSMs (I/L variants split), peptides, proteins, PTM tables in PEAKS CSV format
  8 concordance  optional comparison with a PEAKS reference export
Every stage caches its output under --out and is skipped when the cache exists.
"""
from __future__ import annotations

import argparse
import json
import logging
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .comet_search import read_comet, run_comet
from .config import SearchConfig
from .export_builder import apply_protein_mapping, assemble_export, build_peptides, build_psms
from .fdr import threshold_at_fdr
from .monoiso import PROTON
from .peaks_export import write_export
from .rescore import apply_imitation, best_per_spectrum, minus10lgp, rescore_tda
from .rescore_features import compute_features

log = logging.getLogger("peaksparity.pipeline")


def stage_precursors(mzml: Path, out: Path, features_path: str | None, cfg: SearchConfig):
    from .precursors import build_precursor_list, load_experiment, ms2_metadata, write_expanded_mzml
    ent_pq = out / "precursors.parquet"
    exp_mzml = out / (mzml.stem + ".expanded.mzML")
    if ent_pq.exists() and exp_mzml.exists():
        return pd.read_parquet(ent_pq), exp_mzml
    exp = load_experiment(mzml)
    ms2 = ms2_metadata(exp)
    feats = None
    if features_path == "none":
        feats = None
    elif features_path:
        feats = pd.read_parquet(features_path)
    else:
        from .features import detect_features, load_config
        fcfg = Path(__file__).resolve().parents[1] / "configs" / "features_lumos_120k.yml"
        kw = {k: v for k, v in load_config(fcfg).items() if not k.startswith("ms2_")} if fcfg.exists() else {}
        feats = detect_features(str(mzml), **kw)
        feats.to_parquet(out / "ms1_features.parquet")
    entries = build_precursor_list(exp, feats, ms2)
    entries = entries.merge(ms2[["scan", "rt_min", "precursor_mz", "target_mz"]], on="scan", how="left")
    entries = write_expanded_mzml(exp, entries, exp_mzml)
    entries.to_parquet(ent_pq)
    return entries, exp_mzml


def attach_entries(cands: pd.DataFrame, entries: pd.DataFrame) -> pd.DataFrame:
    """Map search results on the expanded mzML (keyed by query_id in the scan column) back to
    the raw scan and its precursor entry."""
    e = entries.set_index("query_id")
    c = cands.copy()
    q = c["scan"].to_numpy()
    c["query_id"] = q
    c["scan"] = e["scan"].reindex(q).to_numpy().astype(int)
    c["precursor_idx"] = e["precursor_idx"].reindex(q).to_numpy().astype(int)
    c["feature_id"] = e["feature_id"].reindex(q).fillna(0).to_numpy().astype(int)
    c["feature_area"] = e["area"].reindex(q).to_numpy()
    c["entry_rt"] = e["rt"].reindex(q).to_numpy()
    c["entry_source"] = e["source"].reindex(q).to_numpy()
    c["raw_rt"] = (e["rt_min"].reindex(q).to_numpy() if "rt_min" in e else c["rt"].to_numpy())
    c["raw_mz"] = (e["precursor_mz"].reindex(q).to_numpy() if "precursor_mz" in e else np.nan)
    # rank/delta columns are per query already; spectrum grouping must use scan + precursor_idx
    return c


def peaks_filter(top: pd.DataFrame, cfg: SearchConfig, psm_floor: float = 20.0):
    """PEAKS: the -10lgP threshold is set where the FDR over distinct sequences (I = L, best PSM
    each) reaches cfg.peptide_fdr; PSMs of accepted sequences are listed down to -10lgP 20."""
    t = top.copy()
    t["seqkey"] = t["stripped"].str.replace("I", "L")
    best = t.sort_values("minus10lgp", ascending=False).drop_duplicates("seqkey")
    thr = threshold_at_fdr(best["minus10lgp"].to_numpy(), best["is_decoy"].to_numpy(), cfg.peptide_fdr)
    ok = set(best.loc[best["minus10lgp"] >= thr, "seqkey"])
    floor = min(thr, psm_floor)
    acc = t[t["seqkey"].isin(ok) & (t["minus10lgp"] >= floor)]
    # the same peptide can win several precursor entries of one scan; PEAKS lists it once
    acc = acc.sort_values("minus10lgp", ascending=False).drop_duplicates(["source_file", "scan", "peptide_sage", "charge"])
    n_t, n_d = int((~acc["is_decoy"]).sum()), int(acc["is_decoy"].sum())
    stats = dict(peptide_threshold=float(thr), psm_floor=float(floor), n_target_psms=n_t, n_decoy_psms=n_d,
                 psm_fdr=n_d / max(n_t, 1))
    return acc[~acc["is_decoy"]].copy(), stats


def run(cfg: SearchConfig, inputs: list[Path], out: Path, mode: str = "tda", model_path: str | None = None,
        sample_name: str = "Sample 1", processes: int = 8, ref: str | None = None,
        features_path: str | None = None, candidates_per_spectrum: int = 3) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    stats: dict = {}
    fixed = {r: m.delta for m in cfg.fixed_modifications for r in m.residues if len(r) == 1}
    scored_pq = out / "candidates_features.parquet"
    mzmls: list[Path] = []
    for p in inputs:
        if p.suffix.lower() == ".raw":
            from .convert import raw_to_mzml
            p = raw_to_mzml(p, out / "mzml")
        mzmls.append(p)
    if scored_pq.exists():
        feat = pd.read_parquet(scored_pq)
    else:
        parts = []
        for mzml in mzmls:
            rdir = out / mzml.stem
            rdir.mkdir(parents=True, exist_ok=True)
            entries, exp_mzml = stage_precursors(mzml, rdir, features_path, cfg)
            log.info("%s: %d precursor entries (%.0fs)", mzml.name, len(entries), time.time() - t0)
            txt = run_comet(cfg, [exp_mzml], rdir / "comet", num_output_lines=candidates_per_spectrum)[0]
            cands = read_comet(txt, source_file=mzml.name, fixed_mods=fixed)
            cands["proteins"] = cands["proteins"].str.replace("DECOY_", "#DECOY#", regex=False)
            cands = attach_entries(cands, entries)
            cands["rt"] = cands["raw_rt"].fillna(cands["rt"])
            log.info("%s: %d candidates (%.0fs)", mzml.name, len(cands), time.time() - t0)
            model = "HCD2021"
            parts.append(compute_features(cands, mzml, ms2pip_model=model,
                                          ms2_tolerance=cfg.fragment_mass_error_tolerance_da, processes=processes))
        feat = pd.concat(parts, ignore_index=True)
        feat.to_parquet(scored_pq)
    if mode == "imitate":
        with open(model_path, "rb") as fh:
            scored = apply_imitation(pickle.load(fh), feat)
    else:
        scored = rescore_tda(feat)
    top = best_per_spectrum(scored, "score")
    top["minus10lgp"] = minus10lgp(top["score"], top["is_decoy"])
    top.to_parquet(out / "psms_scored.parquet")
    accepted, fstats = peaks_filter(top, cfg)
    stats.update(fstats)
    accepted["rt"] = accepted["entry_rt"].fillna(accepted["rt"]) if "entry_rt" in accepted else accepted["rt"]
    psms_tbl = build_psms(accepted, cfg, sample_name)
    psms_tbl = apply_protein_mapping(psms_tbl, cfg.database)
    peptides_tbl = build_peptides(psms_tbl, sample_name)
    from .proteins import infer_proteins
    proteins_tbl, protpep_tbl = infer_proteins(peptides_tbl, cfg.database, score_threshold=cfg.protein_score_min,
                                               min_unique=cfg.protein_unique_peptides_min, psms_df=psms_tbl,
                                               sample=sample_name)
    ptm_profile = modified = None
    try:
        from .ptm_tables import build_modified_peptides, build_ptm_profile, read_isolation_targets
        ptm_profile = build_ptm_profile(psms_tbl, peptides_tbl, proteins_tbl, protpep_tbl, cfg.database, sample_name)
        iso = {(Path(m).name.replace(".mzML", ".raw"), s): v for m in mzmls for s, v in read_isolation_targets(m).items()}
        modified = build_modified_peptides(psms_tbl, peptides_tbl, proteins_tbl, protpep_tbl, sample_name,
                                           isolation_mz=iso)
    except Exception as e:  # PTM tables are secondary; never fail the run for them
        log.warning("PTM tables skipped: %s", e)
    exp = assemble_export(psms_tbl, peptides_tbl, proteins_tbl, protpep_tbl, ptm_profile=ptm_profile,
                          modified_peptides=modified, sample_name=sample_name, params=cfg.to_dict())
    export_dir = out / "export"
    write_export(exp, str(export_dir))
    stats.update(n_psm_rows=len(psms_tbl), n_peptides=len(peptides_tbl), n_protein_rows=len(proteins_tbl),
                 n_protein_groups=int(proteins_tbl["Protein Group"].nunique()), export_dir=str(export_dir),
                 elapsed_s=time.time() - t0)
    (out / "run_stats.json").write_text(json.dumps(stats, indent=2))
    if ref:
        from .concordance import main as conc_main
        conc_main(["--ref", ref, "--test", str(export_dir), "--out", str(out / "concordance")])
    return stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", nargs="+", required=True, help="mzML or Thermo raw files (one sample)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["tda", "imitate"], default="tda")
    ap.add_argument("--model", help="pickled imitation model (mode imitate)")
    ap.add_argument("--sample-name", default="Sample 1")
    ap.add_argument("--processes", type=int, default=8)
    ap.add_argument("--features", help="precomputed MS1 feature table (parquet), or 'none' to skip features")
    ap.add_argument("--candidates", type=int, default=3, help="candidates kept per spectrum")
    ap.add_argument("--ref", help="PEAKS export dir to compare against")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    cfg = SearchConfig.from_yaml(a.config)
    stats = run(cfg, [Path(m) for m in a.input], Path(a.out), a.mode, a.model, a.sample_name, a.processes,
                a.ref, a.features, a.candidates)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
