"""End-to-end driver: mzML (or Thermo raw) -> PEAKS-format export directory.

    python -m proteoscan.pipeline --config configs/peaks_facility_human.yml \
        --input run.mzML --out work/run1 [--sample-name "Sample 1"] [--ref <PEAKS export dir>]
        [--mode tda|imitate --model model.pkl] [--features features.parquet]

Stages (docs/ARCHITECTURE.md):
  1 convert      raw -> mzML (skipped for mzML input)
  2 precursors   MS1 features + PEAKS-style precursor list -> expanded mzML
  3 search       Comet (low memory) on the expanded mzML
  4 features     fragment-intensity / retention-time / ion-series evidence per candidate
  5 rescore      target-decoy (mokapot) or PEAKS-imitation model
  6 filter       one peptide per (scan, precursor); threshold at the peptide-sequence FDR
  7 export       PSMs (I/L variants split), peptides, proteins, features, PTM tables in PEAKS CSV format
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

from dataclasses import replace as dataclasses_replace

from .comet_search import fixed_mod_deltas, read_comet, run_comet
from .config import SearchConfig
from .export_builder import apply_protein_mapping, assemble_export, build_peptides, build_psms
from .fdr import threshold_summary
from .monoiso import PROTON
from .peaks_export import write_export
from .rescore import apply_imitation, best_per_spectrum, minus10lgp, rescore_tda
from .rescore_features import compute_features

log = logging.getLogger("proteoscan.pipeline")


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
        fcfg = Path(__file__).resolve().parent / "data" / "features_lumos_120k.yml"
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
    c["entry_mz"] = e["mz"].reindex(q).to_numpy()
    c["entry_source"] = e["source"].reindex(q).to_numpy()
    c["raw_rt"] = (e["rt_min"].reindex(q).to_numpy() if "rt_min" in e else c["rt"].to_numpy())
    c["raw_mz"] = (e["precursor_mz"].reindex(q).to_numpy() if "precursor_mz" in e else np.nan)
    # rank/delta columns are per query already; spectrum grouping must use scan + precursor_idx
    return c


def stage_two_pass(feat: pd.DataFrame, cfg: SearchConfig, out: Path, mzmls: list[Path], fixed: dict,
                   processes: int, candidates_per_spectrum: int, fragment_model: str = "HCD2021") -> pd.DataFrame:
    """PEAKS DB-style second pass: rescore pass 1, accept, re-search every query against the
    accepted proteins (plus their decoys), merge candidates and add the decoy-fused protein-support
    feature. Validated by paired entrapment (docs/STATUS.md)."""
    from .decoys import read_fasta
    from .twopass import add_protein_features, merge_passes, protein_support, rerank
    import proteoscan.rescore_features as rf
    pq = out / "candidates_features_2pass.parquet"
    if pq.exists():
        feat2 = pd.read_parquet(pq)
        if "log_prot_support" not in rf.FEATURE_COLUMNS:
            rf.FEATURE_COLUMNS.append("log_prot_support")
        return feat2
    scored1 = rescore_tda(feat)
    top1 = best_per_spectrum(scored1, "score")
    top1["minus10lgp"] = minus10lgp(top1["score"], top1["is_decoy"])
    acc1, st1 = peaks_filter(top1, cfg)
    log.info("pass 1: %d target PSMs accepted", len(acc1))
    ids = set()
    for p in acc1["proteins"]:
        for x in str(p).split(";"):
            if not x.startswith("#DECOY#"):
                ids.add(x.split()[0])
    sub_fasta = out / "pass1_proteins.fasta"
    with open(sub_fasta, "w") as fh:
        for h, seq in read_fasta(cfg.database):
            if h.split()[0] in ids:
                fh.write(f">{h}\n{seq}\n")
    cfg2 = dataclasses_replace(cfg, database=str(sub_fasta))
    parts = []
    for mzml in mzmls:
        rdir = out / mzml.stem
        entries = pd.read_parquet(rdir / "precursors.parquet")
        exp_mzml = rdir / (mzml.stem + ".expanded.mzML")
        txt = run_comet(cfg2, [exp_mzml], rdir / "comet_pass2", num_output_lines=candidates_per_spectrum + 2)[0]
        c2 = read_comet(txt, source_file=mzml.name, fixed_mods=fixed)
        c2["proteins"] = c2["proteins"].str.replace("DECOY_", "#DECOY#", regex=False)
        c2 = attach_entries(c2, entries)
        c2["rt"] = c2["raw_rt"].fillna(c2["rt"])
        c1 = feat[feat["source_file"] == mzml.name]
        merged = merge_passes(c1, c2, keep=candidates_per_spectrum + 2)
        new = merged.loc[merged["from_pass2"] == 1, [c for c in c2.columns if c in merged.columns]].reset_index(drop=True)
        if len(new):
            ref = acc1[acc1["source_file"] == mzml.name]
            new = compute_features(new, mzml, ms2pip_model=fragment_model,
                                   ms2_tolerance=cfg.fragment_mass_error_tolerance_da, processes=processes,
                                   rt_reference=ref)
            new["from_pass2"] = 1
        old = merged[merged["from_pass2"] == 0]
        both = rerank(pd.concat([old, new], ignore_index=True, sort=False), keep=candidates_per_spectrum + 2)
        parts.append(both)
    feat2 = pd.concat(parts, ignore_index=True, sort=False)
    support, seqs = protein_support(acc1)
    feat2 = add_protein_features(feat2, support, seqs)
    if "log_prot_support" not in rf.FEATURE_COLUMNS:
        rf.FEATURE_COLUMNS.append("log_prot_support")
    feat2.to_parquet(pq)
    log.info("pass 2: %d candidates (%d proteins in the restricted database)", len(feat2), len(ids))
    return feat2


def stage_features_table(psms_tbl: pd.DataFrame, mzmls: list[Path], out: Path, features_path: str | None,
                         sample_name: str, denovo_df: pd.DataFrame | None = None) -> pd.DataFrame | None:
    """<Sample>.db.features.csv: the MS1 features of each run that carry a PSM (or de novo result)."""
    from .features_export import build_features_table, features_from_entries
    if features_path == "none":
        return None
    parts = []
    for mzml in mzmls:
        rdir = out / mzml.stem
        src = mzml.name.replace(".mzML", ".raw") if mzml.name.endswith(".mzML") else mzml.name
        fpq = Path(features_path) if features_path else rdir / "ms1_features.parquet"
        if fpq.exists():
            f = pd.read_parquet(fpq)
        elif (rdir / "precursors.parquet").exists():
            log.warning("%s: no feature table at %s; features.csv RT Begin/End fall back to the apex", mzml.name, fpq)
            f = features_from_entries(pd.read_parquet(rdir / "precursors.parquet"))
        else:
            continue
        f = f.copy()
        f["source_file"] = src
        parts.append(f)
    if not parts:
        return None
    return build_features_table(pd.concat(parts, ignore_index=True, sort=False), psms_tbl, denovo_df, sample_name)


MIN_DECOYS_AT_THRESHOLD = 10  # fewer decoy sequences at the threshold: the FDR estimate is not informative


def peaks_filter(top: pd.DataFrame, cfg: SearchConfig, psm_floor: float = 20.0):
    """PEAKS: the -10lgP threshold is set where the FDR over distinct sequences (I = L, best PSM
    each) reaches cfg.peptide_fdr; PSMs of accepted sequences are listed down to -10lgP 20.

    Besides the accepted PSMs, returns stats that say whether the FDR is estimable: the sequence
    counts at the threshold, D/T and the conservative (D+1)/T, the size of the decoy null behind
    minus10lgp and its old ceiling 10*log10(n_decoys+1), whether the threshold lies above every decoy,
    and fdr_estimable (>= MIN_DECOYS_AT_THRESHOLD decoy sequences at the threshold)."""
    t = top.copy()
    t["seqkey"] = t["stripped"].str.replace("I", "L")
    best = t.sort_values("minus10lgp", ascending=False).drop_duplicates("seqkey")
    bs, bd = best["minus10lgp"].to_numpy(), best["is_decoy"].to_numpy(bool)
    summ = threshold_summary(bs, bd, cfg.peptide_fdr)
    thr = summ["threshold"]
    ok = set(best.loc[best["minus10lgp"] >= thr, "seqkey"])
    floor = min(thr, psm_floor)
    acc = t[t["seqkey"].isin(ok) & (t["minus10lgp"] >= floor)]
    # the same peptide can win several precursor entries of one scan; PEAKS lists it once
    acc = acc.sort_values("minus10lgp", ascending=False).drop_duplicates(["source_file", "scan", "peptide_sage", "charge"])
    n_t, n_d = int((~acc["is_decoy"]).sum()), int(acc["is_decoy"].sum())
    n_null = int(t["is_decoy"].sum())  # decoy spectrum winners: the null minus10lgp was computed on
    finite = bool(np.isfinite(thr))
    T, D = (summ["n_target"], summ["n_decoy"]) if finite else (0, 0)
    dec_scores = bs[bd & ~np.isnan(bs)]
    above_all = bool(finite and (len(dec_scores) == 0 or thr > dec_scores.max()))
    stats = dict(peptide_threshold=float(thr) if finite else None, psm_floor=float(floor),
                 n_target_psms=n_t, n_decoy_psms=n_d, psm_fdr=n_d / max(n_t, 1),
                 n_target_seqs_at_threshold=int(T), n_decoy_seqs_at_threshold=int(D),
                 sequence_fdr=(D / T) if T else None, sequence_fdr_conservative=((D + 1) / T) if T else None,
                 decoy_null_size=n_null, null_ceiling=float(10.0 * np.log10(n_null + 1)) if n_null else None,
                 threshold_above_all_decoys=above_all, fdr_estimable=bool(D >= MIN_DECOYS_AT_THRESHOLD))
    if not finite:
        log.warning("peaks_filter: no threshold reaches %.3g sequence FDR; nothing accepted", cfg.peptide_fdr)
    elif above_all:
        log.warning("peaks_filter: threshold %.2f lies above every decoy sequence (%d target sequences accepted); "
                    "the reported FDR rests on 0 decoys and is not estimable (conservative (D+1)/T = %.3g)",
                    thr, T, 1.0 / T if T else float("nan"))
    elif D < MIN_DECOYS_AT_THRESHOLD:
        log.warning("peaks_filter: only %d decoy sequence(s) at threshold %.2f (%d targets); FDR estimate is "
                    "imprecise (D/T = %.3g, (D+1)/T = %.3g)", D, thr, T, D / T, (D + 1) / T)
    return acc[~acc["is_decoy"]].copy(), stats


def run(cfg: SearchConfig, inputs: list[Path], out: Path, mode: str = "tda", model_path: str | None = None,
        sample_name: str = "Sample 1", processes: int = 8, ref: str | None = None,
        features_path: str | None = None, candidates_per_spectrum: int = 3, two_pass: bool = True,
        fragment_model: str = "HCD2021") -> dict:
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    stats: dict = {}
    fixed = fixed_mod_deltas(cfg)
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
            parts.append(compute_features(cands, mzml, ms2pip_model=fragment_model,
                                          ms2_tolerance=cfg.fragment_mass_error_tolerance_da, processes=processes))
        feat = pd.concat(parts, ignore_index=True)
        feat.to_parquet(scored_pq)
    if two_pass:
        feat = stage_two_pass(feat, cfg, out, mzmls, fixed, processes, candidates_per_spectrum, fragment_model)
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
    features_tbl = None
    try:
        features_tbl = stage_features_table(psms_tbl, mzmls, out, features_path, sample_name)
    except Exception as e:  # secondary table; never fail the run for it
        log.warning("features table skipped: %s", e)
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
    exp = assemble_export(psms_tbl, peptides_tbl, proteins_tbl, protpep_tbl, features=features_tbl,
                          ptm_profile=ptm_profile, modified_peptides=modified, sample_name=sample_name,
                          params=cfg.to_dict())
    export_dir = out / "export"
    write_export(exp, str(export_dir))
    fid = psms_tbl["Feature Id"].fillna(0).astype(int)
    stats.update(n_features_exported=0 if features_tbl is None else len(features_tbl),
                 psm_feature_fraction=float((fid > 0).mean()) if len(fid) else None)
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
    ap.add_argument("--single-pass", action="store_true", help="skip the protein-restricted second pass")
    ap.add_argument("--fragment-model", default="HCD2021",
                    help="fragment-intensity predictor: an MS2PIP model name, or peptdeep:<ckpt> / "
                         "peptdeep:hcd=<ckpt>,cid=<ckpt> (docs/IONTRAP_MODEL.md)")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    cfg = SearchConfig.from_yaml(a.config)
    stats = run(cfg, [Path(m) for m in a.input], Path(a.out), a.mode, a.model, a.sample_name, a.processes,
                a.ref, a.features, a.candidates, two_pass=not a.single_pass, fragment_model=a.fragment_model)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
