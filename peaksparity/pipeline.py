"""End-to-end driver for one or more runs: raw/mzML -> PEAKS-format export directory.

    python -m peaksparity.pipeline --config configs/peaks_facility_human.yml \
        --mzml run.mzML --out work/run1 [--chunks 16] [--mode tda|imitate --model model.pkl]

Stages (see docs/ARCHITECTURE.md): search (Sage, chunked) -> features -> rescore -> one PSM per
spectrum -> PEAKS-style FDR thresholds -> protein inference -> export -> (optional) concordance.
"""
from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .config import SearchConfig
from .export_builder import assemble_export, build_peptides, build_psms
from .peaks_export import write_export
from .rescore import best_per_spectrum, minus10lgp, rescore_tda, apply_imitation
from .rescore_features import compute_features
from .sage_search import read_sage, run_sage, run_sage_chunked

log = logging.getLogger("peaksparity.pipeline")


def stage_search(cfg: SearchConfig, mzmls: list[Path], out: Path, chunks: int, report_psms: int) -> pd.DataFrame:
    sdir = out / "search"
    tsv = sdir / "results.sage.tsv"
    if not tsv.exists():
        if chunks > 1:
            run_sage_chunked(cfg, mzmls, sdir, n_chunks=chunks, report_psms=report_psms)
        else:
            run_sage(cfg, mzmls, sdir, report_psms=report_psms)
    return read_sage(tsv)


def stage_features(df: pd.DataFrame, mzmls: list[Path], out: Path, cfg: SearchConfig, processes: int) -> pd.DataFrame:
    pq = out / "features.parquet"
    if pq.exists():
        return pd.read_parquet(pq)
    parts = []
    by_file = {Path(m).name: m for m in mzmls}
    model = "HCD2021" if not cfg.low_res_fragments else "CID"  # ion-trap HCD: low-res model; tuned later
    for src, part in df.groupby("source_file"):
        mz = by_file.get(src)
        if mz is None:
            raise FileNotFoundError(f"no mzML for {src}")
        parts.append(compute_features(part.reset_index(drop=True), mz, ms2pip_model=model,
                                      ms2_tolerance=cfg.fragment_mass_error_tolerance_da, processes=processes))
    feat = pd.concat(parts, ignore_index=True)
    feat.to_parquet(pq)
    return feat


def stage_rescore(feat: pd.DataFrame, mode: str, model_path: str | None) -> pd.DataFrame:
    if mode == "imitate":
        import pickle
        with open(model_path, "rb") as fh:
            model = pickle.load(fh)
        return apply_imitation(model, feat)
    return rescore_tda(feat)


def stage_threshold(psms: pd.DataFrame, cfg: SearchConfig) -> tuple[pd.DataFrame, dict]:
    """PEAKS: choose the -10lgP threshold at the requested peptide-sequence FDR (decoys/targets
    over best PSM per peptidoform), then keep all PSMs (targets) at or above it."""
    from .fdr import threshold_at_fdr  # provided by the proteins/fdr stage
    best = psms.sort_values("minus10lgp", ascending=False).drop_duplicates(["peptide_sage"])
    thr = threshold_at_fdr(best["minus10lgp"].to_numpy(), best["is_decoy"].to_numpy(), cfg.peptide_fdr)
    keep = psms[(psms["minus10lgp"] >= thr)]
    stats = {
        "peptide_threshold": float(thr),
        "psm_fdr": float(keep["is_decoy"].sum() / max((~keep["is_decoy"]).sum(), 1)),
        "n_target_psms": int((~keep["is_decoy"]).sum()),
    }
    return keep[~keep["is_decoy"]].copy(), stats


def run(cfg: SearchConfig, mzmls: list[Path], out: Path, chunks: int = 1, report_psms: int = 5,
        mode: str = "tda", model_path: str | None = None, sample_name: str = "Sample 1",
        processes: int = 8, ref: str | None = None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cands = stage_search(cfg, mzmls, out, chunks, report_psms)
    log.info("search: %d candidate rows in %.0fs", len(cands), time.time() - t0)
    feat = stage_features(cands, mzmls, out, cfg, processes)
    log.info("features: %d rows in %.0fs", len(feat), time.time() - t0)
    scored = stage_rescore(feat, mode, model_path)
    top = best_per_spectrum(scored, "score")
    top["minus10lgp"] = minus10lgp(top["score"], top["is_decoy"])
    top.to_parquet(out / "psms_scored.parquet")
    accepted, stats = stage_threshold(top, cfg)
    log.info("threshold: %s", stats)
    psms_tbl = build_psms(accepted, cfg, sample_name)
    peptides_tbl = build_peptides(psms_tbl, sample_name)
    proteins_tbl = protein_peptides_tbl = None
    try:
        from .proteins import infer_proteins
        proteins_tbl, protein_peptides_tbl = infer_proteins(peptides_tbl, cfg.database,
                                                             score_threshold=cfg.protein_score_min,
                                                             min_unique=cfg.protein_unique_peptides_min,
                                                             psms=psms_tbl, sample_name=sample_name)
    except ImportError:
        log.warning("proteins stage not available yet; exporting PSMs and peptides only")
    exp = assemble_export(psms_tbl, peptides_tbl, proteins_tbl, protein_peptides_tbl, sample_name=sample_name,
                          params=cfg.to_dict())
    export_dir = out / "export"
    write_export(exp, export_dir)
    stats["export_dir"] = str(export_dir)
    stats["elapsed_s"] = time.time() - t0
    (out / "run_stats.json").write_text(json.dumps(stats, indent=2))
    if ref:
        from .concordance import main as conc_main
        conc_main(["--ref", ref, "--test", str(export_dir), "--out", str(out / "concordance")])
    return stats


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--mzml", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunks", type=int, default=1)
    ap.add_argument("--report-psms", type=int, default=5)
    ap.add_argument("--mode", choices=["tda", "imitate"], default="tda")
    ap.add_argument("--model")
    ap.add_argument("--sample-name", default="Sample 1")
    ap.add_argument("--processes", type=int, default=8)
    ap.add_argument("--ref", help="PEAKS export dir to compare against")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    cfg = SearchConfig.from_yaml(a.config)
    stats = run(cfg, [Path(m) for m in a.mzml], Path(a.out), a.chunks, a.report_psms, a.mode, a.model,
                a.sample_name, a.processes, a.ref)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
