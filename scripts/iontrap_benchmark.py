"""Re-score cached pipeline runs with a different fragment-intensity predictor and measure the effect
on ground-truth identification (ProteomeTools pools) - docs/IONTRAP_MODEL.md.

    python scripts/iontrap_benchmark.py --pools pools.txt --source-root /drive_a/drew/peaks-parity-work/proteometools \
        --variant ms2pip=cached --variant ft=peptdeep:hcd=models/hcd.pth,cid=models/cid.pth \
        --out-root /drive_a/drew/peaks-parity-work/iontrap_model/bench --summary bench_summary.csv

For every run <source-root>/<raw stem>/<type> (type in hcd_it, cid_it) with a cached
candidates_features_2pass.parquet: no re-search. The candidate table is copied, its fragment columns
(ms2_corr, ms2_cos, ms2_frac_pred_found, longest_b_obs, longest_y_obs, n_ions_obs, residue_coverage)
are replaced by the variant predictor's values (variant "cached" keeps the MS2PIP HCD2021 values the
pipeline computed), and proteoscan.pipeline.run() is called on the cached tables, which rescored with
rescore_tda, filters with peaks_filter and writes the export. The export is then scored exactly as
scripts/aggregate_proteometools.py scores "ours" (true / false / truncation, measured FDP, true peptides,
true peptides at measured 1% FDP), truth = Supplementary Table 1 + QC peptides.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")  # several runs go in parallel on a shared server

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

FRAG_COLS = ["ms2_corr", "ms2_cos", "ms2_frac_pred_found", "longest_b_obs", "longest_y_obs", "n_ions_obs", "residue_coverage"]
CONFIG = Path(__file__).resolve().parents[1] / "configs" / "peaks_facility_human.yml"
D = Path("/drive_a/drew/public-data/PXD004732")


def prepare(src: Path, out: Path, mzml: Path, spec: str, features_from: Path | None = None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for name in ("candidates_features.parquet", mzml.stem):
        link = out / name
        if not link.exists():
            link.symlink_to(src / name)
    dst = out / "candidates_features_2pass.parquet"
    if dst.exists():
        return
    if features_from is not None and (features_from / "candidates_features_2pass.parquet").exists():
        dst.symlink_to((features_from / "candidates_features_2pass.parquet").resolve())
        return
    if spec == "cached":
        dst.symlink_to(src / "candidates_features_2pass.parquet")
        return
    from proteoscan.rescore_features import peptdeep_correlations
    d = pd.read_parquet(src / "candidates_features_2pass.parquet")
    f = peptdeep_correlations(d, mzml, spec, ms2_tolerance=0.5)
    for c in FRAG_COLS:
        d[c] = f[c].fillna(0.0).to_numpy()
    tmp = dst.with_suffix(".tmp.parquet")
    d.to_parquet(tmp); tmp.rename(dst)


def run_one(job: dict) -> dict:
    import logging
    logging.basicConfig(level=logging.WARNING)
    try:
        import torch
        torch.set_num_threads(2)
    except Exception:
        pass
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", job.get("gpu", "0"))
    from proteoscan.config import SearchConfig
    import proteoscan.pipeline as pl
    from proteoscan.rescore import rescore_tda
    run = pl.run
    seed = int(job.get("seed", 0))
    # mokapot's semi-supervised training depends on its random state: seed > 0 gives replicate rescorings
    pl.rescore_tda = (lambda df, feature_cols=None, seed=seed: rescore_tda(df, feature_cols, seed=seed))
    src, out, mzml = Path(job["src"]), Path(job["out"]), Path(job["mzml"])
    if not (out / "run_stats.json").exists():
        prepare(src, out, mzml, job["spec"], Path(job["features_from"]) if job.get("features_from") else None)
        cfg = SearchConfig.from_yaml(job.get("config", CONFIG))
        run(cfg, [mzml], out, sample_name=f"{job['run']}.{job['type']}", processes=4, two_pass=True,
            ref=job.get("ref"))
    return evaluate(job)


def evaluate(job: dict) -> dict:
    from scripts.aggregate_proteometools import at_measured_fdp
    from scripts.benchmark_proteometools import classify, truth_for_pool
    from proteoscan.seqnorm import strip_mods
    out = Path(job["out"])
    stats = json.loads((out / "run_stats.json").read_text())
    row = dict(pool=job["pool"], run=job["run"], scan_type=job["type"], variant=job["variant"], seed=int(job.get("seed", 0)),
               peptide_threshold=stats.get("peptide_threshold"))
    if job["pool"] is None:
        return row
    truth, src = truth_for_pool(D / "truth", D / "mq", job["pool"])
    idx: dict[str, set[str]] = {}
    for t in truth:
        for i in range(len(t) - 1):
            idx.setdefault(t[i:i + 2], set()).add(t)
    acc = pd.read_csv(out / "export" / "db.psms.csv")
    acc["seq"] = acc["Peptide"].map(strip_mods).str.replace("I", "L")
    acls = acc["seq"].map(lambda q: classify(q, truth, idx)).to_numpy()
    k1, tp1 = at_measured_fdp(acc["-10LgP"].to_numpy(), acls, acc["seq"].to_numpy(), 0.01)
    row.update(psms=len(acc), false=int((acls == "false").sum()),
               fdp=round(float((acls == "false").mean()) if len(acc) else 0, 4),
               true_peptides=int(acc.loc[acls == "true", "seq"].nunique()), truth=len(truth), truth_source=src,
               psms_at_1pct_measured=k1, true_peptides_at_1pct_measured=tp1)
    return row


def jobs_for(pools: list[str], source_root: Path, out_root: Path, variants: dict[str, str], types, gpus, seeds=(0,)):
    jobs, k = [], 0
    for pool in pools:
        runs = sorted(p for p in source_root.glob(f"*-{pool}_01_01-2xIT_2xHCD-*") if p.is_dir())
        for run in runs:
            for t in types:
                src = run / t
                if not (src / "candidates_features_2pass.parquet").exists():
                    continue
                mz = [p for p in (run / "split").glob(f"*.{t}.mzML")] if (run / "split").is_dir() else []
                if not mz:  # validation runs keep the split mzMLs next to the raw conversion
                    mz = list(Path("/drive_a/drew/peaks-parity-work/iontrap_model/runs", run.name, "split").glob(f"*.{t}.mzML"))
                for v, spec in variants.items():
                    for sd in seeds:
                        vd = v if sd == 0 else f"{v}_seed{sd}"
                        jobs.append(dict(pool=pool, run=run.name, type=t, variant=v, spec=spec, src=str(src), mzml=str(mz[0]),
                                         out=str(out_root / vd / run.name / t), gpu=gpus[k % len(gpus)], seed=sd,
                                         features_from=None if sd == 0 else str(out_root / v / run.name / t)))
                        k += 1
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pools", required=True, help="file with pool names")
    ap.add_argument("--source-root", required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--variant", action="append", required=True, help="name=spec; spec 'cached' or 'peptdeep:...'")
    ap.add_argument("--types", nargs="+", default=["hcd_it", "cid_it"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--gpus", default="2,3")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0], help="mokapot seeds (replicate rescorings)")
    a = ap.parse_args()
    variants = dict(v.split("=", 1) for v in a.variant)
    pools = Path(a.pools).read_text().split()
    jobs = jobs_for(pools, Path(a.source_root), Path(a.out_root), variants, a.types, a.gpus.split(","), a.seeds)
    # seed-0 jobs first: the replicates reuse their feature tables
    jobs.sort(key=lambda j: j["seed"])
    print(f"{len(jobs)} jobs", flush=True)
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {}
        for j in jobs:
            env_gpu = j["gpu"]
            futs[ex.submit(_with_gpu, j, env_gpu)] = j
        for f in as_completed(futs):
            j = futs[f]
            try:
                r = f.result(); rows.append(r); print(r, flush=True)
            except Exception as e:
                print("FAILED", j["run"], j["type"], j["variant"], repr(e), flush=True)
    df = pd.DataFrame(rows).sort_values(["scan_type", "pool", "variant"])
    df.to_csv(a.summary, index=False)
    summarize(df)


def _with_gpu(job, gpu):
    os.environ["CUDA_VISIBLE_DEVICES"] = gpu
    return run_one(job)


def summarize(df: pd.DataFrame, exclude=("TUM_first_pool_123",)):
    df = df[~df["pool"].isin(exclude)]
    for t, g in df.groupby("scan_type"):
        print(f"\n## {t}")
        print(g.groupby("variant")[["psms", "fdp", "true_peptides", "psms_at_1pct_measured", "true_peptides_at_1pct_measured"]].mean().round(4).to_string())
        piv = g.pivot_table(index="pool", columns="variant", values="true_peptides_at_1pct_measured")
        vs = list(piv.columns)
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                d = (piv[vs[j]] - piv[vs[i]]).dropna()
                print(f"{vs[j]} - {vs[i]}: mean {d.mean():+.1f}, {int((d > 0).sum())} up / {int((d < 0).sum())} down over {len(d)} pools")


if __name__ == "__main__":
    main()
