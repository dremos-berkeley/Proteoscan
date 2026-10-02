"""Append de novo agreement features (proteoscan.denovo_features) to a candidate table.

    python scripts/add_denovo_features.py --candidates <run>/candidates_features_2pass.parquet \
        --denovo <run>.cand.denovo_all.parquet [<run>.rest.denovo_all.parquet ...] \
        --out <run>/candidates_features_2pass_denovo.parquet

The de novo tables are the output of denovo_features.run_denovo_queries (one row per beam
candidate, key query_id, ``picked`` = precursor-consistent choice). Candidates are joined on
query_id; queries without a de novo call get NaN features.
"""
import argparse, logging, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from proteoscan.denovo_features import DENOVO_FEATURE_COLUMNS, add_denovo_features


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--denovo", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--key", default="query_id")
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--conf-threshold", type=float, default=0.3)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    t = time.time()
    c = pd.read_parquet(a.candidates)
    d = pd.concat([pd.read_parquet(p, columns=[a.key, "sequence", "per_residue_scores", "picked"])
                   for p in a.denovo], ignore_index=True)
    t1 = time.time()
    out = add_denovo_features(c, d, key=a.key, tol=a.tol, conf_threshold=a.conf_threshold,
                              threads=a.threads)
    t2 = time.time()
    out.to_parquet(a.out)
    cov = out["caa"].notna().mean()
    logging.info("%s: %d candidates, %.1f%% with a de novo call; read %.0fs, features %.0fs, write %.0fs",
                 a.out, len(out), 100 * cov, t1 - t, t2 - t1, time.time() - t2)
    logging.info("\n%s", out[DENOVO_FEATURE_COLUMNS].describe().T.to_string())


if __name__ == "__main__":
    main()
