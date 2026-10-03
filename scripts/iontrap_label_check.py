"""Agreement between the targeted labels (iontrap_label.py) and the pipeline's PSMs on the same runs.

    python scripts/iontrap_label_check.py --runs-root .../iontrap_model/runs --pipeline-root .../iontrap_model/val_pipeline \
        --out label_check.csv

Pipeline PSMs used: target PSMs of the export whose stripped peptide (I=L) is a pool member and -10LgP >= 25.
Reported per run and scan type:
  scans_both            scans carrying a label and a confident pipeline PSM
  agree                 of those, same peptide (I=L) and charge
  agree_modseq          same exact modified peptide and charge
  label_recall          share of the confident pipeline (peptide, charge) pairs that have a label
  labels_with_pipeline_psm_any   share of labelled scans with any accepted pipeline PSM
  labels_contradicted   labelled scans where the pipeline accepted a different peptide (any score)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan.seqnorm import canonical, strip_mods  # noqa: E402
from scripts.iontrap_label import pool_peptides  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-root", required=True)
    ap.add_argument("--pipeline-root", required=True)
    ap.add_argument("--min-score", type=float, default=25.0)
    ap.add_argument("--labels-name", default="labels_ev.parquet")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for run in sorted(Path(a.pipeline_root).iterdir()):
        lab_f = Path(a.runs_root) / run.name / a.labels_name
        if not lab_f.exists():
            continue
        lab = pd.read_parquet(lab_f)
        pool = lab["pool"].iloc[0]
        members = {p.replace("I", "L") for p in pool_peptides(pool)}
        for t in ("hcd_it", "cid_it"):
            ex = run / t / "export" / "db.psms.csv"
            if not ex.exists():
                continue
            p = pd.read_csv(ex)
            p["seq"] = p["Peptide"].map(strip_mods).str.replace("I", "L")
            p["canon"] = p["Peptide"].map(canonical)
            conf = p[p["seq"].isin(members) & (p["-10LgP"] >= a.min_score)]
            l = lab[lab["scan_type"] == t].copy()
            l["seq"] = l["sequence"].str.replace("I", "L")
            l["canon"] = l["modseq"].map(canonical)
            m = l.merge(conf[["Scan", "seq", "canon", "z", "-10LgP"]], left_on="scan", right_on="Scan", how="inner",
                        suffixes=("", "_p"))
            # a scan can carry several pipeline PSMs (precursor entries): agreement if any matches
            m["agree"] = (m["seq"] == m["seq_p"]) & (m["charge"] == m["z"])
            m["agree_mod"] = (m["canon"] == m["canon_p"]) & (m["charge"] == m["z"])
            per_scan = m.groupby("scan")[["agree", "agree_mod"]].any()
            pairs = set(zip(conf["seq"], conf["z"]))
            lpairs = set(zip(l["seq"], l["charge"]))
            anyp = p.groupby("Scan")["seq"].apply(set)
            contra = sum(1 for s, q in zip(l["scan"], l["seq"]) if s in anyp.index and q not in anyp[s])
            rows.append(dict(run=run.name, pool=pool, scan_type=t, labels=len(l), pipeline_conf_psms=len(conf),
                             scans_both=len(per_scan), agree=round(per_scan["agree"].mean(), 4) if len(per_scan) else None,
                             agree_modseq=round(per_scan["agree_mod"].mean(), 4) if len(per_scan) else None,
                             label_recall=round(len(pairs & lpairs) / max(len(pairs), 1), 4),
                             labels_with_pipeline_psm_any=round(l["scan"].isin(anyp.index).mean(), 4),
                             labels_contradicted=contra, labels_contradicted_frac=round(contra / max(len(l), 1), 4)))
            print(rows[-1], flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    if len(df):
        print(df.groupby("scan_type")[["labels", "scans_both", "agree", "agree_modseq", "label_recall",
                                       "labels_with_pipeline_psm_any", "labels_contradicted_frac"]].mean().round(4).to_string())


if __name__ == "__main__":
    main()
