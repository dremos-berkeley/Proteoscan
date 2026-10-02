"""Aggregate the ProteomeTools ground-truth benchmark over all pools: our pipeline per scan type
versus PEAKS (DL boost on/off) at each method's own 1% FDR and at matched measured error.

    python scripts/aggregate_proteometools.py --out /drive_a/drew/peaks-parity-work/proteometools/SUMMARY.md
"""
import argparse, glob, json, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
from scripts.benchmark_proteometools import truth_for_pool, classify, peaks_filter
from peaksparity.seqnorm import strip_mods

W = Path("/drive_a/drew/peaks-parity-work/proteometools")
D = Path("/drive_a/drew/public-data/PXD004732")
PEAKS = Path("/drive_a/drew/vendor-exports/PXD004732_peaks")
TYPES = ("hcd_it", "cid_it", "hcd_ft")


def scan_map(split_dir: Path) -> dict[int, str]:
    import pyopenms as oms
    m = {}
    for mz in split_dir.glob("*.mzML"):
        t = mz.stem.split(".")[-1]
        od = oms.OnDiscMSExperiment(); od.openFile(str(mz)); meta = od.getMetaData()
        for i in range(meta.getNrSpectra()):
            if meta[i].getMSLevel() == 2:
                mm = re.search(r"scan=(\d+)", meta[i].getNativeID())
                if mm:
                    m[int(mm.group(1))] = t
    return m


def curve(df: pd.DataFrame, truth, idx):
    """df sorted by descending score with columns seq; returns cumulative false/true arrays."""
    cls = df["seq"].map(lambda q: classify(q, truth, idx)).to_numpy()
    false = np.cumsum(cls == "false"); n = np.arange(1, len(df) + 1)
    return cls, false / n


def true_peptides_at(df, cls, fdp, target):
    k = int(np.sum(fdp <= target))
    sub = df.iloc[:k]
    return k, sub.loc[cls[:k] == "true", "seq"].nunique(), int((cls[:k] == "false").sum())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    rows = []
    for pool in D.joinpath("pools.txt").read_text().split():
        runs = [p for p in W.glob(f"*-{pool}_01_01-2xIT_2xHCD-*") if (p / "split").is_dir()]
        if not runs:
            continue
        run = runs[0]
        truth, src = truth_for_pool(D / "truth", D / "mq", pool)
        idx = {}
        for t in truth:
            for i in range(len(t) - 1):
                idx.setdefault(t[i:i + 2], set()).add(t)
        smap = scan_map(run / "split")
        ours = {}
        for t in TYPES:
            pq = run / t / "psms_scored.parquet"
            if not pq.exists():
                continue
            s = pd.read_parquet(pq, columns=["scan", "peptide_sage", "charge", "stripped", "score", "is_decoy", "minus10lgp"])
            s = s[~s["is_decoy"]].sort_values("score", ascending=False).drop_duplicates(["scan", "peptide_sage", "charge"]).reset_index(drop=True)
            s["seq"] = s["stripped"].str.replace("I", "L")
            cls, fdp = curve(s, truth, idx)
            acc = pd.read_csv(run / t / "export" / "db.psms.csv") if (run / t / "export" / "db.psms.csv").exists() else None
            n_acc = len(acc) if acc is not None else 0
            ours[t] = dict(df=s, cls=cls, fdp=fdp, n_acc=n_acc)
            if acc is not None:  # the pipeline's actual accepted set (PEAKS-style filter)
                acc["seq"] = acc["Peptide"].map(strip_mods).str.replace("I", "L")
                acls = acc["seq"].map(lambda q: classify(q, truth, idx)).to_numpy()
                rows.append(dict(pool=pool, method="ours", scan_type=t, psms=len(acc), false=int((acls == "false").sum()),
                                 fdp=round(float((acls == "false").mean()) if len(acc) else 0, 4),
                                 true_peptides=int(acc.loc[acls == "true", "seq"].nunique()), truth=len(truth), truth_source=src))
        for mode in ("dlboost_on", "dlboost_off"):
            pdir = PEAKS / mode / pool
            if not (pdir / "db.psms.csv").exists():
                continue
            p = pd.read_csv(pdir / "db.psms.csv")
            p = peaks_filter(p) if "decoy" in p else p
            p["seq"] = p["Peptide"].map(strip_mods).str.replace("I", "L")
            p["type"] = p["Scan"].map(smap).fillna("other")
            for t in TYPES:
                g = p[p["type"] == t]
                if not len(g):
                    continue
                cls = g["seq"].map(lambda q: classify(q, truth, idx)).to_numpy()
                fdp = float((cls == "false").mean())
                row = dict(pool=pool, method=f"peaks_{mode}", scan_type=t, psms=len(g), false=int((cls == "false").sum()), fdp=round(fdp, 4),
                           true_peptides=int(g.loc[cls == "true", "seq"].nunique()), truth=len(truth), truth_source=src)
                if t in ours:
                    k, tp, f = true_peptides_at(ours[t]["df"], ours[t]["cls"], ours[t]["fdp"], fdp)
                    row.update(ours_psms_at_matched_fdp=k, ours_true_peptides_at_matched_fdp=tp)
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(Path(a.out).with_suffix(".csv"), index=False)
    lines = ["# ProteomeTools ground-truth benchmark", "", f"Pools with results: {df['pool'].nunique()}", ""]
    for t in TYPES:
        sub = df[df["scan_type"] == t]
        if not len(sub):
            continue
        lines += [f"## {t}", "", "| method | pools | PSMs (mean) | measured FDP (mean) | true peptides (mean) | pool recall (mean) |", "|---|---|---|---|---|---|"]
        for m, g in sub.groupby("method"):
            lines.append(f"| {m} | {len(g)} | {g['psms'].mean():.0f} | {g['fdp'].mean():.4f} | {g['true_peptides'].mean():.0f} | {(g['true_peptides']/g['truth']).mean():.3f} |")
        pk = sub[(sub["method"] == "peaks_dlboost_on") & sub["ours_true_peptides_at_matched_fdp"].notna()] if "ours_true_peptides_at_matched_fdp" in sub else pd.DataFrame()
        if len(pk):
            lines += ["", f"At PEAKS's measured error rate (mean {pk['fdp'].mean():.4f}): PEAKS {pk['true_peptides'].mean():.0f} true peptides, ours {pk['ours_true_peptides_at_matched_fdp'].mean():.0f} (over {len(pk)} pools; ours higher in {(pk['ours_true_peptides_at_matched_fdp'] > pk['true_peptides']).sum()})."]
        lines.append("")
    Path(a.out).write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
