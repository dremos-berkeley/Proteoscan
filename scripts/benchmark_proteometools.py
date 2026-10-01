"""Ground-truth accuracy benchmark on ProteomeTools synthetic peptide pools.

    python scripts/benchmark_proteometools.py --pool TUM_first_pool_30 --truth-dir <dir with pool truth>
        --export <our export dir or PEAKS export dir> [--export ...] --out results.json

Truth per pool: the pool's peptide list (Supplementary Table 1 of Zolg et al. 2017 if present as
truth/pool_peptides.parquet with columns pool, sequence), else the union of MaxQuant identifications
across the pool's ProteomeTools result archives (mq/<pool>_*-tryptic.zip, peptides.txt), which were
searched against the pool-specific database and therefore enumerate most of the pool.

An accepted PSM is classified as:
  true        stripped sequence (I=L) is a pool member
  truncation  a sub-sequence of a pool member (synthesis truncation products are real molecules)
  false       anything else (including peptides from other pools: all are human sequences)
Reported per export: PSMs and distinct peptides accepted, counts by class, false-discovery
proportion FDP = false / accepted, and the same restricted to each MS2 scan type when the export's
Source File names carry the split suffix (.hcd_it, .cid_it, .hcd_ft).
"""
import argparse, io, json, re, sys, zipfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from peaksparity.seqnorm import strip_mods


def truth_from_mq(mq_dir: Path, pool: str) -> set[str]:
    seqs = set()
    for z in sorted(mq_dir.glob(f"{pool}_*.zip")):
        with zipfile.ZipFile(z) as zf:
            names = [n for n in zf.namelist() if n.endswith("peptides.txt")]
            if not names:
                continue
            p = pd.read_csv(io.BytesIO(zf.read(names[0])), sep="\t", usecols=["Sequence", "Reverse", "Potential contaminant"], low_memory=False)
            p = p[(p["Reverse"] != "+") & (p["Potential contaminant"] != "+")]
            seqs |= set(p["Sequence"].str.replace("I", "L"))
    return seqs


def truth_for_pool(truth_dir: Path, mq_dir: Path, pool: str) -> tuple[set[str], str]:
    pq = truth_dir / "pool_peptides.parquet"
    if pq.exists():
        t = pd.read_parquet(pq)
        s = set(t.loc[t["pool"] == pool, "sequence"].str.replace("I", "L"))
        if s:
            return s, "supplementary_table"
    return truth_from_mq(mq_dir, pool), "maxquant_union"


def classify(seq: str, truth: set[str], truth_sub: dict[str, set[str]]) -> str:
    if seq in truth:
        return "true"
    if len(seq) >= 6 and any(seq in t for t in truth_sub.get(seq[:2], ())):
        return "truncation"
    return "false"


def evaluate(export_dir: Path, truth: set[str]) -> dict:
    psms = pd.read_csv(export_dir / "db.psms.csv")
    if "decoy" in psms:
        psms = psms[psms["decoy"].astype(str).str.lower() != "true"]
    psms["seq"] = psms["Peptide"].map(strip_mods).str.replace("I", "L")
    idx: dict[str, set[str]] = {}
    for t in truth:
        for i in range(len(t) - 1):
            idx.setdefault(t[i:i + 2], set()).add(t)
    psms["cls"] = psms["seq"].map(lambda s: classify(s, truth, idx))
    def summarise(d):
        n = len(d); c = d["cls"].value_counts().to_dict()
        peps = d.drop_duplicates("seq"); cp = peps["cls"].value_counts().to_dict()
        return dict(psms=n, psm_true=c.get("true", 0), psm_truncation=c.get("truncation", 0), psm_false=c.get("false", 0),
                    psm_fdp=round(c.get("false", 0) / max(n, 1), 4), peptides=len(peps), pep_true=cp.get("true", 0),
                    pep_false=cp.get("false", 0), pool_recall=round(cp.get("true", 0) / max(len(truth), 1), 4))
    out = {"all": summarise(psms)}
    if "Source File" in psms:
        for sf, g in psms.groupby("Source File"):
            m = re.search(r"\.(hcd_it|cid_it|hcd_ft|etd_it|ethcd_ft)", str(sf))
            out[m.group(1) if m else str(sf)] = summarise(g)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True); ap.add_argument("--truth-dir", default="/drive_a/drew/public-data/PXD004732/truth")
    ap.add_argument("--mq-dir", default="/drive_a/drew/public-data/PXD004732/mq")
    ap.add_argument("--export", nargs="+", required=True, help="export dirs; use name=path to label")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    truth, source = truth_for_pool(Path(a.truth_dir), Path(a.mq_dir), a.pool)
    res = {"pool": a.pool, "truth_source": source, "truth_size": len(truth), "exports": {}}
    for e in a.export:
        name, _, path = e.partition("=") if "=" in e else (Path(e).name, "", e)
        res["exports"][name] = evaluate(Path(path), truth)
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
