"""What are the 'false' PSMs of the ProteomeTools benchmark? For each method: share whose sequence belongs to
another ProteomeTools pool (carry-over between consecutively acquired pools = real peptide, not in this pool's truth),
and for ours the split into the best PSM of a scan and further PSMs of the same scan.


    python scripts/proteometools_false_breakdown.py [--work <benchmark dir>] [--old <earlier benchmark dir>]
"""
import argparse, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
from scripts.benchmark_proteometools import truth_for_pool, classify, peaks_filter
from scripts.aggregate_proteometools import scan_map
from proteoscan.seqnorm import strip_mods
D=Path("/drive_a/drew/public-data/PXD004732"); B=Path("/drive_a/drew/proteoscan-work/benchmarks")
ap = argparse.ArgumentParser()
ap.add_argument("--work", default=str(B/"proteometools_2026-11-12")); ap.add_argument("--old", default=str(B/"proteometools"))
ap.add_argument("--exclude", nargs="*", default=["TUM_first_pool_123"])
a = ap.parse_args()
NEW=Path(a.work); OLD=Path(a.old); PEAKS=Path("/drive_a/drew/vendor-exports/PXD004732_peaks")
allp=pd.read_parquet(D/"truth"/"pool_peptides.parquet"); allp["seq"]=allp["sequence"].str.replace("I","L")
print("pools in truth table:", allp["pool"].nunique(), "peptides:", allp["seq"].nunique())
anypool=set(allp["seq"])
rows=[]
for run in sorted(NEW.glob("*-TUM_first_pool_*")):
    pool=re.search(r"(TUM_first_pool_\d+)", run.name).group(1)
    truth,_=truth_for_pool(D/"truth", D/"mq", pool); idx={}
    for t in truth:
        for i in range(len(t)-1): idx.setdefault(t[i:i+2],set()).add(t)
    smap=scan_map(run/"split")
    def stats(p, method, t, extra=None):
        p=p.copy(); p["seq"]=p["Peptide"].map(strip_mods).str.replace("I","L")
        p["cls"]=p["seq"].map(lambda q: classify(q,truth,idx))
        f=p[p["cls"]=="false"]; other=f["seq"].isin(anypool)
        r=dict(pool=pool,method=method,scan_type=t,psms=len(p),false=len(f),false_other_pool=int(other.sum()),
               false_seqs=f["seq"].nunique(),false_seqs_other_pool=f.loc[other,"seq"].nunique(),
               seqs=p["seq"].nunique(), true_seqs=p.loc[p["cls"]=="true","seq"].nunique(),
               deam=float(p["Peptide"].str.contains(".98)",regex=False).mean()))
        if extra: r.update(extra(p))
        rows.append(r)
    for t in ("hcd_it","cid_it","hcd_ft"):
        f=run/t/"export"/"db.psms.csv"
        if f.exists():
            p=pd.read_csv(f)
            p=p.sort_values("-10LgP",ascending=False); p["rk"]=p.groupby("Scan").cumcount()
            def ex(q):
                fa=q["cls"]=="false"; sec=q["rk"]>0
                return dict(first_psms=int((~sec).sum()),first_false=int((fa&~sec).sum()),second_psms=int(sec.sum()),second_false=int((fa&sec).sum()))
            stats(p,"ours",t,ex)
        oldrun=OLD/run.name
        f=oldrun/t/"export"/"db.psms.csv"
        if f.exists(): stats(pd.read_csv(f),"ours_oct",t)
    for mode in ("dlboost_on","dlboost_off"):
        f=PEAKS/mode/pool/"db.psms.csv"
        if not f.exists(): continue
        p=pd.read_csv(f); p=peaks_filter(p) if "decoy" in p else p
        p["type"]=p["Scan"].map(smap).fillna("other")
        for t in ("hcd_it","cid_it","hcd_ft"): stats(p[p["type"]==t],f"peaks_{mode}",t)
df=pd.DataFrame(rows); df.to_csv(NEW/"false_breakdown.csv",index=False)
df=df[~df.pool.isin(a.exclude)]
g=df.groupby(["scan_type","method"]).sum(numeric_only=True)
out=pd.DataFrame({"psms":g.psms,"fdp_pooled%":100*g.false/g.psms,"false_in_other_pool%":100*g.false_other_pool/g.false,
  "fdp_excl_other_pool%":100*(g.false-g.false_other_pool)/g.psms,"seq_fdp%":100*g.false_seqs/g.seqs,"seq_fdp_excl%":100*(g.false_seqs-g.false_seqs_other_pool)/g.seqs,
  "deam%":None,"first_fdp%":100*g.first_false/g.first_psms,"second_fdp%":100*g.second_false/g.second_psms,"second_psms_per_pool":g.second_psms/df["pool"].nunique()})
out["deam%"]=100*df.groupby(["scan_type","method"])["deam"].mean()
pd.set_option("display.width",250); print(out.round(2).to_string())
