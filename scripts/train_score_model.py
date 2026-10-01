"""Train a model of PEAKS's -10lgP from our candidate features, using a decoy-inclusive PEAKS export.

    python scripts/train_score_model.py --peaks-psms <db.psms.csv with decoy column> \
        --candidates <run>/candidates_features_2pass.parquet [...] --out models/score_v1.pkl [--cv]

Labels: for every candidate that PEAKS also lists for the same scan, modified sequence and charge,
the PEAKS -10lgP; for candidates PEAKS does not list (score below its export floor, or not its
choice for the scan) a censored value below the floor. Cross-validation is by scan parity, and
the evaluation reproduces PEAKS's accepted set (its own 1% peptide-sequence FDR threshold).
"""
import argparse, json, pickle, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
import lightgbm as lgb
from peaksparity.fdr import threshold_at_fdr
from peaksparity.rescore import best_per_spectrum
from peaksparity.seqnorm import canonical, strip_mods
import peaksparity.rescore_features as rf

FLOOR = 15.0


def key(scan, modseq, z):
    return (scan.astype(str) + "|" + modseq + "|" + z.astype(str)).str.replace("I", "L")


def load_peaks(path):
    p = pd.read_csv(path)
    p["decoy"] = p["decoy"].astype(str).str.lower().eq("true") if "decoy" in p else False
    p["key"] = key(p["Scan"], p["Peptide"], p["z"])
    p["seq"] = p["Peptide"].map(strip_mods).str.replace("I", "L")
    return p


def peaks_accepted(p):
    """PEAKS's own filter: threshold at 1% FDR over distinct sequences (best PSM each)."""
    out = {}
    for f, g in p.groupby("Source File"):
        best = g.sort_values("-10LgP", ascending=False).drop_duplicates("seq")
        thr = threshold_at_fdr(best["-10LgP"].to_numpy(), best["decoy"].to_numpy(), 0.01)
        out[f] = (thr, set(g.loc[(g["-10LgP"] >= thr) & ~g["decoy"], "key"]))
    return out


def features(df):
    from peaksparity.ion_features import ION_FEATURE_COLUMNS
    extra = ["log_prot_support", "from_pass2", "caa", "caa_pct", "tag_len", "max_caa_any_beam", "denovo_conf_mean", "denovo_exact"]
    cols = [c for c in rf.FEATURE_COLUMNS + extra + ION_FEATURE_COLUMNS if c in df.columns]
    cols = list(dict.fromkeys(cols))
    return df[cols].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0), cols


def fit(df, y):
    X, cols = features(df)
    m = lgb.LGBMRegressor(n_estimators=1500, learning_rate=0.03, num_leaves=63, subsample=0.8,
                          colsample_bytree=0.8, min_child_samples=40, verbose=-1, random_state=0)
    m.fit(X, y)
    m.feature_cols_ = cols
    return m


def predict(m, df):
    X, _ = features(df)
    return m.predict(X.reindex(columns=m.feature_cols_, fill_value=0.0))  # absent features -> 0


def evaluate(df, pred, acc_by_file, peaks_thr_by_file):
    """Accept by PEAKS's rule on predicted scores (threshold from our decoys at 1% sequence FDR);
    also at PEAKS's own numeric threshold. Returns recall/precision vs PEAKS's accepted set."""
    d = df.copy(); d["score"] = pred
    top = best_per_spectrum(d, "score")
    top = top.sort_values("score", ascending=False).drop_duplicates(["source_file", "scan", "peptide_sage", "charge"])
    top["seq"] = top["stripped"].str.replace("I", "L")
    res = {}
    for f, g in top.groupby("source_file"):
        raw = f.replace(".mzML", ".raw")
        if raw not in acc_by_file:
            continue
        ref = acc_by_file[raw]
        best = g.sort_values("score", ascending=False).drop_duplicates("seq")
        thr = threshold_at_fdr(best["score"].to_numpy(), best["is_decoy"].to_numpy(), 0.01)
        for name, t in (("own_fdr_threshold", thr), ("peaks_threshold", peaks_thr_by_file[raw])):
            a = g[(g["score"] >= t) & ~g["is_decoy"]]
            sh = len(set(a["key"]) & ref)
            res[f"{raw}:{name}"] = dict(threshold=round(float(t), 2), n=len(a), n_ref=len(ref), shared=sh,
                                        recall=round(sh / max(len(ref), 1), 4), precision=round(sh / max(len(a), 1), 4))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks-psms", required=True)
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cv", action="store_true")
    ap.add_argument("--floor", type=float, default=FLOOR, help="PEAKS export score floor; unmatched candidates get floor-5")
    a = ap.parse_args()
    p = load_peaks(a.peaks_psms)
    acc = peaks_accepted(p)
    acc_by_file = {f: s for f, (t, s) in acc.items()}
    thr_by_file = {f: t for f, (t, s) in acc.items()}
    score_by_key = p[~p["decoy"]].groupby("key")["-10LgP"].max()
    parts = []
    for c in a.candidates:
        d = pd.read_parquet(c)
        ion = Path(c).with_name("ion_features.parquet")
        if ion.exists():  # optional extra fragment-series evidence (peaksparity.ion_features)
            f = pd.read_parquet(ion)
            if len(f) == len(d):
                for col in f.columns:
                    d[col] = f[col].to_numpy()
        d["modseq"] = d["peptide_sage"].map(canonical)
        d["key"] = key(d["scan"], d["modseq"], d["charge"])
        d["peaks_score"] = d["key"].map(score_by_key)
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    if "log_prot_support" not in rf.FEATURE_COLUMNS:
        rf.FEATURE_COLUMNS.append("log_prot_support")
    matched = df["peaks_score"].notna() & ~df["is_decoy"]
    print(f"candidates {len(df)}, matched to a PEAKS-listed PSM {int(matched.sum())}, PEAKS listed targets {int((~p['decoy']).sum())}")
    y = df["peaks_score"].fillna(a.floor - 5.0).to_numpy()   # censored: below PEAKS's export floor
    y = np.where(df["is_decoy"], a.floor - 5.0, y)
    report = {"n_candidates": int(len(df)), "n_matched": int(matched.sum())}
    if a.cv:
        pred = np.zeros(len(df))
        for fold in (0, 1):
            tr = (df["scan"] % 2 != fold).to_numpy(); te = ~tr
            m = fit(df[tr], y[tr]); pred[te] = predict(m, df[te])
        r = pred[matched.to_numpy()]; t = y[matched.to_numpy()]
        report["cv_matched_pearson"] = round(float(np.corrcoef(r, t)[0, 1]), 4)
        report["cv_matched_mae"] = round(float(np.abs(r - t).mean()), 3)
        report["cv_eval"] = evaluate(df, pred, acc_by_file, thr_by_file)
        print(json.dumps(report, indent=1))
    m = fit(df, y)
    m.peaks_thresholds_ = thr_by_file
    pickle.dump(m, open(a.out, "wb"))
    imp = pd.Series(m.feature_importances_, index=m.feature_cols_).sort_values(ascending=False)
    report["feature_importance"] = imp.head(15).to_dict()
    Path(a.out).with_suffix(".json").write_text(json.dumps(report, indent=1))
    print("saved", a.out)


if __name__ == "__main__":
    main()
