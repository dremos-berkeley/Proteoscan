"""Validate the de novo agreement features against a PEAKS export.

    python scripts/eval_denovo_features.py --peaks-export <export dir with db.psms.csv> \
        --candidates <run>/candidates_features_2pass_denovo.parquet [...] \
        [--denovo <run>.*.denovo_all.parquet ...] --out report.json

Reports
* Spearman between our features and PEAKS's per-PSM columns (Tag Length, CAA (%), mean
  Positional Confidence, -10LgP) on candidates that PEAKS lists (same file, scan, modified
  sequence, charge; I = L);
* feature distributions for PEAKS-listed targets, non-listed targets and decoys;
* the same features computed from PEAKS's own de novo calls (<file>.denovo.csv, best ALC per
  scan and charge), and the agreement of our de novo calls with PEAKS's.
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from peaksparity.denovo_features import DENOVO_FEATURE_COLUMNS, denovo_features
from peaksparity.seqnorm import canonical

COLS = ["source_file", "scan", "query_id", "peptide_sage", "charge", "is_decoy", "rank", "peptide_len"]


def _key(f, scan, modseq, z):
    return (f.astype(str).str.replace(r"\.(mzML|raw)$", "", regex=True) + "|" + scan.astype(str) + "|"
            + modseq.astype(str) + "|" + z.astype(str)).str.replace("I", "L")


def _sp(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return None
    return round(float(spearmanr(a[m], b[m]).statistic), 4)


def summary(df, cols):
    out = {}
    for c in cols:
        v = df[c].astype(float)
        out[c] = {"mean": round(float(v.mean()), 3), "median": round(float(v.median()), 3),
                  "frac_pos": round(float((v > 0).mean()), 4)}
    out["n"] = int(len(df))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks-export", required=True)
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--denovo", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    exp = Path(a.peaks_export)
    p = pd.read_csv(exp / "db.psms.csv")
    p["decoy"] = p["decoy"].astype(str).str.lower().eq("true") if "decoy" in p else False
    p["key"] = _key(p["Source File"], p["Scan"], p["Peptide"], p["z"])
    p["posconf_mean"] = p["Positional Confidence"].astype(str).str.split().map(
        lambda v: float(np.mean(np.asarray(v, dtype=float))) if len(v) else np.nan)
    pt = p[~p["decoy"]].sort_values("-10LgP", ascending=False).drop_duplicates("key").set_index("key")
    rep = {"peaks_rows": int(len(p)), "peaks_targets": int((~p["decoy"]).sum()),
           "peaks_tag_length_equals_length": float((p["Tag Length"] == p["Length"]).mean())}
    # PEAKS's CAA (%) = % residues whose Positional Confidence >= the report threshold (2 %)
    est = p["Positional Confidence"].astype(str).str.split().map(
        lambda v: 100 * float(np.mean(np.asarray(v, dtype=float) >= 2)))
    rep["peaks_caa_pct_equals_posconf_ge2_pct"] = float((np.abs(est.round(1) - p["CAA (%)"]) < 0.06).mean())

    parts = []
    for c in a.candidates:
        cols = COLS + [x for x in DENOVO_FEATURE_COLUMNS]
        parts.append(pd.read_parquet(c, columns=cols))
    d = pd.concat(parts, ignore_index=True)
    d["modseq"] = d["peptide_sage"].map(canonical)
    d["key"] = _key(d["source_file"], d["scan"], d["modseq"], d["charge"])
    d["listed"] = d["key"].isin(pt.index) & ~d["is_decoy"]
    rep["candidates"] = int(len(d))
    rep["with_denovo"] = round(float(d["caa"].notna().mean()), 4)
    rep["peaks_targets_matched"] = int(d.loc[d["listed"], "key"].nunique())

    m = d[d["listed"]].drop_duplicates("key").join(pt[["Tag Length", "CAA (%)", "posconf_mean", "-10LgP", "Length"]], on="key")
    rep["spearman_listed"] = {
        "tag_len~Tag Length": _sp(m["tag_len"].to_numpy(float), m["Tag Length"].to_numpy(float)),
        "caa_pct~CAA (%)": _sp(m["caa_pct"].to_numpy(float), m["CAA (%)"].to_numpy(float)),
        "caa~CAA (%)": _sp(m["caa"].to_numpy(float), m["CAA (%)"].to_numpy(float)),
        "denovo_conf_mean~mean Positional Confidence": _sp(m["denovo_conf_mean"].to_numpy(float), m["posconf_mean"].to_numpy(float)),
        "caa_pct~-10LgP": _sp(m["caa_pct"].to_numpy(float), m["-10LgP"].to_numpy(float)),
        "tag_len~-10LgP": _sp(m["tag_len"].to_numpy(float), m["-10LgP"].to_numpy(float)),
        "max_caa_any_beam~-10LgP": _sp(m["max_caa_any_beam"].to_numpy(float), m["-10LgP"].to_numpy(float)),
        "CAA (%)~-10LgP (PEAKS vs PEAKS)": _sp(m["CAA (%)"].to_numpy(float), m["-10LgP"].to_numpy(float)),
        "n": int(len(m)),
    }
    feats = DENOVO_FEATURE_COLUMNS
    groups = {"peaks_listed_target": d[d["listed"]], "unlisted_target": d[~d["listed"] & ~d["is_decoy"]],
              "decoy": d[d["is_decoy"]]}
    rep["distributions_all_candidates"] = {g: summary(x.dropna(subset=["caa"]), feats) for g, x in groups.items()}
    r1 = d["rank"] == 1
    rep["distributions_rank1"] = {g: summary(x[r1.loc[x.index]].dropna(subset=["caa"]), feats) for g, x in groups.items()}

    # PEAKS's own de novo calls: the same features against PEAKS-listed PSMs, and agreement with ours
    pdn = []
    for f in sorted(exp.glob("*.raw.denovo.csv")):
        x = pd.read_csv(f, usecols=["Scan", "Peptide", "ALC (%)", "z", "local confidence (%)"])
        x["source_file"] = f.name.replace(".raw.denovo.csv", "")
        pdn.append(x)
    if pdn:
        x = pd.concat(pdn, ignore_index=True).sort_values("ALC (%)", ascending=False)
        x["sk"] = x["source_file"] + "|" + x["Scan"].astype(str) + "|" + x["z"].astype(str)
        x = x.drop_duplicates("sk")
        x["sequence"] = x["Peptide"]
        x["per_residue_scores"] = x["local confidence (%)"].astype(str).str.split().map(
            lambda v: np.asarray(v, dtype=float) / 100.0)
        x["picked"] = True
        mm = m.copy()
        mm["sk"] = (mm["source_file"].str.replace(r"\.mzML$", "", regex=True) + "|" + mm["scan"].astype(str)
                    + "|" + mm["charge"].astype(str))
        f2 = denovo_features(mm, x[["sk", "sequence", "per_residue_scores", "picked"]], key="sk")
        mm = mm.join(f2.add_prefix("peaksdn_"))
        has = mm["peaksdn_caa"].notna()
        rep["peaks_denovo_on_listed"] = {
            "with_peaks_denovo": round(float(has.mean()), 4),
            "peaksdn_caa_pct_mean": round(float(mm.loc[has, "peaksdn_caa_pct"].mean()), 2),
            "our_caa_pct_mean_same_rows": round(float(mm.loc[has, "caa_pct"].mean()), 2),
            "peaksdn_exact": round(float(mm.loc[has, "peaksdn_denovo_exact"].mean()), 4),
            "our_exact_same_rows": round(float(mm.loc[has, "denovo_exact"].mean()), 4),
            "spearman_peaksdn_caa_pct~our_caa_pct": _sp(mm.loc[has, "peaksdn_caa_pct"].to_numpy(float), mm.loc[has, "caa_pct"].to_numpy(float)),
            "spearman_peaksdn_caa_pct~CAA (%)": _sp(mm.loc[has, "peaksdn_caa_pct"].to_numpy(float), mm.loc[has, "CAA (%)"].to_numpy(float)),
            "spearman_peaksdn_caa_pct~-10LgP": _sp(mm.loc[has, "peaksdn_caa_pct"].to_numpy(float), mm.loc[has, "-10LgP"].to_numpy(float)),
            "spearman_peaksdn_tag_len~-10LgP": _sp(mm.loc[has, "peaksdn_tag_len"].to_numpy(float), mm.loc[has, "-10LgP"].to_numpy(float)),
        }
        if a.denovo:
            ours = pd.concat([pd.read_parquet(f, columns=["query_id", "scan", "charge", "peptide", "picked"])
                              .assign(source_file=Path(f).name.split(".")[0]) for f in a.denovo], ignore_index=True)
            ours = ours[ours["picked"]]
            ours["sk"] = ours["source_file"] + "|" + ours["scan"].astype(str) + "|" + ours["charge"].astype(str)
            j = ours.merge(x[["sk", "Peptide", "ALC (%)"]], on="sk")
            same = j["peptide"].fillna("").str.replace("I", "L") == j["Peptide"].str.replace("I", "L")
            rep["our_vs_peaks_denovo"] = {
                "queries_compared": int(len(j)), "exact_IL": round(float(same.mean()), 4),
                "exact_IL_peaks_alc_ge80": round(float(same[j["ALC (%)"] >= 80].mean()), 4),
                "n_peaks_alc_ge80": int((j["ALC (%)"] >= 80).sum())}
    Path(a.out).write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
