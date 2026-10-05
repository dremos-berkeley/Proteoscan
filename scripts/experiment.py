"""Run rescoring + export + concordance for an existing Comet result (one run).

    python scripts/experiment.py --comet-txt X.txt --mzml run.mzML --config cfg.yml \
        --ref <PEAKS export dir> --out <dir> [--entries precursors.parquet] [--sample "Sample 4"]

Writes features.parquet (cached), export_tda/, export_imitate_cv/, concordance_*/ and summary.json.
The imitation model is cross-validated by scan parity so it is never evaluated on its training scans.
"""
import argparse, json, logging, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
from proteoscan.comet_search import fixed_mod_deltas, read_comet
from proteoscan.config import SearchConfig
from proteoscan.export_builder import apply_protein_mapping, assemble_export, build_peptides, build_psms
from proteoscan.peaks_export import write_export
from proteoscan.pipeline import attach_entries, peaks_filter
from proteoscan.proteins import infer_proteins
from proteoscan.rescore import apply_imitation, best_per_spectrum, minus10lgp, rescore_tda, train_imitation
from proteoscan.rescore_features import compute_features
from proteoscan.seqnorm import canonical


def export_and_compare(scored, cfg, sample, out, name, ref, fixed_count=None):
    top = best_per_spectrum(scored, "score")
    top["minus10lgp"] = minus10lgp(top["score"], top["is_decoy"])
    if fixed_count:
        t = top[~top["is_decoy"]].sort_values("score", ascending=False).drop_duplicates(["source_file", "scan", "peptide_sage", "charge"]).head(fixed_count)
        acc, st = t.copy(), dict(fixed_count=fixed_count)
    else:
        acc, st = peaks_filter(top, cfg)
    if "entry_rt" in acc:
        acc["rt"] = acc["entry_rt"].fillna(acc["rt"])
    psms = apply_protein_mapping(build_psms(acc, cfg, sample), cfg.database)
    peps = build_peptides(psms, sample)
    prot, pp = infer_proteins(peps, cfg.database, score_threshold=cfg.protein_score_min,
                              min_unique=cfg.protein_unique_peptides_min, psms_df=psms, sample=sample)
    write_export(assemble_export(psms, peps, prot, pp, sample_name=sample, params=cfg.to_dict()), str(out / f"export_{name}"))
    from proteoscan.concordance import main as conc
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        conc(["--ref", ref, "--test", str(out / f"export_{name}"), "--out", str(out / f"concordance_{name}")])
    m = json.load(open(out / f"concordance_{name}" / "metrics.json"))
    lv = m["levels"]
    res = {k: {x: round(lv[k][x], 4) if isinstance(lv[k][x], float) else lv[k][x] for x in ("n_ref", "n_test", "shared", "recall", "precision")}
           for k in ("psm", "peptide", "protein_group")}
    res["filter"] = st
    res["scan_agreement"] = m.get("scan_agreement")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comet-txt", required=True); ap.add_argument("--mzml", required=True)
    ap.add_argument("--config", required=True); ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--entries"); ap.add_argument("--sample", default="Sample 4")
    ap.add_argument("--pass2-txt", help="Comet result of the restricted second-pass search"); ap.add_argument("--pass1-dir", help="experiment dir of pass 1 (features.parquet)")
    ap.add_argument("--candidates", type=int, default=3); ap.add_argument("--processes", type=int, default=24)
    a = ap.parse_args()
    logging.basicConfig(level=logging.WARNING)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cfg = SearchConfig.from_yaml(a.config)
    fixed = fixed_mod_deltas(cfg)
    t0 = time.time()
    if (out / "features.parquet").exists():
        feat = pd.read_parquet(out / "features.parquet")
    else:
        c = read_comet(a.comet_txt, source_file=Path(a.mzml).name, fixed_mods=fixed)
        c = c[c["rank"] <= a.candidates].reset_index(drop=True)
        c["proteins"] = c["proteins"].str.replace("DECOY_", "#DECOY#", regex=False)
        c["query_id"] = c["scan"]
        if a.pass2_txt:
            from proteoscan.twopass import merge_passes
            c2 = read_comet(a.pass2_txt, source_file=Path(a.mzml).name, fixed_mods=fixed)
            c2["proteins"] = c2["proteins"].str.replace("DECOY_", "#DECOY#", regex=False)
            c2["query_id"] = c2["scan"]
            c = merge_passes(c, c2, keep=a.candidates + 2)
        if a.entries:
            c = attach_entries(c, pd.read_parquet(a.entries))
            c["rt"] = pd.Series(c["raw_rt"]).fillna(c["rt"]).to_numpy()
            best = c.groupby(["source_file", "scan", "precursor_idx"])["hyperscore"].transform("max")
            c["delta_best"] = best - c["hyperscore"]
        feat = compute_features(c, a.mzml, ms2pip_model="HCD2021", ms2_tolerance=cfg.fragment_mass_error_tolerance_da, processes=a.processes)
        feat.to_parquet(out / "features.parquet")
    if a.pass1_dir:
        from proteoscan.twopass import add_protein_features, protein_support
        p1 = pd.read_parquet(Path(a.pass1_dir) / "features.parquet")
        p1s = rescore_tda(p1); t1 = best_per_spectrum(p1s, "score"); t1["minus10lgp"] = minus10lgp(t1["score"], t1["is_decoy"])
        acc1, _ = peaks_filter(t1, cfg)
        support, seqs = protein_support(acc1)
        feat = add_protein_features(feat, support, seqs)
        import proteoscan.rescore_features as rf
        if "log_prot_support" not in rf.FEATURE_COLUMNS:
            rf.FEATURE_COLUMNS.append("log_prot_support")
    peaks = pd.read_csv(Path(a.ref) / "db.psms.csv")
    pk = set(peaks["Scan"].astype(int).astype(str) + "|" + peaks["Peptide"] + "|" + peaks["z"].astype(str))
    feat["modseq"] = feat["peptide_sage"].map(canonical)
    # PEAKS lists each I/L variant; label on the I=L collapsed key
    pk_il = {k.replace("I", "L") for k in pk}
    feat["key_il"] = (feat["scan"].astype(str) + "|" + feat["modseq"] + "|" + feat["charge"].astype(str)).str.replace("I", "L")
    feat["peaks_label"] = feat["key_il"].isin(pk_il).astype(int)
    summary = dict(n_candidates=len(feat), peaks_psm_keys=len(pk_il),
                   candidate_ceiling=round(len(set(feat["key_il"]) & pk_il) / len(pk_il), 4),
                   rank1_ceiling=round(len(set(feat.loc[feat["rank"] == 1, "key_il"]) & pk_il) / len(pk_il), 4),
                   features_s=round(time.time() - t0))
    print(json.dumps(summary), flush=True)
    tda = rescore_tda(feat)
    summary["tda"] = export_and_compare(tda, cfg, a.sample, out, "tda", a.ref)
    summary["tda_at_peaks_count"] = export_and_compare(tda, cfg, a.sample, out, "tda_n", a.ref, fixed_count=len(peaks))
    print(json.dumps(summary["tda"]), flush=True)
    parts = []
    for fold in (0, 1):
        tr = feat[feat["scan"] % 2 != fold]; te = feat[feat["scan"] % 2 == fold]
        model = train_imitation(tr, (tr["peaks_label"] == 1) & (~tr["is_decoy"]))
        parts.append(apply_imitation(model, te))
    imit = pd.concat(parts, ignore_index=True)
    summary["imitate_cv"] = export_and_compare(imit, cfg, a.sample, out, "imitate_cv", a.ref)
    summary["imitate_cv_at_peaks_count"] = export_and_compare(imit, cfg, a.sample, out, "imitate_cv_n", a.ref, fixed_count=len(peaks))
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1)); print("DONE")


if __name__ == "__main__":
    main()
