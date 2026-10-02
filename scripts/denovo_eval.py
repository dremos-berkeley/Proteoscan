"""Evaluate Casanovo checkpoints against PEAKS identifications of a test run (Jurkat).

Modes (each writes <out>/<model>.<mode>[.bN].{mztab,json,curve.tsv}):

  peaks    Sequence the test MGF built by denovo_make_mgf.py (PEAKS PSM scans with PEAKS's
           corrected precursor m/z and z; one spectrum per (scan, peptide, z)) and compare with
           the PEAKS DB peptide.
  mzml     Sequence the whole mzML with its own precursor records; compare on the scans PEAKS
           identified (best PEAKS PSM per scan).
  denovo   Sequence PEAKS's denovoOnly spectra (PEAKS m/z and z per row) and measure agreement
           with PEAKS's own de novo peptide (exact and I=L), overall and by PEAKS ALC.

With --top-match > 1 the highest-scoring candidate within 15 ppm of the precursor is chosen.

Example:
  python scripts/denovo_eval.py --model released=~/.cache/casanovo/casanovo_orbitrap_v5-2-0.ckpt \
      --model ft=/drive_a/.../best.ckpt --out /drive_a/.../eval --gpu 3 --modes peaks denovo \
      --test-mgf .../mgf/test.jurkat_OT4.mgf --denovo-only '.../Sample 4.db.denovoOnly.csv' \
      --mzml .../DrewR_Adam_91526_DDA_OT_4.mzML
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan import denovo  # noqa: E402

log = logging.getLogger("denovo_eval")


def _summ(res):
    out = {k: v for k, v in res.items() if k not in ("curve", "per_spectrum")}
    out["curve_points"] = denovo.curve_points(res["curve"]).to_dict(orient="records")
    return out


def _run(mgf_or_mzml, ckpt, mztab, a, **kw):
    return denovo.sequence_mzml(mgf_or_mzml, ckpt, mztab, n_beams=a.n_beams, top_match=a.top_match,
                                gpu=a.gpu, batch_size=a.batch_size, reuse=not a.rerun, config=a.config, **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", action="append", required=True, help="name=ckpt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--modes", nargs="+", default=["peaks"])
    ap.add_argument("--test-mgf")
    ap.add_argument("--test-labels", help="default: <test-mgf stem>.labels.tsv")
    ap.add_argument("--mzml")
    ap.add_argument("--peaks-psms", help="PEAKS db.psms.csv (mode mzml)")
    ap.add_argument("--denovo-only", help="PEAKS <Sample>.db.denovoOnly.csv (mode denovo)")
    ap.add_argument("--mzml-mztab", action="append", default=[],
                    help="name=existing whole-mzML mztab to reuse (mode mzml)")
    ap.add_argument("--n-beams", type=int, default=1)
    ap.add_argument("--top-match", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=1024)
    ap.add_argument("--gpu", default=None)
    ap.add_argument("--rerun", action="store_true")
    ap.add_argument("--config", help="Casanovo YAML (e.g. preprocessing overrides such as min_mz)")
    ap.add_argument("--tag", default="", help="suffix for output names")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = "" if (a.n_beams == 1 and a.top_match == 1) else f".b{a.n_beams}t{a.top_match}"
    tag += a.tag
    existing = dict(x.split("=", 1) for x in a.mzml_mztab)
    for spec in a.model:
        name, ckpt = spec.split("=", 1)
        results = {}
        if "peaks" in a.modes:
            labels = pd.read_csv(a.test_labels or Path(a.test_mgf).with_suffix(".labels.tsv"), sep="\t")
            labels["spectra_ref"] = [f"index={i}" for i in range(len(labels))]
            pred = _run(a.test_mgf, ckpt, out / f"{name}.peaks{tag}.mztab", a)
            res = denovo.evaluate(pred, labels, key="spectra_ref")
            res["curve"].to_csv(out / f"{name}.peaks{tag}.curve.tsv", sep="\t", index=False)
            res["per_spectrum"].to_csv(out / f"{name}.peaks{tag}.per_spectrum.tsv", sep="\t", index=False)
            s = _summ(res)
            s["precursor_ok_rate"] = float(pred["precursor_ok"].mean())
            results["peaks"] = s
        if "mzml" in a.modes:
            psms = pd.read_csv(a.peaks_psms)
            psms = psms.sort_values("-10LgP", ascending=False).drop_duplicates("Scan")
            truth = pd.DataFrame({"scan": psms["Scan"].astype(int), "peptide": psms["Peptide"]})
            truth["seq"] = truth["peptide"].map(denovo.to_casanovo)
            truth = truth.dropna(subset=["seq"])
            if name in existing and not tag:
                pred = denovo.read_mztab(existing[name])
                pred["ppm"] = [denovo.precursor_ppm(s_, m, int(z)) for s_, m, z in
                               zip(pred["sequence"], pred["exp_mz"], pred["charge"])]
                pred["precursor_ok"] = pred["ppm"] <= 15
            else:
                pred = _run(a.mzml, ckpt, out / f"{name}.mzml{tag}.mztab", a)
            pred["scan"] = pred["scan"].astype(int)
            res = denovo.evaluate(pred.drop_duplicates("scan"), truth, key="scan")
            res["curve"].to_csv(out / f"{name}.mzml{tag}.curve.tsv", sep="\t", index=False)
            s = _summ(res)
            s["n_spectra_sequenced"] = int(len(pred))
            results["mzml"] = s
        if "denovo" in a.modes:
            dn = pd.read_csv(a.denovo_only, keep_default_na=False)
            dn["seq"] = dn["Peptide"].map(denovo.to_casanovo)
            ent = pd.DataFrame({"scan": dn["Scan"].astype(int), "charge": dn["z"].astype(int),
                                "mz": dn["m/z"].astype(float)})
            mgf = out / "denovo_only.input.mgf"
            if not mgf.exists():
                n = denovo.write_annotated_mgf(ent, a.mzml, mgf)
                assert n == len(ent), (n, len(ent))
            pred = _run(mgf, ckpt, out / f"{name}.denovo{tag}.mztab", a)
            truth = pd.DataFrame({"spectra_ref": [f"index={i}" for i in range(len(dn))],
                                  "seq": dn["seq"], "alc": dn["ALC (%)"].astype(float)})
            # charges > 4 are skipped by Casanovo (max_charge 4); they stay in the denominator
            truth = truth.dropna(subset=["seq"])
            res = denovo.evaluate(pred, truth, key="spectra_ref")
            ps = res["per_spectrum"].merge(truth[["spectra_ref", "alc"]], on="spectra_ref")
            ps.to_csv(out / f"{name}.denovo{tag}.per_spectrum.tsv", sep="\t", index=False)
            bins = pd.cut(ps["alc"], [50, 60, 70, 80, 90, 100.01], right=False)
            by = ps.groupby(bins, observed=True).agg(n=("exact_il", "size"), agree_exact=("exact", "mean"),
                                                     agree_il=("exact_il", "mean"),
                                                     casanovo_pep=("casanovo_pep", "mean"))
            s = _summ(res)
            s["by_peaks_alc"] = {str(k): v for k, v in by.to_dict(orient="index").items()}
            # our ALC (mean per-residue score) distribution on these spectra
            alc = pred["per_residue_scores"].map(lambda x: float(np.mean(x)) * 100 if len(x) else 0)
            s["our_alc_ge50_rate"] = float((alc >= 50).mean())
            results["denovo"] = s
        (out / f"{name}{tag}.json").write_text(json.dumps(results, indent=1, default=float))
        print(name, json.dumps({m: {k: v for k, v in r.items() if k in ("n", "exact", "exact_il", "casanovo_pep",
                                                                        "aa_precision", "aa_recall")}
                                for m, r in results.items()}, default=float), flush=True)


if __name__ == "__main__":
    main()
