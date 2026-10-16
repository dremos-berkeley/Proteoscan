"""Re-score the cached Jurkat pipeline run (DrewR_Adam_91526_DDA_OT_4) with other fragment predictors and
compare each export with the facility PEAKS export (proteoscan.concordance). docs/IONTRAP_MODEL.md.

    python scripts/iontrap_jurkat.py --variant ms2pip=cached --variant ft=peptdeep:<ckpt> --out-root <dir>
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.iontrap_benchmark import run_one  # noqa: E402

SRC = "/drive_a/drew/proteoscan-work/benchmarks/pipeline_jurkat"
MZML = "/drive_a/drew/openms-data/260921-DrewR_Adam_91526/mzml/DrewR_Adam_91526_DDA_OT_4.mzML"
REF = "/drive_a/drew/vendor-exports/DrewR_Adam_91526_DDA_OT/Analysis 1 - Sample 4"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", action="append", required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0], help="mokapot seeds (replicate rescorings)")
    a = ap.parse_args()
    rows = []
    for v in a.variant:
      for sd in a.seeds:
        name0, spec = v.split("=", 1)
        name = name0 if sd == 0 else f"{name0}_seed{sd}"
        out = Path(a.out_root) / name
        run_one(dict(pool=None, run="DrewR_Adam_91526_DDA_OT_4", type="Sample 4", variant=name, spec=spec, src=SRC,
                     mzml=MZML, out=str(out), ref=REF, seed=sd,
                     features_from=None if sd == 0 else str(Path(a.out_root) / name0)))
        st = json.loads((out / "run_stats.json").read_text())
        m = json.loads((out / "concordance" / "metrics.json").read_text())
        h = m["headline"]
        rows.append(dict(variant=name0, seed=sd, psm_rows=st["n_psm_rows"], peptides=st["n_peptides"],
                         protein_groups=st["n_protein_groups"], peptide_threshold=round(st["peptide_threshold"], 2),
                         **{k: h[k] for k in h if any(x in k for x in ("recall", "precision", "shared", "n_ref", "n_test"))}))
        for lvl in ("peptide", "protein_group"):
            x = m.get(lvl, {})
            rows[-1][f"{lvl}_recall"] = x.get("recall"); rows[-1][f"{lvl}_precision"] = x.get("precision")
        print(rows[-1], flush=True)
    import pandas as pd
    pd.DataFrame(rows).to_csv(Path(a.out_root) / "jurkat_summary.csv", index=False)
    print(pd.DataFrame(rows).T.to_string())


if __name__ == "__main__":
    main()
