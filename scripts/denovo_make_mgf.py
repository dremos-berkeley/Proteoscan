"""Build Casanovo annotated MGFs (SEQ= lines) from PEAKS-format PSM exports.

Example (HeLa PXD004940 ion-trap runs; one whole run held out for validation):

    python scripts/denovo_make_mgf.py --out /drive_a/drew/peaks-parity-work/denovo_finetune/mgf \
        --psms-root /drive_a/drew/peaks-parity-work/hela --mzml-root /drive_a/drew/public-data/PXD004940/mzml \
        --pattern '*_HCD_IT_*' --valid 151118_Q_Hela_03_03_50cm_HCD_IT_EasySpray_90min_500ng

    # test set: PEAKS PSMs of another run, every PSM (no score filter, no peptide dedup)
    python scripts/denovo_make_mgf.py --out ... --test-psms db.psms.csv --test-mzml run.mzML --test-name jurkat

Each run's PSMs: -10LgP >= --min-score, vocabulary-compatible mods, best PSM per (peptide, z).
Precursor m/z and charge are the PSM's (feature-corrected) values, not the raw mzML ones.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaksparity import denovo  # noqa: E402


def _one(args):
    name, psms, mzml, out, min_score, unique = args
    ent = denovo.training_psms(psms, min_score=min_score, unique_peptide_charge=unique)
    n = denovo.write_annotated_mgf(ent, mzml, out, title_prefix=name)
    ent.to_csv(Path(out).with_suffix(".labels.tsv"), sep="\t", index=False)
    return name, len(ent), n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--psms-root", help="dir with <run>/export/db.psms.csv")
    ap.add_argument("--mzml-root", help="dir with <run>.mzML")
    ap.add_argument("--pattern", default="*")
    ap.add_argument("--valid", action="append", default=[], help="run name(s) held out for validation")
    ap.add_argument("--min-score", type=float, default=30.0)
    ap.add_argument("--test-psms")
    ap.add_argument("--test-mzml")
    ap.add_argument("--test-name", default="test")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tasks = []
    if a.psms_root:
        for d in sorted(Path(a.psms_root).glob(a.pattern)):
            psms = d / "export" / "db.psms.csv"
            mzml = Path(a.mzml_root) / f"{d.name}.mzML"
            if not d.is_dir() or not psms.exists() or not mzml.exists():
                continue
            split = "valid" if d.name in a.valid else "train"
            tasks.append((d.name, str(psms), str(mzml), str(out / f"{split}.{d.name}.mgf"), a.min_score, True))
    if a.test_psms:
        tasks.append((a.test_name, a.test_psms, a.test_mzml, str(out / f"test.{a.test_name}.mgf"), -1e9, False))
    summary = {}
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for name, n_lab, n_written in ex.map(_one, tasks):
            summary[name] = {"labels": n_lab, "spectra_written": n_written}
            print(name, n_lab, n_written, flush=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
