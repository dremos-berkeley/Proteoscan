"""Build fragment-intensity training/evaluation tables for the ion-trap MS2 model.

    python scripts/iontrap_dataset.py --labels <runs>/*/labels.parquet --split train --out data/pt_train
    python scripts/iontrap_dataset.py --peaks-psms <export>/db.psms.csv --mzml <run.mzML> --min-score 30 \
        --split train --source hela --out data/hela_<run>

Output: <out>.parquet (one row per spectrum: sequence, mods, mod_sites, charge, nce, instrument,
activation, scan_mode, source, run, scan, split, nAA, frag_start_idx, frag_stop_idx) and <out>.npy
(float32 matrix, one row per fragment position, columns fragpred.DEFAULT_FRAG_TYPES). Intensities are
max-normalised per spectrum over the observable ions; unobservable ions (fragment charge > precursor
charge, or m/z outside the scan window) are -1.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # single-threaded numpy/scipy: many runs go in parallel

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan.fragpred import DEFAULT_FRAG_TYPES, annotate, fragment_mz, parse_peptide, valid_mask  # noqa: E402
from scripts.iontrap_label import ms2_table  # noqa: E402

BENCH = set(Path("/drive_a/drew/public-data/PXD004732/pools.txt").read_text().split())
INSTRUMENT = "Lumos"


def extract(rows: pd.DataFrame, mzml: str, tol: float, frag_types=DEFAULT_FRAG_TYPES):
    """rows: modseq, charge, scan (+ any metadata). Returns (meta df, list of matrices)."""
    od, ms2 = ms2_table(Path(mzml))
    ms2 = ms2.set_index("scan")
    metas, mats = [], []
    for r in rows.itertuples(index=False):
        if r.scan not in ms2.index:
            continue
        info = ms2.loc[r.scan]
        seq, res, mods, sites, ok = parse_peptide(r.modseq)
        if not ok or len(seq) < 5 or len(seq) > 50:
            continue
        mz, it = od.getSpectrum(int(info.idx)).get_peaks()
        mz = np.asarray(mz, float); it = np.asarray(it, float)
        o = np.argsort(mz); mz, it = mz[o], it[o]
        theo = fragment_mz(res, frag_types)
        obs = annotate(theo, mz, it, tol)
        valid = valid_mask(theo, frag_types, int(r.charge), info.lo_mz, info.hi_mz)
        mx = obs[valid].max() if valid.any() else 0.0
        if mx <= 0:
            continue
        mat = np.where(valid, obs / mx, -1.0).astype(np.float32)
        d = r._asdict()
        d.update(sequence=seq, mods=mods, mod_sites=sites, nce=float(info.nce), activation=info.activation,
                 scan_mode=info.scan_mode, instrument=INSTRUMENT, precursor_mz=float(info.precursor_mz),
                 tic=float(it.sum()), mzml=str(mzml))
        metas.append(d); mats.append(mat)
    return pd.DataFrame(metas), mats


def write(meta: pd.DataFrame, mats: list[np.ndarray], out: str):
    meta = meta.reset_index(drop=True)
    meta["nAA"] = meta["sequence"].str.len()
    lens = meta["nAA"].to_numpy() - 1
    stop = np.cumsum(lens); start = stop - lens
    meta["frag_start_idx"] = start; meta["frag_stop_idx"] = stop
    np.save(out + ".npy", np.concatenate(mats) if mats else np.zeros((0, len(DEFAULT_FRAG_TYPES)), np.float32))
    meta.to_parquet(out + ".parquet")
    print(f"{out}: {len(meta)} spectra, {int(stop[-1]) if len(stop) else 0} fragment rows")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="*", help="labels.parquet files from iontrap_label.py")
    ap.add_argument("--peaks-psms", help="PEAKS-format db.psms.csv (ours or PEAKS's)")
    ap.add_argument("--mzml", help="mzML for --peaks-psms")
    ap.add_argument("--min-score", type=float, default=30.0)
    ap.add_argument("--best-per-peptide", action="store_true", help="keep the best PSM per (peptide, charge)")
    ap.add_argument("--source", default="proteometools")
    ap.add_argument("--split", required=True)
    ap.add_argument("--allow-benchmark", action="store_true", help="only for test-set evaluation tables")
    ap.add_argument("--tol", type=float, default=0.4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    metas, mats = [], []
    if a.labels:
        for f in a.labels:
            lab = pd.read_parquet(f)
            if not len(lab):
                continue
            pool = lab["pool"].iloc[0]
            if pool in BENCH and not a.allow_benchmark:
                raise SystemExit(f"refusing benchmark pool {pool} in a training table")
            for mzml, g in lab.groupby("mzml"):
                g = g.rename(columns={"nce": "label_nce"})[["pool", "run", "scan_type", "scan", "modseq", "charge", "score", "n_by"]]
                m, x = extract(g, mzml, a.tol)
                m["source"] = a.source; m["split"] = a.split
                metas.append(m); mats += x
                print(f, mzml.split("/")[-1], len(m), flush=True)
    if a.peaks_psms:
        p = pd.read_csv(a.peaks_psms)
        if "decoy" in p:
            p = p[~p["decoy"].astype(str).str.lower().eq("true")]
        p = p[p["-10LgP"] >= a.min_score]
        if a.best_per_peptide:
            p = p.sort_values("-10LgP", ascending=False).drop_duplicates(["Peptide", "z"])
        g = pd.DataFrame({"run": Path(a.mzml).stem, "scan": p["Scan"].astype(int), "modseq": p["Peptide"],
                          "charge": p["z"].astype(int), "score": p["-10LgP"]}).drop_duplicates(["scan", "modseq", "charge"])
        m, x = extract(g, a.mzml, a.tol)
        m["source"] = a.source; m["split"] = a.split; m["scan_type"] = "hcd_it"
        metas.append(m); mats += x
    write(pd.concat(metas, ignore_index=True), mats, a.out)


if __name__ == "__main__":
    main()
