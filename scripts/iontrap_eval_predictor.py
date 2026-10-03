"""Evaluate fragment-intensity predictors on labelled ion-trap spectra (docs/IONTRAP_MODEL.md).

    python scripts/iontrap_eval_predictor.py --table data/pt_val --name pt_val --activation hcd \
        --predictor ms2pip:HCD2021 --predictor peptdeep:released --predictor peptdeep:models/hcd_it_v1.pth \
        --out eval/pt_val_hcd.csv

Per spectrum: spectral angle (SA = 1 - 2*angle/pi on max-normalised linear intensities) and Pearson
correlation, over the observable b/y 1+ ions (common to all predictors); for peptdeep models also over
every fragment type the model predicts ("all"). Each spectrum is predicted at its own NCE (filter string).
Writes per-spectrum values (<out>.spectra.parquet) and a summary row per predictor (<out>).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan.fragpred import DEFAULT_FRAG_TYPES, PeptdeepPredictor  # noqa: E402
from scripts.iontrap_train import benchmark_sequences, spectrum_metrics  # noqa: E402

BY1 = [DEFAULT_FRAG_TYPES.index("b_z1"), DEFAULT_FRAG_TYPES.index("y_z1")]


def load(prefix):
    m = pd.read_parquet(prefix + ".parquet"); X = np.load(prefix + ".npy")
    return m, X


def proforma(seq, mods, sites):
    names = {}
    if isinstance(mods, str) and mods:
        for nm, st in zip(mods.split(";"), str(sites).split(";")):
            names[int(st)] = nm.split("@")[0]
    s = "".join(a + (f"[{names[i + 1]}]" if i + 1 in names else "") for i, a in enumerate(seq))
    return (f"[{names[0]}]-" if 0 in names else "") + s


def ms2pip_predict(meta, model, processes=16):
    from ms2pip.core import predict_batch
    from psm_utils import PSM, PSMList
    psms = [PSM(peptidoform=f"{proforma(r.sequence, r.mods, r.mod_sites)}/{int(r.charge)}", spectrum_id=str(i))
            for i, r in enumerate(meta.itertuples(index=False))]
    out = [None] * len(psms)
    for st in range(0, len(psms), 100_000):
        for r in predict_batch(PSMList(psm_list=psms[st:st + 100_000]), model=model, processes=processes):
            b = np.clip(2.0 ** np.asarray(r.predicted_intensity["b"], float) - 0.001, 0, None)
            y = np.clip(2.0 ** np.asarray(r.predicted_intensity["y"], float) - 0.001, 0, None)
            mat = np.full((len(b), len(DEFAULT_FRAG_TYPES)), 0.0, np.float32)
            mat[:, BY1[0]] = b; mat[:, BY1[1]] = y[::-1]  # alphabase rows: y_{n-1-i}
            out[st + r.psm_index] = mat
    return out


def peptdeep_predict(meta, ckpt):
    P = PeptdeepPredictor(ckpt)
    preds = P.predict(meta["sequence"], meta["mods"].fillna(""), meta["mod_sites"].fillna(""), meta["charge"],
                      nce=meta["nce"].to_numpy(float))
    cols = [DEFAULT_FRAG_TYPES.index(t) if t in DEFAULT_FRAG_TYPES else None for t in P.frag_types]
    out = []
    for p in preds:
        mat = np.zeros((p.shape[0], len(DEFAULT_FRAG_TYPES)), np.float32)
        for k, c in enumerate(cols):
            if c is not None:
                mat[:, c] = p[:, k]
        out.append(mat)
    return out, [c for c in cols if c is not None]


def score(meta, X, preds, cols):
    sa = np.zeros(len(meta)); pr = np.zeros(len(meta))
    st = meta["frag_start_idx"].to_numpy(); sp = meta["frag_stop_idx"].to_numpy()
    for i in range(len(meta)):
        t = X[st[i]:sp[i]][None]; p = preds[i][None]
        a, b = spectrum_metrics(p, t, cols)
        sa[i], pr[i] = a[0], b[0]
    return sa, pr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", nargs="+", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--activation", choices=["hcd", "cid"], default=None)
    ap.add_argument("--predictor", action="append", required=True)
    ap.add_argument("--max-spectra", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ms, Xs, off = [], [], 0
    for t in a.table:
        m, X = load(t); m["frag_start_idx"] += off; m["frag_stop_idx"] += off; off += len(X); ms.append(m); Xs.append(X)
    meta = pd.concat(ms, ignore_index=True); X = np.concatenate(Xs)
    if a.activation:
        meta = meta[meta["activation"].eq(a.activation)].reset_index(drop=True)
    if a.max_spectra and len(meta) > a.max_spectra:
        meta = meta.sample(a.max_spectra, random_state=0).reset_index(drop=True)
    meta["benchmark_sequence"] = meta["sequence"].str.replace("I", "L").isin(benchmark_sequences())
    per = pd.DataFrame({"sequence": meta["sequence"], "charge": meta["charge"], "nce": meta["nce"],
                        "benchmark_sequence": meta["benchmark_sequence"]})
    rows = []
    for spec in a.predictor:
        kind, _, arg = spec.partition(":")
        if kind == "ms2pip":
            preds = ms2pip_predict(meta, arg); allcols = None
        else:
            preds, allcols = peptdeep_predict(meta, arg)
        sa, pr = score(meta, X, preds, BY1)
        label = spec if kind == "ms2pip" else ("peptdeep:released" if arg == "released" else f"peptdeep:{Path(arg).stem}")
        per[label + "|sa_by1"] = sa; per[label + "|pcc_by1"] = pr
        row = dict(dataset=a.name, predictor=label, n=len(meta), median_sa_by1=np.nanmedian(sa), median_pcc_by1=np.nanmedian(pr),
                   mean_sa_by1=np.nanmean(sa))
        if allcols is not None:
            sa2, pr2 = score(meta, X, preds, sorted(allcols))
            per[label + "|sa_all"] = sa2
            row.update(median_sa_all=np.nanmedian(sa2), median_pcc_all=np.nanmedian(pr2), n_types=len(allcols))
        rows.append(row); print(row, flush=True)
    res = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out, index=False)
    per.to_parquet(a.out.replace(".csv", ".spectra.parquet"))
    print(res.to_string())


if __name__ == "__main__":
    main()
