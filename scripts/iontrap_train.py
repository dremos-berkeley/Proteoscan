"""Fine-tune AlphaPeptDeep's released MS2 model on ion-trap spectra (docs/IONTRAP_MODEL.md).

    python scripts/iontrap_train.py --train data/pt_train data/hela_* --val data/pt_val \
        --activation hcd --out models/hcd_it_v1 [--epochs 30 --lr 1e-4 --patience 4]

Transfer learning: every layer of the released generic model is loaded; the output head is widened
from (b, y) x (1+, 2+) to also predict b/y - H2O and - NH3 (1+), the new rows initialised from the
parent ion's row. Loss: masked spectral angle (Prosit-style) over all observable ions; unobservable
ions (-1 in the tables) are masked. Early stopping on the median spectral angle of the validation
spectra (held-out ProteomeTools pools, never the benchmark pools). Writes <out>.pth + <out>.json
(fragpred.PeptdeepPredictor reads both) and <out>.history.csv.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan.fragpred import DEFAULT_FRAG_TYPES  # noqa: E402

BENCH = set(Path("/drive_a/drew/public-data/PXD004732/pools.txt").read_text().split())
TRUTH = Path("/drive_a/drew/public-data/PXD004732/truth")


def benchmark_sequences() -> set[str]:
    t = pd.read_parquet(TRUTH / "pool_peptides.parquet")
    s = set(t.loc[t["pool"].isin(BENCH), "sequence"].str.replace("I", "L"))
    s |= set(pd.read_parquet(TRUTH / "qc_peptides.parquet")["sequence"].str.replace("I", "L"))
    return s


def load(prefixes: list[str]):
    metas, mats, off = [], [], 0
    for p in prefixes:
        m = pd.read_parquet(p + ".parquet"); x = np.load(p + ".npy")
        m["frag_start_idx"] += off; m["frag_stop_idx"] += off
        off += len(x)
        metas.append(m); mats.append(x)
    meta = pd.concat(metas, ignore_index=True)
    if "pool" in meta and meta["pool"].isin(BENCH).any():
        raise SystemExit("benchmark pool found in a training/validation table")
    return meta, np.concatenate(mats)


def select(meta: pd.DataFrame, activation: str, exclude: set[str] | None):
    keep = meta["activation"].eq(activation)
    if exclude:
        keep &= ~meta["sequence"].str.replace("I", "L").isin(exclude)
    return meta[keep].reset_index(drop=True)


def init_model(frag_types, device):
    from peptdeep.model.ms2 import pDeepModel
    from peptdeep.pretrained_models import ModelManager
    mm = ModelManager(device=device); mm.load_installed_models()
    old = mm.ms2_model.model.state_dict()
    old_types = [t for t in mm.ms2_model.model.supported_charged_frag_types if "modloss" not in t]
    model = pDeepModel(charged_frag_types=frag_types, device=device)
    new = model.model.state_dict()
    for k, v in old.items():
        if k in new and new[k].shape == v.shape:
            new[k] = v.clone()
    W, b = old["output_nn.nn.2.weight"], old["output_nn.nn.2.bias"]
    for i, t in enumerate(model.model.supported_charged_frag_types):
        parent = t if t in old_types else t.split("_")[0] + "_z1"
        j = old_types.index(parent)
        new["output_nn.nn.2.weight"][i] = W[j]
        new["output_nn.nn.2.bias"][i] = b[j] - (0.0 if t in old_types else 0.3)
    model.model.load_state_dict(new)
    return model


def batches(meta: pd.DataFrame, batch_size: int, rng: np.random.Generator | None):
    groups = [g.index.to_numpy() for _, g in meta.groupby("nAA")]
    out = []
    for idx in groups:
        if rng is not None:
            idx = rng.permutation(idx)
        out += [idx[i:i + batch_size] for i in range(0, len(idx), batch_size)]
    if rng is not None:
        out = [out[i] for i in rng.permutation(len(out))]
    return out


def targets(meta: pd.DataFrame, X: np.ndarray, idx: np.ndarray):
    st = meta["frag_start_idx"].to_numpy()[idx]
    L = int(meta["nAA"].to_numpy()[idx[0]]) - 1
    rows = st[:, None] + np.arange(L)[None, :]
    return X[rows]  # (B, L, K)


def sa_loss(pred, tgt):
    mask = (tgt >= 0).float()
    p = pred * mask
    t = tgt.clamp(min=0) * mask
    p = p.flatten(1); t = t.flatten(1)
    cos = (p * t).sum(1) / (p.norm(dim=1) * t.norm(dim=1) + 1e-8)
    cos = cos.clamp(-1 + 1e-7, 1 - 1e-7)
    return (2 * torch.arccos(cos) / np.pi).mean()


def spectrum_metrics(pred: np.ndarray, tgt: np.ndarray, cols=None):
    """Per-spectrum spectral angle (1 - 2 angle/pi, on max-normalised linear intensities) and Pearson,
    over observable ions; cols restricts to fragment-type columns."""
    if cols is not None:
        pred = pred[..., cols]; tgt = tgt[..., cols]
    B = pred.shape[0]
    p = pred.reshape(B, -1).astype(float); t = tgt.reshape(B, -1).astype(float)
    m = t >= 0
    p = np.where(m, np.clip(p, 0, None), 0); t = np.where(m, t, 0)
    cos = (p * t).sum(1) / (np.linalg.norm(p, axis=1) * np.linalg.norm(t, axis=1) + 1e-12)
    sa = 1 - 2 * np.arccos(np.clip(cos, -1, 1)) / np.pi
    pr = np.full(B, np.nan)
    for i in range(B):
        mi = m[i]
        if mi.sum() >= 3 and p[i][mi].std() > 0 and t[i][mi].std() > 0:
            pr[i] = np.corrcoef(p[i][mi], t[i][mi])[0, 1]
    return sa, pr


@torch.no_grad()
def predict_tables(model, meta, X, batch_size=512):
    model.model.eval()
    out = np.zeros_like(X)
    for idx in batches(meta, batch_size, None):
        feats = model._get_features_from_batch_df(meta.iloc[idx])
        pr = model.model(*feats).float().cpu().numpy()
        st = meta["frag_start_idx"].to_numpy()[idx]
        L = pr.shape[1]
        out[(st[:, None] + np.arange(L)[None, :]).ravel()] = pr.reshape(-1, pr.shape[2])
    return out


def evaluate(model, meta, X, cols=None):
    P = predict_tables(model, meta, X)
    sas, prs = [], []
    for idx in batches(meta, 4096, None):
        sa, pr = spectrum_metrics(targets(meta, P, idx), targets(meta, X, idx), cols)
        sas.append(sa); prs.append(pr)
    return np.concatenate(sas), np.concatenate(prs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", nargs="+", required=True)
    ap.add_argument("--val", nargs="+", required=True)
    ap.add_argument("--activation", default="hcd", choices=["hcd", "cid"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--patience", type=int, default=4)
    ap.add_argument("--warmup-epochs", type=int, default=1)
    ap.add_argument("--keep-benchmark-sequences", action="store_true",
                    help="do not drop training spectra of peptides that are in the benchmark pools")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--hashsplit-val", type=int, default=10,
                    help="non-ProteomeTools training rows (HeLa) whose peptide hashes to 0 mod N are held out "
                         "as a peptide-disjoint second validation set (0 = off)")
    ap.add_argument("--other-per-epoch", type=int, default=0,
                    help="subsample non-ProteomeTools training rows to this many per epoch (0 = all)")
    ap.add_argument("--head-epochs", type=int, default=1, help="epochs training only the output heads first")
    a = ap.parse_args()
    torch.manual_seed(a.seed); np.random.seed(a.seed)
    rng = np.random.default_rng(a.seed)
    excl = None if a.keep_benchmark_sequences else benchmark_sequences()
    tm, TX = load(a.train); vm, VX = load(a.val)
    n_before = int(tm["activation"].eq(a.activation).sum())
    tm = select(tm, a.activation, excl); vm = select(vm, a.activation, excl)
    hm = tm.iloc[:0]
    if a.hashsplit_val:
        import zlib
        h = tm["sequence"].str.replace("I", "L").map(lambda q: zlib.crc32(q.encode()) % a.hashsplit_val == 0)
        other = tm["source"].ne("proteometools")
        hm = tm[other & h].reset_index(drop=True)
        tm = tm[~(other & h)].reset_index(drop=True)
        # a peptide held out from the complex-sample data must not be trained on through another source either
        held = set(hm["sequence"].str.replace("I", "L"))
        tm = tm[~tm["sequence"].str.replace("I", "L").isin(held)].reset_index(drop=True)
    print(f"train {len(tm)} spectra ({n_before - len(tm) - len(hm)} dropped: benchmark sequences), val {len(vm)}, "
          f"val (held-out peptides of non-ProteomeTools sources) {len(hm)}", flush=True)
    print(tm.groupby(["source", "nce", "scan_mode"]).size().to_string(), flush=True)
    model = init_model(DEFAULT_FRAG_TYPES, "gpu")
    by_cols = [DEFAULT_FRAG_TYPES.index(c) for c in ("b_z1", "y_z1")]
    sa0, pr0 = evaluate(model, vm, VX)
    sab0, _ = evaluate(model, vm, VX, by_cols)
    print(f"epoch 0 (released weights + copied loss head): val median SA {np.median(sa0):.4f} PCC {np.nanmedian(pr0):.4f} "
          f"SA b/y 1+ {np.median(sab0):.4f}", flush=True)
    is_other = tm["source"].ne("proteometools").to_numpy()

    def epoch_rows():
        if a.other_per_epoch and is_other.sum() > a.other_per_epoch:
            keep = ~is_other
            keep[rng.choice(np.nonzero(is_other)[0], a.other_per_epoch, replace=False)] = True
            return tm[keep]
        return tm
    opt = torch.optim.AdamW(model.model.parameters(), lr=a.lr, weight_decay=1e-4)
    nb = len(batches(epoch_rows(), a.batch_size, rng))
    total = a.epochs * nb; warm = a.warmup_epochs * nb
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / max(warm, 1)) * 0.5 * (1 + np.cos(np.pi * min(s, total) / total)))
    hist, best, bad = [], -1.0, 0
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    for ep in range(1, a.epochs + 1):
        model.model.train(); t0 = time.time(); losses = []
        head_only = ep <= a.head_epochs
        for n_, p_ in model.model.named_parameters():
            p_.requires_grad = (not head_only) or n_.startswith("output_nn")
        for idx in batches(epoch_rows(), a.batch_size, rng):
            feats = model._get_features_from_batch_df(tm.iloc[idx])
            tgt = torch.as_tensor(targets(tm, TX, idx), device=feats[0].device)
            pred = model.model(*feats)
            loss = sa_loss(pred, tgt)
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(model.model.parameters(), 1.0)
            opt.step(); sched.step(); losses.append(loss.item())
        sa, pr = evaluate(model, vm, VX)
        sab, _ = evaluate(model, vm, VX, by_cols)
        h = dict(epoch=ep, head_only=head_only, train_loss=float(np.mean(losses)), val_median_sa=float(np.median(sa)),
                 val_median_pcc=float(np.nanmedian(pr)), val_median_sa_by1=float(np.median(sab)),
                 lr=sched.get_last_lr()[0], sec=time.time() - t0)
        crit = h["val_median_sa"]
        if len(hm):
            sa2, _ = evaluate(model, hm, TX); sab2, _ = evaluate(model, hm, TX, by_cols)
            h.update(val2_median_sa=float(np.median(sa2)), val2_median_sa_by1=float(np.median(sab2)))
            crit = 0.5 * (crit + h["val2_median_sa"])
        h["criterion"] = crit
        hist.append(h); print(h, flush=True)
        if crit > best + 1e-4:
            best, bad = crit, 0
            torch.save(model.model.state_dict(), str(out) + ".pth")
        else:
            bad += 1
            if bad >= a.patience:
                print(f"early stop at epoch {ep}; best val median SA {best:.4f}")
                break
    meta = dict(charged_frag_types=DEFAULT_FRAG_TYPES, nce=float(tm["nce"].median()), instrument="Lumos",
                activation=a.activation, base="peptdeep pretrained_models_v3 generic ms2", loss="masked spectral angle",
                train=a.train, val=a.val, n_train=len(tm), n_val=len(vm), best_val_median_sa=best,
                epochs_run=len(hist), lr=a.lr, batch_size=a.batch_size, seed=a.seed, n_val2=len(hm),
                hashsplit_val=a.hashsplit_val, other_per_epoch=a.other_per_epoch, head_epochs=a.head_epochs,
                criterion="median SA (all fragment types) on validation pools" + (" averaged with held-out-peptide complex-sample SA" if len(hm) else ""),
                excluded_benchmark_sequences=excl is not None, released_val_median_sa=float(np.median(sa0)))
    Path(str(out) + ".json").write_text(json.dumps(meta, indent=1))
    pd.DataFrame(hist).to_csv(str(out) + ".history.csv", index=False)


if __name__ == "__main__":
    main()
