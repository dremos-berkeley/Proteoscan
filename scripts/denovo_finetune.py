"""Fine-tune Casanovo on annotated MGFs within a wall-clock budget, keeping snapshots.

Writes <out>/config.yaml (Casanovo defaults + the overrides below), runs
``casanovo train --model <ckpt> --config config.yaml -p <valid.mgf> <train.mgf...>`` and every
--snapshot-min minutes copies the current best checkpoint (lowest validation CE loss) to
<out>/snapshots/best.<minutes>min.ckpt. Training is stopped (SIGINT, then SIGTERM) after
--hours; <out>/ft.best.ckpt is the selected model.

Example:
  python scripts/denovo_finetune.py --out /drive_a/drew/peaks-parity-work/denovo_finetune/runs/ft1 \
      --model ~/.cache/casanovo/casanovo_orbitrap_v5-2-0.ckpt \
      --train '/drive_a/.../mgf/train.*.mgf' --valid /drive_a/.../mgf/valid.*.mgf --gpu 2 --hours 4
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaksparity import denovo  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", help="checkpoint to start from; omit to train from scratch")
    ap.add_argument("--train", required=True, help="glob of training MGFs")
    ap.add_argument("--valid", required=True, help="glob of validation MGFs")
    ap.add_argument("--gpu", default="2")
    ap.add_argument("--hours", type=float, default=4.0)
    ap.add_argument("--snapshot-min", type=float, default=30.0)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--warmup", type=int, default=1000)
    ap.add_argument("--period", type=int, default=60000, help="cosine half period (steps)")
    ap.add_argument("--val-every", type=int, default=2000, help="validation interval (steps)")
    ap.add_argument("--max-epochs", type=int, default=1000)
    ap.add_argument("--precision", default="32-true")
    ap.add_argument("--dropout", type=float, default=None, help="only used when training from scratch")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--lance-dir", help="keep Casanovo's parsed Lance datasets here; later runs can pass "
                    "<lance-dir>/train.lance and valid_0.lance as --train/--valid to skip MGF parsing")
    a = ap.parse_args()

    out = Path(a.out)
    (out / "snapshots").mkdir(parents=True, exist_ok=True)
    overrides = dict(
        learning_rate=a.lr, warmup_iters=a.warmup, cosine_schedule_period_iters=a.period,
        train_batch_size=a.batch, val_check_interval=a.val_every, max_epochs=a.max_epochs,
        precision=a.precision, calculate_precision=True, log_metrics=True, tb_summarywriter=False,
        # low-resolution MS2 preprocessing (these are also Casanovo's defaults)
        min_mz=50.0, max_peaks=150, max_charge=4, precursor_mass_tol=15, predict_batch_size=512,
        devices=1,
    )
    if a.lance_dir:
        Path(a.lance_dir).mkdir(parents=True, exist_ok=True)
        overrides["lance_dir"] = str(a.lance_dir)
    if a.dropout is not None:
        overrides["dropout"] = a.dropout
    cfg = denovo.write_config(out / "config.yaml", **overrides)
    train = sorted(glob.glob(a.train))
    valid = sorted(glob.glob(a.valid))
    assert train and valid and not (set(train) & set(valid))
    cmd = ["casanovo", "train", "-c", str(cfg), "-d", str(out), "-o", "ft", "-f", "-v", "info"]
    if a.model:
        cmd += ["-m", str(Path(a.model).expanduser())]
    for v in valid:
        cmd += ["-p", v]
    cmd += train
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(a.gpu), OMP_NUM_THREADS=str(a.threads),
               MKL_NUM_THREADS=str(a.threads))
    (out / "command.txt").write_text(" ".join(cmd) + "\n")
    t0 = time.time()
    with open(out / "train.stdout", "w") as log:
        p = subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
        next_snap = a.snapshot_min * 60
        best = out / "ft.best.ckpt"
        while p.poll() is None:
            time.sleep(20)
            el = time.time() - t0
            if el >= next_snap:
                if best.exists():
                    shutil.copy2(best, out / "snapshots" / f"best.{int(el // 60)}min.ckpt")
                next_snap += a.snapshot_min * 60
            if el >= a.hours * 3600:
                p.send_signal(signal.SIGINT)
                try:
                    p.wait(timeout=600)
                except subprocess.TimeoutExpired:
                    p.terminate()
                    p.wait()
                break
    el = time.time() - t0
    if best.exists():
        shutil.copy2(best, out / "snapshots" / f"best.{int(el // 60)}min.final.ckpt")
    (out / "elapsed_seconds.txt").write_text(f"{el:.0f}\n")
    print(f"done in {el / 3600:.2f} h, exit {p.returncode}")


if __name__ == "__main__":
    main()
