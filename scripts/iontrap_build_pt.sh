#!/bin/bash
# Build the ProteomeTools fragment-intensity tables from the labelled runs: one table per pool, split
# train/val by the pool lists (benchmark pools are refused by iontrap_dataset.py).
source /home/drew/miniforge3/etc/profile.d/conda.sh; conda activate peaksparity
cd /home/drew/code/proteoscan
W=/drive_a/drew/proteoscan-work/models/iontrap_model
for split in train val; do
  for pool in $(cat $W/${split}_pools.txt); do
    R=$(ls -d $W/runs/*-${pool}_01_01-2xIT_2xHCD-* 2>/dev/null | head -1)
    [ -n "$R" ] && [ -s $R/labels_ev.parquet ] || continue
    out=$W/data/pt_${split}_${pool}
    [ -s $out.npy ] || echo "python scripts/iontrap_dataset.py --labels $R/labels_ev.parquet --split $split --out $out > $W/logs/ds_$pool.log 2>&1"
  done
done | xargs -P ${1:-6} -I{} bash -c '{}'
