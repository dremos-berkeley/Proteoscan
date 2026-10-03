#!/bin/bash
# Run the standard pipeline (as for the benchmark pools) on the held-out validation pools' ion-trap
# scan types, to (1) validate the targeted labelling against pipeline PSMs and (2) select the
# fragment-model checkpoint / feature variant without touching the benchmark pools.
#   scripts/iontrap_val_pipeline.sh <pool list> [gpu]
source /home/drew/miniforge3/etc/profile.d/conda.sh; conda activate peaksparity
cd /home/drew/code/peaks-parity
W=/drive_a/drew/peaks-parity-work/iontrap_model; BENCH=/drive_a/drew/public-data/PXD004732/pools.txt
LIST=$1; GPU=${2:-3}
for pool in $(cat $LIST); do
  if grep -qx "$pool" $BENCH; then echo "REFUSE benchmark pool $pool"; continue; fi
  while true; do
    R=$(ls -d $W/runs/*-${pool}_01_01-2xIT_2xHCD-* 2>/dev/null | head -1)
    [ -n "$R" ] && [ -f "$R/label.log" ] && break  # label.log appears once the split is complete
    sleep 60
  done
  b=$(basename $R)
  for t in hcd_it cid_it; do
    O=$W/val_pipeline/$b/$t; mkdir -p $O
    [ -f $O/run_stats.json ] && continue
    CUDA_VISIBLE_DEVICES=$GPU python -m proteoscan.pipeline --config configs/peaks_facility_human.yml --input $R/split/$b.$t.mzML --out $O --sample-name "$b.$t" --processes 12 > $O.log 2>&1
    echo "VAL RUN EXIT $? $b $t $(date)" >> $W/logs/val_pipeline.log
  done
done
