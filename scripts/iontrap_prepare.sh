#!/bin/bash
# Convert, split and label the ProteomeTools training/validation pools (docs/IONTRAP_MODEL.md).
# Raw files are SHA1-checked against PRIDE before use (a parallel-download race once corrupted some).
#   scripts/iontrap_prepare.sh <pool list file> <lane> <nlanes>
# Waits for each raw file to finish downloading (no .part left).
source /home/drew/miniforge3/etc/profile.d/conda.sh; conda activate peaksparity
cd /home/drew/code/proteoscan
RAW=/drive_a/drew/public-data/PXD004732/raw_train; W=/drive_a/drew/proteoscan-work/models/iontrap_model
BENCH=/drive_a/drew/public-data/PXD004732/pools.txt
LIST=$1; LANE=${2:-0}; NLANES=${3:-1}; i=0
for pool in $(cat $LIST); do
  i=$((i+1)); [ $((i % NLANES)) -eq $LANE ] || continue
  if grep -qx "$pool" $BENCH; then echo "REFUSE benchmark pool $pool"; continue; fi
  while true; do
    raw=$(ls $RAW/*-${pool}_01_01-2xIT_2xHCD-*.raw 2>/dev/null | head -1)
    [ -n "$raw" ] && break
    sleep 60
  done
  b=$(basename "$raw" .raw); R=$W/runs/$b; mkdir -p $R
  [ -s "$R/$b.mzML" ] || ThermoRawFileParser -i "$raw" -o $R -f 2 > $R/convert.log 2>&1
  [ -d "$R/split" ] || python -m proteoscan.split_mzml "$R/$b.mzML" "$R/split" > $R/split.log 2>&1
  # IT-only labels (comparison) and block labels (used for training; docs/IONTRAP_MODEL.md)
  [ -s "$R/labels.parquet" ] || python scripts/iontrap_label.py --pool $pool --split-dir $R/split --out $R/labels.parquet > $R/label.log 2>&1
  [ -s "$R/labels_ev.parquet" ] || python scripts/iontrap_label.py --pool $pool --mzml $R/$b.mzML --out $R/labels_ev.parquet > $R/labels_ev.log 2>&1
  echo "DONE $pool $b $(date)" >> $W/logs/prepare.log
done
echo "LANE $LANE DONE $(date)" >> $W/logs/prepare.log
