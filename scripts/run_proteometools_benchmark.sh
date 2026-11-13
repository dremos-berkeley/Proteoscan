#!/bin/bash
# ProteomeTools (PXD004732) ground-truth benchmark, one lane: split each pool's 2xIT_2xHCD run by scan type and
# run the default pipeline on every type. Start several lanes with different GPUs and pool lists, then run
# scripts/aggregate_proteometools.py --work $WORK.
#
#   WORK=<benchmark dir> scripts/run_proteometools_benchmark.sh GPU POOL [POOL ...]      (POOL = number, e.g. 6)
#
# Environment: DATA (default /drive_a/drew/public-data/PXD004732, with mzml/), IT_CONFIG and FT_CONFIG (search
# configs for ion-trap and Orbitrap MS2), PROCESSES (default 8), SPLIT_FROM (optional directory holding
# <run>/split from an earlier benchmark, reused through symlinks).
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd)
DATA=${DATA:-/drive_a/drew/public-data/PXD004732}
IT_CONFIG=${IT_CONFIG:-$REPO/configs/peaks_facility_human.yml}
FT_CONFIG=${FT_CONFIG:-$REPO/configs/peaks_hela_hcd_ft.yml}
PROCESSES=${PROCESSES:-8}
: "${WORK:?set WORK to the benchmark output directory}"
GPU=$1; shift
cd "$REPO"; export PYTHONPATH=$REPO
mkdir -p "$WORK"
for n in "$@"; do
  pool=TUM_first_pool_$n
  mz=$(ls "$DATA"/mzml/*-${pool}_01_01-2xIT_2xHCD-*.mzML | head -1); b=$(basename "$mz" .mzML); W=$WORK/$b
  mkdir -p "$W"
  if [ ! -d "$W/split" ]; then
    if [ -n "${SPLIT_FROM:-}" ] && [ -d "$SPLIT_FROM/$b/split" ]; then ln -s "$SPLIT_FROM/$b/split" "$W/split"
    else python -m proteoscan.split_mzml "$mz" "$W/split" >> "$WORK/split.log" 2>&1; fi
  fi
  feat=""
  for t in hcd_it cid_it hcd_ft; do
    m=$W/split/$b.$t.mzML; o=$W/$t
    case $t in *_ft) c=$FT_CONFIG;; *) c=$IT_CONFIG;; esac
    if [ ! -f "$o/run_stats.json" ]; then
      # the three splits share the run's MS1 scans: detect features once and reuse them
      CUDA_VISIBLE_DEVICES=$GPU python -m proteoscan.pipeline --config "$c" --input "$m" --out "$o" \
        --sample-name "$b.$t" --processes "$PROCESSES" --no-summary $feat > "$o.log" 2>&1
      echo "RUN EXIT $? $pool $t $(date)" >> "$WORK/progress.log"
    fi
    [ -z "$feat" ] && f0=$(ls "$o"/*/ms1_features.parquet 2>/dev/null | head -1) && [ -n "$f0" ] && feat="--features $f0"
  done
done
echo "LANE DONE $GPU $(date)" >> "$WORK/progress.log"
