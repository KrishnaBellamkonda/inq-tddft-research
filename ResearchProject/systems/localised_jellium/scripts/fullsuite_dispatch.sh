#!/bin/bash
# Full-suite production dispatcher (D2 re-run). p3wp + p5wp concurrently on the
# two free GPUs (same pre-built fullsuite_wp binary), then p5cl (classical,
# separate build) on GPU0 once a card frees. Runtime needs only the stock-inq
# share paths; p5cl build needs INQ_SOURCE=inq-study.
ROOT=/local/data/public/skcb2/tddft
SHR=$ROOT/inq/install/share
export INQ_SHARE_PATH=$SHR
export PSEUDOPOD_SHARE_PATH=$SHR/pseudopod

WP=$ROOT/ResearchProject/systems/localised_jellium/scripts/fullsuite_wp
CL=$ROOT/ResearchProject/systems/localised_jellium/scripts/fullsuite_classical

cd "$WP" || exit 9
echo "[$(date '+%F %T')] launch p3wp (GPU0, no CAP, 880) + p5wp (GPU1, CAP, 900)"
CUDA_VISIBLE_DEVICES=0 LJ_OUT=p3_wp LJ_CAP=0 ./run > p3wp_run.log 2>&1 &
P3=$!
CUDA_VISIBLE_DEVICES=1 LJ_OUT=p5_wp LJ_CAP=1 ./run > p5wp_run.log 2>&1 &
P5=$!

wait $P3; echo "[$(date '+%F %T')] p3wp finished, exit=$?"
wait $P5; echo "[$(date '+%F %T')] p5wp finished, exit=$?"

cd "$CL" || exit 9
echo "[$(date '+%F %T')] build+run p5cl (classical, GPU0)"
CUDA_VISIBLE_DEVICES=0 INQ_SOURCE=$ROOT/inq-study LJ_OUT=p5_classical LJ_CAP=1 \
  "$ROOT/shared/bin/inq-run" > p5cl_build_run.log 2>&1
echo "[$(date '+%F %T')] p5cl finished, exit=$?"
echo "ALL FULLSUITE PRODUCTION RUNS COMPLETE"
