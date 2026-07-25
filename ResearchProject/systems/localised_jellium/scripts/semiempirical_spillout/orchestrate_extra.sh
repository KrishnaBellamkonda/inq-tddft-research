#!/bin/bash
# Semi-empirical spill-out: finish the p2 matrix + add p3 (PBC) runs for the open-z test.
# 4 pending p2 (lz240 N164 N328 es60) + 3 p3 (p3_lz90 p3_lz160 p3_lz240). Both GPUs.
set -u
ROOT=/local/data/public/skcb2/tddft
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$INQ_SHARE_PATH/pseudopod
export INQ_SOURCE=$ROOT/inq-study
BIN=$ROOT/ResearchProject/systems/localised_jellium/scripts/semiempirical_spillout/gs/run
BASE=$ROOT/ResearchProject/systems/localised_jellium/scripts/semiempirical_spillout/runs
mkdir -p "$BASE"

# tag|LZ|N|EDGE_W|EXTRA_STATES|PERIODICITY  (half=12.5, spacing=0.5, T=0.00862)
GPU0_RUNS=( "N328|160|328|0|20|2" "es60|160|82|0|60|2" "p3_lz90|90|82|0|20|3" "p3_lz160|160|82|0|20|3" )
GPU1_RUNS=( "lz240|240|82|0|20|2" "N164|160|164|0|20|2" "p3_lz240|240|82|0|20|3" )

run_queue () {
  local gpu=$1; shift
  for spec in "$@"; do
    IFS='|' read -r tag LZ N EW ES PER <<< "$spec"
    local d="$BASE/$tag"; mkdir -p "$d"
    echo "[GPU$gpu] START $tag  (Lz=$LZ N=$N w=$EW es=$ES per=$PER)  $(date +%T)"
    ( cd "$d" && CUDA_VISIBLE_DEVICES=$gpu \
        LJ_LX=50 LJ_LY=50 LJ_LZ=$LZ LJ_HALF=12.5 LJ_N=$N LJ_EDGE_W=$EW \
        LJ_PERIODICITY=$PER LJ_SPACING=0.5 LJ_EXTRA_STATES=$ES LJ_TEMP_EV=0.00862 \
        LJ_GS_DIR="$d/checkpoint" LJ_TAG=$tag "$BIN" > "$d/run.log" 2>&1 )
    if grep -q "run_completed = true" "$d/results/run_summary.txt" 2>/dev/null; then
      echo "[GPU$gpu] DONE  $tag  $(date +%T)"
    else
      echo "[GPU$gpu] FAIL  $tag  (see $d/run.log)  $(date +%T)"
    fi
  done
  echo "[GPU$gpu] QUEUE COMPLETE $(date +%T)"
}

run_queue 0 "${GPU0_RUNS[@]}" &
P0=$!
run_queue 1 "${GPU1_RUNS[@]}" &
P1=$!
wait $P0 $P1
echo "ALL EXTRA QUEUES COMPLETE $(date +%T)"
