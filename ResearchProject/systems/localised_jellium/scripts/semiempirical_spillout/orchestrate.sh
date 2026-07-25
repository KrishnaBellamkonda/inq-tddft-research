#!/bin/bash
# Semi-empirical spill-out GS matrix (campaign localised_jellium_parameter_study_2).
# 10 p2 GS SCF runs across BOTH GPUs. Reference: Lz=160, N=82, w=0, es=20.
#   A box/Lz : lz90 lz120 lz160 lz240      (numerical-vs-physical: floor charge vs vacuum)
#   B soft-w : w1 w2 w4                     (does edge softness move the plateau?)
#   C confine: N164 N328                    (does deeper well reduce spill?)
#   D solver : es60                         (is the floor thermal/state-count?)
set -u
ROOT=/local/data/public/skcb2/tddft
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$INQ_SHARE_PATH/pseudopod
export INQ_SOURCE=$ROOT/inq-study
BIN=$ROOT/ResearchProject/systems/localised_jellium/scripts/semiempirical_spillout/gs/run
BASE=$ROOT/ResearchProject/systems/localised_jellium/scripts/semiempirical_spillout/runs
mkdir -p "$BASE"

# tag|LZ|N|EDGE_W|EXTRA_STATES  (all p2, half=12.5, spacing=0.5, T=0.00862)
GPU0_RUNS=( "lz90|90|82|0|20" "w1|160|82|1|20" "w4|160|82|4|20" "N164|160|164|0|20" "lz240|240|82|0|20" )
GPU1_RUNS=( "lz120|120|82|0|20" "lz160|160|82|0|20" "w2|160|82|2|20" "es60|160|82|0|60" "N328|160|328|0|20" )

run_queue () {
  local gpu=$1; shift
  local specs=("$@")
  for spec in "${specs[@]}"; do
    IFS='|' read -r tag LZ N EW ES <<< "$spec"
    local d="$BASE/$tag"; mkdir -p "$d"
    echo "[GPU$gpu] START $tag  (Lz=$LZ N=$N w=$EW es=$ES)  $(date +%T)"
    ( cd "$d" && CUDA_VISIBLE_DEVICES=$gpu \
        LJ_LX=50 LJ_LY=50 LJ_LZ=$LZ LJ_HALF=12.5 LJ_N=$N LJ_EDGE_W=$EW \
        LJ_PERIODICITY=2 LJ_SPACING=0.5 LJ_EXTRA_STATES=$ES LJ_TEMP_EV=0.00862 \
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
echo "ALL QUEUES COMPLETE $(date +%T)"
