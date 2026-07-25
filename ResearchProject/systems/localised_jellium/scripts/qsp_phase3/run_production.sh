#!/usr/bin/env bash
# Phase-3 production launcher: big-box (50x50x90) WP + classical pair, two-sided CAP,
# equidistant launch (-23.75), tau = 100 a.u. Reuses the smoke-built binaries (inq-run
# detects the up-to-date build and just runs). Usage:  bash run_production.sh [DT]
#   DT defaults to 0.04 (set 0.02 if the smoke shows 0.04 is unstable).
set -u
ROOT=/local/data/public/skcb2/tddft
P3=$ROOT/ResearchProject/systems/localised_jellium/scripts/qsp_phase3
DT=${1:-0.04}
TAU=100.0
PY=$ROOT/venv/bin/python3
NSTEPS=$($PY -c "print(int(round($TAU/$DT)))")
WE=$($PY -c "print(max(1, $NSTEPS//300))")     # ~300 density/wavefunction frames
export INQ_SOURCE=$ROOT/inq-study
echo "[production] dt=$DT  n_steps=$NSTEPS  write_every=$WE  (tau=$TAU a.u.), two-sided CAP"

( cd "$P3/wp"        && CUDA_VISIBLE_DEVICES=0 LJ_OUT=p3_wp        LJ_CAP=1 LJ_DT=$DT \
    LJ_N_STEPS=$NSTEPS LJ_LAUNCH_Z=-23.75 LJ_WRITE_EVERY=$WE LJ_WF_EVERY=$WE \
    $ROOT/shared/bin/inq-run > prod_wp.log 2>&1 ) & WP=$!
( cd "$P3/classical" && CUDA_VISIBLE_DEVICES=1 LJ_OUT=p3_classical LJ_CAP=1 LJ_DT=$DT \
    LJ_N_STEPS=$NSTEPS LJ_LAUNCH_Z=-23.75 LJ_WRITE_EVERY=$WE \
    $ROOT/shared/bin/inq-run > prod_cl.log 2>&1 ) & CL=$!

wait $WP; rc_wp=$?; echo "[production] WP finished exit=$rc_wp (log: $P3/wp/prod_wp.log)"
wait $CL; rc_cl=$?; echo "[production] classical finished exit=$rc_cl (log: $P3/classical/prod_cl.log)"
echo "[production] DONE  wp_exit=$rc_wp cl_exit=$rc_cl"
