#!/usr/bin/env bash
# Phase-4 production launcher: energy-matched S(v=2.0 / 54.42 eV) point.
# WP (unchanged phase-3 mechanism, 90-box, two-sided CAP, equidistant launch)
# + redesigned classical (Ehrenfest -> park+remove at |z|>=35 -> projectile-free
# tail). Reuses the smoke-built binaries (inq-run just runs if up to date).
#
# Usage:  bash run_production.sh [WP_TAU] [CL_TAU] [DT]
#   WP_TAU default 100  (WP absorbed by tau=100, as phase-3)
#   CL_TAU default 160  (longer: the 54 eV ion must travel out to |z|=35; sized
#                        from the smoke trajectory — see handover)
#   DT     default 0.04 (smoke-confirmed stable)
set -u
ROOT=/local/data/public/skcb2/tddft
P4=$ROOT/ResearchProject/systems/localised_jellium/scripts/qsp_phase4
WP_TAU=${1:-100.0}
CL_TAU=${2:-160.0}
DT=${3:-0.04}
PY=$ROOT/venv/bin/python3
WP_NSTEPS=$($PY -c "print(int(round($WP_TAU/$DT)))")
CL_NSTEPS=$($PY -c "print(int(round($CL_TAU/$DT)))")
WP_WE=$($PY -c "print(max(1, $WP_NSTEPS//300))")
CL_WE=$($PY -c "print(max(1, $CL_NSTEPS//300))")
export INQ_SOURCE=$ROOT/inq-study
echo "[prod-p4] WP: tau=$WP_TAU n=$WP_NSTEPS we=$WP_WE | CL: tau=$CL_TAU n=$CL_NSTEPS we=$CL_WE chunk=25 park|z|>=35 | dt=$DT"

( cd "$P4/wp"        && CUDA_VISIBLE_DEVICES=0 LJ_OUT=p4_wp        LJ_CAP=1 LJ_DT=$DT \
    LJ_N_STEPS=$WP_NSTEPS LJ_LAUNCH_Z=-23.75 LJ_WRITE_EVERY=$WP_WE LJ_WF_EVERY=$WP_WE \
    $ROOT/shared/bin/inq-run > prod_wp.log 2>&1 ) & WP=$!
( cd "$P4/classical" && CUDA_VISIBLE_DEVICES=1 LJ_OUT=p4_classical LJ_CAP=1 LJ_DT=$DT \
    LJ_N_STEPS=$CL_NSTEPS LJ_LAUNCH_Z=-23.75 LJ_WRITE_EVERY=$CL_WE LJ_CHUNK=25 \
    $ROOT/shared/bin/inq-run > prod_cl.log 2>&1 ) & CL=$!

wait $WP; rc_wp=$?; echo "[prod-p4] WP finished exit=$rc_wp (log: $P4/wp/prod_wp.log)"
wait $CL; rc_cl=$?; echo "[prod-p4] classical finished exit=$rc_cl (log: $P4/classical/prod_cl.log)"
echo "[prod-p4] DONE  wp_exit=$rc_wp cl_exit=$rc_cl"
