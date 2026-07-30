#!/usr/bin/env bash
# ============================================================================
# run_extkin_test.sh — extensive-kinetic (norm-division fix) validation, vacuum,
# DOUBLE-SIDED CAP. Two matched runs on ONE GPU, sequential (clean timing):
#   dcap_extkin   : WP_EXTKIN=1 -> orbital_kinetic_stats.csv (bare per-orbital
#                   KE + norm each step) alongside INQ's reported energies.csv.
#   dcap_baseline : WP_EXTKIN=0 -> per-step wall-time baseline (identical run).
# extkin runs FIRST so any GPU warm-up penalises the instrumented run ->
# conservative (upper-bound) overhead estimate.
#
# Geometry (double-sided): LZ=60 box [-30,30], CAP_L=15 BOTH ends
# (+z [15,30], -z [-30,-15]), launch z=0 -> 5*sigma0=15 Bohr clearance to both
# CAP inner edges (boundary rule). sigma0=3, E=400 eV (k0*sigma0=16.3, ~5%
# transit spread), h=0.4 (cutoff guard PASS, E_cut=839 eV), dt=0.01, eta=-3.5
# (survival exp(-|eta|W/v)~6e-5), NSTEPS=700 (~15 Bohr free flight + full CAP
# transit at v=5.421 + tail).
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
export INQ_SOURCE=$ROOT/inq-study CUDA_VISIBLE_DEVICES=0
export INQ_SHARE_PATH=${INQ_SHARE_PATH:-$ROOT/inq/install/share}
export PSEUDOPOD_SHARE_PATH=${PSEUDOPOD_SHARE_PATH:-$ROOT/inq/install/share/pseudopod}
export PATH="$ROOT/shared/bin:$PATH"
LOG=extkin_test.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

# WF_EVERY=700: VTI frames at t=0 and final step ONLY — /local/data is at 100%
# (2026-07-29, 377 MB free); the comparison is CSV-based, frames not needed.
COMMON="WP_SIGMA=3 WP_K0=5.421 WP_LZ=60 WP_LPERP=30 WP_H=0.4 WP_DT=0.01 \
WP_CAP_L=15 WP_LAUNCH_Z=0 WP_NSTEPS=700 WP_ETA=-3.5 WP_CAP2=1 \
WP_WF_EVERY=700 WP_MOM_EVERY=1 WP_ABS=cap WP_PROP=etrs"

say "=== extkin test: dcap_extkin (observable ON) ==="
env $COMMON WP_EXTKIN=1 WP_EXTKIN_EVERY=1 WP_OUT=dcap_extkin \
    ./run > dcap_extkin.run.log 2>&1 || say "  dcap_extkin FAILED"
say "=== extkin test: dcap_baseline (observable OFF) ==="
env $COMMON WP_EXTKIN=0 WP_OUT=dcap_baseline \
    ./run > dcap_baseline.run.log 2>&1 || say "  dcap_baseline FAILED"
say "=== done ==="
touch EXTKIN_TEST_DONE.txt

# auto-build: regenerate the study notebook from the completed pair
HYP=$ROOT/ResearchProject/systems/vacuum/hypotheses/cap_norm_investigation/extensive_kinetic
PYTHONPATH=$ROOT/inq-stack/python $ROOT/venv/bin/python3 \
    "$HYP/build_extkin_report.py" >> "$LOG" 2>&1 || say "  notebook build FAILED"
