#!/usr/bin/env bash
# ============================================================================
# run_extkin_plateau.sh — extkin_plateau_E100: N=92 / r_s=4.0 / 20.13-Bohr slab,
# 35x35x120 box, WP sigma=1.5 / E=100 eV, TWO-SIDED CAP (15 Bohr, eta=-1.0),
# dt=0.04 x 1500 steps (t=60 a.u.). CAP run ONLY (user decision 2026-07-29:
# no no-CAP twin, no dt gate). First jellium run with the norm-division fix
# IN-RUN (OrbitalKineticStats, all states, every step).
# Plan: docs/plans/norm-corrected-stopping-power.md "Run design (2026-07-29)".
#
# Sequence: GS (skip if present) -> CAP run -> notebook auto-build.
# Fully autonomous; emails each milestone via notify.py. Resumable
# (WP_RESUME=1 + larger WP_N_STEPS extends the finished run).
# Env override: EXTKIN_GPU (default 0).
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
SYS=$ROOT/ResearchProject/systems/localised_jellium
PY=$ROOT/venv/bin/python3
export INQ_SOURCE=$ROOT/inq-study CUDA_VISIBLE_DEVICES=${EXTKIN_GPU:-0}
export INQ_SHARE_PATH=${INQ_SHARE_PATH:-$ROOT/inq/install/share}
export PSEUDOPOD_SHARE_PATH=${PSEUDOPOD_SHARE_PATH:-$ROOT/inq/install/share/pseudopod}
export PATH="$ROOT/shared/bin:$PATH" PYTHONPATH=$ROOT/inq-stack/python
GSDIR=$SYS/shared_gs/slab_n92_L35x35x120_w0p5_h0p5
LOG=$(pwd)/extkin_plateau.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$SYS/scripts/wp_cap_energy_plateau/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }
fail(){ say "FAILED: $1"; notify "extkin_plateau_E100 FAILED: $1" "See $LOG"; exit 1; }

say "=== extkin_plateau_E100 (GPU ${EXTKIN_GPU:-0}) ==="
notify "extkin_plateau_E100 STARTED" "GS -> CAP (1500 steps, dt=0.04); log $LOG"

# --- 1. Ground state (new N=92 box) -----------------------------------------
if [ -f "$GSDIR/run_summary.txt" ]; then
  say "GS already present at $GSDIR — skipping."
else
  say "computing GS (35x35x120, N=92)..."
  ( cd gs && GS_CKPT=$GSDIR inq-run ) >>"$LOG" 2>&1 || fail "ground state"
  say "GS done -> $GSDIR"; notify "extkin_plateau_E100: GS done" "starting CAP run"
fi

# --- 2. CAP run (two-sided eta=-1.0, extkin every step) ----------------------
say "CAP run (eta=-1.0 two-sided, 1500 steps, extkin every step)..."
( cd wp && WP_GS_DIR=$GSDIR WP_OUT=cap inq-run ) >>"$LOG" 2>&1 || fail "CAP run"
say "CAP run done"
touch EXTKIN_PLATEAU_DONE.txt

# --- 3. auto-build the study notebook ----------------------------------------
HYP=$SYS/hypotheses/extkin_plateau_E100
$PY "$HYP/build_extkin_plateau_report.py" >>"$LOG" 2>&1 || say "  notebook build FAILED (rerun manually)"

say "=== extkin_plateau_E100 COMPLETE ==="
notify "extkin_plateau_E100 COMPLETE" "CAP run done. Results: $(pwd)/wp/results/cap; notebook: $HYP/extkin_plateau_E100_study.ipynb"
