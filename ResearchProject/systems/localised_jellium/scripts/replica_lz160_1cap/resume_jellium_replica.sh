#!/usr/bin/env bash
# ============================================================================
# resume_jellium_replica.sh — resume the CAP leg from its rt_ckpt (last_step=3200)
# after the 2026-07-28 disk-full abort at step 4382, then run the no-CAP leg.
# Uses the EXISTING wp/run binary directly (NO inq-run: avoids any rebuild of a
# possibly-drifted run.cpp and the truncate-while-running failure mode).
# Same physics env as run_jellium_replica.sh. GPU 0. Emails each milestone.
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
JCAMP=$ROOT/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau
PY=$ROOT/venv/bin/python3
export INQ_SOURCE=$ROOT/inq-study CUDA_VISIBLE_DEVICES=0
export INQ_SHARE_PATH=${INQ_SHARE_PATH:-$ROOT/inq/install/share}
export PSEUDOPOD_SHARE_PATH=${PSEUDOPOD_SHARE_PATH:-$ROOT/inq/install/share/pseudopod}
export PATH="$ROOT/shared/bin:$PATH" PYTHONPATH=$ROOT/inq-stack/python
GSDIR=$ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n102_L25x25x160_w0p5_h0p5
LOG=$(pwd)/jellium_replica_resume.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }
fail(){ say "FAILED: $1"; notify "jellium replica RESUME FAILED: $1" "See $LOG"; exit 1; }

say "=== jellium replica RESUME (CAP from ckpt, then no-CAP) on GPU 0 ==="
notify "jellium replica RESUMED (GPU 0)" "CAP from step 3200 -> 8000, then no-CAP; log $LOG"

# --- 1. CAP leg: resume from rt_ckpt (last_step=3200) -----------------------
say "CAP resume (WP_RESUME=1, target 8000 steps, ckpt every 1600)..."
( cd wp && WP_GS_DIR=$GSDIR WP_CAP_ETA=-0.7 WP_OUT=cap WP_N_STEPS=8000 \
  WP_CKPT_EVERY=1600 WP_RESUME=1 ./run ) >>"$LOG" 2>&1 || fail "CAP resume"
say "CAP run done"; notify "jellium replica: CAP done (resumed)" "starting no-CAP run"

# --- 2. no-CAP leg: fresh 8000 steps ---------------------------------------
say "no-CAP run (WP_CAP_ETA=0, 8000 steps, ckpt every 1600)..."
( cd wp && WP_GS_DIR=$GSDIR WP_CAP_ETA=0 WP_OUT=nocap WP_N_STEPS=8000 \
  WP_CKPT_EVERY=1600 ./run ) >>"$LOG" 2>&1 || fail "no-CAP run"
say "no-CAP run done"

say "=== jellium replica COMPLETE (via resume) ==="
notify "jellium replica COMPLETE" "CAP (resumed from 3200) + no-CAP done (8000 steps each). Runs in .../replica_lz160_1cap/wp/results/{cap,nocap}. Next: analysis + two-scheme stopping."
say "done."
