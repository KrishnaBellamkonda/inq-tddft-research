#!/usr/bin/env bash
# ============================================================================
# run_jellium_replica.sh — localised-jellium replica (LZ=160, ONE-SIDED CAP 20 Bohr,
# 8000 steps, r_s=3.32, WP sigma=1/E=100eV). GPU 0. inq-study.
# Sequence: GS -> CAP run -> no-CAP run. 5 CHECKPOINTS per WP run (WP_CKPT_EVERY=1600
# -> interior 1600/3200/4800/6400 + final = 5). Resumable (*_RESUME=1). Emails each
# milestone. All decomposed energies tracked (energies.csv); the per-particle
# normalization fix is applied in analysis (E_ext = E_reported - e_kin_ha*(1-norm)).
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
LOG=$(pwd)/jellium_replica.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }
fail(){ say "FAILED: $1"; notify "jellium replica FAILED: $1" "See $LOG"; exit 1; }

say "=== jellium replica (LZ=160, 1-sided CAP, 8000 steps) on GPU 0 ==="
notify "jellium replica STARTED (GPU 0)" "GS -> CAP -> noCAP; 5 checkpoints each; log $LOG"

# --- 1. Ground state (new LZ=160 box) --------------------------------------
if [ -f "$GSDIR/ground_state.gs" ] || [ -d "$GSDIR/kpin0000000000" ]; then
  say "GS already present at $GSDIR — skipping."
else
  say "computing GS (LZ=160)..."
  ( cd gs && GS_CKPT=$GSDIR inq-run ) >>"$LOG" 2>&1 || fail "ground state"
  say "GS done -> $GSDIR"; notify "jellium replica: GS done" "starting CAP run (22h, 5 checkpoints)"
fi

# --- 2. CAP run (8000 steps, 5 checkpoints) --------------------------------
say "CAP run (WP_CAP_ETA=-0.7, 8000 steps, ckpt every 1600)..."
( cd wp && WP_GS_DIR=$GSDIR WP_CAP_ETA=-0.7 WP_OUT=cap WP_N_STEPS=8000 WP_CKPT_EVERY=1600 inq-run ) >>"$LOG" 2>&1 || fail "CAP run"
say "CAP run done"; notify "jellium replica: CAP done" "starting no-CAP run"

# --- 3. no-CAP run (8000 steps, 5 checkpoints) -----------------------------
say "no-CAP run (WP_CAP_ETA=0, 8000 steps, ckpt every 1600)..."
( cd wp && WP_GS_DIR=$GSDIR WP_CAP_ETA=0 WP_OUT=nocap WP_N_STEPS=8000 WP_CKPT_EVERY=1600 ./run ) >>"$LOG" 2>&1 || fail "no-CAP run"
say "no-CAP run done"

say "=== jellium replica COMPLETE ==="
notify "jellium replica COMPLETE" "CAP + no-CAP done (8000 steps each). Runs in .../replica_lz160_1cap/wp/results/{cap,nocap}. Next: analysis + two-scheme stopping."
say "done."
