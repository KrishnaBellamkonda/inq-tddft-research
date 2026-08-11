#!/usr/bin/env bash
# shared/bin/run-lovelace-test.sh — graphene lovelace_test campaign, local Lovelace GPUs
#
# USAGE (detach so it survives session end):
#   setsid nohup bash shared/bin/run-lovelace-test.sh > lovelace_test.log 2>&1 &
#
# Phase 0  GS — bilayer graphene, Lz=90, dx=0.5 Bohr (standard inq, GPU 0, foreground)
# Phase 1  WP (150 eV, GPU 0) + classical (25 eV, GPU 1) — inq-study for CAP, parallel
#
# Env knobs:
#   SKIP_GS=1         skip Phase 0 if GS already complete
#   GPUS="0 1"        GPU ids for [wp, classical]
#   WP_E_EV=150       WP energy in eV (defaults to 150)
#   CL_E_EV=25        classical energy in eV (defaults to 25)
set -uo pipefail

REPO=/local/data/public/skcb2/tddft
GR=$REPO/ResearchProject/systems/graphene
SCR=$GR/scripts/lovelace_test

export PATH="$REPO/shared/bin:$PATH"
export TDDFT_ROOT="$REPO"
export INQ_BUILD_JOBS="${INQ_BUILD_JOBS:-8}"

# INQ config (will be overridden to inq-study for dynamics)
export INQ_SOURCE="$REPO/inq"
export INQ_SHARE_PATH="$REPO/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO/inq/install/share/pseudopod"
# FetchContent cache (avoids network on compute — reuse engine build's deps)
export INQ_DEPS_CACHE="$REPO/inq/build/_deps"

GPUS="${GPUS:-0 1}"
GPU_WP=$(echo "$GPUS" | awk '{print $1}')
GPU_CL=$(echo "$GPUS" | awk '{print $2}')

# k0 = sqrt(2 * E_eV / 27.211386)
eV_to_k0() { python3 -c "import math; print(round(math.sqrt(2*$1/27.211386),4))"; }
WP_E_EV="${WP_E_EV:-150}"; WP_K0=$(eV_to_k0 "$WP_E_EV")
CL_E_EV="${CL_E_EV:-25}";  CL_K0=$(eV_to_k0 "$CL_E_EV")

GS_DIR="$SCR/gs/results/bi_dx0p5_Lz90"
WP_OUT="bi_E${WP_E_EV}_sigma4_dx0p5"
CL_OUT="bi_E${CL_E_EV}_sigma4_dx0p5"

log() { echo "== $(date '+%H:%M:%S') $*"; }

echo "======================================================================"
echo "== lovelace_test START $(date) host=$(hostname)"
echo "== GPUs: WP=$GPU_WP  CL=$GPU_CL"
echo "== WP ${WP_E_EV} eV (k0=${WP_K0})   Classical ${CL_E_EV} eV (k0=${CL_K0})"
echo "== Disk: $(df -h /local/data 2>/dev/null | tail -1 | awk '{print $4}') free on /local/data"
echo "======================================================================"

# ── Phase 0: Ground State ────────────────────────────────────────────────────
if [ "${SKIP_GS:-0}" = "1" ] || \
   { [ -f "$GS_DIR/run_summary.txt" ] && grep -q "run_completed = true" "$GS_DIR/run_summary.txt"; }; then
    log "GS already complete at $GS_DIR — skipping"
else
    log "PHASE 0: building + running GS (GPU $GPU_WP, standard inq/)"
    ( cd "$SCR/gs" && CUDA_VISIBLE_DEVICES="$GPU_WP" inq-run )
    log "GS done"
fi

if ! grep -q "run_completed = true" "$GS_DIR/run_summary.txt" 2>/dev/null; then
    echo "FATAL: GS did not complete — $GS_DIR/run_summary.txt" >&2; exit 1; fi

# ── Phase 1: Build WP + classical against inq-study ──────────────────────────
export INQ_SOURCE="$REPO/inq-study"
export INQ_DEPS_CACHE="$REPO/inq-study/build-gpu/_deps"
log "Switched INQ_SOURCE → inq-study (needed for perturbations::absorbing CAP)"

# Build WP binary (run exits immediately on missing GS, giving us just the compiled binary)
if [ ! -x "$SCR/wp/run" ]; then
    log "Building WP binary..."
    ( cd "$SCR/wp" && \
      LJ_OUT=_build_check LJ_GS_DIR=/tmp/lovelace_nonexistent_gs \
      CUDA_VISIBLE_DEVICES="$GPU_WP" inq-run ) || true
fi
[ -x "$SCR/wp/run" ] || { echo "FATAL: $SCR/wp/run not built"; exit 2; }

# Build classical binary
if [ ! -x "$SCR/classical/run" ]; then
    log "Building classical binary..."
    ( cd "$SCR/classical" && \
      LJ_OUT=_build_check LJ_GS_DIR=/tmp/lovelace_nonexistent_gs \
      CUDA_VISIBLE_DEVICES="$GPU_CL" inq-run ) || true
fi
[ -x "$SCR/classical/run" ] || { echo "FATAL: $SCR/classical/run not built"; exit 2; }

log "Both binaries ready."

# ── Launch WP (GPU $GPU_WP) ───────────────────────────────────────────────────
WP_RESULT="$SCR/wp/results/$WP_OUT"
mkdir -p "$WP_RESULT"
log "Launching WP [$WP_E_EV eV, k0=$WP_K0] on GPU $GPU_WP → $WP_OUT"
( cd "$SCR/wp" && setsid env \
    LJ_OUT="$WP_OUT"    LJ_GS_DIR="$GS_DIR"  \
    LJ_K0="$WP_K0"      LJ_SIGMA=4.0          \
    LJ_LAUNCH_Z=-19     LJ_N_STEPS=500         \
    LJ_DT=0.04          LJ_CAP_ETA=-1.0        \
    LJ_CAP_L=10.0       LJ_SAVE_EVERY=10       \
    INQ_SHARE_PATH="$INQ_SHARE_PATH"           \
    PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
    TDDFT_ROOT="$REPO"  \
    CUDA_VISIBLE_DEVICES="$GPU_WP" \
    ./run > "$WP_RESULT/run.log" 2>&1 ) &
WP_PID=$!
echo "   WP PID $WP_PID   log: $WP_RESULT/run.log"

# ── Launch classical (GPU $GPU_CL) ────────────────────────────────────────────
CL_RESULT="$SCR/classical/results/$CL_OUT"
mkdir -p "$CL_RESULT"
log "Launching classical [$CL_E_EV eV, k0=$CL_K0] on GPU $GPU_CL → $CL_OUT"
( cd "$SCR/classical" && setsid env \
    LJ_OUT="$CL_OUT"    LJ_GS_DIR="$GS_DIR"  \
    LJ_K0="$CL_K0"      LJ_SIGMA=4.0          \
    LJ_LAUNCH_Z=-19     LJ_N_STEPS=500         \
    LJ_DT=0.04          LJ_CAP_ETA=-1.0        \
    LJ_CAP_L=10.0       LJ_SAVE_EVERY=10       \
    INQ_SHARE_PATH="$INQ_SHARE_PATH"           \
    PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
    TDDFT_ROOT="$REPO"  \
    CUDA_VISIBLE_DEVICES="$GPU_CL" \
    ./run > "$CL_RESULT/run.log" 2>&1 ) &
CL_PID=$!
echo "   Classical PID $CL_PID   log: $CL_RESULT/run.log"

log "Both launched. Waiting for completion..."

wait $WP_PID;  WP_ST=$?
wait $CL_PID;  CL_ST=$?

echo "======================================================================"
echo "== lovelace_test DONE $(date)"
printf "== WP        exit=%s\n" "$WP_ST"
printf "== Classical exit=%s\n" "$CL_ST"
[ -f "$WP_RESULT/run_summary.txt" ] && \
    grep -E "run_completed|total_time_au|wall_time" "$WP_RESULT/run_summary.txt" | \
    sed "s/^/   WP: /"
[ -f "$CL_RESULT/run_summary.txt" ] && \
    grep -E "run_completed|total_time_au|wall_time" "$CL_RESULT/run_summary.txt" | \
    sed "s/^/   CL: /"
echo "======================================================================"
