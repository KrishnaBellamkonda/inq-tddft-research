#!/usr/bin/env bash
# ============================================================================
# launch_parallel_wp_classical.sh — run the WP and classical projectile
# simulations concurrently, one per GPU, after the GS checkpoint exists.
#
# Usage:
#   ./launch_parallel_wp_classical.sh [--dry]
#
# Preconditions (script will check):
#   - GS checkpoint dir present:
#       checkpoints/gs_L50_cubic_N162_dx0p248/
#   - Both A30s available (no CUDA processes on either)
#
# What it does:
#   - WP run on GPU 0 (CUDA_VISIBLE_DEVICES=0), output in
#       run_wp_e1500_L50_cubic/full_run.stdout
#   - Classical run on GPU 1 (CUDA_VISIBLE_DEVICES=1), output in
#       run_classical_e1500_L50_cubic/full_run.stdout
#   - Each is launched detached; this script exits as soon as both are
#     started, leaving them to run in background.
#
# Each run is expected to take ~30-60 min wall on a single A30.
# ============================================================================
set -euo pipefail

REPO=/local/data/public/skcb2/tddft
JELLIUM=$REPO/ResearchProject/systems/jellium
GS_DIR=$JELLIUM/checkpoints/gs_L50_cubic_N162_dx0p248
WP_DIR=$JELLIUM/run_wp_e1500_L50_cubic
CL_DIR=$JELLIUM/run_classical_e1500_L50_cubic

DRY=0
[[ "${1:-}" == "--dry" ]] && DRY=1

echo "=== launch_parallel_wp_classical ==="

# 1. Check GS checkpoint exists.
if [[ ! -d "$GS_DIR" ]]; then
    echo "FATAL: GS checkpoint missing: $GS_DIR"
    echo "Run save_gs/gs_L50_cubic_N162_dx0p248/run.cpp first."
    exit 1
fi
echo "[ok] GS checkpoint present: $GS_DIR"

# 2. Check both run dirs have run.cpp.
[[ -f "$WP_DIR/run.cpp" ]] || { echo "FATAL: $WP_DIR/run.cpp missing"; exit 1; }
[[ -f "$CL_DIR/run.cpp" ]] || { echo "FATAL: $CL_DIR/run.cpp missing"; exit 1; }
echo "[ok] Both run.cpp files present"

# 3. Check GPU 0 and GPU 1 are idle (no compute processes).
GPU_PROCS=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | wc -l)
if [[ "$GPU_PROCS" -gt 0 ]]; then
    echo "WARN: GPU(s) currently in use by $GPU_PROCS process(es):"
    nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv 2>/dev/null
    echo "Continue anyway? (Ctrl-C to abort, Enter to proceed)"
    [[ "$DRY" == "0" ]] && read -r
fi

# 4. Launch.
if [[ "$DRY" == "1" ]]; then
    echo "[dry-run] would launch:"
    echo "  cd $WP_DIR && CUDA_VISIBLE_DEVICES=0 inq-run > full_run.stdout 2>&1 &"
    echo "  cd $CL_DIR && CUDA_VISIBLE_DEVICES=1 inq-run > full_run.stdout 2>&1 &"
    exit 0
fi

echo "[launching] WP run on GPU 0 ..."
( cd "$WP_DIR" && CUDA_VISIBLE_DEVICES=0 nohup inq-run > full_run.stdout 2>&1 & )

echo "[launching] Classical run on GPU 1 ..."
( cd "$CL_DIR" && CUDA_VISIBLE_DEVICES=1 nohup inq-run > full_run.stdout 2>&1 & )

sleep 2

echo
echo "=== launched ==="
nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader
echo
echo "Logs:"
echo "  $WP_DIR/full_run.stdout"
echo "  $CL_DIR/full_run.stdout"
echo
echo "Track progress with:"
echo "  tail -f $WP_DIR/full_run.stdout"
echo "  tail -f $CL_DIR/full_run.stdout"
