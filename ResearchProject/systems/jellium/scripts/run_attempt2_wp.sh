#!/usr/bin/env bash
# run_attempt2_wp.sh — attempt-2 relaunch for WP E=50 eV on a single GPU.
#
# Attempt 1 was killed externally at step 347/913 (08:44 on 2026-05-14)
# at the same instant as the classical sibling — see
# docs/handovers/li_extensive_kick.md for diagnostic. The relaunched run
# lives in run_wp_n162_L50_E50_attempt2/ to keep attempt-1's partial
# results untouched for salvage.
#
# GPU is read from CUDA_VISIBLE_DEVICES (default 0).
#
# Usage:
#   CUDA_VISIBLE_DEVICES=0 ./run_attempt2_wp.sh
set -euo pipefail

GPU=${CUDA_VISIBLE_DEVICES:-0}
ROOT=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium
VENV=/local/data/public/skcb2/tddft/venv/bin/activate
RUN=run_wp_n162_L50_E50_attempt2

echo "=== $(date -Iseconds): attempt-2 WP chain starting on GPU $GPU ==="
echo "    run: $RUN"

cd "$ROOT/$RUN"
echo "=== $(date -Iseconds): launching $RUN ==="
CUDA_VISIBLE_DEVICES=$GPU inq-run > full_run.stdout 2> full_run.stderr
echo "=== $(date -Iseconds): $RUN sim done, running analyse ==="
source "$VENV"
python3 analyse.py > analyse.stdout 2> analyse.stderr
echo "=== $(date -Iseconds): $RUN analyse done ==="
echo "=== $(date -Iseconds): attempt-2 WP chain complete ==="
