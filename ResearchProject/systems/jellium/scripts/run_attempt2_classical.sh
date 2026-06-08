#!/usr/bin/env bash
# run_attempt2_classical.sh — attempt-2 relaunch for classical E=50 eV.
#
# Attempt 1 was killed externally at step 484/913 (08:44 on 2026-05-14).
# Partial attempt-1 data preserved under run_classical_n162_L50_E50/.
# This relaunch writes into run_classical_n162_L50_E50_attempt2/.
#
# GPU is read from CUDA_VISIBLE_DEVICES (default 1).
#
# Usage:
#   CUDA_VISIBLE_DEVICES=1 ./run_attempt2_classical.sh
set -euo pipefail

GPU=${CUDA_VISIBLE_DEVICES:-1}
ROOT=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium
VENV=/local/data/public/skcb2/tddft/venv/bin/activate
RUN=run_classical_n162_L50_E50_attempt2

echo "=== $(date -Iseconds): attempt-2 classical chain starting on GPU $GPU ==="
echo "    run: $RUN"

cd "$ROOT/$RUN"
echo "=== $(date -Iseconds): launching $RUN ==="
CUDA_VISIBLE_DEVICES=$GPU inq-run > full_run.stdout 2> full_run.stderr
echo "=== $(date -Iseconds): $RUN sim done, running analyse ==="
source "$VENV"
python3 analyse.py > analyse.stdout 2> analyse.stderr
echo "=== $(date -Iseconds): $RUN analyse done ==="
echo "=== $(date -Iseconds): attempt-2 classical chain complete ==="
