#!/usr/bin/env bash
# run_chain_wp.sh — sequential chain for wave-packet runs.
#
# Runs (in order): WP 300 eV → analyse → WP 50 eV → analyse.
# (WP 25 eV stretch dropped 2026-05-14 to fit a 6-hour deadline —
#  the lowest-velocity sim is the longest tail of the chain.)
#
# GPU is read from CUDA_VISIBLE_DEVICES (default 0 — the GPU assigned to
# WP sims by the regime-classification plan (§8.3)).
#
# Usage:
#   CUDA_VISIBLE_DEVICES=0 ./run_chain_wp.sh
# or for background:
#   CUDA_VISIBLE_DEVICES=0 nohup ./run_chain_wp.sh \
#       > chain_wp.log 2>&1 &
set -euo pipefail

GPU=${CUDA_VISIBLE_DEVICES:-0}
ROOT=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium
VENV=/local/data/public/skcb2/tddft/venv/bin/activate

RUNS=(
    run_wp_n162_L50_E300
    run_wp_n162_L50_E50
)

echo "=== $(date -Iseconds): WP chain starting on GPU $GPU ==="
echo "    runs: ${RUNS[*]}"

for run in "${RUNS[@]}"; do
    cd "$ROOT/$run"
    echo "=== $(date -Iseconds): launching $run ==="
    CUDA_VISIBLE_DEVICES=$GPU inq-run > full_run.stdout 2> full_run.stderr
    echo "=== $(date -Iseconds): $run sim done, running analyse ==="
    source "$VENV"
    python3 analyse.py > analyse.stdout 2> analyse.stderr
    echo "=== $(date -Iseconds): $run analyse done ==="
done

echo "=== $(date -Iseconds): WP chain complete ==="
