#!/usr/bin/env bash
# run_chain_classical.sh — sequential chain for classical-electron runs.
#
# Runs (in order): classical 50 eV → analyse.
# (Classical 300 eV dropped 2026-05-14 to fit a 6-hour deadline —
#  E=50 preserved because it is the Bragg-peak-onset matched pair
#  for WP E=50, the most scientifically informative comparison
#  in this sweep.)
# Each step is: build (inq-run) → propagate → run-local analyse.py.
#
# GPU is read from CUDA_VISIBLE_DEVICES (default 1 — the GPU assigned to
# classical sims by the regime-classification plan
# (docs/plans/jellium-regime-constrained-simulations.md §8.2)).
#
# `set -euo pipefail` stops the chain at the first error: a failed build
# or failed analyse will halt before contaminating subsequent runs.
#
# Usage:
#   CUDA_VISIBLE_DEVICES=1 ./run_chain_classical.sh
# or for background:
#   CUDA_VISIBLE_DEVICES=1 nohup ./run_chain_classical.sh \
#       > chain_classical.log 2>&1 &
set -euo pipefail

GPU=${CUDA_VISIBLE_DEVICES:-1}
ROOT=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium
VENV=/local/data/public/skcb2/tddft/venv/bin/activate

RUNS=(
    run_classical_n162_L50_E50
)

echo "=== $(date -Iseconds): classical chain starting on GPU $GPU ==="
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

echo "=== $(date -Iseconds): classical chain complete ==="
