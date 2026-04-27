#!/usr/bin/env bash
# Run all diagnostic ground-state simulations, up to 2 in parallel when GPUs are free.
# Each run: inq-run (GPU build) followed by python analysis.py.
# Logs per run to <run_dir>/run.log.
#
# Usage: ./run_all_diagnoses.sh [--max-parallel N]
#
# Requires: inq-run, python (with inqview), nvidia-smi (for GPU count detection).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RUNS=(
    run_01_tight_scf
    run_02_benzene
    run_03_coronene_hardcoded
    run_04_graphene
    run_05_quarter_coronene
)

# Default to 2 parallel jobs; override with --max-parallel N
MAX_PARALLEL=2
while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# Detect available GPUs and cap MAX_PARALLEL to what's present
if command -v nvidia-smi &>/dev/null; then
    N_GPU=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l || echo 1)
else
    echo "WARNING: nvidia-smi not found. Assuming 1 GPU and running serially."
    N_GPU=1
fi

if [ "$N_GPU" -lt "$MAX_PARALLEL" ]; then
    echo "[$(date '+%H:%M:%S')] WARNING: $N_GPU GPU(s) found, requested $MAX_PARALLEL parallel jobs."
    echo "[$(date '+%H:%M:%S')]          Running at most $N_GPU job(s) in parallel."
    MAX_PARALLEL=$N_GPU
fi

echo "[$(date '+%H:%M:%S')] GPUs found: $N_GPU  |  MAX_PARALLEL: $MAX_PARALLEL"
echo "[$(date '+%H:%M:%S')] Runs to execute: ${RUNS[*]}"
echo "---"

# Run one simulation in a subshell (background).
# Writes stdout+stderr to <run_dir>/run.log.
run_one() {
    local run_name="$1"
    local run_dir="$SCRIPT_DIR/$run_name"
    local log="$run_dir/run.log"

    echo "[$(date '+%H:%M:%S')] Starting: $run_name  (log: $log)"
    (
        cd "$run_dir"
        {
            echo "=== $run_name started at $(date) ==="
            inq-run
            echo "=== simulation done at $(date) ==="
            echo "--- running analysis.py ---"
            python analysis.py
            echo "=== $run_name complete at $(date) ==="
        } > "$log" 2>&1
    ) &
}

# Job slot tracking: count background children with `jobs -p`
running=0

for run in "${RUNS[@]}"; do
    # Wait for a slot: poll until a background job finishes
    if [ "$running" -ge "$MAX_PARALLEL" ]; then
        # bash 4.3+: wait for any one child
        if wait -n 2>/dev/null; then
            echo "[$(date '+%H:%M:%S')] A job finished (OK)"
        else
            # wait -n not available or a job failed — still a slot freed
            echo "[$(date '+%H:%M:%S')] A job finished"
        fi
        running=$((running - 1))
    fi

    run_one "$run"
    running=$((running + 1))
done

# Drain remaining jobs
echo "[$(date '+%H:%M:%S')] Waiting for remaining $running job(s)..."
wait
echo "[$(date '+%H:%M:%S')] All diagnostic runs complete."
echo ""
echo "Log files:"
for run in "${RUNS[@]}"; do
    echo "  $SCRIPT_DIR/$run/run.log"
done
