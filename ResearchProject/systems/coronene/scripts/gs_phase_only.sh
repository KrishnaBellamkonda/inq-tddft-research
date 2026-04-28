#!/usr/bin/env bash
# Re-run only the gs phase per run, so eigenvalue viz lands in
# results/analysis/ground_state/.
# After Phase 3 retrofit_eigenvalues.py populates the eigenvalue CSVs.

cd "$(dirname "$0")/.."  # systems/coronene/

source /local/data/public/skcb2/tddft/venv/bin/activate

LOG=scripts/gs_phase_only.log
> "$LOG"

if [[ -n "$1" ]]; then
    RUNS=( "$@" )
else
    RUNS=( run_base run_E30 run_E800 run_s0p33 run_s3 run_E800_s0p33
           run_E30_s3 run_b18_35x35x80 run_b6_35x35x80 run_35x35x40 )
fi

echo "[$(date '+%H:%M:%S')] starting gs phase for ${#RUNS[@]} runs" >>"$LOG"
for run in "${RUNS[@]}"; do
    echo "[$(date '+%H:%M:%S')] === $run ===" >>"$LOG"
    if python3 scripts/coronene_postprocess.py run \
            --results "$run/results" --run-name "$run" \
            --phases gs --rebuild \
            >>"$LOG" 2>&1; then
        echo "[$(date '+%H:%M:%S')] $run gs OK" >>"$LOG"
    else
        echo "[$(date '+%H:%M:%S')] $run gs FAILED (exit=$?)" >>"$LOG"
    fi
done

echo "[$(date '+%H:%M:%S')] all done." >>"$LOG"
