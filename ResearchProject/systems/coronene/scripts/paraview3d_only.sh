#!/usr/bin/env bash
# Run only the paraview_3d phase for every run in the framework.
# Used in Phase 3 closeout: auto-postprocess produced every other phase
# but skipped paraview_3d (heavy; needs --with-paraview).

cd "$(dirname "$0")/.."  # systems/coronene/

source /local/data/public/skcb2/tddft/venv/bin/activate

LOG=scripts/paraview3d_only.log
> "$LOG"

RUNS=( run_base run_E30 run_E800 run_s0p33 run_s3 run_E800_s0p33
       run_E30_s3 run_b18_35x35x80 run_b6_35x35x80 run_35x35x40 )

echo "[$(date '+%H:%M:%S')] starting paraview_3d phase for 10 runs" >>"$LOG"
for run in "${RUNS[@]}"; do
    echo "[$(date '+%H:%M:%S')] === $run ===" >>"$LOG"
    if python3 scripts/coronene_postprocess.py run \
            --results "$run/results" --run-name "$run" \
            --phases paraview_3d --rebuild --with-paraview \
            >>"$LOG" 2>&1; then
        echo "[$(date '+%H:%M:%S')] $run paraview_3d OK" >>"$LOG"
    else
        echo "[$(date '+%H:%M:%S')] $run paraview_3d FAILED (exit=$?)" >>"$LOG"
    fi
done

echo "[$(date '+%H:%M:%S')] all done." >>"$LOG"
