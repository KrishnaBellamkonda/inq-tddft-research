#!/usr/bin/env bash
# Re-runs the full postprocess pipeline (all phases, --rebuild) on every
# completed run + every hypothesis comparison. Logs to scripts/repostprocess.log.

cd "$(dirname "$0")/.."  # systems/coronene/

# Activate the project venv directly (don't source ~/.bashrc — it tries to
# read PS1 which is unbound in a non-interactive shell).
source /local/data/public/skcb2/tddft/venv/bin/activate

LOG=scripts/repostprocess.log
> "$LOG"

RUNS=( run_base run_E30 run_E800 run_s0p33 run_s3 run_E800_s0p33
       run_E30_s3 run_b18_35x35x80 run_b6_35x35x80 run_35x35x40 )

echo "[$(date '+%H:%M:%S')] starting per-run postprocess (--with-paraview, --rebuild)" >>"$LOG"
for run in "${RUNS[@]}"; do
    echo "[$(date '+%H:%M:%S')] === $run ===" >>"$LOG"
    if python3 scripts/coronene_postprocess.py run \
            --results "$run/results" --run-name "$run" \
            --rebuild --with-paraview \
            >>"$LOG" 2>&1; then
        echo "[$(date '+%H:%M:%S')] $run OK" >>"$LOG"
    else
        echo "[$(date '+%H:%M:%S')] $run FAILED" >>"$LOG"
    fi
done

echo "[$(date '+%H:%M:%S')] starting hypothesis comparisons (--rebuild)" >>"$LOG"
RP=/local/data/public/skcb2/tddft/ResearchProject/systems/coronene

run_h() {
    local h=$1; shift
    echo "[$(date '+%H:%M:%S')] === $h ===" >>"$LOG"
    python3 scripts/coronene_postprocess.py hypothesis \
        --hypothesis-dir "hypotheses/$h" \
        --runs "$@" --rebuild >>"$LOG" 2>&1
}

run_h 00_base                       "run_base=$RP/run_base"
run_h 01_wp_energy_spread           "run_E30=$RP/run_E30" "run_base=$RP/run_base" "run_E800=$RP/run_E800"
run_h 02_wp_sigma_spread            "run_s0p33=$RP/run_s0p33" "run_base=$RP/run_base" "run_s3=$RP/run_s3"
run_h 03_fast_projectile_classical  "run_E800_s0p33=$RP/run_E800_s0p33" "run_base=$RP/run_base"
run_h 04_electron_capture           "run_E30_s3=$RP/run_E30_s3" "run_base=$RP/run_base"
run_h 05_box_length_and_distance \
    "run_b18_35x35x80=$RP/run_b18_35x35x80" \
    "run_b6_35x35x80=$RP/run_b6_35x35x80" \
    "run_35x35x40=$RP/run_35x35x40" \
    "run_base=$RP/run_base"

echo "[$(date '+%H:%M:%S')] all done." >>"$LOG"
