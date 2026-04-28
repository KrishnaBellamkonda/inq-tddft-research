#!/usr/bin/env bash
# Re-run only the 6 hypothesis comparisons (no per-run postprocess).
# Used in Phase 3 after Branch-3 reruns finished + auto-postprocess populated
# every run's analysis/ tree.

cd "$(dirname "$0")/.."  # systems/coronene/

source /local/data/public/skcb2/tddft/venv/bin/activate

LOG=scripts/hypotheses_only.log
> "$LOG"

RP=/local/data/public/skcb2/tddft/ResearchProject/systems/coronene

run_h() {
    local h=$1; shift
    echo "[$(date '+%H:%M:%S')] === $h ===" >>"$LOG"
    if python3 scripts/coronene_postprocess.py hypothesis \
            --hypothesis-dir "hypotheses/$h" \
            --runs "$@" --rebuild >>"$LOG" 2>&1; then
        echo "[$(date '+%H:%M:%S')] $h OK" >>"$LOG"
    else
        echo "[$(date '+%H:%M:%S')] $h FAILED (exit=$?)" >>"$LOG"
    fi
}

echo "[$(date '+%H:%M:%S')] starting 6 hypothesis comparisons" >>"$LOG"

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
