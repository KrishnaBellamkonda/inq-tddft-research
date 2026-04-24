#!/usr/bin/env bash
# run_all_wp_rt.sh — sequentially run all WP real-time propagation simulations
#
# Order:
#   1. Tutorial/free-propagation-wp-rt/  (7 runs, ~1-2 h each)
#   2. ResearchProject/jellium/jellium-wp-rt/  (6 runs, ~3-5 h each)
#   3. ResearchProject/systems/coronene/coronene-wp-rt/  (6 runs, ~10-14 h each)
#
# Each run's output is tee'd to results/run.log in the run directory.
# A run that fails is recorded; the script continues to the next run.
#
# Usage:
#   ./run_all_wp_rt.sh             # run everything
#   ./run_all_wp_rt.sh free        # only free-propagation runs
#   ./run_all_wp_rt.sh jellium     # only jellium runs
#   ./run_all_wp_rt.sh coronene    # only coronene runs

set -uo pipefail

export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILTER="${1:-all}"

# ── Run lists ─────────────────────────────────────────────────────────────────

FREE_RUNS=(
    "Tutorial/free-propagation-wp-rt/run_01_base"
    "Tutorial/free-propagation-wp-rt/run_02_low_momentum"
    "Tutorial/free-propagation-wp-rt/run_03_high_momentum"
    "Tutorial/free-propagation-wp-rt/run_04_tilted_45"
    "Tutorial/free-propagation-wp-rt/run_05_transverse_x"
    "Tutorial/free-propagation-wp-rt/run_06_wide_sigma"
    "Tutorial/free-propagation-wp-rt/run_07_narrow_sigma"
)

JELLIUM_RUNS=(
    "ResearchProject/jellium/jellium-wp-rt/run_01_base"
    "ResearchProject/jellium/jellium-wp-rt/run_02_low_energy"
    "ResearchProject/jellium/jellium-wp-rt/run_03_high_energy"
    "ResearchProject/jellium/jellium-wp-rt/run_04_tilted_45"
    "ResearchProject/jellium/jellium-wp-rt/run_05_wide_sigma"
    "ResearchProject/jellium/jellium-wp-rt/run_06_narrow_sigma"
)

CORONENE_RUNS=(
    "ResearchProject/systems/coronene/coronene-wp-rt/run_01_d635_base"
    "ResearchProject/systems/coronene/coronene-wp-rt/run_02_d3"
    "ResearchProject/systems/coronene/coronene-wp-rt/run_03_d10"
    "ResearchProject/systems/coronene/coronene-wp-rt/run_04_d15"
    "ResearchProject/systems/coronene/coronene-wp-rt/run_05_d20"
    "ResearchProject/systems/coronene/coronene-wp-rt/run_06_projectile"
)

# ── Select runs based on filter ───────────────────────────────────────────────

RUNS=()
case "$FILTER" in
    free)     RUNS=("${FREE_RUNS[@]}") ;;
    jellium)  RUNS=("${JELLIUM_RUNS[@]}") ;;
    coronene) RUNS=("${CORONENE_RUNS[@]}") ;;
    all)      RUNS=("${FREE_RUNS[@]}" "${JELLIUM_RUNS[@]}" "${CORONENE_RUNS[@]}") ;;
    *)
        echo "Unknown filter: $FILTER  (valid: free | jellium | coronene | all)"
        exit 1
        ;;
esac

N_TOTAL=${#RUNS[@]}
echo "============================================================"
echo " WP RT propagation — $N_TOTAL runs (filter=$FILTER)"
echo " Started: $(date)"
echo "============================================================"

# ── Execute runs ──────────────────────────────────────────────────────────────

PASSED=()
FAILED=()

for i in "${!RUNS[@]}"; do
    run_rel="${RUNS[$i]}"
    run_abs="$SCRIPT_DIR/$run_rel"
    run_num=$((i + 1))
    run_name="$(basename "$run_rel")"

    echo ""
    echo "------------------------------------------------------------"
    echo " [$run_num/$N_TOTAL] $run_rel"
    echo " Started: $(date)"
    echo "------------------------------------------------------------"

    if [[ ! -d "$run_abs" ]]; then
        echo "ERROR: directory not found: $run_abs"
        FAILED+=("$run_rel (directory missing)")
        continue
    fi

    mkdir -p "$run_abs/results"
    LOG="$run_abs/results/run.log"

    # Run inq-run from within the run directory; tee output to log
    if (cd "$run_abs" && inq-run 2>&1 | tee "$LOG"); then
        echo " [$run_num/$N_TOTAL] PASSED: $run_name  ($(date))"
        PASSED+=("$run_rel")
    else
        exit_code=$?
        echo " [$run_num/$N_TOTAL] FAILED (exit $exit_code): $run_name  ($(date))"
        FAILED+=("$run_rel (exit $exit_code)")
    fi
done

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo " Summary — $(date)"
echo "============================================================"
echo " Passed: ${#PASSED[@]} / $N_TOTAL"
echo " Failed: ${#FAILED[@]} / $N_TOTAL"

if [[ ${#PASSED[@]} -gt 0 ]]; then
    echo ""
    echo " PASSED:"
    for r in "${PASSED[@]}"; do echo "   ✓ $r"; done
fi

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    echo " FAILED:"
    for r in "${FAILED[@]}"; do echo "   ✗ $r"; done
    echo ""
    echo " To re-run a failed job:"
    echo "   cd <run_dir> && inq-run"
    exit 1
fi

echo ""
echo " All runs completed successfully."
