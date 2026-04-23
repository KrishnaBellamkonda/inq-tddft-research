#!/usr/bin/env bash
# run_all_wp_rt.sh — run all 7 jellium WP-RT simulations sequentially.
#
# Usage: bash run_all_wp_rt.sh
# Each run is built and executed by inq-run in its own directory.
# Output from each run is tee'd to a per-run log file.
# A summary is printed at the end showing which runs succeeded/failed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS=(
    run_01_base
    run_02_low_energy
    run_03_high_energy
    run_04_tilted_45
    run_05_wide_sigma
    run_06_narrow_sigma
    run_07_open_shell
)

declare -A STATUS

echo "========================================================"
echo "  jellium WP-RT sequential launcher"
echo "  $(date)"
echo "  Runs: ${#RUNS[@]}"
echo "========================================================"

for RUN in "${RUNS[@]}"; do
    RUN_DIR="$SCRIPT_DIR/$RUN"
    LOG="$RUN_DIR/run.log"

    echo ""
    echo "--------------------------------------------------------"
    echo "  Starting: $RUN  ($(date +%H:%M:%S))"
    echo "  Log:      $LOG"
    echo "--------------------------------------------------------"

    if [ ! -d "$RUN_DIR" ]; then
        echo "  ERROR: directory not found: $RUN_DIR"
        STATUS[$RUN]="SKIP (no dir)"
        continue
    fi

    if [ ! -f "$RUN_DIR/run.cpp" ]; then
        echo "  ERROR: run.cpp not found in $RUN_DIR"
        STATUS[$RUN]="SKIP (no run.cpp)"
        continue
    fi

    pushd "$RUN_DIR" > /dev/null
    if inq-run 2>&1 | tee "$LOG"; then
        STATUS[$RUN]="OK"
        echo "  DONE: $RUN  ($(date +%H:%M:%S))"
    else
        STATUS[$RUN]="FAILED (exit $?)"
        echo "  FAILED: $RUN — see $LOG"
    fi
    popd > /dev/null
done

echo ""
echo "========================================================"
echo "  Summary  ($(date))"
echo "========================================================"
for RUN in "${RUNS[@]}"; do
    printf "  %-28s  %s\n" "$RUN" "${STATUS[$RUN]:-UNKNOWN}"
done
echo "========================================================"
