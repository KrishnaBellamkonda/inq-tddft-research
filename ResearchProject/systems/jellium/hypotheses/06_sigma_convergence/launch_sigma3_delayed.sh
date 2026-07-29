#!/usr/bin/env bash
# ============================================================================
# Delayed detached launcher for the σ=3 ladder.
#
# Sleeps DELAY_SECONDS (default 7200 = 2 h) then execs orchestrate_sigma3.sh.
# Detached via setsid+nohup so it survives the session that armed it.
#
# Arm:
#   setsid nohup bash launch_sigma3_delayed.sh </dev/null >launch_sigma3.log 2>&1 &
# ============================================================================
set -u
HERE="/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/hypotheses/06_sigma_convergence"
DELAY_SECONDS="${DELAY_SECONDS:-7200}"
echo "[$(date)] armed; sleeping ${DELAY_SECONDS}s until $(date -d "+${DELAY_SECONDS} seconds")"
sleep "$DELAY_SECONDS"
echo "[$(date)] waking; launching σ=3 orchestrator"
exec bash "$HERE/orchestrate_sigma3.sh"
