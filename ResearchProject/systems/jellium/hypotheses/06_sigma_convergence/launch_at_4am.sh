#!/usr/bin/env bash
# Detached launcher: sleeps until the next 04:00, then runs the sweep.
# Started via:  setsid nohup ./launch_at_4am.sh </dev/null >launch.log 2>&1 &
set -u
HERE="/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/hypotheses/06_sigma_convergence"
TARGET=$(date -d '04:00' +%s)
NOW=$(date +%s)
# if 04:00 already passed today, aim for tomorrow 04:00
if [ "$TARGET" -le "$NOW" ]; then TARGET=$(date -d 'tomorrow 04:00' +%s); fi
SECS=$(( TARGET - NOW ))
echo "[$(date -u +%FT%TZ)] launcher armed; sleeping ${SECS}s until $(date -d @${TARGET})"
sleep "$SECS"
echo "[$(date -u +%FT%TZ)] firing orchestrator"
exec bash "$HERE/orchestrate_sigma_sweep.sh"
