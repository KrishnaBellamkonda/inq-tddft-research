#!/usr/bin/env bash
# Watches dispatch.log for "finished run_* exit=0" lines and runs the
# coronene postprocess CLI on each one as it completes. Writes its own log
# at scripts/auto_postprocess.log.

set -u
cd "$(dirname "$0")/.."  # systems/coronene/

source /local/data/public/skcb2/tddft/venv/bin/activate

LOG=scripts/auto_postprocess.log
DISP=scripts/dispatch.log
DONE_LIST=scripts/auto_postprocess_done.txt
touch "$DONE_LIST"

echo "[$(date '+%H:%M:%S')] auto_postprocess: watching $DISP" >>"$LOG"

while true; do
  # Find every run_* that exit=0'd
  finished=$(grep -E 'finished (run_[A-Za-z0-9_]+) exit=0' "$DISP" 2>/dev/null \
             | sed -E 's/.* finished (run_[A-Za-z0-9_]+) exit=0.*/\1/' \
             | sort -u)
  for run in $finished; do
    if grep -qx "$run" "$DONE_LIST"; then continue; fi
    # GS save runs don't have a results/ tree the postprocess understands;
    # skip them.
    if [ ! -d "$run/results/raw" ]; then
      echo "[$(date '+%H:%M:%S')] $run: no results/raw, skipping" >>"$LOG"
      echo "$run" >>"$DONE_LIST"
      continue
    fi
    echo "[$(date '+%H:%M:%S')] $run: postprocess starting" >>"$LOG"
    if python3 scripts/coronene_postprocess.py run \
         --results "$run/results" --run-name "$run" \
         >>"$LOG" 2>&1
    then
      echo "[$(date '+%H:%M:%S')] $run: postprocess OK" >>"$LOG"
    else
      echo "[$(date '+%H:%M:%S')] $run: postprocess FAILED" >>"$LOG"
    fi
    echo "$run" >>"$DONE_LIST"
  done
  sleep 60
done
