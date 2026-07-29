#!/usr/bin/env bash
# Sequential S(v) ladder dispatcher for ONE GPU.
# Usage: dispatch_ladder.sh <gpu_id> "v0:nsteps:subdir[:psp] v0:nsteps:subdir ..."
set -u
GPU="$1"; shift
QUEUE="$*"
RUN_DIR="/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_sv_sigma0p5"
cd "$RUN_DIR"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
LOG="$RUN_DIR/dispatch_gpu${GPU}.log"
echo "=== GPU $GPU dispatcher start $(date -u +%FT%TZ) queue: $QUEUE ===" >> "$LOG"
for item in $QUEUE; do
  IFS=':' read -r V0 NSTEPS SUB PSP <<< "$item"
  echo "--- [$(date -u +%FT%TZ)] GPU$GPU v0=$V0 nsteps=$NSTEPS sub=$SUB psp=${PSP:-default} ---" >> "$LOG"
  ENVV=(CUDA_VISIBLE_DEVICES="$GPU" PROJ_V0="$V0" SV_N_STEPS="$NSTEPS" SV_WRITE_EVERY=50 SV_OUT_SUBDIR="$SUB")
  if [ -n "${PSP:-}" ]; then
    ENVV+=(SV_PSEUDO="/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/$PSP")
  fi
  env "${ENVV[@]}" ./run >> "$LOG" 2>&1
  echo "--- [$(date -u +%FT%TZ)] GPU$GPU v0=$V0 DONE rc=$? ---" >> "$LOG"
done
echo "=== GPU $GPU dispatcher finished $(date -u +%FT%TZ) ===" >> "$LOG"
