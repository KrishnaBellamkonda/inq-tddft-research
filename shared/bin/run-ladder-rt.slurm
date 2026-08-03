#!/bin/bash
#SBATCH --job-name=ladder-rt
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-ladder-rt.slurm — ONE half of ONE rung of the cylindrical proximity ladder.
# Plan: docs/plans/cylindrical-proximity-ladder.md
#
#   Usage:  sbatch shared/bin/run-ladder-rt.slurm <wp|classical> <rung> [smoke|prod]
#           e.g.  sbatch shared/bin/run-ladder-rt.slurm wp r08 prod
#
# THE BINARY IS RUNG-INDEPENDENT. Geometry is selected at RUNTIME from CJ_RUNG, so
# ONE build per half serves every rung — which is why the build lives in
# run-ladder-buildcheck.slurm and this script refuses to build. That also removes
# the only reason rungs would have to run sequentially: with ./run already present,
# every (half, rung) pair is independent and the chain can fan out across GPUs.
#
# EXPECTED PHYSICS, so a reader does not mistake it for a fault:
#   * The projectile DECELERATES — it is a mass-1 electron carrying 50 eV
#     (.claude/rules/light-projectile-stopping.md). S is the INITIAL DRAG over the
#     early v >= 0.85 v0 window, never a full-run regression. Estimated fractional
#     energy loss over the run climbs from ~13 % (r10) to ~66 % (r00).
#   * The strong rungs are SUPPOSED to look violent. r04 starts with 37 % of the
#     packet already inside the wall and r00 is immersed from step 0.
#
# Resumable per .claude/rules/final-timestep-checkpoint.md: re-submit with
# CH_RESUME=1 and a larger CH_N_STEPS to EXTEND rather than recompute.
set -uo pipefail

HALF="${1:-}"
RUNG="${2:-}"
STAGE="${3:-prod}"
case "$HALF" in wp|classical) ;;
  *) echo "FATAL: arg 1 must be 'wp' or 'classical'; got '$HALF'" >&2; exit 2 ;; esac
case "$RUNG" in r08|r06|r04|r00|r10|r04n160) ;;
  *) echo "FATAL: arg 2 must be a rung (r08 r06 r04 r00 r04n160); got '$RUNG'" >&2; exit 2 ;; esac
case "$STAGE" in smoke|prod) ;;
  *) echo "FATAL: arg 3 must be 'smoke' or 'prod'; got '$STAGE'" >&2; exit 2 ;; esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 1
fi
# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

SYS="$REPO_ROOT/ResearchProject/systems/cylindrical_jellium"
RUN_DIR="$SYS/scripts/proximity_ladder/$HALF"
GS_DIR="$SYS/shared_gs/ladder_rs3_${RUNG}_L60_dx0p5"

if [ "$STAGE" = "smoke" ]; then
  NAME="smoke_${RUNG}"; NSTEPS=20; SAVE=10; WF=0; CKPT=999999
else
  NAME="${RUNG}"
  NSTEPS="${CH_N_STEPS:-1500}"; SAVE="${CH_SAVE_EVERY:-15}"
  WF="${CH_WF_EVERY:-150}";     CKPT="${CH_CKPT_EVERY:-500}"
fi
export CH_MAX_CKPT="${CH_MAX_CKPT:-1}"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  ladder-rt [$HALF/$RUNG/$STAGE]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " rung   : $RUNG   out: $NAME"
echo " steps  : $NSTEPS x dt 0.02 = $(awk "BEGIN{printf \"%.1f\", $NSTEPS*0.02}") a.u."
echo " cadence: density/$SAVE  wavefn/$WF  ckpt/$CKPT  scalars/1"
echo " resume : ${CH_RESUME:-0}"
echo " gs     : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

# --- correctness preconditions (these DO block; cost never does) -------------
if [ ! -d "$GS_DIR" ] || [ -z "$(ls -A "$GS_DIR" 2>/dev/null)" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-ladder-gs.slurm must succeed first." >&2
  exit 2
fi
cd "$RUN_DIR" || exit 1
if [ ! -x ./run ]; then
  echo "FATAL: ./run not built in $RUN_DIR." >&2
  echo "       run-ladder-buildcheck.slurm builds both halves; it must succeed first." >&2
  exit 2
fi

# Disk WARNS, never blocks (.claude/rules/checkpoint-dont-block.md). A killed run
# resumes from its last checkpoint; a run never launched buys nothing.
avail=$(df -BG --output=avail "$REPO_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
echo " disk available: ${avail:-unknown} GB"
if [ -n "${avail:-}" ] && [ "$avail" -lt 10 ] && [ "$STAGE" = "prod" ]; then
  echo " [WARN] under 10 GB free; this half wants ~5-8 GB. Proceeding — a truncated"
  echo "        run is recoverable via CH_RESUME=1, a run never started is not."
fi

export CJ_RUNG="$RUNG"
export CH_OUT="$NAME" CH_N_STEPS="$NSTEPS" CH_SAVE_EVERY="$SAVE" \
       CH_WF_EVERY="$WF" CH_CKPT_EVERY="$CKPT" \
       CH_GS_DIR="$GS_DIR" CH_RESUME="${CH_RESUME:-0}"

t0=$SECONDS
env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run
status=$?

echo "=============================================="
echo " exit status : $status   elapsed: $(( SECONDS - t0 ))s"
echo " finished    : $(date)"
S="results/$NAME/run_summary.txt"
[ -f "$S" ] && { echo " --- run_summary.txt ---"; cat "$S"; }
echo "=============================================="
exit $status
