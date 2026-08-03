#!/bin/bash
#SBATCH --job-name=ladder-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=16:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-ladder-gs.slurm — ground states for the CYLINDRICAL PROXIMITY LADDER.
# Plan: docs/plans/cylindrical-proximity-ladder.md
#
# One bare r_s = 3.000000 tube GS per rung. Geometry, N_e and extra_states come
# from shared/configs/proximity_ladder_rs3.hpp; the binary selects its rung at
# RUNTIME from CJ_RUNG and refuses to default.
#
#   rung   R_in   R_out    N_e   states   shape
#   r08     8.0   14.000   220     143    annulus
#   r06     6.0   13.986   266     173    annulus
#   r04     4.0   14.000   300     195    annulus
#   r00     0.0   13.986   326     212    CYLINDER (filled)
#   r04n160 4.0   10.583   160     104    annulus  (same-N control, opt-in)
#
# r10 (the completed channeling twin) is NOT here — its GS already exists at
# shared_gs/tube_rs3_R10_14_L60_dx0p5 and must not be recomputed.
#
# SEQUENTIAL, NOT AN ARRAY. inq-run builds in the run directory, so concurrent
# array tasks sharing scripts/proximity_ladder/gs would race on the same build
# tree. One job, one build, then a loop. Each rung is idempotent (skips a
# populated checkpoint), so a re-submit after a timeout resumes the chain.
#
#   Usage:  sbatch shared/bin/run-ladder-gs.slurm
#           CJ_RUNGS="r08 r06" sbatch shared/bin/run-ladder-gs.slurm
set -uo pipefail   # NOT -e: one rung failing its gates must not abandon the rest

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2
  echo "       Submit with: cd <repo> && sbatch shared/bin/run-ladder-gs.slurm" >&2
  exit 1
fi

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

SYS="$REPO_ROOT/ResearchProject/systems/cylindrical_jellium"
RUN_DIR="$SYS/scripts/proximity_ladder/gs"
RUNGS="${CJ_RUNGS:-r08 r06 r04 r00}"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  ladder-gs"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " rungs  : $RUNGS"
echo " started: $(date)"
echo "=============================================="

# DISK: WARN, DO NOT BLOCK (.claude/rules/checkpoint-dont-block.md — cost budgets
# belong to the user, correctness gates to the orchestrator; user confirmed
# 2026-08-02 that space will be available). Ground states scale with the STATE
# COUNT: rung 1 (104 states) saved 1.2 GB, so the 143/173/195/212-state rungs cost
# roughly 1.7/2.0/2.3/2.5 GB — about 8.5 GB for all four. The RT twins that follow
# are far larger (~10-16 GB per rung).
#
# The only HARD stop is a filesystem so full that electrons.save() would truncate
# a checkpoint, which corrupts the run rather than merely slowing it.
WANT_GB="${CJ_WANT_GB:-15}"
HARD_GB="${CJ_HARD_GB:-2}"
avail=$(df -BG --output=avail "$REPO_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
echo " disk available: ${avail:-unknown} GB (ladder ground states want ~${WANT_GB} GB)"
if [ -n "${avail:-}" ] && [ "$avail" -lt "$WANT_GB" ]; then
  echo " [WARN] below the ${WANT_GB} GB comfort mark — proceeding anyway."
  echo "        A rung that cannot save fails its own gates and the loop continues."
fi
if [ -n "${avail:-}" ] && [ "$avail" -lt "$HARD_GB" ]; then
  echo "FATAL: only ${avail} GB free — a checkpoint would be truncated, not merely tight." >&2
  exit 3
fi

cd "$RUN_DIR" || exit 1

fails=0
for RUNG in $RUNGS; do
  GS_DIR="$SYS/shared_gs/ladder_rs3_${RUNG}_L60_dx0p5"
  echo
  echo "########## RUNG $RUNG ##########"
  echo " gs dir : $GS_DIR"

  if [ -d "$GS_DIR" ] && [ -n "$(ls -A "$GS_DIR" 2>/dev/null)" ]; then
    echo " [skip] ground state already present — reusing it (idempotent chain)."
    continue
  fi

  # Re-check before EACH rung, but only for the corrupting case (see above).
  avail=$(df -BG --output=avail "$REPO_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
  if [ -n "${avail:-}" ] && [ "$avail" -lt "$HARD_GB" ]; then
    echo "FATAL: only ${avail} GB free — stopping before rung $RUNG rather than "
    echo "       writing a truncated checkpoint." >&2
    fails=$((fails+1)); break
  fi

  t0=$SECONDS
  export CJ_RUNG="$RUNG" GS_DIR GS_TAG="gs_${RUNG}" GS_SPACING="${GS_SPACING:-0.5}"
  inq-run
  rc=$?
  echo " rung $RUNG exit=$rc elapsed=$(( SECONDS - t0 ))s"

  S="results/gs_${RUNG}/run_summary.txt"
  if [ -f "$S" ]; then
    echo " --- run_summary.txt ($RUNG) ---"
    cat "$S"
  fi
  [ "$rc" -ne 0 ] && fails=$((fails+1))
done

echo
echo "=============================================="
echo " rungs with a non-zero exit : $fails"
echo " finished    : $(date)"
echo "=============================================="
exit $fails
