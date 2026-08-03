#!/bin/bash
#SBATCH --job-name=ladder-buildcheck
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-ladder-buildcheck.slurm — compile the proximity-ladder RT binaries and
# verify their rung-selection guard, WITHOUT running any physics.
# Plan: docs/plans/cylindrical-proximity-ladder.md
#
# WHY THIS EXISTS. The ground-state chain runs for hours; discovering a compile
# error in the RT halves only after it finishes wastes a night. inq-run has no
# build-only mode, but both RT binaries validate CJ_RUNG and CH_GS_DIR in the
# first few statements of main() and return 2 — so "build, then exit in
# milliseconds" is available for free, and it exercises the guards as a bonus.
#
# Three cases per binary:
#   1. CJ_RUNG unset      -> exit 2, lists the valid rungs
#   2. CJ_RUNG=nonsense   -> exit 2, names the bad value
#   3. CJ_RUNG=r08, no GS -> exit 2 at the GS guard (so the rung table, the
#                            r_out() solve and the density self-check all ran)
#
# A binary that exits 0 here would mean a guard is missing, which is a FAILURE.
#
#   Usage:  cd <repo> && sbatch shared/bin/run-ladder-buildcheck.slurm
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2
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

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  ladder-buildcheck"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>&1 | head -1)"
echo " started: $(date)"
echo "=============================================="

fails=0
for HALF in wp classical; do
  echo
  echo "########## $HALF ##########"
  cd "$SYS/scripts/proximity_ladder/$HALF" || { fails=$((fails+1)); continue; }

  # --- compile (the run exits at the first guard, so this is a build check) ---
  t0=$SECONDS
  env -u CJ_RUNG -u CH_GS_DIR inq-run > build_and_guard.log 2>&1
  rc=$?
  echo " build+guard1 exit=$rc  elapsed=$(( SECONDS - t0 ))s"

  if grep -qE "error:|Error [0-9]|make.*\*\*\*" build_and_guard.log; then
    echo " [FAIL] COMPILE ERROR — first errors:"
    grep -nE "error:|Error [0-9]|make.*\*\*\*" build_and_guard.log | head -15
    fails=$((fails+1)); continue
  fi
  echo " [PASS] compiled."

  if [ "$rc" -eq 2 ] && grep -q "CJ_RUNG is unset" build_and_guard.log; then
    echo " [PASS] guard 1: unset CJ_RUNG refuses (exit 2)."
  else
    echo " [FAIL] guard 1: expected exit 2 + 'CJ_RUNG is unset', got exit $rc"
    tail -5 build_and_guard.log; fails=$((fails+1))
  fi

  # --- guard 2: an unknown rung must be named, not defaulted -----------------
  # NOTE the flag order: GNU env stops parsing options at the first NAME=VALUE,
  # so `env CJ_RUNG=x -u CH_GS_DIR ./run` tries to EXECUTE `-u` and exits 127.
  # Every -u must precede every assignment.
  env -u CH_GS_DIR CJ_RUNG=notarung ./run > guard2.log 2>&1
  rc=$?
  if [ "$rc" -eq 2 ] && grep -q "notarung" guard2.log; then
    echo " [PASS] guard 2: unknown rung refuses and names it (exit 2)."
  else
    echo " [FAIL] guard 2: expected exit 2 naming 'notarung', got exit $rc"
    tail -5 guard2.log; fails=$((fails+1))
  fi

  # --- guard 3: a VALID rung gets past the table self-check ------------------
  # Reaching the CH_GS_DIR guard proves rung lookup, r_out() and
  # max_density_error() < 1e-12 all executed successfully.
  env -u CH_GS_DIR CJ_RUNG=r08 ./run > guard3.log 2>&1
  rc=$?
  if [ "$rc" -eq 2 ] && grep -q "CH_GS_DIR" guard3.log; then
    echo " [PASS] guard 3: r08 passes the rung table, stops at the GS guard."
  else
    echo " [FAIL] guard 3: expected exit 2 at CH_GS_DIR, got exit $rc"
    tail -8 guard3.log; fails=$((fails+1))
  fi
done

echo
echo "=============================================="
echo " failures : $fails"
echo " finished : $(date)"
echo "=============================================="
exit $fails
