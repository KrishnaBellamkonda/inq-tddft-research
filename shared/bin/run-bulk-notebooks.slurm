#!/bin/bash
#SBATCH --job-name=bulk-nb
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --array=0-11%4
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.out
#
# run-bulk-notebooks.slurm — rebuild every bulk-jellium twin-pair notebook now
# that all four pairs carry interactions.csv on BOTH halves.
#
#   Usage:  sbatch shared/bin/run-bulk-notebooks.slurm
#
# TWELVE NOTEBOOKS, two kinds:
#
#   task  kind   family                       half        product
#   ----  -----  ---------------------------  ----------  ---------------------------
#    0    run    bulk_ks_stopping             wp          <family>_wp.ipynb
#    1    run    bulk_ks_stopping             classical   <family>_classical.ipynb
#    2    run    bulk_ks_stopping_sigma3      wp
#    3    run    bulk_ks_stopping_sigma3      classical
#    4    run    bulk_ks_stopping_rs4         wp
#    5    run    bulk_ks_stopping_rs4         classical
#    6    run    bulk_ks_stopping_rs4_sigma3  wp
#    7    run    bulk_ks_stopping_rs4_sigma3  classical
#    8    phase  bulk_ks_stopping             (pair)      run_pair_phase_<family>.ipynb
#    9    phase  bulk_ks_stopping_sigma3      (pair)
#   10    phase  bulk_ks_stopping_rs4         (pair)
#   11    phase  bulk_ks_stopping_rs4_sigma3  (pair)
#
# WHY THE RUN NOTEBOOKS ARE REBUILT: they were last built 2026-07-30/31, BEFORE
# the interaction-energy re-run, so they carry none of the pairwise decomposition.
# Their builders have since gained a section 7 (E_SS / E_PP / E_PS, the closure
# gate against INQ's own Hartree energy, and the projectile-cloud clipping
# diagnostic); the old takeaway moved to section 8.
#
# WHY THE PHASE NOTEBOOKS ARE NEW: the existing per-pair notebook compares the two
# halves in the ENERGY ledger. The phase notebook compares them in PHASE SPACE --
# the (z, v) portrait, where the trajectories part, and the local S(z) = -dT/dz
# along the path. Both halves are light free-Ehrenfest projectiles and DECELERATE
# (.claude/rules/light-projectile-stopping.md), so a single S per run averages
# over a whole velocity range; the portrait puts that range back on the page.
#
# NO GPU COMPUTE HAPPENS HERE. These are pure post-processing (VTI reads,
# matplotlib, GIF encoding). The GPU is requested only because this project's
# account is associated with the ampere partition alone -- see the account list
# in `sacctmgr show assoc user=skcb2`. Keep --cpus-per-task modest for the same
# reason: the work is single-process Python.
#
# Each task is INDEPENDENT: one failure does not stop the others, and a task can
# be resubmitted alone with --array=<n>.
#
# Plan: docs/plans/bulk-jellium-ks-stopping.md
set -uo pipefail

FAMILIES=(bulk_ks_stopping bulk_ks_stopping bulk_ks_stopping_sigma3 \
          bulk_ks_stopping_sigma3 bulk_ks_stopping_rs4 bulk_ks_stopping_rs4 \
          bulk_ks_stopping_rs4_sigma3 bulk_ks_stopping_rs4_sigma3)
HALVES=(wp classical wp classical wp classical wp classical)
PHASE_FAMILIES=(bulk_ks_stopping bulk_ks_stopping_sigma3 \
                bulk_ks_stopping_rs4 bulk_ks_stopping_rs4_sigma3)

I="${SLURM_ARRAY_TASK_ID:-}"
if [ -z "$I" ] || [ "$I" -lt 0 ] || [ "$I" -gt 11 ]; then
  echo "FATAL: SLURM_ARRAY_TASK_ID must be 0..11 (got '${I:-unset}')" >&2; exit 2
fi

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

HYP="$REPO_ROOT/ResearchProject/systems/jellium/hypotheses"
PY="$REPO_ROOT/venv/bin/python"
[ -x "$PY" ] || { echo "FATAL: $PY not found" >&2; exit 2; }

if [ "$I" -lt 8 ]; then
  KIND="run"; FAM="${FAMILIES[$I]}"; HALF="${HALVES[$I]}"
  WORKDIR="$HYP/$FAM"; SCRIPT="$WORKDIR/build_run_notebook.py"; ARGS=("$HALF")
else
  KIND="phase"; FAM="${PHASE_FAMILIES[$((I-8))]}"; HALF="(pair)"
  WORKDIR="$HYP/bulk_ks_stopping"
  SCRIPT="$HYP/bulk_ks_stopping/build_phase_notebook.py"; ARGS=("$FAM")
fi

echo "=============================================="
echo " job    : ${SLURM_ARRAY_JOB_ID:-none}_${I}  bulk-nb [$KIND: $FAM/$HALF]"
echo " node   : $(hostname)"
echo " script : $SCRIPT"
echo " started: $(date)"
echo "=============================================="

[ -f "$SCRIPT" ] || { echo "FATAL: $SCRIPT missing" >&2; exit 2; }

# The run-notebook builders gained the interaction-energy section; refuse to
# rebuild one that somehow lacks it, rather than silently reproducing the stale
# notebook this job exists to replace.
if [ "$KIND" = "run" ] && ! grep -q "Pairwise interaction energies" "$SCRIPT"; then
  echo "FATAL: $FAM/build_run_notebook.py has no interaction-energy section" >&2
  exit 2
fi

cd "$WORKDIR" || exit 2
"$PY" "$SCRIPT" "${ARGS[@]}"
status=$?

echo "=============================================="
echo " [$KIND: $FAM/$HALF] builder exit : $status"
echo " finished                          : $(date)"
if [ "$KIND" = "run" ]; then
  NB="$HYP/$FAM/${FAM}_${HALF}.ipynb"
else
  NB="$HYP/$FAM/run_pair_phase_${FAM}.ipynb"
fi
if [ -f "$NB" ]; then
  echo " notebook : $NB"
  echo " size     : $(du -h "$NB" | cut -f1)   modified: $(date -r "$NB")"
else
  echo " notebook : MISSING at $NB"
  status=1
fi
echo "=============================================="
exit $status
