#!/bin/bash
#SBATCH --job-name=s56-fin
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=36:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-s56-finalize.slurm — the autonomous tail of the sigma56_sv campaign:
# status -> repair -> build deliverables -> report. Replaces the one-shot
# notebook stage, which would have built a figure out of whatever happened to be
# on disk.
#
#   Usage:  sbatch shared/bin/run-s56-finalize.slurm            # attempt 1 of 2
#           sbatch shared/bin/run-s56-finalize.slurm 2 2        # final attempt
#           sbatch shared/bin/run-s56-finalize.slurm 1 2 status # inspect only
#
# WHAT IT REPAIRS. Any expected run that is MISSING or SHORT of its step target
# is (re)run in place by invoking its OWN dispatcher as plain bash -- no sbatch
# from inside a job, and no duplicated launch logic. Runs with a checkpoint are
# resumed (LJ_RESUME=1); the rest start fresh. That covers the realistic failure
# modes: a walltime kill, a preempted node, a dependency that never fired.
#
# WHY A GPU AND 36 h. The repairs ARE production propagation, so this job needs
# the same resources a sweep task does. The budget is split: repairs get at most
# S56_REPAIR_BUDGET_S (default 30 h) and the rest is reserved for building the
# figures and notebooks, so the job ALWAYS reaches the reporting stage. A
# finalizer that spends its whole allocation waiting is not autonomy -- that
# exact failure cost 9.5 hours on 2026-06-28.
#
# TWO OF THESE ARE CHAINED with afterany. If the first is killed at walltime
# mid-repair, the second resumes from the checkpoints it left behind. Both always
# write ResearchProject/.../hypotheses/sigma56_sv/CAMPAIGN_REPORT.md.
#
# EMAIL is attempted and its failure is logged, never raised: there are no Gmail
# credentials on this machine yet (`python -m inqview.email setup` adds them), and
# a missing credential must not turn a finished campaign into a failed job.
#
# Plan: docs/plans/sigma56-sv-twin.md
set -euo pipefail

ATTEMPT="${1:-1}"
MAXATT="${2:-2}"
MODE="${3:-run}"

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# The repairs invoke the run dispatchers, which do their own toolchain setup.
# This job only needs the venv for the analysis, plus the paths the dispatchers
# inherit.
export PATH="$REPO_ROOT/shared/bin:$PATH"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

HYP="$REPO_ROOT/ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  s56-finalize [attempt $ATTEMPT/$MAXATT, $MODE]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " started: $(date)"
echo "=============================================="

# shellcheck source=/dev/null
. "$REPO_ROOT/venv/bin/activate"
cd "$HYP"

if [ "$MODE" = "status" ]; then
  python finalize.py --attempt "$ATTEMPT" --max "$MAXATT" --status-only
  exit $?
fi

# finalize.py returns 0 whenever it managed to write a report -- a non-zero code
# would only make the NEXT chained attempt look like it should be skipped, and an
# incomplete campaign is exactly when the next attempt matters.
python finalize.py --attempt "$ATTEMPT" --max "$MAXATT"
status=$?

echo ""
echo "=============================================="
echo " finished : $(date)"
echo " exit     : $status"
if [ -f "$HYP/CAMPAIGN_REPORT.md" ]; then
  echo " --- CAMPAIGN_REPORT.md ---"
  cat "$HYP/CAMPAIGN_REPORT.md"
fi
echo "=============================================="
exit $status
