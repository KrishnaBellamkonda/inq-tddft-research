#!/bin/bash
#SBATCH --job-name=lzb-fin
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
# run-lzb-finalize.slurm — the autonomous tail of lz_bulk_sweep:
# status -> repair -> build deliverables -> report.
# Plan: docs/plans/jellium-slab-extend-Lz.md
#
#   Usage:  sbatch shared/bin/run-lzb-finalize.slurm            # attempt 1 of 2
#           sbatch shared/bin/run-lzb-finalize.slurm 2 2        # final attempt
#           sbatch shared/bin/run-lzb-finalize.slurm 1 2 status # inspect only
#
# Repairs re-run MISSING/SHORT runs in place through their own dispatchers as
# plain bash (no sbatch from inside a job), each capped at ~2x its measured cost
# so one wedged run cannot swallow the allocation (the sigma56 lesson, already
# fixed in this clone's finalize.py). Two of these are chained afterany; both
# always write hypotheses/lz_bulk_sweep/CAMPAIGN_REPORT.md.
set -euo pipefail

ATTEMPT="${1:-1}"
MAXATT="${2:-2}"
MODE="${3:-run}"

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

HYP="$REPO_ROOT/ResearchProject/systems/localised_jellium/hypotheses/lz_bulk_sweep"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  lzb-finalize [attempt $ATTEMPT/$MAXATT, $MODE]"
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
