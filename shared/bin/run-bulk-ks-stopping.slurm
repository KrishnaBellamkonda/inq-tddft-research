#!/bin/bash
#SBATCH --job-name=bulk-ks
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-bulk-ks-stopping.slurm — one half of the bulk-jellium KS-stopping twin pair.
#
#   Usage:  sbatch shared/bin/run-bulk-ks-stopping.slurm wp
#           sbatch shared/bin/run-bulk-ks-stopping.slurm classical
#
# Both halves share the ground state produced by run-bulk-gs.slurm, so submit
# them with a dependency on it:
#
#   gs=$(sbatch --parsable shared/bin/run-bulk-gs.slurm)
#   sbatch --dependency=afterok:$gs shared/bin/run-bulk-ks-stopping.slurm wp
#   sbatch --dependency=afterok:$gs shared/bin/run-bulk-ks-stopping.slurm classical
#
# afterok (not afterany) is deliberate: the GS run.cpp exits 3 if any of its
# analytic gates fail (E_kinetic = 7.395907 Ha, E_hartree = 0, gap 0.244 eV,
# occupations 2.0/0.0), so a bad ground state stops the chain instead of
# silently seeding two 24-hour production runs.
#
# STOCK inq — this study constructs no absorbing potential anywhere, so the
# complexified inq-study fork is not needed (contrast run-dispersion.slurm).
#
# Plan: docs/plans/bulk-jellium-ks-stopping.md
set -euo pipefail

HALF="${1:-}"
if [ "$HALF" != "wp" ] && [ "$HALF" != "classical" ]; then
  echo "FATAL: first argument must be 'wp' or 'classical' (got '${HALF}')" >&2
  exit 2
fi

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/jellium/scripts/bulk_ks_stopping/$HALF"
GS_DIR="$REPO_ROOT/ResearchProject/systems/jellium/checkpoints/gs_L46x46x80_orth_N218_dx0p40"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  bulk-ks-stopping [$HALF]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " rundir : $RUN_DIR"
echo " gs     : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-bulk-gs.slurm must succeed first." >&2
  exit 2
fi

cd "$RUN_DIR"
inq-run
status=$?

echo "=============================================="
echo " [$HALF] exit status : $status"
echo " finished            : $(date)"
if [ -f results/run_summary.txt ]; then
  echo " --- run_summary.txt ---"
  cat results/run_summary.txt
fi
for f in results/raw/observables/electron_track.csv \
         results/raw/observables/wp_momentum_stats.csv; do
  if [ -f "$f" ]; then
    echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"
    head -2 "$f"
    tail -1 "$f"
  fi
done
echo "=============================================="
exit $status
