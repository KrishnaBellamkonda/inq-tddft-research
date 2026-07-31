#!/bin/bash
#SBATCH --job-name=bulk-ks-rs4
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-bulk-ks-stopping-rs4.slurm — one half of the DENSITY-REPLICA twin pair.
#
#   Usage:  sbatch shared/bin/run-bulk-ks-stopping-rs4.slurm wp
#           sbatch shared/bin/run-bulk-ks-stopping-rs4.slurm classical
#
# 40 x 40 x 80 Bohr, N = 482 (r_s = 3.987), dx = 0.50, 646 steps at dt = 0.04.
# Projectile IDENTICAL to the r_s = 5.702 pair (sigma = 2, 100 eV, launch
# z = -32) — density is the only variable that moves.
#
# WALL TIME. Projected ~2.4 h (wp) / ~3.2 h (classical) from the 1.57x cost
# factor over the r_s = 5.702 pair. --time=12:00:00 gives ~4x headroom on the
# pessimistic O(n_states^2) projection, so the non-resumable classical half
# (INQ cannot restart moving ions; user decision 2026-07-30 to accept that)
# has room to overrun substantially without being killed.
#
# Both halves share the ground state produced by run-bulk-gs-rs4.slurm, so
# submit them with a dependency on it:
#
#   gs=$(sbatch --parsable shared/bin/run-bulk-gs-rs4.slurm)
#   sbatch --dependency=afterok:$gs shared/bin/run-bulk-ks-stopping-rs4.slurm wp
#   sbatch --dependency=afterok:$gs shared/bin/run-bulk-ks-stopping-rs4.slurm classical
#
# afterok (not afterany) is deliberate: the GS run.cpp exits 3 if any analytic
# gate fails, so a bad ground state stops the chain instead of silently seeding
# two production runs.
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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/jellium/scripts/bulk_ks_stopping_rs4/$HALF"
GS_DIR="$REPO_ROOT/ResearchProject/systems/jellium/checkpoints/gs_L40x40x80_orth_N482_dx0p50"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  bulk-ks-stopping-rs4 [$HALF]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " rundir : $RUN_DIR"
echo " gs     : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-bulk-gs-rs4.slurm must succeed first." >&2
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
