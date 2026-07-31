#!/bin/bash
#SBATCH --job-name=bulk-sig
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
# run-bulk-sigma-sweep.slurm — one run of the sigma sweep.
#
#   Usage:  sbatch shared/bin/run-bulk-sigma-sweep.slurm <variant> <half>
#     variant : sigma1 | sigma3 | rs4_sigma1 | rs4_sigma3
#     half    : wp | classical
#
# Tests the WIDTH axis of the S_classical / S_WP gap, which the density sweep
# left open (2.92x density moved the ratio only 6.49x -> 5.65x).
#
#   sigma1 / sigma3      : r_s = 5.702 bath (46x46x80, N=218, dx=0.40)
#   rs4_sigma1 / _sigma3 : r_s = 3.987 bath (40x40x80, N=482, dx=0.50)
#
# NO NEW GROUND STATE IS NEEDED: the GS depends only on box + density, not on
# the projectile, so these reuse the two existing checkpoints. That is why the
# sweep has no afterok chain in front of it and all eight runs start at once.
#
# Launch position follows the BOUNDARY RULE z0 = -L_z/2 + 4 sigma (user
# decision, 2026-07-31): -36 for sigma=1, -32 for sigma=2, -28 for sigma=3.
# Path lengths therefore DIFFER per sigma; cross-sigma comparison must use the
# common TIME window [4, 9.37] a.u., not raw path. See the config headers.
#
# Plan: docs/plans/bulk-jellium-ks-stopping.md
set -euo pipefail

VARIANT="${1:-}"
HALF="${2:-}"
case "$VARIANT" in
  sigma1|sigma3)          GS_NAME="gs_L46x46x80_orth_N218_dx0p40" ;;
  rs4_sigma1|rs4_sigma3)  GS_NAME="gs_L40x40x80_orth_N482_dx0p50" ;;
  *) echo "FATAL: variant must be sigma1|sigma3|rs4_sigma1|rs4_sigma3 (got '$VARIANT')" >&2; exit 2 ;;
esac
if [ "$HALF" != "wp" ] && [ "$HALF" != "classical" ]; then
  echo "FATAL: second argument must be 'wp' or 'classical' (got '$HALF')" >&2; exit 2
fi

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/jellium/scripts/bulk_ks_stopping_${VARIANT}/$HALF"
GS_DIR="$REPO_ROOT/ResearchProject/systems/jellium/checkpoints/$GS_NAME"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  sigma-sweep [$VARIANT/$HALF]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " rundir : $RUN_DIR"
echo " gs     : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

if [ ! -d "$RUN_DIR" ]; then
  echo "FATAL: run directory missing at $RUN_DIR" >&2; exit 2
fi
if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR" >&2; exit 2
fi

cd "$RUN_DIR"
inq-run
status=$?

echo "=============================================="
echo " [$VARIANT/$HALF] exit status : $status"
echo " finished                     : $(date)"
if [ -f results/run_summary.txt ]; then
  echo " --- run_summary.txt ---"
  cat results/run_summary.txt
fi
echo "=============================================="
exit $status
