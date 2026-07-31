#!/bin/bash
#SBATCH --job-name=slabks-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-slab-ks-gs.slurm — ground state for the CAP-free wrap-around slab
# KS-stopping study (docs/plans/slab-ks-orbital-stopping-wrap.md).
#
# 35 x 35 x 85 Bohr, periodicity(2), dx = 0.40, 25-Bohr slab (half-width 12.5,
# erfc edge 1.0), extra_states 24, T = 0.00862 eV, LDA.
#
#   Usage:  sbatch shared/bin/run-slab-ks-gs.slurm 40     # r_s = 5.67 (NEW)
#           sbatch shared/bin/run-slab-ks-gs.slurm 100    # r_s = 4.18 (already exists)
#
# N is the ONLY thing that changes between the two campaigns — the slab volume
# is compile-time, so n0 = N/30625 follows from N alone and the two systems
# cannot differ in anything else by construction.
#
#   N = 100  ->  n0 = 3.2653e-3, r_s = 4.18, 50 occupied + 24 = 74 states
#   N =  40  ->  n0 = 1.3061e-3, r_s = 5.67, 20 occupied + 24 = 44 states
#
# The N=100 dx=0.40 ground state was already produced by the WP high-density
# campaign (run-wp-hd-gs.slurm) and is REUSED — this script skips it if the
# checkpoint is already there, so re-running the chain is idempotent.
set -euo pipefail

N="${1:-40}"
case "$N" in
  40|100) ;;
  *) echo "FATAL: pass 40 or 100; got '$N'" >&2; exit 2 ;;
esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/slab_ks_wrap/gs"
GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n${N}_L35x35x85_dx0p4_per2"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  slabks-gs  N=$N"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " gs dir : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

if [ -d "$GS_DIR" ] && [ -n "$(ls -A "$GS_DIR" 2>/dev/null)" ]; then
  echo " [skip] ground state already present — reusing it (idempotent chain)."
  exit 0
fi

cd "$RUN_DIR"
export GS_DIR GS_SPACING="0.40" GS_TAG="n${N}_dx0p4" GS_N="$N" GS_PERIODICITY=2
inq-run
status=$?

echo "=============================================="
echo " exit status : $status"
echo " finished    : $(date)"
[ -f "results/n${N}_dx0p4/run_summary.txt" ] && { echo " --- run_summary.txt ---"; cat "results/n${N}_dx0p4/run_summary.txt"; }
echo "=============================================="
exit $status
