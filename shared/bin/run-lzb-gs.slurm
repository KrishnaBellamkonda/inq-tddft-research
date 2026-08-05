#!/bin/bash
#SBATCH --job-name=lzb-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-lzb-gs.slurm — one ground state of the lz_bulk_sweep campaign.
# Plan: docs/plans/jellium-slab-extend-Lz.md
#
#   Usage:  sbatch shared/bin/run-lzb-gs.slurm s0p5_L15         # builds + runs
#           sbatch shared/bin/run-lzb-gs.slurm s5p0_L35 exec    # reuse ./run
#
# FOUR ground states, one per box preset (the sweep grows the slab at fixed n0,
# so every thickness is a new electronic system: N = 60/140, states 45/104).
# The chain submits the FIRST as the builder (inq-run) and the other three as
# `exec` afterok on it, so concurrent jobs never race on one CMake tree — the
# same build-once pattern as the sigma56 sweep arrays.
#
# Gates live in the binary: electron count and r_s are HARD; E_GS is recorded as
# the per-box deposit reference (never compared across boxes). Bulk-likeness of
# the interior density is checked later by pilot_gate.py from the density VTI
# this run writes (WARN, not a gate — checkpoint-dont-block).
set -euo pipefail

CFG="${1:?pass a box preset: s0p5_L15 | s0p5_L35 | s5p0_L15 | s5p0_L35}"
MODE="${2:-build}"

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

export LZB_CFG="$CFG"
# shellcheck source=./lzb-params.sh
. "$REPO_ROOT/shared/bin/lzb-params.sh"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/lz_bulk_sweep/gs"
GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/$GS_TAG"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  lzb-gs [$CFG, $MODE]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " box    : 35 x 35 x $LZ Bohr, slab $(awk "BEGIN{print 2*$HALF}") Bohr thick, N=$N_E, dx=0.40"
echo " gs dir : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"
export GS_DIR

if [ "$MODE" = "exec" ]; then
  if [ ! -x ./run ]; then
    echo "FATAL: ./run not built. Submit the 'build' mode first." >&2
    exit 2
  fi
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run
  status=$?
else
  inq-run
  status=$?
fi

echo "=============================================="
echo " [$CFG] exit status : $status"
echo " finished           : $(date)"
[ -f "results/$CFG/run_summary.txt" ] && { echo " --- run_summary.txt ---"; cat "results/$CFG/run_summary.txt"; }
echo "=============================================="
exit $status
