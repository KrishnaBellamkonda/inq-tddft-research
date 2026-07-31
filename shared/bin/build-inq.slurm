#!/bin/bash
#SBATCH --job-name=inq-build
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
# build-inq.slurm — build the INQ engine on an A100 (ampere) compute node.
#
# Why a compute node rather than the login node: CSD3 caps a login user at 20 GB,
# and each nvcc `cicc` on INQ's template-heavy translation units takes several GB.
# An unbounded parallel build there is OOM-killed ("cicc died due to signal 9").
# The ampere node also has the A100 that ctest needs.
#
# Usage:
#   sbatch shared/bin/build-inq.slurm            # build inq/
#   INQ_TREE=inq-study sbatch shared/bin/build-inq.slurm   # build the CAP engine
#
# Env:
#   INQ_TREE    engine dir relative to repo root (default: inq)
#   BUILD_JOBS  parallel compile jobs (default: cpus-per-task)
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
INQ_TREE="${INQ_TREE:-inq}"
SRC="$REPO_ROOT/$INQ_TREE"
BUILD="$SRC/build"
JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}"

# One shared toolchain definition for configure, build and runs (also strips the
# stale CUDA 11.4 that rhel8/default-amp drags in — see csd3-env.sh).
# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

NVCC="$INQ_NVCC"

echo "=============================================="
echo " node      : $(hostname)"
echo " engine    : $SRC"
echo " nvcc      : $NVCC ($("$NVCC" --version | tail -1))"
echo " jobs      : $JOBS"
echo " gpu       : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>&1 | head -1)"
echo " started   : $(date)"
echo "=============================================="

# CONFIGURE MUST HAPPEN ON THE LOGIN NODE. INQ's gpurun uses CMake FetchContent to
# clone Catch2 from github.com, and CSD3 compute nodes have NO outbound network:
# configuring here dies with "Failed to clone repository: .../Catch2.git".
# Run configure-inq.sh on the login node first (it has network); this job then only
# compiles. CLEAN=1 wipes the build dir, which also discards the fetched deps — do
# that on the login node too, never here.
if [ ! -f "$BUILD/CMakeCache.txt" ]; then
  echo "ERROR: $BUILD is not configured." >&2
  echo "       Compute nodes have no network, so FetchContent cannot run here." >&2
  echo "       On the login node run:  shared/bin/configure-inq.sh $INQ_TREE" >&2
  exit 1
fi

cmake --build "$BUILD" --parallel "$JOBS"
cmake --install "$BUILD"

echo "=============================================="
echo " BUILD OK — finished $(date)"
echo " install prefix: $SRC/install"
echo "=============================================="
