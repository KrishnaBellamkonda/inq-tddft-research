#!/usr/bin/env bash
#
# configure-inq.sh — CMake-configure an INQ engine tree. RUN ON THE LOGIN NODE.
#
# INQ's gpurun pulls Catch2 via CMake FetchContent (github.com). CSD3 compute
# nodes have no outbound network, so configure there fails with
#   "Failed to clone repository: 'https://github.com/catchorg/Catch2.git'".
# The login node does have network, so configure runs here (cheap, single
# process) and the compile then runs on an ampere node via build-inq.slurm.
#
# Usage (from the repo root):
#   shared/bin/configure-inq.sh            # configure inq/
#   shared/bin/configure-inq.sh inq-study  # configure the CAP engine
#   CLEAN=1 shared/bin/configure-inq.sh    # wipe and reconfigure from scratch
#
# Then:  sbatch shared/bin/build-inq.slurm
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INQ_TREE="${1:-inq}"
SRC="$REPO_ROOT/$INQ_TREE"
BUILD="$SRC/build"

# Same toolchain the compile and the runs use, so the cache is consistent.
# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

NVCC="$INQ_NVCC"

if [ "${CLEAN:-0}" = "1" ] && [ -d "$BUILD" ]; then
  echo "==> CLEAN=1: removing $BUILD ($(du -sh "$BUILD" | cut -f1))"
  rm -rf "$BUILD"
fi

echo "==> configuring $SRC  (nvcc=$NVCC, sm_80)"
cmake -S "$SRC" -B "$BUILD" \
  --install-prefix="$SRC/install" \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_CUDA_COMPILER="$NVCC"

echo
echo "Configured. Now compile on a GPU node:"
echo "  sbatch shared/bin/build-inq.slurm"
