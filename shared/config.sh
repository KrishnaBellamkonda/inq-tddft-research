# INQ shared configuration
# Edit this file if INQ source, CUDA version, or Python path changes.
#
# Paths are derived from the repository root (this file's location) so the same
# checkout works on any device. Every value can still be overridden from the
# environment. Migrated from absolute /local/data/public paths on 2026-07-29
# when the repo moved to CSD3 — see docs/plans/csd3-setup-cuda121-build.md.

_INQ_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# INQ_SOURCE may be overridden in the environment to build against a fork
# (e.g. inq-study). Share paths default off INQ_SOURCE but can be pinned
# separately (a fork without its own install/share can reuse inq's).
export INQ_SOURCE="${INQ_SOURCE:-$_INQ_REPO_ROOT/inq}"
export INQ_SHARE_PATH="${INQ_SHARE_PATH:-$INQ_SOURCE/install/share}"
export PSEUDOPOD_SHARE_PATH="${PSEUDOPOD_SHARE_PATH:-$INQ_SOURCE/install/share/pseudopod}"

# CUDA. CSD3's newest toolkit is 12.1 (nothing >= 12.4 exists on the cluster);
# the gpu/reduce.hpp CUB fix in inq-local.patch is written to work on 12.1
# through 12.6+. sm_80 = A100 (ampere partition).
export INQ_CUDA_COMPILER="${INQ_CUDA_COMPILER:-$(command -v nvcc || echo /usr/local/software/cuda/12.1/bin/nvcc)}"
export INQ_CUDA_ARCH="${INQ_CUDA_ARCH:-80}"

# Python used for INQ's cmake Python detection and for inqview post-processing.
export INQ_PYTHON_EXE="${INQ_PYTHON_EXE:-$_INQ_REPO_ROOT/venv/bin/python3}"

# FetchContent source cache. INQ pulls catch2/pybind11/spdlog/spglib from GitHub at
# configure time, but CSD3 COMPUTE NODES HAVE NO OUTBOUND NETWORK, so a fresh
# configure there dies with "Failed to clone repository: .../Catch2.git". The main
# engine configure (run on the login node) already populated these, so inq-run
# points FetchContent at them instead of cloning. See shared/bin/configure-inq.sh.
export INQ_DEPS_CACHE="${INQ_DEPS_CACHE:-$INQ_SOURCE/build/_deps}"

# Parallel compile jobs. Deliberately NOT nproc: CSD3 login nodes cap a user at
# 20 GB and each nvcc `cicc` on INQ's template-heavy TUs takes several GB, which
# is what OOM-killed the first setup run (16x "cicc died due to signal 9").
# On a compute node with a real memory allocation, raise this via the env.
export INQ_BUILD_JOBS="${INQ_BUILD_JOBS:-4}"
