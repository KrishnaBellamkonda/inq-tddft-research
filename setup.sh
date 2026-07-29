#!/usr/bin/env bash
#
# setup.sh — bootstrap the INQ engine + inq-stack on a fresh device.
#
# This repository carries ONLY source: the inqkit/inqview library
# (inq-stack/), the research run machinery (ResearchProject/), and the
# engine deltas needed to build. The heavy engine trees are NOT vendored:
#
#   * inq/       — pristine upstream INQ, pinned to commit ${INQ_COMMIT},
#                  cloned from GitLab and patched with inq-local.patch
#                  (the CUB fix for CUDA 12.5+ and the read-only ham()
#                  accessor used by inqkit KS-energy observables).
#   * inq-study/ — the project-modified engine (complex absorbing potential
#                  support), pulled in as a git submodule. REQUIRED for any
#                  CAP run; stock upstream inq cannot compile one.
#
# Ground states and run outputs are NOT shipped (see .gitignore). Regenerate
# ground states on this device by running the save_gs/*/run.cpp builders.
#
# Edit the CONFIG block for your machine, then: bash setup.sh
set -euo pipefail

# ----------------------------- CONFIG --------------------------------------
INQ_UPSTREAM="https://gitlab.com/npneq/inq.git"
INQ_COMMIT="44f73d9527ab677f38ed2138c2e83a28a5ab6c79"   # pinned upstream
CUDA_ARCH="${CUDA_ARCH:-80}"                            # sm_80 (A100); adjust
NVCC="${NVCC:-$(command -v nvcc || echo /usr/local/cuda/bin/nvcc)}"
INSTALL_PREFIX="${INSTALL_PREFIX:-$HOME}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# ---------------------------------------------------------------------------

echo "==> repo root: $ROOT"

# 1. Modified engine (inq-study) — git submodule
echo "==> [1/4] fetching inq-study submodule (project-modified engine, CAP)"
git -C "$ROOT" submodule update --init --recursive inq-study

# 2. Pristine upstream inq at the pinned commit + local patch
if [ ! -d "$ROOT/inq/.git" ]; then
  echo "==> [2/4] cloning upstream inq @ ${INQ_COMMIT:0:12}"
  git clone "$INQ_UPSTREAM" "$ROOT/inq"
  git -C "$ROOT/inq" checkout "$INQ_COMMIT"
  echo "    applying inq-local.patch (CUB fix + ham() accessor)"
  git -C "$ROOT/inq" apply "$ROOT/inq-local.patch"
else
  echo "==> [2/4] inq/ already present — skipping clone (verify patch manually)"
fi

# 3. Build inq
echo "==> [3/4] building inq (CUDA arch sm_${CUDA_ARCH}, nvcc=$NVCC)"
cmake -S "$ROOT/inq" -B "$ROOT/inq/build" \
  --install-prefix="$INSTALL_PREFIX" \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCH" \
  -DCMAKE_CUDA_COMPILER="$NVCC"
cmake --build "$ROOT/inq/build" --parallel
cmake --install "$ROOT/inq/build"

# 4. inqview (Python post-processing) into the active environment
echo "==> [4/4] installing inqview (pip install -e inq-stack/)"
python3 -m pip install -e "$ROOT/inq-stack/"

cat <<EOF

Done. Add to your shell profile so 'inq-run' works:
  export PATH="$ROOT/shared/bin:\$PATH"
  export INQ_SHARE_PATH=$ROOT/inq/install/share
  export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod

To build inq-study (CAP runs) point your build at inq-study/ instead of inq/.
Regenerate ground states with the ResearchProject/systems/**/save_gs/*/run.cpp
builders before launching runs that load a saved GS.
EOF
