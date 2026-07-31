#!/bin/bash
#SBATCH --job-name=rebuild-pinq
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --output=rebuild-pinq-%j.out

# Rebuild INQ's `_pinq` python extension against the CURRENT venv interpreter,
# then re-run the 4 python-interface ctest cases.
#
# WHY: the venv was rebuilt 3.6.8 -> 3.11.13 on 2026-07-30. The existing
# extension is `_pinq.cpython-36m-x86_64-linux-gnu.so` (built 2026-07-29 18:19
# against the OLD 3.6.8 venv). A cpython-36m module cannot be imported by 3.11,
# so inq::tests::gaas.py and na2+.py failed with
#   ModuleNotFoundError: No module named '_pinq'
# Those two failures are a PACKAGING consequence of the venv swap — not CUDA,
# not physics. (gaas.py / na2+.py start with a bare `import pinq`, no ASE guard,
# so rebuilding the extension is sufficient to fix them.)
#
# NOTE: ENABLE_PYTHON only gates INSTALLING the python API (its cache docstring
# is "Install Python API interface"); the `_pinq` target builds regardless. So the
# fix is to re-point the CACHED python, not to flip that option.
#
# pybind11's classic search caches the interpreter, the library AND the version
# (PYBIND11_PYTHON_VERSION was pinned to 3.6.8, PYTHON_LIBRARY to
# libpython3.6m.so). Those must be UNSET or cmake keeps resolving 3.6 headers.
#
# MUST NOT run while the main ctest job is active — reconfiguring the shared
# inq/build tree under a running ctest can corrupt it.

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2
  exit 1
fi
cd "$REPO_ROOT" || exit 1
. "$REPO_ROOT/shared/bin/csd3-env.sh" || { echo "csd3-env.sh failed"; exit 1; }

BUILD="$REPO_ROOT/inq/build"
PY="$REPO_ROOT/venv/bin/python"

echo "=============================================="
echo " node : $(hostname)"
echo " gpu  : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo " venv : $($PY --version 2>&1)"
echo " start: $(date)"
echo "=============================================="

echo "--- BEFORE ---"
ls -la "$BUILD"/python/_pinq*.so 2>/dev/null || echo "  (none)"

FC_ARGS=()
for d in "$REPO_ROOT"/inq/build/_deps/*-src; do
  [[ -d "$d" ]] || continue
  n=$(basename "$d"); n="${n%-src}"
  FC_ARGS+=("-DFETCHCONTENT_SOURCE_DIR_$(printf '%s' "$n" | tr '[:lower:]' '[:upper:]')=$d")
done

echo
echo "--- reconfigure: pin the CLASSIC pybind11 python variables to 3.11 ---"
#
# ATTEMPT 1 FAILED (job 32396176) and made things WORSE. It merely `-U`'d the stale
# entries and set the MODERN vars (Python_EXECUTABLE / Python3_EXECUTABLE). But this
# build has PYBIND11_FINDPYTHON=OFF, i.e. pybind11 uses the CLASSIC FindPythonLibs
# path, which does NOT read those — both landed in the cache as `UNINITIALIZED` and
# were ignored. With the 3.6 pin removed and nothing to replace it, the classic search
# fell back to the oldest interpreter on the box, Python 2.7.18:
#     PYTHON_EXECUTABLE   = /usr/bin/python2.7
#     PYTHON_INCLUDE_DIRS = /usr/include/python2.7   (PY_VERSION "2.7.18")
#     PYBIND11_PYTHON_VERSION = 2.7
# and the compile died on pybind11/detail/common.h:277
#     error "PYTHON < 3.6 IS UNSUPPORTED."
# i.e. 3.6-wrong became 2.7-rejected.
#
# So set the CLASSIC trio explicitly (all three verified present):
#   PYTHON_EXECUTABLE / PYTHON_INCLUDE_DIR / PYTHON_LIBRARY  + PYBIND11_PYTHON_VERSION
# The venv is built on the system 3.11, so its headers ARE /usr/include/python3.11
# (confirmed via sysconfig.get_paths()['include']).
PY_INC="$(${PY} -c 'import sysconfig;print(sysconfig.get_paths()["include"])')"
PY_LIB="/usr/lib64/libpython3.11.so"
PY_VER="$(${PY} -c 'import sys;print("%d.%d"%sys.version_info[:2])')"
echo "  interpreter : $PY"
echo "  include dir : $PY_INC"
echo "  library     : $PY_LIB"
echo "  version     : $PY_VER"

cmake -S "$REPO_ROOT/inq" -B "$BUILD" \
  -DPYTHON_EXECUTABLE="$PY" \
  -DPYTHON_INCLUDE_DIR="$PY_INC" \
  -DPYTHON_LIBRARY="$PY_LIB" \
  -DPYBIND11_PYTHON_VERSION="$PY_VER" \
  -DPython_EXECUTABLE="$PY" \
  -DPython3_EXECUTABLE="$PY" \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_CUDA_COMPILER="$INQ_NVCC" \
  "${FC_ARGS[@]}" 2>&1 | grep -iE 'python|pybind|error|Configuring|Generating' | tail -20
rc_cfg=${PIPESTATUS[0]}
echo "configure exit=$rc_cfg"

echo
echo "--- POST-CONFIGURE GUARD: what did pybind11 actually resolve? ---"
# Fail fast. Attempt 1 spent 4 minutes compiling only to die on a header mismatch.
# The cache is the authority on what nvcc will be handed, so check it, not our intent.
RES_INC=$(grep -m1 '^PYTHON_INCLUDE_DIRS:INTERNAL=' "$BUILD/CMakeCache.txt" | cut -d= -f2-)
RES_VER=$(grep -m1 '^PYTHON_VERSION:INTERNAL='      "$BUILD/CMakeCache.txt" | cut -d= -f2-)
RES_LIB=$(grep -m1 '^PYTHON_LIBRARIES:INTERNAL='    "$BUILD/CMakeCache.txt" | cut -d= -f2-)
echo "  resolved include : ${RES_INC:-<unset>}"
echo "  resolved version : ${RES_VER:-<unset>}"
echo "  resolved library : ${RES_LIB:-<unset>}"
case "${RES_VER:-}" in
  3.1[0-9]*|3.[6-9]*)
    echo "  OK: pybind11 resolved a supported Python (>=3.6)." ;;
  *)
    echo "ERROR: pybind11 resolved Python '${RES_VER:-<unset>}' (include: ${RES_INC:-<unset>})." >&2
    echo "       Refusing to compile — this is the 2.7 fallback that killed attempt 1." >&2
    exit 1 ;;
esac

echo
echo "--- build target _pinq (bounded -j4: heavy CUDA TU) ---"
t0=$SECONDS
cmake --build "$BUILD" --target _pinq -j4 > "$BUILD/pinq-build.log" 2>&1
rc_b=$?
echo "build exit=$rc_b elapsed=$(( SECONDS - t0 ))s"
if [[ $rc_b -ne 0 ]]; then
  # Report the COMPILER diagnostic, not just gmake's exit lines. Attempt 1's output
  # showed only "gmake: *** [_pinq] Error 2" because the grep matched 'Error ' first
  # and head-20 truncated before the real cause (a pybind11 #error) appeared.
  echo "--- build errors (compiler diagnostics first) ---"
  grep -nE '(^|[^A-Za-z])error:|#error|fatal error' "$BUILD/pinq-build.log" | head -15
  echo "--- last 15 lines of build log ---"
  tail -15 "$BUILD/pinq-build.log"
fi

echo
echo "--- AFTER ---"
ls -la "$BUILD"/python/_pinq*.so 2>/dev/null || echo "  (none)"

echo
echo "--- import smoke test (the exact failure mode) ---"
PYTHONPATH="$BUILD/python:$REPO_ROOT/inq/python${PYTHONPATH:+:$PYTHONPATH}" \
  "$PY" -c "import _pinq; print('  _pinq OK'); import pinq; print('  pinq  OK')" 2>&1 | tail -5

echo
echo "--- rerun the 4 python ctest cases ---"
cd "$BUILD" || exit 1
ctest -R '\.py' --output-on-failure --timeout 2000 \
      --output-junit "$BUILD/pinq-retest.junit.xml" 2>&1 | tail -30

echo
echo "=============================================="
echo " DONE — $(date)"
echo "=============================================="
