#!/bin/bash
#SBATCH --job-name=inq-ctest
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=inq-ctest-%j.out

# Run INQ's OWN test suite (249 tests) on an A100.
#
# This is the broad correctness check on the local CUDA patch in
# inq/external_libs/gpurun/include/gpu/reduce.hpp (`fixed_return`, replacing
# cuda::proclaim_return_type which is broken on CUDA 12.1). That patch sits under
# a large part of INQ, and neither the inqkit suite (32/32) nor the vacuum
# free-Gaussian sweep exercises Hartree / XC / pseudopotentials / PAW.
#
# NO BUILD STEP: build job 32355421 already produced src/all_unit_tests (98 MB,
# runs all 231 unit tests via Catch2 tag filters) and the 11 executables in
# tests/, plus scripts/run_test.sh. Verified all 37 ctest-referenced executables
# exist before submitting. Runtime is the only cost here.
#
# QOS: mphil-...-sl2-gpu => gpu1, priority 5000, 36 h cap (user asked for 5000).
# SL2 has ~364 GPU-hours; this job is charged against that pot.
#
# SERIAL on purpose (no `ctest -j`): the tests share ONE GPU, and concurrent GPU
# tests risk memory contention producing FALSE failures. A slow honest result
# beats a fast ambiguous one.

set -uo pipefail   # NOT -e: a failing test must not abort before results collect.

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root (no shared/bin/csd3-env.sh)." >&2
  echo "       Submit with: cd <repo> && sbatch shared/bin/run-inq-ctest.slurm" >&2
  exit 1
fi
cd "$REPO_ROOT" || exit 1

. "$REPO_ROOT/shared/bin/csd3-env.sh" || { echo "csd3-env.sh failed"; exit 1; }

BUILD="$REPO_ROOT/inq/build"
TMO="${INQ_CTEST_TIMEOUT:-2000}"   # per-test, matches CLAUDE.md guidance

echo "=============================================="
echo " node    : $(hostname)"
echo " gpu     : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -1)"
echo " nvcc    : $($INQ_NVCC --version | tail -1)"
echo " timeout : ${TMO}s per test"
echo " started : $(date)"
echo "=============================================="

# Guard: confirm the engine is linked against the 12.x toolkit, not the 11.4 that
# rhel8/default-amp drags in. A run against a mis-linked binary would produce
# meaningless failures (fault #6 in the handover).
echo "--- linkage guard ---"
objdump -p "$BUILD/src/all_unit_tests" 2>/dev/null | grep NEEDED \
  | grep -iE 'cufft|cudart|cublas|cusolver' || echo "  (no cuda NEEDED lines?!)"
if objdump -p "$BUILD/src/all_unit_tests" 2>/dev/null | grep -q 'libcufft\.so\.10'; then
  echo "ERROR: all_unit_tests links CUDA 11.x cuFFT (libcufft.so.10). Aborting." >&2
  exit 1
fi

cd "$BUILD" || exit 1
echo
echo "--- ctest: ALL $(ctest -N 2>/dev/null | awk '/Total Tests:/{print $3}') tests, serial ---"
t0=$SECONDS
ctest --output-on-failure --timeout "$TMO" \
      --output-junit "$BUILD/inq-ctest.junit.xml" > "$BUILD/inq-ctest.log" 2>&1
rc=$?
echo "ctest exit=$rc  elapsed=$(( SECONDS - t0 ))s"
echo

echo "=== SUMMARY ==="
grep -E 'tests passed|tests failed|Total Test time' "$BUILD/inq-ctest.log" | tail -5
echo
echo "=== FAILED / TIMEOUT / NOT RUN (if any) ==="
grep -nE '\*\*\*(Failed|Timeout|Not Run|Exception)' "$BUILD/inq-ctest.log" | head -60
echo
echo "=== slowest 10 ==="
grep -oE 'Passed +[0-9.]+ sec' "$BUILD/inq-ctest.log" >/dev/null 2>&1 && \
  grep -E '^ *[0-9]+/[0-9]+ Test' "$BUILD/inq-ctest.log" \
    | sed -E 's/.*Test +#([0-9]+): +([^ ]+) +\.*.* ([0-9.]+) sec.*/\3 \2/' \
    | sort -rn | head -10
echo
echo "=============================================="
echo " DONE — $(date)"
echo " log  : $BUILD/inq-ctest.log"
echo " junit: $BUILD/inq-ctest.junit.xml"
echo "=============================================="
