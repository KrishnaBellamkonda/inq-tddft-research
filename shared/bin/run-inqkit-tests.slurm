#!/bin/bash
#SBATCH --job-name=inqkit-tests
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=inqkit-tests-%j.out

# Build + run the inqkit Catch2 test suite (both tiers) on an A100.
#
# PREREQUISITE — configure on the LOGIN node first (compute nodes have no
# outbound network; see docs/handovers/csd3-device-setup.md fault #3):
#   cmake -S inq-stack/tests/include        -B inq-stack/tests/include/build ...
#   cmake -S inq-stack/tests/include/engine -B inq-stack/tests/include/engine/build ...
# This job hard-errors if either CMakeCache.txt is missing rather than trying
# to configure here and dying on a git clone.
#
# Tiers (ADR 0001):
#   pure   - 12 tests, no INQ link, CPU, seconds
#   engine - 22 tests, links INQ + CUDA sm_80, GPU; the expensive half
#
# Emits JUnit XML per tier so results can be tabulated exactly rather than
# scraped from console text.

set -uo pipefail   # NOT -e: a failing test must not abort the job before the
                   # results are collected and summarised.

# SLURM copies the batch script to /var/spool/slurm/slurmd/job<ID>/slurm_script, so
# ${BASH_SOURCE[0]} does NOT point into the repo — deriving REPO_ROOT from it yields
# /var/spool/slurm and every path breaks (killed job 32392850 in 3 s). Use
# SLURM_SUBMIT_DIR, matching build-inq.slurm and run-dispersion.slurm. This means the
# job MUST be submitted from the repo root; the marker check below enforces that
# rather than failing later with a confusing missing-file error.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root" >&2
  echo "       (no shared/bin/csd3-env.sh). Submit with: cd <repo> && sbatch shared/bin/run-inqkit-tests.slurm" >&2
  exit 1
fi
cd "$REPO_ROOT" || exit 1

. "$REPO_ROOT/shared/bin/csd3-env.sh" || { echo "csd3-env.sh failed"; exit 1; }

PURE_DIR="$REPO_ROOT/inq-stack/tests/include/build"
ENG_DIR="$REPO_ROOT/inq-stack/tests/include/engine/build"
JOBS="${INQ_BUILD_JOBS:-12}"

# The per-test ENVIRONMENT in the CMakeLists already points these at
# INQ_SOURCE/install/share, but export them too so a manually-run binary works.
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

echo "=============================================="
echo " node    : $(hostname)"
echo " gpu     : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo " nvcc    : $($INQ_NVCC --version | tail -1)"
echo " g++     : $(g++ --version | head -1)"
echo " jobs    : $JOBS"
echo " started : $(date)"
echo "=============================================="

for d in "$PURE_DIR" "$ENG_DIR"; do
  if [[ ! -f "$d/CMakeCache.txt" ]]; then
    echo "ERROR: $d not configured. Run cmake on the LOGIN node first." >&2
    exit 1
  fi
done

run_tier () {          # run_tier <label> <build-dir> <ctest-label> <timeout>
  local label="$1" dir="$2" ctlabel="$3" tmo="$4"
  echo
  echo "########## TIER: $label ##########"
  local t0=$SECONDS
  echo "--- building ($label) ---"
  cmake --build "$dir" -j"$JOBS" > "$dir/build.log" 2>&1
  local brc=$?
  local bt=$(( SECONDS - t0 ))
  echo "build exit=$brc  elapsed=${bt}s  log=$dir/build.log"
  if [[ $brc -ne 0 ]]; then
    echo "--- build FAILED; first errors: ---"
    grep -nE 'error:|Error ' "$dir/build.log" | head -25
    # Still attempt ctest: targets that DID build can still be exercised.
  fi

  echo "--- ctest ($label) ---"
  local t1=$SECONDS
  ( cd "$dir" && ctest -L "$ctlabel" --output-on-failure --timeout "$tmo" \
      --output-junit "$dir/results.junit.xml" ) > "$dir/ctest.log" 2>&1
  local crc=$?
  echo "ctest exit=$crc  elapsed=$(( SECONDS - t1 ))s  log=$dir/ctest.log"
  tail -12 "$dir/ctest.log"
  echo "BUILD_RC_$label=$brc CTEST_RC_$label=$crc BUILD_SECONDS_$label=$bt"
}

# Pure tier first: cheap, gives fast signal even if the engine build is long.
run_tier pure   "$PURE_DIR" pure   120
run_tier engine "$ENG_DIR"  engine 900

echo
echo "=============================================="
echo " ALL TIERS DONE — $(date)"
echo " junit: $PURE_DIR/results.junit.xml"
echo "        $ENG_DIR/results.junit.xml"
echo "=============================================="
