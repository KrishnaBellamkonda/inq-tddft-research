#!/bin/bash
#SBATCH --job-name=bulk-ie-build
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-bulk-ie-build.slurm — PREFLIGHT for the interaction-energy re-run.
#
# Compiles all SIX run definitions that were just wired with the P/S/B pairwise
# decomposition (.claude/rules/decomposed-interaction-energies.md) and FAILS THE
# WHOLE JOB if any one of them does not build.
#
# WHY A SEPARATE BUILD STAGE. The production array is chained on this with
# afterok, so a compile error in ANY half stops ALL six before they consume GPU
# hours -- and, more importantly, prevents a PARTIAL set of twin pairs, which
# would be useless for the classical-vs-WP comparison this re-run exists to make.
#
# The build dirs are already cmake-configured from previous runs and only
# run.cpp changed, so `cmake --build` here is exactly the build step inq-run
# performs (inq-run has no build-only mode -- it execs the binary).
#
# Usage:  sbatch shared/bin/run-bulk-ie-build.slurm
#
# Plan: docs/plans/bulk-jellium-ks-stopping.md
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
JOBS="${SLURM_CPUS_PER_TASK:-8}"
SCRIPTS="$REPO_ROOT/ResearchProject/systems/jellium/scripts"

TARGETS=(
  "bulk_ks_stopping/wp"
  "bulk_ks_stopping/classical"
  "bulk_ks_stopping_rs4/wp"
  "bulk_ks_stopping_rs4/classical"
  "bulk_ks_stopping_rs4_sigma3/wp"
  "bulk_ks_stopping_rs4_sigma3/classical"
)

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  bulk-ie PREFLIGHT BUILD"
echo " node   : $(hostname)"
echo " engine : $INQ_SOURCE"
echo " jobs   : $JOBS"
echo " started: $(date)"
echo "=============================================="

rc=0
for t in "${TARGETS[@]}"; do
  d="$SCRIPTS/$t"
  echo ""
  echo "--- building $t ---"
  if [ ! -f "$d/run.cpp" ]; then
    echo "FATAL: $d/run.cpp missing" >&2; rc=1; continue
  fi
  if ! grep -q "interaction_energies.hpp" "$d/run.cpp"; then
    echo "FATAL: $t is NOT wired with the pairwise decomposition — refusing to" >&2
    echo "       build it, the re-run would produce no interactions.csv." >&2
    rc=1; continue
  fi
  if [ ! -f "$d/build/CMakeLists.txt" ]; then
    echo "FATAL: $d/build is not cmake-configured; run inq-run there once first." >&2
    rc=1; continue
  fi
  if ( cd "$d/build" && cmake --build . --target run -j"$JOBS" ); then
    echo "=== BUILD OK: $t ==="
    ls -la "$d/run" 2>/dev/null
  else
    echo "=== BUILD FAILED: $t ===" >&2
    rc=1
  fi
done

echo ""
echo "=============================================="
echo " preflight build exit : $rc"
echo " finished             : $(date)"
if [ "$rc" -ne 0 ]; then
  echo " ONE OR MORE BUILDS FAILED — the production array will NOT start."
fi
echo "=============================================="
exit $rc
