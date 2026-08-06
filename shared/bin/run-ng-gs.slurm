#!/bin/bash
#SBATCH --job-name=ng-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-ng-gs.slurm — ground state for the Nazarov-Gross mass ladder.
# Plan: docs/plans/nazarov-gross-slab-mass-ladder.md  (Phase 0, step 2)
#
#   Usage:  sbatch shared/bin/run-ng-gs.slurm build   # compile only, no GS needed
#           sbatch shared/bin/run-ng-gs.slurm         # full SCF + gates + save
#
# System: 30 x 30 x 120 Bohr, periodicity(2), dx 0.50, 15-Bohr slab,
#         N = 206 -> r_s = 2.5011, v_F = 0.767, E_F = 8.01 eV, LDA.
# The binary's own gates (electron count, r_s, finite E) refuse to SAVE a bad GS,
# so a nonzero exit here must stop the whole chain — the orchestrator chains on it.
set -uo pipefail

STAGE="${1:-prod}"
case "$STAGE" in build|prod) ;; *) echo "FATAL: stage must be 'build' or 'prod'; got '$STAGE'" >&2; exit 2 ;; esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 1
fi

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

SYS="$REPO_ROOT/ResearchProject/systems/localised_jellium"
RUN_DIR="$SYS/scripts/ng_mass_ladder/gs"
export GS_DIR="${GS_DIR:-$SYS/shared_gs/slab_n206_L30x30x120_rs2p5_dx0p5_per2}"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  ng-gs [$STAGE]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " system : slab r_s=2.5011, 30x30x15 Bohr slab in a 30x30x120 box, N=206"
echo " grid   : dx 0.50 -> 60 x 60 x 240 = 864k points, 123 states"
echo " target : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

if [ "$STAGE" = "build" ]; then
  # IDEMPOTENCE MATTERS HERE MORE THAN FOR THE RT HALVES. `inq-run` builds and
  # then EXECS, and this binary has no prerequisite to stop it, so a build job
  # runs the whole SCF (~48 min of GPU on job 32877570 — which is how the ground
  # state came to exist before P1 ever ran). Re-submitting the build after an
  # unrelated failure would recompute it for nothing, so skip when the artefact
  # is already newer than its source.
  if [ -x ./run ] && [ ./run -nt run.cpp ]; then
    echo " [build] ./run is already newer than run.cpp — skipping (no recompile, no SCF)"
    status=0
  else
    inq-run || true
    if [ -x ./run ] && [ ./run -nt run.cpp ]; then
      echo " [build] ./run compiled successfully ($(stat -c%s ./run) bytes)"; status=0
    else
      echo " [build] FAILED: ./run was not produced from run.cpp." >&2; status=1
    fi
  fi
else
  if [ -d "$GS_DIR" ] && [ -f results/gs/run_summary.txt ] \
     && grep -q "run_completed = true" results/gs/run_summary.txt; then
    echo " [skip] ground state already present and complete at $GS_DIR"
    status=0
  else
    if [ -x ./run ]; then
      env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run
    else
      inq-run
    fi
    status=$?
  fi
fi

echo "=============================================="
echo " [ng-gs/$STAGE] exit status : $status"
echo " finished                   : $(date)"
[ -f results/gs/run_summary.txt ] && { echo " --- run_summary.txt ---"; cat results/gs/run_summary.txt; }
echo "=============================================="
exit $status
