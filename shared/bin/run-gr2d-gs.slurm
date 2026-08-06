#!/bin/bash
#SBATCH --job-name=gr2d-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=10:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-gr2d-gs.slurm — Phase 0 GS battery for the two-definitions campaign
# (bilayer + monolayer graphene). Plan: docs/plans/real-material-stopping-comparison.md
#
#   Usage:  sbatch shared/bin/run-gr2d-gs.slurm build     # compile only
#           sbatch shared/bin/run-gr2d-gs.slurm           # full battery
#
# Battery (sequential, one GPU job): bilayer cutoff probes {35,45} Ha, then
# production bi@50 + mono@50 (checkpoints saved), then bi@50 Lz=100 vacuum probe.
# Each variant writes results/<tag>/run_summary.txt; completed variants are
# skipped on resubmit (idempotent resume).
set -uo pipefail

STAGE="${1:-prod}"
case "$STAGE" in build|prod) ;; *) echo "FATAL: stage 'build'|'prod', got '$STAGE'" >&2; exit 2 ;; esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 1
fi
. "$REPO_ROOT/shared/bin/csd3-env.sh"
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
. "$REPO_ROOT/shared/config.sh"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/graphene/scripts/twodef_sv/gs"
cd "$RUN_DIR"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  gr2d-gs [$STAGE]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " started: $(date)"
echo "=============================================="

if [ "$STAGE" = "build" ]; then
  inq-run || true
  if [ -x ./run ] && [ ./run -nt run.cpp ]; then
    echo " [build] ./run compiled ($(stat -c%s ./run) bytes)"; exit 0
  else
    echo " [build] FAILED: ./run not produced." >&2; exit 1
  fi
fi

# Compile once if needed (inq-run builds then execs the DEFAULT variant, which
# is battery item 3 — bi@50 probe tag; harmless duplicate, skipped below).
if [ ! -x ./run ] || [ run.cpp -nt ./run ]; then
  echo " [battery] binary missing/stale — building via inq-run (runs default variant)"
  GR_VARIANT=bi GR_CUTOFF_HA=50 GR_LZ_BOHR=80 GR_SAVE=0 inq-run || {
    [ -x ./run ] || { echo " [battery] build failed" >&2; exit 1; }
  }
fi

run_variant () {  # args: VARIANT CUTOFF LZ SAVE
  local v="$1" c="$2" lz="$3" save="$4"
  local tag="${v}_cut${c}_Lz${lz}"
  # skip only if complete AND (for save variants) the checkpoint actually exists
  # (2026-08-05 bug: probe run marked bi_cut50_Lz80 complete, save was skipped)
  local ckpt="$REPO_ROOT/ResearchProject/systems/graphene/shared_gs/gs_3x2_${tag}"
  if [ -f "results/$tag/run_summary.txt" ] \
     && grep -q "run_completed = true" "results/$tag/run_summary.txt" \
     && { [ "$save" != "1" ] || [ -d "$ckpt" ]; }; then
    echo " [skip] $tag already complete"; return 0
  fi
  echo " [run ] $tag (save=$save)  $(date)"
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
      GR_VARIANT="$v" GR_CUTOFF_HA="$c" GR_LZ_BOHR="$lz" GR_SAVE="$save" ./run
  local st=$?
  echo " [done] $tag exit=$st  $(date)"
  return $st
}

overall=0
run_variant bi   35  80 0 || overall=$?
run_variant bi   45  80 0 || overall=$?
run_variant bi   50  80 1 || overall=$?    # production bilayer GS (checkpoint)
run_variant mono 50  80 1 || overall=$?    # production monolayer GS (checkpoint)
run_variant bi   50 100 0 || overall=$?    # vacuum-length probe

echo "=============================================="
echo " [gr2d-gs] battery exit : $overall   finished $(date)"
for f in results/*/run_summary.txt; do
  echo " --- $f ---"; grep -E "run =|cutoff|energy_per_atom|ground_state_energy|electron_count|run_completed" "$f"
done
echo "=============================================="
exit $overall
