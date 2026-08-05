#!/bin/bash
#SBATCH --job-name=lzb-smoke
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
# run-lzb-smoke.slurm — builds BOTH production binaries (wp, classical) and runs
# a 20-step smoke with full t=0 gates for EVERY box preset x half: 8 smokes in
# one job. Exits non-zero if any build or gate fails, which afterok-blocks the
# whole pilot — a compile error or a bad packet costs one short job, never a
# sweep.
# Plan: docs/plans/jellium-slab-extend-Lz.md
#
# Per-box t=0 gates matter here because the launch clearance differs per box
# (per-family standoff) and sigma = 0.5 has only ~1.25 grid points per sigma —
# the relative-tolerance gates are calibrated for exactly that.
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

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

LJ_BASE="$REPO_ROOT/ResearchProject/systems/localised_jellium"
CFGS=(s0p5_L15 s0p5_L35 s5p0_L15 s5p0_L35)
fails=0

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  lzb-smoke (build + 8 x t=0 gates)"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " started: $(date)"
echo "=============================================="

for HALF_DIR in wp classical; do
  cd "$LJ_BASE/scripts/lz_bulk_sweep/$HALF_DIR"
  first=1
  for CFG in "${CFGS[@]}"; do
    export LZB_CFG="$CFG"
    unset V NSTEPS SAVE WF CKPT
    VIDX=2   # v = 3.0, the pilot velocity
    # shellcheck source=./lzb-params.sh
    . "$REPO_ROOT/shared/bin/lzb-params.sh"
    GS_DIR="$LJ_BASE/shared_gs/$GS_TAG"
    if [ ! -d "$GS_DIR" ]; then
      echo "FAIL [$HALF_DIR/$CFG]: ground state missing at $GS_DIR"
      fails=$((fails+1)); continue
    fi
    echo ""
    echo "### smoke $HALF_DIR/$CFG — 20 steps at v=$V, sigma=$SIGMA ###"
    export LJ_OUT="${CFG}_smoke" LJ_K0="$V" LJ_N_STEPS=20 LJ_SIGMA="$SIGMA" \
           LJ_SAVE_EVERY=5 LJ_WF_EVERY=10 LJ_CKPT_EVERY=10 \
           LJ_GS_DIR="$GS_DIR" LJ_RESUME=0
    if [ "$first" = "1" ]; then
      # inq-run builds (or incrementally rebuilds) then runs — never trust a
      # stale ./run for the first variant (the sigma56 vac lesson).
      if ! inq-run; then echo "FAIL [$HALF_DIR/$CFG]: build/smoke failed"; fails=$((fails+1)); fi
    else
      if ! env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run; then
        echo "FAIL [$HALF_DIR/$CFG]: smoke failed (t=0 gates or runtime)"; fails=$((fails+1))
      fi
    fi
    first=0
  done
done

echo ""
echo "=============================================="
echo " finished : $(date)"
echo " failures : $fails / 8 smokes"
echo "=============================================="
[ "$fails" -eq 0 ]
