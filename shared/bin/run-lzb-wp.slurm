#!/bin/bash
#SBATCH --job-name=lzb-wp
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=%x-%A_%a-%j.out
#SBATCH --error=%x-%A_%a-%j.out
#
# run-lzb-wp.slurm — one WAVEPACKET point of the lz_bulk_sweep campaign.
# Plan: docs/plans/jellium-slab-extend-Lz.md
#
# THREE MODES:
#   smoke   sbatch --export=ALL,LZB_CFG=s0p5_L15 run-lzb-wp.slurm smoke
#           20 steps; BUILDS the binary (used standalone by the finalizer; the
#           chain's own build+gate pass is run-lzb-smoke.slurm).
#   pilot / manual / repair
#           sbatch --export=ALL,LZB_CFG=s5p0_L35 run-lzb-wp.slurm 2   # v = 3.0
#           velocity INDEX 0-3 for the box named by LZB_CFG. Resume with
#           LJ_RESUME=1 to extend.
#   production array
#           sbatch --array=0-11%4 run-lzb-wp.slurm
#           index = 3*box + slot, boxes cheap-first (s0p5_L15, s5p0_L15,
#           s0p5_L35, s5p0_L35), slot -> v in {2.0, 2.5, 3.5} — v = 3.0 is the
#           PILOT velocity and already done by the time this array is released.
#           %4 caps GPU pressure at 4 concurrent wp tasks.
#
# Box geometry, sigma family, step counts and cadences all come from
# shared/bin/lzb-params.sh (the bash mirror of lzb_boxes.hpp). Run names carry
# the full preset: results/s0p5_L15_v3p0 etc.
set -euo pipefail

STAGE="${1:-}"
ARRAY_ID="${SLURM_ARRAY_TASK_ID:-}"

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
RUN_DIR="$LJ_BASE/scripts/lz_bulk_sweep/wp"

CFG_ORDER=(s0p5_L15 s5p0_L15 s0p5_L35 s5p0_L35)   # cheap-first for %4 arrays
PROD_VIDX=(0 1 3)                                  # v = 2.0, 2.5, 3.5

SMOKE=0
if [ "$STAGE" = "smoke" ]; then
  SMOKE=1
  export LZB_CFG="${LZB_CFG:-s0p5_L15}"
  VIDX=2
elif [ -n "$STAGE" ]; then
  case "$STAGE" in (*[!0-9]*|"") echo "FATAL: bad velocity index '$STAGE'" >&2; exit 2 ;; esac
  [ "$STAGE" -ge 0 ] && [ "$STAGE" -le 3 ] || { echo "FATAL: index out of range: $STAGE" >&2; exit 2; }
  : "${LZB_CFG:?FATAL: velocity-index mode needs LZB_CFG}"
  VIDX="$STAGE"
elif [ -n "$ARRAY_ID" ]; then
  [ "$ARRAY_ID" -ge 0 ] && [ "$ARRAY_ID" -le 11 ] || { echo "FATAL: array index out of range: $ARRAY_ID" >&2; exit 2; }
  export LZB_CFG="${CFG_ORDER[$((ARRAY_ID / 3))]}"
  VIDX="${PROD_VIDX[$((ARRAY_ID % 3))]}"
else
  echo "FATAL: pass 'smoke' or a velocity index 0-3, or submit as --array=0-11" >&2
  exit 2
fi

# shellcheck source=./lzb-params.sh
. "$REPO_ROOT/shared/bin/lzb-params.sh"
GS_DIR="$LJ_BASE/shared_gs/$GS_TAG"

if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-lzb-gs.slurm $LZB_CFG must succeed first." >&2
  exit 2
fi

if [ "$SMOKE" = "1" ]; then
  NAME="${LZB_CFG}_smoke"; NSTEPS=20; SAVE=5; WF=10; CKPT=10
else
  NAME="${LZB_CFG}_v$(echo "$V" | tr '.' 'p')"
fi

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  lzb-wp [$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " box    : 35 x 35 x $LZ Bohr, slab $(awk "BEGIN{print 2*$HALF}") Bohr, N=$N_E, dx=0.40"
echo " sigma  : $SIGMA Bohr (density std $(awk "BEGIN{printf \"%.4f\", $SIGMA/sqrt(2)}"))"
echo " launch : z = $LAUNCH  (standoff $(awk "BEGIN{printf \"%.1f\", -($HALF)-($LAUNCH)}") Bohr from the face at -$HALF)"
echo " v = k0 : $V     N_steps : $NSTEPS   (t = $(awk "BEGIN{printf \"%.2f\", $NSTEPS*0.04}") a.u.)"
echo " cadence: density/$SAVE  wavefn/$WF  ckpt/$CKPT  stats/1"
echo " resume : ${LJ_RESUME:-0}"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export LJ_OUT="$NAME" LJ_K0="$V" LJ_N_STEPS="$NSTEPS" LJ_SIGMA="$SIGMA" \
       LJ_SAVE_EVERY="$SAVE" LJ_WF_EVERY="$WF" LJ_CKPT_EVERY="$CKPT" \
       LJ_GS_DIR="$GS_DIR" LJ_RESUME="${LJ_RESUME:-0}" LJ_LAUNCH_Z="$LAUNCH"

# A resume must not collide with frames written PAST the checkpoint it loads.
if [ "${LJ_RESUME:-0}" = "1" ]; then
  bash "$REPO_ROOT/shared/bin/prune_orphan_frames.sh" "results/$NAME"
fi

if [ "$SMOKE" = "1" ]; then
  inq-run                       # builds the binary, then exec's it
  status=$?
else
  if [ ! -x ./run ]; then
    echo "FATAL: ./run not built. Run the smoke stage first (it builds)." >&2
    exit 2
  fi
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run
  status=$?
fi

echo "=============================================="
echo " [$NAME] exit status : $status"
echo " finished            : $(date)"
S="results/$NAME/run_summary.txt"
if [ -f "$S" ]; then echo " --- run_summary.txt ---"; cat "$S"; fi
for f in "results/$NAME/raw/observables/wp_momentum_stats.csv" \
         "results/$NAME/raw/observables/wp_real_space_stats.csv" \
         "results/$NAME/raw/observables/interactions.csv"; do
  [ -f "$f" ] && { echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"; head -2 "$f"; tail -1 "$f"; }
done
echo " --- retained checkpoints ---"
ls -d "results/$NAME"/ckpt_step* 2>/dev/null || echo "   (none yet)"
echo "=============================================="
exit $status
