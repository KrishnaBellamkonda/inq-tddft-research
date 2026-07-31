#!/bin/bash
#SBATCH --job-name=slabks-cl
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
# run-slab-ks-classical.slurm — the WRAPPED CLASSICAL TWIN of one (density,
# velocity) point of the CAP-free wrap-around slab KS-stopping study.
# Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
#
#   Usage:  sbatch shared/bin/run-slab-ks-classical.slurm smoke
#           sbatch --array=0-7 shared/bin/run-slab-ks-classical.slurm
#
# Index layout, run names, step counts and cadences are IDENTICAL to
# run-slab-ks-wp.slurm — that is what makes the two halves comparable point for
# point. The only physical difference is the projectile representation:
#
#   WP half        : a Gaussian electron WAVEPACKET, an occupied KS orbital
#   classical half : a moving Gaussian CHARGE (perturbation), velocity-Verlet
#                    Ehrenfest, mass 1, charge -1
#
# Both at sigma_WP = 2.0. The classical charge std is the derived
# sigma_pot = sigma_WP/sqrt2 = 1.41421 (.claude/rules/sigma-wp-convention.md) —
# it lives inside the binary and is never used as a label.
#
# WRAPPED (LJ_WRAP_Z=1): the projectile re-enters at the opposite face, so it
# makes the same repeated slab crossings the wavepacket makes, and its charge is
# built from the MINIMUM-IMAGE displacement so it wraps smoothly instead of being
# clipped by the face (which would cost 54 % of its charge at the crossing --
# measured, test_gaussian_minimum_image_engine).
#
# CADENCES REDUCED 2026-07-31 (user instruction): density frames ~119/run,
# retained checkpoints 3 (N/3), final step always step-stamped.
#
# EXPECTED, NOT A BUG: a mass-1 electron at v = 2.0 carries 54 eV and deposits
# ~27 eV per slab crossing, so it STOPS after about two crossings and then sits
# (.claude/rules/light-projectile-stopping.md). S comes from the INITIAL DRAG.
set -euo pipefail

STAGE="${1:-${SLURM_ARRAY_TASK_ID:-}}"
if [ -z "$STAGE" ]; then
  echo "FATAL: pass 'smoke' or an index 0-7 (or submit as --array=0-7)" >&2
  exit 2
fi

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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/slab_ks_wrap/classical"
GS_ROOT="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs"

V_LIST=(    2.0    2.5    3.0    3.5 )
N_LIST=(   4529   3623   3020   2588 )
SAVE_LIST=(  38     30     25     22 )
CKPT_LIST=(1509   1207   1006    862 )
NE_LIST=(  100     40 )

if [ "$STAGE" = "smoke" ]; then
  NAME="smoke"; NE=100; V="2.0"; NSTEPS=20; SAVE=5; CKPT=10
else
  case "$STAGE" in (*[!0-9]*|"") echo "FATAL: bad index '$STAGE'" >&2; exit 2 ;; esac
  [ "$STAGE" -ge 0 ] && [ "$STAGE" -le 7 ] || { echo "FATAL: index out of range: $STAGE" >&2; exit 2; }
  DB=$(( STAGE / 4 )); VI=$(( STAGE % 4 ))
  NE="${NE_LIST[$DB]}"
  V="${V_LIST[$VI]}"; NSTEPS="${N_LIST[$VI]}"
  SAVE="${SAVE_LIST[$VI]}"; CKPT="${CKPT_LIST[$VI]}"
  NAME="n${NE}_v$(echo "$V" | tr '.' 'p')"
fi

GS_DIR="$GS_ROOT/slab_n${NE}_L35x35x85_dx0p4_per2"
if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-slab-ks-gs.slurm must succeed first." >&2
  exit 2
fi

RS=$(awk "BEGIN{n=$NE/30625.0; printf \"%.3f\", (3.0/(4.0*3.14159265358979*n))^(1.0/3.0)}")

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  slabks-cl [$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " system : N=$NE  r_s=$RS   sigma_WP=2.0 (sigma_pot=1.41421)"
echo " v = k0 : $V     N_steps : $NSTEPS   (t = $(awk "BEGIN{printf \"%.2f\", $NSTEPS*0.04}") a.u.)"
echo " wrap   : ON (minimum-image Gaussian charge)"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export LJ_OUT="$NAME" LJ_K0="$V" LJ_N_STEPS="$NSTEPS" LJ_SIGMA="${LJ_SIGMA:-2.0}" \
       LJ_N="$NE" LJ_LX=35 LJ_LY=35 LJ_LZ=85 LJ_HALF=12.5 LJ_EDGE_W=1.0 \
       LJ_PERIODICITY=2 LJ_SPACING="${LJ_SPACING:-0.40}" LJ_EXTRA_STATES=24 \
       LJ_LAUNCH_Z="${LJ_LAUNCH_Z:--24.0}" LJ_DT="${LJ_DT:-0.04}" LJ_MASS=1.0 \
       LJ_WRAP_Z="${LJ_WRAP_Z:-1}" LJ_CONST_V=0 \
       LJ_SAVE_EVERY="$SAVE" LJ_CKPT_EVERY="$CKPT" \
       LJ_GS_DIR="$GS_DIR" LJ_RESUME="${LJ_RESUME:-0}"

if [ "$STAGE" = "smoke" ]; then
  inq-run
  status=$?
else
  if [ ! -x ./run ]; then
    echo "FATAL: ./run not built. Submit the 'smoke' stage first (it builds)." >&2
    exit 2
  fi
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run
  status=$?
fi

echo "=============================================="
echo " [$NAME] exit status : $status"
echo " finished            : $(date)"
S="results/$NAME/run_summary.txt"
[ -f "$S" ] && { echo " --- run_summary.txt ---"; cat "$S"; }
for f in "results/$NAME/raw/observables/projectile.csv" \
         "results/$NAME/raw/observables/observables.csv"; do
  [ -f "$f" ] && { echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"; head -3 "$f"; tail -1 "$f"; }
done
echo " --- wrap events ---"
grep -c "\[wrap\]" "slabks-cl-${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}-${SLURM_JOB_ID:-0}.out" 2>/dev/null || true
echo "=============================================="
exit $status
