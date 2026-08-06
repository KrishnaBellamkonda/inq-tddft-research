#!/bin/bash
#SBATCH --job-name=s56-cl
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=%x-%A_%a-%j.out
#SBATCH --error=%x-%A_%a-%j.out
#
# run-s56-cl.slurm — one CLASSICAL velocity point of the sigma_WP = 5/6 twin
# campaign: the width-matched twin of the run-s56-wp.slurm point at the same
# index and the same LJ_SIGMA.
#
#   Usage:  sbatch --export=ALL,LJ_SIGMA=6.0 shared/bin/run-s56-cl.slurm smoke
#           sbatch --export=ALL,LJ_SIGMA=6.0 --array=0-3 shared/bin/run-s56-cl.slurm
#           sbatch --export=ALL,LJ_SIGMA=6.0,LJ_CAP_ETA=0 \
#                  shared/bin/run-s56-cl.slurm 2      # CAP-FREE control at v=3.0
#
# THE TWIN CONTRACT. Identical to the wavepacket half in EVERY parameter except
# the projectile representation: same GS, box, dx, dt, step counts, cadences,
# CAP, launch z, velocity and sigma LABEL. LJ_SIGMA is sigma_WP; the binary
# derives sigma_pot = sigma_WP/sqrt(2) itself (3.53553 at sigma=5, 4.24264 at
# sigma=6), which is the width at which the classical charge cloud equals the
# packet's t=0 density. Never pass sigma_pot here.
#
# PROJECTILE = DIRECT erf/r POTENTIAL (moving_gaussian_projectile_potential), not
# the Poisson charge perturbation (user decision 2026-08-02). No charge inside the
# cell means the projectile may leave the box with the in-box field simply
# becoming its free-space tail — and no periodic G=0 sheet, hence no exit
# transient in the ledger.
#
# CAP ON BY DEFAULT (user decision 2026-08-02) so E_absorbed/L_slab is the SAME
# estimator on both halves. LJ_CAP_ETA=0 makes the absorbing bands an exact no-op
# and produces the CAP-free control — same binary, so the control cannot drift
# away from production. Control runs are named cl_nocap_*.
#
#   idx  v = k0   N_steps   t (a.u.)   save/  ckpt/
#     0    2.0      4360      174.4      14     872
#     1    2.5      3488      139.5      12     698
#     2    3.0      2907      116.3      10     581
#     3    3.5      2491       99.6       8     498
#
# LIGHT PROJECTILE (.claude/rules/light-projectile-stopping.md): a mass-1 electron
# under free Ehrenfest DECELERATES strongly. That is physics. Nothing here gates
# on velocity drift, and the final v is reported, not checked.
#
# ENGINE: inq-study — the CAP needs the complexified scalar potential.
#
# Plan: docs/plans/sigma56-sv-twin.md
set -euo pipefail

STAGE="${1:-${SLURM_ARRAY_TASK_ID:-}}"
if [ -z "$STAGE" ]; then
  echo "FATAL: pass 'smoke' or a velocity index 0-3 (or submit as --array=0-3)" >&2
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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/sigma56_sv/classical"
SPACING="${LJ_SPACING:-0.40}"
GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x105_dx0p4_per2"

if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-s56-gs.slurm must succeed first." >&2
  exit 2
fi

#            idx:   0      1      2      3
V_LIST=(   2.0    2.5    3.0    3.5 )
N_LIST=(   4360   3488   2907   2491 )
SAVE_LIST=(  14     12     10      8 )
CKPT_LIST=( 872    698    581    498 )

SIGMA="${LJ_SIGMA:-6.0}"
STAG="s$(echo "$SIGMA" | awk '{printf "%.1f", $1}' | tr '.' 'p')_"

# CAP-free control gets its own name prefix, so a control can never be mistaken
# for production in a summary table (they differ ONLY by eta).
CAP_ETA="${LJ_CAP_ETA:--1.0}"
if [ "$(awk "BEGIN{print ($CAP_ETA == 0)}")" = "1" ]; then
  TAG="cl_nocap_${STAG}"
else
  TAG="cl_${STAG}"
fi

LAUNCH_Z="${LJ_LAUNCH_Z:--27.5}"

if [ "$STAGE" = "smoke" ]; then
  NAME="${TAG}smoke";  V="2.0"; NSTEPS=20; SAVE=5; CKPT=10
else
  case "$STAGE" in (*[!0-9]*|"") echo "FATAL: bad index '$STAGE'" >&2; exit 2 ;; esac
  [ "$STAGE" -ge 0 ] && [ "$STAGE" -le 3 ] || { echo "FATAL: index out of range: $STAGE" >&2; exit 2; }
  V="${V_LIST[$STAGE]}"; NSTEPS="${N_LIST[$STAGE]}"
  SAVE="${SAVE_LIST[$STAGE]}"; CKPT="${CKPT_LIST[$STAGE]}"
  NAME="${TAG}v$(echo "$V" | tr '.' 'p')"
fi

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  s56-cl [$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " box    : 35 x 35 x 105 Bohr, dx = $SPACING"
echo " sigma  : sigma_WP = $SIGMA Bohr  ->  sigma_pot = $(awk "BEGIN{printf \"%.5f\", $SIGMA/sqrt(2)}") Bohr (derived in the binary)"
echo " launch : z = $LAUNCH_Z  (standoff $(awk "BEGIN{printf \"%.1f\", -12.5-($LAUNCH_Z)}") Bohr from the face at -12.5)"
echo " CAP    : eta = $CAP_ETA Ha  ($([ "$(awk "BEGIN{print ($CAP_ETA == 0)}")" = "1" ] && echo 'OFF — CAP-free control' || echo 'ON, 12.5 Bohr/face, bands +/-[40, 52.5]'))"
echo " v = k0 : $V     N_steps : $NSTEPS   (t = $(awk "BEGIN{printf \"%.2f\", $NSTEPS*0.04}") a.u.)"
echo " cadence: density/$SAVE  ckpt/$CKPT  stats/1"
echo " resume : ${LJ_RESUME:-0}"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export LJ_OUT="$NAME" LJ_K0="$V" LJ_N_STEPS="$NSTEPS" LJ_SIGMA="$SIGMA" \
       LJ_SAVE_EVERY="$SAVE" LJ_CKPT_EVERY="$CKPT" LJ_CAP_ETA="$CAP_ETA" \
       LJ_GS_DIR="$GS_DIR" LJ_RESUME="${LJ_RESUME:-0}" LJ_SPACING="$SPACING" \
       LJ_LAUNCH_Z="$LAUNCH_Z"

# A resume must not collide with frames written PAST the checkpoint it loads
# (VTIImageDataWriter is overwrite=false by design). Those frames are recomputed
# deterministically by the resumed segment, so they are dropped first.
if [ "${LJ_RESUME:-0}" = "1" ]; then
  bash "$REPO_ROOT/shared/bin/prune_orphan_frames.sh" "results/$NAME"
fi

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
if [ -f "$S" ]; then echo " --- run_summary.txt ---"; cat "$S"; fi
for f in "results/$NAME/raw/observables/projectile.csv" \
         "results/$NAME/raw/observables/interactions.csv" \
         "results/$NAME/raw/observables/observables.csv"; do
  [ -f "$f" ] && { echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"; head -2 "$f"; tail -1 "$f"; }
done
echo " --- density frames ---"
ls "results/$NAME/raw/vti/density_total" 2>/dev/null | wc -l
echo " --- retained checkpoints ---"
ls -d "results/$NAME"/ckpt_step* 2>/dev/null || echo "   (none yet)"
echo "=============================================="
exit $status
