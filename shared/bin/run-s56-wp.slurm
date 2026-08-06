#!/bin/bash
#SBATCH --job-name=s56-wp
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
# run-s56-wp.slurm — one WAVEPACKET velocity point of the sigma_WP = 5/6 twin
# campaign.
#
#   Usage:  sbatch --export=ALL,LJ_SIGMA=6.0 shared/bin/run-s56-wp.slurm smoke
#           sbatch --export=ALL,LJ_SIGMA=6.0 --array=0-3 shared/bin/run-s56-wp.slurm
#           sbatch --export=ALL,LJ_SIGMA=5.0 shared/bin/run-s56-wp.slurm 2   # v=3.0
#           LJ_RESUME=1 sbatch --export=ALL,LJ_SIGMA=6.0,LJ_RESUME=1 \
#                  shared/bin/run-s56-wp.slurm 0    # extend v=2.0
#
# RUN NAMES are always prefixed s<sigma>_ (s5p0_, s6p0_) — unlike the 85-Bohr
# campaigns there is no bare-name legacy to preserve here, so every run in this
# sweep carries its sigma explicitly.
#
# WHAT CHANGES BETWEEN THE TWO CAMPAIGNS: sigma_WP only. Same GS (sigma does not
# enter the bath), same dx = 0.40, same launch z = -27.5, same CAP, same dt and
# step counts, same four velocities. That is what makes the traces comparable
# point for point — and comparable to the classical twins, which take the same
# LJ_SIGMA and derive sigma_pot = sigma_WP/sqrt(2) internally.
#
# VELOCITY GRID / STEP COUNTS. Same convention as the 85-Bohr campaigns,
# re-evaluated for the longer path (launch -27.5, L_z = 105):
#     N_STEPS = round( 4.36 * (|launch_z| + L_z/2) / (v * dt) ),  dt = 0.04
# which reproduces the recorded 3623 at v=2.0/z0=-24/L_z=85 to within one step.
#
#   idx  v = k0   N_steps   t (a.u.)   save/  wf/   ckpt/
#     0    2.0      4360      174.4      14    43     872
#     1    2.5      3488      139.5      12    35     698
#     2    3.0      2907      116.3      10    29     581
#     3    3.5      2491       99.6       8    25     498
#
# The packet clears the slab by t = 20 a.u. and reaches the CAP inner edge (+40)
# by t ~ 34 a.u. at the slowest velocity; the rest is plateau time for the
# deposit estimator S = (E_total(t_f) - E_GS)/25.
#
# NO ALIASING CAVEAT HERE. sigma_p = 1/(sqrt2 sigma) = 0.141 (sigma=5) / 0.118
# (sigma=6) against k_Nyq = pi/0.40 = 7.85, so every momentum moment is exact to
# ~1e-14 %. v = 4.0/4.5 would be perfectly recoverable in this campaign; the grid
# is held at four points for comparability with the existing sigma traces.
#
# ckpt/ gives 5 RETAINED numbered checkpoints per run plus the rolling
# results/<name>/checkpoint that LJ_RESUME=1 loads.
#
# The smoke stage BUILDS the binary; array tasks exec ./run directly, so
# concurrent tasks never race on one CMake tree. Chain with --dependency=afterok.
#
# ENGINE: inq-study — the two CAPs need the complexified scalar potential.
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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/sigma56_sv/wp"
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
WF_LIST=(    43     35     29     25 )
CKPT_LIST=( 872    698    581    498 )

SIGMA="${LJ_SIGMA:-6.0}"
TAG="s$(echo "$SIGMA" | awk '{printf "%.1f", $1}' | tr '.' 'p')_"

# Launch z is COMMON to both halves and both sigma (plan §4). It is exposed as a
# knob only so a future campaign can move it; changing it here without changing
# the classical dispatcher would silently break the twin pairing.
LAUNCH_Z="${LJ_LAUNCH_Z:--27.5}"

if [ "$STAGE" = "smoke" ]; then
  NAME="${TAG}smoke";  V="2.0"; NSTEPS=20; SAVE=5; WF=10; CKPT=10
else
  case "$STAGE" in (*[!0-9]*|"") echo "FATAL: bad index '$STAGE'" >&2; exit 2 ;; esac
  [ "$STAGE" -ge 0 ] && [ "$STAGE" -le 3 ] || { echo "FATAL: index out of range: $STAGE" >&2; exit 2; }
  V="${V_LIST[$STAGE]}"; NSTEPS="${N_LIST[$STAGE]}"
  SAVE="${SAVE_LIST[$STAGE]}"; WF="${WF_LIST[$STAGE]}"; CKPT="${CKPT_LIST[$STAGE]}"
  NAME="${TAG}v$(echo "$V" | tr '.' 'p')"
fi

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  s56-wp [$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " box    : 35 x 35 x 105 Bohr, dx = $SPACING, CAP bands +/-[40, 52.5]"
echo " sigma  : $SIGMA Bohr (density std $(awk "BEGIN{printf \"%.4f\", $SIGMA/sqrt(2)}"), spread rate $(awk "BEGIN{printf \"%.3f\", 1/(sqrt(2)*$SIGMA)}") Bohr/a.u.)"
echo " launch : z = $LAUNCH_Z  (standoff $(awk "BEGIN{printf \"%.1f\", -12.5-($LAUNCH_Z)}") Bohr from the face at -12.5; $(awk "BEGIN{printf \"%.1f\", -40-($LAUNCH_Z)}") Bohr from the CAP at -40)"
echo " v = k0 : $V     N_steps : $NSTEPS   (t = $(awk "BEGIN{printf \"%.2f\", $NSTEPS*0.04}") a.u.)"
echo " cadence: density/$SAVE  wavefn/$WF  ckpt/$CKPT  stats/1"
echo " resume : ${LJ_RESUME:-0}"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export LJ_OUT="$NAME" LJ_K0="$V" LJ_N_STEPS="$NSTEPS" LJ_SIGMA="$SIGMA" \
       LJ_SAVE_EVERY="$SAVE" LJ_WF_EVERY="$WF" LJ_CKPT_EVERY="$CKPT" \
       LJ_GS_DIR="$GS_DIR" LJ_RESUME="${LJ_RESUME:-0}" LJ_SPACING="$SPACING" \
       LJ_LAUNCH_Z="$LAUNCH_Z"

# A resume must not collide with frames written PAST the checkpoint it loads
# (VTIImageDataWriter is overwrite=false by design). Those frames are recomputed
# deterministically by the resumed segment, so they are dropped first.
if [ "${LJ_RESUME:-0}" = "1" ]; then
  bash "$REPO_ROOT/shared/bin/prune_orphan_frames.sh" "results/$NAME"
fi

if [ "$STAGE" = "smoke" ]; then
  inq-run                       # builds the binary, then exec's it
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
for f in "results/$NAME/raw/observables/wp_momentum_stats.csv" \
         "results/$NAME/raw/observables/wp_real_space_stats.csv" \
         "results/$NAME/raw/observables/interactions.csv"; do
  [ -f "$f" ] && { echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"; head -2 "$f"; tail -1 "$f"; }
done
echo " --- retained checkpoints ---"
ls -d "results/$NAME"/ckpt_step* 2>/dev/null || echo "   (none yet)"
echo "=============================================="
exit $status
