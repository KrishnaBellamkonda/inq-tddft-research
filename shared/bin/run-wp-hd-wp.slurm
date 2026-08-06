#!/bin/bash
#SBATCH --job-name=wp-hd
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
# run-wp-hd-wp.slurm — one wavepacket velocity point of the WP twin of the
# high-density classical S(v) benchmark (campaign classical-highdensity-sv).
#
#   Usage:  sbatch shared/bin/run-wp-hd-wp.slurm smoke      # build + 20-step gate
#           sbatch --array=0-3 shared/bin/run-wp-hd-wp.slurm
#           sbatch shared/bin/run-wp-hd-wp.slurm 3          # single point, v=3.5
#           LJ_RESUME=1 sbatch shared/bin/run-wp-hd-wp.slurm 0   # extend v=2.0
#           sbatch --export=ALL,LJ_SIGMA=2.0 --array=0-3 \
#                  shared/bin/run-wp-hd-wp.slurm             # sigma = 2 campaign
#
# SIGMA CAMPAIGNS (user instruction 2026-07-31). LJ_SIGMA selects sigma_WP; the
# run name is PREFIXED with s<sigma> for every campaign except the original
# sigma = 0.5, whose names stay bare (v2p0, ...) so the completed runs, their
# notebooks and wp_S_summary.csv keep working untouched.
#
#   LJ_SIGMA=0.5  ->  results/v2p0        (legacy, already run)
#   LJ_SIGMA=2.0  ->  results/s2p0_v2p0
#   LJ_SIGMA=3.0  ->  results/s3p0_v2p0
#
# sigma is the ONLY thing that changes: same GS (sigma does not enter the bath),
# same dx, same launch z = -24, same CAPs, same dt/N_steps/cadences. The velocity
# grid and step counts below are therefore reused verbatim, which is what makes
# the three campaigns directly comparable point for point.
#
# NOTE ON THE ALIASING COLUMN BELOW: it is a sigma = 0.5 property. sigma_p =
# 1/(sqrt2 sigma), so sigma = 2 and 3 have sigma_p = 0.354 and 0.236 against
# k_Nyq = pi/0.40 = 7.85 — every moment is exact to <1e-4 %, and v = 4.0/4.5
# would be recoverable in those campaigns (not requested; grid kept at 4 points).
#
# NOTE ON THE -z CAP AT sigma = 3. Launch is z = -24 and the CAP inner edge is
# z = -30, i.e. 6 Bohr of clearance. The t=0 density std is sigma/sqrt2, so the
# clearance is 4.2 std at sigma = 2 (tail in CAP 1e-5, negligible) but only 2.8
# std at sigma = 3, where 0.23 % of the packet starts INSIDE the absorbing band
# and is removed over the first few a.u. That is a small one-off norm loss, it is
# reproduced by the matching vacuum control, and it is reported in the notebooks
# rather than silently absorbed into S.
#
# VELOCITY GRID — the classical campaign's six points, at 1.5x its step count
# (user instruction 2026-07-30). N_cl = ceil(1.4*69/(0.5 v dt)); N_wp = 1.5 N_cl.
#
#   idx  v = k0   N_steps   t (a.u.)   save/  wf/   ckpt/   sigma_pz^2 err @dx=0.40
#     0    2.0      3623     144.9      12    36     724      +0.05 %
#     1    2.5      2898     115.9      10    30     579      +0.26 %
#     2    3.0      2415      96.6       8    24     483      +1.24 %
#     3    3.5      2070      82.8       7    21     414      +5.06 %
#   ------------------------------------------------------------------------
#     4    4.0      1812      72.5       6    18     362     +17.9 %  DROPPED
#     5    4.5      1611      64.4       5    15     322     +55.1 %  DROPPED
#
# INDICES 4 AND 5 ARE NOT PART OF THE PRODUCTION SWEEP (user decision 2026-07-30).
# sigma_WP = 0.5 fixes the packet's momentum width at sigma_p = 1/(sqrt2 sigma) =
# 1.414, so its k-distribution folds at k_Nyq = pi/dx. Measured at dx = 0.5 and
# v = 4.5: <p_z> = 3.44 vs 4.5 (-24 %), sigma_pz^2 = 9.05 vs 2.0 (+353 %). A fold
# model reproduces that and predicts the table above. sigma_pz^2 sets T1 - T2, a
# headline observable, so the two fastest points are excluded rather than
# caveated. They are still runnable (pass the index explicitly) if dx is reduced
# to 0.30, where every moment is accurate to <= 0.11 %.
#
# ckpt/ gives 5 RETAINED numbered checkpoints per run (>= the 4 the user asked
# for) plus the rolling results/<name>/checkpoint that LJ_RESUME=1 loads.
#
# The smoke stage BUILDS the binary; the array tasks then exec ./run directly, so
# six concurrent tasks never race on one CMake tree. Chain them with
# --dependency=afterok (see submit-wp-highdensity-sv.sh).
#
# ENGINE: inq-study — the two CAPs need the complexified scalar potential.
#
# Plan: docs/plans/wavepacket-highdensity-sv-twin.md
set -euo pipefail

STAGE="${1:-${SLURM_ARRAY_TASK_ID:-}}"
if [ -z "$STAGE" ]; then
  echo "FATAL: pass 'smoke' or a velocity index 0-5 (or submit as --array=0-5)" >&2
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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/wp"
SPACING="${LJ_SPACING:-0.40}"          # production grid (see run.cpp header)
GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x85_dx0p4_per2"

if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-wp-hd-gs.slurm must succeed first." >&2
  exit 2
fi

#            idx:   0      1      2      3      4      5
V_LIST=(   2.0    2.5    3.0    3.5    4.0    4.5 )
N_LIST=(   3623   2898   2415   2070   1812   1611 )
SAVE_LIST=(  12     10      8      7      6      5 )
WF_LIST=(    36     30     24     21     18     15 )
CKPT_LIST=( 724    579    483    414    362    322 )

# sigma_WP and the run-name prefix it implies. Bare names for the original 0.5
# campaign; s<sigma>_ prefix otherwise (see the header).
SIGMA="${LJ_SIGMA:-0.5}"
if [ "$(awk "BEGIN{print ($SIGMA == 0.5)}")" = "1" ]; then
  TAG=""
else
  TAG="s$(echo "$SIGMA" | awk '{printf "%.1f", $1}' | tr '.' 'p')_"
fi

# LAUNCH DISTANCE (effective-sigma hypothesis, docs/plans/effective-sigma-near-launch.md).
# The far-launch campaigns all start at z = -24 and the packet is 4.7-8.1 Bohr wide
# by the time it reaches the slab, so their sigma label describes the packet at
# LAUNCH, not during the physics. The near-launch runs start just outside the slab
# face so the packet arrives essentially undispersed at sigma/sqrt2 = 0.354 Bohr.
# The exact z comes from the injection scan (run-wp-hd-scan.slurm), which picks the
# smallest standoff removing < 3 % of the packet to orthogonalisation.
#
# Default -24.0 keeps EVERY existing name and result resolving unchanged; any other
# value prefixes the run name with nl_ (near-launch).
#
# LJ_LAUNCH_Z=auto reads the accepted z from the scan's own artifact, so the
# sweep can be chained behind the scan with --dependency=afterok and no human
# relays the number between the two jobs.
LAUNCH_Z="${LJ_LAUNCH_Z:--24.0}"
if [ "$LAUNCH_Z" = "auto" ]; then
  SCAN_Z_FILE="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/inject_scan/results/scan/accepted_launch_z.txt"
  if [ ! -f "$SCAN_Z_FILE" ]; then
    echo "FATAL: LJ_LAUNCH_Z=auto but the scan produced no accepted launch z." >&2
    echo "       expected: $SCAN_Z_FILE" >&2
    echo "       Run shared/bin/run-wp-hd-scan.slurm first." >&2
    exit 2
  fi
  LAUNCH_Z="$(tr -d '[:space:]' < "$SCAN_Z_FILE")"
  echo "LJ_LAUNCH_Z=auto -> $LAUNCH_Z (from $SCAN_Z_FILE)"
fi
if [ "$(awk "BEGIN{print ($LAUNCH_Z == -24.0)}")" = "1" ]; then
  NLTAG=""
else
  NLTAG="nl_"
fi
TAG="${NLTAG}${TAG}"

if [ "$STAGE" = "smoke" ]; then
  NAME="${TAG}smoke";  V="2.0"; NSTEPS=20; SAVE=5; WF=10; CKPT=10
else
  case "$STAGE" in (*[!0-9]*|"") echo "FATAL: bad index '$STAGE'" >&2; exit 2 ;; esac
  [ "$STAGE" -ge 0 ] && [ "$STAGE" -le 5 ] || { echo "FATAL: index out of range: $STAGE" >&2; exit 2; }
  V="${V_LIST[$STAGE]}"; NSTEPS="${N_LIST[$STAGE]}"
  SAVE="${SAVE_LIST[$STAGE]}"; WF="${WF_LIST[$STAGE]}"; CKPT="${CKPT_LIST[$STAGE]}"
  NAME="${TAG}v$(echo "$V" | tr '.' 'p')"
fi

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  wp-hd [$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " sigma  : $SIGMA Bohr (density std $(awk "BEGIN{printf \"%.4f\", $SIGMA/sqrt(2)}"), spread rate $(awk "BEGIN{printf \"%.3f\", 1/(sqrt(2)*$SIGMA)}") Bohr/a.u.)"
echo " launch : z = $LAUNCH_Z  (standoff $(awk "BEGIN{printf \"%.1f\", -12.5-($LAUNCH_Z)}") Bohr from the face at -12.5)"
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

if [ "$STAGE" = "smoke" ]; then
  # inq-run builds the binary and then exec's it — this is the one-off build.
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
if [ -f "$S" ]; then
  echo " --- run_summary.txt ---"; cat "$S"
fi
for f in "results/$NAME/raw/observables/wp_momentum_stats.csv" \
         "results/$NAME/raw/observables/wp_real_space_stats.csv" \
         "results/$NAME/raw/observables/interactions.csv"; do
  [ -f "$f" ] && { echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"; head -2 "$f"; tail -1 "$f"; }
done
echo " --- retained checkpoints ---"
ls -d "results/$NAME"/ckpt_step* 2>/dev/null || echo "   (none yet)"
echo "=============================================="
exit $status
