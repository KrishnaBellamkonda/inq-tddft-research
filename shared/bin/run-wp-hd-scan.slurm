#!/bin/bash
#SBATCH --job-name=wp-hd-scan
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-wp-hd-scan.slurm — launch-distance scan for the effective-sigma hypothesis.
#
# Finds the SMALLEST standoff from the jellium slab face at which injecting the
# sigma_WP = 0.5 wavepacket removes less than 3 % of its weight to
# orthogonalisation against the occupied bath (user criterion, 2026-08-01).
#
# SEARCH RULE (user decisions, 2026-08-01):
#   * start at launch_z = -14.0  (1.5 Bohr outside the face at -12.5, i.e. just
#     beyond the 1 Bohr erfc softening) — this is a FLOOR, never go closer;
#   * on reject, RETREAT in 0.5 Bohr steps (-14.5, -15.0, ...);
#   * scan at v = 2.0 ONLY. The overlap with the occupied manifold is suppressed
#     as exp(-sigma^2 k0^2), so the SLOWEST velocity is the worst case; the
#     accepted z is then used for all four production velocities so that the
#     four runs share one geometry.
#   * Gaussianity of the momentum profile is measured and reported but NEVER
#     vetoes — the 3 % number alone decides.
#
# Each trial is a GS load + one injection. NOTHING is propagated: inqkit
# orthonormalises exactly once, at injection, and the ETRS propagator these runs
# use never re-orthogonalises. So a trial costs ~a minute, not a run.
#
# The scan program's EXIT CODE is the decision: 0 = accept, 3 = reject, other =
# error (a closure failure in the loss bookkeeping aborts the whole scan).
#
# Usage:  cd <repo> && sbatch shared/bin/run-wp-hd-scan.slurm
#         sbatch --export=ALL,LJ_K0=3.5 shared/bin/run-wp-hd-scan.slurm
#
# Plan: docs/plans/effective-sigma-near-launch.md
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root. Submit from the repo root." >&2
  exit 1
fi
cd "$REPO_ROOT" || exit 1

. "$REPO_ROOT/shared/bin/csd3-env.sh"

export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
. "$REPO_ROOT/shared/config.sh"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/inject_scan"
GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x85_dx0p4_per2"

if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR" >&2
  exit 2
fi

# Search parameters
Z_START="${LJ_Z_START:--14.0}"     # floor: 1.5 Bohr outside the face
Z_STEP="${LJ_Z_STEP:-0.5}"         # retreat step, Bohr
Z_MAX_TRIALS="${LJ_Z_MAX_TRIALS:-14}"   # -14.0 ... -20.5; far-launch is -24
K0="${LJ_K0:-2.0}"                 # worst-case velocity
SIGMA="${LJ_SIGMA:-0.5}"
TOL="${LJ_REMOVED_TOL_PC:-3.0}"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  wp-hd-scan"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " search : z from $Z_START, step -$Z_STEP, max $Z_MAX_TRIALS trials"
echo " wp     : sigma=$SIGMA  k0=$K0   criterion: removed < $TOL %"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export LJ_GS_DIR="$GS_DIR" LJ_SIGMA="$SIGMA" LJ_K0="$K0" \
       LJ_SPACING="${LJ_SPACING:-0.40}" LJ_REMOVED_TOL_PC="$TOL"

# ---- build once ------------------------------------------------------------
# inq-run builds AND runs; do the build with a throwaway trial at the far-launch
# point z = -24, which doubles as a REGRESSION CHECK: that is the geometry the
# completed campaign used, and its recorded max_overlap is 3.69e-4
# (wp/results/v2p0/raw/observables/wp_config.txt). If the build is sane the scan
# program must reproduce that number.
echo ""
echo "--- build + far-launch regression trial (z = -24, expect max_overlap ~3.69e-4) ---"
LJ_LAUNCH_Z=-24.0 LJ_TRIAL="regression_z24" inq-run
build_status=$?
if [ ! -x ./run ]; then
  echo "FATAL: ./run was not built (inq-run exit $build_status)" >&2
  exit 2
fi

# ---- the scan --------------------------------------------------------------
Z="$Z_START"
ACCEPTED=""
for ((i=0; i<Z_MAX_TRIALS; i++)); do
  echo ""
  echo "--- trial $((i+1)): launch_z = $Z ---"
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
      LJ_LAUNCH_Z="$Z" ./run
  rc=$?
  case $rc in
    0) ACCEPTED="$Z"; echo ">>> ACCEPTED at launch_z = $Z"; break ;;
    3) echo ">>> rejected at launch_z = $Z — retreating ${Z_STEP} Bohr" ;;
    5) echo "FATAL: orthogonalisation-loss closure check FAILED — bookkeeping is wrong." >&2
       exit 5 ;;
    *) echo "FATAL: trial exited $rc" >&2; exit "$rc" ;;
  esac
  Z=$(awk "BEGIN{printf \"%.1f\", $Z - $Z_STEP}")
done

echo ""
echo "=============================================="
if [ -n "$ACCEPTED" ]; then
  STANDOFF=$(awk "BEGIN{printf \"%.1f\", -12.5 - ($ACCEPTED)}")
  echo " SCAN RESULT : ACCEPTED launch_z = $ACCEPTED  (standoff $STANDOFF Bohr)"
  # The single artifact the sweep dispatcher reads.
  echo "$ACCEPTED" > results/scan/accepted_launch_z.txt
  echo " written     : $RUN_DIR/results/scan/accepted_launch_z.txt"
else
  echo " SCAN RESULT : NO launch_z satisfied removed < $TOL % within $Z_MAX_TRIALS trials"
  echo "               (searched $Z_START down to $Z)"
fi
echo " --- scan ledger ---"
[ -f results/scan/scan.csv ] && cat results/scan/scan.csv
echo " finished : $(date)"
echo "=============================================="

[ -n "$ACCEPTED" ] || exit 3
exit 0
