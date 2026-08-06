#!/bin/bash
#SBATCH --job-name=slabks-wp
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
# run-slab-ks-wp.slurm — one (density, velocity) point of the CAP-free
# wrap-around slab KS-orbital stopping study.
# Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
#
#   Usage:  sbatch shared/bin/run-slab-ks-wp.slurm smoke     # build + t=0 gates
#           sbatch --array=0-7 shared/bin/run-slab-ks-wp.slurm
#           sbatch shared/bin/run-slab-ks-wp.slurm 5         # single point
#           LJ_RESUME=1 sbatch shared/bin/run-slab-ks-wp.slurm 0   # extend
#
# ARRAY INDEX = density_block * 4 + velocity_index:
#
#   idx  N    r_s    v      N_steps   t (a.u.)  save/  wf/   ckpt/   run name
#     0  100  4.18   2.0     4529      181.2     38   228    1509    n100_v2p0
#     1  100  4.18   2.5     3623      144.9     30   180    1207    n100_v2p5
#     2  100  4.18   3.0     3020      120.8     25   150    1006    n100_v3p0
#     3  100  4.18   3.5     2588      103.5     22   132     862    n100_v3p5
#     4   40  5.67   2.0     4529      181.2     38   228    1509    n40_v2p0
#     5   40  5.67   2.5     3623      144.9     30   180    1207    n40_v2p5
#     6   40  5.67   3.0     3020      120.8     25   150    1006    n40_v3p0
#     7   40  5.67   3.5     2588      103.5     22   132     862    n40_v3p5
#
# CADENCES REDUCED 2026-07-31 (user instruction) after /rds hit 100 % and killed
# 11 of 16 runs mid-flight with a VTI write failure. Density frames 301 -> ~119
# per run, wavefunctions 100 -> ~19, retained checkpoints 5 -> 3 (ckpt = N/3, and
# the FINAL step is always written as a step-stamped snapshot that pruning can
# never remove). Footprint per run drops roughly 3x.
#
# N_steps = 1.5 x the CAP-free classical campaign's step count (user instruction
# 2026-07-31), so every velocity covers the SAME 362.3 Bohr of path = 4.26 box
# lengths ~ 14.5 slab crossings. At r_s = 4.18 (T_plasmon = 31.0 a.u.) the
# v = 2.0 run spans 5.8 plasma periods — the first time in this line of work that
# the 5-period wake criterion is actually reachable.
#
# NO CAP (LJ_CAP_ETA=0). The packet wraps on the FFT grid, which it always did;
# the CAP existed only to stop it. With the CAP gone the Hamiltonian is Hermitian
# and time-independent, so energy_total is CONSERVED — the run gates on that.
#
# sigma_WP = 2.0 throughout, matching the bulk KS study exactly. At that width
# sigma_p = 1/(sqrt2 sigma) = 0.354 against k_Nyq = pi/0.40 = 7.85, so momentum
# aliasing is nil at every velocity here (and would be at 4.0/4.5 too — those two
# are dropped by user instruction, not by a numerical limit).
#
# The smoke stage BUILDS the binary; array tasks then exec ./run directly, so
# eight concurrent tasks never race on one CMake tree. Chain with afterok.
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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/slab_ks_wrap/wp"
GS_ROOT="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs"

#            vi:   0      1      2      3
V_LIST=(    2.0    2.5    3.0    3.5 )
N_LIST=(   4529   3623   3020   2588 )
SAVE_LIST=(  38     30     25     22 )
WF_LIST=(   228    180    150    132 )
CKPT_LIST=(1509   1207   1006    862 )
NE_LIST=(  100     40 )               # density blocks

if [ "$STAGE" = "smoke" ]; then
  NAME="smoke"; NE=100; V="2.0"; NSTEPS=20; SAVE=5; WF=10; CKPT=10
else
  case "$STAGE" in (*[!0-9]*|"") echo "FATAL: bad index '$STAGE'" >&2; exit 2 ;; esac
  [ "$STAGE" -ge 0 ] && [ "$STAGE" -le 7 ] || { echo "FATAL: index out of range: $STAGE" >&2; exit 2; }
  DB=$(( STAGE / 4 )); VI=$(( STAGE % 4 ))
  NE="${NE_LIST[$DB]}"
  V="${V_LIST[$VI]}"; NSTEPS="${N_LIST[$VI]}"
  SAVE="${SAVE_LIST[$VI]}"; WF="${WF_LIST[$VI]}"; CKPT="${CKPT_LIST[$VI]}"
  NAME="n${NE}_v$(echo "$V" | tr '.' 'p')"
fi

# CAP'd variants get a _cap suffix so they cannot overwrite this study's CAP-FREE
# results, which are the whole point of slab_ks_wrap and are already complete.
# (User decision 2026-08-01: the density S(v) plot must contain ONLY CAP'd runs.
# E_absorbed = E_total(t_f) - E_GS is uninformative without an absorber: with the
# packet still in the box it just returns the projectile's own energy, identically
# 58.82 eV at BOTH densities, carrying no stopping information whatsoever.)
if [ "${LJ_CAP_ETA:-0}" != "0" ]; then
  NAME="${NAME}_cap"
fi

GS_DIR="$GS_ROOT/slab_n${NE}_L35x35x85_dx0p4_per2"
if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-slab-ks-gs.slurm must succeed first." >&2
  exit 2
fi

RS=$(awk "BEGIN{n=$NE/30625.0; printf \"%.3f\", (3.0/(4.0*3.14159265358979*n))^(1.0/3.0)}")

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  slabks-wp [$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " system : N=$NE  n0=$(awk "BEGIN{printf \"%.6e\", $NE/30625.0}")  r_s=$RS"
echo " sigma  : 2.0 Bohr (density std 1.4142, spread rate 0.354 Bohr/a.u.)"
echo " v = k0 : $V     N_steps : $NSTEPS   (t = $(awk "BEGIN{printf \"%.2f\", $NSTEPS*0.04}") a.u.,"
echo "          path = $(awk "BEGIN{printf \"%.1f\", $V*$NSTEPS*0.04}") Bohr = $(awk "BEGIN{printf \"%.2f\", $V*$NSTEPS*0.04/85.0}") box lengths)"
echo " CAP    : OFF  -> packet wraps; energy_total must be conserved"
echo " cadence: density/$SAVE  wavefn/$WF  ckpt/$CKPT  stats/1"
echo " resume : ${LJ_RESUME:-0}"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export LJ_OUT="$NAME" LJ_K0="$V" LJ_N_STEPS="$NSTEPS" LJ_SIGMA="${LJ_SIGMA:-2.0}" \
       LJ_N="$NE" LJ_CAP_ETA="${LJ_CAP_ETA:-0}" LJ_SPACING="${LJ_SPACING:-0.40}" \
       LJ_SAVE_EVERY="$SAVE" LJ_WF_EVERY="$WF" LJ_CKPT_EVERY="$CKPT" \
       LJ_LAUNCH_Z="${LJ_LAUNCH_Z:--24.0}" LJ_DT="${LJ_DT:-0.04}" \
       LJ_GS_DIR="$GS_DIR" LJ_RESUME="${LJ_RESUME:-0}"

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
[ -f "$S" ] && { echo " --- run_summary.txt ---"; cat "$S"; }
for f in "results/$NAME/raw/observables/wp_momentum_stats.csv" \
         "results/$NAME/raw/observables/wp_real_space_stats.csv" \
         "results/$NAME/raw/observables/wp_slab_occupancy.csv" \
         "results/$NAME/raw/observables/interactions.csv"; do
  [ -f "$f" ] && { echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"; head -3 "$f"; tail -1 "$f"; }
done
echo " --- retained checkpoints ---"
ls -d "results/$NAME"/ckpt_step* 2>/dev/null || echo "   (none yet)"
echo "=============================================="
exit $status
