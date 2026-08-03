#!/bin/bash
#SBATCH --job-name=chan-twin
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-chan-twin.slurm — ONE HALF of the annular-tube channeling twin.
# Plan: docs/plans/cylindrical-channeling-ks-stopping.md
#
#   Usage:  sbatch shared/bin/run-chan-twin.slurm wp smoke        # build + t=0 gates
#           sbatch shared/bin/run-chan-twin.slurm classical smoke
#           sbatch shared/bin/run-chan-twin.slurm wp
#           sbatch shared/bin/run-chan-twin.slurm classical
#           CH_RESUME=1 CH_N_STEPS=3000 sbatch shared/bin/run-chan-twin.slurm wp
#
# The two halves are physically independent, so both can hold a GPU at once; the
# submit script queues them together.
#
# THE PAIR IS THE MEASUREMENT. Every parameter below is shared by both halves and
# comes from shared/configs/channeling_tube_rs3.hpp, so the ONLY difference
# between the runs is that the projectile is an external Gaussian charge in one
# and an occupied KS orbital in the other:
#
#   tube      R_in 10, R_out 14, L_z 60, edge 0.5, fully periodic, dx 0.5
#   bath      N = 160  ->  n0 = 8.84194e-3,  r_s = 3.000000,  hbar w_p = 9.07 eV
#   proj      50 eV  ->  v = k0 = 1.917 a.u.,  v/v_F = 3.00
#   sigma     sigma_WP = 4 Bohr (classical charge std sigma_pot = 2.8284)
#   launch    (0, 0, -28), +z
#   time      1500 x 0.02 = 30 a.u.  ->  57.5 Bohr of path, ONE traversal
#
# SMOKE FIRST. The smoke stage BUILDS the binary and runs 20 steps; the
# production stage exec's ./run directly, so the two halves never race on one
# CMake tree. Chain production afterok the smoke.
set -euo pipefail

HALF="${1:-}"
STAGE="${2:-prod}"
case "$HALF" in
  wp|classical) ;;
  *) echo "FATAL: first argument must be 'wp' or 'classical'; got '$HALF'" >&2; exit 2 ;;
esac
case "$STAGE" in
  build|smoke|prod) ;;
  *) echo "FATAL: second argument must be 'build', 'smoke' or 'prod'; got '$STAGE'" >&2; exit 2 ;;
esac

# STAGES, and why there are three rather than two.
#   build  COMPILE ONLY. Needs no ground state, so it can run the moment the
#          library tests pass and a compile error surfaces in minutes instead of
#          after an 8-hour SCF. `inq-run` builds and then execs the binary, which
#          exits 2 at its own CH_GS_DIR check — expected, and ignored here: this
#          stage's success criterion is "./run exists and is newer than run.cpp".
#   smoke  20 real steps + the t=0 analytic gates. NEEDS the ground state.
#   prod   the full run.
# The first submission collapsed 'build' and 'smoke' into one stage that ran
# before the ground state existed; the script's own GS guard then killed it in 5
# seconds WITHOUT COMPILING, defeating the early-compile-check it was supposed to
# provide, and its afterok cancelled both production halves (2026-08-01).

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

SYS="$REPO_ROOT/ResearchProject/systems/cylindrical_jellium"
RUN_DIR="$SYS/scripts/channeling_twin/$HALF"
GS_DIR="$SYS/shared_gs/tube_rs3_R10_14_L60_dx0p5"

if [ "$STAGE" = "smoke" ] || [ "$STAGE" = "build" ]; then
  # CKPT > NSTEPS so no INTERIOR checkpoint is written. A 20-step gate run is
  # not something anyone will resume, and each checkpoint of 104 states on this
  # grid is 1.2 GB — the first version left 5.1 GB of throwaway checkpoints
  # across the two smokes, which was 30 % of the whole study's footprint on a
  # filesystem with 25 GB free (2026-08-01). The FINAL checkpoint is still
  # written by the binary (.claude/rules/final-timestep-checkpoint.md) and is
  # pruned below once the gates have passed.
  NAME="smoke"; NSTEPS=20; SAVE=10; WF=0; CKPT=999999
else
  # CADENCES SIZED TO THE DISK, NOT TO TASTE (measured 2026-08-01).
  # A VTI frame on this 80x80x120 grid costs 8.7 MB and an RT checkpoint of 104
  # states costs 1.28 GB. At the original SAVE_EVERY=5 / WF=75 / MAX_CKPT=3 the
  # pair needs 24.9 GB against 25 GB free — it would have filled /rds mid-run and
  # killed both halves, which is exactly what happened to 11 of 16 runs on
  # 2026-07-31. At 15/150/1 it needs 11.0 GB.
  # Nothing scientific is lost: every scalar observable is still written EVERY
  # step, and 101 frames is well above the 30 the GIF builders actually sample.
  NAME="${CH_OUT:-$HALF}"; NSTEPS="${CH_N_STEPS:-1500}"
  SAVE="${CH_SAVE_EVERY:-15}"; WF="${CH_WF_EVERY:-150}"; CKPT="${CH_CKPT_EVERY:-500}"
fi
export CH_MAX_CKPT="${CH_MAX_CKPT:-1}"   # rolling `checkpoint` + 1 numbered snapshot

# PRE-FLIGHT DISK GUARD. The notebook job has always had one; the RUN jobs did
# not, so a full filesystem was discovered only after the GPU hours were spent.
# ~7 GB is the larger half's own footprint plus headroom for its checkpoints.
NEED_GB="${CH_NEED_GB:-9}"
hd=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); printf "%.0f", $3-$2; exit}')
if [ -n "${hd:-}" ] && [ "$STAGE" != "smoke" ] && [ "$hd" -lt "$NEED_GB" ]; then
  echo "FATAL: only ${hd} GB free on /rds; this half needs ~${NEED_GB} GB." >&2
  echo "       Free space, or lower CH_SAVE_EVERY/CH_MAX_CKPT, before resubmitting." >&2
  exit 3
fi
echo " disk headroom: ${hd:-unknown} GB (need ~${NEED_GB} GB)"

if [ "$STAGE" != "build" ] && [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR — run-chan-gs.slurm must succeed first." >&2
  exit 2
fi

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  chan-twin [$HALF/$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " half   : $HALF   stage: $STAGE"
echo " system : annular tube r_s=3, R_in=10 R_out=14 L_z=60, N=160"
echo " proj   : 50 eV (v=1.917, v/v_F=3.00), sigma_WP=4 (sigma_pot=2.8284), z0=-28"
echo " steps  : $NSTEPS x dt 0.02 = $(awk "BEGIN{printf \"%.1f\", $NSTEPS*0.02}") a.u."
echo "          path = $(awk "BEGIN{printf \"%.1f\", 1.91701127*$NSTEPS*0.02}") Bohr"
echo " cadence: density/$SAVE  wavefn/$WF  ckpt/$CKPT  scalars/1"
echo " resume : ${CH_RESUME:-0}"
echo " gs     : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export CH_OUT="$NAME" CH_N_STEPS="$NSTEPS" CH_SAVE_EVERY="$SAVE" \
       CH_WF_EVERY="$WF" CH_CKPT_EVERY="$CKPT" \
       CH_GS_DIR="$GS_DIR" CH_RESUME="${CH_RESUME:-0}"
# Everything else defaults to shared/configs/channeling_tube_rs3.hpp; overriding
# any of CH_N / CH_V0 / CH_SIGMA / CH_LAUNCH_Z / CH_DT / CH_SPACING breaks twin
# parity unless it is done for BOTH halves (check_twin.py will catch it).

if [ "$STAGE" = "build" ]; then
  # inq-run builds, then execs the binary, which exits 2 on the missing GS.
  # That runtime failure is EXPECTED here; the compile is what is being tested.
  inq-run || true
  if [ -x ./run ] && [ ./run -nt run.cpp ]; then
    echo " [build] ./run compiled successfully ($(stat -c%s ./run) bytes)"
    status=0
  else
    echo " [build] FAILED: ./run was not produced from run.cpp — see the log above." >&2
    status=1
  fi
elif [ "$STAGE" = "smoke" ]; then
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

# A PASSED smoke has served its purpose: keep the observables and the gate log as
# provenance, drop the multi-GB checkpoints nobody will resume from.
if [ "$STAGE" = "smoke" ] && [ "$status" -eq 0 ]; then
  freed=$(du -sm results/"$NAME"/checkpoint results/"$NAME"/ckpt_step* 2>/dev/null | awk '{s+=$1} END{print s+0}')
  rm -rf results/"$NAME"/checkpoint results/"$NAME"/ckpt_step* 2>/dev/null
  echo " [cleanup] smoke passed — removed ${freed:-0} MB of throwaway checkpoints"
fi

echo "=============================================="
echo " [$HALF/$NAME] exit status : $status"
echo " finished                  : $(date)"
S="results/$NAME/run_summary.txt"
[ -f "$S" ] && { echo " --- run_summary.txt ---"; cat "$S"; }
for f in results/"$NAME"/raw/observables/*.csv; do
  [ -f "$f" ] || continue
  echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"
  head -2 "$f"; tail -1 "$f"
done
echo " --- retained checkpoints ---"
ls -d "results/$NAME"/ckpt_step* 2>/dev/null || echo "   (none yet)"
echo "=============================================="
exit $status
