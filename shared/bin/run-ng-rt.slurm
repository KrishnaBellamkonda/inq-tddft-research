#!/bin/bash
#SBATCH --job-name=ng-rt
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
# run-ng-rt.slurm — ONE run of the Nazarov-Gross mass ladder (either half).
# Plan: docs/plans/nazarov-gross-slab-mass-ladder.md
#
#   Usage:  sbatch shared/bin/run-ng-rt.slurm wp build          # compile only
#           sbatch shared/bin/run-ng-rt.slurm classical build
#           NG_OUT=wp_m1 NG_MASS=1.0 NG_N_STEPS=2560 \
#             sbatch shared/bin/run-ng-rt.slurm wp
#
# EVERY PHYSICAL PARAMETER DEFAULTS TO THE CONFIG HEADER
# (shared/configs/slab_n206_L30x30x120_rs2p5.hpp). The orchestrator sets only
# what legitimately varies between rungs: NG_OUT, NG_MASS, NG_SIGMA_WP,
# NG_N_STEPS, NG_DT. Overriding anything else breaks twin parity.
#
# THE TIMESTEP IS MASS-DEPENDENT, and this is the single easiest thing to get
# wrong: dt <= 0.08 * M * h^2 (= 0.02*M at h=0.5), calibrated on p3 (M=1,
# dt=0.02) and sigma1_masspair (M=2, dt=0.04), both of which sat exactly on the
# ceiling. Cost therefore scales as 1/M — the LIGHT rungs are the expensive ones.
# If NG_DT is unset the binary computes the ceiling itself; if it is set too
# large the binary REFUSES to start. Do not "fix" that by raising the ceiling.
set -uo pipefail

HALF="${1:-}"
STAGE="${2:-prod}"
case "$HALF"  in wp|classical) ;; *) echo "FATAL: first arg must be 'wp' or 'classical'; got '$HALF'" >&2; exit 2 ;; esac
case "$STAGE" in build|prod)   ;; *) echo "FATAL: second arg must be 'build' or 'prod'; got '$STAGE'" >&2; exit 2 ;; esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 1
fi

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

SYS="$REPO_ROOT/ResearchProject/systems/localised_jellium"
RUN_DIR="$SYS/scripts/ng_mass_ladder/$HALF"
export NG_GS_DIR="${NG_GS_DIR:-$SYS/shared_gs/slab_n206_L30x30x120_rs2p5_dx0p5_per2}"

export NG_OUT="${NG_OUT:-$HALF}"
export NG_MASS="${NG_MASS:-1.0}"
export NG_N_STEPS="${NG_N_STEPS:-2560}"
export NG_RESUME="${NG_RESUME:-0}"
export NG_CAP="${NG_CAP:-1}"
# NG_SAVE_EVERY / NG_WF_EVERY / NG_CKPT_EVERY are DELIBERATELY not defaulted
# here: the binary derives them from N_STEPS (density/24, wavefn/4, ckpt/3) so
# every rung produces the same NUMBER of artefacts regardless of its step count.
# A fixed cadence would give the M=0.5 run (5120 steps at dt=0.01) six times the
# frames of the M=3 run (853 steps at dt=0.06) for identical physics.
# Measured artefact sizes on this grid: density VTI 6.9 MB, wavefunction VTI
# 13.8 MB, RT checkpoint 1.7 GB (rolling — one copy on disk, not one per write).
# Budget per run: ~390 MB frames (WP) or ~170 MB (classical) + 1.7 GB final ckpt.

# PRE-FLIGHT DISK GUARD. A full /rds killed 11 of 16 runs on 2026-07-31; the RT
# jobs had no guard until then, so the loss surfaced only after the GPU hours
# were spent. This box is small (864k points, 124 states) so the footprint is
# modest, but the guard is cheap and the failure mode is total.
NEED_GB="${NG_NEED_GB:-6}"
hd=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); printf "%.0f", $3-$2; exit}')
if [ -n "${hd:-}" ] && [ "$STAGE" = "prod" ] && [ "$hd" -lt "$NEED_GB" ]; then
  echo "FATAL: only ${hd} GB free on /rds; this run needs ~${NEED_GB} GB." >&2
  exit 3
fi

if [ "$STAGE" = "prod" ] && [ ! -d "$NG_GS_DIR" ]; then
  echo "FATAL: ground state missing at $NG_GS_DIR — run-ng-gs.slurm must succeed first." >&2
  exit 2
fi

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  ng-rt [$HALF/$NG_OUT] stage=$STAGE"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " system : slab r_s=2.5011 (v_F=0.767, E_F=8.01 eV), 15 Bohr thick, N=206"
echo " proj   : mass=$NG_MASS  sigma_WP=${NG_SIGMA_WP:-4.0 (config)}  v=${NG_V:-1.0742685 (config, 1.40 v_F)}"
echo " steps  : $NG_N_STEPS   dt=${NG_DT:-auto 0.08*M*h^2}   cap=$NG_CAP"
# These three are intentionally NOT exported above (the binary derives them from
# N_STEPS), so under `set -u` they must be referenced with defaults or the script
# dies before it ever reaches the compiler. That killed builds 32877571/2.
echo " cadence: density/${NG_SAVE_EVERY:-auto N/24}  wavefn/${NG_WF_EVERY:-auto N/4}" \
     " ckpt/${NG_CKPT_EVERY:-auto N/3}  scalars+interactions/1"
echo " resume : $NG_RESUME"
echo " gs     : $NG_GS_DIR"
echo " disk   : ${hd:-unknown} GB free (need ~${NEED_GB})"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

if [ "$STAGE" = "build" ]; then
  if [ -x ./run ] && [ ./run -nt run.cpp ]; then
    echo " [build] ./run is already newer than run.cpp — skipping the recompile"
    status=0
  else
    # COMPILE ONLY. `inq-run` builds and then EXECS, and there is no build-only
    # flag, so the ONLY way to stop the binary after compiling is to make its own
    # NG_GS_DIR guard fire. That used to happen for free because the ground state
    # did not exist yet; once it did, jobs 32881739/40 sailed past the guard and
    # started full 2560-step propagations inside a 45-minute build slot. Point the
    # guard at a path that cannot exist and the binary exits 2 immediately after
    # the compile — which is the outcome this stage is actually testing.
    NG_GS_DIR="/nonexistent/ng-build-compile-only" inq-run || true
    if [ -x ./run ] && [ ./run -nt run.cpp ]; then
      echo " [build] ./run compiled successfully ($(stat -c%s ./run) bytes)"; status=0
    else
      echo " [build] FAILED: ./run was not produced from run.cpp." >&2; status=1
    fi
  fi
else
  if [ -f "results/$NG_OUT/run_summary.txt" ] \
     && grep -q "run_completed = true" "results/$NG_OUT/run_summary.txt" \
     && [ "$NG_RESUME" = "0" ]; then
    echo " [skip] $NG_OUT already complete."
    status=0
  elif [ -x ./run ]; then
    env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run
    status=$?
  else
    echo "FATAL: ./run not built. Submit the 'build' stage first." >&2
    exit 2
  fi
fi

# PILOT CLEANUP. A passed pilot has served its purpose — nobody resumes a
# 600-step gate run — but its final checkpoint is 1.7 GB. Three pilots would
# hold 5 GB hostage for no scientific return. The observables, interactions and
# run_summary stay as provenance; only the checkpoint goes. Production runs keep
# theirs (rule final-timestep-checkpoint.md: every run must be extensible).
if [ "$STAGE" = "prod" ] && [ "$status" -eq 0 ] && [ "${NG_PRUNE_CKPT:-0}" = "1" ]; then
  freed=$(du -sm "results/$NG_OUT/rt_ckpt" "results/$NG_OUT/checkpoint" 2>/dev/null | awk '{s+=$1} END{print s+0}')
  rm -rf "results/$NG_OUT/rt_ckpt" "results/$NG_OUT/checkpoint" 2>/dev/null
  echo " [cleanup] NG_PRUNE_CKPT=1 — removed ${freed:-0} MB of pilot checkpoint"
fi

echo "=============================================="
echo " [$HALF/$NG_OUT] exit status : $status"
echo " finished                    : $(date)"
S="results/$NG_OUT/run_summary.txt"
[ -f "$S" ] && { echo " --- run_summary.txt ---"; cat "$S"; }
for f in results/"$NG_OUT"/raw/observables/*.csv; do
  [ -f "$f" ] || continue
  echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"
  head -2 "$f"; tail -1 "$f"
done
echo "=============================================="
exit $status
