#!/bin/bash
#SBATCH --job-name=chan-sic
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
# run-chan-sic.slurm — the SIC-corrected WP half of the channeling comparison.
# Plan: docs/plans/wp-self-interaction-correction.md (reviewed 2026-08-02).
#
#   sbatch shared/bin/run-chan-sic.slurm tierb    # build + 200-step bath gates
#   sbatch shared/bin/run-chan-sic.slurm prod     # the 1500-step production run
#   CH_RESUME=1 CH_N_STEPS=3000 sbatch shared/bin/run-chan-sic.slurm prod
#
# Chain (submit-channeling-sic.sh): tests -> tierV(vacuum) -> tierb -> prod.
# The classical reference and the uncorrected WP run are NOT re-run — they are
# the completed channeling_twin halves; this run is the third leg of the
# three-way comparison and shares their GS, grid, dt and cadences exactly.
#
# TIER B (200 steps of the production system) exists because the vacuum tier
# never exercises the Q projection (no occupied bath there — plan §0/D6). Its
# binary exits 4 on hard gate failure (cum_norm_removed, max_overlap_pre,
# E_corrected drift), which fails this job and afterok blocks production.
set -euo pipefail

STAGE="${1:-prod}"
case "$STAGE" in
  tierb|prod) ;;
  *) echo "FATAL: argument must be 'tierb' or 'prod'; got '$STAGE'" >&2; exit 2 ;;
esac

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
RUN_DIR="$SYS/scripts/channeling_sic/wp"
GS_DIR="$SYS/shared_gs/tube_rs3_R10_14_L60_dx0p5"

export CH_SIC="${CH_SIC:-pzrun}"

if [ "$STAGE" = "tierb" ]; then
  # No interior checkpoint (CKPT > NSTEPS): a 200-step gate run is never
  # resumed, and checkpoints here are 1.2 GB each. Final ckpt still written by
  # the binary and pruned below once the gates pass.
  NAME="tierb"; NSTEPS=200; SAVE=100; WF=0; CKPT=999999
  export CH_SIC_TIER=b
else
  # Cadences identical to the channeling twin production halves (disk-sized,
  # measured 2026-08-01): scalars every step, density/15, wavefn/150, 1 ckpt.
  NAME="${CH_OUT:-wp_sic}"; NSTEPS="${CH_N_STEPS:-1500}"
  SAVE="${CH_SAVE_EVERY:-15}"; WF="${CH_WF_EVERY:-150}"; CKPT="${CH_CKPT_EVERY:-500}"
  export CH_SIC_TIER=prod
fi
export CH_MAX_CKPT="${CH_MAX_CKPT:-1}"

NEED_GB="${CH_NEED_GB:-9}"
hd=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); printf "%.0f", $3-$2; exit}')
if [ -n "${hd:-}" ] && [ "$STAGE" = "prod" ] && [ "$hd" -lt "$NEED_GB" ]; then
  echo "FATAL: only ${hd} GB free on /rds; this run needs ~${NEED_GB} GB." >&2
  exit 3
fi
echo " disk headroom: ${hd:-unknown} GB (need ~${NEED_GB} GB)"

if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: ground state missing at $GS_DIR (the channeling_twin GS is reused)." >&2
  exit 2
fi

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  chan-sic [$STAGE/$NAME]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " sic    : $CH_SIC   tier: $CH_SIC_TIER"
echo " system : annular tube r_s=3, R_in=10 R_out=14 L_z=60, N=160 (twin GS)"
echo " proj   : 50 eV (v=1.917, v/v_F=3.00), sigma_WP=4, z0=-28"
echo " steps  : $NSTEPS x dt 0.02 = $(awk "BEGIN{printf \"%.1f\", $NSTEPS*0.02}") a.u."
echo " cadence: density/$SAVE  wavefn/$WF  ckpt/$CKPT  scalars/1"
echo " resume : ${CH_RESUME:-0}"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

export CH_OUT="$NAME" CH_N_STEPS="$NSTEPS" CH_SAVE_EVERY="$SAVE" \
       CH_WF_EVERY="$WF" CH_CKPT_EVERY="$CKPT" \
       CH_GS_DIR="$GS_DIR" CH_RESUME="${CH_RESUME:-0}"

if [ "$STAGE" = "tierb" ]; then
  inq-run                       # builds the binary, then exec's it
  status=$?
else
  if [ ! -x ./run ]; then
    echo "FATAL: ./run not built. Submit the 'tierb' stage first (it builds)." >&2
    exit 2
  fi
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" ./run
  status=$?
fi

# A passed Tier B has served its purpose: keep observables + gate log, drop the
# multi-GB checkpoints nobody will resume from.
if [ "$STAGE" = "tierb" ] && [ "$status" -eq 0 ]; then
  freed=$(du -sm results/"$NAME"/checkpoint results/"$NAME"/ckpt_step* 2>/dev/null | awk '{s+=$1} END{print s+0}')
  rm -rf results/"$NAME"/checkpoint results/"$NAME"/ckpt_step* 2>/dev/null
  echo " [cleanup] tierb passed — removed ${freed:-0} MB of throwaway checkpoints"
fi

echo "=============================================="
echo " [$STAGE/$NAME] exit status : $status"
echo " finished                   : $(date)"
S="results/$NAME/run_summary.txt"
[ -f "$S" ] && { echo " --- run_summary.txt ---"; cat "$S"; }
for f in results/"$NAME"/raw/observables/*.csv; do
  [ -f "$f" ] || continue
  echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"
  head -2 "$f"; tail -1 "$f"
done
echo "=============================================="
exit $status
