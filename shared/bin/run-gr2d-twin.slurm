#!/bin/bash
#SBATCH --job-name=gr2d-twin
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
# run-gr2d-twin.slurm — Phase 2 classical+WP twin pairs, two-definitions campaign.
# Plan: docs/plans/real-material-stopping-comparison.md
#
#   sbatch --export=ALL,GR_CAP_ETA=-2.0,GR_CAP_L=16 --array=0-5 \
#          shared/bin/run-gr2d-twin.slurm
#
# Array map: 0/1/2 = CLASSICAL at E=25/100/300 eV; 3/4/5 = WP at same energies.
# k0 = sqrt(2E/27.211386): 1.3555 / 2.7110 / 4.6955.
# CAP (GR_CAP_ETA, GR_CAP_L) comes from the Phase 1 sigma=2.0 capscan verdict —
# REQUIRED, no default. Vacuum baselines = the capscan runs at the same k.
# Idempotent: a member whose run_summary says run_completed=true is skipped.
# Resume: --export=ALL,...,LJ_RESUME=1 re-invokes the same member.
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
[[ -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]] || { echo "bad REPO_ROOT" >&2; exit 1; }
. "$REPO_ROOT/shared/bin/csd3-env.sh"
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
. "$REPO_ROOT/shared/config.sh"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

: "${GR_CAP_ETA:?GR_CAP_ETA required (from capscan verdict, <0)}"
: "${GR_CAP_L:?GR_CAP_L required (from capscan verdict, Bohr/face)}"

SYS="$REPO_ROOT/ResearchProject/systems/graphene"
GS_DIR="$SYS/shared_gs/gs_3x2_bi_cut50_Lz80"
UPF="$SYS/shared/pseudopotentials/electron_gaussian_wpsigma2p0_He.upf"

IDX="${SLURM_ARRAY_TASK_ID:?run as --array=0-5}"
EVS=(25 100 300 25 100 300)
K0S=(1.3555 2.7110 4.6955 1.3555 2.7110 4.6955)
HALF=$([ "$IDX" -lt 3 ] && echo classical || echo wp)
EV="${EVS[$IDX]}"; K0="${K0S[$IDX]}"
OUT="E${EV}"
RUN_DIR="$SYS/scripts/twodef_sv/$HALF"
cd "$RUN_DIR"

echo "== gr2d-twin[$IDX] $HALF E=$EV eV k0=$K0 cap(eta=$GR_CAP_ETA,W=$GR_CAP_L) =="
echo "   node $(hostname)  job ${SLURM_JOB_ID:-none}  $(date)"

if [ -f "results/$OUT/run_summary.txt" ] \
   && grep -q "run_completed = true" "results/$OUT/run_summary.txt" \
   && [ "${LJ_RESUME:-0}" != "1" ]; then
  echo " [skip] $HALF/$OUT already complete"; exit 0
fi

[ -x ./run ] || { echo "FATAL: ./run missing — run the compile check first" >&2; exit 1; }

ENV_COMMON=(GR_VARIANT=bi LJ_SIGMA=2.0 LJ_K0="$K0" LJ_LAUNCH_Z=-12
            LJ_DT=0.025 LJ_CAP_ETA="$GR_CAP_ETA" LJ_CAP_L="$GR_CAP_L"
            LJ_OUT="$OUT" LJ_GS_DIR="$GS_DIR" LJ_RESUME="${LJ_RESUME:-0}")
EXTRA_ENV=()
[ "$HALF" = "classical" ] && EXTRA_ENV=(LJ_PROJ_UPF="$UPF")

env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
    "${ENV_COMMON[@]}" "${EXTRA_ENV[@]}" ./run
st=$?

echo "== gr2d-twin[$IDX] $HALF/$OUT exit=$st  $(date) =="
[ -f "results/$OUT/run_summary.txt" ] && grep -E "run =|run_completed|proj_final|wall_time" "results/$OUT/run_summary.txt"
exit $st
