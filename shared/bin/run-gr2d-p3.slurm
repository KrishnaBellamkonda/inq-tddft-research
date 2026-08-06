#!/bin/bash
#SBATCH --job-name=gr2d-p3
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
# run-gr2d-p3.slurm — Phase 3 WP production sweep, two-definitions campaign.
# Plan: docs/plans/real-material-stopping-comparison.md
#
#   sbatch --array=0-8 shared/bin/run-gr2d-p3.slurm
#
# Array map (WP only, hex-centre site, CAP verdict eta=1/W=20 baked in):
#   0-5: MONO  E = 15/25/50/100/200/300 eV
#   6-8: BI    E = 15/50/200 eV          (E25/100/300 done in Phase 2)
# E15 carries the DOCUMENTED ~3% low-k CAP residue (plan Stage D) — the
# Geelen-overlap point; analysis quotes the residue, does not hide it.
# Idempotent; resume via --export=ALL,LJ_RESUME=1 + same array index.
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
[[ -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]] || { echo "bad REPO_ROOT" >&2; exit 1; }
. "$REPO_ROOT/shared/bin/csd3-env.sh"
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
. "$REPO_ROOT/shared/config.sh"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"; mkdir -p "$TMPDIR"

SYS="$REPO_ROOT/ResearchProject/systems/graphene"
RUN_DIR="$SYS/scripts/twodef_sv/wp"
cd "$RUN_DIR"
[ -x ./run ] || { echo "FATAL: wp ./run missing" >&2; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:?run as --array=0-8}"
VARIANTS=(mono mono mono mono mono mono bi bi bi)
EVS=(15 25 50 100 200 300 15 50 200)
K0S=(1.0500 1.3555 1.9170 2.7110 3.8340 4.6955 1.0500 1.9170 3.8340)
V="${VARIANTS[$IDX]}"; EV="${EVS[$IDX]}"; K0="${K0S[$IDX]}"
OUT="${V}_E${EV}"
GS_DIR="$SYS/shared_gs/gs_3x2_${V}_cut50_Lz80"

echo "== gr2d-p3[$IDX] wp $V E=$EV eV k0=$K0  node $(hostname)  $(date) =="

if [ -f "results/$OUT/run_summary.txt" ] \
   && grep -q "run_completed = true" "results/$OUT/run_summary.txt" \
   && [ "${LJ_RESUME:-0}" != "1" ]; then
  echo " [skip] $OUT already complete"; exit 0
fi

env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
    GR_VARIANT="$V" LJ_SIGMA=2.0 LJ_K0="$K0" LJ_LAUNCH_Z=-12 LJ_DT=0.025 \
    LJ_CAP_ETA=-1.0 LJ_CAP_L=20 LJ_OUT="$OUT" LJ_GS_DIR="$GS_DIR" \
    LJ_RESUME="${LJ_RESUME:-0}" ./run
st=$?
echo "== gr2d-p3[$IDX] $OUT exit=$st $(date) =="
[ -f "results/$OUT/run_summary.txt" ] && grep -E "run =|run_completed|wall_time" "results/$OUT/run_summary.txt"
exit $st
