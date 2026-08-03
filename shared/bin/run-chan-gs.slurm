#!/bin/bash
#SBATCH --job-name=chan-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-chan-gs.slurm — the SHARED ground state of the annular-tube channeling
# twin. Plan: docs/plans/cylindrical-channeling-ks-stopping.md
#
# Bare r_s = 3 tube: 40 x 40 x 60 Bohr fully periodic, dx = 0.5 (80x80x120),
# background between R_in = 10 and R_out = 14 with a 0.5-Bohr erfc edge,
# N = 160 electrons -> n0 = 8.84194e-3, r_s = 3.000000, 80 + 24 = 104 states.
#
# BOTH twins load THIS checkpoint. A different ground state under either half
# silently destroys the comparison, which is why there is one script and one
# directory rather than a per-run GS.
#
# IDEMPOTENT: if the checkpoint is already populated the job exits 0 immediately,
# so re-running the chain is safe.
#
#   Usage:  sbatch shared/bin/run-chan-gs.slurm
set -euo pipefail

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
RUN_DIR="$SYS/scripts/channeling_twin/gs"
GS_DIR="$SYS/shared_gs/tube_rs3_R10_14_L60_dx0p5"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  chan-gs"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " gs dir : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

# PRE-FLIGHT DISK GUARD. The ground state saves 104 states on a 768k-point grid
# = 1.3 GB, and the two RT halves that follow need ~11 GB between them. Checking
# here means a full filesystem costs seconds rather than an 8-hour SCF.
NEED_GB="${CH_NEED_GB:-4}"
hd=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); printf "%.0f", $3-$2; exit}')
echo " disk headroom: ${hd:-unknown} GB (need ~${NEED_GB} GB for the GS alone)"
if [ -n "${hd:-}" ] && [ "$hd" -lt "$NEED_GB" ]; then
  echo "FATAL: only ${hd} GB free on /rds — refusing to start an 8-hour SCF that cannot save." >&2
  exit 3
fi

if [ -d "$GS_DIR" ] && [ -n "$(ls -A "$GS_DIR" 2>/dev/null)" ]; then
  echo " [skip] ground state already present — reusing it (idempotent chain)."
  exit 0
fi

cd "$RUN_DIR"
export GS_DIR GS_TAG="gs" GS_N="${GS_N:-160}" GS_SPACING="${GS_SPACING:-0.5}"
inq-run
status=$?

echo "=============================================="
echo " exit status : $status"
echo " finished    : $(date)"
S="results/gs/run_summary.txt"
[ -f "$S" ] && { echo " --- run_summary.txt ---"; cat "$S"; }
echo "=============================================="
exit $status
