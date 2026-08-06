#!/bin/bash
#SBATCH --job-name=s56-gs
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
# run-s56-gs.slurm — ground state for the sigma_WP = 5 and 6 twin campaign.
#
# 35 x 35 x 105 Bohr, periodicity(2), dx = 0.40, 25-Bohr slab (half 12.5, erfc
# edge 1.0), N = 100 (r_s = 4.183), extra_states 24, T = 0.00862 eV, LDA.
#
# ONE ground state serves ALL 18 production runs: sigma does not enter the bath,
# and both halves (wavepacket and classical) load the same checkpoint.
#
# The 85-Bohr campaigns' E_GS = 207.18322 Ha is REPORTED but NOT gated — a
# different box is a different calculation. What IS gated is the electron count
# and r_s, because those would only move if the background wiring were broken.
#
#   Usage:  sbatch shared/bin/run-s56-gs.slurm          # dx = 0.40 (production)
#           sbatch shared/bin/run-s56-gs.slurm 0.5      # coarse cross-check
#
# Plan: docs/plans/sigma56-sv-twin.md
set -euo pipefail

DX="${1:-0.4}"
case "$DX" in
  0.5|0.50) DXTAG="dx0p5" ;;
  0.4|0.40) DXTAG="dx0p4" ;;
  *) echo "FATAL: pass a spacing of 0.4 (production) or 0.5 (cross-check); got '$DX'" >&2; exit 2 ;;
esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

# inq-study: the RT binaries need it for the CAP, so the GS is built the same way.
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/sigma56_sv/gs"
GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x105_${DXTAG}_per2"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  s56-gs"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " box    : 35 x 35 x 105 Bohr, dx = $DX"
echo " gs dir : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"
export GS_DIR GS_SPACING="$DX" GS_TAG="$DXTAG"
inq-run
status=$?

echo "=============================================="
echo " exit status : $status"
echo " finished    : $(date)"
[ -f "results/$DXTAG/run_summary.txt" ] && { echo " --- run_summary.txt ---"; cat "results/$DXTAG/run_summary.txt"; }
echo "=============================================="
exit $status
