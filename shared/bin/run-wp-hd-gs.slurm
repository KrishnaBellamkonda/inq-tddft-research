#!/bin/bash
#SBATCH --job-name=wp-hd-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-wp-hd-gs.slurm — ground state for the WP twin of classical-highdensity-sv.
#
# 35 x 35 x 85 Bohr, periodicity(2), dx = 0.5, 25-Bohr slab (half-width 12.5,
# erfc edge 1.0), N = 100 (r_s = 4.18), extra_states 24, T = 0.00862 eV, LDA.
#
# The original GS lived on the decommissioned /local/data/public machine, so it
# is recomputed here. run.cpp ABORTS (exit 3) unless it reproduces
# E_GS = 207.18322156141 Ha and integral n dV = 100 — reproducing that number is
# what licenses the six classical S values as the reference curve. Downstream
# jobs therefore chain with afterok, not afterany.
#
# TWO ground states are produced (user decision 2026-07-30):
#   0.5  FIDELITY CHECK — must reproduce E_GS = 207.18322156141 Ha, proving the
#        CSD3 rebuild really is the classical campaign's system. run.cpp exits 3
#        if it does not, so downstream jobs chain with afterok, not afterany.
#   0.4  PRODUCTION grid for the WP runs. sigma_WP = 0.5 is ONE grid point per
#        sigma at dx = 0.5 and aliases 10.4 % of its z-momentum past Nyquist at
#        v = 4.5; at 0.40 that falls to <= 0.9 %. The E_GS gate does NOT apply
#        here (a finer grid legitimately shifts the energy) — run.cpp reports it.
#
#   Usage:  sbatch shared/bin/run-wp-hd-gs.slurm 0.5     # fidelity check
#           sbatch shared/bin/run-wp-hd-gs.slurm 0.4     # production GS
#
# Plan: docs/plans/wavepacket-highdensity-sv-twin.md
set -euo pipefail

DX="${1:-0.4}"
case "$DX" in
  0.5|0.50) DXTAG="dx0p5" ;;
  0.4|0.40) DXTAG="dx0p4" ;;
  *) echo "FATAL: pass a spacing of 0.5 (fidelity) or 0.4 (production); got '$DX'" >&2; exit 2 ;;
esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

# inq-study: the WP production binary needs it for the CAP, and the classical
# campaign recorded engine = inq-study, so the GS is built the same way.
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/gs"
GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x85_${DXTAG}_per2"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  wp-hd-gs"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
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
