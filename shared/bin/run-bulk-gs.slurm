#!/bin/bash
#SBATCH --job-name=bulk-gs
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-bulk-gs.slurm — ground state for the bulk-jellium KS-stopping twin pair.
#
#   46 x 46 x 80 Bohr orthorhombic periodic, N = 218 electrons, dx = 0.40, LDA.
#   Plan: docs/plans/bulk-jellium-ks-stopping.md
#
# Builds against STOCK inq (no CAP anywhere in this study, so inq-study is not
# needed — contrast run-dispersion.slurm, whose run.cpp always instantiates
# perturbations::absorbing and therefore requires the complexified fork).
#
# The GS run.cpp self-validates against the analytic plane-wave result
# (E_kinetic = 7.395907 Ha, E_hartree = 0, gap 0.244 eV, occupations 2.0/0.0)
# and EXITS 3 if any gate fails, so a dependent production job chained with
# afterok will not start on a bad ground state.
#
# Usage:  sbatch shared/bin/run-bulk-gs.slurm
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/jellium/save_gs/gs_L46x46x80_orth_N218_dx0p40"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  bulk-jellium GS"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " rundir : $RUN_DIR"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"
inq-run

status=$?
echo "=============================================="
echo " GS exit status : $status"
echo " finished       : $(date)"
if [ -f results/run_summary.txt ]; then
  echo " --- run_summary.txt ---"
  cat results/run_summary.txt
fi
echo "=============================================="
exit $status
