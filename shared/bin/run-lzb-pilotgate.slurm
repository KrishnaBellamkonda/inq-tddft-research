#!/bin/bash
#SBATCH --job-name=lzb-gate
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-lzb-pilotgate.slurm — the AUTOMATED CHECK between the pilot (v = 3.0 at
# all four boxes, both halves) and the ~70 GPU-h of production velocities.
# Plan: docs/plans/jellium-slab-extend-Lz.md (user, 2026-08-05: "do one or two
# runs of velocity for all the Lz. Check if everything is alright before
# committing the massive number of GPU hours.")
#
# HARD gates (exit non-zero -> the production arrays' afterok NEVER fires):
#   completeness, ledger closure, finite/positive S. Correctness only.
# WARN (reported, never gates): GS bulk-likeness n(z=0) vs n0, S(L)
# monotonicity, cost re-projection from the measured s/step
# (checkpoint-dont-block: cost belongs to the user, correctness to the gate).
#
# Writes hypotheses/lz_bulk_sweep/PILOT_REPORT.md on EVERY path and attempts an
# email (best-effort). A gate failure therefore never fails silently: the
# production jobs sit DependencyNeverSatisfied and the report says why.
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
export PATH="$REPO_ROOT/shared/bin:$PATH"

HYP="$REPO_ROOT/ResearchProject/systems/localised_jellium/hypotheses/lz_bulk_sweep"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  lzb-pilotgate"
echo " node   : $(hostname)"
echo " started: $(date)"
echo "=============================================="

# shellcheck source=/dev/null
. "$REPO_ROOT/venv/bin/activate"
cd "$HYP"

python pilot_gate.py
status=$?

echo ""
echo "=============================================="
echo " finished : $(date)"
echo " exit     : $status  (0 = production released, non-zero = production BLOCKED)"
if [ -f "$HYP/PILOT_REPORT.md" ]; then
  echo " --- PILOT_REPORT.md ---"
  cat "$HYP/PILOT_REPORT.md"
fi
echo "=============================================="
exit $status
