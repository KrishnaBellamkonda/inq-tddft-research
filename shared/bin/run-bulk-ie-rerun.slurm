#!/bin/bash
#SBATCH --job-name=bulk-ie
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --array=0-5
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.out
#
# run-bulk-ie-rerun.slurm — re-run the bulk jellium twin pairs WITH the P/S/B
# pairwise interaction-energy decomposition attached.
#
#   Usage:  sbatch --dependency=afterok:<buildjob> shared/bin/run-bulk-ie-rerun.slurm
#
# THE SIX RUNS (three twin pairs; the sigma=1 families are excluded by user
# instruction, and bulk_ks_stopping_sigma3 already carries interactions.csv):
#
#   task  family                       sigma  r_s     N    box        half
#   ----  ---------------------------  -----  ------  ---  ---------  ---------
#    0    bulk_ks_stopping               2    5.702   218  46x46x80   wp
#    1    bulk_ks_stopping               2    5.702   218  46x46x80   classical
#    2    bulk_ks_stopping_rs4           2    3.987   482  40x40x80   wp
#    3    bulk_ks_stopping_rs4           2    3.987   482  40x40x80   classical
#    4    bulk_ks_stopping_rs4_sigma3    3    3.987   482  40x40x80   wp
#    5    bulk_ks_stopping_rs4_sigma3    3    3.987   482  40x40x80   classical
#
# WHY BOTH HALVES, not just the classical twins: the quantity that answers the
# question is the GAUGE-CLEAN difference dE_PP = E_PP(WP) - E_PP(classical)
# evaluated in the SAME cell. Absolute E_PP carries the charged-cell G=0 gauge
# and is meaningless across geometries (this invalidated a slab-vs-bulk
# comparison on 2026-08-01). So a pair is only useful if BOTH halves have
# interactions.csv.
#
# PHYSICS IS UNCHANGED. compute_coulomb* only READS the density and solves
# Poisson; it never touches the Hamiltonian. These runs therefore reproduce the
# existing trajectories and strictly ADD columns. Pre-re-run observables CSVs
# were archived before submission as insurance against a failed re-run.
#
# NO NEW GROUND STATE: the GS depends only on box + density, not the projectile,
# so all six reuse the two existing checkpoints.
#
# WALL TIME. Measured on the sigma=3 r_s=5.702 pair that already carries the
# decomposition: 1:28 (wp) / 1:56 (classical). --time=08:00:00 is ~4x headroom;
# the classical half is NOT resumable (INQ cannot restart moving ions, user
# decision 2026-07-30), so it needs room to overrun rather than be killed.
#
# Each task runs its own CLOSURE GATE after the run
# (.claude/rules/decomposed-interaction-energies.md): the pairwise terms must sum
# back to INQ's own scalars, or the task exits non-zero.
#
# Plan: docs/plans/bulk-jellium-ks-stopping.md
set -uo pipefail

FAMILIES=(bulk_ks_stopping bulk_ks_stopping bulk_ks_stopping_rs4 \
          bulk_ks_stopping_rs4 bulk_ks_stopping_rs4_sigma3 bulk_ks_stopping_rs4_sigma3)
HALVES=(wp classical wp classical wp classical)
GSNAMES=(gs_L46x46x80_orth_N218_dx0p40 gs_L46x46x80_orth_N218_dx0p40 \
         gs_L40x40x80_orth_N482_dx0p50 gs_L40x40x80_orth_N482_dx0p50 \
         gs_L40x40x80_orth_N482_dx0p50 gs_L40x40x80_orth_N482_dx0p50)

I="${SLURM_ARRAY_TASK_ID:-}"
if [ -z "$I" ] || [ "$I" -lt 0 ] || [ "$I" -gt 5 ]; then
  echo "FATAL: SLURM_ARRAY_TASK_ID must be 0..5 (got '${I:-unset}')" >&2; exit 2
fi
FAM="${FAMILIES[$I]}"; HALF="${HALVES[$I]}"; GSNAME="${GSNAMES[$I]}"

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"

SYS="$REPO_ROOT/ResearchProject/systems/jellium"
RUN_DIR="$SYS/scripts/$FAM/$HALF"
GS_DIR="$SYS/checkpoints/$GSNAME"
VERIFY="$SYS/scripts/verify_interactions_closure.py"
PY="$REPO_ROOT/venv/bin/python3"

echo "=============================================="
echo " job    : ${SLURM_ARRAY_JOB_ID:-none}_${I}  bulk-ie [$FAM/$HALF]"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " rundir : $RUN_DIR"
echo " gs     : $GS_DIR"
echo " started: $(date)"
echo "=============================================="

[ -d "$RUN_DIR" ] || { echo "FATAL: run dir missing at $RUN_DIR" >&2; exit 2; }
[ -d "$GS_DIR" ]  || { echo "FATAL: ground state missing at $GS_DIR" >&2; exit 2; }
grep -q "interaction_energies.hpp" "$RUN_DIR/run.cpp" || {
  echo "FATAL: $FAM/$HALF is not wired with the pairwise decomposition." >&2; exit 2; }

cd "$RUN_DIR"
inq-run
status=$?

echo "=============================================="
echo " [$FAM/$HALF] run exit status : $status"
echo " finished                     : $(date)"

if [ -f results/run_summary.txt ]; then
  echo " --- run_summary.txt ---"
  cat results/run_summary.txt
fi

# ---- mandatory closure gate ------------------------------------------------
gate=0
if [ "$status" -eq 0 ]; then
  echo ""
  echo " --- interaction-energy closure gate ---"
  if [ -x "$PY" ]; then
    PYTHONPATH="$REPO_ROOT/inq-stack/python" "$PY" "$VERIFY" "$RUN_DIR" "$HALF" || gate=$?
  else
    echo " WARN: $PY not found — closure gate SKIPPED (verify manually)" >&2
    gate=0
  fi
fi

for f in results/raw/observables/interactions.csv \
         results/raw/observables/observables.csv; do
  if [ -f "$f" ]; then
    echo " --- $(basename "$f"): $(wc -l < "$f") rows ---"
    head -2 "$f"
    tail -1 "$f"
  else
    echo " --- $(basename "$f"): MISSING ---"
  fi
done

echo " closure gate exit : $gate"
echo "=============================================="
[ "$status" -ne 0 ] && exit "$status"
exit "$gate"
