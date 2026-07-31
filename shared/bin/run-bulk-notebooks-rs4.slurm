#!/bin/bash
#SBATCH --job-name=bulk-nb-rs4
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-bulk-notebooks-rs4.slurm — autonomous final stage of the r_s = 3.99 chain.
#
# Builds and EXECUTES both run notebooks once the two production halves have
# finished. Chained with afterok on BOTH halves:
#
#   sbatch --dependency=afterok:$wp:$cl shared/bin/run-bulk-notebooks-rs4.slurm
#
# Runs on ampere with the project GPU account, matching every other script here
# (run-wp-hd-notebooks.slurm does the same). The stage itself is matplotlib +
# nbclient and does not use the card; there is no CPU allocation on this account.
#
# Plan: docs/plans/bulk-jellium-ks-stopping.md
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
HYP="$REPO_ROOT/ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping_rs4"
PY="$REPO_ROOT/venv/bin/python3"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  bulk-ks-stopping-rs4 notebooks"
echo " node   : $(hostname)"
echo " started: $(date)"
echo "=============================================="

cd "$HYP"
rc=0
for half in wp classical; do
  echo "--- building notebook: $half ---"
  if PYTHONPATH="$REPO_ROOT/inq-stack/python" "$PY" build_run_notebook.py "$half"; then
    echo "=== $half notebook OK ==="
  else
    echo "=== $half notebook FAILED (exit $?) ===" >&2
    rc=1
  fi
done

echo "=============================================="
echo " notebooks exit : $rc"
echo " finished       : $(date)"
ls -la "$HYP"/*.ipynb 2>/dev/null || true
echo "=============================================="
exit $rc
