#!/bin/bash
#SBATCH --job-name=wp-hd-nb
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-wp-hd-notebooks.slurm — build and EXECUTE the run notebooks once the sweep
# finishes, so the campaign is autonomous end to end (user requirement
# 2026-07-30: "Once the simulations are done, run notebooks for each of them must
# be produced").
#
# Produces, in hypotheses/wp_highdensity_sv/:
#   run_v2p0.ipynb  run_v2p5.ipynb  run_v3p0.ipynb  run_v3p5.ipynb
#   synthesis.ipynb + wp_S_summary.csv + wp_vs_classical_Sv.png
#
# Each run notebook opens with the density-matrix GIFs displayed inline
# (.claude/rules/notebook-density-gif.md), then derives S step by step, shows
# every KE and position term separately, the CAP/vacuum-control diagnostics, the
# full energy decomposition with its ledger-closure gates, and the comparison
# against the classical benchmark.
#
#   Usage:  sbatch --dependency=afterok:<sweep>,<vac> shared/bin/run-wp-hd-notebooks.slurm
#           sbatch shared/bin/run-wp-hd-notebooks.slurm 3.0    # single velocity
#
# A GPU is requested only because inqview pulls CUDA-linked VTK/numpy paths; the
# work itself is post-processing.
#
# Plan: docs/plans/wavepacket-highdensity-sv-twin.md
set -uo pipefail

TARGETS=("${@:-all}")        # one or more: all | synthesis | a velocity
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
HYP="$REPO_ROOT/ResearchProject/systems/localised_jellium/hypotheses/wp_highdensity_sv"
PY="$REPO_ROOT/venv/bin/python3"

# matplotlib must not try to write a font cache into a read-only or shared HOME
export MPLCONFIGDIR="${MPLCONFIGDIR:-/rds/user/skcb2/hpc-work/tmp/mplcache}"
mkdir -p "$MPLCONFIGDIR"
export PYTHONPATH="$REPO_ROOT/inq-stack/python:${PYTHONPATH:-}"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  wp-hd-notebooks [${TARGETS[*]}]"
echo " node   : $(hostname)"
echo " hyp    : $HYP"
echo " python : $PY"
echo " started: $(date)"
echo "=============================================="

cd "$HYP"

# The CAP-validation figure is cheap and independent of the sweep — refresh it so
# the campaign's evidence is all regenerated from the same code path.
"$PY" cap_check/build_cap_check_figure.py || echo "WARN: cap_check figure failed (continuing)"

status=0
for T in "${TARGETS[@]}"; do
  echo ""; echo "--- building notebook target: $T ---"
  "$PY" build_run_notebooks.py "$T" || status=$?
done

echo "=============================================="
echo " exit status : $status  (1 = some cells errored; notebooks still written)"
echo " finished    : $(date)"
ls -la "$HYP"/*.ipynb "$HYP"/*.csv "$HYP"/*.png 2>/dev/null
echo "=============================================="
exit $status
