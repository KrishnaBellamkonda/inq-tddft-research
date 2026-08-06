#!/bin/bash
#SBATCH --job-name=s56-nb
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-s56-notebooks.slurm — post-processing for the sigma56_sv twin campaign:
# the S(v) summary table, the combined S(v) figure, and the per-run + synthesis
# notebooks (each carrying the mandatory density-matrix GIF,
# .claude/rules/notebook-density-gif.md).
#
#   Usage:  sbatch shared/bin/run-s56-notebooks.slurm            # everything
#           sbatch shared/bin/run-s56-notebooks.slurm summary    # tables + figure only
#
# The GPU is requested because the density-GIF battery is heavy on array work,
# not because anything here needs CUDA; it also keeps this job on the same node
# type as the runs so the filesystem view is identical.
#
# NOT gated on every run existing. Runs that are missing or incomplete are listed
# and SKIPPED, so one failed velocity still yields the rest. `complete` from
# s56_stopping.deposit_stopping() is checked per point — a still-propagating run
# once returned a perfectly plausible S from 86 of 3623 steps.
#
# Plan: docs/plans/sigma56-sv-twin.md
set -euo pipefail

MODE="${1:-all}"
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

HYP="$REPO_ROOT/ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  s56-notebooks [$MODE]"
echo " node   : $(hostname)"
echo " started: $(date)"
echo "=============================================="

# shellcheck source=/dev/null
. "$REPO_ROOT/venv/bin/activate"
cd "$HYP"

status=0

echo ""
echo "### S(v) summary tables + combined figure ###"
python build_sv_figure.py || { echo "WARN: figure stage failed"; status=1; }

echo ""
echo "### effective-width figure + per-run width table ###"
# <sigma_r> is averaged over t=0 -> 1 % norm loss and the legend quotes the MAX
# across velocities (the slowest run). Both live in s56_stopping.sigma_r_window
# so these two stages and the notebooks share one definition.
python build_sv_effective_width_s6.py || { echo "WARN: eff-width figure failed"; status=1; }
python make_width_table.py            || { echo "WARN: width table failed"; status=1; }
python make_sigma_r_traces.py         || { echo "WARN: sigma_r traces failed"; status=1; }
python make_sigma6_summary.py         || { echo "WARN: sigma6 summary failed"; status=1; }

if [ "$MODE" = "all" ]; then
  echo ""
  echo "### per-run + synthesis notebooks ###"
  python build_run_notebooks.py || { echo "WARN: notebook stage failed"; status=1; }
fi

echo ""
echo "=============================================="
echo " finished : $(date)"
echo " --- artefacts ---"
ls -la "$HYP"/*.csv "$HYP"/*.png "$HYP"/*.ipynb 2>/dev/null || echo "   (none)"
echo " exit status : $status"
echo "=============================================="
exit $status
