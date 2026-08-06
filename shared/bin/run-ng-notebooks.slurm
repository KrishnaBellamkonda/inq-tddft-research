#!/bin/bash
#SBATCH --job-name=ng-notebooks
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-ng-notebooks.slurm — figures, run notebooks and phase notebooks for the
# Nazarov-Gross mass ladder. Plan: docs/plans/nazarov-gross-slab-mass-ladder.md
#
# Pure post-processing: it reads results/ and writes hypotheses/ng_mass_ladder/.
# It never touches a run directory, so it is safe to re-run at any point and is
# the last phase the orchestrator calls.
#
# NOTEBOOKS ARE EXECUTED, NOT JUST WRITTEN. A notebook whose cells have not run
# is a script with markdown in it: the reader sees no numbers and no density GIF.
# nbconvert --execute is therefore the deliverable step, and a cell that raises
# is surfaced here rather than silently shipped.
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 1
fi
# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

PY="$REPO_ROOT/venv/bin/python3"
HYP="$REPO_ROOT/ResearchProject/systems/localised_jellium/hypotheses/ng_mass_ladder"
export PYTHONPATH="$REPO_ROOT/inq-stack/python:$HYP:${PYTHONPATH:-}"
export MPLBACKEND=Agg

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  ng-notebooks"
echo " node   : $(hostname)"
echo " target : $HYP"
echo " started: $(date)"
echo "=============================================="

cd "$HYP"
status=0

echo; echo "--- [1/3] figures ---"
"$PY" make_figures.py || { echo "figures FAILED"; status=1; }

echo; echo "--- [2/3] notebooks (write) ---"
"$PY" build_notebooks.py || { echo "notebook build FAILED"; status=1; }

echo; echo "--- [3/3] notebooks (execute) ---"
# Each notebook is executed independently so one bad run's notebook cannot take
# the whole set down. Timeout is generous: the density-GIF cell decodes tens of
# VTI frames.
for nb in phase_*.ipynb run_*.ipynb; do
  [ -f "$nb" ] || continue
  echo "  executing $nb"
  if "$PY" -m nbconvert --to notebook --execute --inplace \
       --ExecutePreprocessor.timeout=1800 --ExecutePreprocessor.kernel_name=python3 \
       "$nb" >/dev/null 2>"exec_$nb.log"; then
    echo "    OK"
  else
    echo "    FAILED — see $HYP/exec_$nb.log"
    tail -5 "exec_$nb.log" | sed 's/^/      /'
    status=1
  fi
done

echo; echo "--- inventory ---"
ls -1 figures/*.png 2>/dev/null | wc -l | xargs echo " figures:"
ls -1 *.ipynb        2>/dev/null | wc -l | xargs echo " notebooks:"
du -sh gifs 2>/dev/null || true
[ -f ng_ladder_table.csv ] && { echo " --- ng_ladder_table.csv ---"; cat ng_ladder_table.csv; }

echo "=============================================="
echo " exit status : $status"
echo " finished    : $(date)"
echo "=============================================="
# A failed notebook is NOT a reason to fail the campaign — the runs and their
# data are intact and the build is idempotent. Report, do not propagate.
exit 0
