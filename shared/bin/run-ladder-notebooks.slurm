#!/bin/bash
#SBATCH --job-name=ladder-nb
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
# run-ladder-notebooks.slurm — figures, per-rung notebooks and the cross-rung
# comparison for the cylindrical proximity ladder.
# Plan: docs/plans/cylindrical-proximity-ladder.md
#
# CHAINED WITH afterany, SO IT MUST DEGRADE RATHER THAN DIE. The production runs
# it depends on may have partly failed; build_ladder_figures.py discovers which
# rungs have BOTH halves complete, reports the rest, and writes up what exists.
# This job therefore exits 0 when it successfully analysed a subset — a non-zero
# exit here means the ANALYSIS broke, not that a rung was missing.
#
#   Usage:  cd <repo> && sbatch shared/bin/run-ladder-notebooks.slurm
#           LADDER_RUNGS="r10,r08" sbatch shared/bin/run-ladder-notebooks.slurm
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
[[ -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]] || {
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 1; }
# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export MPLBACKEND=Agg

# The venv carries inqview (matplotlib/VTK) — the analysis needs it, the runs do not.
if [ -f "$REPO_ROOT/venv/bin/activate" ]; then
  # shellcheck source=/dev/null
  . "$REPO_ROOT/venv/bin/activate"
fi

HYP="$REPO_ROOT/ResearchProject/systems/cylindrical_jellium/hypotheses/proximity_ladder"
RUNGS="${LADDER_RUNGS:-r10,r08,r06,r04,r00}"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  ladder-nb"
echo " node   : $(hostname)"
echo " python : $(python3 --version 2>&1)  ($(command -v python3))"
echo " rungs  : $RUNGS"
echo " started: $(date)"
echo "=============================================="

cd "$HYP" || exit 1
fails=0

echo
echo "########## figures + cross-rung comparison ##########"
python3 build_ladder_figures.py --rungs "$RUNGS" 2>&1
rc=$?
echo " build_ladder_figures exit=$rc"
[ "$rc" -ne 0 ] && fails=$((fails+1))

if [ -f build_ladder_notebooks.py ]; then
  echo
  echo "########## per-rung + comparison notebooks ##########"
  python3 build_ladder_notebooks.py --rungs "$RUNGS" 2>&1
  rc=$?
  echo " build_ladder_notebooks exit=$rc"
  [ "$rc" -ne 0 ] && fails=$((fails+1))
fi

echo
echo "=============================================="
echo " summary table:"
[ -f figures/ladder_summary.csv ] && cat figures/ladder_summary.csv
echo " failures : $fails"
echo " finished : $(date)"
echo "=============================================="
exit $fails
