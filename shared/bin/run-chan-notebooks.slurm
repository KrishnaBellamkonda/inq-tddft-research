#!/bin/bash
#SBATCH --job-name=chan-nb
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
# run-chan-notebooks.slurm — the whole post-processing stage of the annular-tube
# channeling twin, in dependency order.
# Plan: docs/plans/cylindrical-channeling-ks-stopping.md
#
#   1. check_twin.py     GATE. Verifies the pair really is a twin (same GS, cell,
#                        N, spacing, sigma, launch_z; projectile representation
#                        the ONLY difference) and that both halves carry the full
#                        energy decomposition. Writes twin_manifest.json. A
#                        failing pair is NOT analysed — a silent parity break
#                        would show up as a bogus quantum effect.
#   2. build_run_notebooks.py       the two per-run deep dives (classical + WP),
#                                   each opening with the density-GIF battery
#                                   (.claude/rules/notebook-density-gif.md).
#   3. build_comparison_notebook.py THE PHASE NOTEBOOK: classical vs WP side by
#                                   side, the four KS-stopping definitions
#                                   against the classical S, the var(p)-freeze
#                                   and f_bore channeling diagnostics, the
#                                   pairwise energy ledger with E_PP (the WP
#                                   self-Hartree, which has no classical
#                                   counterpart), and the explicit AIM-MET
#                                   verdict.
#
# Stage 1 is a hard gate; 2 and 3 run only if it passes. Stage 3 is the
# deliverable, so it runs even if a per-run notebook fails (they are independent).
#
#   Usage:  sbatch shared/bin/run-chan-notebooks.slurm
#           sbatch --dependency=afterany:<wp>:<cl> shared/bin/run-chan-notebooks.slurm
#           sbatch shared/bin/run-chan-notebooks.slurm --skip-gate   # re-render only
#
# A GPU is requested only because inqview pulls CUDA-linked VTK/numpy paths; the
# work itself is post-processing.
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
SYS="$REPO_ROOT/ResearchProject/systems/cylindrical_jellium"
HYP="$SYS/hypotheses/channeling_twin"
SCR="$SYS/scripts/channeling_twin"
PY="$REPO_ROOT/venv/bin/python3"
SKIP_GATE=0
[ "${1:-}" = "--skip-gate" ] && { SKIP_GATE=1; shift; }

export MPLCONFIGDIR="${MPLCONFIGDIR:-/rds/user/skcb2/hpc-work/tmp/mplcache}"
mkdir -p "$MPLCONFIGDIR"
export PYTHONPATH="$REPO_ROOT/inq-stack/python:${PYTHONPATH:-}"

echo "=============================================="
echo " job     : ${SLURM_JOB_ID:-none}  chan-nb"
echo " node    : $(hostname)"
echo " hyp dir : $HYP"
echo " started : $(date)"
echo "=============================================="

# Refuse to build into a full filesystem — GIFs are not small.
#
# THRESHOLD SIZED TO THIS STUDY (2026-08-01). The 30 GB inherited from
# run-slab-ks-notebooks.slurm is for a SIXTEEN-run campaign on a 1.65M-point
# grid; it blocked this two-run study at 26 GB free even though its notebooks
# need well under 1 GB. Two runs x (a ~30-frame GIF battery on 80x80 slices +
# ~20 figures) is a few hundred MB, so 3 GB is generous and 30 GB was simply the
# wrong number carried across. Raise it again if the run count grows.
NEED_GB="${CH_NB_NEED_GB:-3}"
hd=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); printf "%.0f", $3-$2; exit}')
echo " disk headroom: ${hd:-unknown} GB (need ~${NEED_GB} GB)"
if [ -n "${hd:-}" ] && [ "$hd" -lt "$NEED_GB" ]; then
  echo "FATAL: under ${NEED_GB} GB free — free space before building notebooks." >&2
  exit 3
fi

rc=0

if [ "$SKIP_GATE" -eq 0 ]; then
  echo
  echo "########## 1. twin parity gate ##########"
  mkdir -p "$SCR/pair"
  # --dynamic: the classical half is a MOVING projectile, so U_proj_bg is not a
  # constant and the per-step trajectory columns (proj_z, energy_proj_ke,
  # energy_proj_bg_ideal) are required, not optional.
  "$PY" "$REPO_ROOT/.claude/skills/twin-run-generation/check_twin.py" \
      --wp        "$SCR/wp/results/wp" \
      --classical "$SCR/classical/results/classical" \
      --pair-dir  "$SCR/pair" \
      --dynamic
  gate=$?
  echo " check_twin exit: $gate"
  if [ "$gate" -ne 0 ]; then
    echo "FATAL: the pair is NOT a valid twin — refusing to analyse it." >&2
    echo "       See $SCR/pair/twin_manifest.json for the failing fields." >&2
    exit "$gate"
  fi
else
  echo " [skip] twin parity gate (--skip-gate)"
fi

cd "$HYP" || exit 1

echo
echo "########## 2. per-run notebooks ##########"
"$PY" build_run_notebooks.py "$@" || { echo "  per-run notebooks FAILED"; rc=1; }

echo
echo "########## 3. comparison (phase) notebook ##########"
"$PY" build_comparison_notebook.py || { echo "  comparison notebook FAILED"; rc=1; }

echo "=============================================="
echo " exit status : $rc"
echo " finished    : $(date)"
echo " --- notebooks present ---"
ls -la "$HYP"/*.ipynb 2>/dev/null || echo "   (none)"
echo " --- summary table ---"
[ -f "$HYP/stopping_summary.csv" ] && cat "$HYP/stopping_summary.csv"
echo "=============================================="
exit $rc
