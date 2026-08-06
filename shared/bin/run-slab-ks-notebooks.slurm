#!/bin/bash
#SBATCH --job-name=slabks-nb
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
# run-slab-ks-notebooks.slurm — build and EXECUTE the single-run notebooks for
# the CAP-free wrap-around slab KS-stopping study.
# Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
# Skill: .claude/skills/run-notebook/ (the builder); this job is the dispatcher.
#
# Produces, in hypotheses/slab_ks_wrap/:
#   wp_n100_v2p0.ipynb … wp_n40_v3p5.ipynb          (8 wavepacket runs)
#   classical_n100_v2p0.ipynb … classical_n40_v3p5.ipynb   (8 wrapped twins)
#
# Each notebook opens with the density-matrix GIF battery displayed inline
# (.claude/rules/notebook-density-gif.md), then carries the full standard
# battery: energetics with the ΔE_total/N(t) 2-panel, the projectile/transport
# section, KS eigen-energy GIFs, collective response with the FFT-pipeline panel,
# momentum, KL metric, physical anchors and the loss function.
#
# ONLY COMPLETE RUNS ARE BUILT. build_run_notebooks.py gates on
# `run_completed = true` in run_summary.txt; a still-propagating run yields
# plausible-looking but WRONG numbers. Incomplete runs are listed and skipped, so
# re-running this job later picks them up — it is idempotent and safe to chain.
#
#   Usage:  sbatch shared/bin/run-slab-ks-notebooks.slurm
#           sbatch shared/bin/run-slab-ks-notebooks.slurm --only wp
#           sbatch --dependency=afterany:<wp_array>,<cl_array> shared/bin/run-slab-ks-notebooks.slurm
#
# A GPU is requested only because inqview pulls CUDA-linked VTK/numpy paths; the
# work itself is post-processing.
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
HYP="$REPO_ROOT/ResearchProject/systems/localised_jellium/hypotheses/slab_ks_wrap"
PY="$REPO_ROOT/venv/bin/python3"

# matplotlib must not write its font cache into a shared/read-only HOME
export MPLCONFIGDIR="${MPLCONFIGDIR:-/rds/user/skcb2/hpc-work/tmp/mplcache}"
mkdir -p "$MPLCONFIGDIR"
export PYTHONPATH="$REPO_ROOT/inq-stack/python:${PYTHONPATH:-}"

echo "=============================================="
echo " job     : ${SLURM_JOB_ID:-none}  slabks-nb"
echo " node    : $(hostname)"
echo " args    : ${*:-<all runs>}"
echo " hyp dir : $HYP"
echo " started : $(date)"
echo "=============================================="

# Refuse to build into a full filesystem — GIFs are not small.
hd=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); printf "%.0f", $3-$2; exit}')
echo " disk headroom: ${hd:-unknown} GB"
if [ -n "${hd:-}" ] && [ "$hd" -lt 30 ]; then
  echo "FATAL: under 30 GB free — free space before building notebooks." >&2
  exit 3
fi

cd "$HYP" || exit 1
"$PY" build_run_notebooks.py "$@"
status=$?

echo "=============================================="
echo " exit status : $status"
echo " finished    : $(date)"
echo " --- notebooks present ---"
ls -la "$HYP"/*.ipynb 2>/dev/null || echo "   (none)"
echo "=============================================="
exit $status
