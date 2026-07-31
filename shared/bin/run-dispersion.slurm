#!/bin/bash
#SBATCH --job-name=wp-dispersion
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
# run-dispersion.slurm — free-Gaussian dispersion validation of the rebuilt engine.
#
# Runs the sigma_WP x energy vacuum wavepacket sweep (9 runs, no CAP, stock inq)
# and checks the density width against sigma_d(t) = sqrt(sigma0^2/2 + t^2/(2 sigma0^2)).
# See ResearchProject/systems/vacuum/scripts/wp_traversal_energy/dispatch_dispersion.py
#
# Usage:  sbatch shared/bin/run-dispersion.slurm
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

# MUST build against inq-study, not stock inq. wp_traversal_energy/run.cpp always
# CONSTRUCTS perturbations::absorbing and passes it to real_time::propagate, so the
# CAP code path is INSTANTIATED even when WP_ETA=0 makes it numerically inert.
# Stock inq keeps the scalar potential real (field<real_space,double>) and the CAP's
# `vk[...] += complex(0.0, ...)` fails to compile; inq-study complexifies it
# (self_consistency.hpp / ks_hamiltonian.hpp). Setting INQ_SOURCE before sourcing
# config.sh also redirects INQ_DEPS_CACHE and the share paths to that engine.
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"

echo "=============================================="
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " share  : $INQ_SHARE_PATH"
echo " jobs   : $INQ_BUILD_JOBS"
echo " started: $(date)"
echo "=============================================="

cd "$REPO_ROOT/ResearchProject/systems/vacuum/scripts/wp_traversal_energy"
"$INQ_PYTHON_EXE" --version 2>/dev/null || true

# python3.11 for the dispatcher (the repo venv is 3.6 and cannot run it).
/usr/bin/python3.11 dispatch_dispersion.py

echo "=============================================="
echo " SWEEP DONE — $(date)"
echo "=============================================="
