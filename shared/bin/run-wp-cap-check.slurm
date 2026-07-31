#!/bin/bash
#SBATCH --job-name=wp-cap-check
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-wp-cap-check.slurm — CAP DEFINITION VALIDATION REPLICA.
#
# Free Gaussian wavepacket (sigma_WP = 0.5, launch z = -24) in an EMPTY box with
# the EXACT production geometry (35 x 35 x 85 Bohr, periodicity(2), dx = 0.5) and
# the EXACT production CAP (two sin^2 bands, 12.5 Bohr per z face, |eta| = 1 Ha).
# Proves the CAP is defined and placed correctly BEFORE six 1.5x-length
# production runs are launched.
#
#   Usage:  sbatch shared/bin/run-wp-cap-check.slurm
#
# Three variants, one build:
#   cap_v2p0    CAP on,  k0 = 2.0  (slowest production velocity, longest transit)
#   cap_v4p5    CAP on,  k0 = 4.5  (fastest -> weakest absorption, worst case:
#                                   analytic survival exp(-12.5/v) = 6.2e-2 per hump)
#   nocap_v2p0  CAP off, k0 = 2.0  (CONTROL: must show the packet wrapping the +z
#                                   face and re-entering at -z, which is what the
#                                   CAP exists to prevent)
#
# ENGINE: inq-study, NOT stock inq. run.cpp CONSTRUCTS perturbations::absorbing,
# whose `vk[...] += complex(0.0, ...)` does not compile against stock inq's REAL
# scalar potential. inq-study complexifies it (self_consistency.hpp:176). Setting
# INQ_SOURCE before sourcing config.sh also redirects INQ_DEPS_CACHE, so the share
# paths are pinned back to the stock install explicitly.
#
# Plan: docs/plans/wavepacket-highdensity-sv-twin.md
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/cap_check"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  wp-cap-check"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " rundir : $RUN_DIR"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

# ---- variant 1: CAP on at v = 2.0 (also performs the one-off build) ---------
echo ""
echo "### [1/3] cap_v2p0 — CAP on, k0 = 2.0 (builds the binary) ###"
WPC_OUT=cap_v2p0 WPC_ETA=-1.0 WPC_K0=2.0 WPC_NSTEPS=1200 \
  inq-run || { echo "FATAL: build/run of cap_v2p0 failed" >&2; exit 1; }

# inq-run exec's ./run after building, so the binary now exists and the remaining
# variants are plain invocations with different env.
run_variant() {
  local name="$1" eta="$2" k0="$3" nsteps="$4"
  echo ""
  echo "### $name — eta=$eta k0=$k0 nsteps=$nsteps ###"
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
      WPC_OUT="$name" WPC_ETA="$eta" WPC_K0="$k0" WPC_NSTEPS="$nsteps" \
      ./run || echo "WARN: $name failed (continuing so the others still ship)"
}

# ---- variant 2: worst-case absorption, fastest production velocity ----------
# k0=4.5 reaches the +z CAP inner edge (+30) at t=(30+24)/4.5=12 a.u. = step 300;
# 900 steps is ample to see full absorption and confirm no re-entry.
run_variant cap_v4p5 -1.0 4.5 900

# ---- variant 3: no-CAP CONTROL ---------------------------------------------
# k0=2.0 crosses the +z face (+42.5) at t=33.25 a.u. = step 831, so 1200 steps
# runs ~370 steps PAST the wrap — long enough for the re-entry to be unmistakable.
run_variant nocap_v2p0 0.0 2.0 1200

echo ""
echo "=============================================="
echo " finished : $(date)"
for d in cap_v2p0 cap_v4p5 nocap_v2p0; do
  f="results/$d/run_summary.txt"
  if [ -f "$f" ]; then
    echo " --- $d : run_completed = $(grep -c 'run_completed = true' "$f") ---"
    m="results/$d/raw/observables/wp_momentum_stats.csv"
    [ -f "$m" ] && echo "     wp_momentum_stats.csv : $(wc -l < "$m") rows"
  else
    echo " --- $d : MISSING run_summary.txt ---"
  fi
done
echo "=============================================="
