#!/bin/bash
#SBATCH --job-name=wp-hd-vac
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
# run-wp-hd-vac.slurm — VACUUM CONTROLS for the WP twin sweep.
#
# WHY THESE EXIST. The cap_check replica (job 32417361) showed that in an EMPTY
# box with NO bath and NO forces, the CAP ALONE drags the surviving packet's
# <p_z> from 2.00 down to 0.61 over 48 a.u., while the CAP-off control holds it
# flat at 1.985. Because sigma_WP = 0.5 spreads at 1.414 Bohr/a.u., the packet's
# LEADING edge reaches the +z absorbing band first and is preferentially removed,
# so the mean momentum of what survives falls for purely non-physical reasons.
#
# Any S = -dT/ds fitted where the CAP is active would therefore be measuring the
# CAP, not the jellium. These runs give the CAP-ONLY baseline T(t), <p_z>(t) and
# norm(t) at each production velocity, on the SAME grid (dx = 0.40) and for the
# SAME number of steps, so the bath's true contribution is the DIFFERENCE between
# each slab run and its vacuum twin, step for step.
#
# They are cheap: 1 state instead of 74, non_interacting, so ~50x less work than
# a production point. All four run sequentially inside one job.
#
#   Usage:  sbatch shared/bin/run-wp-hd-vac.slurm
#           sbatch --export=ALL,LJ_SIGMA=2.0 shared/bin/run-wp-hd-vac.slurm
#
# LJ_SIGMA selects the campaign, matching run-wp-hd-wp.slurm: sigma = 0.5 keeps
# the bare names (vac_v2p0, already run), other sigmas get an s<sigma>_ prefix
# (vac_s2p0_v2p0). Each production campaign needs its OWN controls because the
# CAP's effect depends on how fast the packet spreads into the absorbing band,
# which is exactly what sigma sets: the spreading rate is 1/(sqrt2 sigma), so the
# sigma = 0.5 baseline (1.414 Bohr/a.u.) says nothing about sigma = 3 (0.236).
#
# Depends on the smoke stage only for the cap_check binary, which it builds
# itself on first use (separate run dir, separate CMake tree).
#
# Plan: docs/plans/wavepacket-highdensity-sv-twin.md (section 6b, item 5)
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
DX=0.40                      # production grid — must match run-wp-hd-wp.slurm
SIGMA="${LJ_SIGMA:-0.5}"     # campaign selector — must match run-wp-hd-wp.slurm
if [ "$(awk "BEGIN{print ($SIGMA == 0.5)}")" = "1" ]; then
  TAG=""
else
  TAG="s$(echo "$SIGMA" | awk '{printf "%.1f", $1}' | tr '.' 'p')_"
fi

# LAUNCH DISTANCE — must MATCH the slab runs these baselines are subtracted from.
# The control's whole job is to supply the CAP-only <p_z>(t) that is differenced
# out of the slab run step for step; if the packet starts somewhere else it
# reaches the +z CAP at a different time and the subtraction is meaningless.
# Default -24.0 leaves every existing vac_* name resolving unchanged.
# LJ_LAUNCH_Z=auto reads the scan's accepted z, so this can be chained behind
# the sweep with no human relaying the number.
LAUNCH_Z="${LJ_LAUNCH_Z:--24.0}"
if [ "$LAUNCH_Z" = "auto" ]; then
  SCAN_Z_FILE="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/inject_scan/results/scan/accepted_launch_z.txt"
  if [ ! -f "$SCAN_Z_FILE" ]; then
    echo "FATAL: LJ_LAUNCH_Z=auto but no accepted launch z at $SCAN_Z_FILE" >&2
    exit 2
  fi
  LAUNCH_Z="$(tr -d '[:space:]' < "$SCAN_Z_FILE")"
  echo "LJ_LAUNCH_Z=auto -> $LAUNCH_Z"
fi
if [ "$(awk "BEGIN{print ($LAUNCH_Z == -24.0)}")" = "1" ]; then
  NLTAG=""
else
  NLTAG="nl_"
fi
TAG="${NLTAG}${TAG}"

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  wp-hd-vac (CAP-only baselines)"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " dx     : $DX  (matches production)"
echo " sigma  : $SIGMA Bohr   name prefix: '${TAG}'"
echo " launch : z = $LAUNCH_Z  (must match the slab runs)"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

#          v      N_steps (identical to the slab runs)
V_LIST=( 2.0    2.5    3.0    3.5 )
N_LIST=( 3623   2898   2415   2070 )

first=1
for i in 0 1 2 3; do
  V="${V_LIST[$i]}"; N="${N_LIST[$i]}"
  NAME="vac_${TAG}v$(echo "$V" | tr '.' 'p')"
  echo ""
  echo "### $NAME — k0=$V  nsteps=$N  dx=$DX ###"
  if [ "$first" = "1" ]; then
    # ALWAYS go through inq-run on the first variant, even when ./run already
    # exists. Guarding on [ ! -x ./run ] silently reused a STALE binary in job
    # 32418538: cap_check/run.cpp had gained the WPC_H knob after that binary was
    # built, so all four controls ran at the HARDCODED dx = 0.50 instead of the
    # production 0.40. The t=0 gates caught it at v = 3.0 and 3.5 (sigma_pz^2
    # aliases +31 % at dx = 0.5) and aborted, but v = 2.0 and 2.5 completed on the
    # wrong grid and were unusable as CAP baselines. cmake rebuilds only what
    # changed, so this costs nothing when the source is already current.
    WPC_OUT="$NAME" WPC_K0="$V" WPC_NSTEPS="$N" WPC_H="$DX" WPC_ETA=-1.0 \
      WPC_SIGMA="$SIGMA" WPC_LAUNCH_Z="$LAUNCH_Z" \
      WPC_SAVE_EVERY=0 inq-run || echo "WARN: $NAME failed (continuing)"
  else
    env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
        WPC_OUT="$NAME" WPC_K0="$V" WPC_NSTEPS="$N" WPC_H="$DX" WPC_ETA=-1.0 \
      WPC_SIGMA="$SIGMA" WPC_LAUNCH_Z="$LAUNCH_Z" \
        WPC_SAVE_EVERY=0 ./run || echo "WARN: $NAME failed (continuing)"
  fi
  first=0
done

echo ""
echo "=============================================="
echo " finished : $(date)"
for i in 0 1 2 3; do
  NAME="vac_${TAG}v$(echo "${V_LIST[$i]}" | tr '.' 'p')"
  f="results/$NAME/raw/observables/wp_momentum_stats.csv"
  if [ -f "$f" ]; then
    echo " $NAME : $(wc -l < "$f") rows"
  else
    echo " $NAME : MISSING"
  fi
done
echo "=============================================="
