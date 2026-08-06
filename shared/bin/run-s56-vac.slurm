#!/bin/bash
#SBATCH --job-name=s56-vac
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
# run-s56-vac.slurm — the four CAP-only vacuum baselines for one sigma campaign
# of sigma56_sv. A free wavepacket in an EMPTY 35 x 35 x 105 box with the exact
# production CAP: no slab, no background, one electron. Any norm loss or <p_z>
# drift here is unambiguously the ABSORBER, not the physics.
#
#   Usage:  sbatch --export=ALL,LJ_SIGMA=6.0 shared/bin/run-s56-vac.slurm
#           sbatch --export=ALL,LJ_SIGMA=5.0 shared/bin/run-s56-vac.slurm
#
# EACH sigma NEEDS ITS OWN CONTROLS. The CAP's effect depends on how fast the
# packet spreads into the absorbing band, and that rate is 1/(sqrt2 sigma) —
# exactly what sigma sets. At these widths it is small (0.141 / 0.118 Bohr per
# a.u.), which is the point: these runs QUANTIFY the predicted 0.020 % (sigma=5)
# / 0.16 % (sigma=6) of the packet that starts inside the -z band at t = 0, so
# the slab runs can quote that one-off loss instead of absorbing it into S.
#
# Box, dx, launch z, dt and step counts are IDENTICAL to the slab runs, because a
# baseline is only subtractable step for step if the packet meets the absorber at
# the same moment. Density frames are ON at the production cadence (user, storage
# is not a constraint), so each control carries its own density GIF and the
# free-vs-slab packet can be compared frame for frame, not just scalar for scalar.
#
# ENGINE: inq-study — the CAP needs the complexified scalar potential.
# Builds its own binary on first use (separate run dir, separate CMake tree).
#
# Plan: docs/plans/sigma56-sv-twin.md
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

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/sigma56_sv/vac"
DX=0.40                       # must match run-s56-wp.slurm
LZ=105.0                      # must match the production box
SIGMA="${LJ_SIGMA:-6.0}"
TAG="s$(echo "$SIGMA" | awk '{printf "%.1f", $1}' | tr '.' 'p')_"
LAUNCH_Z="${LJ_LAUNCH_Z:--27.5}"   # must match the slab runs

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  s56-vac (CAP-only baselines)"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " box    : 35 x 35 x $LZ Bohr, dx = $DX  (matches production)"
echo " sigma  : $SIGMA Bohr   name prefix: 'vac_${TAG}'"
echo " launch : z = $LAUNCH_Z  (must match the slab runs)"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

#             v      N_steps, save cadence (identical to the slab runs)
V_LIST=(    2.0    2.5    3.0    3.5 )
N_LIST=(    4360   3488   2907   2491 )
SAVE_LIST=(   14     12     10      8 )

first=1
for i in 0 1 2 3; do
  V="${V_LIST[$i]}"; N="${N_LIST[$i]}"; SAVE="${SAVE_LIST[$i]}"
  NAME="vac_${TAG}v$(echo "$V" | tr '.' 'p')"
  echo ""
  echo "### $NAME — k0=$V  nsteps=$N  dx=$DX  Lz=$LZ ###"
  if [ "$first" = "1" ]; then
    # ALWAYS go through inq-run on the first variant, even when ./run exists.
    # Guarding on [ ! -x ./run ] silently reused a STALE binary once (job
    # 32418538): the source had gained a knob after the binary was built, so the
    # controls ran at hard-coded values and were unusable as baselines. cmake
    # rebuilds only what changed, so this costs nothing when already current.
    WPC_OUT="$NAME" WPC_K0="$V" WPC_NSTEPS="$N" WPC_H="$DX" WPC_LZ="$LZ" \
      WPC_ETA=-1.0 WPC_SIGMA="$SIGMA" WPC_LAUNCH_Z="$LAUNCH_Z" \
      WPC_SAVE_EVERY="$SAVE" inq-run || echo "WARN: $NAME failed (continuing)"
  else
    env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
        WPC_OUT="$NAME" WPC_K0="$V" WPC_NSTEPS="$N" WPC_H="$DX" WPC_LZ="$LZ" \
        WPC_ETA=-1.0 WPC_SIGMA="$SIGMA" WPC_LAUNCH_Z="$LAUNCH_Z" \
        WPC_SAVE_EVERY="$SAVE" ./run || echo "WARN: $NAME failed (continuing)"
  fi
  first=0
done

echo ""
echo "=============================================="
echo " finished : $(date)"
for i in 0 1 2 3; do
  NAME="vac_${TAG}v$(echo "${V_LIST[$i]}" | tr '.' 'p')"
  f="results/$NAME/raw/observables/wp_momentum_stats.csv"
  if [ -f "$f" ]; then echo " $NAME : $(wc -l < "$f") rows"; else echo " $NAME : MISSING"; fi
done
echo "=============================================="
