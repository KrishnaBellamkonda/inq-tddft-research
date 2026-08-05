#!/bin/bash
#SBATCH --job-name=lzb-vac
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
# run-lzb-vac.slurm — the four CAP-only vacuum baselines for ONE box preset of
# lz_bulk_sweep: a free wavepacket in the empty box with the exact production
# CAP. Any norm loss or <p_z> drift here is unambiguously the absorber.
# Plan: docs/plans/jellium-slab-extend-Lz.md
#
#   Usage:  sbatch --export=ALL,LZB_CFG=s5p0_L35 shared/bin/run-lzb-vac.slurm
#
# Box, dx, launch z, dt and step counts are IDENTICAL to the slab runs of the
# same preset (lzb-params.sh), because a baseline is only subtractable step for
# step if the packet meets the absorber at the same moment.
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

: "${LZB_CFG:?FATAL: pass the box via --export=ALL,LZB_CFG=<preset>}"
# shellcheck source=./lzb-params.sh
. "$REPO_ROOT/shared/bin/lzb-params.sh"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/scripts/lz_bulk_sweep/vac"
DX=0.40

echo "=============================================="
echo " job    : ${SLURM_JOB_ID:-none}  lzb-vac [$LZB_CFG] (CAP-only baselines)"
echo " node   : $(hostname)"
echo " gpu    : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine : $INQ_SOURCE"
echo " box    : 35 x 35 x $LZ Bohr, dx = $DX  (matches production)"
echo " sigma  : $SIGMA Bohr   launch z = $LAUNCH"
echo " started: $(date)"
echo "=============================================="

cd "$RUN_DIR"

V_TABLE=(2.0 2.5 3.0 3.5)
first=1
for i in 0 1 2 3; do
  VIDX="$i"
  # shellcheck source=./lzb-params.sh
  . "$REPO_ROOT/shared/bin/lzb-params.sh"
  NAME="vac_${LZB_CFG}_v$(echo "$V" | tr '.' 'p')"
  echo ""
  echo "### $NAME — k0=$V  nsteps=$NSTEPS  dx=$DX  Lz=$LZ ###"
  if [ "$first" = "1" ]; then
    # ALWAYS inq-run on the first variant — never trust a stale ./run.
    WPC_OUT="$NAME" WPC_K0="$V" WPC_NSTEPS="$NSTEPS" WPC_H="$DX" WPC_LZ="$LZ" \
      WPC_ETA=-1.0 WPC_SIGMA="$SIGMA" WPC_LAUNCH_Z="$LAUNCH" \
      WPC_SAVE_EVERY="$SAVE" inq-run || echo "WARN: $NAME failed (continuing)"
  else
    env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
        WPC_OUT="$NAME" WPC_K0="$V" WPC_NSTEPS="$NSTEPS" WPC_H="$DX" WPC_LZ="$LZ" \
        WPC_ETA=-1.0 WPC_SIGMA="$SIGMA" WPC_LAUNCH_Z="$LAUNCH" \
        WPC_SAVE_EVERY="$SAVE" ./run || echo "WARN: $NAME failed (continuing)"
  fi
  first=0
done

echo ""
echo "=============================================="
echo " finished : $(date)"
for i in 0 1 2 3; do
  NAME="vac_${LZB_CFG}_v$(echo "${V_TABLE[$i]}" | tr '.' 'p')"
  f="results/$NAME/raw/observables/wp_momentum_stats.csv"
  if [ -f "$f" ]; then echo " $NAME : $(wc -l < "$f") rows"; else echo " $NAME : MISSING"; fi
done
echo "=============================================="
