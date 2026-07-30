#!/usr/bin/env bash
# v=4.5 replica of the fastest previous sweep run (dyn/results/v4p5) but with the
# DIRECT-potential fix (direct erf/r perturbation + direct force + direct ledger).
# Exact old-v4p5 env; only the projectile representation differs. GPU 1 (GPU 0 busy).
set -euo pipefail
ROOT=/local/data/public/skcb2/tddft
SYS=$ROOT/ResearchProject/systems/localised_jellium
DIR=$SYS/scripts/classical_highdensity_sv/dyn_direct
cd "$DIR"
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod
export CUDA_VISIBLE_DEVICES=1
export LJ_LX=35 LJ_LY=35 LJ_LZ=85 LJ_HALF=12.5 LJ_N=100 LJ_EDGE_W=1.0
export LJ_PERIODICITY=2 LJ_SPACING=0.5 LJ_SIGMA=0.5 LJ_MASS=1.0 LJ_DELTA=0.1
export LJ_DT=0.04 LJ_CONST_V=0 LJ_LAUNCH_Z=-24.0 LJ_K0=4.5
export LJ_N_STEPS=1074 LJ_SAVE_EVERY=4
export LJ_GS_DIR=$SYS/shared_gs/slab_n100_L35x35x85_dx0p5_per2
export LJ_OUT=v4p5_direct LJ_RESUME="${LJ_RESUME:-0}"
echo "launch v4p5_direct on GPU $CUDA_VISIBLE_DEVICES  n_steps=$LJ_N_STEPS  $(date)"
exec ./run
