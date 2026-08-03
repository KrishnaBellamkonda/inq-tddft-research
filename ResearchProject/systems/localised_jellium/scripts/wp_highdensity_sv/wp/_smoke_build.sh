#!/usr/bin/env bash
set -e
export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
export PATH=/local/data/public/skcb2/tddft/shared/bin:$PATH
export CUDA_VISIBLE_DEVICES=0
export LJ_SIGMA=4 LJ_K0=3.8340 LJ_LAUNCH_Z=-24.5 LJ_SPACING=0.5 LJ_DT=0.04
export LJ_N_STEPS=5 LJ_CAP_ETA=-1.0 LJ_CAP_L=12.5 LJ_SAVE_EVERY=0 LJ_WF_EVERY=0 LJ_CKPT_EVERY=0
export LJ_GS_DIR=/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x85_dx0p5_per2
export LJ_OUT=smoke
echo "=== build+smoke start $(date) ==="
inq-run run.cpp
echo "=== smoke exit code $? $(date) ==="
