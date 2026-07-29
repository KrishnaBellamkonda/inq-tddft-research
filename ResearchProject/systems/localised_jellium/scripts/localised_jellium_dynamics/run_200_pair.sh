#!/bin/bash
set -e
export TMPDIR=/local/data/public/skcb2/tddft/.buildtmp CUDA_VISIBLE_DEVICES=0
export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
GS=/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/campaign_autorun/runs/h2/gs_p2_lz120/checkpoint
BASE=/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/localised_jellium_dynamics
RUN=/local/data/public/skcb2/tddft/shared/bin/inq-run
C="LJ_GS_DIR=$GS LJ_LAUNCH_Z=-24.5 LJ_SIGMA=0.5 LJ_PERIODICITY=2 LJ_LZ=120 LJ_HALF=12.5 LJ_N=82 LJ_SPACING=0.5 LJ_K0=1.0 LJ_N_STEPS=200 LJ_DT=0.05"
echo "### classical 0->200 $(date)"; cd $BASE/proj_dyn
env $C LJ_OUT=pdyn_k1_200 LJ_MASS=1.0 LJ_DELTA=0.1 "$RUN" > run_200.log 2>&1
echo "### classical done $(grep run_completed results/pdyn_k1_200/run_summary.txt)"
echo "### WP 0->200 $(date)"; cd $BASE/phase5_wp
env $C LJ_OUT=wp_k1_200 LJ_SAVE_EVERY=25 "$RUN" > run_200.log 2>&1
echo "### WP done $(grep run_completed results/wp_k1_200/run_summary.txt)  $(date)"
