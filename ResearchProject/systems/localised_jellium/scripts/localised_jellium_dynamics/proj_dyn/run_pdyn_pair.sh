#!/bin/bash
# P-dyn twin pair (Rung 2): moving classical projectile (proj_dyn, fresh build) +
# WP twin (phase5_wp, k0 matched). One GPU (GPU0), sequential. 50 steps, dt=0.05.
set -e
export TMPDIR=/local/data/public/skcb2/tddft/.buildtmp
export CUDA_VISIBLE_DEVICES=0
export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
GS=/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/campaign_autorun/runs/h2/gs_p2_lz120/checkpoint
BASE=/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/localised_jellium_dynamics
RUN=/local/data/public/skcb2/tddft/shared/bin/inq-run
COMMON="LJ_GS_DIR=$GS LJ_LAUNCH_Z=-24.5 LJ_SIGMA=0.5 LJ_PERIODICITY=2 LJ_LZ=120 LJ_HALF=12.5 LJ_N=82 LJ_SPACING=0.5 LJ_N_STEPS=50 LJ_DT=0.05"

echo "### [1/2] classical proj_dyn (build + run) $(date)"
cd "$BASE/proj_dyn"
env $COMMON LJ_OUT=pdyn_k1 LJ_K0=1.0 LJ_MASS=1.0 LJ_DELTA=0.1 "$RUN" > build_run.log 2>&1
echo "### classical exit=$?  proj summary:"; grep -E "proj_z_final|run_completed" results/pdyn_k1/run_summary.txt 2>/dev/null || echo "  (no summary — check build_run.log)"

echo "### [2/2] WP twin phase5_wp k0=1.0 (run) $(date)"
cd "$BASE/phase5_wp"
env $COMMON LJ_OUT=wp_k1_dyn LJ_K0=1.0 LJ_SAVE_EVERY=25 "$RUN" > wp_run.log 2>&1
echo "### wp exit=$?  wp summary:"; grep -E "run_completed|wp_norm_after" results/wp_k1_dyn/run_summary.txt 2>/dev/null || echo "  (no summary — check wp_run.log)"
echo "### PAIR DONE $(date)"
