#!/usr/bin/env bash
set -e
cd /local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/effmass_12h/quantum
export CUDA_VISIBLE_DEVICES=0
export EM_GS_DIR=/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/shared_gs/slab_n42_L36x36x80_dx0p40
export EM_DT=0.05 EM_WRITE_EVERY=1 EM_CAP=1
echo "### SEG1: fresh 0->3 (ckpt final at 3) $(date)"
EM_OUT=smoke_seg EM_N_STEPS=3 EM_CKPT_EVERY=100 EM_RESUME=0 \
  EM_RT_CKPT_DIR=results/smoke_seg/rt_ckpt ./run 2>&1 | grep -E "step|RESUME|FRESH|ckpt|error|nan|inf" | tail -8
echo "### SEG2: resume 3->6 $(date)"
EM_OUT=smoke_seg EM_N_STEPS=6 EM_RESUME=1 \
  EM_RT_CKPT_DIR=results/smoke_seg/rt_ckpt ./run 2>&1 | grep -E "step|RESUME|start_step|ckpt|error|nan|inf" | tail -8
echo "### REF: fresh 0->6 continuous $(date)"
EM_OUT=smoke_ref EM_N_STEPS=6 EM_CKPT_EVERY=100 EM_RESUME=0 \
  EM_RT_CKPT_DIR=results/smoke_ref/rt_ckpt ./run 2>&1 | grep -E "step|FRESH|error|nan|inf" | tail -8
echo "### SMOKE DONE $(date)"
