#!/bin/bash
# args: gpu mode subdir label   (mode = b2|b3)
gpu=$1; mode=$2; subdir=$3; label=$4
cd /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/scripts/cap_baselines
CUDA_VISIBLE_DEVICES=$gpu CAP_MODE=$mode CAP_N_STEPS=7000 CAP_WRITE_EVERY=23 \
  CAP_ETA=-0.5 CAP_OUT_SUBDIR=$subdir ./run > ${subdir}.log 2>&1
status=$?
/local/data/public/skcb2/tddft/venv/bin/python3 email_on_done.py "$subdir" "$label" "$status" >> ${subdir}.log 2>&1
