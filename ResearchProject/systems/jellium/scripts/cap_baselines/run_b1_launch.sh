#!/bin/bash
# args: gpu eta subdir label
gpu=$1; eta=$2; subdir=$3; label=$4
cd /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/scripts/cap_baselines
CUDA_VISIBLE_DEVICES=$gpu CAP_MODE=b1 CAP_N_STEPS=7000 CAP_WRITE_EVERY=23 \
  CAP_ETA=$eta CAP_OUT_SUBDIR=$subdir ./run > ${subdir}.log 2>&1
status=$?
/local/data/public/skcb2/tddft/venv/bin/python3 email_on_done.py "$subdir" "$label" "$status" >> ${subdir}.log 2>&1
