#!/bin/bash
BIN=./run_pilot_validated
run_one() { # gpu eta
  CUDA_VISIBLE_DEVICES=$1 CAP_MODE=b1 CAP_N_STEPS=100 CAP_WRITE_EVERY=1 \
  CAP_ETA=$2 CAP_OUT_SUBDIR="eta_compare/eta_${2#-}" $BIN \
  > "eta_${2#-}.log" 2>&1
}
( run_one 0 -0.05; run_one 0 -0.20 ) &
( run_one 1 -0.10; run_one 1 -0.30 ) &
wait
echo "ETA_COMPARE_DONE"
