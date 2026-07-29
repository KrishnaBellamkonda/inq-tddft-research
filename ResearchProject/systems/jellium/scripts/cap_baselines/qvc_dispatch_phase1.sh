#!/usr/bin/env bash
# Autonomous dispatcher — Study A (WITH-CAP) quantum-vs-classical, sigma_wp=3, 300 eV.
# Phase 1: R4 baseline (b1), R1 WP (b3), R2 matched classical (b2). One GPU, sequential.
# Config (locked w/ user): 50^3 box, N=162, reuse GS; CAP L=20 (10/side) eta=-1.0;
#   launch z0=-6.5; V0=4.696 (300 eV); N_STEPS=336 (tau=6.72 au = rigid full-exit);
#   WRITE_EVERY=5. Matched classical UPF = electron_gaussian_wpsigma3p0.upf (charge std 2.121).
# Emails per run + a final all-done email. PROVISIONAL until Task #7.
set -u
cd /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/scripts/cap_baselines
export CUDA_VISIBLE_DEVICES=1
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
PY=/local/data/public/skcb2/tddft/venv/bin/python3
WPUPF=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma3p0.upf
COMMON="CAP_N_STEPS=336 CAP_WRITE_EVERY=5 CAP_ETA=-1.0 CAP_WIDTH_BOHR=10 CAP_LAUNCH_Z=-6.5 CAP_V0=4.696 CAP_WP_SIGMA=3.0"
FAM="qvc-cap-s3-E300"

run () {  # $1=mode $2=subdir $3=label $4=extra-env
  echo "[$(date '+%F %T')] START $3"
  env $COMMON CAP_MODE=$1 CAP_OUT_SUBDIR=$2 $4 ./run > logs_$2.log 2>&1
  local st=$?
  echo "[$(date '+%F %T')] DONE $3 exit=$st"
  $PY email_run.py "$2" "$3" "$st" "$FAM" || true
}

run b1 qvc_b1_s3      "Study A R4 baseline (CAP on, no projectile)"            ""
run b3 qvc_wp_s3_E300 "Study A R1 wavepacket (CAP, sigma_wp=3, 300 eV)"        ""
run b2 qvc_cl_s3_E300 "Study A R2 matched classical (CAP, charge std 2.121, 300 eV)" "CAP_PROJ_PSEUDO=$WPUPF"

$PY email_run.py qvc_cl_s3_E300 "Study A WITH-CAP — ALL 3 runs complete (R4,R1,R2)" 0 "$FAM" \
   "Next: vacuum-WP SIE control (R3) + no-CAP Study B pending user confirmation." || true
echo "[$(date '+%F %T')] Phase 1 complete."
