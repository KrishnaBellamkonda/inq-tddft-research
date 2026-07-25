#!/usr/bin/env bash
# Autonomous dispatcher — Phase 2: classical S(v) benchmark at the MATCHED width
# (charge std 2.121 = unified sigma_wp=3), PURE jellium, NO CAP (CAP_ETA=0), to compare
# the classical S(v) shape against analytical Lindhard. Auto-follows Phase 1.
# Energies 150/300/450/600 eV; launch z0=-16 (3 sigma_wp from box edge); ~1.5 box
# periods each; WRITE_EVERY=5. One GPU (1), sequential. Emails per run + all-done.
set -u
cd /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/scripts/cap_baselines
export CUDA_VISIBLE_DEVICES=1
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
PY=/local/data/public/skcb2/tddft/venv/bin/python3
WPUPF=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma3p0.upf
FAM="qvc-nocap-bench-s3"

# --- wait for Phase 1 (R2 classical) to complete ---
echo "[$(date '+%F %T')] Phase 2 watcher: waiting for Phase 1 completion..."
until [ -f results/qvc_cl_s3_E300/run_summary.txt ] && grep -q "run_completed  = true" results/qvc_cl_s3_E300/run_summary.txt 2>/dev/null; do
  sleep 60
done
echo "[$(date '+%F %T')] Phase 1 done -> starting Phase 2 benchmark."

run_E () {  # $1=E_eV
  local E=$1
  local V0=$($PY -c "import math;print(f'{math.sqrt(2*$E/27.2114):.5f}')")
  local NS=$($PY -c "import math;print(int(math.ceil(1.5*50/math.sqrt(2*$E/27.2114)/0.02)))")
  local SUB=qvc_bench_cl_s3_E$E
  echo "[$(date '+%F %T')] START benchmark E=$E eV  V0=$V0  N_STEPS=$NS"
  env CAP_MODE=b2 CAP_ETA=0 CAP_WIDTH_BOHR=10 CAP_LAUNCH_Z=-16 CAP_V0=$V0 \
      CAP_WP_SIGMA=3.0 CAP_PROJ_PSEUDO=$WPUPF CAP_N_STEPS=$NS CAP_WRITE_EVERY=5 \
      CAP_OUT_SUBDIR=$SUB ./run > logs_$SUB.log 2>&1
  local st=$?
  echo "[$(date '+%F %T')] DONE benchmark E=$E exit=$st"
  $PY email_run.py "$SUB" "no-CAP classical S(v) benchmark, matched sigma, E=$E eV" "$st" "$FAM" || true
}

for E in 150 300 450 600; do run_E $E; done

$PY email_run.py qvc_bench_cl_s3_E600 "Phase 2 no-CAP S(v) benchmark COMPLETE (150-600 eV) — ALL TASKS DONE" 0 "$FAM" \
   "Study A (CAP) + Phase 2 (no-CAP S(v) benchmark) finished. Study B (no-CAP 300 eV twin) is a prompt awaiting your command." || true
echo "[$(date '+%F %T')] Phase 2 complete. ALL DONE."
