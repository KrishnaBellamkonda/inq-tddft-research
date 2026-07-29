#!/bin/bash
# Launch the 2026-05-31 sigma-sweep WF reruns. NVML/nvidia-smi is broken; GPUs
# confirmed free via cuInit probe (GPU0/GPU1 ~25GB each). Stagger: sigma3 on GPU0,
# sigma0p5 on GPU1; sigma8 waits for GPU1 to free (smallest N_STEPS=240 for 0p5).
export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
JB=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium

run_one () {
  local dir=$1 gpu=$2
  cd "$JB/$dir" || exit 1
  echo "[$(date +%H:%M:%S)] building+launching $dir on GPU$gpu" > _run.log
  CUDA_VISIBLE_DEVICES=$gpu inq-run >> _run.log 2>&1
  echo "[$(date +%H:%M:%S)] $dir EXIT=$?" >> _run.log
}

# GPU0: sigma3 (N_STEPS~461)
run_one run_wp_n162_L50_E100_sigma3_wf 0 &
PID3=$!
# GPU1: sigma0p5 (N_STEPS=240, fastest) then sigma8 after it frees
( run_one run_wp_n162_L50_E100_sigma0p5_wf 1
  run_one run_wp_n162_L50_E100_sigma8_wf 1 ) &
PID18=$!
echo "launched: sigma3 pid=$PID3 (GPU0); sigma0p5->sigma8 chain pid=$PID18 (GPU1)"
wait
echo "ALL DONE $(date)"
