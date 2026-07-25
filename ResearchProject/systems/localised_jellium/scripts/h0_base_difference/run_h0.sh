#!/bin/bash
# ============================================================================
# H0 dispatcher — base WP-vs-classical E_total(0) gap, localised jellium slab.
# Builds against inq-study; runs on GPU $GPU (default 1; GPU0 is busy with the
# qsp_phase5 run). Sequence: GS (L_z=120) -> wp{r4,r40} -> classical{r4,r40}.
# Each projectile builds once (inq-run) then reruns the binary (./run) at the
# second distance. r measured from the near slab face (z=-12.5):
#   r=4  -> z=-16.5 ;  r=40 -> z=-52.5.
# ============================================================================
set -eo pipefail
ROOT=/local/data/public/skcb2/tddft
SHR=$ROOT/inq/install/share
export INQ_SHARE_PATH=$SHR
export PSEUDOPOD_SHARE_PATH=$SHR/pseudopod
export INQ_SOURCE=$ROOT/inq-study
GPU=${GPU:-1}
INQRUN=$ROOT/shared/bin/inq-run
H0=$ROOT/ResearchProject/systems/localised_jellium/scripts/h0_base_difference

echo "[$(date '+%F %T')] H0 start (GPU=$GPU, INQ_SOURCE=inq-study)"

# --- 1. GS (build + run) ----------------------------------------------------
cd "$H0/gs"
echo "[$(date '+%F %T')] build+run H0 GS (slab_n82_L50x50x120)"
CUDA_VISIBLE_DEVICES=$GPU "$INQRUN"
GSE=$(grep ground_state_energy_ha results/run_summary.txt | awk '{print $3}')
echo "[$(date '+%F %T')] GS energy = $GSE Ha   (90-box reference ~ -160.99)"

# --- 2. WP runs (build once @ r4, rerun @ r40) ------------------------------
cd "$H0/wp"
echo "[$(date '+%F %T')] build+run H0 wp r4 (z=-16.5)"
CUDA_VISIBLE_DEVICES=$GPU LJ_OUT=h0_wp_r4  LJ_LAUNCH_Z=-16.5 "$INQRUN"
echo "[$(date '+%F %T')] run H0 wp r40 (z=-52.5)"
CUDA_VISIBLE_DEVICES=$GPU LJ_OUT=h0_wp_r40 LJ_LAUNCH_Z=-52.5 ./run

# --- 3. Classical runs (build once @ r4, rerun @ r40) -----------------------
cd "$H0/classical"
echo "[$(date '+%F %T')] build+run H0 classical r4 (z=-16.5)"
CUDA_VISIBLE_DEVICES=$GPU LJ_OUT=h0_cl_r4  LJ_LAUNCH_Z=-16.5 "$INQRUN"
echo "[$(date '+%F %T')] run H0 classical r40 (z=-52.5)"
CUDA_VISIBLE_DEVICES=$GPU LJ_OUT=h0_cl_r40 LJ_LAUNCH_Z=-52.5 ./run

# --- 4. Collect E_total(0) --------------------------------------------------
echo "[$(date '+%F %T')] === H0 RESULTS: E_total(step 0) in Ha ==="
for tag in wp/results/h0_wp_r4 wp/results/h0_wp_r40 classical/results/h0_cl_r4 classical/results/h0_cl_r40; do
  f="$H0/$tag/raw/observables/observables.csv"
  if [ -f "$f" ]; then
    e0=$(awk -F, 'NR==2{print $3}' "$f")
    echo "  $tag : E_total(0) = $e0 Ha"
  else
    echo "  $tag : MISSING ($f)"
  fi
done
echo "[$(date '+%F %T')] H0 ALL DONE"
