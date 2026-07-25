#!/bin/bash
# ============================================================================
# quantum-stopping-power Phase-2 P2.1 — AUTONOMOUS dispatcher.
# WP (GPU0) + matched classical (GPU1) convergence+CAP test, 40 a.u. (2000 steps),
# CAP on (eta -0.7, 10 Bohr/side) => build against inq-study. Then post-process:
# analysis + notebook, fully unattended. Logs/markers left for the agent to read.
# ============================================================================
set -u
ROOT=/local/data/public/skcb2/tddft
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$INQ_SHARE_PATH/pseudopod
export PATH="$ROOT/shared/bin:$PATH"
export INQ_SOURCE=$ROOT/inq-study          # CAP (perturbations::absorbing) needs inq-study
P2=$ROOT/ResearchProject/systems/localised_jellium/scripts/qsp_phase2
HYP=$ROOT/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase2
VENV=$ROOT/venv/bin/python3
mkdir -p "$HYP"
log(){ echo "[$(date '+%F %T')] $*"; }

log "BUILD+RUN: WP (GPU0) + classical (GPU1), inq-study, CAP on, 2000 steps (40 a.u.)"
( cd "$P2/wp"        && CUDA_VISIBLE_DEVICES=0 LJ_OUT=p2_wp        LJ_CAP=1 inq-run > wp_run.log 2>&1 ) & WP=$!
( cd "$P2/classical" && CUDA_VISIBLE_DEVICES=1 LJ_OUT=p2_classical LJ_CAP=1 inq-run > cl_run.log 2>&1 ) & CL=$!
wait $WP; rc_wp=$?; log "WP finished exit=$rc_wp (log: $P2/wp/wp_run.log)"
wait $CL; rc_cl=$?; log "classical finished exit=$rc_cl (log: $P2/classical/cl_run.log)"

export PYTHONPATH=$ROOT/inq-stack/python
RNB=$ROOT/.claude/skills/run-notebook/run_notebook_builder.py

log "ANALYSIS -> $HYP/analyse_phase2.py"
$VENV "$HYP/analyse_phase2.py" > "$HYP/analyse_phase2.log" 2>&1 && log "analysis OK" || log "analysis ERROR (see analyse_phase2.log)"

log "STUDY NOTEBOOK -> $HYP/build_phase2_notebook.py"
$VENV "$HYP/build_phase2_notebook.py" > "$HYP/build_phase2_notebook.log" 2>&1 && log "study notebook OK" || log "study notebook ERROR (see build_phase2_notebook.log)"

# per-run deep-dive notebooks (run-notebook skill builder); CPU-only post-processing
log "RUN NOTEBOOK (WP) -> p2wp_run_notebook.ipynb"
CUDA_VISIBLE_DEVICES="" $VENV "$RNB" "$P2/wp/results/p2_wp" "$HYP/p2wp_run_notebook.ipynb" \
  --run-cpp "$P2/wp/run.cpp" --cap-inner 25 --rs 5.666 --launch-z -22 --v0 2.711 --lindhard both \
  > "$HYP/p2wp_runnb.log" 2>&1 && log "WP run notebook OK" || log "WP run notebook ERROR (see p2wp_runnb.log)"
log "RUN NOTEBOOK (classical) -> p2cl_run_notebook.ipynb"
CUDA_VISIBLE_DEVICES="" $VENV "$RNB" "$P2/classical/results/p2_classical" "$HYP/p2cl_run_notebook.ipynb" \
  --run-cpp "$P2/classical/run.cpp" --cap-inner 25 --rs 5.666 --launch-z -22 --v0 2.711 --lindhard both \
  --measured-s 0.018632 --measured-v 2.711 \
  > "$HYP/p2cl_runnb.log" 2>&1 && log "classical run notebook OK" || log "classical run notebook ERROR (see p2cl_runnb.log)"

log "P2.1 PIPELINE COMPLETE (wp_exit=$rc_wp cl_exit=$rc_cl)"
echo "wp_exit=$rc_wp cl_exit=$rc_cl done=$(date '+%F %T')" > "$HYP/P2_1_DONE"
