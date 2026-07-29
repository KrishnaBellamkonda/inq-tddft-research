#!/usr/bin/env bash
# ============================================================================
# effmass_sigma1 LEAN orchestrator — fully detached (survives SSH logout).
# Plan: docs/plans/effmass-sigma1-lean-rerun.md
#
# Chain (idempotent):
#   1. Ensure the lean GS exists (wait for a running GS; launch it if absent).
#   2. Launch WP (GPU 0) + classical twin (GPU 1) CONCURRENTLY.
#   3. Wait for both; write ORCHESTRATOR_STATUS.txt at every transition.
#
# Launch:  setsid nohup ./orchestrate.sh > orchestrate.log 2>&1 &
# ============================================================================
set -u
BASE="/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/effmass_sigma1"
GSCK="/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/shared_gs/slab_n52_L40x40x80_dx0p333"
STATUS="$BASE/ORCHESTRATOR_STATUS.txt"

# Environment (bashrc not sourced in a detached shell)
export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod

log(){ echo "[$(date +%F' '%T)] $*" | tee -a "$STATUS"; }

log "orchestrator started (pid $$, pgid $(ps -o pgid= -p $$ | tr -d ' '))"

# ---- Stage 1: ensure GS -----------------------------------------------------
gs_done(){ [[ -f "$GSCK/spin_density" ]] && grep -q "run_completed = true" "$BASE/gs/results/run_summary.txt" 2>/dev/null; }

if ! gs_done; then
  log "GS not complete yet — waiting up to 3 h (relaunching if it died)"
  for i in $(seq 1 180); do
    gs_done && break
    # if no GS binary is running AND checkpoint absent -> (re)launch GS ourselves
    if ! pgrep -f "$BASE/gs/run" >/dev/null 2>&1 && ! gs_done; then
      # small grace: the harness-launched GS may be between build and run
      sleep 90
      if ! pgrep -f "$BASE/gs/run" >/dev/null 2>&1 && ! gs_done; then
        log "no GS process found — (re)launching GS on GPU 0"
        ( cd "$BASE/gs" && CUDA_VISIBLE_DEVICES=0 inq-run >> gs_build_run.log 2>&1 )
        log "GS relaunch finished with exit $?"
      fi
    fi
    sleep 60
  done
fi

if ! gs_done; then
  log "FATAL: GS still not complete after wait — aborting orchestration"
  exit 2
fi
log "GS COMPLETE: $(grep -E 'ground_state_energy_ha|r_s' "$BASE/gs/results/run_summary.txt" | tr '\n' ' ')"

# ---- Stage 2: launch WP (GPU 0) + classical (GPU 1) concurrently ------------
launch_wp(){
  cd "$BASE/wp" || return 2
  export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
  export CUDA_VISIBLE_DEVICES=0
  export EM_OUT=sigma1 EM_DT=0.04 EM_N_STEPS=900 EM_CKPT_EVERY=300 \
         EM_WRITE_EVERY=20 EM_CAP=1 EM_RESUME=0 EM_FOCUS_DIST=4.0
  log "WP launching on GPU 0 (900 steps, ckpt every 300)"
  inq-run > rt_run_lean.log 2>&1
  local rc=$?
  log "WP finished, exit=$rc  $(grep -E 'run_completed|wall_time_s' results/sigma1/run_summary*.txt 2>/dev/null | tr '\n' ' ')"
  return $rc
}
launch_cl(){
  cd "$BASE/classical" || return 2
  export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
  export CUDA_VISIBLE_DEVICES=1
  export EM_OUT=classical EM_DT=0.04 EM_N_STEPS=900 EM_CKPT_EVERY=300 \
         EM_WRITE_EVERY=20 EM_CAP=1 EM_RESUME=0
  log "classical twin launching on GPU 1 (900 steps, ckpt every 300)"
  inq-run > rt_run_lean.log 2>&1
  local rc=$?
  log "classical finished, exit=$rc  $(grep -E 'run_completed|wall_time_s|park_' results/classical/run_summary.txt 2>/dev/null | tr '\n' ' ')"
  return $rc
}

launch_wp & WP_PID=$!
launch_cl & CL_PID=$!
log "runs dispatched (wp shell $WP_PID, classical shell $CL_PID) — waiting"

wait $WP_PID; WP_RC=$?
wait $CL_PID; CL_RC=$?

log "ALL DONE: wp_exit=$WP_RC classical_exit=$CL_RC"
log "results: $BASE/wp/results/sigma1  +  $BASE/classical/results/classical"
