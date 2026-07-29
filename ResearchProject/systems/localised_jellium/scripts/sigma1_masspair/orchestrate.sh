#!/usr/bin/env bash
# sigma1_masspair serial chain: wp_m2_k4p5 -> wp_m3_k4p5 on ONE GPU.
# Plan: docs/plans/sigma1-masspair-decay-runs.md
#
# Guarantees (user requirement 2026-07-15):
#  * idempotent — completed runs (run_completed=true) are skipped; a partial
#    run (rt_ckpt/rt_state.txt present) relaunches itself with LJ_RESUME=1,
#    losing at most LJ_CKPT_EVERY=200 steps (~15 min).
#  * stall watchdog — if a run's log stops growing for STALL_SECS the hung
#    process is killed and the chain exits 42 (distinct from crash codes) so
#    the supervising session can dispatch an intervention agent.
#  * no build here — requires a pre-built, up-to-date ./run (cap_fix lesson:
#    never cmake per-run). Exit 3 if missing/stale.
set -u
ROOT=/local/data/public/skcb2/tddft
WP=$ROOT/ResearchProject/systems/localised_jellium/scripts/sigma1_masspair/wp
GPU=${SMP_GPU:-1}
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod
STALL_SECS=${SMP_STALL_SECS:-900}
POLL=60

[ -x "$WP/run" ] || { echo "[chain] FATAL: $WP/run missing — build first (inq-run)"; exit 3; }
[ "$WP/run" -nt "$WP/run.cpp" ] || { echo "[chain] FATAL: binary older than run.cpp — rebuild"; exit 3; }

run_one() {  # <name> <inv_mass> <k0>
  local name=$1 inv=$2 k0=$3
  local res=$WP/results/$name log=$WP/${name}.log
  local summary=$res/run_summary.txt state=$res/rt_ckpt/rt_state.txt
  if [ -f "$summary" ] && grep -q "run_completed = true" "$summary"; then
      echo "[chain] $name already complete, skipping"; return 0
  fi
  local resume=0
  if [ -f "$state" ]; then resume=1; echo "[chain] $name resuming ($(grep last_step "$state"))"; fi
  echo "[chain] launching $name (resume=$resume, gpu=$GPU) at $(date '+%F %T')"
  # SIGMA: 1.4142 = wavefunction-width param giving DENSITY std 1.0 Bohr — the
  # packet the plan's spreading/aliasing numbers were computed for (house label
  # sigma_WP=1.41; see plan "sigma convention correction", 2026-07-15).
  ( cd "$WP" && CUDA_VISIBLE_DEVICES=$GPU LJ_OUT=$name LJ_INV_MASS=$inv LJ_K0=$k0 \
      LJ_SIGMA_WP=1.4142135623730951 LJ_LAUNCH_Z=-16.5 LJ_N_STEPS=2500 LJ_DT=0.04 \
      LJ_WRITE_EVERY=8 LJ_WF_EVERY=8 LJ_CKPT_EVERY=200 LJ_RESUME=$resume \
      ./run >> "$log" 2>&1 ) &
  local pid=$!
  local last_size=0 stall=0
  while kill -0 $pid 2>/dev/null; do
      sleep $POLL
      local size; size=$(stat -c%s "$log" 2>/dev/null || echo 0)
      if [ "$size" -gt "$last_size" ]; then last_size=$size; stall=0
      else stall=$((stall+POLL)); fi
      if [ "$stall" -ge "$STALL_SECS" ]; then
          echo "[watchdog] $name: no log growth for ${STALL_SECS}s — killing pid $pid"
          kill $pid; sleep 10; kill -9 $pid 2>/dev/null
          return 42
      fi
  done
  wait $pid; local rc=$?
  echo "[chain] $name exited rc=$rc at $(date '+%F %T')"
  return $rc
}

run_one wp_m2_k4p5 0.5 4.5 || exit $?
run_one wp_m3_k4p5 0.33333333333333331 4.5 || exit $?
echo "[chain] ALL DONE at $(date '+%F %T')"
