#!/usr/bin/env bash
# ============================================================================
# qsp_phase5 — autonomous WP velocity sweep for the quantum stopping power S(E).
#
# Hands-off, no Claude in the loop. Flow:
#   1. SMOKE GATE: build (INQ_SOURCE=inq-study) + a short v=6.0 run; verify it
#      propagated and wrote observables. On failure -> email abort + exit.
#   2. PRODUCTION, value-first on 2 GPUs (clean+cheap+converged first; the long
#      marginal v=1.3 alone on GPU1). The dispatcher OWNS each ./run in a
#      foreground subshell, so completion = the binary returning (no log polling).
#   3. PER-RUN CHAIN (the moment a run returns 0): analyse_phase5 -> append
#      se_state.csv -> build_se_plot --email -> per-run notebook. State-mutating
#      part is flock-serialised so the two lanes can't corrupt se_state.csv.
#   4. ALL DONE: study notebook + POSTPROC_DONE.
#
# Usage:  nohup bash run_sweep.sh > sweep.log 2>&1 &
# ============================================================================
set -u
ROOT=/local/data/public/skcb2/tddft
LJ=$ROOT/ResearchProject/systems/localised_jellium
P5=$LJ/scripts/qsp_phase5
WP=$P5/wp
HYP=$LJ/hypotheses/qsp_phase5
PY=$ROOT/venv/bin/python3
RNB=$ROOT/.claude/skills/run-notebook/run_notebook_builder.py
EGS=-70.22568216820937
DT=0.04
LOCK=$HYP/.se_state.lock
STAMP() { date '+%F %T'; }

export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod
export PYTHONPATH=$ROOT/inq-stack/python
mkdir -p "$HYP/figs"

email_abort() {  # $1 = reason
  $PY - "$1" <<'PYE' 2>/dev/null || true
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
try:
    from inqview.email import send_run_email
    send_run_email("[lj-wp-se-sweep] ABORTED — smoke gate failed",
                   "qsp_phase5 sweep aborted before production.\n\nReason: " + sys.argv[1],
                   to="chiddukanna@gmail.com")
    print("abort email sent")
except Exception as e:
    print("abort email failed:", e)
PYE
}

# nsteps / write_every for a given tau (echo "NSTEPS WE")
plan() { $PY -c "
tau=min(200.0,200.0/$1); n=int(round(tau/$DT)); we=max(1,n//300)
print(n, we)"; }

# ---------------------------------------------------------------- 1. SMOKE GATE
echo "[sweep $(STAMP)] === SMOKE GATE (build + v=6.0 short run) ==="
SMOKE_N=80; SMOKE_WE=10
( cd "$WP" && CUDA_VISIBLE_DEVICES=0 INQ_SOURCE=$ROOT/inq-study \
    LJ_OUT=p5_smoke LJ_CAP=1 LJ_K0=6.0 LJ_DT=$DT LJ_N_STEPS=$SMOKE_N \
    LJ_WRITE_EVERY=$SMOKE_WE LJ_WF_EVERY=$SMOKE_WE LJ_LAUNCH_Z=-23.75 \
    $ROOT/shared/bin/inq-run > smoke.log 2>&1 )
SR=$WP/results/p5_smoke
ok=1
grep -q "run_completed = true" "$SR/run_summary.txt" 2>/dev/null || ok=0
obs=$($PY -c "import pandas as p;print(len(p.read_csv('$SR/raw/observables/observables.csv')))" 2>/dev/null || echo 0)
[ "${obs:-0}" -ge 4 ] || ok=0
[ -f "$SR/raw/observables/wp_real_space_stats.csv" ] || ok=0
crash=$(grep -ciE "terminate called|segmentation fault|nan_seen = true|what\(\):" "$WP/smoke.log" 2>/dev/null || echo 0)
[ "${crash:-0}" -ge 1 ] && ok=0
if [ "$ok" != "1" ]; then
  echo "[sweep $(STAMP)] SMOKE FAILED (run_completed/obs=$obs/crash=$crash) — see $WP/smoke.log"
  email_abort "smoke gate: run_completed or observables missing (obs rows=$obs, crash=$crash). See $WP/smoke.log"
  exit 5
fi
if [ ! -x "$WP/run" ]; then echo "[sweep $(STAMP)] SMOKE built no executable — abort"; email_abort "no ./run executable after build"; exit 6; fi
echo "[sweep $(STAMP)] SMOKE OK (obs rows=$obs, ./run built). Removing smoke output."
rm -rf "$SR"

# ---------------------------------------------------------------- per-run chain
postproc() {  # $1=tag  $2=v  $3=results_dir
  local tag=$1 v=$2 rdir=$3
  echo "[pp $(STAMP)] $tag (v=$v) — analyse + S(E) + email + notebook"
  ( flock 9
    CUDA_VISIBLE_DEVICES="" $PY "$HYP/analyse_phase5.py" "$rdir" "$tag" "$v" \
      && CUDA_VISIBLE_DEVICES="" $PY "$HYP/build_se_plot.py" --email --note "v=$v ($tag) done"
  ) 9>"$LOCK"
  CUDA_VISIBLE_DEVICES="" $PY "$RNB" "$rdir" "$HYP/${tag}_run_notebook.ipynb" \
    --run-cpp "$WP/run.cpp" --cap-inner 35 --rs 5.666 --launch-z -23.75 --v0 "$v" \
    --e-gs-ha "$EGS" --l-slab 25 --lindhard point \
    && echo "[pp $(STAMP)] $tag notebook OK" || echo "[pp $(STAMP)] $tag notebook ERROR (non-fatal)"
}

# run one velocity to completion on a given GPU, then postproc
run_one() {  # $1=gpu  $2=tag  $3=v
  local gpu=$1 tag=$2 v=$3
  read -r n we < <(plan "$v")
  echo "[sweep $(STAMP)] START $tag v=$v on GPU$gpu (n=$n we=$we)"
  ( cd "$WP" && CUDA_VISIBLE_DEVICES=$gpu \
      env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
      LJ_OUT=$tag LJ_CAP=1 LJ_K0=$v LJ_DT=$DT LJ_N_STEPS=$n \
      LJ_WRITE_EVERY=$we LJ_WF_EVERY=$we LJ_LAUNCH_Z=-23.75 \
      ./run > "$WP/${tag}.log" 2>&1 )
  local rc=$?
  if [ $rc -ne 0 ] || ! grep -q "run_completed = true" "$WP/results/$tag/run_summary.txt" 2>/dev/null; then
    echo "[sweep $(STAMP)] $tag FAILED (rc=$rc) — skipping, continuing sweep"
    return 1
  fi
  echo "[sweep $(STAMP)] DONE $tag (v=$v)"
  postproc "$tag" "$v" "$WP/results/$tag"
}

# ---------------------------------------------------------------- 2. PRODUCTION
echo "[sweep $(STAMP)] === PRODUCTION (2 GPUs, value-first) ==="
# Lane 0 (GPU0): clean/cheap/converged first → high-quality points land early.
( run_one 0 p5_wp_v6p0 6.0
  run_one 0 p5_wp_v5p0 5.0
  run_one 0 p5_wp_v4p0 4.0
  run_one 0 p5_wp_v3p0 3.0 ) & L0=$!
# Lane 1 (GPU1): the long marginal point, alone.
( run_one 1 p5_wp_v1p3 1.3 ) & L1=$!

wait $L0; echo "[sweep $(STAMP)] lane0 done"
wait $L1; echo "[sweep $(STAMP)] lane1 done"

# ---------------------------------------------------------------- 3. STUDY NB
echo "[sweep $(STAMP)] === study notebook ==="
( cd "$HYP" && CUDA_VISIBLE_DEVICES="" $PY build_phase5_notebook.py ) \
  && echo "[sweep $(STAMP)] study notebook OK" || echo "[sweep $(STAMP)] study notebook ERROR"

echo "QSP_PHASE5_SWEEP_DONE $(STAMP)" > "$HYP/POSTPROC_DONE"
echo "[sweep $(STAMP)] ALL DONE — S(E) sweep complete"
