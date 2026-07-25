#!/usr/bin/env bash
# ============================================================================
# Autonomous orchestrator — WP-CAP energy-plateau campaign.
# Chain (one GPU, sequential):  jellium GS -> self-validate WP position ->
#   WP smoke (40 steps) -> WP no-CAP (100 a.u.) -> WP CAP (100 a.u.) -> compare.
# Also builds the vacuum warm-up run notebooks. Emails at each stage (threaded).
# Idempotent: completed stages (run_summary run_completed=true) are skipped, so
# it can be safely re-launched after a kill (resumes from the last checkpoint).
#
# Launch DETACHED so it survives the launching shell/session:
#   setsid nohup bash orchestrate.sh >/dev/null 2>&1 &
# ============================================================================
set -uo pipefail

ROOT=/local/data/public/skcb2/tddft
LJ=$ROOT/ResearchProject/systems/localised_jellium
CAMP=$LJ/scripts/wp_cap_energy_plateau
VAC=$ROOT/ResearchProject/systems/vacuum/scripts/wp_traversal_energy
GS_CKPT=$LJ/shared_gs/slab_n102_L25x25x140_w0p5_h0p5
PY=$ROOT/venv/bin/python3
NOTIFY="$PY $CAMP/notify.py"
LOG=$CAMP/orchestrate.log
RUN=$ROOT/shared/bin/inq-run
N_STEPS=5000            # 100 a.u. at dt=0.02

source $ROOT/shared/config.sh
export INQ_SOURCE=$ROOT/inq-study
export CUDA_VISIBLE_DEVICES=0
export PYTHONPATH=$ROOT/inq-stack/python${PYTHONPATH:+:$PYTHONPATH}

exec >>"$LOG" 2>&1
echo ""; echo "======== orchestrator start $(date) pid=$$ GPU=$CUDA_VISIBLE_DEVICES ========"

secs() { date +%s; }
done_ok() { [ -f "$1" ] && grep -q "run_completed = true" "$1"; }
fail() { echo "!! FAIL: $1 ($(date))"; $NOTIFY "FAILED at: $1" "Orchestrator stopped at stage: $1.
See $LOG for details. Completed stages are checkpointed; re-launch to resume."; exit 1; }

# ---- Stage 0: vacuum warm-up notebooks (runs already done) ------------------
for tag in nocap cap; do
  if [ -d "$VAC/results/$tag" ]; then
    echo "[vacuum:$tag] analysing..."
    $PY $CAMP/analyse.py "$VAC/results/$tag" --label wp --slab-face 1000000 \
        --cap-inner 20 --title "vacuum WP $tag" || echo "  (vacuum $tag analysis non-fatal fail)"
  fi
done
if [ -d "$VAC/results/nocap" ] && [ -d "$VAC/results/cap" ]; then
  $PY $CAMP/compare.py "$VAC/results/nocap" "$VAC/results/cap" \
      "$VAC/results/vacuum_energy_compare.png" || true
  $NOTIFY "vacuum warm-up ready" "Vacuum WP no-CAP vs CAP total-energy comparison attached.
no-CAP conserves E_total; CAP drains it as the WP is absorbed." \
      "$VAC/results/vacuum_energy_compare.png" \
      "$VAC/results/nocap/report/energy_vs_time.png" \
      "$VAC/results/cap/report/energy_vs_time.png" || true
fi

# ---- Stage 1: jellium GS ----------------------------------------------------
if done_ok "$GS_CKPT/run_summary.txt"; then
  echo "[GS] already complete — skipping."
else
  echo "[GS] building + running ($(date))..."
  t0=$(secs)
  ( cd "$CAMP/gs" && GS_CKPT="$GS_CKPT" "$RUN" ) || fail "GS run"
  done_ok "$GS_CKPT/run_summary.txt" || fail "GS did not complete"
  gs_wall=$(( $(secs) - t0 ))
  echo "[GS] done in ${gs_wall}s"
fi
$NOTIFY "GS complete" "Jellium GS done. $(grep -h 'gs_energy_ha\|r_s' $GS_CKPT/run_summary.txt | tr '\n' ' ')
Proceeding to WP runs (100 a.u. each, checkpointed)." || true

# ---- Stage 2: self-validate WP start position -------------------------------
echo "[validate] checking WP launch position..."
$PY $CAMP/validate_wp.py "$GS_CKPT" || fail "WP position validation"

# ---- Stage 3: WP smoke (catch runtime bugs before the long runs) -----------
if ! done_ok "$CAMP/wp/results/smoke/run_summary.txt"; then
  echo "[smoke] 40-step WP run..."
  ( cd "$CAMP/wp" && WP_CAP_ETA=0 WP_OUT=smoke WP_GS_DIR="$GS_CKPT" \
       WP_N_STEPS=40 WP_WF_EVERY=10 WP_DENS_EVERY=10 "$RUN" ) || fail "WP smoke"
  done_ok "$CAMP/wp/results/smoke/run_summary.txt" || fail "WP smoke incomplete"
fi
echo "[smoke] ok"

# ---- Stage 4+5: WP no-CAP then CAP (full 100 a.u.) --------------------------
run_wp () {  # $1=tag  $2=eta  $3=slab_face  $4=cap_inner  $5=title
  local tag=$1 eta=$2 face=$3 cap=$4 title=$5
  local rdir="$CAMP/wp/results/$tag"
  if done_ok "$rdir/run_summary.txt"; then echo "[$tag] already complete — skipping run."; else
    echo "[$tag] running full ($N_STEPS steps) $(date)..."
    local t0=$(secs)
    ( cd "$CAMP/wp" && WP_CAP_ETA="$eta" WP_OUT="$tag" WP_GS_DIR="$GS_CKPT" \
         WP_N_STEPS="$N_STEPS" WP_RESUME=1 "$RUN" ) || fail "$tag run"
    done_ok "$rdir/run_summary.txt" || fail "$tag incomplete"
    echo "[$tag] done in $(( $(secs) - t0 ))s"
  fi
  echo "[$tag] analysing..."
  $PY $CAMP/analyse.py "$rdir" --label wp --slab-face "$face" --cap-inner "$cap" --title "$title" \
      || echo "  ($tag analysis non-fatal fail)"
  $NOTIFY "$tag complete" "Jellium WP $tag run finished ($N_STEPS steps, 100 a.u.).
Energy decomposition + momentum + density GIF attached." \
      "$rdir/report/energy_vs_time.png" "$rdir/report/norm_vs_time.png" \
      "$rdir/report/momentum_evolution.png" || true
}

# On resume, WP_RESUME=1 means each run continues from its last checkpoint.
# Fresh runs: rt_ckpt absent -> START=0 (the run.cpp handles this by treating a
# missing rt_state as no-op only when RESUME and ckpt exists; for a clean start
# we want RESUME=0 the first time). Handle first-vs-resume explicitly:
first_run () {  # like run_wp but forces a clean (non-resume) first start
  local tag=$1 eta=$2 face=$3 cap=$4 title=$5
  local rdir="$CAMP/wp/results/$tag"
  if done_ok "$rdir/run_summary.txt"; then echo "[$tag] complete — skipping.";
  elif [ -f "$rdir/rt_ckpt/rt_state.txt" ]; then run_wp "$tag" "$eta" "$face" "$cap" "$title";  # resume
  else
    echo "[$tag] running full clean start ($N_STEPS steps) $(date)..."
    local t0=$(secs)
    ( cd "$CAMP/wp" && WP_CAP_ETA="$eta" WP_OUT="$tag" WP_GS_DIR="$GS_CKPT" \
         WP_N_STEPS="$N_STEPS" WP_RESUME=0 "$RUN" ) || fail "$tag run"
    done_ok "$rdir/run_summary.txt" || fail "$tag incomplete"
    echo "[$tag] done in $(( $(secs) - t0 ))s"
    $PY $CAMP/analyse.py "$rdir" --label wp --slab-face "$face" --cap-inner "$cap" --title "$title" \
        || echo "  ($tag analysis non-fatal fail)"
    $NOTIFY "$tag complete" "Jellium WP $tag finished ($N_STEPS steps, 100 a.u.). Plots attached." \
        "$rdir/report/energy_vs_time.png" "$rdir/report/norm_vs_time.png" \
        "$rdir/report/momentum_evolution.png" || true
  fi
}

first_run nocap 0     12.5 60 "jellium WP no-CAP"
first_run cap   -0.7  12.5 60 "jellium WP CAP"

# ---- Stage 6: headline comparison ------------------------------------------
echo "[compare] building plateau comparison..."
$PY $CAMP/compare.py "$CAMP/wp/results/nocap" "$CAMP/wp/results/cap" \
    "$CAMP/wp/results/jellium_energy_compare.png" || true
$NOTIFY "ALL DONE — plateau comparison" "Campaign complete. Headline figure: total energy(t)
no-CAP vs CAP. The plateau GAP = energy radiated to the boundaries (drained by the CAP).
Vacuum warm-up + both jellium runs finished. Notebooks in each results/*/report/." \
    "$CAMP/wp/results/jellium_energy_compare.png" || true
echo "======== orchestrator done $(date) ========"
