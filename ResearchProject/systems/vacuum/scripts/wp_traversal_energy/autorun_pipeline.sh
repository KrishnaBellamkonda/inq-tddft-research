#!/usr/bin/env bash
# ============================================================================
# autorun_pipeline.sh — fully autonomous, session-surviving pipeline on GPU 1.
#
# Runs unattended (launch with:  setsid bash autorun_pipeline.sh &  )
#   1. vacuum WP reruns (TRUE vacuum, 30x30x80, one-sided +z CAP) on GPU 1:
#         nocap (WP_ETA=0) then cap (WP_ETA=-0.7), same prebuilt binary.
#   2. regenerate vacuum notebooks (per-run + setup figure + comparison).
#   3. regenerate jellium wp_cap_energy_plateau notebooks (cap, nocap) with the
#      NEW ΔE component + pairwise E_ss/E_ps/E_pp/E_sb/E_pb decomposition plots.
#   4. email the user at start / success / failure (notify.py).
#
# GPU 1 chosen because GPU 0 is occupied by another task (probe: GPU0 52% free,
# GPU1 99% free, 2026-07-27). NVML/nvidia-smi is broken (cosmetic) — compute works.
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
VDIR=$ROOT/ResearchProject/systems/vacuum/scripts/wp_traversal_energy
JCAMP=$ROOT/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau
JWP=$JCAMP/wp/results
PY=$ROOT/venv/bin/python3
export INQ_SOURCE=$ROOT/inq-study
export CUDA_VISIBLE_DEVICES=1
export INQ_SHARE_PATH=${INQ_SHARE_PATH:-$ROOT/inq/install/share}
export PSEUDOPOD_SHARE_PATH=${PSEUDOPOD_SHARE_PATH:-$ROOT/inq/install/share/pseudopod}
export PATH="$ROOT/shared/bin:$PATH"
export PYTHONPATH=$ROOT/inq-stack/python

LOG=$VDIR/autorun_pipeline.log
: > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" ${3:+"$3"} >>"$LOG" 2>&1 || true; }
fail(){ say "FAILED: $1"; notify "pipeline FAILED at: $1" "See $LOG"; exit 1; }

say "=== autonomous pipeline start (GPU $CUDA_VISIBLE_DEVICES) ==="
notify "pipeline started (GPU 1)" "vacuum true-vacuum reruns + all notebook regen; log $LOG"

# --- 1. vacuum runs (binary already built by smoke) -------------------------
cd "$VDIR"
rm -rf results/nocap results/cap results/comparison results/smoke
say "vacuum no-CAP run..."
WP_OUT=nocap WP_ETA=0    ./run >>"$LOG" 2>&1 || fail "vacuum nocap run"
say "vacuum CAP run..."
WP_OUT=cap   WP_ETA=-0.7 ./run >>"$LOG" 2>&1 || fail "vacuum cap run"
grep -H "cell_bohr\|launch\|cap_z\|run_completed" results/nocap/run_summary.txt results/cap/run_summary.txt | tee -a "$LOG"

# --- 2. vacuum notebooks ----------------------------------------------------
say "vacuum notebooks (per-run + setup + comparison)..."
bash regen_notebooks.sh >>"$LOG" 2>&1 || fail "vacuum notebook regen"

# --- 3. jellium notebooks with ΔE decomposition -----------------------------
for tag in cap nocap; do
  say "jellium $tag notebook (+ ΔE component & pairwise plots)..."
  $PY "$JCAMP/analyse.py" "$JWP/$tag" --label wp --dt 0.02 \
      --slab-face 12.5 --cap-inner 60 --cap-lines 60,70 \
      --title "jellium WP-CAP / $tag" >>"$LOG" 2>&1 || fail "jellium $tag notebook"
done

# --- 4. done ----------------------------------------------------------------
say "=== pipeline COMPLETE ==="
{
  echo "Autonomous pipeline complete (GPU 1)."
  echo ""
  echo "Vacuum (true vacuum, 30x30x80, one-sided +z CAP, launch z=-30):"
  grep -h "run_completed\|cell_bohr" "$VDIR"/results/{nocap,cap}/run_summary.txt 2>/dev/null
  echo ""
  echo "Notebooks:"
  echo "  $VDIR/results/comparison/nocap_vs_cap_comparison.ipynb"
  echo "  $VDIR/results/{nocap,cap}/report/run_report.ipynb"
  echo "  $JWP/{cap,nocap}/report/run_report.ipynb  (now with ΔE component + pairwise E_ss/E_ps/E_pp plots)"
} > "$VDIR/PIPELINE_DONE.txt"
# NB: do NOT attach the comparison notebook — its embedded GIFs (~25 MB) exceed
# Gmail's attachment cap and bounce the mail. Text summary only; paths inside.
notify "pipeline COMPLETE" "$(cat "$VDIR/PIPELINE_DONE.txt")"
say "emailed completion; artifacts listed in PIPELINE_DONE.txt"
