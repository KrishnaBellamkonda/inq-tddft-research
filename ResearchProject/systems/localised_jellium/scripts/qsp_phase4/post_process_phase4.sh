#!/usr/bin/env bash
# Autonomous post-processing watcher for the Phase-4 production pair (54 eV / v=2.0).
# Waits for BOTH runs to print "done." then: full analysis (GIF battery) -> shrink
# GIFs -> study notebook -> the two per-run notebooks. Fully hands-off (polls).
# If only the WP ran (classical halted on smoke-fail), set LJ_WP_ONLY=1 to skip the
# classical-done wait.
set -u
ROOT=/local/data/public/skcb2/tddft
LJ=$ROOT/ResearchProject/systems/localised_jellium
P4=$LJ/scripts/qsp_phase4
HYP=$LJ/hypotheses/qsp_phase4
PY=$ROOT/venv/bin/python3
export PYTHONPATH=$ROOT/inq-stack/python
WPLOG=$P4/wp/prod_wp.log
CLLOG=$P4/classical/prod_cl.log
WP_ONLY=${LJ_WP_ONLY:-0}
mkdir -p "$HYP"
STAMP() { date '+%F %T'; }

# --- build-success precheck: both runs must START propagating within 40 min,
#     else a compile failed (silent) -> abort instead of hanging for 18 h. ---
echo "[pp4 $(STAMP)] precheck: waiting for runs to start propagating ..."
SECONDS=0
while true; do
  wp_go=$(grep -ac "starting real-time propagation" "$WPLOG" 2>/dev/null); wp_go=${wp_go:-0}
  cl_go=$(grep -ac "starting real-time propagation" "$CLLOG" 2>/dev/null); cl_go=${cl_go:-0}
  [ "$WP_ONLY" = "1" ] && cl_go=1
  if [ "${wp_go:-0}" -ge 1 ] && [ "${cl_go:-0}" -ge 1 ]; then echo "[pp4 $(STAMP)] both runs propagating — precheck OK"; break; fi
  if [ "$SECONDS" -gt 2400 ]; then echo "[pp4 $(STAMP)] PRECHECK FAIL (40 min, a build likely failed) wp_go=$wp_go cl_go=$cl_go — aborting"; exit 5; fi
  sleep 30
done

echo "[pp4 $(STAMP)] waiting for production runs to print 'done.' (wp_only=$WP_ONLY) ..."
SECONDS=0
while true; do
  wp_done=$(grep -c "done\. wall" "$WPLOG" 2>/dev/null); wp_done=${wp_done:-0}
  cl_done=$(grep -c "done\. parked" "$CLLOG" 2>/dev/null); cl_done=${cl_done:-0}   # specific marker (build logs contain git 'done.')
  [ "$WP_ONLY" = "1" ] && cl_done=1
  crash=$(cat "$WPLOG" "$CLLOG" 2>/dev/null | grep -ciE "terminate called|segmentation fault|nan_seen = true|what\(\):"); crash=${crash:-0}
  if [ "${wp_done:-0}" -ge 1 ] && [ "${cl_done:-0}" -ge 1 ]; then echo "[pp4 $(STAMP)] runs DONE"; break; fi
  if [ "${crash:-0}" -ge 1 ]; then echo "[pp4 $(STAMP)] CRASH in a prod log — aborting"; exit 3; fi
  if [ "$SECONDS" -gt 64800 ]; then echo "[pp4 $(STAMP)] TIMEOUT (18 h) — aborting"; exit 4; fi
  sleep 60
done

echo "[pp4 $(STAMP)] WP summary:"; grep -E "run_completed|wall_time_s|n_steps|wp_norm_after|launch_z" "$P4/wp/results/p4_wp/run_summary.txt" 2>/dev/null
echo "[pp4 $(STAMP)] CL summary:"; grep -E "run_completed|wall_time_s|park_|steps_run|final_z" "$P4/classical/results/p4_classical/run_summary.txt" 2>/dev/null

cd "$HYP"
echo "[pp4 $(STAMP)] === FULL analysis (with GIF battery) ==="
P4_TAG=p4 LJ_WP_ONLY=$WP_ONLY $PY analyse_phase4.py && echo "[pp4] analysis OK" || echo "[pp4] analysis ERROR (continuing)"

echo "[pp4 $(STAMP)] === shrink study GIFs (470px, 64-colour) ==="
$PY - <<'PYG'
import glob, os
from PIL import Image, ImageSequence
TW = 470
for fn in sorted(glob.glob("figs/*.gif")):
    im = Image.open(fn); fr = [f.copy() for f in ImageSequence.Iterator(im)]
    du = [f.info.get("duration", 100) for f in fr]; w, h = fr[0].size
    if w > TW:
        s = TW / w; fr = [f.convert("RGB").resize((int(w*s), int(h*s)), Image.LANCZOS) for f in fr]
    fr = [f.convert("RGB").quantize(colors=64, method=Image.FASTOCTREE) for f in fr]
    fr[0].save(fn, save_all=True, append_images=fr[1:], loop=0, duration=du, optimize=True, disposal=2)
g = glob.glob('figs/*.gif')
if g: print("study gifs:", round(sum(os.path.getsize(f) for f in g)/1048576, 1), "MB total")
PYG

echo "[pp4 $(STAMP)] === study notebook ==="
$PY build_phase4_notebook.py && echo "[pp4] study notebook OK" || echo "[pp4] study notebook ERROR"

echo "[pp4 $(STAMP)] === per-run notebooks ==="
RNB=$ROOT/.claude/skills/run-notebook/run_notebook_builder.py
CUDA_VISIBLE_DEVICES="" $PY "$RNB" "$P4/wp/results/p4_wp" "$HYP/p4wp_run_notebook.ipynb" \
  --run-cpp "$P4/wp/run.cpp" --cap-inner 35 --rs 5.666 --launch-z -23.75 --v0 2.0 --lindhard both \
  && echo "[pp4] WP run-notebook OK" || echo "[pp4] WP run-notebook ERROR"
if [ "$WP_ONLY" != "1" ]; then
CUDA_VISIBLE_DEVICES="" $PY "$RNB" "$P4/classical/results/p4_classical" "$HYP/p4cl_run_notebook.ipynb" \
  --run-cpp "$P4/classical/run.cpp" --cap-inner 35 --rs 5.666 --launch-z -23.75 --v0 2.0 --lindhard both \
  && echo "[pp4] CL run-notebook OK" || echo "[pp4] CL run-notebook ERROR"
fi

echo "[pp4 $(STAMP)] === shrink run-notebook GIFs ==="
$PY - <<'PYG'
import glob, os
from PIL import Image, ImageSequence
TW = 470
for d in ("p4wp_run_notebook_figs", "p4cl_run_notebook_figs"):
    for fn in sorted(glob.glob(os.path.join(d, "*.gif"))):
        im = Image.open(fn); fr = [f.copy() for f in ImageSequence.Iterator(im)]
        du = [f.info.get("duration", 100) for f in fr]; w, h = fr[0].size
        if w > TW:
            s = TW / w; fr = [f.convert("RGB").resize((int(w*s), int(h*s)), Image.LANCZOS) for f in fr]
        fr = [f.convert("RGB").quantize(colors=64, method=Image.FASTOCTREE) for f in fr]
        fr[0].save(fn, save_all=True, append_images=fr[1:], loop=0, duration=du, optimize=True, disposal=2)
    g = glob.glob(os.path.join(d, "*.gif"))
    if g: print(f"{d}:", round(sum(os.path.getsize(f) for f in g)/1048576, 1), "MB")
PYG

echo "PHASE4_POSTPROC_DONE $(STAMP)" > "$HYP/POSTPROC_DONE"
echo "[pp4 $(STAMP)] ALL DONE — study + run-notebook(s) built"
