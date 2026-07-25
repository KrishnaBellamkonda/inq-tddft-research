#!/usr/bin/env bash
# Autonomous post-processing watcher for the Phase-3 production pair.
# Waits for BOTH production runs to finish (their run.cpp print "done." at the end),
# then: full analysis (with GIF battery) -> shrink GIFs -> study notebook -> the two
# per-run notebooks. Launch in the BACKGROUND; it polls, so it is fully hands-off.
set -u
ROOT=/local/data/public/skcb2/tddft
LJ=$ROOT/ResearchProject/systems/localised_jellium
P3=$LJ/scripts/qsp_phase3
HYP=$LJ/hypotheses/qsp_phase3
PY=$ROOT/venv/bin/python3
export PYTHONPATH=$ROOT/inq-stack/python
WPLOG=$P3/wp/prod_wp.log
CLLOG=$P3/classical/prod_cl.log
STAMP() { date '+%F %T'; }

echo "[postproc $(STAMP)] waiting for both production runs to print 'done.' ..."
SECONDS=0
while true; do
  wp_done=$(grep -c "done\. wall" "$WPLOG" 2>/dev/null); wp_done=${wp_done:-0}
  cl_done=$(grep -c "done\. " "$CLLOG" 2>/dev/null); cl_done=${cl_done:-0}
  crash=$(cat "$WPLOG" "$CLLOG" 2>/dev/null | grep -ciE "terminate called|segmentation fault|nan_seen = true|what\(\):"); crash=${crash:-0}
  if [ "${wp_done:-0}" -ge 1 ] && [ "${cl_done:-0}" -ge 1 ]; then echo "[postproc $(STAMP)] both runs DONE"; break; fi
  if [ "${crash:-0}" -ge 1 ]; then echo "[postproc $(STAMP)] CRASH detected in a prod log — aborting post-proc"; exit 3; fi
  if [ "$SECONDS" -gt 43200 ]; then echo "[postproc $(STAMP)] TIMEOUT (12 h) — aborting"; exit 4; fi
  sleep 60
done

echo "[postproc $(STAMP)] WP summary:"; grep -E "run_completed|wall_time_s|n_steps|wp_norm_after|launch_z" "$P3/wp/results/p3_wp/run_summary.txt" 2>/dev/null
echo "[postproc $(STAMP)] CL summary:"; grep -E "run_completed|wall_time_s|n_steps" "$P3/classical/results/p3_classical/run_summary.txt" 2>/dev/null

cd "$HYP"
echo "[postproc $(STAMP)] === FULL analysis (with GIF battery) ==="
P3_TAG=p3 $PY analyse_phase3.py && echo "[postproc] analysis OK" || echo "[postproc] analysis ERROR (continuing)"

echo "[postproc $(STAMP)] === shrink study GIFs (470px, 64-colour) ==="
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
print("study gifs:", round(sum(os.path.getsize(f) for f in glob.glob('figs/*.gif'))/1048576, 1), "MB total")
PYG

echo "[postproc $(STAMP)] === study notebook ==="
$PY build_phase3_notebook.py && echo "[postproc] study notebook OK" || echo "[postproc] study notebook ERROR"

echo "[postproc $(STAMP)] === per-run notebooks ==="
RNB=$ROOT/.claude/skills/run-notebook/run_notebook_builder.py
CUDA_VISIBLE_DEVICES="" $PY "$RNB" "$P3/wp/results/p3_wp" "$HYP/p3wp_run_notebook.ipynb" \
  --run-cpp "$P3/wp/run.cpp" --cap-inner 35 --rs 5.666 --launch-z -23.75 --v0 2.711 --lindhard both \
  && echo "[postproc] WP run-notebook OK" || echo "[postproc] WP run-notebook ERROR"
CUDA_VISIBLE_DEVICES="" $PY "$RNB" "$P3/classical/results/p3_classical" "$HYP/p3cl_run_notebook.ipynb" \
  --run-cpp "$P3/classical/run.cpp" --cap-inner 35 --rs 5.666 --launch-z -23.75 --v0 2.711 --lindhard both \
  && echo "[postproc] CL run-notebook OK" || echo "[postproc] CL run-notebook ERROR"

echo "[postproc $(STAMP)] === shrink run-notebook GIFs ==="
$PY - <<'PYG'
import glob, os
from PIL import Image, ImageSequence
TW = 470
for d in ("p3wp_run_notebook_figs", "p3cl_run_notebook_figs"):
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

echo "PHASE3_POSTPROC_DONE $(STAMP)" > "$HYP/POSTPROC_DONE"
echo "[postproc $(STAMP)] ALL DONE — study + 2 run-notebooks built"
