#!/usr/bin/env bash
# ============================================================================
# qsp_phase5 / rerun_v5_h035 — re-run the v=5 (340 eV) WP on a finer grid
# (h=0.35, E_cut=1096 eV) to remove the h=0.5 high-k aliasing. Autonomous:
#   1. GS at h=0.35 (new checkpoint; the 0.5 GS can't load into a 0.35 grid)
#   2. stability smoke (20 steps) — confirm dt=0.04 is stable on the finer grid
#   3. full WP run (k0=5, tau=40, n=1000)
#   4. chain: analyse_phase5 (NEW E_GS anchor) -> S(E) plot -> email -> run notebook
#
# Usage:  nohup bash run_rerun.sh > rerun.log 2>&1 &
# ============================================================================
set -u
ROOT=/local/data/public/skcb2/tddft
LJ=$ROOT/ResearchProject/systems/localised_jellium
RR=$LJ/scripts/qsp_phase5/rerun_v5_h035
HYP=$LJ/hypotheses/qsp_phase5
PY=$ROOT/venv/bin/python3
RNB=$ROOT/.claude/skills/run-notebook/run_notebook_builder.py
GSDIR=$LJ/shared_gs/slab_n82_L50x50x90_h0p35
GPU=0
STAMP(){ date '+%F %T'; }

export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod
export PYTHONPATH=$ROOT/inq-stack/python

email() {  # $1 subject  $2 body  ($3 attachment optional)
  $PY - "$1" "$2" "${3:-}" <<'PYE' 2>/dev/null || true
import sys
sys.path.insert(0,"/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.email import send_run_email
sub,body,att=sys.argv[1],sys.argv[2],sys.argv[3]
kw={}
import os
tf="/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/email_thread.txt"
if os.path.exists(tf):
    r=open(tf).read().strip()
    if r: kw=dict(in_reply_to=r, references=[r])
try:
    send_run_email(sub, body, attachments=[att] if att else None, to="chiddukanna@gmail.com", **kw)
    print("emailed")
except Exception as e: print("email failed:", e)
PYE
}

# -------------------------------------------------------------- 1. GS at h=0.35
if grep -q "run_completed = true" "$RR/gs/results/run_summary.txt" 2>/dev/null && [ -d "$GSDIR" ]; then
  echo "[rerun $(STAMP)] GS already present — skipping GS build/run"
else
  echo "[rerun $(STAMP)] === GS at h=0.35 (SCF, new checkpoint) ==="
  ( cd "$RR/gs" && CUDA_VISIBLE_DEVICES=$GPU INQ_SOURCE=$ROOT/inq-study \
      LJ_SPACING=0.35 LJ_GS_DIR=$GSDIR $ROOT/shared/bin/inq-run > gs.log 2>&1 )
  if ! grep -q "run_completed = true" "$RR/gs/results/run_summary.txt" 2>/dev/null || [ ! -d "$GSDIR" ]; then
    echo "[rerun $(STAMP)] GS FAILED — see $RR/gs/gs.log"
    email "[lj-wp-se-sweep] v5 re-run ABORTED — GS (h=0.35) failed" "The finer-grid ground state did not converge/save. See $RR/gs/gs.log"
    exit 3
  fi
fi
GSE=$(grep -oE "ground_state_energy_ha = [-0-9.eE+]+" "$RR/gs/results/run_summary.txt" | grep -oE "[-0-9.eE+]+$" | head -1)
echo "[rerun $(STAMP)] GS(h=0.35) energy = $GSE Ha  (vs h=0.5: -70.2257)"
if [ -z "$GSE" ]; then echo "[rerun $(STAMP)] could not parse GS energy — abort"; email "[lj-wp-se-sweep] v5 re-run ABORTED — GS energy unparseable" "No ground_state_energy_ha in GS summary."; exit 4; fi

# -------------------------------------------------------------- 2. STABILITY SMOKE
if [ "${SKIP_SMOKE:-0}" = "1" ]; then
  echo "[rerun $(STAMP)] smoke SKIPPED (dt=0.04 stability already proven at h=0.35)"
else
echo "[rerun $(STAMP)] === stability smoke (20 steps, h=0.35) ==="
( cd "$RR/wp" && CUDA_VISIBLE_DEVICES=$GPU INQ_SOURCE=$ROOT/inq-study \
    LJ_OUT=p5_v5_h035_smoke LJ_CAP=1 LJ_K0=5.0 LJ_DT=0.04 LJ_N_STEPS=20 \
    LJ_WRITE_EVERY=5 LJ_WF_EVERY=5 LJ_LAUNCH_Z=-23.75 LJ_SPACING=0.35 LJ_GS_DIR=$GSDIR \
    $ROOT/shared/bin/inq-run > smoke.log 2>&1 )
SR=$RR/wp/results/p5_v5_h035_smoke
ok=1
grep -q "run_completed = true" "$SR/run_summary.txt" 2>/dev/null || ok=0
nan=$(grep -ciE "nan|terminate called|segmentation fault|what\(\):" "$RR/wp/smoke.log" 2>/dev/null); nan=${nan:-0}
[ "$nan" -ge 1 ] && ok=0
[ -x "$RR/wp/run" ] || ok=0
if [ "$ok" != "1" ]; then
  echo "[rerun $(STAMP)] SMOKE FAILED (nan/instability at h=0.35, dt=0.04) — see $RR/wp/smoke.log"
  email "[lj-wp-se-sweep] v5 re-run ABORTED — smoke unstable at h=0.35" "dt=0.04 may be unstable on the finer grid (nan=$nan). Consider dt=0.02. See $RR/wp/smoke.log"
  exit 5
fi
echo "[rerun $(STAMP)] SMOKE OK — dt=0.04 stable at h=0.35. Removing smoke."
rm -rf "$SR"
fi

# -------------------------------------------------------------- 3. FULL WP RUN
# LEAN config: at h=0.35 each VTI frame is 60 MB and compute is ~60 s/step, so VTIs
# are the cost driver and NONE are needed for the energy-method S. tau=28 (n=700,
# enough for v=5 transit+deposit), observables every 25 steps (28 pts), and the
# wavefunction/wp-density VTIs disabled (LJ_WF_EVERY huge) → ~15 h instead of ~56 h.
echo "[rerun $(STAMP)] === full WP run: v=5.0 (340 eV), h=0.35, tau=28 (n=700, we=25, lean VTI) ==="
( cd "$RR/wp" && CUDA_VISIBLE_DEVICES=$GPU \
    env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
    LJ_OUT=p5_wp_v5p0_h035 LJ_CAP=1 LJ_K0=5.0 LJ_DT=0.04 LJ_N_STEPS=700 \
    LJ_WRITE_EVERY=25 LJ_WF_EVERY=100000 LJ_LAUNCH_Z=-23.75 LJ_SPACING=0.35 LJ_GS_DIR=$GSDIR \
    ./run > prod.log 2>&1 )
RC=$?
RD=$RR/wp/results/p5_wp_v5p0_h035
if [ $RC -ne 0 ] || ! grep -q "run_completed = true" "$RD/run_summary.txt" 2>/dev/null; then
  echo "[rerun $(STAMP)] WP run FAILED (rc=$RC) — see $RR/wp/prod.log"
  email "[lj-wp-se-sweep] v5 re-run FAILED — WP run (h=0.35)" "rc=$RC, no run_completed. See $RR/wp/prod.log"
  exit 6
fi
echo "[rerun $(STAMP)] WP run DONE"

# -------------------------------------------------------------- 4. CHAIN
echo "[rerun $(STAMP)] === analyse (E_GS=$GSE) -> S(E) plot -> email -> notebook ==="
P5_EGS=$GSE CUDA_VISIBLE_DEVICES="" $PY "$HYP/analyse_phase5.py" "$RD" p5_wp_v5p0 5.0
CUDA_VISIBLE_DEVICES="" $PY "$HYP/build_se_plot.py" --email \
  --note "v=5.0 (340 eV) RE-RUN at h=0.35 — aliasing fixed (E_cut 537->1096 eV)"
CUDA_VISIBLE_DEVICES="" $PY "$RNB" "$RD" "$HYP/p5_wp_v5p0_h035_run_notebook.ipynb" \
  --run-cpp "$RR/wp/run.cpp" --cap-inner 35 --rs 5.666 --launch-z -23.75 --v0 5.0 \
  --e-gs-ha "$GSE" --l-slab 25 --lindhard point \
  && echo "[rerun $(STAMP)] run-notebook OK" || echo "[rerun $(STAMP)] run-notebook ERROR (non-fatal)"

echo "V5_H035_RERUN_DONE $(STAMP)" > "$RR/RERUN_DONE"
echo "[rerun $(STAMP)] ALL DONE — corrected 340 eV point in se_state.csv + emailed"
