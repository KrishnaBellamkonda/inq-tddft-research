#!/usr/bin/env bash
# ============================================================================
# rerun_compact_dualgpu.sh — COMPACT non-dispersing vacuum WP reruns, both GPUs.
#
# Autonomous + session-surviving (launch: setsid bash rerun_compact_dualgpu.sh &).
# Params (user-chosen 2026-07-27): sigma0=3, E=100 eV (k0*sigma0=8.1 -> ~4.6%
# transit expansion), box 30x30x40, one-sided +z CAP z in [10,20], launch z=-5,
# 600 steps. Defaults now baked into run.cpp.
#   1. build once (GPU 1 smoke) with the new defaults.
#   2. run no-CAP on GPU 0 and CAP on GPU 1 CONCURRENTLY (same binary).
#   3. run + phase notebooks (per-run reports, setup figure, comparison).
#   4. quick dispersion/geometry verification; email summary.
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
JCAMP=$ROOT/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau
PY=$ROOT/venv/bin/python3
export INQ_SOURCE=$ROOT/inq-study
export INQ_SHARE_PATH=${INQ_SHARE_PATH:-$ROOT/inq/install/share}
export PSEUDOPOD_SHARE_PATH=${PSEUDOPOD_SHARE_PATH:-$ROOT/inq/install/share/pseudopod}
export PATH="$ROOT/shared/bin:$PATH"
export PYTHONPATH=$ROOT/inq-stack/python
LOG=rerun_compact.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }
fail(){ say "FAILED: $1"; notify "vacuum compact rerun FAILED: $1" "See $(pwd)/$LOG"; exit 1; }

say "=== compact vacuum rerun (sigma0=3, E=100eV, 30x30x40) — dual GPU ==="
notify "vacuum compact rerun started (GPU0=nocap, GPU1=cap)" "sigma0=3 E=100eV box 30x30x40 CAP[10,20] launch -5; log $(pwd)/$LOG"

rm -rf results/nocap results/cap results/comparison results/smoke

# --- 1. build once (compile new defaults) via a 2-step smoke on GPU 1 --------
say "building binary (2-step smoke on GPU 1)..."
CUDA_VISIBLE_DEVICES=1 WP_OUT=smoke WP_NSTEPS=2 WP_WF_EVERY=1 inq-run >>"$LOG" 2>&1 || fail "build/smoke"
grep -q "WP injected" "$LOG" || fail "smoke did not inject WP"
rm -rf results/smoke

# --- 2. concurrent runs: nocap on GPU0, cap on GPU1 --------------------------
say "launching no-CAP (GPU0) + CAP (GPU1) concurrently..."
CUDA_VISIBLE_DEVICES=0 WP_OUT=nocap WP_ETA=0    ./run > nocap.run.log 2>&1 &
P0=$!
CUDA_VISIBLE_DEVICES=1 WP_OUT=cap   WP_ETA=-0.7 ./run > cap.run.log   2>&1 &
P1=$!
say "  nocap pid=$P0 (GPU0), cap pid=$P1 (GPU1); waiting..."
wait $P0; R0=$?; wait $P1; R1=$?
[ $R0 -eq 0 ] || fail "no-CAP run (exit $R0; see nocap.run.log)"
[ $R1 -eq 0 ] || fail "CAP run (exit $R1; see cap.run.log)"
grep -H "run_completed\|cell_bohr\|launch\|cap_z" results/nocap/run_summary.txt results/cap/run_summary.txt | tee -a "$LOG"

# --- 3. run + phase notebooks -----------------------------------------------
say "regenerating run + phase notebooks..."
bash regen_notebooks.sh >>"$LOG" 2>&1 || fail "notebook regen"

# --- 4. dispersion / geometry verification ----------------------------------
say "verifying dispersion + injection..."
$PY - <<'PYV' 2>&1 | tee -a "$LOG"
import numpy as np, glob, re
from inqview import load_vti
sig0=3.0
for tag in ("nocap","cap"):
    fs=sorted(glob.glob(f"results/{tag}/raw/vti/density_wp/density_wp_t*.vti"))
    v0=load_vti(fs[0],expect_centered_axis=None); dV=(v0.x[1]-v0.x[0])*(v0.y[1]-v0.y[0])*(v0.z[1]-v0.z[0])
    def sig_z(v):
        z=v.z; nz=np.asarray(v.data).sum(axis=(0,1)); w=nz/nz.sum(); c=(w*z).sum(); return float(np.sqrt((w*(z-c)**2).sum()))
    def sig_x(v):
        x=v.x; nx=np.asarray(v.data).sum(axis=(1,2)); w=nx/nx.sum(); c=(w*x).sum(); return float(np.sqrt((w*(x-c)**2).sum()))
    n0=np.asarray(v0.data).sum()*dV
    # transit frame ~ when peak reaches z~10 (CAP inner) for nocap; use mid-run
    sx0=sig_x(v0)*np.sqrt(2); sz0=sig_z(v0)*np.sqrt(2)   # -> wavefunction sigma
    vT=load_vti(fs[len(fs)//2],expect_centered_axis=None); sxT=sig_x(vT)*np.sqrt(2)
    nT=np.asarray(load_vti(fs[-1],expect_centered_axis=None).data).sum()*dV
    print(f"{tag}: N(t0)={n0:.4f} sigma_wf(t0)~{sx0:.2f} (target {sig0}) sigma_x(mid)~{sxT:.2f} expand={100*(sxT/sx0-1):+.0f}% norm(tF)={nT:.3f}")
PYV

# --- 5. done ----------------------------------------------------------------
say "=== compact vacuum rerun COMPLETE ==="
notify "vacuum compact rerun COMPLETE" "sigma0=3 E=100eV 30x30x40. Notebooks:
  results/{nocap,cap}/report/run_report.ipynb
  results/comparison/nocap_vs_cap_comparison.ipynb
Verification + geometry in $(pwd)/$LOG"
say "done."
