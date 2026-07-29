#!/usr/bin/env bash
# ============================================================================
# rerun_lowspread_dualgpu.sh — LOW-SPREADING vacuum WP reruns, both GPUs.
#
# Autonomous + session-surviving (setsid bash rerun_lowspread_dualgpu.sh &).
# sigma0=3, E=400 eV (k0*sigma0=16.3 -> ~5% transit spread), grid h=0.4 (cutoff
# guard PASS), dt=0.01. Box 30x30x45, one-sided +z CAP z in [7.5,22.5] (W=15,
# eta=-1.0), launch z=-7.5. CAP run 800 steps (full absorption, no wrap remnant);
# no-CAP control 350 steps (stops in CAP region BEFORE it can wrap). Defaults in
# run.cpp; here we override NSTEPS/ETA per run.
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
# eta=-3.5: survival exp(-|eta|W/v)=exp(-3.5*15/5.42)~6e-5, BELOW the log-GIF floor
# (gmax*1e-4) so the wrapped remnant is invisible. (sin^2 CAP averages to 1/2 -> the
# effective absorption is exp(-|eta|W/v); validated: eta=-1.0->0.063, eta=-2.5->0.0010.)
# Reflection stayed 0.000 at eta=-1.0 and -2.5, so |eta|=3.5 is safe (adiabatic W=15).
NOCAP_STEPS=350; CAP_STEPS=800; CAP_ETA=-3.5
LOG=rerun_lowspread.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }
fail(){ say "FAILED: $1"; notify "vacuum low-spread rerun FAILED: $1" "See $(pwd)/$LOG"; exit 1; }

say "=== low-spread vacuum rerun (sigma0=3, E=400eV, h=0.4) — dual GPU ==="
notify "vacuum low-spread rerun started" "sigma0=3 E=400eV h=0.4 box 30x30x45 CAP[7.5,22.5] eta=$CAP_ETA; nocap=$NOCAP_STEPS cap=$CAP_STEPS steps; log $(pwd)/$LOG"

rm -rf results/nocap results/cap results/comparison results/smoke

# --- 1. build once (2-step smoke, GPU1) -------------------------------------
say "building (2-step smoke, GPU1)..."
CUDA_VISIBLE_DEVICES=1 WP_OUT=smoke WP_NSTEPS=2 WP_WF_EVERY=1 inq-run >>"$LOG" 2>&1 || fail "build/smoke"
grep -q "WP injected" "$LOG" || fail "smoke did not inject WP"
rm -rf results/smoke

# --- 2. concurrent: nocap(short) GPU0, cap(full) GPU1 -----------------------
say "launching no-CAP (GPU0, $NOCAP_STEPS steps) + CAP (GPU1, $CAP_STEPS steps)..."
CUDA_VISIBLE_DEVICES=0 WP_OUT=nocap WP_ETA=0        WP_NSTEPS=$NOCAP_STEPS ./run > nocap.run.log 2>&1 &
P0=$!
CUDA_VISIBLE_DEVICES=1 WP_OUT=cap   WP_ETA=$CAP_ETA WP_NSTEPS=$CAP_STEPS   ./run > cap.run.log   2>&1 &
P1=$!
wait $P0; R0=$?; wait $P1; R1=$?
[ $R0 -eq 0 ] || fail "no-CAP run (exit $R0; nocap.run.log)"
[ $R1 -eq 0 ] || fail "CAP run (exit $R1; cap.run.log)"

# --- 3. verification: spread, transverse wrap, absorption, reflection, z-wrap
say "verifying (spread / transverse / absorption / reflection / z-wrap)..."
$PY - <<'PYV' 2>&1 | tee -a "$LOG"
import numpy as np, glob, re
from inqview import load_vti
sig0=3.0; halfz=22.5; halfx=15.0
def load(tag,st):
    return load_vti(sorted(glob.glob(f"results/{tag}/raw/vti/density_wp/density_wp_t{st:06d}.vti"))[0],expect_centered_axis=None) if glob.glob(f"results/{tag}/raw/vti/density_wp/density_wp_t{st:06d}.vti") else None
for tag in ("nocap","cap"):
    fs=sorted(glob.glob(f"results/{tag}/raw/vti/density_wp/density_wp_t*.vti"))
    v0=load_vti(fs[0],expect_centered_axis=None); vT=load_vti(fs[-1],expect_centered_axis=None)
    x,y,z=v0.x,v0.y,v0.z; dV=(x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0])
    def sx(v):
        nx=np.asarray(v.data).sum(axis=(1,2)); w=nx/nx.sum(); c=(w*x).sum(); return float(np.sqrt((w*(x-c)**2).sum()))
    def pz(v):
        return float(z[np.asarray(v.data).sum(axis=(0,1)).argmax()])
    def edge(v):
        nx=np.asarray(v.data).sum(axis=(1,2)); nx=nx/nx.sum(); return float(nx[np.abs(x)>halfx-1].sum())
    # spread at ~transit (t=2.8 -> step 280): frame nearest
    steps=[int(re.search(r'_t(\d+)',f).group(1)) for f in fs]
    it=min(range(len(steps)),key=lambda i:abs(steps[i]*0.01-2.8))
    vt=load_vti(fs[it],expect_centered_axis=None)
    R=sx(vt)/sx(v0)
    normT=np.asarray(vT.data).sum()*dV
    # z-wrap: any frame whose peak_z jumps to negative (wrapped) after going positive?
    wrapped=False; prev=-99
    for f in fs:
        p=pz(load_vti(f,expect_centered_axis=None))
        if prev>10 and p<0: wrapped=True
        prev=p
    print(f"{tag}: N0={np.asarray(v0.data).sum()*dV:.4f} spread@transit={100*(R-1):+.0f}% "
          f"edge_x(tF)={edge(vT):.1e} peak_z(tF)={pz(vT):+.1f} norm(tF)={normT:.4f} z_wrapped={wrapped}")
# reflection check (cap): momentum distribution should have NO -k0 peak
import pandas as pd
mp=sorted(glob.glob("results/cap/raw/observables/momentum_distribution*.csv"))
if mp:
    df=pd.read_csv(mp[0],comment='#'); last=df[df.step==df.step.max()]
    kcol=[c for c in df.columns if 'k' in c.lower()][0]; ncol=[c for c in df.columns if 'n' in c.lower()][-1]
    neg=last[last[kcol]<0][ncol].sum(); pos=last[last[kcol]>0][ncol].sum()
    print(f"cap reflection: neg-k weight / pos-k weight = {neg/max(pos,1e-9):.3f} (small => no reflection)")
PYV

# --- 4. notebooks -----------------------------------------------------------
say "regenerating run + phase notebooks..."
bash regen_notebooks.sh >>"$LOG" 2>&1 || fail "notebook regen"

say "=== low-spread vacuum rerun COMPLETE ==="
notify "vacuum low-spread rerun COMPLETE" "sigma0=3 E=400eV. Verification + geometry in $(pwd)/$LOG. Notebooks:
  results/{nocap,cap}/report/run_report.ipynb
  results/comparison/nocap_vs_cap_comparison.ipynb"
say "done."
