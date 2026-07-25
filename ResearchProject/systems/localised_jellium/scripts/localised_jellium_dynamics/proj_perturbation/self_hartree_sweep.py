#!/usr/bin/env python3
"""Run the EMPIRICAL wavepacket self-Hartree diagnostic (self_hartree binary) across
boundary condition (p2 open-z vs p3 full-PBC), grid spacing, box length Lz, and σ.

E_self = 1/2 ∫ n_WP · poisson(n_WP) in the ACTUAL run cell — INQ picks the kernel from
periodicity (p2 → Rozzi 2D cutoff = open-z; p3 → periodic FFT). Compares to the measured
INQ residuals (p2=20.81, p3=21.49 at dx=0.5) so the analytic self-Hartree can be replaced
by this boundary-matched empirical value. GPU0 only.

Writes hypotheses/perturbation_method/self_hartree_empirical.csv.
"""
import os, re, subprocess, csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIN  = HERE/"self_hartree"
OUT  = HERE.parent.parent.parent/"hypotheses/perturbation_method/self_hartree_empirical.csv"
HALF = 12.5

def run(per, Lz, dx, sigma, r=12):
    z = -(HALF + r)
    env = dict(os.environ,
        CUDA_VISIBLE_DEVICES="0",
        LJ_LX="50", LJ_LY="50", LJ_LZ=str(Lz), LJ_N="82",
        LJ_PERIODICITY=str(per), LJ_SPACING=str(dx), LJ_SIGMA=str(sigma), LJ_LAUNCH_Z=str(z))
    out = subprocess.run([str(BIN)], env=env, capture_output=True, text=True).stdout
    m = re.search(r"E_self_ev=([-\d.]+)", out)
    nm = re.search(r"norm=([-\d.]+)", out)
    if not m:
        print("  FAILED per=%s Lz=%s dx=%s sigma=%s\n%s" % (per,Lz,dx,sigma,out[-400:])); return None
    return {"per":per,"Lz":Lz,"dx":dx,"sigma":sigma,"r":r,
            "norm":float(nm.group(1)) if nm else float("nan"),
            "E_self_ev":float(m.group(1))}

rows = []
# boundary contrast at grid-matched dx=0.5 and converged dx=0.3
for per in (2, 3):
    for dx in (0.5, 0.3):
        r = run(per, 120, dx, 0.5); rows.append(r) if r else None
        if r: print("per=%d Lz=120 dx=%.1f  E_self=%.3f eV" % (per, dx, r["E_self_ev"]))
# Lz sweep, open-z (p2), dx=0.5
for Lz in (90, 160, 240):
    r = run(2, Lz, 0.5, 0.5); rows.append(r) if r else None
    if r: print("per=2 Lz=%d dx=0.5  E_self=%.3f eV" % (Lz, r["E_self_ev"]))
# sigma sweep, open-z (p2), dx=0.5
for sg in (0.35, 0.7, 1.0):
    r = run(2, 120, 0.5, sg); rows.append(r) if r else None
    if r: print("per=2 sigma=%.2f dx=0.5  E_self=%.3f eV" % (sg, r["E_self_ev"]))

rows = [r for r in rows if r]
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nwrote", OUT, "(%d rows)" % len(rows))
