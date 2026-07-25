#!/usr/bin/env python3
"""Q1(iii): cumulative absorbed norm (CAP overflow) vs time — WP vs classical.

Localised-jellium cap_stopping runs (p5_wp, p5_classical), 100 eV, two-sided sin2 CAP.

Framing (user-agreed 2026-06-23): the jellium BATH overflows into the CAP in BOTH
runs, computed identically by the integrated-charge approach (electron count). The WP
run absorbs an EXTRA component — the WP electron itself (it IS electron density). The
classical projectile is an external 1-unit Gaussian layer added ON TOP of the
integrated bath charge; it TRANSMITS (Ehrenfest ion, not subject to the CAP) and so
contributes ZERO absorbed norm.

Sources (per step / per frame, no invented values):
  N_total(t): raw/observables/electron_number.csv  (per step, authoritative count).
  N_wp(t):    integral of density_wp VTIs (91 frames) = surviving WP norm.
  N_bath = N_total - N_wp  (run-independent bath, CONTEXT density decomposition).
Validation: abs_wpself + abs_bath == abs_combined (sum check).
Outputs: qa_iii_absorbed_norm.png + .csv in this directory.
"""
import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inqview import load_vti
from inqview.visualisation import style

style.apply_theme()
DT = 0.02
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"

# --- N_total(t), per step (authoritative electron count) ---
en_wp = pd.read_csv(f"{WP}/observables/electron_number.csv")
en_cl = pd.read_csv(f"{CL}/observables/electron_number.csv")

# --- N_wp(t) by integrating the density_wp VTIs ---
def ftime(p):
    return int(re.search(r"_t(\d+)\.vti", p).group(1)) * DT

def integrate(path):
    v = load_vti(path)
    return float(v.data.sum() * v.spacing[0] * v.spacing[1] * v.spacing[2])

wpf = sorted(glob.glob(f"{WP}/vti/density_wp/density_t*.vti"), key=ftime)
t_f = np.array([ftime(p) for p in wpf])
N_wp = np.array([integrate(p) for p in wpf])

# N_total at frame times -> bath at frame times
Ntot_f = np.interp(t_f, en_wp.time_au, en_wp.N_total)
N_bath = Ntot_f - N_wp

# --- cumulative absorbed (overflow) components ---
abs_comb_step = en_wp.N_total.iloc[0] - en_wp.N_total          # WP run, per step
abs_comb_f = Ntot_f[0] - Ntot_f                                # WP run, at frames
abs_wpself = N_wp[0] - N_wp                                    # WP electron absorbed
abs_bath_wp = N_bath[0] - N_bath                               # bath overflow (WP run)
abs_bath_cl = en_cl.N_total.iloc[0] - en_cl.N_total           # bath overflow (classical)

# sum check (validation)
resid = np.max(np.abs((abs_wpself + abs_bath_wp) - abs_comb_f))

# --- plot ---
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(en_wp.time_au, abs_comb_step, "k-", lw=1.9,
        label=f"WP run — combined: {abs_comb_f[-1]:.3f}")
ax.plot(t_f, abs_wpself, "C0-o", ms=3, lw=1.4,
        label=f"WP run — WP electron (self): {abs_wpself[-1]:.3f}")
ax.plot(t_f, abs_bath_wp, "C1-", lw=1.5,
        label=f"WP run — bath overflow: {abs_bath_wp[-1]:.3f}")
ax.plot(en_cl.time_au, abs_bath_cl, "C1--", lw=1.7,
        label=f"classical — bath overflow: {abs_bath_cl.iloc[-1]:.3f}")
ax.plot(en_cl.time_au, 0.0 * en_cl.time_au, "C2:", lw=1.6,
        label="classical — projectile (transmits): 0")
ax.set_xlabel("time (a.u.)")
ax.set_ylabel("cumulative absorbed norm  (electrons)")
ax.set_title("CAP overflow vs time — WP vs classical (localised slab, 100 eV)",
             fontsize=9)
ax.legend(fontsize=7, frameon=False, loc="upper left")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{OUT}/qa_iii_absorbed_norm.png", dpi=200)
plt.close(fig)

pd.DataFrame({
    "time_au": t_f, "N_wp": N_wp, "N_total_wprun": Ntot_f, "N_bath_wprun": N_bath,
    "abs_wpself": abs_wpself, "abs_bath_wprun": abs_bath_wp,
}).to_csv(f"{OUT}/qa_iii_absorbed_norm.csv", index=False)

print(f"WP run : combined={abs_comb_f[-1]:.3f}  WP-self={abs_wpself[-1]:.3f}  "
      f"bath={abs_bath_wp[-1]:.3f}   (sum-check resid={resid:.2e})")
print(f"classical: bath overflow={abs_bath_cl.iloc[-1]:.3f}  projectile absorbed=0 (transmits)")
print("wrote qa_iii_absorbed_norm.png + .csv")
