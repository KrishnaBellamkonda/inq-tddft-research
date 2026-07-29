#!/usr/bin/env python3
"""Q1(vi) / todo-3: total norm of the simulation — WP vs classical, easy to compare.

To compare like-for-like, the classical run's reconstructed 1-unit projectile is ADDED to
its electron count, so BOTH runs start at 235 charge units:
  WP run    : N_total electrons (235 = 234 bath + 1 WP).
  classical : N_total electrons (234 bath) + reconstructed projectile (~1, transmits).
The gap that opens = the WP electron being absorbed (the classical projectile is not).
Outputs: qa_vi_total_norm_compare.png + .csv.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inqview.visualisation import style

style.apply_theme()
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"
KEYS = ["mzCAP", "leftfree", "slab", "rightfree", "pzCAP"]

en_wp = pd.read_csv(f"{WP}/observables/electron_number.csv")
en_cl = pd.read_csv(f"{CL}/observables/electron_number.csv")
qi = pd.read_csv(f"{OUT}/qa_i_region_densities.csv")
proj_norm = qi[[f"cl_{k}" for k in KEYS]].sum(axis=1).values      # reconstructed projectile (~1)
t_f = qi.time_au.values

N_cl_e = np.interp(t_f, en_cl.time_au, en_cl.N_total)             # classical electrons (bath)
N_cl_total = N_cl_e + proj_norm                                   # + projectile on top

fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.plot(en_wp.time_au, en_wp.N_total, "C0-", lw=1.9,
        label=f"WP run: total charge (235 → {en_wp.N_total.iloc[-1]:.3f}); loses WP+bath")
ax.plot(t_f, N_cl_total, "C1-", lw=1.9,
        label=f"classical: electrons + projectile (235 → {N_cl_total[-1]:.3f}); loses bath only")
ax.plot(t_f, N_cl_e, "C1--", lw=1.1,
        label=f"classical: electrons only (234 → {N_cl_e[-1]:.3f}) — projectile not counted")
ax.set_xlabel("time (a.u.)"); ax.set_ylabel("total norm / charge units")
ax.set_title("Total norm of the simulation — WP vs classical (both start at 235)", fontsize=9)
ax.legend(fontsize=7, frameon=False, loc="lower left"); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/qa_vi_total_norm_compare.png", dpi=200); plt.close(fig)

pd.DataFrame({"time_au": t_f, "N_total_wp": np.interp(t_f, en_wp.time_au, en_wp.N_total),
              "N_cl_electrons": N_cl_e, "proj_norm": proj_norm,
              "N_cl_total_with_proj": N_cl_total}).to_csv(
    f"{OUT}/qa_vi_total_norm_compare.csv", index=False)
print(f"WP total: 235 -> {en_wp.N_total.iloc[-1]:.3f}  (loses {235-en_wp.N_total.iloc[-1]:.3f})")
print(f"classical total (+proj): 235 -> {N_cl_total[-1]:.3f}  (loses {235-N_cl_total[-1]:.3f}, bath only)")
print(f"projectile norm: start {proj_norm[0]:.3f}  end {proj_norm[-1]:.3f}  (conserved, transmits)")
print("wrote qa_vi_total_norm_compare.png + .csv")
