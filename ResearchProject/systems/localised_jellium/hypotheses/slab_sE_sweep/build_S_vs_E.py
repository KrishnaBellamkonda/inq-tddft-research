#!/usr/bin/env python3
"""S(E) — bulk (L_z→∞) electronic stopping power vs projectile energy, classical
and wavepacket traces.

Reads sE_results.csv (one row per energy, written by build_S_vs_invLz.py as each
energy's 5-length sweep completes and is extrapolated to 1/L_z=0). Each point is
the quadratic-extrapolated bulk limit. The shaded band between the traces is the
quantum correction S_wp − S_cl.

Output: hypotheses/slab_sE_sweep/S_vs_E.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import inqview.visualisation.style as st

st.apply_theme()
LJ = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
SE_CSV = LJ / "scripts/slab_sE_sweep/sE_results.csv"
OUT = Path(__file__).resolve().parent / "S_vs_E.png"

df = pd.read_csv(SE_CSV).sort_values("E_eV")
E = df["E_eV"].to_numpy(float)
S_cl = df["S_cl_bulk"].to_numpy(float)
S_wp = df["S_wp_bulk"].to_numpy(float)
print("S(E) bulk-extrapolated points:")
for e, c, w in zip(E, S_cl, S_wp):
    print(f"  E={e:.0f} eV: S_cl={c:.3f}  S_wp={w:.3f}  (quantum {w-c:+.3f}) eV/Bohr")

C_CL, C_WP = st.CLASSICAL, st.WP
fig, ax = plt.subplots(figsize=(5.6, 4.4), layout="constrained")

mcl, mwp = np.isfinite(S_cl), np.isfinite(S_wp)
common = mcl & mwp
if common.any():
    ax.fill_between(E[common], S_cl[common], S_wp[common], color="0.5", alpha=0.18,
                    zorder=1, label="quantum correction")
ax.plot(E[mcl], S_cl[mcl], "-s", color=C_CL, lw=1.8, ms=8, zorder=5, label="classical")
ax.plot(E[mwp], S_wp[mwp], "-o", color=C_WP, lw=1.8, ms=8, zorder=5, label="wavepacket")

ax.set_xlabel(r"projectile energy  $E$  (eV)", fontsize=14)
ax.set_ylabel(r"bulk stopping power  $S$  (eV/Bohr)", fontsize=14)
ax.tick_params(direction="in", which="both", top=True, right=True, length=4)
handles = [
    Line2D([], [], color=C_CL, marker="s", ls="-", ms=8, label="classical"),
    Line2D([], [], color=C_WP, marker="o", ls="-", ms=8, label="wavepacket"),
    Patch(facecolor="0.5", alpha=0.18, label="quantum correction"),
]
ax.legend(handles=handles, loc="best", frameon=True, fancybox=False, edgecolor="0.7",
          framealpha=1.0, fontsize=11)
fig.savefig(OUT, dpi=200)
print(f"wrote {OUT}")
