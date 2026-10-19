#!/usr/bin/env python3
"""S vs L_z — the DIRECTLY-SIMULATED finite-slab stopping power at each slab
length, per projectile energy. Shows how S approaches the bulk limit as the slab
thickens (the convergence the 1/L_z extrapolation rests on) and where a LONGER
L_z would land.

Data: scripts/slab_sE_sweep/sweep_results.csv (per-cell S at L_z=15,30,45, both
halves). The NEW 2-GPU classical point (E276, L_z=60) is injected here, marked
distinctly — it is not yet merged into the campaign CSV (different run method).

Output: hypotheses/slab_sE_sweep/S_vs_Lz.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import inqview.visualisation.style as st

st.apply_theme()
LJ = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
CSV = LJ / "scripts/slab_sE_sweep/sweep_results.csv"
OUT = Path(__file__).resolve().parent / "S_vs_Lz.png"

df = pd.read_csv(CSV)
energies = sorted(df["E_eV"].unique())            # 31,54,100,167,276
# colour by energy (viridis, low E dark -> high E bright)
cmap = matplotlib.cm.get_cmap("viridis")
cnorm = {e: cmap(i / (len(energies) - 1)) for i, e in enumerate(energies)}

# NEW 2-GPU point (measured this session; slab-stopping-power skill)
NEW = {"E_eV": 276, "L_z": 60, "S_classical": 0.0928}

fig, ax = plt.subplots(figsize=(6.2, 4.6), layout="constrained")

for e in energies:
    sub = df[df["E_eV"] == e].sort_values("L_z")
    c = cnorm[e]
    # classical (filled circles, solid)
    m = sub["S_classical"].notna()
    ax.plot(sub["L_z"][m], sub["S_classical"][m], "-o", color=c, lw=1.6, ms=7, zorder=4)
    # wavepacket (open squares, dashed)
    mw = sub["S_wp"].notna()
    ax.plot(sub["L_z"][mw], sub["S_wp"][mw], "--s", color=c, lw=1.2, ms=6,
            mfc="white", mew=1.3, zorder=3)

# the new 2-GPU classical L_z=60 point for E276, connected to its 45-Bohr neighbour
e276 = df[df["E_eV"] == 276].sort_values("L_z")
x_prev, y_prev = 45, float(e276[e276["L_z"] == 45]["S_classical"].iloc[0])
ax.plot([x_prev, NEW["L_z"]], [y_prev, NEW["S_classical"]], "-", color=cnorm[276],
        lw=1.6, zorder=4)
ax.plot(NEW["L_z"], NEW["S_classical"], "*", color=cnorm[276], ms=20,
        mec="k", mew=1.2, zorder=6)
ax.annotate("L_z=60\n(2-GPU, new)", (NEW["L_z"], NEW["S_classical"]),
            textcoords="offset points", xytext=(8, -4), fontsize=9, va="top")

ax.set_xlabel(r"slab thickness  $L_z$  (Bohr)", fontsize=14)
ax.set_ylabel(r"stopping power  $S$  (eV/Bohr)", fontsize=14)
ax.set_xlim(10, 66)
ax.tick_params(direction="in", which="both", top=True, right=True, length=4)

# legends: energies (colour) + representation (marker/line)
en_handles = [Line2D([], [], color=cnorm[e], marker="o", ls="-",
                     label=f"{int(e)} eV") for e in energies]
rep_handles = [
    Line2D([], [], color="0.3", marker="o", ls="-", label="classical"),
    Line2D([], [], color="0.3", marker="s", ls="--", mfc="white", label="wavepacket"),
    Line2D([], [], color="0.3", marker="*", ls="none", ms=14, mec="k", label="new (2-GPU)"),
]
leg1 = ax.legend(handles=en_handles, loc="upper left", title="energy", frameon=True,
                 edgecolor="0.7", framealpha=1.0, fontsize=10, title_fontsize=10)
ax.add_artist(leg1)
ax.legend(handles=rep_handles, loc="lower right", frameon=True, edgecolor="0.7",
          framealpha=1.0, fontsize=10)

fig.savefig(OUT, dpi=200)
print(f"wrote {OUT}")
# echo the table
for e in energies:
    sub = df[df["E_eV"] == e].sort_values("L_z")
    row = "  ".join(f"L{int(l)}={s:.3f}" for l, s in zip(sub["L_z"], sub["S_classical"]) if pd.notna(s))
    extra = "  L60=0.093(NEW)" if e == 276 else ""
    print(f"E={int(e):>3} eV classical:  {row}{extra}")
