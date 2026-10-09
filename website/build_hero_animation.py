#!/usr/bin/env python3
"""Hero animation for the project website.

Animates the corrected total-energy traces E_tot(t) - E_GS for the E=100 eV,
L_z=45 Bohr jellium-slab pair (classical vs wavepacket), classical climbing
visibly above the wavepacket -- the site's headline, visualised.

Physics is IDENTICAL to the vetted static energy_vs_time.png: it reuses the
slab_stopping skill loaders (extensive_kinetic_total + norm-correction, +E_PB
for the classical half) and the same per-box E_GS. No re-derivation.

Run:  /local/data/public/skcb2/tddft/venv/bin/python3 website/build_hero_animation.py
Out:  website/assets/hero/hero_energy_E100_Lz45.gif  (+ .mp4 if ffmpeg present)
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, "/local/data/public/skcb2/tddft/.claude/skills/slab-stopping-power")
import slab_stopping as SS  # noqa: E402
from inqview.visualisation import style  # noqa: E402

HA = SS.HA_TO_EV
LJ = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
RUNS = LJ / "scripts/slab_sE_sweep/runs"
OUT = Path(__file__).resolve().parent / "assets/hero"
OUT.mkdir(parents=True, exist_ok=True)

LZ = 45
E_GS_LZ45 = 386.538854026  # Ha, per box (from build_results_gallery.py)
C_CL, C_WP = style.CLASSICAL, style.WP


def corrected_trace(half):
    """Replicate build_results_gallery.load_energy_trace for E100/Lz45."""
    sub = "classical" if half == "classical" else "wp"
    od = RUNS / f"E100_Lz{LZ}" / "results" / sub / "raw/observables"
    ext = SS.extensive_kinetic_total(od)          # step,time_au,energy_total,e_total_corr_ha,...
    if half == "classical":
        ix = SS.concat_obs(od, "interactions")
        m = pd.merge(ext, ix[["step", "e_pb"]], on="step")
        e_corr = m["e_total_corr_ha"] + m["e_pb"]
        t = m["time_au"].to_numpy(float)
    else:
        e_corr = ext["e_total_corr_ha"]
        t = ext["time_au"].to_numpy(float)
    y = (e_corr.to_numpy(float) - E_GS_LZ45) * HA
    return t, y


tc, yc = corrected_trace("classical")
tw, yw = corrected_trace("wavepacket")

S_cl, S_wp = 0.384, 0.214  # bulk-extrapolated, from sE_results.csv (annotation only)

tmax = max(tc.max(), tw.max())
ymax = max(yc.max(), yw.max())
ymin = min(0.0, yc.min(), yw.min())
pad = 0.08 * (ymax - ymin)

plt.rcParams.update({"font.size": 13, "axes.linewidth": 1.1})
fig, ax = plt.subplots(figsize=(7.2, 4.5), dpi=120)
ax.set_xlim(0, tmax)
ax.set_ylim(ymin - pad, ymax + pad)
ax.set_xlabel("time  (a.u.)")
ax.set_ylabel(r"$E_{\mathrm{tot}}(t) - E_{\mathrm{GS}}$  (eV)")
ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
ax.tick_params(direction="in", which="both", top=True, right=True, length=4)
ax.set_title("Energy deposited in the jellium slab  ·  E = 100 eV", fontsize=13)

(line_cl,) = ax.plot([], [], color=C_CL, lw=2.4, label="classical point charge")
(line_wp,) = ax.plot([], [], color=C_WP, lw=2.4, label="quantum wavepacket")
dot_cl = ax.plot([], [], "o", color=C_CL, ms=7)[0]
dot_wp = ax.plot([], [], "o", color=C_WP, ms=7)[0]
ax.legend(loc="lower right", frameon=True, edgecolor="0.7", framealpha=0.95)
fig.tight_layout()

N = 130


def frac_slice(t, y, f):
    k = max(2, int(np.ceil(f * len(t))))
    return t[:k], y[:k]


def update(i):
    f = (i + 1) / N
    xc, vc = frac_slice(tc, yc, f)
    xw, vw = frac_slice(tw, yw, f)
    line_cl.set_data(xc, vc)
    line_wp.set_data(xw, vw)
    dot_cl.set_data([xc[-1]], [vc[-1]])
    dot_wp.set_data([xw[-1]], [vw[-1]])
    return line_cl, line_wp, dot_cl, dot_wp


anim = FuncAnimation(fig, update, frames=N, interval=45, blit=True)

gif = OUT / "hero_energy_E100_Lz45.gif"
anim.save(gif, writer=PillowWriter(fps=24))
print(f"wrote {gif}")

try:
    mp4 = OUT / "hero_energy_E100_Lz45.mp4"
    anim.save(mp4, writer="ffmpeg", fps=24, dpi=120,
              savefig_kwargs={"facecolor": "white"})
    print(f"wrote {mp4}")
except Exception as e:
    print(f"mp4 skipped ({type(e).__name__}: {e})")

# static poster frame (final state) for the <video>/fallback
fig2, ax2 = plt.subplots(figsize=(7.2, 4.5), dpi=120)
ax2.plot(tc, yc, color=C_CL, lw=2.4, label="classical point charge")
ax2.plot(tw, yw, color=C_WP, lw=2.4, label="quantum wavepacket")
ax2.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
ax2.set_xlabel("time  (a.u.)")
ax2.set_ylabel(r"$E_{\mathrm{tot}}(t) - E_{\mathrm{GS}}$  (eV)")
ax2.set_title("Energy deposited in the jellium slab  ·  E = 100 eV", fontsize=13)
ax2.tick_params(direction="in", which="both", top=True, right=True, length=4)
ax2.legend(loc="lower right", frameon=True, edgecolor="0.7", framealpha=0.95)
fig2.tight_layout()
poster = OUT / "hero_energy_E100_Lz45_poster.png"
fig2.savefig(poster, dpi=120)
print(f"wrote {poster}")
