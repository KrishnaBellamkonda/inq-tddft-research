#!/usr/bin/env python3
"""Q1(iv): norm overview + bath region densities + slab<->CAP reflection diagnostic.

Adds to qa_i/ii/iii the pieces the user asked for (2026-06-23):
  - total norm N(t) of the system overlaid with per-region norms (the trend);
  - the density NOT yet absorbed but sitting BETWEEN the slab and the CAP
    (left-free [-17.5,-12.5] and right-free [12.5,17.5]) — where reflection hides;
  - the non-zero t=0 baseline in those regions (bath spill-out + the WP launched at
    z=-15.5 in left-free) made explicit.

Bands (Bohr): [-25,-17.5] -zCAP | [-17.5,-12.5] left-free | [-12.5,12.5] slab |
[12.5,17.5] right-free | [17.5,25] +zCAP.  Free regions = left-free + right-free.

WP run: total/bath/wp per band (bath = total - wp). Classical: total per band.
Outputs: qa_iv_norms.png, qa_iv_bath_bands.png, qa_iv_reflection_freeregion.png + .csv.
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
EDGES = np.array([-25.0, -17.5, -12.5, 12.5, 17.5, 25.0])
KEYS = ["mzCAP", "leftfree", "slab", "rightfree", "pzCAP"]
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"


def ftime(p):
    return int(re.search(r"_t(\d+)\.vti", p).group(1)) * DT


def band_sums(path):
    v = load_vti(path)
    dx, dy, dz = v.spacing
    nz = v.data.sum(axis=(0, 1)) * dx * dy           # 1D linear density n(z)
    return np.array([float(nz[(v.z >= lo) & (v.z < hi)].sum() * dz)
                     for lo, hi in zip(EDGES[:-1], EDGES[1:])])


# --- WP run: total & wp bands -> bath bands ---
totf = sorted(glob.glob(f"{WP}/vti/density_total/density_t*.vti"), key=ftime)
t_f = np.array([ftime(p) for p in totf])
tot_b = np.array([band_sums(p) for p in totf])
wp_b = pd.read_csv(f"{OUT}/qa_i_region_densities.csv")[[f"wp_{k}" for k in KEYS]].values
bath_b = tot_b - wp_b

# --- classical run: total (=bath) bands ---
totf_cl = sorted(glob.glob(f"{CL}/vti/density_total/density_t*.vti"), key=ftime)
t_cl = np.array([ftime(p) for p in totf_cl])
tot_cl = np.array([band_sums(p) for p in totf_cl])

# --- norms ---
en_wp = pd.read_csv(f"{WP}/observables/electron_number.csv")
en_cl = pd.read_csv(f"{CL}/observables/electron_number.csv")
N_wp = wp_b.sum(axis=1)            # WP norm in box
N_tot_wp = np.interp(t_f, en_wp.time_au, en_wp.N_total)
N_bath_wp = N_tot_wp - N_wp

# free-region (slab<->CAP) un-absorbed density
i_lf, i_rf = KEYS.index("leftfree"), KEYS.index("rightfree")
wp_free = wp_b[:, i_lf] + wp_b[:, i_rf]
bath_free = bath_b[:, i_lf] + bath_b[:, i_rf]

# ============================ figures ============================
cols = ["C3", "C1", "C2", "C0", "C4"]
labs = ["-z CAP", "left-free", "slab", "right-free", "+z CAP"]
t_wrap = (25.0 - (-15.5)) / 2.711      # ballistic time WP first reaches +z box edge

# (1) norm overview — twin axis (bath/total ~234 left; WP 0..1 right)
fig, axL = plt.subplots(figsize=(6.6, 4.2))
axR = axL.twinx()
l1, = axL.plot(en_wp.time_au, en_wp.N_total, "k-", lw=1.8, label="WP run: total N(t)")
l2, = axL.plot(t_f, N_bath_wp, "C1-", lw=1.3, label="WP run: bath norm")
l3, = axL.plot(en_cl.time_au, en_cl.N_total, "C2--", lw=1.6, label="classical: total N (=bath)")
l4, = axR.plot(t_f, N_wp, "C0-o", ms=3, lw=1.4, label="WP run: WP norm (right axis)")
axL.set_xlabel("time (a.u.)"); axL.set_ylabel("bath / total norm (electrons)")
axR.set_ylabel("WP norm (electrons)", color="C0"); axR.tick_params(axis="y", colors="C0")
axR.set_ylim(0, 1.05)
axL.set_title("Total & component norms vs time", fontsize=9)
axL.legend(handles=[l1, l2, l3, l4], fontsize=7, frameon=False, loc="center right")
axL.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/qa_iv_norms.png", dpi=200); plt.close(fig)

# (2) bath per band — change from t=0 baseline (drainage / redistribution)
fig, ax = plt.subplots(figsize=(6.8, 4.3))
for k in range(5):
    ax.plot(t_f, bath_b[:, k] - bath_b[0, k], "-", color=cols[k], lw=1.5,
            label=f"{labs[k]} (t0={bath_b[0,k]:.1f})")
ax.axhline(0, color="grey", lw=0.6)
ax.set_xlabel("time (a.u.)"); ax.set_ylabel("bath density change  n(t) - n(0)  (electrons)")
ax.set_title("Bath redistribution / drainage per z-band — WP run (Δ from t=0)", fontsize=9)
ax.legend(fontsize=6.8, frameon=False, ncol=2); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/qa_iv_bath_bands.png", dpi=200); plt.close(fig)

# (3) reflection diagnostic — two panels (full range + zoom)
fig, (a0, a1) = plt.subplots(2, 1, figsize=(6.6, 6.2), sharex=True)
for ax in (a0, a1):
    ax.plot(t_f, wp_b[:, i_lf], "C1-", lw=1.7, label="WP LEFT-free [-17.5,-12.5]  (reflection side)")
    ax.plot(t_f, wp_b[:, i_rf], "C0-", lw=1.7, label="WP RIGHT-free [12.5,17.5]  (transmission side)")
    ax.plot(t_f, wp_free, "k-", lw=1.3, label="WP total between slab & CAPs (un-absorbed)")
    ax.axvline(t_wrap, color="purple", ls="--", lw=1.0)
    ax.grid(alpha=0.25)
a0.annotate("WP launched in left-free (z=-15.5) -> starts at 1.0", xy=(0.3, 1.0),
            fontsize=6.5, color="grey", va="top")
a0.annotate("periodic wrap\nbeyond here", xy=(t_wrap, 0.6), fontsize=6.5, color="purple",
            ha="right", rotation=90, va="top")
a0.set_ylabel("density (electrons)"); a0.set_title(
    "Reflection diagnostic — un-absorbed WP between slab & CAPs", fontsize=9)
a0.legend(fontsize=6.8, frameon=False, loc="upper right")
a1.set_ylim(0, 0.12)
a1.set_xlabel("time (a.u.)"); a1.set_ylabel("density (electrons) — ZOOM")
a1.annotate("a refill of LEFT-free after it empties = reflection toward -z CAP",
            xy=(2, 0.11), fontsize=6.5, color="C1", va="top")
fig.tight_layout(); fig.savefig(f"{OUT}/qa_iv_reflection_freeregion.png", dpi=200); plt.close(fig)

# --- csv ---
df = pd.DataFrame({"time_au": t_f, "N_total_wp": N_tot_wp, "N_wp": N_wp, "N_bath_wp": N_bath_wp,
                   "wp_leftfree": wp_b[:, i_lf], "wp_rightfree": wp_b[:, i_rf],
                   "wp_free_total": wp_free, "bath_free_total": bath_free})
for k in range(5):
    df[f"bath_{KEYS[k]}"] = bath_b[:, k]
df.to_csv(f"{OUT}/qa_iv_region_detail.csv", index=False)

print(f"t=0 free-region (slab<->CAP): WP={wp_free[0]:.3f}  bath={bath_free[0]:.3f}")
print(f"peak WP between slab & CAPs: {wp_free.max():.3f} at t={t_f[wp_free.argmax()]:.1f}")
print(f"WP RIGHT-free (transmission side): peak={wp_b[:,i_rf].max():.3f} "
      f"at t={t_f[wp_b[:,i_rf].argmax()]:.1f}  end={wp_b[-1,i_rf]:.3f}")
post = t_f > 3.0   # after the launch transient leaves left-free
print(f"WP LEFT-free: t0={wp_b[0,i_lf]:.3f}  min(after launch)={wp_b[post,i_lf].min():.3f}  "
      f"end={wp_b[-1,i_lf]:.3f}  max(after launch)={wp_b[post,i_lf].max():.3f}")
print(f"  (a rise of LEFT-free after its min = reflection back toward -z CAP)")
print(f"WP total between slab&CAPs: max(after launch)={wp_free[post].max():.3f} "
      f"at t={t_f[post][wp_free[post].argmax()]:.1f}")
print(f"WP reaches +z box edge (wrap) at t~{t_wrap:.1f} a.u.")
print("wrote qa_iv_norms.png, qa_iv_bath_bands.png, qa_iv_reflection_freeregion.png + .csv")
