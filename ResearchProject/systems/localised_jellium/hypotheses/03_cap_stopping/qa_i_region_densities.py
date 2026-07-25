#!/usr/bin/env python3
"""Q1(i): projectile charge in 5 z-bands vs time — WP vs classical.

Bands along z (Bohr): [-25,-17.5] -z CAP | [-17.5,-12.5] left-free |
[-12.5,12.5] slab | [12.5,17.5] right-free | [17.5,25] +z CAP.

Projectile charge per band:
  WP run     -> integrate density_wp VTIs (|psi_WP|^2) over each band.
  classical  -> reconstruct the 1-unit Gaussian (charge std 0.350) at the tracked
                ion z(t) (periodic-wrapped) and integrate its z-marginal per band
                (added ON TOP of the integrated quantum bath, per user 2026-06-23).
Both projectiles carry exactly 1 unit; this shows transit / transmission (+z bands)
/ reflection (-z bands) directly and comparably.
Outputs: qa_i_region_densities.png + .csv in this directory.
"""
import glob
import re
import numpy as np
import pandas as pd
from scipy.special import erf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inqview import load_vti
from inqview.visualisation import style

style.apply_theme()
DT = 0.02
SIG_CHARGE = 0.350           # classical projectile charge std (legacy UPF, as-run)
EDGES = np.array([-25.0, -17.5, -12.5, 12.5, 17.5, 25.0])
LABELS = ["-z CAP", "left-free", "slab", "right-free", "+z CAP"]
KEYS = ["mzCAP", "leftfree", "slab", "rightfree", "pzCAP"]
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"


def ftime(p):
    return int(re.search(r"_t(\d+)\.vti", p).group(1)) * DT


def zmarginal(vf):
    """1D linear density n(z) = integral over x,y, and the z axis."""
    dx, dy, _ = vf.spacing
    return vf.z, vf.data.sum(axis=(0, 1)) * dx * dy


def band_sum(z, nz, edges):
    dz = z[1] - z[0]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (z >= lo) & (z < hi)
        out.append(float(nz[m].sum() * dz))
    return np.array(out)


def gauss_band_frac(edges, z0, s=SIG_CHARGE):
    """Fraction of a unit 1D Gaussian(mean z0, std s) inside each band."""
    c = erf((edges - z0) / (s * np.sqrt(2.0)))
    return 0.5 * (c[1:] - c[:-1])


# --- WP projectile: density_wp per band ---
wpf = sorted(glob.glob(f"{WP}/vti/density_wp/density_t*.vti"), key=ftime)
t_f = np.array([ftime(p) for p in wpf])
wp_bands = np.zeros((len(wpf), 5))
for i, p in enumerate(wpf):
    z, nz = zmarginal(load_vti(p))
    wp_bands[i] = band_sum(z, nz, EDGES)

# --- classical projectile: reconstructed Gaussian per band, at tracked z(t) ---
trk = pd.read_csv(f"{CL}/observables/electron_track.csv").drop_duplicates("step")
z_ion_f = np.interp(t_f, trk.time_au, trk.z)
z_ion_wrapped = ((z_ion_f + 25.0) % 50.0) - 25.0        # periodic fold into [-25,25]
cl_bands = np.array([gauss_band_frac(EDGES, z0) for z0 in z_ion_wrapped])

# --- plot ---
fig, ax = plt.subplots(figsize=(7.0, 4.4))
cols = ["C3", "C1", "C2", "C0", "C4"]
for k in range(5):
    ax.plot(t_f, wp_bands[:, k], "-", color=cols[k], lw=1.6, label=f"WP: {LABELS[k]}")
    ax.plot(t_f, cl_bands[:, k], "--", color=cols[k], lw=1.3, alpha=0.8)
ax.plot([], [], "k-", label="WP (solid)")
ax.plot([], [], "k--", label="classical (dashed)")
ax.set_xlabel("time (a.u.)")
ax.set_ylabel("projectile charge in band (electrons)")
ax.set_title("Projectile charge per z-band vs time — WP vs classical (slab, 100 eV)",
             fontsize=9)
ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{OUT}/qa_i_region_densities.png", dpi=200)
plt.close(fig)

df = pd.DataFrame({"time_au": t_f, "z_ion_wrapped": z_ion_wrapped})
for k in range(5):
    df[f"wp_{KEYS[k]}"] = wp_bands[:, k]
    df[f"cl_{KEYS[k]}"] = cl_bands[:, k]
df.to_csv(f"{OUT}/qa_i_region_densities.csv", index=False)

print("end-state projectile charge per band [-zCAP, left, slab, right, +zCAP]:")
print("  WP       :", np.round(wp_bands[-1], 3), " sum=%.3f" % wp_bands[-1].sum())
print("  classical:", np.round(cl_bands[-1], 3), " sum=%.3f (z_ion=%.1f)"
      % (cl_bands[-1].sum(), z_ion_wrapped[-1]))
print("  WP totals in-box (sum bands):", np.round(wp_bands.sum(axis=1)[[0, -1]], 3))
print("wrote qa_i_region_densities.png + .csv")
