#!/usr/bin/env python3
"""Phase-1 GS validation figure: the converged slab density (xz slice, log).
Confirms the background localises the electrons in |z|<12.5 at the new
50×50×70 geometry, with the new region layout marked. Loads via the canonical
inqview.load_vti (physical order — NEVER fftshift a VTI)."""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview import load_vti
from inqview.visualisation import style

style.apply_theme()
HERE = os.path.dirname(os.path.abspath(__file__))
GSDIR = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
         "scripts/qsp_phase1/gs/results/density_gs_system")

vti = sorted(glob.glob(f"{GSDIR}/*.vti"))[0]
d = load_vti(vti, expect_centered_axis="z")          # asserts physical order; slab at centre
n, x, y, z = d.data, d.x, d.y, d.z
iy = len(y) // 2
sl = n[:, iy, :]                                     # (x, z)

# interior-flat check: mean density in |z|<10, |x|<20 vs its std
Zc, Xc = np.meshgrid(z, x)
interior = (np.abs(Zc) < 10) & (np.abs(Xc) < 20)
n_int = sl[interior]
print(f"interior density mean={n_int.mean():.5f}  std/mean={n_int.std()/n_int.mean():.3%}  "
      f"(flat interior = small)")

cmap = style.cmap_for("sequential")
ext = [x[0], x[-1], z[0], z[-1]]


def mark(ax):
    for zz in (12.5, -12.5):
        ax.axhline(zz, ls="--", lw=1.1, color="cyan")
    for zz in (25, -25):
        ax.axhline(zz, ls="--", lw=1.0, color="lime")
    for zz in (35, -35):
        ax.axhline(zz, ls=":", lw=1.0, color="0.5")


fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.7), sharey=True)

# (left) log — small free-region density is visually exaggerated
imL = axL.imshow(np.log10(np.maximum(sl.T, 1e-8)), origin="lower", aspect="auto",
                 extent=ext, cmap=cmap)
mark(axL)
axL.set_xlabel("x (Bohr)"); axL.set_ylabel("z (Bohr)")
axL.set_title("log₁₀ n — exaggerates the (near-empty) free region", fontsize=9)
fig.colorbar(imL, ax=axL, label="log₁₀ n")

# (right) linear — true contrast: density is confined to the slab
imR = axR.imshow(sl.T, origin="lower", aspect="auto", extent=ext, cmap=cmap,
                 vmin=0.0, vmax=float(sl.max()))
mark(axR)
axR.set_xlabel("x (Bohr)")
axR.set_title("linear n — density confined to |z|<12.5, free region ≈ 0", fontsize=9)
fig.colorbar(imR, ax=axR, label="n (a₀⁻³)")

fig.suptitle("GS slab density (xz) — 82 e, r$_s$≈5.67, box 50×50×70  ·  "
             "cyan=slab faces ±12.5 · lime=CAP edges ±25 · grey=box ±35", fontsize=9)
fig.tight_layout()
fig.savefig(f"{HERE}/gs_density_xz.png", dpi=200)
plt.close(fig)
print(f"wrote {HERE}/gs_density_xz.png  (log | linear)")
print(f"free-region (12.5<|z|<25) mean density = "
      f"{sl[:, (np.abs(z)>12.5)&(np.abs(z)<25)].mean():.2e} (cf interior {n_int.mean():.2e})")
