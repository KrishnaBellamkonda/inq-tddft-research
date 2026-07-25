#!/usr/bin/env python3
"""T1 validation for the localised jellium slab GS.

Loads the GS electron density VTI, computes the planar-averaged profile
n(z) = <n(x,y,z)>_{x,y}, and checks:
  * density peaks INSIDE the slab (|z| < half_width), ~0 in vacuum,
  * interior flat to a few % of n0 away from the surface (the "increase R_cl"
    gate, worksheet VC-6),
  * Friedel wavelength near the surface ~ pi/k_F.
Writes n_of_z.png for the notebook and prints a verdict.
"""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VTI = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/01_slab_validation/gs_slab/results/density_gs_system/"
       "density_gs_system.vti")
N0 = 234.0 / (50.0 * 50.0 * 25.0)        # 3.744e-3
HALF = 12.5
RS = (3.0 / (4.0 * np.pi * N0)) ** (1.0 / 3.0)
KF = (9.0 * np.pi / 4.0) ** (1.0 / 3.0) / RS
LAMBDA_FRIEDEL = np.pi / KF

r = vtk.vtkXMLImageDataReader(); r.SetFileName(VTI); r.Update()
img = r.GetOutput()
nx, ny, nz = img.GetDimensions()
sx, sy, sz = img.GetSpacing()
ox, oy, oz = img.GetOrigin()
arr = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)  # z slow

n_of_z = arr.mean(axis=(1, 2))           # planar average
z = oz + sz * np.arange(nz)
# INQ FFT-natural order: if origin is 0, fold to symmetric range about box centre.
if abs(oz) < 1e-9:
    L = sz * nz
    z = ((z + L / 2) % L) - L / 2
    order = np.argsort(z); z = z[order]; n_of_z = n_of_z[order]

interior = np.abs(z) < (HALF - LAMBDA_FRIEDEL)     # away from surface
vac = np.abs(z) > (HALF + 3.0)
mean_int = n_of_z[interior].mean()
dev_int = np.abs(n_of_z[interior] / N0 - 1.0).max()
peak_in_slab = np.abs(z[np.argmax(n_of_z)]) < HALF
vac_level = np.abs(n_of_z[vac]).max() / N0

print(f"n0(target)      = {N0:.5e}  r_s={RS:.3f}  k_F={KF:.3f}  "
      f"lambda_Friedel={LAMBDA_FRIEDEL:.2f} Bohr")
print(f"interior <n>    = {mean_int:.5e}  ({mean_int/N0*100:.1f}% of n0)")
print(f"interior max|dev| = {dev_int*100:.1f}%   (gate: a few %)")
print(f"peak inside slab  = {peak_in_slab}")
print(f"vacuum level    = {vac_level*100:.2f}% of n0")
verdict = (peak_in_slab and dev_int < 0.10 and vac_level < 0.05)
print(f"T1 INTERIOR VERDICT = {'PASS' if verdict else 'REVIEW'}")

# Plot in explicit 10^-3 a0^-3 units so the y-axis carries its own scale and no
# ambiguous offset (1e-5) label is needed (cell-23 TODO).
SCALE = 1.0e-3
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(z, n_of_z / SCALE, lw=1.5, color="C0")
ax.axhline(N0 / SCALE, ls="--", color="k", lw=0.8, label=r"$n_0$")
for s in (-HALF, HALF):
    ax.axvline(s, ls=":", color="C3", lw=0.8)
ax.axvspan(-HALF, HALF, color="C3", alpha=0.07, label="slab background")
ax.ticklabel_format(axis="y", style="plain")
ax.set_xlabel("z (Bohr)")
ax.set_ylabel(r"$\langle n\rangle_{xy}(z)$  ($10^{-3}\,a_0^{-3}$)")
ax.set_title(f"Localised jellium slab GS: planar density  (r$_s$={RS:.2f}, N=234)")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
out = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "hypotheses/01_slab_validation/n_of_z.png")
fig.savefig(out, dpi=140)
print("wrote", out)
