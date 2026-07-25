#!/usr/bin/env python3
"""Phase-5 classical+CAP post-processing → stopping power S = ΔKE_ion / x.

The clean classical estimate: the Ehrenfest ion decelerates via electronic
drag while crossing the slab. ΔKE_ion between the entrance (z=−12.5) and exit
(z=+12.5) faces, divided by the traversal length x=25 Bohr, is the stopping
power. Cross-checks the WP bath-energy estimate.
"""
import glob, os, re
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

BASE = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/03_cap_stopping/classical_cap/results")
OUT = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "hypotheses/03_cap_stopping")
HALF, DT, WE = 12.5, 0.02, 10


def load_vti(p):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(p); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    return np.fft.fftshift(vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx))


d = np.genfromtxt(os.path.join(BASE, "classical_trace.csv"), delimiter=",", names=True)
t, z, ke = d["time_au"], d["ion_z"], d["ke_ion_ha"]

# KE at the slab entrance/exit faces (interpolate vs ion z, monotonic in z).
ke_in = np.interp(-HALF, z, ke)
ke_out = np.interp(+HALF, z, ke)
dKE = ke_in - ke_out
S = dKE / (2 * HALF)
print(f"classical: KE_ion in={ke_in:.4f} out={ke_out:.4f} Ha  "
      f"ΔKE={dKE:.4f} Ha over x={2*HALF} Bohr")
print(f"STOPPING POWER S = ΔKE/x = {S:.5f} Ha/Bohr = {S*27.211:.4f} eV/Bohr")
print(f"bath over-drain check: num_e {d['num_electrons'][0]:.2f}->{d['num_electrons'][-1]:.2f}")

fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
a1.plot(t, ke, lw=1.4); a1.set_ylabel(r"KE$_{ion}$ (Ha)")
a1.set_title(f"Classical projectile deceleration  (S={S*27.211:.3f} eV/Bohr)")
a2.plot(t, z, lw=1.4, color="C2"); a2.set_ylabel("ion z (Bohr)"); a2.set_xlabel("time (a.u.)")
for s in (-HALF, HALF): a2.axhline(s, ls=":", color="C3", lw=0.7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "classical_stopping.png"), dpi=140)
print("wrote classical_stopping.png")

frames = sorted(glob.glob(os.path.join(BASE, "density_frames", "density_*.vti")),
                key=lambda p: int(re.search(r"(\d+)\.vti", p).group(1)))
if frames:
    vols = [load_vti(f) for f in frames]; nz, ny, nx = vols[0].shape; L = 50.0
    xz = [v[:, ny // 2, :] for v in vols]; vmax = np.percentile(np.stack(xz), 99.5)
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(xz[0], origin="lower", extent=[-L/2, L/2, -L/2, L/2],
                   vmin=0, vmax=vmax, cmap="inferno", aspect="equal")
    ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
    fig.colorbar(im, ax=ax, fraction=0.046).set_label(r"n (a$_0^{-3}$)")
    for s in (-HALF, HALF): ax.axhline(s, ls=":", color="cyan", lw=0.7)
    ttl = ax.set_title("classical+CAP  t=0.00")
    def upd(i):
        im.set_data(xz[i]); ttl.set_text(f"classical+CAP  t={i*WE*DT:.2f}"); return im, ttl
    FuncAnimation(fig, upd, frames=len(xz), blit=False).save(
        os.path.join(OUT, "classical_xz_density.gif"), writer=PillowWriter(fps=12))
    print(f"wrote classical_xz_density.gif ({len(xz)} frames)")
