#!/usr/bin/env python3
"""Post-process the Phase-3 WP-projectile run: xz density gif + response traces."""
import glob, os, re
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

BASE = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/02_projectile_slab/wp_slab/results")
OUT = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "hypotheses/02_projectile_slab")
os.makedirs(OUT, exist_ok=True)
HALF = 12.5
DT, WE = 0.02, 10


def load_vti(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return np.fft.fftshift(a)


# ---- response traces (dipole / energy) ------------------------------------
csv = os.path.join(BASE, "energy_dipole_vs_time.csv")
if os.path.exists(csv):
    d = np.genfromtxt(csv, delimiter=",", names=True)
    t = d["time_au"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    a1.plot(t, (d["total_ha"] - d["total_ha"][0]) * 1e3, lw=1.3)
    a1.set_ylabel(r"$E(t)-E(0)$ (mHa)"); a1.set_title("WP through slab: system response")
    a2.plot(t, d["dipole_z"], lw=1.3, color="C1")
    a2.set_ylabel(r"dipole$_z$ (a.u.)"); a2.set_xlabel("time (a.u.)")
    for ax in (a1, a2): ax.axhline(0, ls=":", color="k", lw=0.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "wp_response.png"), dpi=140)
    print("wrote wp_response.png  | final dE =",
          f"{(d['total_ha'][-1]-d['total_ha'][0])*1e3:.3f} mHa")

# ---- xz density gif --------------------------------------------------------
frames = sorted(glob.glob(os.path.join(BASE, "density_frames", "density_*.vti")),
                key=lambda p: int(re.search(r"(\d+)\.vti", p).group(1)))
if frames:
    vols = [load_vti(f) for f in frames]
    nz, ny, nx = vols[0].shape; L = 50.0
    xz = [v[:, ny // 2, :] for v in vols]
    vmax = np.percentile(np.stack(xz), 99.5)           # WP blob visible vs bath
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(xz[0], origin="lower", extent=[-L/2, L/2, -L/2, L/2],
                   vmin=0.0, vmax=vmax, cmap="inferno", aspect="equal")
    ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
    fig.colorbar(im, ax=ax, fraction=0.046).set_label(r"n (a$_0^{-3}$)")
    ttl = ax.set_title("WP→slab  t=0.00 a.u.")
    for s in (-HALF, HALF): ax.axhline(s, ls=":", color="cyan", lw=0.7)

    def upd(i):
        im.set_data(xz[i]); ttl.set_text(f"WP→slab  t={i*WE*DT:.2f} a.u."); return im, ttl
    FuncAnimation(fig, upd, frames=len(xz), blit=False).save(
        os.path.join(OUT, "wp_xz_density.gif"), writer=PillowWriter(fps=12))
    print("wrote wp_xz_density.gif", f"({len(xz)} frames)")
