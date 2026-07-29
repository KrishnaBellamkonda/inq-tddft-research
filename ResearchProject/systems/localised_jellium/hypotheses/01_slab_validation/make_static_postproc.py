#!/usr/bin/env python3
"""Post-process the Phase-2 static 2 au run.

  * energy_conservation.png + verdict (T3.4: static background ⇒ |ΔE| ≈ 0)
  * static_xz_density.gif — xz density slice vs time, fixed colorbar (shared clim)
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
        "scripts/01_slab_validation/static_2au/results")
OUT = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "hypotheses/01_slab_validation")
N0 = 234.0 / (50.0 * 50.0 * 25.0)


def load_vti(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update()
    img = r.GetOutput()
    nx, ny, nz = img.GetDimensions()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return np.fft.fftshift(a)          # INQ FFT-natural order → centre the box


# ---- energy conservation (T3.4) -------------------------------------------
ev = os.path.join(BASE, "energy_vs_time.csv")
if os.path.exists(ev):
    d = np.genfromtxt(ev, delimiter=",", names=True)
    t, E = np.atleast_1d(d["time_au"]), np.atleast_1d(d["total_ha"])
    dE = E - E[0]
    drift = np.abs(dE).max()
    print(f"T3.4 energy drift over 2 au: max|dE| = {drift:.3e} Ha "
          f"({drift/abs(E[0])*1e6:.2f} ppm)")
    print(f"T3.4 VERDICT = {'PASS' if drift < 1e-3 else 'REVIEW'}")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, dE * 1e3, lw=1.5)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$E(t)-E(0)$ (mHa)")
    ax.set_title("Static run: total-energy conservation (T3.4)")
    ax.axhline(0, ls=":", color="k", lw=0.6)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "energy_conservation.png"), dpi=140)
    print("wrote energy_conservation.png")

# ---- xz density gif (fixed clim) ------------------------------------------
frames = sorted(glob.glob(os.path.join(BASE, "density_frames", "density_*.vti")),
                key=lambda p: int(re.search(r"(\d+)\.vti", p).group(1)))
if frames:
    vols = [load_vti(f) for f in frames]
    nz, ny, nx = vols[0].shape
    L = 50.0
    xz = [v[:, ny // 2, :] for v in vols]      # slice at y=0 → (z, x)
    vmax = max(s.max() for s in xz)
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(xz[0], origin="lower", extent=[-L/2, L/2, -L/2, L/2],
                   vmin=0.0, vmax=vmax, cmap="inferno", aspect="equal")
    ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
    cb = fig.colorbar(im, ax=ax, fraction=0.046); cb.set_label(r"n (a$_0^{-3}$)")
    ttl = ax.set_title("static  t=0.00 a.u.")
    for s in (-12.5, 12.5):
        ax.axhline(s, ls=":", color="cyan", lw=0.7)

    def upd(i):
        im.set_data(xz[i]); ttl.set_text(f"static  t={i*2*0.02:.2f} a.u."); return im, ttl
    anim = FuncAnimation(fig, upd, frames=len(xz), blit=False)
    gif = os.path.join(OUT, "static_xz_density.gif")
    anim.save(gif, writer=PillowWriter(fps=12)); print("wrote", gif, f"({len(xz)} frames)")
