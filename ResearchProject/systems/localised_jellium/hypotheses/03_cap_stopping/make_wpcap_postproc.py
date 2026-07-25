#!/usr/bin/env python3
"""Phase-5 WP+CAP post-processing: absorbed-norm + energy traces, xz gif."""
import glob, os, re
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

BASE = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/03_cap_stopping/wp_cap/results")
OUT = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "hypotheses/03_cap_stopping")
os.makedirs(OUT, exist_ok=True)
HALF, DT, WE = 12.5, 0.02, 10
N_BATH = 234.0


def load_vti(p):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(p); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    return np.fft.fftshift(vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx))


d = np.genfromtxt(os.path.join(BASE, "cap_trace.csv"), delimiter=",", names=True)
t, E, ne = d["time_au"], d["total_ha"], d["num_electrons"]
absorbed = ne[0] - ne                                   # cumulative absorbed norm
print(f"WP+CAP: norm {ne[0]:.3f}->{ne[-1]:.3f}  absorbed={absorbed[-1]:.3f} "
      f"of WP(1.0); bath intact={'yes' if ne[-1] > N_BATH - 0.05 else 'NO (over-drain!)'}")
print(f"  energy removed by CAP = {E[0]-E[-1]:.3f} Ha")

fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
a1.plot(t, absorbed, lw=1.4); a1.set_ylabel("cumulative absorbed norm")
a1.set_title("WP+CAP: CAP absorption & energy")
a1.axhline(1.0, ls=":", color="k", lw=0.6, label="full WP")
a1.legend(frameon=False, fontsize=8)
a2.plot(t, E - E[0], lw=1.4, color="C3"); a2.set_ylabel(r"$E(t)-E(0)$ (Ha)")
a2.set_xlabel("time (a.u.)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "wpcap_traces.png"), dpi=140)
print("wrote wpcap_traces.png")

frames = sorted(glob.glob(os.path.join(BASE, "density_frames", "density_*.vti")),
                key=lambda p: int(re.search(r"(\d+)\.vti", p).group(1)))
vols = [load_vti(f) for f in frames]
nz, ny, nx = vols[0].shape; L = 50.0
xz = [v[:, ny // 2, :] for v in vols]
vmax = np.percentile(np.stack(xz), 99.5)
fig, ax = plt.subplots(figsize=(5, 5))
im = ax.imshow(xz[0], origin="lower", extent=[-L/2, L/2, -L/2, L/2],
               vmin=0, vmax=vmax, cmap="inferno", aspect="equal")
ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
fig.colorbar(im, ax=ax, fraction=0.046).set_label(r"n (a$_0^{-3}$)")
for s in (-HALF, HALF): ax.axhline(s, ls=":", color="cyan", lw=0.7)
for s in (-17.5, 17.5): ax.axhline(s, ls="--", color="lime", lw=0.8)   # CAP edges
ttl = ax.set_title("WP+CAP  t=0.00")

def upd(i):
    im.set_data(xz[i]); ttl.set_text(f"WP+CAP  t={i*WE*DT:.2f} a.u."); return im, ttl
FuncAnimation(fig, upd, frames=len(xz), blit=False).save(
    os.path.join(OUT, "wpcap_xz_density.gif"), writer=PillowWriter(fps=12))
print(f"wrote wpcap_xz_density.gif ({len(xz)} frames)")
