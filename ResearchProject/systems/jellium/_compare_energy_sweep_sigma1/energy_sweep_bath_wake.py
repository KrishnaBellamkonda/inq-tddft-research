#!/usr/bin/env python3
"""T-e: energy-sweep bath-only wake, σ=1, E={20,25,50,100,200,300} eV.

Bath line density = (density_total − density_wp)(z) summed over x,y (e/Bohr); the WP
orbital is removed so the profile is the pure jellium response (same total−wp method as
tb_wake_comparability_E20.py, verified there). Induced wake = bath(z,t) − bath(z,0).

Outputs (this dir):
  energy_sweep_bath_wake_overlay.png  static: induced bath z-profile at matched times,
                                      one colour per energy, SHARED y-axis.
  energy_sweep_bath_wake.gif          animated: induced bath z-profiles vs z over a
                                      common time grid; y-limits FIXED ONCE (global
                                      min/max over all frames+energies) — this is the
                                      1D analogue of the colour-scale fix (the buggy
                                      density_difference_compact.gif used per-frame
                                      vmin/vmax → an animated scale).

Known-case (printed): induced bath at t0 == 0 for every energy.
"""
import glob, re, sys, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.pipeline import wake as _wake   # exact n_total-n_wp subtraction
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path

JB = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
THIS = Path(JB) / "_compare_energy_sweep_sigma1"
ENERGIES = [20, 25, 50, 100, 200, 300]
DT = 0.01
T_MAX_COMMON = 7.0                 # limited by E300 (last wp frame t=7.0)
N_GRID = 36                        # animation frames on the common time grid
TARGETS = [1.5, 3.0, 4.5, 6.0]     # static-overlay matched times


def zline(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    oz, sz = img.GetOrigin()[2], img.GetSpacing()[2]
    sx, sy = img.GetSpacing()[0], img.GetSpacing()[1]
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return oz + sz * np.arange(nz), a.sum(axis=(1, 2)) * sx * sy


def frames(run, sub):
    out = []
    for f in glob.glob(f"{run}/results/raw/vti/{sub}/*.vti"):
        m = re.search(r"_t(\d+)\.vti$", f)
        if m:
            out.append((int(m.group(1)), f))
    return sorted(out)


def bath_z(run, t):
    """bath line density = total - wp, EXACT-step subtraction (no moving-WP
    residual; snaps to a frame with an exact density_wp partner)."""
    return _wake.bath_line_z(run, t)


# baselines + known-case
base = {}
for E in ENERGIES:
    run = f"{JB}/run_wp_n162_L50_E{E}_sigma1_v2"
    z, b0, _ = bath_z(run, 0.0)
    _, b0b, _ = bath_z(run, 0.0)
    print(f"[known-case] E{E}: induced bath at t0 max|.|={np.abs(b0b - b0).max():.2e} (==0)")
    base[E] = (z, b0, run)

colors = plt.cm.plasma(np.linspace(0.05, 0.9, len(ENERGIES)))

# per-E WP centroid markers (user request 2026-06-01); _wake imported at top.

# ---- static overlay: induced bath z-profile at matched times (shared y) ----
fig, axs = plt.subplots(1, len(TARGETS), figsize=(4 * len(TARGETS), 4), sharey=True)
ymax = 0.0
panel = {t: [] for t in TARGETS}; cent = {t: {} for t in TARGETS}
for E, c in zip(ENERGIES, colors):
    z, b0, run = base[E]
    for t in TARGETS:
        zz, bw, tt = bath_z(run, t)
        ind = bw - b0
        panel[t].append((E, zz, ind, c))
        cent[t][E] = _wake.wp_centroid_z(run, t)
        ymax = max(ymax, np.abs(ind).max())
for ax, t in zip(axs, TARGETS):
    for E, zz, ind, c in panel[t]:
        ax.plot(zz, ind, color=c, lw=1.4, label=f"E={E}")
        cz = cent[t].get(E)
        if cz is not None and np.isfinite(cz):
            ax.axvline(cz, color=c, ls=":", lw=0.8)        # per-E WP centroid
    ax.axhline(0, color="0.6", lw=0.6); ax.set_title(f"t≈{t:.1f} a.u.")
    ax.set_xlabel("z (Bohr)"); ax.grid(alpha=0.3); ax.set_ylim(-1.05 * ymax, 1.05 * ymax)
axs[0].set_ylabel("induced bath Δn(z) (e/Bohr)")
axs[-1].legend(fontsize=7, title="σ=1 (dotted=centroid)")
fig.suptitle("T-e: energy-sweep bath-only induced wake (total−wp), σ=1 — shared axes")
fig.tight_layout()
fig.savefig(THIS / "energy_sweep_bath_wake_overlay.png", dpi=150)
plt.close(fig)
print(f"wrote energy_sweep_bath_wake_overlay.png (shared ymax={ymax:.3e})")

# ---- animation: common time grid, GLOBAL fixed y-limits ----
tgrid = np.linspace(0.3, T_MAX_COMMON, N_GRID)
curves = {E: [] for E in ENERGIES}      # curves[E][k] = (z, induced)
gmax = 0.0
for E in ENERGIES:
    z, b0, run = base[E]
    for t in tgrid:
        zz, bw, _ = bath_z(run, t)
        ind = bw - b0
        curves[E].append((zz, ind))
        gmax = max(gmax, np.abs(ind).max())
ylim = 1.05 * gmax                       # computed ONCE over all frames + energies
print(f"animation global fixed ylim = ±{ylim:.3e} (NOT per-frame)")

figA, axA = plt.subplots(figsize=(7, 5))
lines = {}
for E, c in zip(ENERGIES, colors):
    (ln,) = axA.plot([], [], color=c, lw=1.6, label=f"E={E}")
    lines[E] = ln
axA.axhline(0, color="0.6", lw=0.6)
axA.set_xlim(base[20][0][0], base[20][0][-1]); axA.set_ylim(-ylim, ylim)
axA.set_xlabel("z (Bohr)"); axA.set_ylabel("induced bath Δn(z) (e/Bohr)")
axA.legend(fontsize=8, title="σ=1", loc="upper left"); axA.grid(alpha=0.3)
ttl = axA.set_title("")


def update(k):
    for E in ENERGIES:
        zz, ind = curves[E][k]
        lines[E].set_data(zz, ind)
    ttl.set_text(f"T-e energy-sweep bath wake  t={tgrid[k]:.2f} a.u.  (fixed scale)")
    return list(lines.values()) + [ttl]


anim = animation.FuncAnimation(figA, update, frames=N_GRID, interval=120, blit=False)
anim.save(THIS / "energy_sweep_bath_wake.gif", writer="pillow", dpi=100)
plt.close(figA)
print("wrote energy_sweep_bath_wake.gif")
