#!/usr/bin/env python3
"""T-d: σ-sweep bath-only wake at E=100 eV, σ={0.5,1,3,8}.

Bath line density = (density_total − density_wp)(z) summed over x,y (e/Bohr) — WP orbital
removed (same total−wp method as tb_wake_comparability / T-e). Induced wake = bath(z,t) −
bath(z,0). Same v₀ across all σ (same energy), so this isolates the WP-WIDTH dependence
of the jellium wake.

Outputs (this dir):
  sigma_sweep_bath_wake_overlay.png  static: induced bath z-profile at matched times,
                                     one colour per σ, SHARED y-axis.
  sigma_sweep_bath_wake.gif          animated: induced bath z-profiles over a common time
                                     grid; y-limits FIXED ONCE (global min/max) — fixes
                                     the per-frame-scale bug of density_difference_compact.gif.

Known-case (printed): induced bath at t0 == 0 for every σ.
"""
import glob, re, sys, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.postprocess import wake as _wake   # exact n_total-n_wp subtraction

JB = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
THIS = Path(JB) / "_compare_sigma_sweep"
RUNS = {
    0.5: f"{JB}/run_wp_n162_L50_E100_sigma0p5_wf",
    1.0: f"{JB}/run_wp_n162_L50_E100_sigma1_v2",
    3.0: f"{JB}/run_wp_n162_L50_E100_sigma3_wf",
    8.0: f"{JB}/run_wp_n162_L50_E100_sigma8_wf",
}
TARGETS = [1.2, 2.4, 3.6, 4.5]   # within σ0.5's 4.8 a.u. cap (N_STEPS=240 self-spread cap)
N_GRID = 36


def dt_of(run):
    txt = (Path(run) / "results/run_summary.txt").read_text()
    m = re.search(r"dt_au\s*=\s*([\d.]+)", txt)
    return float(m.group(1)) if m else 0.02


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


def bath_z(run, dt, t):
    # EXACT n_total-n_wp subtraction (snaps to a frame with an exact density_wp
    # partner) -> no moving-WP residual. Returns (z, bath_line, t_au).
    return _wake.bath_line_z(run, t)


sigmas = sorted(RUNS)
dts = {s: dt_of(RUNS[s]) for s in sigmas}
# common time range = min over σ of max available (wp-limited) time
tmax = {}
for s in sigmas:
    fw = frames(RUNS[s], "density_wp"); tmax[s] = max(st for st, _ in fw) * dts[s]
T_COMMON = min(tmax.values())
print(f"per-σ max wp time: {{ {', '.join(f'{s}:{tmax[s]:.1f}' for s in sigmas)} }} → common {T_COMMON:.1f} a.u.")

base = {}
for s in sigmas:
    z, b0, _ = bath_z(RUNS[s], dts[s], 0.0)
    _, b0b, _ = bath_z(RUNS[s], dts[s], 0.0)
    print(f"[known-case] σ={s}: induced bath at t0 max|.|={np.abs(b0b - b0).max():.2e} (==0)")
    base[s] = (z, b0)

colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(sigmas)))

# classical reference (E=100 eV) + per-σ WP centroid markers (user request 2026-06-01)
CLASSICAL = f"{JB}/run_classical_n162_L50_E100_v2"
_clz, _clz0, _ = _wake.bath_line_z(CLASSICAL, 0.0)


def classical_induced(t):
    z, l, _ = _wake.bath_line_z(CLASSICAL, t)
    return z, l - _clz0


# static overlay
fig, axs = plt.subplots(1, len(TARGETS), figsize=(4 * len(TARGETS), 4), sharey=True)
ymax = 0.0; panel = {t: [] for t in TARGETS}; cent = {t: {} for t in TARGETS}
for s, c in zip(sigmas, colors):
    z, b0 = base[s]
    for t in TARGETS:
        if t > T_COMMON:
            continue
        zz, bw, tt = bath_z(RUNS[s], dts[s], t)
        ind = bw - b0; panel[t].append((s, zz, ind, c)); ymax = max(ymax, np.abs(ind).max())
        cent[t][s] = _wake.wp_centroid_z(RUNS[s], t)
for ax, t in zip(axs, TARGETS):
    for s, zz, ind, c in panel[t]:
        ax.plot(zz, ind, color=c, lw=1.4, label=f"σ={s}")
        cz = cent[t].get(s)
        if cz is not None and np.isfinite(cz):
            ax.axvline(cz, color=c, ls=":", lw=0.9)           # per-σ WP centroid
    if t <= T_COMMON:
        zc, indc = classical_induced(t)
        ax.plot(zc, indc, color="k", lw=1.6, ls="--", label="classical")  # classical ref
    ax.axhline(0, color="0.6", lw=0.6); ax.set_title(f"t≈{t:.1f} a.u.")
    ax.set_xlabel("z (Bohr)"); ax.grid(alpha=0.3); ax.set_ylim(-1.05 * ymax, 1.05 * ymax)
axs[0].set_ylabel("induced bath Δn(z) (e/Bohr)")
axs[-1].legend(fontsize=7, title="E=100 eV (dotted = WP centroid)")
fig.suptitle("T-d: σ-sweep bath-only induced wake (total−wp), E=100 eV — shared axes")
fig.tight_layout(); fig.savefig(THIS / "sigma_sweep_bath_wake_overlay.png", dpi=150); plt.close(fig)
print(f"wrote sigma_sweep_bath_wake_overlay.png (shared ymax={ymax:.3e})")

# animation with global fixed y-limits
tgrid = np.linspace(0.3, T_COMMON, N_GRID)
curves = {s: [] for s in sigmas}; cent_g = {s: [] for s in sigmas}; gmax = 0.0
cl_curves = []
for s in sigmas:
    z, b0 = base[s]
    for t in tgrid:
        zz, bw, _ = bath_z(RUNS[s], dts[s], t); ind = bw - b0
        curves[s].append((zz, ind)); gmax = max(gmax, np.abs(ind).max())
        cent_g[s].append(_wake.wp_centroid_z(RUNS[s], t))
for t in tgrid:
    cl_curves.append(classical_induced(t))
ylim = 1.05 * gmax
print(f"animation global fixed ylim = ±{ylim:.3e} (NOT per-frame)")

figA, axA = plt.subplots(figsize=(7, 5)); lines = {}; cmarks = {}
for s, c in zip(sigmas, colors):
    (ln,) = axA.plot([], [], color=c, lw=1.8, label=f"σ={s}"); lines[s] = ln
    cmarks[s] = axA.axvline(np.nan, color=c, ls=":", lw=0.9)
(clln,) = axA.plot([], [], color="k", ls="--", lw=1.6, label="classical")
axA.axhline(0, color="0.6", lw=0.6)
axA.set_xlim(base[1.0][0][0], base[1.0][0][-1]); axA.set_ylim(-ylim, ylim)
axA.set_xlabel("z (Bohr)"); axA.set_ylabel("induced bath Δn(z) (e/Bohr)")
axA.legend(fontsize=8, title="E=100 eV (dotted=centroid)", loc="upper left"); axA.grid(alpha=0.3)
ttl = axA.set_title("")


def update(k):
    for s in sigmas:
        zz, ind = curves[s][k]; lines[s].set_data(zz, ind)
        cz = cent_g[s][k]
        if cz is not None and np.isfinite(cz):
            cmarks[s].set_xdata([cz, cz])
    zc, indc = cl_curves[k]; clln.set_data(zc, indc)
    ttl.set_text(f"T-d σ-sweep bath wake (+classical)  t={tgrid[k]:.2f} a.u.  (fixed scale)")
    return list(lines.values()) + list(cmarks.values()) + [clln, ttl]


anim = animation.FuncAnimation(figA, update, frames=N_GRID, interval=120, blit=False)
anim.save(THIS / "sigma_sweep_bath_wake.gif", writer="pillow", dpi=100); plt.close(figA)
print("wrote sigma_sweep_bath_wake.gif")
