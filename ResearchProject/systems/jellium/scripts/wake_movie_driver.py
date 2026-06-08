#!/usr/bin/env python3
"""Reusable WP-vs-classical induced-bath-density wake movies (M7).

For each case (a WP run + the matched classical run) produces, over the WP
run's FULL duration:
  <tag>_wake_2d.gif        3-panel xz Δn maps: [WP | classical | WP-classical]
                           panels 1&2 SHARE one colorbar (shared_clim rule);
                           panel 3 (difference) has its OWN colorbar. linear.
  <tag>_wake_2d_log.gif    same, symlog colour scale.
  <tag>_wake_1d.gif        z-profile Δn(z,t): WP & classical overlaid (top) +
                           WP-classical (bottom). WP centroid marked.
  <tag>_wake_static.png    multi-time static z-profile overlay (WP vs classical).

Bath density n_system = n_total - n_wp (classical: = n_total). Induced
Δn = n_system(t) - n_system(t0). See inqview.postprocess.wake for the
canonical definition + the shared-colorbar rule.

Known-case (printed): Δn at t0 == 0; WP centroid monotonic.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import SymLogNorm
from inqview.postprocess import wake

JB = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
OUT = Path("/local/data/public/skcb2/tddft/docs/presentations/storyline/tasks/batch2_figures")
N_FRAMES = 60           # even sampling over full duration (smooth, sane I/O)
N_STATIC = 4
CLIP_PCT = 99.5         # suppress the WP self-spike / boundary outliers in clim


def _times(run, n):
    """Sample times for a WP run at EXACT density_wp frame times (so the
    n_total - n_wp subtraction is exact — no moving-WP residual). If more than
    n WP frames exist, subsample evenly to n."""
    dt = wake.dt_of(run)
    wpt = wake.wp_frame_times(run, dt)
    if wpt:
        wpt = np.array(sorted(wpt))
        if len(wpt) > n:
            wpt = wpt[np.linspace(0, len(wpt) - 1, n).astype(int)]
        return wpt, float(wpt[-1])
    ft = wake.frames(run, "density_total")
    tmax = max(s for s, _ in ft) * dt
    return np.linspace(0.0, tmax, n), tmax


def _slab_line(run, t):
    """One VTI read -> (xz central-y slab, z-line, centroid, t_au, x, z)."""
    nsys, origin, spacing, t_au, _ = wake.bath_volume(run, t)
    nz, ny, nx = nsys.shape
    slab = nsys[:, ny // 2, :]
    line = nsys.sum(axis=(1, 2)) * spacing[0] * spacing[1]
    x = origin[0] + spacing[0] * np.arange(nx)
    z = origin[2] + spacing[2] * np.arange(nz)
    cent = wake.wp_centroid_z(run, t)
    return slab, line, cent, t_au, x, z


def cache_case(wp_run, cl_run, n=N_FRAMES):
    """Read each VTI ONCE; return dict of cached induced slabs/profiles + axes."""
    tw, tmax = _times(wp_run, n)
    s0, l0w, _, _, x, z = _slab_line(wp_run, 0.0)
    c0, l0c, _, _, _, _ = _slab_line(cl_run, 0.0)
    wp2d, cl2d, wp1d, cl1d, cent, tt = [], [], [], [], [], []
    for t in tw:
        sw, lw, cw, ta, _, _ = _slab_line(wp_run, t)
        sc, lc, _, _, _, _ = _slab_line(cl_run, t)
        wp2d.append(sw - s0); cl2d.append(sc - c0)
        wp1d.append(lw - l0w); cl1d.append(lc - l0c)
        cent.append(cw); tt.append(ta)
    return dict(x=x, z=z, z1=z, tt=np.array(tt), cent=np.array(cent, dtype=float),
                wp2d=np.array(wp2d), cl2d=np.array(cl2d),
                wp1d=np.array(wp1d), cl1d=np.array(cl1d), tmax=tmax)


def animate_2d(C, tag, title, log=False):
    x, z = C["x"], C["z"]
    ext = [z[0], z[-1], x[0], x[-1]]              # imshow: rows=x? we transpose -> z horiz
    # panels 1,2 share clim; panel 3 own clim
    vmn, vmx = wake.shared_clim(C["wp2d"], C["cl2d"], pct=CLIP_PCT)
    dmn, dmx = wake.shared_clim(C["wp2d"] - C["cl2d"], pct=CLIP_PCT)
    if log:
        lin = max(vmx * 1e-3, 1e-6); dlin = max(dmx * 1e-3, 1e-6)
        norm12 = SymLogNorm(linthresh=lin, vmin=vmn, vmax=vmx, base=10)
        norm3 = SymLogNorm(linthresh=dlin, vmin=dmn, vmax=dmx, base=10)
    else:
        norm12 = norm3 = None
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)

    def imshow(ax, slab, norm, vmn, vmx, ttl):
        im = ax.imshow(slab.T, origin="lower", extent=ext, aspect="auto",
                       cmap="RdBu_r", norm=norm,
                       vmin=None if norm else vmn, vmax=None if norm else vmx)
        ax.set_title(ttl); ax.set_xlabel("z (Bohr)")
        return im

    im0 = imshow(axs[0], C["wp2d"][0], norm12, vmn, vmx, "WP  Δn_system")
    im1 = imshow(axs[1], C["cl2d"][0], norm12, vmn, vmx, "classical  Δn_system")
    im2 = imshow(axs[2], C["wp2d"][0] - C["cl2d"][0], norm3, dmn, dmx, "WP − classical")
    axs[0].set_ylabel("x (Bohr)")
    fig.colorbar(im1, ax=axs[:2], shrink=0.8, label="Δn (e/Bohr³)  [shared]")
    fig.colorbar(im2, ax=axs[2], shrink=0.8, label="Δn diff")
    cmark = axs[0].axvline(C["cent"][0], color="k", ls="--", lw=1.0) if np.isfinite(C["cent"][0]) else None
    sup = fig.suptitle("")

    def upd(k):
        im0.set_data(C["wp2d"][k].T); im1.set_data(C["cl2d"][k].T)
        im2.set_data((C["wp2d"][k] - C["cl2d"][k]).T)
        if cmark is not None and np.isfinite(C["cent"][k]):
            cmark.set_xdata([C["cent"][k], C["cent"][k]])
        sup.set_text(f"{title}  —  t={C['tt'][k]:.2f} a.u.{'  [symlog]' if log else ''}")
        return [im0, im1, im2]

    anim = animation.FuncAnimation(fig, upd, frames=len(C["tt"]), interval=120, blit=False)
    name = f"{tag}_wake_2d{'_log' if log else ''}.gif"
    anim.save(OUT / name, writer="pillow", dpi=90); plt.close(fig)
    print(f"  wrote {name}  (shared clim ±{vmx:.2e}, diff ±{dmx:.2e})")


def animate_1d(C, tag, title):
    z = C["z1"]
    ymax = 1.05 * max(np.abs(C["wp1d"]).max(), np.abs(C["cl1d"]).max())
    dmax = 1.05 * np.abs(C["wp1d"] - C["cl1d"]).max()
    fig, (axt, axb) = plt.subplots(2, 1, figsize=(7.5, 6.2), sharex=True,
                                   gridspec_kw=dict(height_ratios=[2, 1]))
    (lw,) = axt.plot([], [], color="C0", lw=1.8, label="WP (total−wp)")
    (lc,) = axt.plot([], [], color="C3", lw=1.8, label="classical")
    (ld,) = axb.plot([], [], color="C2", lw=1.8, label="WP − classical")
    cmark = axt.axvline(np.nan, color="k", ls="--", lw=1.0, label="WP centroid")
    for ax in (axt, axb):
        ax.axhline(0, color="0.6", lw=0.6); ax.grid(alpha=0.3)
        ax.set_xlim(z[0], z[-1])
    axt.set_ylim(-ymax, ymax); axb.set_ylim(-dmax, dmax)
    axt.set_ylabel("induced Δn(z) (e/Bohr)"); axb.set_ylabel("difference")
    axb.set_xlabel("z (Bohr)"); axt.legend(fontsize=8, loc="upper left")
    ttl = axt.set_title("")

    def upd(k):
        lw.set_data(z, C["wp1d"][k]); lc.set_data(z, C["cl1d"][k])
        ld.set_data(z, C["wp1d"][k] - C["cl1d"][k])
        if np.isfinite(C["cent"][k]):
            cmark.set_xdata([C["cent"][k], C["cent"][k]])
        ttl.set_text(f"{title}  —  t={C['tt'][k]:.2f} a.u.  (fixed scale)")
        return [lw, lc, ld, cmark]

    anim = animation.FuncAnimation(fig, upd, frames=len(C["tt"]), interval=120, blit=False)
    anim.save(OUT / f"{tag}_wake_1d.gif", writer="pillow", dpi=100); plt.close(fig)
    print(f"  wrote {tag}_wake_1d.gif  (z-profile ±{ymax:.2e}, diff ±{dmax:.2e})")


def static_overlay(C, tag, title):
    z = C["z1"]; idx = np.linspace(0, len(C["tt"]) - 1, N_STATIC).astype(int)
    ymax = 1.05 * max(np.abs(C["wp1d"][idx]).max(), np.abs(C["cl1d"][idx]).max())
    fig, axs = plt.subplots(1, N_STATIC, figsize=(4 * N_STATIC, 4), sharey=True)
    for ax, k in zip(axs, idx):
        ax.plot(z, C["wp1d"][k], color="C0", lw=1.5, label="WP")
        ax.plot(z, C["cl1d"][k], color="C3", lw=1.5, label="classical")
        if np.isfinite(C["cent"][k]):
            ax.axvline(C["cent"][k], color="k", ls="--", lw=0.9)
        ax.axhline(0, color="0.6", lw=0.6); ax.grid(alpha=0.3)
        ax.set_title(f"t≈{C['tt'][k]:.2f} a.u."); ax.set_xlabel("z (Bohr)")
        ax.set_ylim(-ymax, ymax)
    axs[0].set_ylabel("induced Δn(z) (e/Bohr)"); axs[-1].legend(fontsize=8)
    fig.suptitle(title); fig.tight_layout()
    fig.savefig(OUT / f"{tag}_wake_static.png", dpi=150); plt.close(fig)
    print(f"  wrote {tag}_wake_static.png")


def run_case(wp_run, cl_run, tag, title):
    print(f"[{tag}] caching {Path(wp_run).name} vs {Path(cl_run).name} ...")
    C = cache_case(wp_run, cl_run)
    print(f"  [KC] Δn(t0) max|.| wp={np.abs(C['wp1d'][0]).max():.2e} cl={np.abs(C['cl1d'][0]).max():.2e} (==0); "
          f"centroid {C['cent'][0]:.1f}->{C['cent'][-1]:.1f}; tmax={C['tmax']:.1f}")
    animate_2d(C, tag, title, log=False)
    animate_2d(C, tag, title, log=True)
    animate_1d(C, tag, title)
    static_overlay(C, tag, title)


SIGMA = {
    "0p5": f"{JB}/run_wp_n162_L50_E100_sigma0p5_wf",
    "1":   f"{JB}/run_wp_n162_L50_E100_sigma1_v2",
    "3":   f"{JB}/run_wp_n162_L50_E100_sigma3_wf",
    "8":   f"{JB}/run_wp_n162_L50_E100_sigma8_wf",
}
CL_E100 = f"{JB}/run_classical_n162_L50_E100_v2"

# Energy sweep @ sigma=1 (v2 WP) vs classical at matched energy.
# Classical L50 runs exist at 20/25/50/100/300 (NONE at 200 -> omitted).
ENERGY = {
    "20":  (f"{JB}/run_wp_n162_L50_E20_sigma1_v2",  f"{JB}/run_classical_n162_L50_E20"),
    "25":  (f"{JB}/run_wp_n162_L50_E25_sigma1_v2",  f"{JB}/run_classical_n162_L50_E25"),
    "50":  (f"{JB}/run_wp_n162_L50_E50_sigma1_v2",  f"{JB}/run_classical_n162_L50_E50_v2"),
    "100": (f"{JB}/run_wp_n162_L50_E100_sigma1_v2", f"{JB}/run_classical_n162_L50_E100_v2"),
    "300": (f"{JB}/run_wp_n162_L50_E300_sigma1_v2", f"{JB}/run_classical_n162_L50_E300_v2"),
}

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sigma"
    if mode in ("sigma", "all"):
        sel = sys.argv[2] if len(sys.argv) > 2 else "all"
        cases = SIGMA.items() if sel == "all" else [(sel, SIGMA[sel])]
        for s, wp_run in cases:
            run_case(wp_run, CL_E100, f"sigma{s}_E100",
                     f"M7 σ={s} @E100: induced bath wake  WP vs classical")
    elif mode == "energy":
        sel = sys.argv[2] if len(sys.argv) > 2 else "all"
        items = ENERGY.items() if sel == "all" else [(sel, ENERGY[sel])]
        for e, (wp_run, cl_run) in items:
            run_case(wp_run, cl_run, f"E{e}_sigma1",
                     f"M7 E={e} eV @σ=1: induced bath wake  WP vs classical")
    elif mode in SIGMA:               # back-compat: bare sigma key
        run_case(SIGMA[mode], CL_E100, f"sigma{mode}_E100",
                 f"M7 σ={mode} @E100: induced bath wake  WP vs classical")
