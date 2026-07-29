#!/usr/bin/env python3
"""Shared density-view renderer for the localised-jellium notebooks.

Two renderers, both through the canonical `inqview.load_vti` loader (physical
order, NO fftshift — fixes the "slab-at-the-edges" bug):

* `render_total_views` — the 3 TOTAL-density views (n, Δ-first, Δ-prev). Used by
  the current-data (minimal-run) figures, which only had `density::total`.
* `render_decomposition_views` — the full THREE-WAY decomposition the full-suite
  re-run unlocked: {total, bath, wp} × {n(t), Δ-vs-first, Δ-vs-prev}. The bath is
  derived run-independently as `bath(t) = n_total(t) − |ψ_WP(t)|²` (CONTEXT.md
  "Density decomposition"). Classical runs (no `density_wp`) fall back to total.
"""
from __future__ import annotations

import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.animation import FuncAnimation, PillowWriter

from inqview.visualisation import load_vti  # canonical loader


def _readable_cbar(fig, im, ax, label, *, sci=True):
    """Colorbar with READABLE ticks: scientific offset (×10^n shown once at top) +
    2 significant figures, so long decimals don't get clipped at the figure edge.
    `sci=False` for log/symlog axes (matplotlib's own 10^n ticks are already short)."""
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    if sci:
        fmt = mticker.ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((0, 0))                 # force the shared ×10^n offset
        cb.ax.yaxis.set_major_formatter(fmt)
        cb.ax.yaxis.set_major_locator(mticker.MaxNLocator(5))
        cb.ax.yaxis.get_offset_text().set_fontsize(7)
    cb.ax.tick_params(labelsize=7)
    cb.set_label(label, fontsize=8)
    return cb


def _sorted_frames(frames_dir):
    paths = glob.glob(os.path.join(frames_dir, "*.vti"))
    return sorted(paths, key=lambda p: int(re.search(r"(\d+)\.vti$", p).group(1)))


def render_total_views(frames_dir, out_dir, prefix, *, dt, write_every,
                       slab_half=None, expect_centered=False, y_plane=0.0,
                       fps=12, label="total"):
    """Render the 3 total-density xz GIFs. Returns {view: path}."""
    paths = _sorted_frames(frames_dir)
    if not paths:
        return {}
    # First frame self-checks the mapping when the feature is centred (GS slab).
    f0 = load_vti(paths[0], expect_centered_axis="z" if expect_centered else None)
    fields = [f0] + [load_vti(p) for p in paths[1:]]
    x, z = f0.x, f0.z
    extent = [x[0], x[-1], z[0], z[-1]]                 # physical coords from axes
    xz = [f.xz_slice(y_plane) for f in fields]          # (nz, nx) rows=z

    series = {
        "total": ("inferno", xz),
        "dfirst": ("RdBu_r", [s - xz[0] for s in xz]),
        "dprev": ("RdBu_r", [np.zeros_like(xz[0])] +
                  [xz[i] - xz[i - 1] for i in range(1, len(xz))]),
    }
    titles = {
        "total": f"{label}  n(t)",
        "dfirst": f"{label}  Δn = n(t)−n(0)",
        "dprev": f"{label}  Δn = n(t)−n(t−Δt)",
    }
    out = {}
    for view, (cmap, frames) in series.items():
        stack = np.stack(frames)
        if view == "total":
            vmin, vmax = 0.0, float(np.percentile(stack, 99.5))
        else:                                            # symmetric diff scale
            a = float(np.percentile(np.abs(stack), 99.5)) or 1e-12
            vmin, vmax = -a, a
        fig, ax = plt.subplots(figsize=(5, 5))
        im = ax.imshow(frames[0], origin="lower", extent=extent,
                       vmin=vmin, vmax=vmax, cmap=cmap, aspect="equal")
        ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
        fig.colorbar(im, ax=ax, fraction=0.046).set_label(r"n (a$_0^{-3}$)")
        if slab_half is not None:
            for s in (-slab_half, slab_half):
                ax.axhline(s, ls=":", color="cyan", lw=0.7)
        ttl = ax.set_title(f"{titles[view]}  t=0.00")

        def upd(i, frames=frames, view=view):
            im.set_data(frames[i])
            ttl.set_text(f"{titles[view]}  t={i*write_every*dt:.2f}")
            return im, ttl

        path = os.path.join(out_dir, f"{prefix}_{view}.gif")
        FuncAnimation(fig, upd, frames=len(frames), blit=False).save(
            path, writer=PillowWriter(fps=fps))
        plt.close(fig)
        out[view] = path
        print(f"wrote {path} ({len(frames)} frames, clim [{vmin:.2e},{vmax:.2e}])")
    return out


# ---------------------------------------------------------------------------
# Full three-way decomposition (total / bath = total − wp / wp) × 3 views.
# ---------------------------------------------------------------------------
def _load_xz_series(frames_dir, y_plane=0.0):
    """Load a VTI series as memory-light xz slices. Returns (steps, x, z, [slice])."""
    paths = _sorted_frames(frames_dir)
    steps, slices, x, z = [], [], None, None
    for p in paths:
        f = load_vti(p)                      # physical order, no fftshift
        if x is None:
            x, z = f.x, f.z
        iy = int(np.argmin(np.abs(f.y - y_plane)))
        slices.append(f.data[:, iy, :].T)    # (nz, nx): rows z, cols x
        steps.append(int(re.search(r"(\d+)\.vti$", p).group(1)))
    return np.array(steps), x, z, slices


def _gif(frames, extent, cmap, vmin, vmax, title_base, steps, dt, write_every,
         out, slab_half, fps, ylabel=r"n (a$_0^{-3}$)", cap_inner=None):
    """One animated figure with TWO panels side by side — LINEAR (left) and LOG
    (right) — sharing fixed colour scales across all frames (the gradations never
    move). Signed fields (Δn, cmap=RdBu_r) use a symmetric-log right panel; positive
    fields (n, inferno) use a true-log right panel."""
    signed = (cmap == "RdBu_r")
    if signed:
        a = max(abs(vmin), abs(vmax)) or 1e-12
        lin_norm = mcolors.Normalize(-a, a)
        log_norm = mcolors.SymLogNorm(linthresh=a / 100.0, vmin=-a, vmax=a, base=10)
    else:
        a = max(vmax, 1e-12)
        lin_norm = mcolors.Normalize(0.0, a)
        log_norm = mcolors.LogNorm(vmin=a / 1e3, vmax=a)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.4, 4.7))
    imL = axL.imshow(frames[0], origin="lower", extent=extent, norm=lin_norm,
                     cmap=cmap, aspect="equal")
    imR = axR.imshow(frames[0], origin="lower", extent=extent, norm=log_norm,
                     cmap=cmap, aspect="equal")
    for ax, pl in ((axL, "linear"), (axR, "log")):
        ax.set_xlabel("x (Bohr)")
        ax.text(0.03, 0.97, pl, transform=ax.transAxes, va="top", ha="left",
                fontsize=8, color="w",
                bbox=dict(fc="k", alpha=0.4, lw=0, pad=1.5))
        if slab_half is not None:
            for s in (-slab_half, slab_half):
                ax.axhline(s, ls=":", color="cyan", lw=0.7)
        if cap_inner is not None:                   # CAP boundaries (dashed)
            for s in (-abs(cap_inner), abs(cap_inner)):
                ax.axhline(s, ls="--", color="lime", lw=1.0)
    axL.set_ylabel("z (Bohr)")
    _readable_cbar(fig, imL, axL, ylabel, sci=True)     # ×10^n offset, 2 s.f.
    _readable_cbar(fig, imR, axR, ylabel, sci=False)    # log decades
    sup = fig.suptitle(f"{title_base}  t=0.00", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    def upd(i):
        imL.set_data(frames[i]); imR.set_data(frames[i])
        sup.set_text(f"{title_base}  t={i*write_every*dt:.2f}")
        return imL, imR, sup

    FuncAnimation(fig, upd, frames=len(frames), blit=False).save(
        out, writer=PillowWriter(fps=fps))
    plt.close(fig)


def render_decomposition_views(total_dir, wp_dir, out_dir, prefix, *, dt,
                               write_every, slab_half=None, y_plane=0.0, fps=12,
                               cap_inner=None):
    """Render {total, bath, wp} × {n, dfirst, dprev} xz GIFs from full-suite data.

    bath(t) = n_total(t) − |psi_WP(t)|^2 (frame-aligned). If `wp_dir` is absent
    (classical run), only the total row is produced. Returns {(system,view): path}.
    """
    steps, x, z, tot = _load_xz_series(total_dir, y_plane)
    have_wp = bool(wp_dir) and os.path.isdir(wp_dir) and bool(_sorted_frames(wp_dir))
    if have_wp:
        _, _, _, wp = _load_xz_series(wp_dir, y_plane)
        n = min(len(tot), len(wp)); steps = steps[:n]
        tot = tot[:n]; wp = wp[:n]
        bath = [tot[i] - wp[i] for i in range(n)]
        systems = [("total", tot, "total"), ("bath", bath, "bath (gas)"),
                   ("wp", wp, "wavepacket")]
    else:
        systems = [("total", tot, "total")]
    extent = [x[0], x[-1], z[0], z[-1]]
    out = {}
    for key, series, lab in systems:
        views = {
            "total": ("inferno", series),
            "dfirst": ("RdBu_r", [s - series[0] for s in series]),
            "dprev": ("RdBu_r", [np.zeros_like(series[0])] +
                      [series[i] - series[i - 1] for i in range(1, len(series))]),
        }
        labels = {"total": f"{lab}  n(t)",
                  "dfirst": f"{lab}  Δn = n(t)−n(0)",
                  "dprev": f"{lab}  Δn = n(t)−n(t−Δt)"}
        for view, (cmap, frames) in views.items():
            stk = np.stack(frames)
            if view == "total":
                vmin, vmax = 0.0, float(np.percentile(stk, 99.5)) or 1e-12
            else:
                a = float(np.percentile(np.abs(stk), 99.5)) or 1e-12
                vmin, vmax = -a, a
            path = os.path.join(out_dir, f"{prefix}_{key}_{view}.gif")
            _gif(frames, extent, cmap, vmin, vmax, labels[view], steps, dt,
                 write_every, path, slab_half, fps, cap_inner=cap_inner)
            out[(key, view)] = path
            print(f"wrote {path} ({len(frames)} frames)")
    return out


# ---------------------------------------------------------------------------
# Classical run: total e-density (3) + projectile Gaussian charge (3) = 6 GIFs.
# The projectile is an external Ehrenfest ion (no density_wp); its charge is the
# exact UPF Gaussian rigidly centred on the tracked position R_ion(t).
# ---------------------------------------------------------------------------
def render_classical_views(total_dir, track_csv, out_dir, prefix, *, dt,
                           write_every, sigma_pot=0.35, slab_half=None,
                           y_plane=0.0, fps=12, cap_inner=None):
    """Render {total electrons, projectile Gaussian} × {n, dfirst, dprev} = 6 GIFs.

    The projectile slice at y=`y_plane` is the analytic Gaussian
    ρ(x,z;t) = (2πσ²)^{-3/2} · exp(−[(x−x_i)² + (y_p−y_i)² + (z−z_i)²]/2σ²),
    with (x_i,y_i,z_i) = R_ion at the frame's step, read from `track_csv`.
    Returns {(system,view): path}.
    """
    steps, x, z, tot = _load_xz_series(total_dir, y_plane)
    # ion position vs step from the track
    trk = np.genfromtxt(track_csv, delimiter=",", names=True)
    tstep = trk["step"].astype(int)
    def ion_at(s):
        i = int(np.argmin(np.abs(tstep - s)))
        return float(trk["x"][i]), float(trk["y"][i]), float(trk["z"][i])
    ZZ, XX = np.meshgrid(z, x, indexing="ij")          # (nz, nx) matching slices
    s2 = 2.0 * sigma_pot * sigma_pot
    norm = (2.0 * np.pi * sigma_pot * sigma_pot) ** (-1.5)
    proj = []
    for s in steps:
        xi, yi, zi = ion_at(s)
        g = norm * np.exp(-(((XX - xi) ** 2) + (yi - y_plane) ** 2 + ((ZZ - zi) ** 2)) / s2)
        proj.append(g)
    extent = [x[0], x[-1], z[0], z[-1]]
    out = {}
    for key, series, lab in [("total", tot, "total e-density"),
                             ("proj", proj, "projectile (Gaussian)")]:
        views = {
            "total": ("inferno", series),
            "dfirst": ("RdBu_r", [s - series[0] for s in series]),
            "dprev": ("RdBu_r", [np.zeros_like(series[0])] +
                      [series[i] - series[i - 1] for i in range(1, len(series))]),
        }
        labels = {"total": f"{lab}  n(t)",
                  "dfirst": f"{lab}  Δn = n(t)−n(0)",
                  "dprev": f"{lab}  Δn = n(t)−n(t−Δt)"}
        for view, (cmap, frames) in views.items():
            stk = np.stack(frames)
            if view == "total":
                vmin, vmax = 0.0, float(np.percentile(stk, 99.5)) or 1e-12
            else:
                a = float(np.percentile(np.abs(stk), 99.5)) or 1e-12
                vmin, vmax = -a, a
            path = os.path.join(out_dir, f"{prefix}_{key}_{view}.gif")
            _gif(frames, extent, cmap, vmin, vmax, labels[view], steps, dt,
                 write_every, path, slab_half, fps, cap_inner=cap_inner)
            out[(key, view)] = path
            print(f"wrote {path} ({len(frames)} frames)")
    return out
