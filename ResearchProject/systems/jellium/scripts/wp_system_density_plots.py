#!/usr/bin/env python3
"""Absolute WP-only and system-only (bath) density plots for the matched run.

Companion to wp_vs_classical_matched.py.  For the WP run
`run_wp_n162_L50_E100_sigma1_v2` this writes the RAW (not t0-subtracted)
component densities so the n_total - n_wp split can be checked by eye:

  n_wp     = |psi_wp|^2                 (the moving Gaussian electron)
  n_system = n_total - n_wp             (the bath alone)

Both are sampled at the WP's exact density_wp frame times (so n_system is an
exact same-step subtraction).  Sequential colormap (inferno) since these are
non-negative densities; the induced/difference views live in the Dn plots.

Outputs -> docs/presentations/storyline/tasks/wp_vs_classical_matched/
  wp_density_2d.gif        n_wp xz central-y slab, moving Gaussian + centroid
  wp_density_zprofile.gif  n_wp(z) line density vs t
  system_density_2d.gif    n_system xz slab (bath; near-uniform jellium + hole)
  system_density_zprofile.gif  n_system(z) line density vs t
  fig_wp_density_2d.png / fig_system_density_2d.png   report-standard statics
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _gifutil import save_gif_fixed_palette

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.pipeline import wake
from applications.report1 import apply_style, panel_label
from applications.report1._shared_style import ONE_COL_IN, STYLE_CONFIG

JB = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
WP_RUN = f"{JB}/run_wp_n162_L50_E100_sigma1_v2"
OUT = Path("/local/data/public/skcb2/tddft/docs/presentations/storyline/tasks/wp_vs_classical_matched")
DPI = STYLE_CONFIG["save_dpi"]


def load_components():
    """Read n_total and n_wp at every exact (paired) frame; return slabs +
    z-profiles for n_wp and n_system, plus axes, times, centroid."""
    dt = wake.dt_of(WP_RUN)
    wp_frames = {s: p for s, p in wake.frames(WP_RUN, "density_wp")}
    tot_frames = {s: p for s, p in wake.frames(WP_RUN, "density_total")}
    steps = sorted(s for s in wp_frames if s in tot_frames)
    print(f"paired exact frames: {len(steps)}  steps {steps[0]}..{steps[-1]}")
    wp2d, sys2d, wp1d, sys1d, tt, cent = [], [], [], [], [], []
    x = z = None
    for s in steps:
        ntot, origin, spacing = wake._read_vti(tot_frames[s])
        nwp, _, _ = wake._read_vti(wp_frames[s])
        nsys = ntot - nwp
        ny = ntot.shape[1]
        if x is None:
            x = origin[0] + spacing[0] * np.arange(ntot.shape[2])
            z = origin[2] + spacing[2] * np.arange(ntot.shape[0])
            dxdy = spacing[0] * spacing[1]
        wp2d.append(nwp[:, ny // 2, :]); sys2d.append(nsys[:, ny // 2, :])
        wp1d.append(nwp.sum(axis=(1, 2)) * dxdy)
        sys1d.append(nsys.sum(axis=(1, 2)) * dxdy)
        tt.append(s * dt); cent.append(wake.wp_centroid_z(WP_RUN, s * dt))
    return dict(x=x, z=z, dxdy=dxdy,
                wp2d=np.array(wp2d), sys2d=np.array(sys2d),
                wp1d=np.array(wp1d), sys1d=np.array(sys1d),
                tt=np.array(tt), cent=np.array(cent, float))


def _centroid_marker(ax, cent, color="cyan", sym="v"):
    """Dashed vertical line at z=cent + a small marker above the panel."""
    c = cent if np.isfinite(cent) else np.nan
    line = ax.axvline(c, color=color, ls="--", lw=1.0)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    mk, = ax.plot([c], [1.02], marker=sym, ms=7, color=color, mec=color,
                  transform=tr, clip_on=False, zorder=6)
    return (line, mk)


def _set_centroid(handle, cent):
    line, mk = handle
    if np.isfinite(cent):
        line.set_xdata([cent, cent]); mk.set_xdata([cent])


def anim_2d(C, key, fname, title, clo_pct=0.5, chi_pct=99.7):
    x, z = C["x"], C["z"]; data = C[key]
    ext = [z[0], z[-1], x[0], x[-1]]
    vmn = float(np.percentile(data, clo_pct)); vmx = float(np.percentile(data, chi_pct))
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    im = ax.imshow(data[0].T, origin="lower", extent=ext, aspect="auto",
                   cmap="inferno", vmin=vmn, vmax=vmx)
    cm = _centroid_marker(ax, C["cent"][0])
    ax.set_xlabel("z (Bohr)"); ax.set_ylabel("x (Bohr)")
    cb = fig.colorbar(im, ax=ax, pad=0.02); cb.set_label(r"density (e/Bohr$^3$)")
    ttl = ax.set_title("")
    fig.tight_layout()

    def upd(k):
        im.set_data(data[k].T)
        im.set_clim(vmn, vmx)              # FROZEN colour scale (every frame)
        _set_centroid(cm, C["cent"][k])
        ttl.set_text(f"{title}  —  t={C['tt'][k]:.2f} a.u.")
        return [im]

    save_gif_fixed_palette(fig, upd, len(C["tt"]), OUT / fname, duration_ms=140, dpi=90)
    plt.close(fig)
    print(f"  wrote {fname}  (clim {vmn:.2e}..{vmx:.2e}; global palette)")


def anim_zprofile(C, key, fname, title):
    z = C["z"]; data = C[key]
    ymn, ymx = float(data.min()), float(data.max()); pad = 0.05 * (ymx - ymn + 1e-30)
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    (ln,) = ax.plot([], [], color="C0", lw=1.8)
    cm = _centroid_marker(ax, C["cent"][0], color="k")
    cm[0].set_label("WP centroid")
    ax.set_xlim(z[0], z[-1]); ax.set_ylim(ymn - pad, ymx + pad)
    ax.set_xlabel("z (Bohr)"); ax.set_ylabel("line density (e/Bohr)")
    ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper right")
    ttl = ax.set_title("")
    fig.tight_layout()

    def upd(k):
        ln.set_data(z, data[k])
        _set_centroid(cm, C["cent"][k])
        ttl.set_text(f"{title}  —  t={C['tt'][k]:.2f} a.u.")
        return [ln]

    save_gif_fixed_palette(fig, upd, len(C["tt"]), OUT / fname, duration_ms=140, dpi=100)
    plt.close(fig)
    print(f"  wrote {fname}  (global palette)")


def static_2d(C, key, fname, label, kstar, clo_pct=0.5, chi_pct=99.7):
    x, z = C["x"], C["z"]; data = C[key]
    vmn = float(np.percentile(data[kstar], clo_pct)); vmx = float(np.percentile(data[kstar], chi_pct))
    fig, ax = plt.subplots(figsize=ONE_COL_IN)
    ext = [z[0], z[-1], x[0], x[-1]]
    im = ax.imshow(data[kstar].T, origin="lower", extent=ext, aspect="auto",
                   cmap="inferno", vmin=vmn, vmax=vmx)
    if np.isfinite(C["cent"][kstar]):
        ax.axvline(C["cent"][kstar], color="cyan", ls="--", lw=0.9)
    ax.set_xlabel(r"$z$ (Bohr)"); ax.set_ylabel(r"$x$ (Bohr)")
    panel_label(ax, label)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046); cb.set_label(r"density (e/Bohr$^3$)")
    fig.subplots_adjust(left=0.16, right=0.86, bottom=0.16, top=0.97)
    fig.savefig(OUT / fname, dpi=DPI); plt.close(fig)
    print(f"  wrote {fname}  (clim {vmn:.2e}..{vmx:.2e})")


if __name__ == "__main__":
    apply_style()
    C = load_components()
    # sanity: integrals (per-frame) should be ~1 (wp) and ~162 (system)
    dz = C["z"][1] - C["z"][0]
    print(f"  integral n_wp:     t0={C['wp1d'][0].sum()*dz:.3f}  end={C['wp1d'][-1].sum()*dz:.3f}  (expect 1.0)")
    print(f"  integral n_system: t0={C['sys1d'][0].sum()*dz:.3f}  end={C['sys1d'][-1].sum()*dz:.3f}  (expect 162)")
    print("[wp-only]")
    anim_2d(C, "wp2d", "wp_density_2d.gif", "WP density  $n_{wp}=|\\psi_{wp}|^2$")
    anim_zprofile(C, "wp1d", "wp_density_zprofile.gif", "WP line density $n_{wp}(z)$")
    print("[system-only / bath]")
    anim_2d(C, "sys2d", "system_density_2d.gif", "system (bath) density  $n_{sys}=n_{tot}-n_{wp}$")
    anim_zprofile(C, "sys1d", "system_density_zprofile.gif", "system line density $n_{sys}(z)$")
    print("[statics]")
    kstar = int(0.6 * (len(C["tt"]) - 1))
    static_2d(C, "wp2d", "fig_wp_density_2d.png", "(a)", kstar)
    static_2d(C, "sys2d", "fig_system_density_2d.png", "(b)", kstar)
    print(f"DONE -> {OUT}")
