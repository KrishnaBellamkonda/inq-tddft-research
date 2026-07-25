#!/usr/bin/env python3
"""make_setup_figure.py — system-setup figure for the vacuum WP-CAP run.

Per the scientific-figures rule §4, a system-design figure is built from the RUN'S
REAL DENSITY, never a hand-drawn cartoon: we plot the actual t=0 total-density xz
slice (so the wavepacket's starting position is read off the data, not asserted)
and overlay ONLY dashed lines marking the two CAP bands and the WP launch plane.

Vacuum (bulk, no slab): total density + dashed CAP extents + WP launch line.

Usage:
  make_setup_figure.py <results_dir> <out.png> [--cap-inner 20] [--cap-outer 30]
                       [--wp-launch -26]
"""
from __future__ import annotations
import argparse, glob
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from inqview.visualisation.field_io import load_vti
try:
    from inqview.visualisation.style import apply_theme, save_presentation
    apply_theme()
except Exception:
    def save_presentation(fig, p): fig.savefig(p, dpi=600, bbox_inches="tight")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("out")
    ap.add_argument("--cap-inner", type=float, default=30.0)
    ap.add_argument("--cap-outer", type=float, default=40.0)
    ap.add_argument("--wp-launch", type=float, default=-30.0)
    ap.add_argument("--title", default="Vacuum WP-CAP run — setup (t = 0)")
    a = ap.parse_args()

    t0 = sorted(glob.glob(str(Path(a.results_dir) /
                "raw/vti/density_total/*.vti")))[0]
    f = load_vti(t0)
    sl = f.xz_slice(0.0)                      # (nz, nx), rows=z
    ext = [f.x[0], f.x[-1], f.z[0], f.z[-1]]
    vmax = float(np.nanmax(sl))

    fig, ax = plt.subplots(figsize=(4.6, 6.2))
    im = ax.imshow(np.clip(sl, vmax * 1e-4, None), origin="lower", aspect="auto",
                   extent=ext, cmap="viridis",
                   norm=LogNorm(vmin=vmax * 1e-4, vmax=vmax))

    # --- ONE-SIDED +z CAP band + dashed inner boundary ----------------------
    # (a single perturbations::absorbing band; the -z wall is the CAP's wrapped
    #  outer edge via periodicity, marked lightly for context.)
    ax.axhspan(a.cap_inner, a.cap_outer, color="crimson", alpha=0.16, zorder=2)
    ax.axhline(a.cap_inner, ls="--", lw=1.5, color="crimson", zorder=3)
    ax.text(f.x[-1] * 0.92, (a.cap_inner + a.cap_outer) / 2, "CAP\n(+z)",
            color="crimson", ha="right", va="center", fontsize=9, weight="bold")
    # -z wall == wrapped CAP outer edge (periodic): light dotted marker only
    ax.axhline(f.z[0], ls=":", lw=1.0, color="crimson", alpha=0.6, zorder=3)

    # --- WP launch plane ----------------------------------------------------
    ax.axhline(a.wp_launch, ls="--", lw=1.5, color="white", zorder=3)
    ax.text(f.x[0] * 0.9, a.wp_launch + 1.5, f"WP launch  z = {a.wp_launch:g}",
            color="white", ha="left", va="bottom", fontsize=8.5)

    ax.set_xlabel("x (Bohr)")
    ax.set_ylabel("z (Bohr)")
    ax.set_title(a.title, fontsize=10)

    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="5%", pad=0.08)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("n (a₀⁻³, log)", fontsize=8)

    save_presentation(fig, a.out)
    plt.close(fig)
    print(f"[setup] wrote {a.out}  (WP peak read from data, CAP |z|∈"
          f"[{a.cap_inner:g},{a.cap_outer:g}])")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
