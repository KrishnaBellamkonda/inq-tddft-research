"""Bath density difference — figures-2 version (2.948 × 2.948 in).

Slot: 0.48\\textwidth = 2.948 in (square).
Produces annotated-bath-density-difference.png (symlog, more detail).

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 build_bath_diff.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
sys.path.insert(0, str(REPO / "inq-stack/python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator, ScalarFormatter

from inqview.visualisation import style
from inqview import load_vti

style.apply_theme()

HERE  = Path(__file__).parent
RUNS  = REPO / "ResearchProject/systems/jellium"

WP_VTI = RUNS / "run_wp_n162_L50_E100_sigma3_wf" / "results" / "raw" / "vti"
CL_VTI = RUNS / "run_classical_n162_L50_E100"    / "results" / "raw" / "vti"

WP_CZ  = -13.0
WP_K0  = 2.71106
STEP   = 100
DT_STEP = 0.02
T_PHYS = STEP * DT_STEP
Z_WP   = WP_CZ + WP_K0 * T_PHYS
CMAP   = "RdBu_r"


def _fmt_cb(cb) -> None:
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    cb.ax.yaxis.set_major_formatter(fmt)
    cb.ax.yaxis.set_major_locator(MaxNLocator(5))


def _make_panel(delta_zx, xs, zs, norm, *, out: Path) -> None:
    fig = plt.figure(figsize=(2.948, 2.948))
    ax  = fig.add_axes([0.17, 0.14, 0.80, 0.80])

    img = ax.imshow(delta_zx, origin="lower",
                    extent=[zs[0], zs[-1], xs[0], xs[-1]],
                    aspect="equal", norm=norm, cmap=CMAP,
                    interpolation="nearest", rasterized=True)
    cax = ax.inset_axes([0.85, 0.04, 0.04, 0.88])
    cb  = fig.colorbar(img, cax=cax)
    _fmt_cb(cb)
    cb.ax.yaxis.set_tick_params(color="black", labelcolor="black", labelsize=7)
    cb.outline.set_edgecolor("black")
    cb.ax.yaxis.get_offset_text().set_color("black")

    ax.axvline(Z_WP, color="black", lw=0.9, ls="--", label="projectile center")
    ax.legend(loc="upper left", handlelength=1.5, fontsize=7)
    ax.set_xlabel(r"$z\ (\mathrm{Bohr})$")
    ax.set_ylabel(r"$x\ (\mathrm{Bohr})$")

    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {out}")


def main() -> None:
    def _vti(parent: Path, subdir: str) -> Path:
        p = parent / subdir / f"density_t{STEP:06d}.vti"
        if not p.exists():
            raise FileNotFoundError(f"VTI not found: {p}")
        return p

    f_total_wp = load_vti(_vti(WP_VTI, "density_total"))
    f_wp_orb   = load_vti(_vti(WP_VTI, "density_wp"))
    f_total_cl = load_vti(_vti(CL_VTI, "density_total"))

    bath_wp = f_total_wp.data - f_wp_orb.data
    bath_cl = f_total_cl.data
    delta   = bath_cl - bath_wp

    xs, ys, zs = f_total_wp.x, f_total_wp.y, f_total_wp.z
    iy_mid     = int(np.argmin(np.abs(ys)))
    delta_zx   = delta[:, iy_mid, :]

    vmax = float(np.max(np.abs(delta_zx)))
    print(f"|delta|_max = {vmax:.3e} Bohr^-3")

    linthresh = vmax / 100.0
    norm_log  = mcolors.SymLogNorm(linthresh=linthresh, linscale=0.5,
                                   vmin=-vmax, vmax=vmax)
    _make_panel(delta_zx, xs, zs, norm_log,
                out=HERE / "annotated-bath-density-difference.png")


if __name__ == "__main__":
    main()
