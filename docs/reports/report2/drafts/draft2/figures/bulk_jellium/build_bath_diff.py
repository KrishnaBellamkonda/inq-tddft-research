"""Draft-2 bath density difference figure (fig 7 base).

n_bath_classical - n_bath_WP in the xz mid-plane at step 100 (t = 2.0 a.u.).
Proxy: sigma_WP=3, r_s=5.69, v=2.711.

Outputs (600 DPI, bbox_inches=None):
  bath_density_diff_linear.png   — linear diverging norm
  bath_density_diff_symlog.png   — symlog diverging norm

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 build_bath_diff.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator, ScalarFormatter

REPO = Path(__file__).resolve().parents[7]
sys.path.insert(0, str(REPO / "inq-stack/python"))

from inqview.visualisation import style
from inqview import load_vti

HERE = Path(__file__).parent
RUNS = REPO / "ResearchProject/systems/jellium"

WP_VTI = RUNS / "run_wp_n162_L50_E100_sigma3_wf" / "results" / "raw" / "vti"
CL_VTI = RUNS / "run_classical_n162_L50_E100"    / "results" / "raw" / "vti"

WP_CZ = -13.0
WP_K0 = 2.71106
STEP  = 100
DT_STEP = 0.02
T_PHYS = STEP * DT_STEP
Z_WP = WP_CZ + WP_K0 * T_PHYS   # ≈ -7.58 Bohr

CMAP = "RdBu_r"


def _fmt_cb(cb) -> None:
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    cb.ax.yaxis.set_major_formatter(fmt)
    cb.ax.yaxis.set_major_locator(MaxNLocator(5))


def _make_panel(delta_zx, xs, zs, norm, *, suffix: str) -> Path:
    style.apply_theme()
    fig = plt.figure(figsize=(3.5, 3.5))
    ax  = fig.add_axes([0.14, 0.13, 0.82, 0.82])

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
    ax.legend(loc="upper left", handlelength=1.5)
    ax.set_xlabel(r"$z\ (\mathrm{Bohr})$")
    ax.set_ylabel(r"$x\ (\mathrm{Bohr})$")

    out = HERE / f"bath_density_diff_{suffix}.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {out}")
    return out


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
    delta_zx   = delta[:, iy_mid, :]   # (nx, nz)

    vmax = float(np.max(np.abs(delta_zx)))
    print(f"Proxy: sigma_WP=3, r_s=5.69, v=2.71, t={T_PHYS:.1f} a.u. (step {STEP})")
    print(f"|delta|_max = {vmax:.3e} Bohr^-3")

    norm_lin = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    _make_panel(delta_zx, xs, zs, norm_lin, suffix="linear")

    linthresh = vmax / 100.0
    norm_log  = mcolors.SymLogNorm(linthresh=linthresh, linscale=0.5,
                                   vmin=-vmax, vmax=vmax)
    _make_panel(delta_zx, xs, zs, norm_log, suffix="symlog")


if __name__ == "__main__":
    main()
