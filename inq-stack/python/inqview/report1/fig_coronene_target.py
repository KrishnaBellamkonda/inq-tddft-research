"""fig_coronene_target — Coronene C24H12 ground-state xy density with atom overlay.

Shows the z-integrated GS electron density n_2D(x,y) as a heatmap with:
  - Carbon atom positions (larger dots)
  - Hydrogen atom positions (smaller dots)
  - WP impact positions (center + CC-bond midpoint) as cross/plus marks
  - 1-sigma Gaussian width circle for the WP

Run:
    python -m inqview.report1.fig_coronene_target
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from pathlib import Path

from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
    palette_sweep5,
)

# --- Paths -------------------------------------------------------------------
GS_VTI = Path(
    "ResearchProject/systems/coronene/run_save_gs_paper_replica/"
    "results/density_gs/density_t000000.vti"
)
XYZ_FILE = Path(
    "ResearchProject/systems/coronene/configurations/"
    "tsubonoya_2014_paper_replica/coronene_centred.xyz"
)
OUT_LOG = Path("docs/reports/report1/figures/fig_coronene_target.png")
OUT_LIN = Path("docs/reports/report1/figures/fig_coronene_target_linear.png")

# --- Config from tsubonoya_2014_coronene.hpp ---------------------------------
ANG_TO_BOHR = 1.8897259886
WP_SIGMA_BOHR = 0.53 * ANG_TO_BOHR   # ~1.0015 Bohr
WP_CX_CENTER = 0.0                     # center impact
WP_CY_CENTER = 0.0
WP_CX_CCBOND = 2.1315 * ANG_TO_BOHR   # ~4.028 Bohr — CC bond midpoint


def load_vti_density(path: str) -> tuple[np.ndarray, tuple, tuple]:
    """Load a VTI density file. Returns (3D array [ix,iy,iz], origin, spacing)."""
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()

    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    origin = data.GetOrigin()

    arr = data.GetPointData().GetArray(0)
    n_pts = arr.GetNumberOfTuples()
    flat = np.array([arr.GetValue(i) for i in range(n_pts)])

    # VTK uses Fortran ordering (x fastest)
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    return rho, origin, spacing


def parse_xyz(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Parse XYZ file. Returns (positions_ang, positions_bohr, species)."""
    lines = path.read_text().strip().split("\n")
    n_atoms = int(lines[0])
    species = []
    positions_ang = []
    for line in lines[2: 2 + n_atoms]:
        parts = line.split()
        species.append(parts[0])
        positions_ang.append([float(parts[1]), float(parts[2]), float(parts[3])])
    pos_ang = np.array(positions_ang)
    pos_bohr = pos_ang * ANG_TO_BOHR
    return pos_ang, pos_bohr, species


def make_figure(n2d_crop, x_crop, y_crop, c_pos, h_pos, *, use_log: bool, out: Path):
    W = column_widths_in["full"]
    fig, ax = plt.subplots(figsize=(W, W * 0.85))

    vmax = n2d_crop.max()
    if use_log:
        vmin = max(1e-3 * vmax, n2d_crop[n2d_crop > 0].min())
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None

    im = ax.pcolormesh(
        y_crop, x_crop, n2d_crop,
        cmap="inferno",
        norm=norm,
        shading="auto",
        rasterized=True,
    )
    ax.set_aspect("equal")

    ax.scatter(c_pos[:, 1], c_pos[:, 0], s=28, c="#00CC66",
               edgecolors="white", linewidths=0.4, zorder=5, label=r"C")
    ax.scatter(h_pos[:, 1], h_pos[:, 0], s=12, c="#66CCFF",
               edgecolors="white", linewidths=0.3, zorder=5, label=r"H")

    ax.plot(WP_CY_CENTER, WP_CX_CENTER, "+", color="white",
            markersize=10, markeredgewidth=1.5, zorder=6, label=r"Target centre")
    ax.plot(0.0, WP_CX_CCBOND, "x", color="#FF6666",
            markersize=8, markeredgewidth=1.5, zorder=6, label=r"Target C--C bond")

    sigma_circle = Circle((WP_CY_CENTER, WP_CX_CENTER), WP_SIGMA_BOHR,
                           fill=False, edgecolor="white", linewidth=0.8,
                           linestyle="--", zorder=6)
    ax.add_patch(sigma_circle)
    ax.annotate(rf"$\sigma = {WP_SIGMA_BOHR:.2f}$ Bohr",
                xy=(WP_CY_CENTER + WP_SIGMA_BOHR * 0.71,
                    WP_CX_CENTER + WP_SIGMA_BOHR * 0.71),
                fontsize=7, color="white", ha="left", va="bottom")

    sigma_circle_cc = Circle((0.0, WP_CX_CCBOND), WP_SIGMA_BOHR,
                              fill=False, edgecolor="#FF6666", linewidth=0.8,
                              linestyle="--", zorder=6)
    ax.add_patch(sigma_circle_cc)

    ax.set_xlabel(r"$y$ (Bohr)")
    ax.set_ylabel(r"$x$ (Bohr)")

    cbar = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
    cbar.set_label(r"$n_{2\mathrm{D}}(x,y)$ (e/Bohr$^2$)")

    leg = ax.legend(loc="upper left", fontsize=7, framealpha=0.85,
                    edgecolor="#404040", handletextpad=0.3,
                    borderpad=0.3, labelspacing=0.25)
    leg.get_frame().set_facecolor("black")
    for text in leg.get_texts():
        text.set_color("white")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic ({out.name}): {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(str(out), dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {out}")
    plt.close(fig)


def main() -> None:
    apply_style()

    print("Loading GS density VTI...")
    rho, origin, spacing = load_vti_density(str(GS_VTI))
    dx, dy, dz = spacing
    nx, ny, nz = rho.shape
    n2d = rho.sum(axis=2) * dz

    x = origin[0] + np.arange(nx) * dx
    y = origin[1] + np.arange(ny) * dy

    pos_ang, pos_bohr, species = parse_xyz(XYZ_FILE)
    c_mask = np.array([s == "C" for s in species])
    h_mask = np.array([s == "H" for s in species])
    c_pos = pos_bohr[c_mask]
    h_pos = pos_bohr[h_mask]

    margin = 4.0
    x_max_atom = max(np.abs(pos_bohr[:, 0]).max(), np.abs(pos_bohr[:, 1]).max())
    crop = x_max_atom + margin

    ix_lo = np.searchsorted(x, -crop)
    ix_hi = np.searchsorted(x, +crop)
    iy_lo = np.searchsorted(y, -crop)
    iy_hi = np.searchsorted(y, +crop)

    n2d_crop = n2d[ix_lo:ix_hi, iy_lo:iy_hi]
    x_crop = x[ix_lo:ix_hi]
    y_crop = y[iy_lo:iy_hi]

    make_figure(n2d_crop, x_crop, y_crop, c_pos, h_pos,
                use_log=True, out=OUT_LOG)
    make_figure(n2d_crop, x_crop, y_crop, c_pos, h_pos,
                use_log=False, out=OUT_LIN)


if __name__ == "__main__":
    main()
