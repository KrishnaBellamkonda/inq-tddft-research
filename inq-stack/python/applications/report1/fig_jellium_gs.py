"""fig_jellium_gs — Jellium ground-state density slices.

Two-panel: (a) N=162 closed shell (uniform), (b) N=138 partial shell
(broken-symmetry artefact from half-filled |G|²=6 shell).

Run:
    python -m applications.report1.fig_jellium_gs
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
)

VTI_162 = (
    "ResearchProject/systems/jellium/save_gs/"
    "gs_L50_cubic_N162_dx0p40/results/density_gs_system/density_gs_system.vti"
)
VTI_138 = (
    "ResearchProject/systems/jellium/save_gs/"
    "gs_L60_cubic_N138_dx0p55/results/density_gs_system/density_gs_system.vti"
)


def load_density_slice(vti_path: str) -> tuple[np.ndarray, np.ndarray, float]:
    """Load VTI, return z=0 plane slice, coordinates, and mean density."""
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(vti_path)
    reader.Update()
    data = reader.GetOutput()
    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    origin = data.GetOrigin()
    arr = data.GetPointData().GetArray(0)
    flat = np.array([arr.GetValue(i) for i in range(arr.GetNumberOfTuples())])
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    iz_mid = dims[2] // 2
    slice_xy = rho[:, :, iz_mid]
    x = origin[0] + np.arange(dims[0]) * spacing[0]
    return slice_xy, x, rho.mean()


def main() -> None:
    apply_style()

    print("Loading N=162 GS density...")
    slice_162, x_162, mean_162 = load_density_slice(VTI_162)
    print("Loading N=138 GS density...")
    slice_138, x_138, mean_138 = load_density_slice(VTI_138)

    norm_162 = slice_162 / mean_162
    norm_138 = slice_138 / mean_138

    mod_162 = (norm_162.max() - norm_162.min()) / 2 * 100
    mod_138 = (norm_138.max() - norm_138.min()) / 2 * 100
    print(f"N=162 modulation: {mod_162:.1f}%, N=138 modulation: {mod_138:.1f}%")

    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(W, W * 0.42))

    # Use tight colorbar centered on 1.0 to highlight contrast
    vmin, vmax = 0.97, 1.03

    im_a = ax_a.imshow(norm_162.T, origin="lower", cmap="RdBu_r",
                       vmin=vmin, vmax=vmax, aspect="equal",
                       extent=[x_162[0], x_162[-1], x_162[0], x_162[-1]])
    ax_a.set_xlabel(r"$x$ (Bohr)")
    ax_a.set_ylabel(r"$y$ (Bohr)")
    panel_label(ax_a, "(a)")
    ax_a.text(0.5, 0.97, f"$N=162$, mod.~${mod_162:.1f}$\\%",
              transform=ax_a.transAxes, fontsize=6, ha="center", va="top",
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

    im_b = ax_b.imshow(norm_138.T, origin="lower", cmap="RdBu_r",
                       vmin=vmin, vmax=vmax, aspect="equal",
                       extent=[x_138[0], x_138[-1], x_138[0], x_138[-1]])
    ax_b.set_xlabel(r"$x$ (Bohr)")
    ax_b.set_yticklabels([])
    panel_label(ax_b, "(b)")
    ax_b.text(0.5, 0.97, f"$N=138$, mod.~${mod_138:.1f}$\\%",
              transform=ax_b.transAxes, fontsize=6, ha="center", va="top",
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

    cbar = fig.colorbar(im_b, ax=[ax_a, ax_b], shrink=0.85, pad=0.02)
    cbar.set_label(r"$n / \bar{n}$", fontsize=8)
    cbar.ax.tick_params(labelsize=6)

    out = "docs/reports/report1/figures/fig_jellium_gs.png"
    with TufteCritic.disabled():
        fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
