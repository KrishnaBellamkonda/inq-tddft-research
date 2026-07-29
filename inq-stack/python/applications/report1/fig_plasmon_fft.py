"""fig_plasmon_fft — Plasmon detection in jellium (two-panel).

(a) Fourier power spectrum of density-FT modes n_q(omega) for m=1,2,3
    with Bohm-Gross prediction.
(b) Real-space Delta n(x,z) density-difference at mid-propagation,
    integrated over y.

Data: run_plasmon_n162_L50_E15 (dedicated long-time plasmon run).

Run:
    python -m applications.report1.fig_plasmon_fft
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.colors import SymLogNorm

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
)

CSV = (
    "ResearchProject/systems/jellium/run_plasmon_n162_L50_E15/"
    "results/analysis/observables/n_q_spectrum.csv"
)
VTI_DIR = (
    "ResearchProject/systems/jellium/run_plasmon_n162_L50_E15/"
    "results/raw/vti/density_rt_delta"
)

HA_TO_EV = 27.2114
L_BOHR = 50.0
N_EL = 1.296e-3  # density a.u.
OMEGA_P_HA = np.sqrt(4 * np.pi * N_EL)
OMEGA_P_EV = OMEGA_P_HA * HA_TO_EV
V_F = 0.337  # Fermi velocity a.u.


def bohm_gross(m: int) -> float:
    """Bohm-Gross dispersion: omega^2 = omega_p^2 + (3/5)v_F^2 q^2."""
    q = 2 * np.pi * m / L_BOHR
    return np.sqrt(OMEGA_P_HA**2 + 3/5 * V_F**2 * q**2) * HA_TO_EV


def load_vti_density_delta(path: str) -> np.ndarray:
    """Load a single VTI file and return 3D numpy array."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(path)
    reader.Update()
    img = reader.GetOutput()
    dims = img.GetDimensions()
    arr = vtk_to_numpy(img.GetPointData().GetArray(0))
    # VTK uses Fortran ordering: x varies fastest
    return arr.reshape(dims[2], dims[1], dims[0])


def main() -> None:
    apply_style()

    df = pd.read_csv(CSV)

    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(W, W * 0.38),
                                      gridspec_kw={"width_ratios": [1.2, 1],
                                                   "wspace": 0.35})

    # ── Panel (a): FFT spectrum, m=1,2,3 only ──
    colors = [palette_sweep5[0], palette_sweep5[2], palette_sweep5[4]]
    omega_max = 12.0

    modes_to_plot = [1, 2, 3]
    for i, m in enumerate(modes_to_plot):
        sub = df[df["m"] == m]
        mask = sub["omega_eV"] <= omega_max
        omega = sub.loc[mask, "omega_eV"].values
        amp = sub.loc[mask, "abs_FFT_n_q"].values
        ax_a.semilogy(omega, amp, linewidth=0.8, color=colors[i],
                      label=f"$m={m}$", alpha=0.85)

    # Bohm-Gross predictions for m=1,2,3
    for i, m in enumerate(modes_to_plot):
        bg = bohm_gross(m)
        ax_a.axvline(bg, color=colors[i], ls=":", lw=0.5, alpha=0.5, zorder=1)

    # Mark Bohm-Gross label once
    bg_m1 = bohm_gross(1)
    ax_a.text(bg_m1 + 0.2, ax_a.get_ylim()[1] * 0.15 if ax_a.get_ylim()[1] > 1 else 1e1,
              f"BG $m{{=}}1$: {bg_m1:.2f} eV",
              fontsize=5, color="grey", va="top")

    ax_a.set_xlabel(r"$\omega$ (eV)")
    ax_a.set_ylabel(r"$|\mathrm{FFT}[n_{q_m}(t)]|$")
    ax_a.set_xlim(0, omega_max)
    ax_a.legend(fontsize=6, loc="upper right", frameon=False)
    panel_label(ax_a, "(a)")

    # ── Panel (b): Real-space Delta n(x,z) at mid-propagation ──
    import os
    vti_path = os.path.join(VTI_DIR, "density_delta_t050000.vti")
    dn_xz = None
    if os.path.exists(vti_path):
        dn_3d = load_vti_density_delta(vti_path)
        dn_xz = dn_3d.sum(axis=1)  # sum over y → (nz, nx)

    _finish_panel_b(fig, ax_a, ax_b, dn_xz, use_log=False,
                    out="docs/reports/report1/figures/fig_plasmon_fft.png")

    # ── Log-scale version ──
    fig2, (ax_a2, ax_b2) = plt.subplots(1, 2, figsize=(W, W * 0.38),
                                         gridspec_kw={"width_ratios": [1.2, 1],
                                                      "wspace": 0.35})
    # Re-draw panel (a) identically
    for i, m in enumerate(modes_to_plot):
        sub = df[df["m"] == m]
        mask = sub["omega_eV"] <= omega_max
        omega = sub.loc[mask, "omega_eV"].values
        amp = sub.loc[mask, "abs_FFT_n_q"].values
        ax_a2.semilogy(omega, amp, linewidth=0.8, color=colors[i],
                       label=f"$m={m}$", alpha=0.85)
    for i, m in enumerate(modes_to_plot):
        bg = bohm_gross(m)
        ax_a2.axvline(bg, color=colors[i], ls=":", lw=0.5, alpha=0.5, zorder=1)
    ax_a2.text(bg_m1 + 0.2, ax_a2.get_ylim()[1] * 0.15 if ax_a2.get_ylim()[1] > 1 else 1e1,
               f"BG $m{{=}}1$: {bg_m1:.2f} eV", fontsize=5, color="grey", va="top")
    ax_a2.set_xlabel(r"$\omega$ (eV)")
    ax_a2.set_ylabel(r"$|\mathrm{FFT}[n_{q_m}(t)]|$")
    ax_a2.set_xlim(0, omega_max)
    ax_a2.legend(fontsize=6, loc="upper right", frameon=False)
    panel_label(ax_a2, "(a)")

    _finish_panel_b(fig2, ax_a2, ax_b2, dn_xz, use_log=True,
                    out="docs/reports/report1/figures/fig_plasmon_fft_log.png")


def _finish_panel_b(fig, ax_a, ax_b, dn_xz, *, use_log: bool, out: str):
    if dn_xz is not None:
        nx, nz = dn_xz.shape[1], dn_xz.shape[0]
        x_edges = np.linspace(-L_BOHR/2, L_BOHR/2, nx + 1)
        z_edges = np.linspace(-L_BOHR/2, L_BOHR/2, nz + 1)

        vmax = np.max(np.abs(dn_xz)) * 0.8
        if use_log:
            linthresh = 1e-4 * vmax
            norm = SymLogNorm(linthresh=linthresh, linscale=1.0,
                              vmin=-vmax, vmax=vmax)
        else:
            norm = plt.Normalize(vmin=-vmax, vmax=vmax)

        im = ax_b.pcolormesh(z_edges, x_edges, dn_xz.T,
                             cmap="RdBu_r", norm=norm,
                             rasterized=True)
        cbar = fig.colorbar(im, ax=ax_b, fraction=0.046, pad=0.04)
        cbar.set_label(r"$\sum_y \Delta n$ (a.u.)", fontsize=7)
        cbar.ax.tick_params(labelsize=6)
    else:
        ax_b.text(0.5, 0.5, "VTI data not found",
                  transform=ax_b.transAxes, ha="center", va="center",
                  fontsize=8, color="#808080")

    ax_b.set_xlabel(r"$z$ (Bohr)")
    ax_b.set_ylabel(r"$x$ (Bohr)")
    ax_b.set_aspect("equal")
    panel_label(ax_b, "(b)")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic ({out}): {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
