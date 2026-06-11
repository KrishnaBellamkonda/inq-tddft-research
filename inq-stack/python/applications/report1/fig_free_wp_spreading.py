"""fig_free_wp_spreading — Free wave-packet spreading verification.

Two-panel: (a) σ_r(t) from analytical and INQ;
(b) residual from analytical, interference-free zone only.

Data: run_free_wp_L50_E100 (σ₀=5 Bohr, E=100 eV).

Run:
    python -m applications.report1.fig_free_wp_spreading
"""

from __future__ import annotations

import numpy as np
import pandas as pd
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

BASE = "ResearchProject/systems/jellium/run_free_wp_L50_E100/results/raw/observables"
INQ_CSV = f"{BASE}/wp_real_space_stats.csv"

SIGMA0 = 5.0
K0 = 2.711063340102429
Z0 = -5.0            # launch position (Bohr)
L_BOHR = 50.0
T_MAX = 7.0           # show full data including interference onset


def sigma_analytic(t: np.ndarray, sigma0: float = SIGMA0) -> np.ndarray:
    """Correct σ_density(t) = √(σ₀²/2 + t²/(2σ₀²))."""
    return np.sqrt(sigma0**2 / 2 + t**2 / (2 * sigma0**2))


def compute_t_if() -> float:
    """Time when centroid + 4σ_r reaches the box boundary."""
    from scipy.optimize import brentq
    z_boundary = L_BOHR / 2.0
    def leading_edge(t):
        return Z0 + K0 * t + 4.0 * sigma_analytic(np.array([t]))[0] - z_boundary
    return brentq(leading_edge, 0, 20)


def main() -> None:
    apply_style()

    T_IF = compute_t_if()

    # INQ data
    inq = pd.read_csv(INQ_CSV, comment="#")
    mask_inq = inq["time_au"] <= T_MAX
    t_inq = inq.loc[mask_inq, "time_au"].values
    sigma_inq = np.sqrt(inq.loc[mask_inq, "sigma_z2"].values)

    # Analytical
    sigma_ana_inq = sigma_analytic(t_inq)

    # Residuals (Bohr)
    resid = sigma_inq - sigma_ana_inq

    c_ana = "black"
    c_inq = palette_sweep5[0]

    W = column_widths_in["single"]
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(W, W * 1.3),
                                      gridspec_kw={"height_ratios": [2, 1],
                                                   "hspace": 0.08})

    # (a) σ(t) with interference-free zone shading
    t_dense = np.linspace(0, T_MAX, 500)
    ax_a.axvspan(0, T_IF, alpha=0.10, color="#188048", zorder=0)
    ax_a.axvline(T_IF, color="#188048", linewidth=0.6, linestyle=":",
                 zorder=1, alpha=0.6)
    ax_a.plot(t_dense, sigma_analytic(t_dense), "-", color=c_ana,
              linewidth=1.0, label="Analytical", zorder=3)
    ax_a.plot(t_inq[::3], sigma_inq[::3], "-o", color=c_inq,
              markersize=2.5, linewidth=0.8, label="INQ", zorder=4)

    ax_a.set_ylabel(r"$\sigma_r(t)$ (Bohr)")
    ax_a.set_xticklabels([])
    ax_a.legend(fontsize=6, loc="center left", frameon=False)
    panel_label(ax_a, "(a)")

    yhi = ax_a.get_ylim()[1]
    ax_a.text(T_IF * 0.45, yhi - 0.002,
              "interference-free domain", fontsize=5.5, ha="center", va="top",
              color="#188048", fontstyle="italic")

    # (b) Residual (Bohr) — full range with IFW shading
    ax_b.axvspan(0, T_IF, alpha=0.10, color="#188048", zorder=0)
    ax_b.axvline(T_IF, color="#188048", linewidth=0.6, linestyle=":",
                 zorder=1, alpha=0.6)
    ax_b.plot(t_inq, resid, "-o", color=c_inq, linewidth=0.9,
              markersize=1.8)
    ax_b.axhline(0, color="#b0b0b0", linewidth=0.4, zorder=0)

    mask_if = t_inq <= T_IF
    max_resid_if = np.max(np.abs(resid[mask_if]))
    ax_b.text(0.97, 0.95,
              f"IFW max: {max_resid_if:.5f} Bohr",
              transform=ax_b.transAxes, fontsize=5.5, va="top", ha="right",
              bbox=dict(facecolor="white", edgecolor="#b0b0b0",
                        linewidth=0.3, pad=2, alpha=0.9))

    ax_b.set_xlabel(r"$t$ (a.u.)")
    ax_b.set_ylabel(r"Residual (Bohr)")
    panel_label(ax_b, "(b)")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig_free_wp_spreading.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
