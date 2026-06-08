"""fig_free_wp_momentum — Momentum-space verification for free WP.

Two-panel: (a) Analytical |ψ̃(p_z)|² Gaussian with INQ ⟨p_z⟩ ± σ_pz
at selected times; (b) ⟨p_z⟩(t) and σ_pz(t) vs analytical constants.

For free propagation the momentum distribution is time-independent:
the WP spreads in real space but its momentum content is frozen.

Data: run_free_wp_L50_E100 (σ₀=5 Bohr, E=100 eV).

Run:
    python -m inqview.report1.fig_free_wp_momentum
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
)

BASE = "ResearchProject/systems/jellium/run_free_wp_L50_E100/results/raw/observables"
MOM_CSV = f"{BASE}/wp_momentum_stats.csv"

SIGMA0 = 5.0
K0 = 2.711063340102429
T_IF = 5.0

SIGMA_DENSITY_0 = SIGMA0 / np.sqrt(2)
SIGMA_P_ANA = 1.0 / (2.0 * SIGMA_DENSITY_0)


def momentum_prob(p: np.ndarray, p0: float, sigma_p: float) -> np.ndarray:
    return (1.0 / (sigma_p * np.sqrt(2 * np.pi))) * np.exp(
        -0.5 * ((p - p0) / sigma_p) ** 2
    )


def main() -> None:
    apply_style()

    mom = pd.read_csv(MOM_CSV, comment="#")
    mask = mom["time_au"] <= T_IF
    t = mom.loc[mask, "time_au"].values
    pz_mean = mom.loc[mask, "pz_mean"].values
    sigma_pz = np.sqrt(mom.loc[mask, "sigma_pz2"].values)

    c_ana = "black"
    c_inq = palette_sweep5[0]
    c_sigma = palette_sweep5[3]

    W = column_widths_in["single"]
    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(W, W * 1.3),
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.25},
    )

    # (a) Analytical momentum probability with INQ markers at selected times
    p_range = np.linspace(K0 - 6 * SIGMA_P_ANA, K0 + 6 * SIGMA_P_ANA, 500)
    prob_ana = momentum_prob(p_range, K0, SIGMA_P_ANA)
    ax_a.plot(p_range, prob_ana, "-", color=c_ana, linewidth=1.0,
              label="Analytical")
    ax_a.fill_between(p_range, prob_ana, alpha=0.06, color=c_ana)

    sample_times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    markers = ["o", "s", "^", "v", "D", "p"]
    for i, t_sample in enumerate(sample_times):
        idx = np.argmin(np.abs(t - t_sample))
        p_mean_i = pz_mean[idx]
        sig_i = sigma_pz[idx]
        prob_at_mean = momentum_prob(np.array([p_mean_i]), K0, SIGMA_P_ANA)[0]
        label = f"$t = {t[idx]:.1f}$ a.u." if i in (0, len(sample_times) - 1) else None
        ax_a.plot(p_mean_i, prob_at_mean, markers[i % len(markers)],
                  color=c_inq, markersize=4, zorder=5, label=label)

    ax_a.set_xlabel(r"$p_z$ (a.u.)")
    ax_a.set_ylabel(r"$|\tilde{\psi}(p_z)|^2$")
    ax_a.legend(fontsize=6, loc="upper right", frameon=False)
    panel_label(ax_a, "(a)")

    # (b) Relative deviation from analytical constants
    rel_pz = (pz_mean - K0) / K0
    rel_sig = (sigma_pz - SIGMA_P_ANA) / SIGMA_P_ANA

    ax_b.axvspan(0, T_IF, alpha=0.10, color="#188048", zorder=0)
    ax_b.axhline(0, color="#b0b0b0", linewidth=0.4, zorder=0)
    ax_b.plot(t, rel_pz * 100, "-o", color=c_inq, markersize=1.5,
              linewidth=0.8, label=r"$\langle p_z \rangle / k_0 - 1$")
    ax_b.plot(t, rel_sig * 100, "-o", color=c_sigma, markersize=1.5,
              linewidth=0.8, label=r"$\sigma_{p_z} / \sigma_{p,\mathrm{ana}} - 1$")

    ax_b.set_xlabel(r"$t$ (a.u.)")
    ax_b.set_ylabel(r"Relative deviation (\%)")
    ax_b.legend(fontsize=6, loc="lower left", frameon=False)
    panel_label(ax_b, "(b)")

    dpz = np.max(np.abs(rel_pz)) * 100
    dsig = np.max(np.abs(rel_sig)) * 100
    ax_b.text(
        0.97, 0.95,
        (f"$\\langle p_z\\rangle$ max dev: {dpz:.2e}\\%\n"
         f"$\\sigma_{{p_z}}$ max dev: {dsig:.2e}\\%"),
        transform=ax_b.transAxes, fontsize=5.5, va="top", ha="right",
        bbox=dict(facecolor="white", edgecolor="#b0b0b0",
                  linewidth=0.3, pad=2, alpha=0.9),
    )

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig_free_wp_momentum.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
