"""fig_momentum_before_after — WP momentum distribution before/after collision.

Two-panel:
(a) |ψ̃_WP(|k|)|² vs |k| on linear scale — initial Gaussian and post-collision
(b) Placeholder for 2D Δ|ψ̃(k_z, k_⊥)|² difference map (awaiting orbital output)

Data: E=100 eV, σ=1 Bohr, N=162, L=50 (high-density jellium, r_s ≈ 5.69).

Run:
    python -m applications.report1.fig_momentum_before_after
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from applications.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
    palette_sweep5,
)

MOM_CSV = (
    "ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1/"
    "results/raw/observables/momentum_distribution.csv"
)
OUT = "docs/reports/report1/figures/fig_momentum_before_after.png"

HA_TO_EV = 27.21138625
AU_TO_FS = 1.0 / 41.341374575751
E_KIN_EV = 100.0
K0 = np.sqrt(2.0 * E_KIN_EV / HA_TO_EV)


def main() -> None:
    apply_style()

    df = pd.read_csv(MOM_CSV, comment="#")
    steps = sorted(df["step"].unique())
    step_0 = steps[0]
    step_final = steps[-1]

    t_0_fs = df[df["step"] == step_0]["time_au"].iloc[0] * AU_TO_FS
    t_f_fs = df[df["step"] == step_final]["time_au"].iloc[0] * AU_TO_FS

    df_0 = df[df["step"] == step_0].sort_values("k_bohr_inv")
    df_f = df[df["step"] == step_final].sort_values("k_bohr_inv")

    k = df_0["k_bohr_inv"].values
    n_wp_0 = df_0["n_wp"].values
    n_wp_f = df_f["n_wp"].values

    c_before = palette_sweep5[4]
    c_after = palette_sweep5[0]

    W = column_widths_in["single"]
    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(W, W * 1.4),
        gridspec_kw={"height_ratios": [2, 1.5], "hspace": 0.30},
    )

    # (a) WP momentum distribution — linear scale
    ax_a.plot(k, n_wp_0, color=c_before, linewidth=1.2,
              label=rf"$t = {t_0_fs:.0f}$ fs")
    ax_a.plot(k, n_wp_f, color=c_after, linewidth=1.2,
              label=rf"$t = {t_f_fs:.3f}$ fs")
    ax_a.fill_between(k, n_wp_0, alpha=0.08, color=c_before)

    ax_a.axvline(K0, color="#808080", linewidth=0.7, linestyle="--", alpha=0.7)
    ax_a.text(K0 + 0.1, n_wp_0.max() * 0.85,
              r"$k_0$", fontsize=8, color="#606060", ha="left", va="center")

    ax_a.set_xlabel(r"$|k|$ (Bohr$^{-1}$)")
    ax_a.set_ylabel(r"$|\tilde\psi_{\mathrm{WP}}(k)|^2$")
    ax_a.set_xlim(0, k.max())
    ax_a.set_ylim(bottom=0)
    ax_a.legend(fontsize=7, loc="upper right", frameon=False)
    panel_label(ax_a, "(a)")

    # (b) Placeholder for 2D Δ|ψ̃(k_z, k_⊥)|² difference map
    ax_b.set_xlim(-5, 5)
    ax_b.set_ylim(-3, 3)
    ax_b.set_facecolor("#F5F5F0")
    ax_b.set_xlabel(r"$k_z$ (Bohr$^{-1}$)")
    ax_b.set_ylabel(r"$k_\perp$ (Bohr$^{-1}$)")

    ax_b.text(
        0.5, 0.5,
        (r"$\Delta|\tilde{\psi}_{\mathrm{WP}}(k_z, k_\perp)|^2$"
         "\n\nAwaiting 3D orbital output\n"
         "from updated simulation"),
        transform=ax_b.transAxes, fontsize=8, ha="center", va="center",
        color="#808080", fontstyle="italic",
    )

    # Sketch the expected features as annotations
    ax_b.annotate("", xy=(-K0, 0), xytext=(K0, 0),
                  arrowprops=dict(arrowstyle="<->", color="#C0C0C0", lw=0.8))
    ax_b.text(0, -0.4, r"$2k_0$", fontsize=6, ha="center", color="#C0C0C0")

    panel_label(ax_b, "(b)")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
