"""fig_momentum_before_after_highdens — Momentum distribution for high-density jellium.

Two panels: (a) WP orbital |ψ(k)|², (b) total system n(k), both before/after.

Run:
    python -m inqview.report1.fig_momentum_before_after_highdens
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.report1._shared_style import (
    apply_style, palette_sweep5, column_widths_in, panel_label, TufteCritic,
)

RUN = "ResearchProject/systems/jellium/run_wp_n162_L30_E100_highdens_sigma1"
MOM_CSV = f"{RUN}/results/raw/observables/momentum_distribution.csv"
OUT = "docs/reports/report1/figures/fig_momentum_before_after_highdens.png"

HA_TO_EV = 27.2114
AU_TO_FS = 0.02418884
K0 = np.sqrt(2.0 * 100.0 / HA_TO_EV)


def main() -> None:
    apply_style()

    df = pd.read_csv(MOM_CSV, comment="#")
    steps = sorted(df["step"].unique())
    step_0, step_f = steps[0], steps[-1]
    t_0_fs = df[df["step"] == step_0]["time_au"].iloc[0] * AU_TO_FS
    t_f_fs = df[df["step"] == step_f]["time_au"].iloc[0] * AU_TO_FS

    df_0 = df[df["step"] == step_0].sort_values("k_bohr_inv")
    df_f = df[df["step"] == step_f].sort_values("k_bohr_inv")

    k = df_0["k_bohr_inv"].values
    c_before = palette_sweep5[4]
    c_after = palette_sweep5[0]

    W = column_widths_in["single"]
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(W, W * 1.1), sharex=True,
                                      gridspec_kw={"hspace": 0.08})

    # (a) WP orbital
    ax_a.semilogy(k, df_0["n_wp"].values, color=c_before, linewidth=1.0,
                  label=rf"$t = {t_0_fs:.0f}$ fs")
    ax_a.semilogy(k, df_f["n_wp"].values, color=c_after, linewidth=1.0,
                  label=rf"$t = {t_f_fs:.3f}$ fs")
    ax_a.axvline(K0, color="#808080", linewidth=0.7, linestyle="--", alpha=0.7)
    ax_a.text(K0 + 0.1, 1e-4, r"$k_0$", fontsize=7, color="#606060")
    ax_a.set_ylabel(r"$|\tilde\psi_{\mathrm{WP}}(k)|^2$")
    ax_a.legend(fontsize=6, loc="upper right", frameon=False)
    wp_max = max(df_0["n_wp"].max(), df_f["n_wp"].max())
    ax_a.set_ylim(1e-15, wp_max * 5)
    panel_label(ax_a, "(a)")

    # (b) Total
    ax_b.semilogy(k, df_0["n_total"].values, color=c_before, linewidth=1.0,
                  label=rf"$t = {t_0_fs:.0f}$ fs")
    ax_b.semilogy(k, df_f["n_total"].values, color=c_after, linewidth=1.0,
                  label=rf"$t = {t_f_fs:.3f}$ fs")
    ax_b.axvline(K0, color="#808080", linewidth=0.7, linestyle="--", alpha=0.7)
    ax_b.set_xlabel(r"$|k|$ (Bohr$^{-1}$)")
    ax_b.set_ylabel(r"$n_{\mathrm{total}}(k)$")
    ax_b.legend(fontsize=6, loc="upper right", frameon=False)
    panel_label(ax_b, "(b)")

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
