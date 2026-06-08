"""fig_energy_decomposition_highdens — Energy decomposition for high-density jellium.

Two-panel:
  (a) Total electronic system: ΔE_total, ΔE_kinetic, ΔE_hartree, ΔE_xc
  (b) Wave-packet orbital only: ⟨p⟩²/2m, σ_p²/2m, their sum

Run:
    python -m inqview.report1.fig_energy_decomposition_highdens
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
OBS_CSV = f"{RUN}/results/raw/observables/observables.csv"
MOM_CSV = f"{RUN}/results/raw/observables/wp_momentum_stats.csv"
OUT = "docs/reports/report1/figures/fig_energy_decomposition_highdens.png"

HA_TO_EV = 27.2114
AU_TO_FS = 0.02418884


def main() -> None:
    apply_style()

    obs = pd.read_csv(OBS_CSV)
    mom = pd.read_csv(MOM_CSV, comment="#")

    t_obs_fs = obs["time_au"].values * AU_TO_FS
    t_mom_fs = mom["time_au"].values * AU_TO_FS

    W = column_widths_in["single"]
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(W, W * 1.1), sharex=True,
                                      gridspec_kw={"hspace": 0.08})

    # Panel (a): Total electronic system
    for col, label, color in [
        ("energy_total", r"$\Delta E_{\mathrm{total}}$", "black"),
        ("energy_kinetic", r"$\Delta E_{\mathrm{kin}}$", palette_sweep5[0]),
        ("energy_hartree", r"$\Delta E_H$", palette_sweep5[2]),
        ("energy_ion_xc", r"$\Delta E_{\mathrm{xc}}$", palette_sweep5[4]),
    ]:
        if col in obs.columns:
            dE = (obs[col].values - obs[col].iloc[0]) * HA_TO_EV
            ax_a.plot(t_obs_fs, dE, "-", color=color, linewidth=0.8, label=label)

    ax_a.axhline(0, color="#b0b0b0", linewidth=0.3)
    ax_a.set_ylabel(r"$\Delta E$ (eV)")
    ax_a.legend(fontsize=5.5, loc="best", frameon=False)
    ax_a.set_title("Total electronic system", fontsize=8, pad=3)
    panel_label(ax_a, "(a)")

    # Panel (b): WP orbital only
    pz = mom["pz_mean"].values
    sig_pz2 = mom["sigma_pz2"].values
    E_centroid = pz**2 / 2 * HA_TO_EV
    E_spread = sig_pz2 / 2 * HA_TO_EV
    E_total_wp = E_centroid + E_spread

    ax_b.plot(t_mom_fs, E_centroid, "-", color=palette_sweep5[0], linewidth=0.8,
              label=r"$\langle p\rangle^2 / 2m$")
    ax_b.plot(t_mom_fs, E_spread, "-", color=palette_sweep5[2], linewidth=0.8,
              label=r"$\sigma_p^2 / 2m$")
    ax_b.plot(t_mom_fs, E_total_wp, "-", color="black", linewidth=1.0,
              label=r"Sum")

    ax_b.set_xlabel(r"Time (fs)")
    ax_b.set_ylabel(r"Energy (eV)")
    ax_b.legend(fontsize=5.5, loc="best", frameon=False)
    ax_b.set_title("Wave-packet orbital only", fontsize=8, pad=3)
    panel_label(ax_b, "(b)")

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
