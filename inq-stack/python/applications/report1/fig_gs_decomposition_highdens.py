"""fig_gs_decomposition_highdens — GS-basis decomposition bar plot.

Bar chart of δn_i^GS(t_end) for the high-density jellium WP run,
showing depletion of occupied states and excitation of virtual states.

Uses the existing effective_gs_occupations_t_end.csv from analyse.py.

Run:
    python -m applications.report1.fig_gs_decomposition_highdens
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from applications.report1._shared_style import (
    apply_style, palette_sweep5, column_widths_in, panel_label, TufteCritic,
)

RUN = Path("ResearchProject/systems/jellium/run_wp_n162_L30_E100_highdens_sigma1_v2")
GS_OCC_CSV = RUN / "results/analysis/observables/gs_projected_occupations/effective_gs_occupations_t_end.csv"
OUT = "docs/reports/report1/figures/fig_gs_decomposition_highdens.png"


def main() -> None:
    apply_style()

    df = pd.read_csv(GS_OCC_CSV)
    print(f"Loaded: {df.shape}, columns: {list(df.columns)}")

    # Handle various column naming conventions
    if "delta_n" in df.columns:
        state_idx = df["state_index"].values
        delta_n = df["delta_n"].values
    elif "delta" in df.columns:
        state_idx = df.get("state_index", df.get("gs_orbital_index")).values
        delta_n = df["delta"].values
    elif "n_gs_final" in df.columns:
        state_idx = df["state_index"].values
        delta_n = df["n_gs_final"].values - df["n_gs_initial"].values
    else:
        print(f"ERROR: unexpected columns: {list(df.columns)}")
        return

    n_occ = 81  # N=162, spin-degenerate

    W = column_widths_in["full"]
    fig, ax = plt.subplots(figsize=(W, W * 0.35))

    colors = np.where(state_idx < n_occ, palette_sweep5[0], palette_sweep5[3])
    ax.bar(state_idx, delta_n * 1e3, color=colors, width=0.8,
           edgecolor="none", zorder=3)

    ax.axhline(0, color="#b0b0b0", linewidth=0.4)
    ax.axvline(n_occ - 0.5, color="#C07020", linewidth=0.8, linestyle="--",
               label=f"Fermi level ($n_{{occ}}={n_occ}$)")

    # Annotate top depleted and excited
    sorted_idx = np.argsort(delta_n)
    for j in sorted_idx[:3]:
        if delta_n[j] < -1e-5:
            ax.text(state_idx[j], delta_n[j] * 1e3 - 0.3,
                    str(state_idx[j]), fontsize=4.5, ha="center", va="top",
                    color=palette_sweep5[0])
    for j in sorted_idx[-3:]:
        if delta_n[j] > 1e-5:
            ax.text(state_idx[j], delta_n[j] * 1e3 + 0.3,
                    str(state_idx[j]), fontsize=4.5, ha="center", va="bottom",
                    color=palette_sweep5[3])

    # Summary annotation
    sum_occ = delta_n[state_idx < n_occ].sum()
    sum_virt = delta_n[state_idx >= n_occ].sum()
    ax.text(0.5, 0.02,
            (rf"$\Sigma_{{\mathrm{{occ}}}}\delta n_i = {sum_occ*1e3:.1f}\times10^{{-3}}$"
             rf"$\quad\Sigma_{{\mathrm{{virt}}}}\delta n_i = {sum_virt*1e3:+.1f}\times10^{{-3}}$"),
            transform=ax.transAxes, fontsize=6, ha="center", va="bottom",
            bbox=dict(facecolor="white", edgecolor="#b0b0b0", linewidth=0.3, pad=2, alpha=0.9))

    ax.set_xlabel(r"GS orbital index $i$")
    ax.set_ylabel(r"$\delta n_i^{\mathrm{GS}}(t_{\mathrm{end}})$ [$\times 10^{-3}$]")
    ax.legend(fontsize=6, loc="upper left", frameon=False)

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
