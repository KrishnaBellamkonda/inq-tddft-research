"""fig05 — Momentum-distribution waterfall (heatmap).

2D heatmap of the WP momentum distribution n_wp(k, t) from the jellium
σ=1 E=100 eV run. Shows momentum-space evolution: initial peak at k₀,
broadening and shifting after interaction with the bath.

Run:
    python -m inqview.report1.fig05_momentum_waterfall
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    TufteCritic,
)

CSV_PATH = (
    "ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1/"
    "results/raw/observables/momentum_distribution.csv"
)


def main() -> None:
    apply_style()

    df = pd.read_csv(CSV_PATH, comment="#")

    steps = sorted(df["step"].unique())
    k_vals = df[df["step"] == steps[0]]["k_bohr_inv"].values
    times = df[df["step"] == steps[0]]["time_au"].values
    # each step has one time value; get unique time per step
    time_per_step = df.groupby("step")["time_au"].first().values

    # build 2D array: n_wp(k, t)
    n_k = len(k_vals)
    n_t = len(steps)
    carpet = np.zeros((n_t, n_k))
    for i, step in enumerate(steps):
        mask = df["step"] == step
        carpet[i, :] = df.loc[mask, "n_wp"].values

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 1.0))

    vmax = carpet.max()
    vmin = max(1e-4 * vmax, carpet[carpet > 0].min()) if (carpet > 0).any() else 1e-10

    im = ax.pcolormesh(k_vals, time_per_step, carpet, cmap="viridis",
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading="auto", rasterized=True)

    # initial momentum reference
    k0_idx = np.argmax(carpet[0, :])
    k0 = k_vals[k0_idx]
    ax.axvline(k0, color="white", linewidth=0.5, linestyle=":", alpha=0.7)
    ax.text(k0 + 0.15, time_per_step[-1] * 0.95,
            fr"$k_0 = {k0:.1f}$ Bohr$^{{-1}}$",
            fontsize=6, color="white", va="top")

    ax.set_xlabel(r"$k$ (Bohr$^{-1}$)")
    ax.set_ylabel(r"$t$ (a.u.)")

    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(r"$n_{\mathrm{WP}}(k, t)$", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig05_momentum_waterfall.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
