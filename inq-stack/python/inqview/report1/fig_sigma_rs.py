"""fig_sigma_rs — The σ/r_s scaling rule.

|ΔE_WP| / |ΔE_cl| vs σ/r_s for two densities (r_s = 5.69 and 3.41).
Shows the quantum-to-classical crossover controlled by σ/r_s.

Uses the stopping_power_data pipeline for interference-free windows.

Run:
    python -m inqview.report1.fig_sigma_rs
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    TufteCritic,
)
from inqview.report1.stopping_power_data import (
    collect_L50_data,
    collect_L30_data,
    MasterData,
)

RS_569 = 5.69  # L=50, N=162
RS_341 = 3.41  # L=30, N=162

OUT = "docs/reports/report1/figures/fig_sigma_rs.png"


def _extract_ratio_points(data: MasterData, r_s: float, E_target: float = 100.0):
    """Extract |ΔE_WP|/|ΔE_cl| vs σ/r_s at a given energy.

    Uses S₂ (KS orbital) definition for both WP and classical.
    """
    # Find classical reference at E_target
    dE_cl = None
    for r in data.classical:
        if abs(r.energy_eV - E_target) < 1 and np.isfinite(r.S_classical):
            dE_cl = abs(r.S_classical * r.dz_bohr) if np.isfinite(r.dz_bohr) and abs(r.dz_bohr) > 0.1 else None
            # Actually use S_classical * matched_dz... but we want |ΔE|
            # Better: use the stopping power ratio directly: S_WP / S_cl
            S_cl = r.S_classical
            break

    if dE_cl is None and S_cl is None:
        return [], [], [], []

    sigma_over_rs = []
    ratios = []
    is_v2_list = []
    is_compromised_list = []

    # Collect all WP runs at E_target across all sigma values
    all_wp = list(data.wp_sigma1) + list(data.wp_sigma5) + list(data.wp_supplementary)
    for r in all_wp:
        if abs(r.energy_eV - E_target) > 1:
            continue
        S_wp = r.S2_eV_per_bohr
        if not np.isfinite(S_wp) or S_wp <= 0:
            continue

        ratio = S_wp / S_cl
        sigma_over_rs.append(r.sigma / r_s)
        ratios.append(ratio)
        is_v2_list.append(r.is_v2)
        is_compromised_list.append(r.is_boundary_compromised)

        tag = "v2" if r.is_v2 else "v1"
        flag = " [COMPROMISED]" if r.is_boundary_compromised else ""
        print(f"    σ={r.sigma:.1f} → σ/r_s={r.sigma/r_s:.3f}, "
              f"S_WP/S_cl={ratio:.4f} ({tag}){flag}")

    return sigma_over_rs, ratios, is_v2_list, is_compromised_list


def main() -> None:
    apply_style()

    print("Collecting L=50 data...")
    data_L50 = collect_L50_data()
    print("\nCollecting L=30 data...")
    data_L30 = collect_L30_data()

    print(f"\n  σ/r_s scaling — r_s = {RS_569}:")
    sr_569, rat_569, v2_569, comp_569 = _extract_ratio_points(data_L50, RS_569)

    print(f"\n  σ/r_s scaling — r_s = {RS_341}:")
    sr_341, rat_341, v2_341, comp_341 = _extract_ratio_points(data_L30, RS_341)

    # Plot
    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.78))

    c1 = palette_sweep5[0]
    c2 = palette_sweep5[3]

    # r_s = 5.69 points
    if sr_569:
        order = np.argsort(sr_569)
        sr_arr = np.array(sr_569)[order]
        rat_arr = np.array(rat_569)[order]
        v2_arr = np.array(v2_569)[order]
        comp_arr = np.array(comp_569)[order]

        ax.plot(sr_arr, rat_arr, "-", color=c1, linewidth=0.9, zorder=2)

        # v2 filled, v1 open
        for i in range(len(sr_arr)):
            if comp_arr[i]:
                ax.plot(sr_arr[i], rat_arr[i], "o", color=c1, markersize=4,
                        markerfacecolor="white", markeredgecolor=c1,
                        markeredgewidth=0.8, zorder=4)
                ax.plot(sr_arr[i], rat_arr[i], "x", color="red",
                        markersize=5, markeredgewidth=0.8, zorder=5)
            elif v2_arr[i]:
                ax.plot(sr_arr[i], rat_arr[i], "o", color=c1, markersize=4,
                        markeredgecolor="white", markeredgewidth=0.3, zorder=4)
            else:
                ax.plot(sr_arr[i], rat_arr[i], "o", color=c1, markersize=4,
                        markerfacecolor="white", markeredgecolor=c1,
                        markeredgewidth=0.8, zorder=4)

        # Single legend entry
        ax.plot([], [], "o-", color=c1, markersize=4, linewidth=0.9,
                label=f"$r_s = {RS_569}$")

    # r_s = 3.41 points
    if sr_341:
        order = np.argsort(sr_341)
        sr_arr = np.array(sr_341)[order]
        rat_arr = np.array(rat_341)[order]
        v2_arr = np.array(v2_341)[order]

        if len(sr_arr) > 1:
            ax.plot(sr_arr, rat_arr, "-", color=c2, linewidth=0.9, zorder=2)

        for i in range(len(sr_arr)):
            if v2_arr[i]:
                ax.plot(sr_arr[i], rat_arr[i], "s", color=c2, markersize=4,
                        markeredgecolor="white", markeredgewidth=0.3, zorder=4)
            else:
                ax.plot(sr_arr[i], rat_arr[i], "s", color=c2, markersize=4,
                        markerfacecolor="white", markeredgecolor=c2,
                        markeredgewidth=0.8, zorder=4)

        ax.plot([], [], "s-", color=c2, markersize=4, linewidth=0.9,
                label=f"$r_s = {RS_341}$")

    ax.axhline(1.0, color="#b0b0b0", linewidth=0.5, ls="--", zorder=1)
    ax.text(0.05, 1.03, "classical agreement", fontsize=5.5, color="#808080")

    ax.set_xlabel(r"$\sigma / r_s$")
    ax.set_ylabel(r"$S_{\mathrm{WP}} \,/\, S_{\mathrm{cl}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=6, loc="upper right", frameon=False)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"\n  TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"    {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"\n  Saved → {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
