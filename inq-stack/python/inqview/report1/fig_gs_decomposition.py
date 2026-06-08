"""fig_gs_decomposition — GS-basis decomposition bar plot with conservation sidebar.

Bar chart of δn_i^GS(t_end) showing depletion of occupied and excitation
of virtual states. Right sidebar shows charge conservation:
  Σ_occ (depletion), Σ_virt (excitation), unaccounted (high-energy states).

Generates both L=50 (standard) and L=30 (high density) versions.

Run:
    python -m inqview.report1.fig_gs_decomposition
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

from inqview.report1._shared_style import (
    apply_style, palette_sweep5, column_widths_in, TufteCritic,
)

CONFIGS = [
    {
        "label": "L=50 standard density ($r_s = 5.69$)",
        "run": Path("ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1"),
        "out": "docs/reports/report1/figures/fig_gs_decomposition.png",
        "L": 50.0, "dx": 0.40, "N_el": 162, "n_occ": 81,
    },
    {
        "label": "L=30 high density ($r_s = 3.41$)",
        "run": Path("ResearchProject/systems/jellium/run_wp_n162_L30_E100_highdens_sigma1_v2"),
        "out": "docs/reports/report1/figures/fig_gs_decomposition_highdens.png",
        "L": 30.0, "dx": 0.40, "N_el": 162, "n_occ": 81,
    },
]


def count_pw_states(L: float, dx: float) -> int:
    """Count plane-wave states within the grid energy cutoff."""
    k_cut = np.pi / dx
    E_cut = k_cut**2 / 2
    k_unit = 2 * np.pi / L
    n_max = int(k_cut / k_unit) + 1
    count = 0
    for nx in range(-n_max, n_max + 1):
        for ny in range(-n_max, n_max + 1):
            for nz in range(-n_max, n_max + 1):
                if (k_unit**2 * (nx**2 + ny**2 + nz**2)) / 2 <= E_cut:
                    count += 1
    return count


def make_plot(cfg: dict) -> None:
    apply_style()

    run = cfg["run"]
    csv_path = run / "results/analysis/observables/gs_projected_occupations/effective_gs_occupations_t_end.csv"
    if not csv_path.exists():
        print(f"  SKIP {cfg['out']}: no CSV at {csv_path}")
        return

    df = pd.read_csv(csv_path)
    state_idx = df.get("state_index", df.get("gs_orbital_index")).values
    delta_n = df["delta"].values
    n_occ = cfg["n_occ"]

    # Conservation accounting
    sum_occ = delta_n[state_idx < n_occ].sum()
    sum_virt = delta_n[state_idx >= n_occ].sum()
    sum_total = delta_n.sum()
    unaccounted = -sum_total

    n_tracked = len(state_idx)
    n_pw_total = count_pw_states(cfg["L"], cfg["dx"])

    print(f"  {cfg['label']}:")
    print(f"    Σ_occ = {sum_occ*1e3:.2f}×10⁻³,  Σ_virt = {sum_virt*1e3:+.2f}×10⁻³")
    print(f"    total δn = {sum_total*1e3:+.2f}×10⁻³,  unaccounted = {unaccounted*1e3:+.2f}×10⁻³")
    print(f"    Tracked: {n_tracked} / {n_pw_total} PW states ({n_tracked/n_pw_total*100:.2f}%)")

    # Figure with GridSpec: main bar + charge balance sidebar
    W = column_widths_in["full"]
    fig = plt.figure(figsize=(W, W * 0.38))
    gs = GridSpec(1, 2, width_ratios=[5.5, 1], wspace=0.12, figure=fig)
    ax = fig.add_subplot(gs[0])
    ax_side = fig.add_subplot(gs[1])

    # ── Main bar chart ──
    # Extend x-axis to show untracked region
    # Compute how many eigenstates exist at each shell up to cutoff
    k_cut = np.pi / cfg["dx"]
    E_cut_eV = (k_cut**2 / 2) * 27.2114
    k_unit = 2 * np.pi / cfg["L"]

    # Build shell-count table for annotation
    shell_counts = {}
    for nx in range(-int(k_cut / k_unit) - 1, int(k_cut / k_unit) + 2):
        for ny in range(-int(k_cut / k_unit) - 1, int(k_cut / k_unit) + 2):
            for nz in range(-int(k_cut / k_unit) - 1, int(k_cut / k_unit) + 2):
                g2 = nx**2 + ny**2 + nz**2
                E = 0.5 * g2 * k_unit**2 * 27.2114
                if E <= E_cut_eV:
                    shell_counts[g2] = shell_counts.get(g2, 0) + 1

    colors = np.where(state_idx < n_occ, palette_sweep5[0], palette_sweep5[3])
    ax.bar(state_idx, delta_n * 1e3, color=colors, width=0.8,
           edgecolor="none", zorder=3)
    ax.axhline(0, color="#b0b0b0", linewidth=0.4)
    ax.axvline(n_occ - 0.5, color="#C07020", linewidth=0.8, linestyle="--",
               label=f"Fermi level ($n_{{\\mathrm{{occ}}}}={n_occ}$)")

    # Shade the untracked region beyond the last tracked state
    last_tracked = int(state_idx.max())
    ax.axvspan(last_tracked + 0.5, last_tracked + 25, alpha=0.06,
               color="#808080", zorder=0)
    ax.annotate(
        rf"{n_pw_total - n_tracked:,} untracked states" "\n"
        rf"($E < E_{{\mathrm{{cut}}}} = {E_cut_eV:.0f}$ eV)",
        xy=(last_tracked + 1, 0), xytext=(last_tracked + 5, delta_n.min() * 1e3 * 0.5),
        fontsize=5, color="#808080", ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color="#808080", lw=0.5),
    )

    # Annotate top 3 depleted and excited
    sorted_idx = np.argsort(delta_n)
    for j in sorted_idx[:3]:
        if delta_n[j] < -1e-5:
            ax.text(state_idx[j], delta_n[j] * 1e3 - 0.3,
                    str(state_idx[j]), fontsize=4, ha="center", va="top",
                    color=palette_sweep5[0])
    for j in sorted_idx[-3:]:
        if delta_n[j] > 1e-5:
            ax.text(state_idx[j], delta_n[j] * 1e3 + 0.3,
                    str(state_idx[j]), fontsize=4, ha="center", va="bottom",
                    color=palette_sweep5[3])

    ax.set_xlim(-1, last_tracked + 26)
    ax.set_xlabel(r"GS orbital index $i$")
    ax.set_ylabel(r"$\delta n_i^{\mathrm{GS}}(t_{\mathrm{end}})$ [$\times 10^{-3}$]")
    ax.legend(fontsize=5.5, loc="lower left", frameon=False)

    # ── Charge balance sidebar ──
    bar_vals = [sum_occ * 1e3, sum_virt * 1e3, unaccounted * 1e3]
    bar_colors = [palette_sweep5[0], palette_sweep5[3], "#808080"]
    bar_labels = [r"$\Sigma_{\mathrm{occ}}$", r"$\Sigma_{\mathrm{virt}}$", "untracked"]
    x_bar = [0, 1, 2]

    bars = ax_side.bar(x_bar, bar_vals, color=bar_colors, width=0.65,
                       edgecolor="white", linewidth=0.3)
    ax_side.axhline(0, color="#b0b0b0", linewidth=0.4)

    for xi, val in zip(x_bar, bar_vals):
        va = "bottom" if val >= 0 else "top"
        offset = abs(max(bar_vals) - min(bar_vals)) * 0.05
        y_pos = val + offset if val >= 0 else val - offset
        ax_side.text(xi, y_pos, f"{val:+.1f}", fontsize=6, fontweight="bold",
                     ha="center", va=va)

    ax_side.set_xticks(x_bar)
    ax_side.set_xticklabels(bar_labels, fontsize=5, rotation=35, ha="right")
    ax_side.set_ylabel(r"$\Sigma\,\delta n_i$ [$\times 10^{-3}$]", fontsize=7)
    ax_side.set_title("Charge balance", fontsize=7, pad=8)
    ax_side.yaxis.tick_right()
    ax_side.yaxis.set_label_position("right")

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"    TufteCritic: {iss}")

    fig.savefig(cfg["out"], dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"  Saved -> {cfg['out']}")
    plt.close(fig)


def main() -> None:
    for cfg in CONFIGS:
        make_plot(cfg)


if __name__ == "__main__":
    main()
