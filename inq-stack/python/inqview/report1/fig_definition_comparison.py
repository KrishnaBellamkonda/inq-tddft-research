"""fig_definition_comparison — Stopping-power definition comparison.

Two-panel bar chart: (a) σ=1 Bohr, (b) σ=5 Bohr, showing WP stopping
power definitions vs classical reference at E=100 eV.

Uses the stopping_power_data pipeline for interference-free windows
and matched time periods.

Generates both L=50 (standard) and L=30 (high density) versions.

Run:
    python -m inqview.report1.fig_definition_comparison
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from inqview.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
)
from inqview.report1.stopping_power_data import (
    collect_L50_data,
    collect_L30_data,
    MasterData,
    HA_TO_EV,
)

OUT_L50 = "docs/reports/report1/figures/fig_definition_comparison.png"
OUT_L30 = "docs/reports/report1/figures/fig_definition_comparison_highdens.png"


def _find_at_E100(results, sigma_target: float):
    """Find the StoppingResult at E=100 eV closest to sigma_target."""
    for r in results:
        if abs(r.energy_eV - 100) < 1 and abs(r.sigma - sigma_target) < 0.1:
            return r
    return None


def _find_classical_E100(results):
    for r in results:
        if abs(r.energy_eV - 100) < 1:
            return r
    return None


def _extract_definitions(r) -> dict:
    """Extract all available stopping power definitions from a StoppingResult."""
    defs = {}
    if r is None:
        return defs
    if np.isfinite(r.S1_eV_per_bohr):
        defs["Momentum"] = r.S1_eV_per_bohr
    if np.isfinite(r.S2_eV_per_bohr):
        defs["KS orbital"] = r.S2_eV_per_bohr
    return defs


def _make_figure(data: MasterData, *, density_label: str, out: str):
    apply_style()

    # Find E=100 data points
    r_s1 = _find_at_E100(data.wp_sigma1, 1.0)
    r_s5 = _find_at_E100(data.wp_sigma5, 5.0)
    # Also check supplementary for σ=0.5, 3, 8
    if r_s1 is None:
        r_s1 = _find_at_E100(data.wp_supplementary, 1.0)
    if r_s5 is None:
        r_s5 = _find_at_E100(data.wp_supplementary, 5.0)

    r_cl = _find_classical_E100(data.classical)

    S_cl = r_cl.S_classical if r_cl and np.isfinite(r_cl.S_classical) else np.nan
    S_cl_err = r_cl.S_classical_err if r_cl and np.isfinite(r_cl.S_classical_err) else 0

    defs_s1 = _extract_definitions(r_s1)
    defs_s5 = _extract_definitions(r_s5)

    if not defs_s1 and not defs_s5:
        print(f"  No E=100 WP data found for {density_label} — skipping")
        return

    # Print provenance
    print(f"\n  Definition comparison ({density_label}):")
    if r_s1:
        tag = "v2" if r_s1.is_v2 else "v1"
        print(f"    σ=1: {r_s1.run_dir} ({tag}), window=[0, {r_s1.window.t_end:.2f}] a.u.")
        for k, v in defs_s1.items():
            print(f"      {k}: S = {v:.5f} eV/Bohr")
    if r_s5:
        tag = "v2" if r_s5.is_v2 else "v1"
        print(f"    σ=5: {r_s5.run_dir} ({tag}), window=[0, {r_s5.window.t_end:.2f}] a.u.")
        for k, v in defs_s5.items():
            print(f"      {k}: S = {v:.5f} eV/Bohr")
    if r_cl:
        tag = "v2" if r_cl.is_v2 else "v1"
        print(f"    Classical: {r_cl.run_dir} ({tag})")
        print(f"      S = {S_cl:.5f} ± {S_cl_err:.5f} eV/Bohr")

    # Build bar chart
    all_labels = sorted(set(list(defs_s1.keys()) + list(defs_s5.keys())))
    all_labels.append("Classical")
    colors_map = {
        "Momentum": palette_sweep5[0],
        "KS orbital": palette_sweep5[2],
        "Classical": "black",
    }
    colors = [colors_map.get(l, palette_sweep5[3]) for l in all_labels]

    has_s5 = bool(defs_s5)
    n_panels = 2 if has_s5 else 1

    W = column_widths_in["full"]
    if n_panels == 2:
        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(W, W * 0.42), sharey=True)
        panels = [(ax_a, defs_s1, r"(a) $\sigma = 1$ Bohr"),
                  (ax_b, defs_s5, r"(b) $\sigma = 5$ Bohr")]
    else:
        fig, ax_a = plt.subplots(1, 1, figsize=(W * 0.55, W * 0.42))
        panels = [(ax_a, defs_s1, r"(a) $\sigma = 1$ Bohr")]

    for ax, defs, plabel in panels:
        vals = [defs.get(l, np.nan) for l in all_labels[:-1]] + [S_cl]
        errs = [0] * (len(all_labels) - 1) + [S_cl_err]
        x = np.arange(len(all_labels))

        bar_vals = [v if np.isfinite(v) else 0.0 for v in vals]
        bar_errs = [e if np.isfinite(e) else 0.0 for e in errs]

        bars = ax.bar(x, bar_vals, yerr=bar_errs, color=colors, width=0.7,
                      edgecolor="white", linewidth=0.3, zorder=3,
                      capsize=2, error_kw={"elinewidth": 0.6, "capthick": 0.5})

        for j, v in enumerate(vals):
            if not np.isfinite(v):
                bars[j].set_hatch("///")
                bars[j].set_facecolor("#d0d0d0")
                bars[j].set_edgecolor("#808080")

        ax.set_xticks(x)
        ax.set_xticklabels(all_labels, fontsize=6, rotation=25, ha="right")
        ax.set_xlabel(r"SP definition", fontsize=7)
        if np.isfinite(S_cl):
            ax.axhline(S_cl, color="grey", ls="--", lw=0.5, zorder=1)
        panel_label(ax, plabel)

    panels[0][0].set_ylabel(r"$S$ (eV/Bohr)")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"  TufteCritic ({out}): {len(issues)} issue(s)")
        for iss in issues:
            print(f"    {iss}")

    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"  Saved → {out}")
    plt.close(fig)


def main() -> None:
    print("Collecting L=50 data...")
    data_L50 = collect_L50_data()
    _make_figure(data_L50, density_label="L=50, r_s=5.69", out=OUT_L50)

    print("\nCollecting L=30 data...")
    data_L30 = collect_L30_data()
    _make_figure(data_L30, density_label="L=30, r_s=3.41", out=OUT_L30)


if __name__ == "__main__":
    main()
