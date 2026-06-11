"""fig_master_stopping — Master S(E) stopping power plot (both densities).

Combined figure with r_s=5.69 (L=50) solid and r_s=3.41 (L=30) dashed,
plus WP data points and Bethe asymptotes for both densities.

Produces:
  - fig_master_stopping.png       Combined, S_classical + WP S₂

Data: stopping_power_data pipeline (collect_L50_data, collect_L30_data).

Run:
    python -m applications.report1.fig_master_stopping
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    palette_regime3,
    column_widths_in,
    panel_label,
    references,
    TufteCritic,
)
from applications.report1.stopping_power_data import (
    MasterData,
    StoppingResult,
    bethe_curve,
    collect_L50_data,
    collect_L30_data,
)

HA_TO_EV = 27.2114

N_EL_L50 = 162 / 50.0**3
N_EL_L30 = 162 / 30.0**3

OUT_DIR = Path("docs/reports/report1/figures")


def _collect_classical(data: MasterData):
    """Extract (E, S, err) arrays from classical results, sorted by energy."""
    E_v2, S_v2, err_v2 = [], [], []
    E_v1, S_v1, err_v1 = [], [], []
    for r in sorted(data.classical, key=lambda x: x.energy_eV):
        S = r.S_classical
        if np.isfinite(S) and S > 0:
            if r.is_v2:
                E_v2.append(r.energy_eV)
                S_v2.append(S)
                err_v2.append(r.S_classical_err)
            else:
                E_v1.append(r.energy_eV)
                S_v1.append(S)
                err_v1.append(r.S_classical_err)
    return E_v2, S_v2, err_v2, E_v1, S_v1, err_v1


def _collect_wp(results: list[StoppingResult], S_attr: str = "S2_eV_per_bohr"):
    """Extract (E_v2, S_v2, E_v1, S_v1) from WP results."""
    v2_E, v2_S = [], []
    v1_E, v1_S = [], []
    for r in sorted(results, key=lambda x: x.energy_eV):
        S = getattr(r, S_attr)
        if not np.isfinite(S) or S <= 0:
            continue
        if r.is_v2:
            v2_E.append(r.energy_eV)
            v2_S.append(S)
        else:
            v1_E.append(r.energy_eV)
            v1_S.append(S)
    return v2_E, v2_S, v1_E, v1_S


def main() -> None:
    apply_style()

    print("Collecting L=50 data...")
    data_L50 = collect_L50_data()
    print("\nCollecting L=30 data...")
    data_L30 = collect_L30_data()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Colour assignments — L=50 (low density) in warm, L=30 (high density) in cool
    c_cl_L50 = "black"
    c_cl_L30 = "#606060"
    c_wp_s1_L50 = palette_sweep5[0]   # wine red
    c_wp_s5_L50 = palette_sweep5[3]   # steel blue
    c_wp_s1_L30 = palette_regime3[1]  # deep navy
    c_bethe_L50 = "#808080"
    c_bethe_L30 = "#A0A0A0"

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.65))

    S_attr = "S2_eV_per_bohr"

    # ═══════════════════════════════════════════════════════════════════
    # L=50 (r_s = 5.69) — SOLID lines
    # ═══════════════════════════════════════════════════════════════════

    # Classical
    cl_E_v2, cl_S_v2, cl_err_v2, cl_E_v1, cl_S_v1, cl_err_v1 = _collect_classical(data_L50)
    all_cl_E = cl_E_v2 + cl_E_v1
    all_cl_S = cl_S_v2 + cl_S_v1
    if all_cl_E:
        order = np.argsort(all_cl_E)
        ax.plot(np.array(all_cl_E)[order], np.array(all_cl_S)[order],
                "-", color=c_cl_L50, linewidth=0.9, zorder=3)
    if cl_E_v2:
        ax.errorbar(cl_E_v2, cl_S_v2, yerr=cl_err_v2,
                    fmt="s", color=c_cl_L50, markersize=4.5,
                    markeredgecolor="white", markeredgewidth=0.3,
                    capsize=2, capthick=0.6, elinewidth=0.6,
                    label=r"Classical $r_s{=}5.69$", zorder=4)
    if cl_E_v1:
        ax.plot(cl_E_v1, cl_S_v1, "s", color=c_cl_L50, markersize=4.5,
                markerfacecolor="white", markeredgecolor=c_cl_L50,
                markeredgewidth=0.8, zorder=4)

    # WP sigma=1
    s1_v2E, s1_v2S, s1_v1E, s1_v1S = _collect_wp(data_L50.wp_sigma1, S_attr)
    all_s1_E = s1_v2E + s1_v1E
    all_s1_S = s1_v2S + s1_v1S
    if all_s1_E:
        order = np.argsort(all_s1_E)
        ax.plot(np.array(all_s1_E)[order], np.array(all_s1_S)[order],
                "-", color=c_wp_s1_L50, linewidth=0.8, zorder=3)
    if s1_v2E:
        ax.plot(s1_v2E, s1_v2S, "o", color=c_wp_s1_L50, markersize=5,
                markeredgecolor="white", markeredgewidth=0.3,
                label=r"WP $\sigma{=}1$, $r_s{=}5.69$", zorder=5)
    if s1_v1E:
        ax.plot(s1_v1E, s1_v1S, "o", color=c_wp_s1_L50, markersize=5,
                markerfacecolor="white", markeredgecolor=c_wp_s1_L50,
                markeredgewidth=0.8, zorder=5)

    # WP sigma=5
    s5_v2E, s5_v2S, s5_v1E, s5_v1S = _collect_wp(data_L50.wp_sigma5, S_attr)
    all_s5_E = s5_v2E + s5_v1E
    all_s5_S = s5_v2S + s5_v1S
    if all_s5_E:
        order = np.argsort(all_s5_E)
        ax.plot(np.array(all_s5_E)[order], np.array(all_s5_S)[order],
                "-", color=c_wp_s5_L50, linewidth=0.8, zorder=3)
    if s5_v2E:
        ax.plot(s5_v2E, s5_v2S, "^", color=c_wp_s5_L50, markersize=5,
                markeredgecolor="white", markeredgewidth=0.3,
                label=r"WP $\sigma{=}5$, $r_s{=}5.69$", zorder=5)
    if s5_v1E:
        ax.plot(s5_v1E, s5_v1S, "^", color=c_wp_s5_L50, markersize=5,
                markerfacecolor="white", markeredgecolor=c_wp_s5_L50,
                markeredgewidth=0.8, zorder=5)

    # Bethe L=50
    E_grid = np.logspace(np.log10(15), np.log10(3000), 200)
    S_bethe_50 = bethe_curve(E_grid, N_EL_L50)
    mask50 = S_bethe_50 > 0
    ax.plot(E_grid[mask50], S_bethe_50[mask50],
            color=c_bethe_L50, linestyle="--", linewidth=0.9,
            label=r"Bethe $r_s{=}5.69$", zorder=2)

    # ═══════════════════════════════════════════════════════════════════
    # L=30 (r_s = 3.41) — DASHED lines
    # ═══════════════════════════════════════════════════════════════════

    # Classical
    hd_cl_E_v2, hd_cl_S_v2, hd_cl_err_v2, hd_cl_E_v1, hd_cl_S_v1, hd_cl_err_v1 = _collect_classical(data_L30)
    all_hd_cl_E = hd_cl_E_v2 + hd_cl_E_v1
    all_hd_cl_S = hd_cl_S_v2 + hd_cl_S_v1
    if all_hd_cl_E:
        order = np.argsort(all_hd_cl_E)
        ax.plot(np.array(all_hd_cl_E)[order], np.array(all_hd_cl_S)[order],
                "--", color=c_cl_L30, linewidth=0.9, zorder=3)
    if hd_cl_E_v2:
        ax.errorbar(hd_cl_E_v2, hd_cl_S_v2, yerr=hd_cl_err_v2,
                    fmt="D", color=c_cl_L30, markersize=4,
                    markeredgecolor="white", markeredgewidth=0.3,
                    capsize=2, capthick=0.6, elinewidth=0.6,
                    label=r"Classical $r_s{=}3.41$", zorder=4)
    if hd_cl_E_v1:
        ax.errorbar(hd_cl_E_v1, hd_cl_S_v1, yerr=hd_cl_err_v1,
                    fmt="D", color=c_cl_L30, markersize=4,
                    markerfacecolor="white", markeredgecolor=c_cl_L30,
                    markeredgewidth=0.8,
                    capsize=2, capthick=0.6, elinewidth=0.6,
                    label=r"Classical $r_s{=}3.41$" if not hd_cl_E_v2 else None,
                    zorder=4)

    # WP sigma=1 (high density)
    hd_s1_v2E, hd_s1_v2S, hd_s1_v1E, hd_s1_v1S = _collect_wp(data_L30.wp_sigma1, S_attr)
    all_hd_s1_E = hd_s1_v2E + hd_s1_v1E
    all_hd_s1_S = hd_s1_v2S + hd_s1_v1S
    if all_hd_s1_E:
        order = np.argsort(all_hd_s1_E)
        ax.plot(np.array(all_hd_s1_E)[order], np.array(all_hd_s1_S)[order],
                "--", color=c_wp_s1_L30, linewidth=0.8, zorder=3)
    if hd_s1_v2E:
        ax.plot(hd_s1_v2E, hd_s1_v2S, "o", color=c_wp_s1_L30, markersize=5,
                markeredgecolor="white", markeredgewidth=0.3,
                label=r"WP $\sigma{=}1$, $r_s{=}3.41$", zorder=5)
    if hd_s1_v1E:
        ax.plot(hd_s1_v1E, hd_s1_v1S, "o", color=c_wp_s1_L30, markersize=5,
                markerfacecolor="white", markeredgecolor=c_wp_s1_L30,
                markeredgewidth=0.8, zorder=5)

    # WP supplementary (high density)
    supp_styles = {
        0.5: ("D", palette_sweep5[2], r"WP $\sigma{=}0.5$, $r_s{=}3.41$"),
    }
    for r in data_L30.wp_supplementary:
        S = getattr(r, S_attr)
        if not np.isfinite(S) or S <= 0:
            continue
        marker, color, label = supp_styles.get(
            r.sigma, ("x", "gray", r"WP $\sigma{=}" + f"{r.sigma}" + r"$, $r_s{=}3.41$")
        )
        ax.plot(r.energy_eV, S, marker, color=color, markersize=5.5,
                markeredgecolor="white", markeredgewidth=0.3,
                label=label, zorder=5)

    # Bethe L=30
    S_bethe_30 = bethe_curve(E_grid, N_EL_L30)
    mask30 = S_bethe_30 > 0
    ax.plot(E_grid[mask30], S_bethe_30[mask30],
            color=c_bethe_L30, linestyle=":", linewidth=0.9,
            label=r"Bethe $r_s{=}3.41$", zorder=2)

    # ═══════════════════════════════════════════════════════════════════
    # Axes formatting
    # ═══════════════════════════════════════════════════════════════════

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Projectile energy $E$ (eV)")
    ax.set_ylabel(r"Stopping power $S$ (eV/Bohr)")

    ax.legend(fontsize=6, loc="upper right", frameon=True,
              edgecolor="#b0b0b0", fancybox=False, ncol=1)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        for iss in issues:
            print(f"  TufteCritic: {iss}")

    out = OUT_DIR / "fig_master_stopping.png"
    fig.savefig(str(out), dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"\nSaved -> {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
