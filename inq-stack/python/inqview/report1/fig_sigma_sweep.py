"""fig_sigma_sweep — σ-dependence of WP stopping at E=100 eV.

(a) |ΔE_KS| vs σ on log-log — absolute energy loss of the WP orbital.
(b) Δσ_pz²/σ_pz²(0) vs σ — fractional momentum-variance change,
    a measure of momentum-space coupling / heating by the bath.

All quantities windowed to the interference-free interval from the
stopping_power_data pipeline.

Run:
    python -m inqview.report1.fig_sigma_sweep
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
    MasterData,
    StoppingResult,
    parse_run_summary,
    compute_time_window,
    HA_TO_EV,
)

OUT = "docs/reports/report1/figures/fig_sigma_sweep.png"


def _momentum_variance_change(run_dir: str, t_end: float) -> tuple[float, float, float]:
    """Compute windowed Δσ_pz²/σ_pz²(0) from wp_momentum_stats.csv.

    Returns (sigma_pz2_0, sigma_pz2_end, fractional_change).
    """
    path = Path(run_dir) / "results" / "raw" / "observables" / "wp_momentum_stats.csv"
    if not path.exists():
        return (np.nan, np.nan, np.nan)

    df = pd.read_csv(path, comment="#")
    mask = df["time_au"] <= t_end
    df_win = df[mask]
    if len(df_win) < 2:
        return (np.nan, np.nan, np.nan)

    sp0 = df_win["sigma_pz2"].iloc[0]
    spf = df_win["sigma_pz2"].iloc[-1]
    frac = (spf - sp0) / sp0 if sp0 > 1e-10 else np.nan
    return (sp0, spf, frac)


def main() -> None:
    apply_style()

    print("Collecting L=50 data...")
    data = collect_L50_data()

    # Collect all E=100 WP runs across all σ values
    all_wp = list(data.wp_sigma1) + list(data.wp_sigma5) + list(data.wp_supplementary)
    e100_runs: list[StoppingResult] = [
        r for r in all_wp if abs(r.energy_eV - 100) < 1
    ]
    e100_runs.sort(key=lambda r: r.sigma)

    # Extract panel (a) data: |ΔE_KS| and panel (b): Δσ_pz²/σ_pz²(0)
    sigmas = []
    dE_abs = []
    frac_dsigma = []
    is_v2 = []
    is_compromised = []

    print("\n  σ-sweep at E=100 eV:")
    for r in e100_runs:
        if r.window is None:
            continue
        dE = abs(r.dE2_eV) if np.isfinite(r.dE2_eV) else np.nan
        if not np.isfinite(dE):
            continue

        sp0, spf, frac = _momentum_variance_change(r.run_dir, r.window.t_end)

        sigmas.append(r.sigma)
        dE_abs.append(dE)
        frac_dsigma.append(frac)
        is_v2.append(r.is_v2)
        is_compromised.append(r.is_boundary_compromised)

        tag = "v2" if r.is_v2 else "v1"
        flag = " [BOUNDARY]" if r.is_boundary_compromised else ""
        print(
            f"    σ={r.sigma:.1f}: |ΔE_KS|={dE:.4f} eV, "
            f"Δσ_pz²/σ_pz²(0)={frac:.4f} ({tag}){flag}"
        )

    sigmas = np.array(sigmas)
    dE_abs = np.array(dE_abs)
    frac_dsigma = np.array(frac_dsigma)
    is_v2 = np.array(is_v2)
    is_compromised = np.array(is_compromised)

    # ── Plot ──
    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(W, W * 0.42),
        gridspec_kw={"wspace": 0.35},
    )

    c = palette_sweep5[0]

    # Panel (a): |ΔE_KS| vs σ
    ax_a.plot(sigmas, dE_abs, "-", color=c, linewidth=0.9, zorder=2)
    for i in range(len(sigmas)):
        kw = dict(markersize=5, zorder=4)
        if is_compromised[i]:
            kw.update(markerfacecolor="white", markeredgecolor=c, markeredgewidth=0.8)
            ax_a.plot(sigmas[i], dE_abs[i], "o", color=c, **kw)
            ax_a.plot(sigmas[i], dE_abs[i], "x", color="red",
                      markersize=6, markeredgewidth=0.8, zorder=5)
        elif is_v2[i]:
            kw.update(markeredgecolor="white", markeredgewidth=0.3)
            ax_a.plot(sigmas[i], dE_abs[i], "o", color=c, **kw)
        else:
            kw.update(markerfacecolor="white", markeredgecolor=c, markeredgewidth=0.8)
            ax_a.plot(sigmas[i], dE_abs[i], "o", color=c, **kw)

    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlabel(r"$\sigma$ (Bohr)")
    ax_a.set_ylabel(r"$|\Delta E_{\mathrm{KS}}|$ (eV)")
    panel_label(ax_a, "(a)")

    # Panel (b): Δσ_pz²/σ_pz²(0) vs σ
    valid = np.isfinite(frac_dsigma)
    ax_b.plot(sigmas[valid], frac_dsigma[valid], "-", color=c, linewidth=0.9, zorder=2)
    for i in range(len(sigmas)):
        if not np.isfinite(frac_dsigma[i]):
            continue
        kw = dict(markersize=5, zorder=4)
        if is_compromised[i]:
            kw.update(markerfacecolor="white", markeredgecolor=c, markeredgewidth=0.8)
            ax_b.plot(sigmas[i], frac_dsigma[i], "o", color=c, **kw)
            ax_b.plot(sigmas[i], frac_dsigma[i], "x", color="red",
                      markersize=6, markeredgewidth=0.8, zorder=5)
        elif is_v2[i]:
            kw.update(markeredgecolor="white", markeredgewidth=0.3)
            ax_b.plot(sigmas[i], frac_dsigma[i], "o", color=c, **kw)
        else:
            kw.update(markerfacecolor="white", markeredgecolor=c, markeredgewidth=0.8)
            ax_b.plot(sigmas[i], frac_dsigma[i], "o", color=c, **kw)

    ax_b.set_xscale("log")
    ax_b.set_xlabel(r"$\sigma$ (Bohr)")
    ax_b.set_ylabel(r"$\Delta\sigma_{p_z}^2 \,/\, \sigma_{p_z}^2(0)$")
    panel_label(ax_b, "(b)")

    # ── Save ──
    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        for iss in issues:
            print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"\n  Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
