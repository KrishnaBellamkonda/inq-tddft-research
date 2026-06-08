"""fig13 — Scenario A vs B schematic.

Two-panel conceptual schematic: elastic backscattering (A) vs inelastic
forward transmission (B).  Each panel has real-space and momentum-space
sub-views.  Pure matplotlib line art, no data.

Run:
    python -m inqview.report1.fig13_scenario_ab
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from inqview.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
)


def _gaussian(x, mu, sigma, amp=1.0):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _draw_scenario(ax_real, ax_mom, scenario: str) -> None:
    """Draw one scenario (A or B) across its real-space and momentum-space axes."""
    c_before = palette_sweep5[4]  # deep navy
    c_after = palette_sweep5[0]   # wine red

    p0 = 3.0
    sigma_p = 0.5
    p = np.linspace(-5, 7, 500)
    x = np.linspace(-8, 8, 500)

    # ── real-space sub-panel ────────────────────────────────────────
    ax = ax_real

    # initial WP
    wp_init = _gaussian(x, -3, 1.2, 0.8)
    ax.fill_between(x, wp_init, color=c_before, alpha=0.25)
    ax.plot(x, wp_init, color=c_before, linewidth=0.8)
    ax.annotate("", xy=(-1, 0.45), xytext=(-3, 0.45),
                arrowprops=dict(arrowstyle="-|>", color=c_before, lw=1.0))
    ax.text(-2, 0.55, r"$+p_0$", fontsize=6, ha="center", color=c_before)

    if scenario == "A":
        # after: split into forward + backward
        wp_fwd = _gaussian(x, 3, 1.2, 0.4)
        wp_bwd = _gaussian(x, -3, 1.2, 0.4)
        ax.fill_between(x, wp_fwd, color=c_after, alpha=0.2)
        ax.plot(x, wp_fwd, color=c_after, linewidth=0.8, linestyle="--")
        ax.fill_between(x, wp_bwd, color=c_after, alpha=0.2)
        ax.plot(x, wp_bwd, color=c_after, linewidth=0.8, linestyle="--")
        ax.annotate("", xy=(5, 0.22), xytext=(3.5, 0.22),
                    arrowprops=dict(arrowstyle="-|>", color=c_after, lw=0.8))
        ax.text(4.3, 0.28, r"$+p_0$", fontsize=5, ha="center", color=c_after)
        ax.annotate("", xy=(-5, 0.22), xytext=(-3.5, 0.22),
                    arrowprops=dict(arrowstyle="-|>", color=c_after, lw=0.8))
        ax.text(-4.3, 0.28, r"$-p_0$", fontsize=5, ha="center", color=c_after)
        ax.text(0, -0.12, r"\textit{after}", fontsize=5, ha="center",
                color=c_after)
    else:
        # after: forward with reduced momentum
        wp_fwd = _gaussian(x, 4, 1.2, 0.8)
        ax.fill_between(x, wp_fwd, color=c_after, alpha=0.25)
        ax.plot(x, wp_fwd, color=c_after, linewidth=0.8, linestyle="--")
        ax.annotate("", xy=(6, 0.45), xytext=(4.5, 0.45),
                    arrowprops=dict(arrowstyle="-|>", color=c_after, lw=0.8))
        ax.text(5.3, 0.55, r"$+p_1$", fontsize=6, ha="center", color=c_after)

    # target line
    ax.axvline(0, color="#404040", linewidth=1.0, linestyle="-", zorder=0)
    ax.text(0, 0.88, "target", fontsize=5, ha="center", color="#606060")

    ax.set_xlim(-7, 7)
    ax.set_ylim(-0.15, 0.95)
    ax.axis("off")

    # ── momentum-space sub-panel ────────────────────────────────────
    ax = ax_mom

    # before
    g_before = _gaussian(p, p0, sigma_p)
    ax.plot(p, g_before, color=c_before, linewidth=0.9, label="before")
    ax.fill_between(p, g_before, color=c_before, alpha=0.15)

    if scenario == "A":
        g_fwd = _gaussian(p, p0, sigma_p, 0.5)
        g_bwd = _gaussian(p, -p0, sigma_p, 0.5)
        g_after = g_fwd + g_bwd
        ax.plot(p, g_after, color=c_after, linewidth=0.9, linestyle="--",
                label="after")
        ax.fill_between(p, g_after, color=c_after, alpha=0.12)
        # annotation
        ax.text(0, 0.65,
                (r"$\langle p \rangle = 0$" "\n"
                 r"$\sigma_p^2 = p_0^2 + \sigma_0^2$" "\n"
                 r"$\Delta E_{\mathrm{kin}} = 0$"),
                fontsize=5, ha="center", va="top", color="#404040",
                bbox=dict(facecolor="white", edgecolor="#c0c0c0",
                          linewidth=0.3, pad=2, alpha=0.9))
    else:
        p1 = 2.0
        g_after = _gaussian(p, p1, sigma_p)
        ax.plot(p, g_after, color=c_after, linewidth=0.9, linestyle="--",
                label="after")
        ax.fill_between(p, g_after, color=c_after, alpha=0.12)
        ax.text(3.5, 0.65,
                (r"$\langle p \rangle = p_1$" "\n"
                 r"$\sigma_p^2 = \sigma_0^2$" "\n"
                 r"$\Delta E_{\mathrm{kin}} < 0$"),
                fontsize=5, ha="center", va="top", color="#404040",
                bbox=dict(facecolor="white", edgecolor="#c0c0c0",
                          linewidth=0.3, pad=2, alpha=0.9))

    ax.set_xlabel(r"$p_z$ (a.u.)", fontsize=7)
    ax.set_ylabel(r"$|\psi(p)|^2$", fontsize=7)
    ax.set_xlim(-5, 6)
    ax.set_ylim(-0.05, 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    apply_style()

    W = column_widths_in["full"]
    fig = plt.figure(figsize=(W, W * 0.42))

    outer = gridspec.GridSpec(1, 2, wspace=0.35)

    for col, (scenario, title) in enumerate([
        ("A", r"Scenario A: elastic redistribution"),
        ("B", r"Scenario B: inelastic transmission"),
    ]):
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[1, 1.2], hspace=0.15)
        ax_real = fig.add_subplot(inner[0])
        ax_mom = fig.add_subplot(inner[1])

        _draw_scenario(ax_real, ax_mom, scenario)

        lbl = "(a)" if col == 0 else "(b)"
        ax_real.text(0.02, 0.95, lbl, transform=ax_real.transAxes,
                     fontsize=8, va="top", ha="left")
        ax_real.text(0.5, 1.05, title, transform=ax_real.transAxes,
                     fontsize=6, ha="center", va="bottom", color="#404040")

    # bottom annotation
    fig.text(0.5, -0.02,
             (r"In jellium, $V_{\mathrm{ion}}(q \neq 0) = 0$ suppresses "
              r"Scenario A $\Rightarrow$ $\Delta E_{\mathrm{kin}}$ "
              r"is dominated by Scenario B."),
             fontsize=5.5, ha="center", va="top", color="#404040",
             bbox=dict(facecolor="#f8f8f8", edgecolor="#c0c0c0",
                       linewidth=0.3, pad=3))

    # stamp
    fig.text(0.99, -0.04, r"\textit{conceptual schematic}",
             fontsize=5, ha="right", va="top", color="#a0a0a0")

    out = "docs/reports/report1/figures/fig13_scenario_ab.png"
    with TufteCritic.disabled():
        fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.04)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
