"""fig03 — Jellium-specific stopping-power methods (tabular figure).

9 methods × 4 columns with colour-coded approximation-class badges.
Pure table rendered as a matplotlib figure — no axes, no data curves.

Run:
    python -m applications.report1.fig03_jellium_methods
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from applications.report1._shared_style import (
    apply_style,
    palette_regime3,
    column_widths_in,
    TufteCritic,
)

# ── table data (from T3_schematic_data.md Fig 3c) ───────────────────
# (method, observable, validity, approx_class)
# approx_class: "linear" | "non-linear" | "non-perturbative"
ROWS = [
    ("Linear Lindhard RPA",
     r"$-\mathrm{Im}[1/\varepsilon(q,\omega)]$",
     "All $v$; linear in $Z_1$",
     "linear"),
    ("Mermin (1970)",
     r"Damped $\varepsilon$ with $\tau$",
     "All $v$; number-conserving",
     "linear"),
    (r"Beyond-RPA $G(q)$ (Corradini 1998)",
     r"Static $G(q)$",
     "Static; xc in $q$",
     "linear"),
    ("GS-DFT (Echenique 1981)",
     r"$\mathrm{d}E/\mathrm{d}x$ (phase shifts)",
     r"$v \ll v_F$; self-consistent",
     "non-linear"),
    ("Non-linear DFT (Echenique 1986)",
     r"$Q(v_F, Z_1)$ friction",
     r"$v \ll v_F$; bound states",
     "non-linear"),
    ("OF-DFT (Moldabekov 2023)",
     r"Density response (TF+)",
     "Requires UEG-exact KE",
     "non-linear"),
    ("Vignale--Kohn TDCDFT (1996)",
     r"Current response",
     "Linear; memory via $f_{xc}^L$",
     "linear"),
    ("rt-TDDFT in jellium (project)",
     r"Full $S(v)$ curve",
     "All $v$; ALDA xc",
     "non-perturbative"),
    ("Nazarov--Gross (2025)",
     r"$S(v)$, $Q_1$ friction",
     "Quantum projectile; exact fact.",
     "non-perturbative"),
]

BADGE_COLORS = {
    "linear":           palette_regime3[1],  # deep navy
    "non-linear":       "#C07020",           # warm orange
    "non-perturbative": palette_regime3[2],  # forest green
}

HEADER = ("Method", "Observable", "Validity / approximation", "Class")
PROJECT_ROWS = {7, 8}  # 0-indexed: rt-TDDFT and Nazarov-Gross


def main() -> None:
    apply_style()

    n_rows = len(ROWS)
    n_cols = len(HEADER)

    W = column_widths_in["full"]
    row_h = 0.30
    header_h = 0.35
    fig_h = header_h + n_rows * row_h + 0.15
    fig, ax = plt.subplots(figsize=(W, fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, fig_h)
    ax.axis("off")

    col_x = [0.0, 0.32, 0.58, 0.84]
    col_w = [0.31, 0.25, 0.25, 0.15]

    y_top = fig_h - 0.08

    # header
    for j, (hdr, xj) in enumerate(zip(HEADER, col_x)):
        ax.text(xj + 0.01, y_top, r"\textbf{" + hdr + "}",
                fontsize=6.5, va="top", ha="left", color="#1a1a1a")
    y_line = y_top - header_h + 0.08
    ax.plot([0, 1], [y_line, y_line], color="#1a1a1a", linewidth=0.6,
            clip_on=False)

    # rows
    for i, (method, obs, validity, approx) in enumerate(ROWS):
        y = y_line - (i + 0.5) * row_h

        # alternating background
        if i % 2 == 1:
            ax.add_patch(plt.Rectangle(
                (0, y - row_h * 0.45), 1, row_h * 0.9,
                facecolor="#f5f5f5", edgecolor="none", zorder=0))

        # project-relevant highlight
        if i in PROJECT_ROWS:
            ax.add_patch(plt.Rectangle(
                (0, y - row_h * 0.45), 1, row_h * 0.9,
                facecolor="#FFFDE0", edgecolor="none", zorder=0))

        fs = 5.8
        ax.text(col_x[0] + 0.01, y, method, fontsize=fs, va="center")
        ax.text(col_x[1] + 0.01, y, obs, fontsize=fs, va="center")
        ax.text(col_x[2] + 0.01, y, validity, fontsize=fs, va="center")

        # approximation-class badge
        badge_color = BADGE_COLORS[approx]
        badge_label = approx.replace("-", "-\n") if len(approx) > 10 else approx
        bx = col_x[3] + 0.01
        badge = FancyBboxPatch(
            (bx, y - 0.08), 0.13, 0.16,
            boxstyle="round,pad=0.02",
            facecolor=badge_color, edgecolor="none", alpha=0.85,
            zorder=2, clip_on=False,
        )
        ax.add_patch(badge)
        ax.text(bx + 0.065, y, approx,
                fontsize=4.5, va="center", ha="center",
                color="white", fontweight="bold", zorder=3)

    # schematic stamp
    ax.text(0.99, 0.02, r"\textit{SCHEMATIC}",
            transform=ax.transAxes, fontsize=5, color="#a0a0a0",
            ha="right", va="bottom")

    critic = TufteCritic()
    with critic.disabled():
        pass

    out = "docs/reports/report1/figures/fig03_jellium_methods.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
