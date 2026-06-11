"""fig_gantt — Report 2 plan timeline.

Gantt-style chart showing relative time-blocks for planned work.
No specific dates — just relative durations.

Run:
    python -m applications.report1.fig_gantt
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    TufteCritic,
)


def main() -> None:
    apply_style()

    tasks = [
        (r"Goal 1: $\sigma$-sweep at multiple $E$", 0, 4, palette_sweep5[0]),
        ("Goal 1: error-bar campaign", 2, 3, palette_sweep5[1]),
        ("Goal 2: induced-field-work definition", 3, 4, palette_sweep5[2]),
        ("Goal 2: channel decomposition", 5, 3, palette_sweep5[2]),
        (r"Goal 3: $\sigma \times v$ 2D map", 4, 5, palette_sweep5[3]),
        ("Coronene quantitative comparison", 1, 3, palette_sweep5[4]),
        ("Report 2 writing", 7, 3, "#808080"),
    ]

    W = column_widths_in["full"]
    fig, ax = plt.subplots(figsize=(W, W * 0.32))

    for i, (label, start, duration, color) in enumerate(tasks):
        ax.barh(i, duration, left=start, height=0.6, color=color,
                edgecolor="white", linewidth=0.3, alpha=0.85, zorder=2)
        ax.text(start + duration / 2, i, label, ha="center", va="center",
                fontsize=5, color="white", fontweight="bold", zorder=3)

    ax.set_yticks([])
    ax.set_xlabel("Relative time (weeks)", fontsize=8)
    ax.set_xlim(-0.5, 11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.invert_yaxis()

    out = "docs/reports/report1/figures/fig_gantt.png"
    with TufteCritic.disabled():
        fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
