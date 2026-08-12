"""S(v) relabelled by effective width — figures-2 version (6.142 × 2.63 in).

Slot: 1\\textwidth = 6.142 in.  Legend inside the axes (ncol=3, upper right).
Title-case label: "Projectile Energy $E$ (eV)".

Reuses B.snapshot() from the source module for data; writes its own draw()
so the figure size, axes rect, and labels are controlled independently.

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_sv_effective_width.py
"""
from __future__ import annotations

import os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
SRC  = REPO / "ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"
PANEL_DIR = REPO / "docs/reports/report2/drafts/draft1/figures"

sys.path.insert(0, str(SRC))
sys.path.insert(0, str(PANEL_DIR))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Apply theme BEFORE importing source module so its module-level apply_theme()
# picks up the updated rcParams (font.serif = Computer Modern Roman).
from inqview.visualisation import style
style.apply_theme()

import build_sv_effective_width_s6 as B   # noqa: E402  (reads CSVs at import)

HERE = Path(__file__).parent
HA   = 27.211386

# Axes rect for 6.142 × 2.63 in: 0.09*6.142=0.55 in left; 0.16*2.63=0.42 in bottom
_RECT = (0.09, 0.16, 0.88, 0.80)


def draw(snap: dict) -> None:
    fig = plt.figure(figsize=(6.142, 2.63))
    ax  = fig.add_axes(_RECT)

    # Classical: hollow, dashed
    for sigma in (0.5, 2.0, 3.0):
        d = snap["legacy_cl"][sigma]
        if d.empty:
            continue
        ax.plot(0.5 * d["v"]**2 * HA, d["S_B_Eabs"],
                ls="--", lw=1.0, color=B.COLOUR[sigma], marker="o",
                mfc="none", ms=4.5, label=B.cl_label(sigma))

    for sigma in B.NEW_SIGMAS:
        d = snap["new_cl"][sigma]
        if d.empty:
            continue
        ax.plot(0.5 * d["v"]**2 * HA, d["S_B"],
                ls="--", lw=1.0, color=B.COLOUR[sigma], marker="s",
                mfc="none", ms=4.5, label=B.cl_label(sigma))

    # WP: filled, solid
    for sigma in (0.5, 2.0, 3.0):
        d = snap["legacy_wp"][sigma]
        if d.empty:
            continue
        ax.plot(0.5 * d["v"]**2 * HA, d["S_deposit_corrected"],
                ls="-", lw=1.0, color=B.COLOUR[sigma], marker="o", ms=4.5,
                label=B.wp_label(sigma, snap["eff_legacy"][sigma]))

    for sigma in B.NEW_SIGMAS:
        d = snap["new_wp"][sigma]
        if d.empty:
            continue
        ax.plot(0.5 * d["v"]**2 * HA, d["S_B"],
                ls="-", lw=1.0, color=B.COLOUR[sigma], marker="s", ms=4.5,
                label=B.wp_label(sigma, snap["eff_new"][sigma]))

    ax.set_xlabel(r"Projectile Energy $E$ (eV)")
    ax.set_ylabel(style.axis_label("stopping_power", symbol="$S_B$"))

    # Reserve headroom above data for the inside legend.
    ax.relim(); ax.autoscale_view()
    _, ymax = ax.get_ylim()
    ax.set_ylim(top=ymax * 1.42)

    ax.legend(fontsize=6.5, frameon=False, ncol=3, loc="upper right",
              columnspacing=1.0, handlelength=1.7, handletextpad=0.4)

    out = HERE / "slab_sv_effective_width_s56.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    snap = B.snapshot()
    draw(snap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
