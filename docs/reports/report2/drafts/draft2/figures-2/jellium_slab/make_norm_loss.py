"""Bath norm lost to CAP vs sigma_WP — figures-2 version (6.142 × 3.00 in).

Slot: 1\\textwidth = 6.142 in.
Title-case labels: "Initial Wavepacket Width ...", "$N_\\mathrm{lost}$ (Electrons)".
Legend inside axes (ncol=4, upper right).

Reuses N.collect() from the source module; writes its own draw().

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_norm_loss.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
SRC  = REPO / "ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"

sys.path.insert(0, str(SRC))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.visualisation import style
style.apply_theme()

import make_norm_loss as N   # noqa: E402  (loads data helpers at import)

HERE = Path(__file__).parent

# Axes rect for 6.142 × 3.00 in: 0.09*6.142=0.55 in left; 0.15*3.0=0.45 in bottom
_RECT = (0.09, 0.15, 0.88, 0.81)


def draw(df) -> None:
    fig = plt.figure(figsize=(6.142, 3.0))
    ax  = fig.add_axes(_RECT)

    for v in N.VS:
        w = df[(df.half == "wp") & (df.v == v)].sort_values("sigma_WP")
        if not w.empty:
            ax.plot(w.sigma_WP, w.excess_norm_lost, "-o", ms=4.5, lw=1.0,
                    color=N.COLOUR[v], label=rf"WP, $v$={v:g}")
        c = df[(df.half == "classical") & (df.v == v)].sort_values("sigma_WP")
        if not c.empty:
            ax.plot(c.sigma_WP, c.excess_norm_lost, "--s", ms=4.5, lw=1.0,
                    mfc="none", color=N.COLOUR[v], label=rf"classical, $v$={v:g}")

    # Box-size change marker
    ax.axvspan(3.0, 5.0, color="0.5", alpha=0.12, lw=0)

    ax.relim(); ax.autoscale_view()
    _, ymax = ax.get_ylim()
    ax.set_ylim(top=ymax * 1.30)

    ax.text(4.0, ymax * 1.30 * 0.97, r"$L_z$: 85 $\to$ 105",
            ha="center", va="top", fontsize=6, color="0.35")

    ax.set_xlabel(r"Initial Wavepacket Width $\sigma_\mathrm{WP}$ (Bohr)")
    ax.set_ylabel(r"$N_\mathrm{lost}$ (Electrons)")

    ax.legend(fontsize=6.5, frameon=False, ncol=4, loc="upper right",
              columnspacing=1.0, handlelength=1.7, handletextpad=0.4)

    out = HERE / "slab_norm_loss_vs_sigma.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    df = N.collect()
    df = df.sort_values(["half", "sigma_WP", "v"]).reset_index(drop=True)
    draw(df)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
