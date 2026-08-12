"""S vs 1/L slab — figures-2, sigma=5 only panel (3.071 × 3.00 in).

Slot: 0.5\\textwidth minipage = 3.071 in.
Shows the sigma_WP = 5 data with linear fit lines and the 1/L extrapolation.
Output: S_of_invL_sigma5_line.png

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_lzb_report_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
LZB_HYP = REPO / "ResearchProject/systems/localised_jellium/hypotheses/lz_bulk_sweep"

sys.path.insert(0, str(LZB_HYP))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.visualisation import style
style.apply_theme()

import build_lzb_figures as BL  # noqa: E402  (imports lzb_stopping as L)
import lzb_stopping as L        # noqa: E402

HERE = Path(__file__).parent

V_COLOR = {2.0: "C0", 2.5: "C1", 3.0: "C2", 3.5: "C3"}

# Axes rect for 3.071 × 3.00 in: 0.19*3.071=0.58 in left
_RECT = (0.19, 0.16, 0.78, 0.80)


def draw_sigma5(t: pd.DataFrame, a: pd.DataFrame,
                fits: pd.DataFrame) -> None:
    sigma = 5.0
    ok = t[t.complete] if not t.empty else t

    fig = plt.figure(figsize=(3.071, 3.0))
    ax  = fig.add_axes(_RECT)

    for half, filled in (("classical", False), ("wp", True)):
        for v in L.VELOCITIES:
            xs, ys = [], []
            if not ok.empty:
                m = ok[(ok.sigma_wp == sigma) & (ok.half == half) & (ok.v == v)]
                xs += list(m.inv_L); ys += list(m.S_deposit_eV_per_Bohr)
            if not a.empty:
                m = a[(a.sigma_wp == sigma) & (a.half == half) & (a.v == v)]
                xs += list(m.inv_L); ys += list(m.S)
            if not xs:
                continue
            col = V_COLOR[v]
            ax.plot(xs, ys,
                    "o" if filled else "s",
                    color=col, mfc=col if filled else "none", ms=5,
                    label=f"{'WP' if filled else 'cl.'} v={v}")
            if not fits.empty:
                f = fits[(fits.sigma_wp == sigma) & (fits.half == half)
                         & (fits.v == v)]
                if len(f):
                    c, sb = float(f.c_slope.iloc[0]), float(f.S_bulk.iloc[0])
                    xx = np.linspace(0.0, 1.0 / 15.0, 50)
                    ax.plot(xx, sb + c * xx,
                            "-" if filled else "--",
                            color=col, lw=0.8, alpha=0.6)
                    ax.plot([0.0], [sb],
                            "*" if filled else "P",
                            color=col,
                            mfc=col if filled else "none", ms=7)

    for lslab in (15.0, 25.0, 35.0):
        ax.axvline(1.0 / lslab, color="0.85", lw=0.5, zorder=0)

    ax.set_xlim(-0.003, 1.0 / 15.0 + 0.004)
    ax.set_xlabel(r"$1/L_{\mathrm{slab}}$ (Bohr$^{-1}$)")
    ax.set_ylabel(r"$S$ (eV/Bohr)")
    ax.annotate(rf"$\sigma_{{WP}} = {sigma:g}$",
                xy=(0.03, 0.96), xycoords="axes fraction", va="top", fontsize=8)
    ax.legend(fontsize=6, ncol=2, loc="upper right", frameon=False)

    out = HERE / "S_of_invL_sigma5_line.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    t, a = BL.assemble()
    if t.empty and a.empty:
        print("no data yet")
        return 1
    fits = BL.fit_rows(t, a)
    draw_sigma5(t, a, fits)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
