"""
Figures + summary CSVs for the slab->bulk L_slab sweep.

    python build_lzb_figures.py        # or build_lzb_figures.main() from finalize

Deliverables (all in this folder, and every figure ALSO written to the report-2
figure source with a `slab_` prefix — the feedback rule of 2026-08-05: wired
into the builder, never a manual copy):

    lzb_S_summary.csv     every new (box, v, half) point, evidence attached
    lzb_anchors.csv       the L = 25 anchor points actually used
    S_of_invL.png         S_deposit vs 1/L_slab per sigma family, velocity as
                          colour, classical open / WP filled, linear 1/L fits
                          with the S_bulk intercept marked at 1/L = 0
    lzb_fits.csv          per-(sigma, half, v) fit: slope c, intercept S_bulk

House standard for the report copy (report2 CLAUDE.md §1): no on-canvas title,
600 dpi, fixed axes rect (bbox_inches=None), legend inside the axes.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import lzb_stopping as L   # noqa: E402

REPORT_DIR = (L.REPO / "docs/reports/report2/drafts/draft1/figures/jellium_slab")

import matplotlib                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt         # noqa: E402

try:
    from inqview.visualisation import style as ivstyle
    ivstyle.apply_theme()
except Exception as _e:                                          # noqa: BLE001
    print(f"  (canonical theme unavailable: {_e}; default rcParams)")

V_COLOR = {2.0: "C0", 2.5: "C1", 3.0: "C2", 3.5: "C3"}


def save_both(fig, name: str, dpi: int = 600) -> None:
    """One render to the hypotheses folder AND the report-2 source folder."""
    local = HERE / name
    fig.savefig(local, dpi=dpi, bbox_inches=None)
    print(f"  wrote {local}")
    try:
        REPORT_DIR.mkdir(parents=True, exist_ok=True)
        rep = REPORT_DIR / f"slab_{name}"
        fig.savefig(rep, dpi=dpi, bbox_inches=None)
        print(f"  wrote {rep}")
    except Exception as e:                                       # noqa: BLE001
        print(f"  REPORT COPY FAILED ({type(e).__name__}: {e}) — local copy stands")


def assemble() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(points, anchors): the new measurements and the L = 25 anchors, both on
    the corrected-deposit estimator."""
    t = L.table()
    a = L.anchors()
    if not t.empty:
        t.to_csv(HERE / "lzb_S_summary.csv", index=False)
        print(f"  wrote lzb_S_summary.csv ({len(t)} rows)")
    if not a.empty:
        a.to_csv(HERE / "lzb_anchors.csv", index=False)
        print(f"  wrote lzb_anchors.csv ({len(a)} rows)")
    return t, a


def fit_rows(t: pd.DataFrame, a: pd.DataFrame) -> pd.DataFrame:
    """Linear fits S = S_bulk + c * (1/L) per (sigma, half, v).

    Requires >= 2 usable points; with all 3 the residual of the middle point is
    the LINEARITY check the 3-thickness design exists for. sigma = 0.5 WP rows
    are fitted too but flagged qualitative (the packet's in-slab width grows
    with L, so its slope mixes surface and width-growth terms — plan caveat).
    """
    rows = []
    ok = t[t.complete] if not t.empty else t
    for sigma in (0.5, 5.0):
        for half in ("wp", "classical"):
            for v in L.VELOCITIES:
                xs, ys, srcs = [], [], []
                if not ok.empty:
                    m = ok[(ok.sigma_wp == sigma) & (ok.half == half) & (ok.v == v)]
                    xs += list(m.inv_L); ys += list(m.S_deposit_eV_per_Bohr)
                    srcs += list(m.cfg)
                if not a.empty:
                    m = a[(a.sigma_wp == sigma) & (a.half == half) & (a.v == v)]
                    xs += list(m.inv_L); ys += list(m.S)
                    srcs += list(m.source)
                if len(xs) < 2:
                    continue
                x, y = np.asarray(xs), np.asarray(ys)
                c, s_bulk = np.polyfit(x, y, 1)
                resid = float(np.max(np.abs(np.polyval([c, s_bulk], x) - y))) \
                    if len(x) == 3 else np.nan
                rows.append({"sigma_wp": sigma, "half": half, "v": v,
                             "n_points": len(x), "S_bulk": s_bulk, "c_slope": c,
                             "max_resid_eV_per_Bohr": resid,
                             "qualitative": (sigma == 0.5 and half == "wp"),
                             "points": ";".join(map(str, srcs))})
    return pd.DataFrame(rows)


def figure(t: pd.DataFrame, a: pd.DataFrame, fits: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.14, top=0.92, wspace=0.25)
    ok = t[t.complete] if not t.empty else t

    for ax, sigma in zip(axes, (0.5, 5.0)):
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
                ax.plot(xs, ys, "o" if filled else "s", color=col,
                        mfc=col if filled else "none", ms=5,
                        label=f"{'WP' if filled else 'cl.'} v={v}")
                if not fits.empty:
                    f = fits[(fits.sigma_wp == sigma) & (fits.half == half)
                             & (fits.v == v)]
                    if len(f):
                        c, sb = float(f.c_slope.iloc[0]), float(f.S_bulk.iloc[0])
                        xx = np.linspace(0.0, 1.0 / 15.0, 50)
                        ax.plot(xx, sb + c * xx, "-" if filled else "--",
                                color=col, lw=0.8, alpha=0.6)
                        ax.plot([0.0], [sb], "*" if filled else "P", color=col,
                                mfc=col if filled else "none", ms=7)
        for lslab in (15.0, 25.0, 35.0):
            ax.axvline(1.0 / lslab, color="0.85", lw=0.5, zorder=0)
        ax.set_xlim(-0.003, 1.0 / 15.0 + 0.004)
        ax.set_xlabel(r"$1/L_{\mathrm{slab}}$ (Bohr$^{-1}$)")
        ax.set_ylabel(r"$S$ (eV/Bohr)")
        ax.annotate(rf"$\sigma_{{WP}} = {sigma:g}$" +
                    ("  (WP trace qualitative)" if sigma == 0.5 else ""),
                    xy=(0.03, 0.96), xycoords="axes fraction", va="top", fontsize=8)
        ax.legend(fontsize=6, ncol=2, loc="upper right", frameon=False)
    save_both(fig, "S_of_invL.png")
    plt.close(fig)


def main() -> int:
    print("=== build_lzb_figures ===")
    t, a = assemble()
    if t.empty and a.empty:
        print("  nothing measurable yet — no figure written")
        return 1
    fits = fit_rows(t, a)
    if not fits.empty:
        fits.to_csv(HERE / "lzb_fits.csv", index=False)
        print(f"  wrote lzb_fits.csv ({len(fits)} fits)")
        print(fits.round(3).to_string(index=False))
    figure(t, a, fits)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
