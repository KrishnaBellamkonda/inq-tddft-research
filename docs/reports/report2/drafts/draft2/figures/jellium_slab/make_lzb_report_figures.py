"""Draft-2 S(1/L) Lz-sweep figures (fig 17 = sigma_WP=5 panel).

Runs family2() from the original script but writes output to this directory.
  Output: slab_S_of_invL_sigma5.png (600 DPI, bbox_inches=None)
  Also writes slab_S_of_invL_sigma0p5.png (not used in report but produced
  alongside for consistency).

Okabe-Ito per-velocity colour scheme — correct for this figure type.
No WP/classical binary: this figure shows WP-only S(1/L) extrapolation.

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_lzb_report_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
HERE = Path(__file__).parent
LZ_SRC = (REPO / "docs/reports/report2/drafts/draft1/figures"
          / "jellium_slab/Lz_sweep")

sys.path.insert(0, str(LZ_SRC))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import matplotlib
matplotlib.use("Agg")

# Import data + helpers from the original script but keep our own save path
import make_lzb_report_figures as M
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def save_here(fig, stem: str) -> None:
    """Save PNG to THIS directory at 600 DPI with bbox_inches=None."""
    p = HERE / f"{stem}.png"
    fig.savefig(p, dpi=600, bbox_inches=None)
    print(f"    {p}")


def family2_draft2(df: pd.DataFrame) -> None:
    """Reproduce family2 but save to draft-2 directory."""
    ok = df[df["complete"] & (df["half"] == "wp")].copy()
    ok["inv_L"] = 1.0 / ok["L_slab"]

    # Pre-compute fits (identical logic to M.family2)
    fit_cache: dict = {}
    print()
    for sigma in (0.5, 5.0):
        sub_s = ok[ok["sigma_wp"] == sigma]
        for v in M.V_LIST if hasattr(M, "V_LIST") else [2.0, 2.5, 3.0, 3.5]:
            rows = sub_s[sub_s["v_au"] == v]
            if rows.empty:
                continue
            x   = rows["inv_L"].values
            y   = rows["S_eV_per_Bohr"].values
            err = rows["S_err_eV_per_Bohr"].values
            finite = ~np.isnan(err)
            x_f, y_f, e_f = x[finite], y[finite], err[finite]
            if len(x_f) < 2:
                continue
            w = 1.0 / e_f ** 2
            a, b, a_err, b_err, chi2_red, dof = M.weighted_lstsq(x_f, y_f, w)
            fit_cache[(sigma, v)] = (a, b, a_err, b_err, chi2_red, dof)

    y_raw = list(ok["S_eV_per_Bohr"].values)
    y_int = [v[0] for v in fit_cache.values()]
    y_arr = np.asarray(y_raw + y_int)
    ymin, ymax = y_arr.min(), y_arr.max()
    yr = ymax - ymin
    ylim = [ymin - 0.06 * yr, ymax + 0.06 * yr]
    xlim = [-0.03 / 15.0, 1.10 / 15.0]

    V_LIST = [2.0, 2.5, 3.0, 3.5]

    for sigma in (0.5, 5.0):
        sub_s = ok[ok["sigma_wp"] == sigma]
        fig, ax = M._figure_one_col()

        ax.axvline(0.0, color="#AAAAAA", linewidth=0.6, zorder=0)
        leg_handles = []

        for v in V_LIST:
            col    = M.V_COLOR[v]
            marker = M.V_MARKER[v]
            rows   = sub_s[sub_s["v_au"] == v].sort_values("inv_L")
            if rows.empty:
                continue

            x   = rows["inv_L"].values
            y   = rows["S_eV_per_Bohr"].values
            err = rows["S_err_eV_per_Bohr"].values
            finite = ~np.isnan(err)

            if finite.any():
                ax.errorbar(x[finite], y[finite], yerr=err[finite],
                            fmt="none", ecolor=col, elinewidth=1, capsize=2, zorder=2)
            ax.plot(x, y, marker=marker, color=col,
                    mfc=col, mec=col, ms=4.5, ls="none", zorder=3)

            if (sigma, v) in fit_cache:
                S_bulk, c, S_bulk_err, c_err, chi2_red, dof = fit_cache[(sigma, v)]
                xx = np.linspace(0.0, 1.10 / 15.0, 100)
                ax.plot(xx, S_bulk + c * xx, ls="--", color=col, lw=0.9, zorder=1, alpha=0.85)
                ax.plot([0.0], [S_bulk], marker=marker, color=col,
                        mfc="white", mec=col, ms=5.5, ls="none", zorder=4)
                label = (f"$E = {M.E_LABEL[v]}$ eV:  "
                         f"$S_{{\\rm bulk}} = {M.compact_val_err(S_bulk, S_bulk_err)}$")
            else:
                label = f"$E = {M.E_LABEL[v]}$ eV: insuff. data"

            leg_handles.append(
                plt.Line2D([0], [0], marker=marker, color=col,
                           mfc=col, mec=col, ms=4.5, ls="--", lw=0.9, label=label)
            )

        if sigma == 0.5:
            ax.text(0.97, 0.97, r"$\sigma_{\rm WP} = 0.5$",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9)
            leg_loc = "upper left"
        else:
            ax.text(0.03, 0.03, r"$\sigma_{\rm WP} = 5$",
                    transform=ax.transAxes, ha="left", va="bottom", fontsize=9)
            leg_loc = "upper right"

        ax.legend(handles=leg_handles, loc=leg_loc, frameon=False, fontsize=6.5)
        ax.set_xlabel(r"$1/L_z$ (Bohr$^{-1}$)")
        ax.set_ylabel("$S$ (eV/Bohr)")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        sigma_str = "0p5" if sigma == 0.5 else "5"
        stem = f"slab_S_of_invL_sigma{sigma_str}"
        save_here(fig, stem)
        plt.close(fig)
        print(f"  sigma={sigma}: done")


def main() -> int:
    print("=== draft-2 make_lzb_report_figures ===")
    df = M.load_data()
    print("\n--- S(1/L), two figures ---")
    family2_draft2(df)
    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
