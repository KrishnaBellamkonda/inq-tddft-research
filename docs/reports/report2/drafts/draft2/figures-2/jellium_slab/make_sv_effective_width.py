"""S(v) relabelled by effective width — figures-2 version (one-column, 3.5 × 3.0 in).

Legend: 2-column inside the axes (Classical | Wavepacket), header row shows
hollow=classical convention; labels carry σ_WP only, no "WP"/"classical" prefix.
Error bars: ±5 % of S at the highest-energy point (last velocity) per trace.
Legacy classical σ=0.5 capped at E ≤ 250 eV (v=4.5 / 275 eV point removed).

Reuses B.snapshot() from the source module for data; writes its own draw()
so the figure size, axes rect, and labels are controlled independently.

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_sv_effective_width.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
SRC  = REPO / "ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"
PANEL_DIR = REPO / "docs/reports/report2/drafts/draft1/figures"

sys.path.insert(0, str(SRC))
sys.path.insert(0, str(PANEL_DIR))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import matplotlib as mpl
from inqview.visualisation import style
style.apply_theme()
# Bump every text size +1 pt relative to the theme defaults.
for _key in ("font.size", "axes.labelsize", "axes.titlesize",
             "xtick.labelsize", "ytick.labelsize", "legend.fontsize"):
    mpl.rcParams[_key] += 1

import build_sv_effective_width_s6 as B   # noqa: E402  (reads CSVs at import)

HERE = Path(__file__).parent
HA   = 27.211386

_ERR_FRAC  = 0.05     # ±5 % of S at the highest-energy point per trace
_GREY      = "#888888"
_E_MAX_EV  = 250.0    # cap legacy classical σ=0.5 (removes the v=4.5 / 275 eV point)


def _add_errorbar(ax, df, x_col, y_col, colour):
    """Overlay a ±5 % error bar at the last (highest-energy) point."""
    if df.empty:
        return
    x = 0.5 * df[x_col].iloc[-1] ** 2 * HA
    y = df[y_col].iloc[-1]
    ax.errorbar(x, y, yerr=_ERR_FRAC * abs(y), fmt="none",
                color=colour, capsize=2.5, elinewidth=0.8, capthick=0.8)


def _cl_label(sigma):
    return rf"$\sigma_\mathrm{{WP}}$={sigma:g}"


def _wp_label(sigma, eff):
    if eff is None:
        return ""
    return rf"($\langle\sigma_r\rangle$={eff:.1f})"


def draw(snap: dict) -> None:
    fig, ax = style.figure_one_col()

    cl_handles, wp_handles = [], []

    # Classical: hollow markers, dashed line
    for sigma in (0.5, 2.0, 3.0):
        d = snap["legacy_cl"][sigma]
        if d.empty:
            continue
        # Remove the >250 eV point for σ=0.5 (v=4.5 → 275 eV)
        d = d[0.5 * d["v"] ** 2 * HA <= _E_MAX_EV]
        h, = ax.plot(0.5 * d["v"]**2 * HA, d["S_B_Eabs"],
                     ls="--", lw=1.0, color=B.COLOUR[sigma], marker="o",
                     mfc="none", ms=4.5, label=_cl_label(sigma))
        _add_errorbar(ax, d, "v", "S_B_Eabs", B.COLOUR[sigma])
        cl_handles.append(h)

    for sigma in B.NEW_SIGMAS:
        d = snap["new_cl"][sigma]
        if d.empty:
            continue
        h, = ax.plot(0.5 * d["v"]**2 * HA, d["S_B"],
                     ls="--", lw=1.0, color=B.COLOUR[sigma], marker="s",
                     mfc="none", ms=4.5, label=_cl_label(sigma))
        _add_errorbar(ax, d, "v", "S_B", B.COLOUR[sigma])
        cl_handles.append(h)

    # Wavepacket: filled markers, solid line
    for sigma in (0.5, 2.0, 3.0):
        d = snap["legacy_wp"][sigma]
        if d.empty:
            continue
        h, = ax.plot(0.5 * d["v"]**2 * HA, d["S_deposit_corrected"],
                     ls="-", lw=1.0, color=B.COLOUR[sigma], marker="o",
                     ms=4.5, label=_wp_label(sigma, snap["eff_legacy"][sigma]))
        _add_errorbar(ax, d, "v", "S_deposit_corrected", B.COLOUR[sigma])
        wp_handles.append(h)

    for sigma in B.NEW_SIGMAS:
        d = snap["new_wp"][sigma]
        if d.empty:
            continue
        h, = ax.plot(0.5 * d["v"]**2 * HA, d["S_B"],
                     ls="-", lw=1.0, color=B.COLOUR[sigma], marker="s",
                     ms=4.5, label=_wp_label(sigma, snap["eff_new"][sigma]))
        _add_errorbar(ax, d, "v", "S_B", B.COLOUR[sigma])
        wp_handles.append(h)

    ax.set_xlabel(r"Projectile Energy $E$ (eV)")
    ax.set_ylabel(style.axis_label("stopping_power", symbol="$S$"))

    ax.set_ylim(top=1.50)

    # matplotlib legend fills column-major (top→bottom, then next column).
    # With ncol=2 and 12 handles, first 6 → left column, next 6 → right column.
    # Layout: [hdr_cl, cl_0…cl_6] | [hdr_wp, wp_0…wp_6] gives paired rows.
    hdr_cl = Line2D([], [], ls="--", lw=1.0, color=_GREY, marker="o",
                    mfc="none", ms=4.5, label="Classical")
    hdr_wp = Line2D([], [], ls="-",  lw=1.0, color=_GREY, marker="o",
                    ms=4.5, label="Wavepacket")
    handles = [hdr_cl] + cl_handles + [hdr_wp] + wp_handles

    ax.legend(handles, [h.get_label() for h in handles],
              fontsize=7.0, frameon=False, ncol=2, loc="upper right",
              handlelength=1.7, handletextpad=0.4, columnspacing=1.0)

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
