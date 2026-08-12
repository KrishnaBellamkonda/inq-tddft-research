#!/usr/bin/env python
"""Bath norm lost to the CAP — the "secondary emission" channel — vs sigma_WP.

THE QUANTITY (user request, 2026-08-05)

    excess_norm_lost = 100 - N_bath(t_final)
                     = 100 - [norm_total(t_f) - norm_wp(t_f)]

i.e. electrons removed from the SLAB by the absorbing potential, over and above
the projectile's own charge. Every WP run starts at norm_total = 101.0 exactly
(100 bath + 1 wavepacket); the wavepacket's own single electron is subtracted
off via norm_wp, so what is left is bath density that was excited hard enough to
reach the CAP and leave. That is the electron-emission (secondary-electron)
channel of the collision.

The classical twins are included as the CONTROL: same Gaussian charge, same
velocity, no wavepacket in the ledger at all. Their `norm_slab` column is the
same quantity with no subtraction needed, so the difference WP-minus-classical
isolates how much extra emission the quantum representation of the projectile
produces.

CROSS-CAMPAIGN CAVEAT — READ BEFORE COMPARING ACROSS THE WHOLE sigma AXIS
    sigma_WP = 0.5/2/3 ran at L_z = 85 (CAP inner edge at z = +-30);
    sigma_WP = 5/6    ran at L_z = 105 (CAP inner edge at z = +-40).
The legacy box holds its absorber 10 Bohr CLOSER to the slab on each face, so a
given excited electron reaches it sooner and with a lower escape threshold. Some
of any legacy-vs-new difference in this quantity is therefore geometric, not
physical. Within each campaign the comparison is clean.

Runs also stop at different times (t_f = 99.6-174.4 a.u. here), and emission
accumulates, so the table reports t_final alongside every value; a common-time
column is emitted too.

Writes, into this directory:
    norm_loss_table.csv        one row per run
    norm_loss_vs_sigma.png     excess norm lost vs sigma_WP, one line per velocity
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
sys.path.insert(0, str(REPO / "docs/reports/report2/drafts/draft1/figures"))
# Mirrored into the report draft on every build — see build_sv_effective_width_s6.py.
REPORT_FIG = (REPO / "docs/reports/report2/drafts/draft1/figures/jellium_slab"
              / "slab_norm_loss_vs_sigma.png")
sys.path.insert(0, str(HERE))
import s56_stopping as S                                        # noqa: E402
import build_sv_effective_width_s6 as B                         # noqa: E402
from _panel import panel_mode, slot_figure                      # noqa: E402
from inqview.visualisation import style                         # noqa: E402

style.apply_theme()

VS = [2.0, 2.5, 3.0, 3.5]
LEGACY = (0.5, 2.0, 3.0)
NEW = (5.0, 6.0)
COLOUR = {2.0: "#1f5fb4", 2.5: "#d1600a", 3.0: "#2ca02c", 3.5: "#7b3ba8"}
T_COMMON = 99.6          # a.u. — the shortest run in the set


def wp_frame(sigma: float, v: float) -> pd.DataFrame:
    if sigma in LEGACY:
        d = B.WPH / "sweep_data" / f"{B.LEGACY_WP_PREFIX[sigma]}{B.VTAG[v]}"
        return B._load_concat(d, "interactions")
    return S._concat(S.run_dir(sigma, v, "wp") / "raw" / "observables",
                     "interactions")


def cl_frame(sigma: float, v: float) -> pd.DataFrame:
    return S._concat(S.run_dir(sigma, v, "classical") / "raw" / "observables",
                     "interactions")


def collect() -> pd.DataFrame:
    rows = []
    for sigma in LEGACY + NEW:
        for v in VS:
            # ---- wavepacket half -------------------------------------------
            try:
                d = wp_frame(sigma, v)
            except Exception:
                d = pd.DataFrame()
            if not d.empty and {"norm_total", "norm_wp"} <= set(d.columns):
                bath = d.norm_total - d.norm_wp
                at_c = (100.0 - bath[(d.time_au - T_COMMON).abs().idxmin()]
                        if d.time_au.max() >= T_COMMON else np.nan)
                rows.append({"sigma_WP": sigma, "v": v, "half": "wp",
                             "box_Lz": 85 if sigma in LEGACY else 105,
                             "t_final_au": float(d.time_au.iloc[-1]),
                             "norm_wp_final": float(d.norm_wp.iloc[-1]),
                             "excess_norm_lost": float(100.0 - bath.iloc[-1]),
                             "excess_at_t99.6": float(at_c)})
            # ---- classical control (this campaign only) ---------------------
            if sigma in NEW:
                try:
                    c = cl_frame(sigma, v)
                except Exception:
                    continue
                if "norm_slab" not in c:
                    continue
                at_c = (100.0 - c.norm_slab[(c.time_au - T_COMMON).abs().idxmin()]
                        if c.time_au.max() >= T_COMMON else np.nan)
                rows.append({"sigma_WP": sigma, "v": v, "half": "classical",
                             "box_Lz": 105, "t_final_au": float(c.time_au.iloc[-1]),
                             "norm_wp_final": np.nan,
                             "excess_norm_lost": float(100.0 - c.norm_slab.iloc[-1]),
                             "excess_at_t99.6": float(at_c)})
    return pd.DataFrame(rows)


def draw(df: pd.DataFrame, out: Path) -> None:
    # PANEL=1 re-authors at the FULL-width slot so the 8-entry legend fits INSIDE
    # the axes; the house standard forbids bbox_inches="tight", so the standalone
    # figure's below-axes legend cannot survive in a panel.
    fig, ax = slot_figure("full") if panel_mode() else style.figure_one_col()
    for v in VS:
        w = df[(df.half == "wp") & (df.v == v)].sort_values("sigma_WP")
        if not w.empty:
            ax.plot(w.sigma_WP, w.excess_norm_lost, "-o", ms=4.5, lw=1.0,
                    color=COLOUR[v], label=rf"WP, $v$={v:g}")
        c = df[(df.half == "classical") & (df.v == v)].sort_values("sigma_WP")
        if not c.empty:
            ax.plot(c.sigma_WP, c.excess_norm_lost, "--s", ms=4.5, lw=1.0,
                    mfc="none", color=COLOUR[v], label=rf"classical, $v$={v:g}")
    # The box changes between sigma=3 and sigma=5 -- mark it, do not hide it.
    ax.axvspan(3.0, 5.0, color="0.5", alpha=0.12, lw=0)
    ax.text(4.0, ax.get_ylim()[1] * 0.97, r"$L_z$: 85 $\to$ 105",
            ha="center", va="top", fontsize=6, color="0.35")
    ax.set_xlabel(r"initial wavepacket width $\sigma_\mathrm{WP}$ (Bohr)")
    ax.set_ylabel("bath norm lost to CAP (electrons)")
    if panel_mode():
        # The full slot is only 2.63 in tall and the descriptive label overruns
        # it (the trailing "(electrons)" clipped on the first build). Panels get
        # the compact symbol; the word version lives in the caption.
        ax.set_ylabel(r"$N_\mathrm{lost}$ (electrons)")
        ax.set_ylim(top=ax.get_ylim()[1] * 1.30)     # room for the legend
        ax.legend(fontsize=6.5, frameon=False, ncol=4, loc="upper right",
                  columnspacing=1.0, handlelength=1.7, handletextpad=0.4)
        p = REPORT_FIG.parent / "slab_panel" / REPORT_FIG.name
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=600)
        plt.close(fig)
        print(f"wrote {p}")
        return
    ax.set_title("Secondary emission: slab charge absorbed,\n"
                 "excluding the projectile's own electron", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.17), columnspacing=1.1, handlelength=1.9)
    fig.savefig(out, dpi=600, bbox_inches=None)
    REPORT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(REPORT_FIG, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"wrote {REPORT_FIG}")


def main() -> int:
    df = collect()
    df = df.sort_values(["half", "sigma_WP", "v"]).reset_index(drop=True)
    df.to_csv(HERE / "norm_loss_table.csv", index=False)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nwrote {HERE/'norm_loss_table.csv'}")
    draw(df, HERE / "norm_loss_vs_sigma.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
