#!/usr/bin/env python
"""S(v) with every curve labelled by its TIME-AVERAGED width — sigma = 5/6 added.

This is the v3 successor to
`hypotheses/classical_highdensity_sv/dyn_direct/S_of_v_v2_timeavg_sigmar.png`
(that file and its builder live only on the user's other device; this script is a
reconstruction of the design from its description, not a copy of it).

WHY A TIME-AVERAGED LABEL AT ALL. A classical projectile has a FIXED width, so
its time-averaged sigma IS its sigma_WP. A wavepacket spreads, so over the slab
transit it behaves like a WIDER packet than its nominal label:

    sigma_d(t) = sqrt( sigma^2/2 + t^2/(2 sigma^2) )
    sigma_eq   = sqrt(2) * <sigma_d>       averaged over the in-slab window

`sigma_eq` is already expressed in the sigma_WP convention
(.claude/rules/sigma-wp-convention.md): it is the sigma_WP label a CONSTANT-width
packet would need to carry to present the same average width. So a WP curve
tagged sigma_bar = 6.2 is directly comparable with a classical curve at 6.2.

That is the whole point of this figure: a nominal-sigma plot silently compares a
spreading object with a rigid one, and the two are only the same object at t = 0.

STOPPING POWER. S = [E_total(t_f) - E_GS - E_PS(t_f)] / L_slab, L_slab = 25 Bohr.
The E_PS term is NOT optional on the classical half -- see s56_stopping.measure()
and the handover. Omitting it inflates classical S by ~4x and inverts the
WP/classical ordering.

Outputs (this directory):
    S_of_v_v3_timeavg_sigma.png
    s56_timeavg_sigma.csv        the per-point sigma_bar behind the labels
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
sys.path.insert(0, str(HERE))

import s56_stopping as S                                    # noqa: E402
import build_sv_figure as BF                                # noqa: E402
from inqview.visualisation import style                     # noqa: E402

# One colour per NOMINAL sigma so a twin pair shares a hue; WP solid + filled,
# classical dashed + hollow, exactly as in the v2 design.
COLOUR = {0.5: "#4C72B0", 2.0: "#DD8452", 3.0: "#55A868",
          5.0: "#C44E52", 6.0: "#8172B3"}
MARKER = {0.5: "o", 2.0: "s", 3.0: "v", 5.0: "D", 6.0: "^"}


def label(sigma: float, half: str, sbar: float) -> str:
    """Legend text. The time-averaged width leads because it is the quantity the
    two halves are being compared AT; the nominal sigma follows in parentheses so
    a run is still identifiable from the legend."""
    if half == "classical":
        # A rigid projectile does not spread: sigma_bar == sigma_WP by
        # construction. Saying so is more honest than printing a fitted number.
        return rf"$\bar\sigma={sigma:.2g}$ (classical, fixed)"
    return rf"$\bar\sigma={sbar:.3g}$ ($\sigma_{{\mathrm{{WP}}}}={sigma:g}$, WP)"


def main() -> int:
    new = S.table()
    new["S_eV_per_Bohr"] = new["S_deposit_eV_per_Bohr"]      # E_PS-corrected
    prod = new[new["complete"] & new["cap"]].copy()

    old_wp = BF.legacy_wp()
    old_cl = BF.legacy_classical()

    rows = []
    fig, ax = style.figure_two_col(height_in=3.6)

    # ---- new sigma = 5/6 twins ------------------------------------------
    for sigma in sorted(prod.sigma_wp.unique()):
        for half, ls, fill in (("wp", "-", True), ("classical", "--", False)):
            d = prod[(prod.sigma_wp == sigma) & (prod.half == half)].sort_values("v")
            if d.empty:
                continue
            sbar = float(d.sigma_eq.mean())
            c = COLOUR.get(sigma, "0.4")
            ax.plot(d["v"], d["S_eV_per_Bohr"], ls=ls, marker=MARKER.get(sigma, "o"),
                    color=c, mfc=c if fill else "none", ms=5.2, lw=1.5,
                    label=label(sigma, half, sbar))
            for r in d.itertuples():
                # sigma_eq is a DISPERSION quantity -- it only means something for
                # a spreading packet. The classical projectile is rigid, so its
                # time-averaged width is its sigma_WP identically. Carrying the
                # WP's sigma_eq onto the classical row would fabricate a spread
                # the run does not have.
                sbar_pt = r.sigma_eq if half == "wp" else sigma
                rows.append({"sigma_wp": sigma, "v": r.v, "half": half,
                             "sigma_eq": sbar_pt, "S_eV_per_Bohr": r.S_eV_per_Bohr,
                             "set": "sigma56_sv (L_z=105)"})

    # ---- legacy sigma = 0.5/2/3 wavepackets (L_z = 85) -------------------
    for sigma in sorted(old_wp.sigma_wp.unique()) if not old_wp.empty else []:
        d = old_wp[old_wp.sigma_wp == sigma].sort_values("v")
        if d.empty:
            continue
        c = COLOUR.get(sigma, "0.4")
        ax.plot(d["v"], d["S_eV_per_Bohr"], ls=":", marker=MARKER.get(sigma, "o"),
                color=c, ms=3.8, lw=1.1, alpha=0.8,
                label=label(sigma, "wp", float(d.sigma_eq.mean())) + r" $L_z$=85")
        for r in d.itertuples():
            rows.append({"sigma_wp": sigma, "v": r.v, "half": "wp",
                         "sigma_eq": r.sigma_eq, "S_eV_per_Bohr": r.S_eV_per_Bohr,
                         "set": "wp_highdensity_sv (L_z=85)"})

    # ---- legacy sigma = 0.5 classical -----------------------------------
    # NOTE: this curve was scored by PROJECTILE KE LOSS (S_of_v_direct.csv carries
    # v_final / v_mean_slab / deposit_eV), not by the field-side deposit used for
    # every other curve here. It is therefore not contaminated by the E_PS tail,
    # but it is also not the same measurement -- hence the distinct grey styling.
    if not old_cl.empty:
        ax.plot(old_cl["v"], old_cl["S_eV_per_Bohr"], ls="-.", marker="x",
                color="0.45", ms=4.2, lw=1.1,
                label=r"$\bar\sigma=0.5$ (classical, KE-loss)")

    ax.set_xlabel("$v$ (a.u.)")
    ax.set_ylabel(style.axis_label("stopping_power", "$S$"))
    ax.set_xticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
    ax.legend(frameon=False, ncol=2, fontsize="x-small", handlelength=2.0,
              labelspacing=0.28, borderaxespad=0.2)

    out = HERE / "S_of_v_v3_timeavg_sigma.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    t = pd.DataFrame(rows)
    t.to_csv(HERE / "s56_timeavg_sigma.csv", index=False)
    print(f"wrote {HERE/'s56_timeavg_sigma.csv'}  ({len(t)} rows)")

    print("\ntime-averaged width per curve (sigma_bar, sigma_WP convention):")
    for (sig, half), g in t.groupby(["sigma_wp", "half"]):
        print(f"  sigma_WP={sig:>4g} {half:<9} "
              f"sigma_bar = {g.sigma_eq.min():.3g}-{g.sigma_eq.max():.3g} "
              f"(mean {g.sigma_eq.mean():.3g})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
