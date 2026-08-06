"""
The sigma56_sv deliverable: S(v) with every classical and wavepacket run on it.

Writes, into this directory:
    s56_S_summary.csv          every point in the sweep, with its evidence columns
    s56_cap_cost.csv           CAP-on minus CAP-free on the classical half (v = 3.0)
    S_of_v_sigma56.png         S vs v -- the main figure
    S_of_sigma_eq.png          S vs TIME-AVERAGED sigma -- the collapse test

--------------------------------------------------------------------------------
DESIGN PROVENANCE -- READ THIS BEFORE CALLING THE FIGURE FINAL
--------------------------------------------------------------------------------
The user asked for the same design as
    hypotheses/classical_highdensity_sv/dyn_direct/S_of_v_v2_timeavg_sigmar.png
That file, and the script that produced it, are NOT on this machine (the
directory holds only S_of_v_direct.csv and two notebook builders; the classical
results tree is empty -- the raw data went with the decommissioned
/local/data/public machine). What is written here therefore follows the PROJECT
standard (canonical theme, .claude/rules/report-figures + number-rounding) and the
one thing the reference's filename states unambiguously: points are placed on a
TIME-AVERAGED sigma axis. Reconcile with the reference when it is transferred.

--------------------------------------------------------------------------------
WHY TWO FIGURES
--------------------------------------------------------------------------------
S vs v is the physics figure: four twin pairs (sigma = 5 and 6, classical and
wavepacket) plus the existing sigma = 0.5/2/3 wavepacket traces for context.

S vs sigma_eq is the TEST. A wavepacket has no single width -- sigma_d(t) grows --
so a classical twin at fixed sigma_pot is only a fair comparison if the packet's
label agrees with its time-average. sigma_eq = sqrt(2)*<sigma_d> over the in-slab
transit is that time-average expressed as a sigma_WP label. On this axis:
  * the new sigma = 5/6 points barely move from their labels (5.3-5.7, 6.2-6.5);
  * the old sigma = 2 points spread across 4.0-6.4 depending on velocity;
  * a sigma = 6 point at v = 2.0 lands at sigma_eq = 6.45, essentially ON TOP of
    the existing sigma = 2, v = 2.0 point (sigma_eq = 6.35).
That last coincidence is a free, decisive check: if a packet that IS 4.5 Bohr wide
throughout gives the same S as one that AVERAGES 4.5 Bohr while sweeping 2.5->6.6,
then time-averaged sigma is a valid collapse variable. If not, it is not.

--------------------------------------------------------------------------------
HETEROGENEITY THAT MUST REACH THE CAPTION (user decision 2026-08-02)
--------------------------------------------------------------------------------
  * sigma = 5/6 ran at L_z = 105 with launch z = -27.5; sigma = 0.5/2/3 at
    L_z = 85, launch z = -24. The slab, r_s, dx, dt and CAP are identical.
  * sigma = 5/6 have CAP-on classical twins. sigma = 0.5/2/3 have none: the only
    classical reference there is the CAP-FREE sigma = 0.5 benchmark, which is a
    DIFFERENT estimator (the medium's gain directly, not its retained excitation).
    It is drawn dashed and grey for exactly that reason.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import s56_stopping as S
from inqview.visualisation import style

HERE = Path(__file__).resolve().parent
LJ_HYP = HERE.parent
HD = LJ_HYP / "wp_highdensity_sv"
CL = LJ_HYP / "classical_highdensity_sv"

# sigma_WP -> (colour, marker) for the NEW twin pairs; the legacy traces are the
# cooler end of the same family so the new work reads as the foreground.
NEW_STYLE = {5.0: ("#d62728", "D"), 6.0: ("#9467bd", "^")}
OLD_STYLE = {0.5: ("#1f77b4", "o"), 2.0: ("#ff7f0e", "s"), 3.0: ("#2ca02c", "v")}


def legacy_wp() -> pd.DataFrame:
    """The completed sigma = 0.5/2/3 wavepacket deposit points (L_z = 85).

    Read from sigma_sweep_S_deposit.csv, which carries `complete` and the same
    S_deposit_corrected column this campaign's estimator produces.
    """
    f = HD / "sigma_sweep_S_deposit.csv"
    if not f.exists():
        print(f"  NOTE: no legacy sweep at {f} — plotting the new twins only")
        return pd.DataFrame()
    d = pd.read_csv(f)
    d = d[d["complete"]].copy()
    # sigma_eq on the LEGACY geometry: launch z = -24, not -27.5. Computing it
    # with this campaign's launch would misplace every legacy point.
    d["sigma_eq"] = [
        float(np.sqrt(2) * _mean_sd(r.sigma, r.v, launch_z=-24.0))
        for r in d.itertuples()
    ]
    return d.rename(columns={"S_deposit_corrected": "S_eV_per_Bohr",
                             "sigma": "sigma_wp"})


def _mean_sd(sigma: float, v: float, launch_z: float, n: int = 4001) -> float:
    ti = (abs(launch_z) - S.SLAB_HALF) / v
    to = (abs(launch_z) + S.SLAB_HALF) / v
    t = np.linspace(ti, to, n)
    return float(np.trapezoid(S.sigma_d(t, sigma), t) / (to - ti))


def legacy_classical() -> pd.DataFrame:
    """The sigma = 0.5 CAP-FREE classical benchmark (6 velocities).

    A DIFFERENT estimator from everything else on the plot — the projectile was an
    external perturbation that never entered the electronic ledger, and with no
    absorber `plateau - E_GS` is the slab's gain directly rather than its retained
    excitation. Drawn dashed grey and captioned as such; it is orientation, not a
    twin.
    """
    f = CL / "sv_sweep" / "S_summary.csv"
    if not f.exists():
        print(f"  NOTE: no classical benchmark at {f}")
        return pd.DataFrame()
    return pd.read_csv(f)


def main() -> int:
    style.apply_theme()

    # ---- data ------------------------------------------------------------
    try:
        e_gs = S.e_gs_ha()
        print(f"E_GS (L_z = 105, dx = 0.40) = {e_gs:.9f} Ha")
    except FileNotFoundError as exc:
        print(f"FATAL: {exc}")
        return 2

    new = S.table()
    if new.empty:
        print("FATAL: no sigma56_sv runs have produced observables yet.")
        return 2
    new.to_csv(HERE / "s56_S_summary.csv", index=False)
    print(f"wrote {HERE/'s56_S_summary.csv'}  ({len(new)} rows)")

    cost = S.cap_cost()
    if not cost.empty:
        cost.to_csv(HERE / "s56_cap_cost.csv", index=False)
        print(f"wrote {HERE/'s56_cap_cost.csv'}")
        print(cost.to_string(index=False))

    # NEVER plot an incomplete point. A still-propagating run returns a perfectly
    # plausible S; the only defence is this filter. Dropped points are listed.
    dropped = new[~new["complete"]]
    for r in dropped.itertuples():
        print(f"  EXCLUDED (incomplete) {r.run}: {r.steps_done}/{r.steps_target} steps")
    unsettled = new[new["complete"] & ~new["settled"]]
    for r in unsettled.itertuples():
        print(f"  WARN (plateau still drifting) {r.run}: "
              f"drift {r.plateau_drift_eV:.3f} eV on E_abs {r.E_absorbed_eV:.2f} eV")
    ok = new[new["complete"]].copy()

    # PLOT THE MEASURED DEPOSIT AS-IS: S = [E_total(t_f) - E_GS]/L_slab.
    # classical -> raw E_absorbed/25 (the published S_B_Eabs convention)
    # WP        -> the norm-corrected deposit (published S_deposit_corrected)
    #
    # No further correction is applied (user decision, 2026-08-03). An earlier
    # revision plotted `S_deposit_eV_per_Bohr`, which subtracts an N_e/z monopole
    # term from the classical half; that changes measured values and must not be
    # the default. The column is still computed and lives in s56_S_summary.csv for
    # anyone who wants to look at it, but it is NOT what the figures show.
    prod = ok[ok["cap"]]                       # controls are not plotted

    old_wp = legacy_wp()
    old_cl = legacy_classical()

    # ---- figure 1: S vs v ------------------------------------------------
    fig, ax = style.figure_two_col(height_in=3.4)

    for sigma, (c, m) in NEW_STYLE.items():
        for half, ls, mfc in (("wp", "-", c), ("classical", "--", "none")):
            d = prod[(prod.sigma_wp == sigma) & (prod.half == half)].sort_values("v")
            if d.empty:
                continue
            ax.plot(d["v"], d["S_eV_per_Bohr"], ls=ls, marker=m, color=c,
                    mfc=mfc, ms=5.0, lw=1.4,
                    label=rf"$\sigma_{{\mathrm{{WP}}}}={sigma:g}$, "
                          + ("WP" if half == "wp" else "classical"))

    if not old_wp.empty:
        for sigma, (c, m) in OLD_STYLE.items():
            d = old_wp[old_wp.sigma_wp == sigma].sort_values("v")
            if d.empty:
                continue
            ax.plot(d["v"], d["S_eV_per_Bohr"], ls=":", marker=m, color=c,
                    ms=3.5, lw=1.0, alpha=0.75,
                    label=rf"$\sigma_{{\mathrm{{WP}}}}={sigma:g}$, WP ($L_z$=85)")

    if not old_cl.empty:
        ax.plot(old_cl["v"], old_cl["S_eV_per_Bohr"], ls="--", marker="x",
                color="0.45", ms=4.0, lw=1.0,
                label=r"$\sigma_{\mathrm{WP}}=0.5$, classical (CAP-free)")

    ax.set_xlabel("$v$ (a.u.)")
    ax.set_ylabel(style.axis_label("stopping_power", "$S$"))
    ax.set_xticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
    ax.legend(frameon=False, ncol=2, fontsize="small", handlelength=1.8,
              labelspacing=0.25, borderaxespad=0.2)
    out1 = HERE / "S_of_v_sigma56.png"
    fig.savefig(out1, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out1}")

    # ---- figure 2: S vs time-averaged sigma ------------------------------
    # The collapse test. Each point sits at the width the packet ACTUALLY had,
    # averaged over the transit, rather than at its launch label.
    fig, ax = style.figure_two_col(height_in=3.4)

    for sigma, (c, m) in NEW_STYLE.items():
        for half, mfc in (("wp", c), ("classical", "none")):
            d = prod[(prod.sigma_wp == sigma) & (prod.half == half)].sort_values("v")
            if d.empty:
                continue
            ax.plot(d["sigma_eq"], d["S_eV_per_Bohr"], ls="none", marker=m,
                    color=c, mfc=mfc, ms=5.5,
                    label=rf"$\sigma_{{\mathrm{{WP}}}}={sigma:g}$, "
                          + ("WP" if half == "wp" else "classical"))

    if not old_wp.empty:
        for sigma, (c, m) in OLD_STYLE.items():
            d = old_wp[old_wp.sigma_wp == sigma].sort_values("v")
            if d.empty:
                continue
            ax.plot(d["sigma_eq"], d["S_eV_per_Bohr"], ls="none", marker=m,
                    color=c, ms=4.0, alpha=0.75,
                    label=rf"$\sigma_{{\mathrm{{WP}}}}={sigma:g}$, WP ($L_z$=85)")

    # Mark the equality line: where a point's time-average equals its label, the
    # classical twin at that label is a fair comparison. Everything to the right
    # of its own label is a packet that outgrew its name.
    for sigma, (c, _) in {**OLD_STYLE, **NEW_STYLE}.items():
        ax.axvline(sigma, color=c, lw=0.6, ls="-", alpha=0.30)

    ax.set_xlabel(r"time-averaged width over the transit,  "
                  r"$\sqrt{2}\,\langle\sigma_d\rangle$ (Bohr)")
    ax.set_ylabel(style.axis_label("stopping_power", "$S$"))
    ax.legend(frameon=False, ncol=2, fontsize="small", handlelength=1.0,
              labelspacing=0.25, borderaxespad=0.2)
    out2 = HERE / "S_of_sigma_eq.png"
    fig.savefig(out2, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out2}")

    # ---- headline table (2 s.f., .claude/rules/number-rounding.md) --------
    print("\nS (eV/Bohr), CAP-on production points:")
    piv = prod.pivot_table(index=["sigma_wp", "half"], columns="v",
                           values="S_eV_per_Bohr")
    print(piv.round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
