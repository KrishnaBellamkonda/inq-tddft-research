#!/usr/bin/env python
"""S(v) relabelled by effective width — the upstream figure, with sigma_WP = 6 added.

Reproduces the design of
`docs/misc/misc-figures/S_of_v_effective_width.png`
(built by wp_highdensity_sv/build_momentum_width_notebook.py, merged 2026-08-03)
and adds this campaign's sigma_WP = 6 twin, plus sigma_WP = 5 as its runs land.

DESIGN, FROM UPSTREAM
  classical : hollow marker, dashed
  WP        : filled marker, solid
  one colour per twin pair; S_B on y, v on x; legend ncol=2 upper right.

LEGEND CONVENTION (user decision, 2026-08-05 -- DEVIATES from upstream)
  BOTH halves of a twin pair are labelled by the SAME sigma_WP, per
  `.claude/rules/sigma-wp-convention.md`: a classical run matched to a WP run is
  reported at its sigma_WP, never at the derived sigma_pot = sigma_WP/sqrt(2)
  that only exists inside the UPF generator. Upstream labelled the classical
  curves "sigma_pot = 3.54" while their WP twins said "5", which put the two
  halves of one pair on axes differing by sqrt(2) and made the pairing
  unreadable. The WP entries additionally carry <sigma_r> in brackets, which is
  what this figure exists to show:

      classical  sigma_WP = 5
      WP         sigma_WP = 5  (<sigma_r> = 7.8)

  METHODS FOOTNOTE (the only place sigma_pot belongs): the classical projectile
  is a Gaussian pseudopotential of width sigma_pot = sigma_WP/sqrt(2), chosen so
  its charge standard deviation equals that of the |psi|^2 of the WP twin. The
  sigma = 5 / 6 pairs therefore use sigma_pot = 3.54 / 4.24 Bohr.

<sigma_r> DEFINITION (user decision, 2026-08-03):
    sigma_r(t) = sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2)      (3-D radial spread)
    <sigma_r>  = time average from t=0 until the norm has fallen by 1 %

This REPLACES the earlier full-run average, which ran far past the point where the
CAP destroys the packet and so measured a smeared remnant (sigma_r spiking to
40-48 Bohr, then oscillating on a residue of norm ~1e-9). Implementation is
s56_stopping.sigma_r_window, shared by every consumer so one definition applies to
the legacy L_z=85 sweeps and this campaign alike.
It is NOT the `sigma_eq` used elsewhere in this campaign (1-D, transit window,
scaled by sqrt(2)). The two answer different questions and their numbers differ by
~4x, so they must never be mixed on one axis. Computed here the upstream way so
the sigma = 6 point is commensurate with the sigma = 0.5/2/3 points beside it.

ONE FIGURE, ONE ESTIMATOR: S_B = [E_tot(t_f) - E_GS]/L_slab, exactly as measured.

  classical -> E_absorbed/25, the published S_B_Eabs column, UNMODIFIED
  WP        -> the norm-corrected deposit, the published S_deposit_corrected

NO further correction is applied to either half (user decision, 2026-08-03). An
earlier version of this script also wrote a variant with an N_e/z monopole term
subtracted from the classical series; that altered measured values and has been
removed. Do not reintroduce it here.
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
import s56_stopping as S                                      # noqa: E402
sys.path.insert(0, str(HERE.parents[4]
                       / "docs/reports/report2/drafts/draft1/figures"))
from _panel import panel_mode as _panel_mode                  # noqa: E402
from _panel import slot_figure as _slot_figure                # noqa: E402
from inqview.visualisation import style                       # noqa: E402

style.apply_theme()

REPO = HERE.parents[4]
HYP = REPO / "ResearchProject/systems/localised_jellium/hypotheses"

# The report draft consumes this figure. Written on every build so the draft can
# never silently hold a stale copy (the runs are live; this figure changed twice
# on 2026-08-05 alone).
REPORT_FIG = (REPO / "docs/reports/report2/drafts/draft1/figures/jellium_slab"
              / "slab_sv_effective_width_s56.png")
WPH = HYP / "wp_highdensity_sv"
CL = HYP / "classical_highdensity_sv/dyn_direct"

HA = 27.211386
VS = [2.0, 2.5, 3.0, 3.5]
VTAG = {2.0: "v2p0", 2.5: "v2p5", 3.0: "v3p0", 3.5: "v3p5"}

# Upstream colours kept identical so the two figures overlay cleanly.
COLOUR = {0.5: "#1f5fb4", 2.0: "#d1600a", 3.0: "#2ca02c",
          5.0: "#c0392b", 6.0: "#7b3ba8"}
LEGACY_WP_SUMMARY = {0.5: WPH / "wp_S_summary.csv",
                     2.0: WPH / "wp_S_summary_s2p0.csv",
                     3.0: WPH / "wp_S_summary_s3p0.csv"}
LEGACY_WP_PREFIX = {0.5: "", 2.0: "s2p0_", 3.0: "s3p0_"}

# Which of THIS campaign's widths to overlay.
NEW_SIGMAS = (5.0, 6.0)

# Velocities EXCLUDED from the <sigma_r> MAX, per sigma_WP (user, 2026-08-05).
# The label is the largest <sigma_r> across the velocity runs of a set, so one
# suspect run silently sets the number for the whole curve -- hence the exclusion
# is a table here rather than a hand-edited constant.
#   sigma_WP = 6, v = 2.0 : user flagged its <sigma_r> as not to be used.
#     Its value (8.45) was the set MAX, so it alone defined the sigma = 6 label;
#     dropping it hands the label to v = 2.5 (8.15).
# Per-velocity values are printed by main() so any future change to this table is
# checkable against the numbers it moves.
EFF_EXCLUDE_V: dict[float, set[float]] = {6.0: {2.0}}

# sigma_WP -> {v: <sigma_r>} actually considered for the MAX; filled as the
# widths are computed and printed by main().
EFF_DETAIL: dict[float, dict[float, float]] = {}


def _load_concat(run_dir: Path, base: str) -> pd.DataFrame:
    """base.csv + every segment-suffixed base.fromNNNN.csv, in step order.

    MUST concatenate. Reading only the base file silently truncates any RESUMED
    run at its first kill point -- e.g. s3p0_v2p0 carries .from2172 and .from2896
    -- and because the packet is widest LATE, dropping the tail biases <sigma_r>
    LOW. That is exactly how this script first labelled the sigma_WP = 3 curve
    <sigma_r> ~ 16 when the published value is 17. Mirrors upstream's load_concat
    in build_momentum_width_notebook.py.
    """
    files = sorted(list(run_dir.glob(f"{base}.csv"))
                   + list(run_dir.glob(f"{base}.from*.csv")))
    if not files:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(f, comment="#") for f in files])
    return df.drop_duplicates("step").sort_values("step").reset_index(drop=True)


def mean_sigma_r_legacy(sigma: float) -> float | None:
    """LARGEST <sigma_r> across the velocity runs of a legacy (L_z = 85) set.

    MAX, not mean (user decision, 2026-08-03): the label should quote the widest
    the packet ever effectively is, which is invariably the SLOWEST run -- it
    spends longest in flight before the 1 % norm-loss cutoff, so it disperses
    most. A mean across velocities would understate the object the slab sees at
    the low-v end of the very curve being plotted.

    Runs that never reach the 1 % threshold are skipped: their window is the whole
    (too-short) run, so their <sigma_r> is not comparable.
    """
    vals = {}
    for v in VS:
        if v in EFF_EXCLUDE_V.get(sigma, set()):
            continue
        d = WPH / "sweep_data" / f"{LEGACY_WP_PREFIX[sigma]}{VTAG[v]}"
        rs = _load_concat(d, "wp_real_space_stats")
        if rs.empty or not {"sigma_x2", "sigma_y2", "sigma_z2"} <= set(rs.columns):
            continue
        m, _t, _n, reached = S.sigma_r_window(rs)
        if reached:
            vals[v] = m
    EFF_DETAIL[sigma] = vals
    return float(np.max(list(vals.values()))) if vals else None


def mean_sigma_r_new(sigma: float) -> float | None:
    """LARGEST <sigma_r> across this campaign's velocity runs (see the legacy twin).

    INCOMPLETE runs still count. The 1 % window closes at t ~ 17-28 a.u., far
    inside even a 45 %-finished run, so the width is fully determined long before
    the run reaches its step target -- unlike S, which needs the whole run.
    """
    vals = {}
    for v in VS:
        if v in EFF_EXCLUDE_V.get(sigma, set()):
            continue
        try:
            rs = S._concat(S.run_dir(sigma, v, "wp") / "raw" / "observables",
                           "wp_real_space_stats")
        except Exception:
            continue
        if not {"sigma_x2", "sigma_y2", "sigma_z2"} <= set(rs.columns):
            continue
        m, _t, _n, reached = S.sigma_r_window(rs)
        if reached:
            vals[v] = m
    EFF_DETAIL[sigma] = vals
    return float(np.max(list(vals.values()))) if vals else None


def legacy_classical(sigma: float) -> pd.DataFrame:
    if sigma == 0.5:
        d = pd.read_csv(CL / "S_of_v_cap.csv")
    else:
        d = pd.read_csv(CL / "S_of_v_cap_sigma.csv")
        d = d[np.isclose(d.sigma_WP, sigma)]
    return d.sort_values("v")


def legacy_wp(sigma: float) -> pd.DataFrame:
    f = LEGACY_WP_SUMMARY[sigma]
    if not f.exists():
        return pd.DataFrame()
    d = pd.read_csv(f)
    if "complete" in d:
        d = d[d["complete"]]
    return d.sort_values("v")


# Per-(sigma, half, v) audit of what this build actually plotted. The figure must
# carry EVERY completed sigma = 5/6 run (user, 2026-08-05) and no partial one, so
# the skip reason is recorded rather than swallowed -- a silently dropped run is
# indistinguishable from a run that was never launched.
COVERAGE: list[dict] = []


def new_series(sigma: float, half: str) -> pd.DataFrame:
    rows = []
    for v in VS:
        try:
            p = S.measure(sigma, v, half)
        except Exception as exc:
            COVERAGE.append({"sigma_WP": sigma, "half": half, "v": v,
                             "used": False,
                             "why": f"{type(exc).__name__}: {exc}"})
            continue
        if not p.complete:
            # Partial run: E_total(t_f) is a mid-transit excitation, NOT a
            # stopping value. Resume it (rt_state.txt / *_RESUME=1) rather than
            # plotting it.
            COVERAGE.append({"sigma_WP": sigma, "half": half, "v": v,
                             "used": False,
                             "why": f"incomplete {p.steps_done}/{p.steps_target}"
                                    f" steps ({100*p.steps_done/p.steps_target:.0f}%)"})
            continue
        COVERAGE.append({"sigma_WP": sigma, "half": half, "v": v, "used": True,
                         "why": f"complete {p.steps_done}/{p.steps_target}"})
        rows.append({"v": v, "S_B": p.S_eV_per_Bohr,
                     "S_corr": p.S_deposit_eV_per_Bohr})
    return pd.DataFrame(rows)


def snapshot() -> dict:
    """Read every series ONCE.

    The runs are live: on the first build, sigma = 5 v3.5 WP finished BETWEEN the
    raw and the corrected draw, so the two panels disagreed about which points
    existed. Two figures that are meant to differ only by a correction must be
    drawn from one read of the data, not two.
    """
    return {
        "legacy_cl": {s: legacy_classical(s) for s in (0.5, 2.0, 3.0)},
        "legacy_wp": {s: legacy_wp(s) for s in (0.5, 2.0, 3.0)},
        "new_cl": {s: new_series(s, "classical") for s in NEW_SIGMAS},
        "new_wp": {s: new_series(s, "wp") for s in NEW_SIGMAS},
        "eff_legacy": {s: mean_sigma_r_legacy(s) for s in (0.5, 2.0, 3.0)},
        "eff_new": {s: mean_sigma_r_new(s) for s in NEW_SIGMAS},
    }


def cl_label(sigma: float) -> str:
    """Classical entry -- sigma_WP, the SHARED twin label (never sigma_pot)."""
    return rf"classical $\sigma_\mathrm{{WP}}$={sigma:g}"


def wp_label(sigma: float, eff: float | None) -> str:
    """WP entry -- the same sigma_WP as its classical twin, plus <sigma_r>."""
    if eff is None:
        return rf"WP $\sigma_\mathrm{{WP}}$={sigma:g}"
    return (rf"WP $\sigma_\mathrm{{WP}}$={sigma:g} "
            rf"($\langle\sigma_r\rangle$={eff:.1f})")


def draw(corrected: bool, out: Path, snap: dict) -> None:
    # PANEL=1 re-authors at the FULL-width slot (6.142 x 2.63 in) rather than
    # scaling the 3.5 in standalone: `\includegraphics` scaling would render the
    # 10 pt labels at 17.6 pt. The extra width is also what lets the 10-entry
    # legend sit INSIDE the axes, which the house standard requires (no
    # bbox_inches="tight").
    if _panel_mode():
        fig, ax = _slot_figure("full")
    else:
        fig, ax = style.figure_one_col()

    # ---- classical: hollow, dashed ---------------------------------------
    for sigma in (0.5, 2.0, 3.0):
        d = snap["legacy_cl"][sigma]
        if d.empty:
            continue
        y = d["S_B_Eabs"].to_numpy()
        if corrected:
            # remove the N_e/z monopole still in the ledger at t_final
            y = (d["E_absorbed_eV"].to_numpy()
                 - 100.0 / d["z_final"].to_numpy() * HA) / 25.0
        ax.plot(d["v"], y, ls="--", lw=1.0, color=COLOUR[sigma], marker="o",
                mfc="none", ms=4.5, label=cl_label(sigma))

    for sigma in NEW_SIGMAS:
        d = snap["new_cl"][sigma]
        if d.empty:
            continue
        y = d["S_corr"] if corrected else d["S_B"]
        ax.plot(d["v"], y, ls="--", lw=1.0, color=COLOUR[sigma], marker="s",
                mfc="none", ms=4.5, label=cl_label(sigma))

    # ---- WP: filled, solid ------------------------------------------------
    for sigma in (0.5, 2.0, 3.0):
        d = snap["legacy_wp"][sigma]
        if d.empty:
            continue
        ax.plot(d["v"], d["S_deposit_corrected"], ls="-", lw=1.0,
                color=COLOUR[sigma], marker="o", ms=4.5,
                label=wp_label(sigma, snap["eff_legacy"][sigma]))

    for sigma in NEW_SIGMAS:
        d = snap["new_wp"][sigma]
        if d.empty:
            continue
        ax.plot(d["v"], d["S_B"], ls="-", lw=1.0, color=COLOUR[sigma],
                marker="s", ms=4.5,
                label=wp_label(sigma, snap["eff_new"][sigma]))

    ax.set_xlabel(r"projectile velocity $v$ (a.u.)")
    ax.set_ylabel(style.axis_label("stopping_power", symbol="$S_B$"))
    # Legend BELOW the axes, not upper-right as upstream: the sigma_WP entries
    # are ~2x wider than the old sigma_pot ones and the 2-column box then covers
    # the sigma_WP = 0.5 classical point at v = 2.0 (S_B = 1.3), the highest
    # point on the figure. `bbox_inches="tight"` at save grows the canvas to fit.
    if _panel_mode():
        ax.set_ylim(top=ax.get_ylim()[1] * 1.42)     # reserve room for the legend
        ax.legend(fontsize=6.5, frameon=False, ncol=3, loc="upper right",
                  columnspacing=1.0, handlelength=1.7, handletextpad=0.4)
    else:
        ax.legend(fontsize=7, frameon=False, ncol=2, loc="upper center",
                  bbox_to_anchor=(0.5, -0.17), columnspacing=1.1,
                  handlelength=1.9, handletextpad=0.5, borderaxespad=0.0)
    if not _panel_mode():
        # Report figures carry NO on-canvas title (house standard §1) -- it also
        # clips, since a panel figure may not use bbox_inches="tight".
        ax.set_title("E-absorbed stopping, WP relabelled by effective width"
                     + ("\n(classical monopole tail removed)" if corrected else ""),
                     fontsize=9)
    if _panel_mode():
        p = REPORT_FIG.parent / "slab_panel" / REPORT_FIG.name
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=600)                       # never 'tight' in a panel
        print(f"wrote {p}")
    else:
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"wrote {out}")
        if not corrected:
            REPORT_FIG.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(REPORT_FIG, dpi=300, bbox_inches="tight")
            print(f"wrote {REPORT_FIG}")
    plt.close(fig)


def main() -> int:
    snap = snapshot()          # ONE read; both figures are drawn from it

    print("sigma = 5/6 run coverage (every completed run is plotted):")
    used = sum(c["used"] for c in COVERAGE)
    for c in sorted(COVERAGE, key=lambda c: (c["sigma_WP"], c["half"], c["v"])):
        flag = "USED " if c["used"] else "SKIP "
        print(f"  {flag} sigma_WP={c['sigma_WP']:<4g} {c['half']:<9s}"
              f" v={c['v']:<4g} {c['why']}")
    print(f"  -> {used}/{len(COVERAGE)} sigma = 5/6 runs on the figure")

    print("\n<sigma_r> = MAX over velocities of the 1%-norm-loss-window mean:")
    for sigma in (0.5, 2.0, 3.0) + NEW_SIGMAS:
        e = (snap["eff_legacy"] if sigma < 5 else snap["eff_new"])[sigma]
        per_v = EFF_DETAIL.get(sigma, {})
        won = max(per_v, key=per_v.get) if per_v else None
        dropped = sorted(EFF_EXCLUDE_V.get(sigma, set()))
        line = (f"  sigma_WP={sigma:<4g} <sigma_r> = "
                + (f"{e:.2f}" if e else "n/a")
                + (f"  (from v={won:g})" if won is not None else "")
                + f"   [classical twin sigma_pot = {sigma/np.sqrt(2):.2f}]")
        if dropped:
            line += f"   EXCLUDED v={', '.join(f'{v:g}' for v in dropped)}"
        print(line)
        print("      per-v: "
              + "  ".join(f"v={v:g}:{m:.2f}" for v, m in sorted(per_v.items())))

    draw(False, HERE / "S_of_v_effective_width_s56.png", snap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
