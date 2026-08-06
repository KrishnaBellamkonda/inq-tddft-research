#!/usr/bin/env python3
"""READABLE RESULTS for the cylindrical proximity ladder: S(T1) and S(T2) vs coupling.

Plan:     docs/plans/cylindrical-proximity-ladder.md
Handover: docs/handovers/cylindrical-proximity-ladder.md

WHAT THIS ADDS OVER `build_ladder_figures.py`
---------------------------------------------
`build_ladder_figures.py` is the exhaustive per-rung report (24 panels x 5 rungs)
driven by the channeling twin's drawing engine. It answers "what happened in
every run". This module answers the much narrower question the user actually
asked: *how does the stopping power, defined via T1 and via T2, trend as the
interaction is turned up?* — five figures, one table, no per-rung detail.

THE FIT WINDOW, AND WHY THE OLD ONE WAS WRONG (corrected 2026-08-03)
-------------------------------------------------------------------
The inherited windows were `T1 9-25`, `T2 21-30`, `T2 5-20` (channeling twin).
Scanning the LOCAL slope -dE/ds(t) at every rung shows those T2 windows are
unusable, and that this is a property of the window, not of T2:

  * S is ~0 at t=2 and peaks at t~8-11 in EVERY estimator and every rung. That
    is WAKE BUILD-UP: at r_s = 3, n = 3/(4 pi r_s^3) = 8.84e-3, so
    omega_p = sqrt(4 pi n) = 0.333 a.u. and a quarter plasma period is 4.7 a.u.
    The bath cannot polarise faster than that. No steady-state S exists before
    t ~ 10, so a window starting at t=5 averages over the transient.
  * S(T2) is NEGATIVE at early times for r10-r04 and crosses zero at
    t = 5.1 / 14.1 / 11.5 / 5.2 (r10/r08/r06/r04); at r00 it never goes negative.
    `T2 5-20` straddles that sign change at every rung and `T2 21-30` sits in the
    late decay, so both produced numbers that jumped around with no trend. The
    old "T2 is erratic" verdict was an artefact of these windows.
  * The upper edge is set by the light-projectile criterion
    (.claude/rules/light-projectile-stopping.md): the CLASSICAL half falls below
    0.85 v0 at t = 20.6 (r00) and t = 26.0 (r04); r10/r06/r08 never do.

So ONE window, t in [11, 20] a.u., is used for every estimator and every rung:
after wake build-up, before the velocity criterion fails anywhere, and clear of
the T2 zero crossing at all rungs but r08 (whose crossing at 14.1 is inside it —
flagged in the table by its low r^2, not hidden).

Applying one window to all three estimators is the point: T1, T2 and the
classical KE are then differenced over the SAME path interval at the SAME
velocity, so the ratios are comparable by construction.

THE COUPLING COORDINATE
-----------------------
x is the MEASURED occupancy `f_wall` — the fraction of |psi|^2 inside the
jellium — averaged over the fit window, with its min/max drawn as horizontal
error bars. Nominal R_in/sigma is a secondary readability axis only: the bore is
not empty (r10's ground state already has 16 of 160 electrons inside it), so
"distance to the background edge" is not "distance to the electrons".
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SYSTEM = HERE.parents[1]                     # .../systems/cylindrical_jellium
TWIN_HYP = SYSTEM / "hypotheses/channeling_twin"
for p in (str(HERE), str(TWIN_HYP)):
    if p not in sys.path:
        sys.path.insert(0, p)

import build_ladder_figures as L              # noqa: E402  rung -> results paths
import build_report_figures as BR             # noqa: E402  load_refined, local_slope

from inqview.visualisation import style       # noqa: E402

OUT = HERE / "figures/results"

# --- the one window, justified in the module docstring ---------------------
WIN = (11.0, 20.0)
V_FLOOR = 0.85          # light-projectile criterion, checked and reported

RUNGS = ["r10", "r08", "r06", "r04", "r00"]
RIN = {"r10": 10.0, "r08": 8.0, "r06": 6.0, "r04": 4.0, "r00": 0.0}
SIGMA_WP = 4.0
# Human labels: the ladder is a coupling ladder, so say the coupling, not the tag.
NICE = {"r10": r"$R_{\rm in}=2.5\,\sigma$", "r08": r"$R_{\rm in}=2.0\,\sigma$",
        "r06": r"$R_{\rm in}=1.5\,\sigma$", "r04": r"$R_{\rm in}=1.0\,\sigma$",
        "r00": "filled tube"}
#: bare tick-style label for annotating points, where the full mathtext is noise
SHORT = {"r10": r"2.5$\sigma$", "r08": r"2.0$\sigma$", "r06": r"1.5$\sigma$",
         "r04": r"1.0$\sigma$", "r00": "filled"}

# Rung colour = interaction strength, light -> dark, so "stronger coupling" is
# legible from the colour alone and never has to be looked up in the legend.
RUNG_COLOR = {r: c for r, c in zip(
    RUNGS, plt.get_cmap("viridis")(np.linspace(0.82, 0.06, len(RUNGS))))}

# Estimator identity is fixed across every figure in this module.
EST = {
    "T1": dict(color="tab:red",  label=r"$T_1=|\langle p\rangle|^2/2m$  (drift)"),
    "T2": dict(color="tab:green", label=r"$T_2=\langle p^2\rangle/2m$  (total)"),
    "CL": dict(color="tab:blue", label=r"classical $\frac{1}{2}mv^2$"),
}
COL = {"T1": "T1_drift_ev", "T2": "T2_total_ev"}


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def collect() -> tuple[dict, pd.DataFrame]:
    """Load every present rung; return per-rung frames and the summary table.

    Degrades like the ladder builder: a rung whose runs are missing or did not
    complete is REPORTED and skipped, never fatal.
    """
    data, rows = {}, []
    for tag, cfg in L.rung_run_sets(RUNGS).items():
        ok, why = L.rung_is_present(tag, cfg)
        if not ok:
            print(f"  [skip] {tag}: {why}")
            continue
        R = BR.load_refined(cfg["wp_results"], cfg["cl_results"])
        wp, cl = R.wp_frame(cfg["wp_name"]), R.cl_frame(cfg["cl_name"])

        d = dict(wp=wp, cl=cl,
                 s_wp=wp.s_pintegral.to_numpy() - wp.s_pintegral.iloc[0],
                 s_cl=cl.z_unwrapped.to_numpy() - cl.z_unwrapped.iloc[0])
        # local -dE/ds, centred +-1 a.u. stencil (same estimator as the report)
        d["S_T1"] = BR.local_slope(wp.T1_drift_ev.to_numpy(),
                                   wp.s_pintegral.to_numpy(), wp.t.to_numpy())
        d["S_T2"] = BR.local_slope(wp.T2_total_ev.to_numpy(),
                                   wp.s_pintegral.to_numpy(), wp.t.to_numpy())
        d["S_CL"] = BR.local_slope(cl.ke_ev.to_numpy(),
                                   cl.z_unwrapped.to_numpy(), cl.t.to_numpy())
        data[tag] = d

        # --- window fits: identical window for all three estimators ---------
        f = {e: R.fit_in_window(wp.s_pintegral.to_numpy(), wp[COL[e]].to_numpy(),
                                wp.t.to_numpy(), *WIN) for e in ("T1", "T2")}
        f["CL"] = R.fit_in_window(cl.z_unwrapped.to_numpy(), cl.ke_ev.to_numpy(),
                                  cl.t.to_numpy(), *WIN)

        mw = (wp.t >= WIN[0]) & (wp.t <= WIN[1])
        mc = (cl.t >= WIN[0]) & (cl.t <= WIN[1])
        fw = wp.f_wall[mw]
        row = {
            "rung": tag, "R_in": RIN[tag], "R_in_over_sigma": RIN[tag] / SIGMA_WP,
            "f_wall_mean": float(fw.mean()), "f_wall_lo": float(fw.min()),
            "f_wall_hi": float(fw.max()),
            # velocity criterion, reported not assumed
            "v_over_v0_min_cl": float(np.sqrt(cl.ke_ev[mc] / cl.ke_ev.iloc[0]).min()),
            "v_over_v0_min_wp": float(np.sqrt(wp.T1_drift_ev[mw]
                                              / wp.T1_drift_ev.iloc[0]).min()),
            # whole-run losses, for the cumulative-loss figure's annotation
            "dE_T1_ev": float(wp.T1_drift_ev.iloc[0] - wp.T1_drift_ev.iloc[-1]),
            "dE_T2_ev": float(wp.T2_total_ev.iloc[0] - wp.T2_total_ev.iloc[-1]),
            "dE_CL_ev": float(cl.ke_ev.iloc[0] - cl.ke_ev.iloc[-1]),
            "var_growth_ev": float(wp.var_term_ev.iloc[-1] - wp.var_term_ev.iloc[0]),
        }
        for e in ("T1", "T2", "CL"):
            row[f"S_{e}"] = float(f[e]["S"])
            row[f"S_{e}_err"] = float(f[e]["sigma"])
            row[f"r2_{e}"] = float(f[e]["r2"])
        row["ratio_T1"] = row["S_T1"] / row["S_CL"]
        row["ratio_T2"] = row["S_T2"] / row["S_CL"]
        # where -dT2/ds turns positive (NaN if it never was negative)
        t, S2 = wp.t.to_numpy(), d["S_T2"]
        g = ~np.isnan(S2)
        pos = np.where(g & (S2 > 0))[0]
        row["t_zero_cross_T2"] = (float(t[pos[0]])
                                  if pos.size and np.nanmin(S2[g]) < 0 else np.nan)
        rows.append(row)

    df = pd.DataFrame(rows)
    if len(df):
        df = df.set_index("rung").loc[[r for r in RUNGS if r in set(df.rung)]] \
               .reset_index() if "rung" in df.columns else df
    return data, df


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

# quarter plasma period at r_s = 3: omega_p = sqrt(4 pi n) = 0.333 a.u.
OMEGA_P = 0.3334
T_WAKE = 0.5 * np.pi / OMEGA_P          # ~4.7 a.u.


def _win_shade(ax):
    """Grey = the fit window; pink = the wake build-up no-go zone.

    Both are drawn on every panel of F2 so the reader can check the window
    against the raw slope rather than trust the caption.
    """
    ax.axvspan(0, T_WAKE, color="#f6dede", zorder=0, lw=0)
    ax.axvspan(*WIN, color="0.85", zorder=0, lw=0)


def f1_cumulative_loss(data, out, dpi):
    """Energy actually lost, vs distance travelled — the raw evidence.

    Three panels share a y-axis on purpose: the whole point is that the SAME
    physical run yields three different loss curves depending on which kinetic
    energy you call the projectile's.
    """
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    for ax, est in zip(axes, ("T1", "T2", "CL")):
        for tag, d in data.items():
            if est == "CL":
                s, e = d["s_cl"], d["cl"].ke_ev.to_numpy()
            else:
                s, e = d["s_wp"], d["wp"][COL[est]].to_numpy()
            ax.plot(s, e[0] - e, lw=1.4, color=RUNG_COLOR[tag], label=NICE[tag])
        ax.axhline(0, color="0.5", lw=0.7, ls=":")
        ax.set_xlabel(style.axis_label("length", "path travelled"))
        ax.set_title(EST[est]["label"], fontsize=9)
    axes[0].set_ylabel("energy lost (eV)")
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper left",
                   title="coupling", title_fontsize=7.5)
    fig.suptitle("Energy lost by the projectile, by definition of its kinetic energy",
                 fontsize=10)
    fig.tight_layout()
    style.save_presentation(fig, out / "F1_energy_loss_vs_path.png", dpi=dpi)


def f2_local_stopping(data, out, dpi):
    """The instantaneous stopping power -dE/ds(t): where a fit window is legal.

    The shaded band is the window every quoted S uses. Drawing it on the raw
    slope makes the window auditable instead of a number in a caption: a reader
    sees directly that it sits after the wake build-up peak and inside the flat
    part of the classical curve.
    """
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    for ax, est in zip(axes, ("T1", "T2", "CL")):
        for tag, d in data.items():
            fr = d["cl"] if est == "CL" else d["wp"]
            ax.plot(fr.t, d[f"S_{est}"], lw=1.4, color=RUNG_COLOR[tag],
                    label=NICE[tag])
        _win_shade(ax)
        ax.axhline(0, color="0.4", lw=0.8)
        ax.set_xlabel("time (a.u.)")
        ax.set_title(EST[est]["label"], fontsize=9)
    axes[0].set_ylabel(style.axis_label("stopping_power", r"$-\,\mathrm{d}E/\mathrm{d}s$"))
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper right")
    axes[0].annotate("pink: wake still\nbuilding up\n" r"($t<\pi/2\omega_p$)",
                     xy=(0.6, -0.085), fontsize=7, color="#a04040")
    axes[1].annotate(r"$S(T_2)<0$ here — the packet gains"
                     "\nmomentum spread faster than\nit loses drift",
                     xy=(6.2, -0.088), fontsize=7, color="0.25")
    fig.suptitle("Instantaneous stopping power; grey band = the fit window used for every quoted $S$",
                 fontsize=10)
    fig.tight_layout()
    style.save_presentation(fig, out / "F2_local_stopping_vs_time.png", dpi=dpi)


def f3_stopping_vs_coupling(df, out, dpi):
    """THE HEADLINE: S(T1) and S(T2) against measured coupling, classical on both.

    Two panels rather than one axis, because the user asked for the definitions
    separately — and because overlaying them hides that they approach the
    classical curve from OPPOSITE sides.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharey=True)
    x = df.f_wall_mean.to_numpy()
    xerr = np.vstack([x - df.f_wall_lo, df.f_wall_hi - x])
    for ax, est in zip(axes, ("T1", "T2")):
        ax.errorbar(x, df.S_CL, xerr=xerr, yerr=df.S_CL_err, marker="s", ms=5,
                    lw=1.4, color=EST["CL"]["color"], label=EST["CL"]["label"],
                    capsize=2)
        ax.errorbar(x, df[f"S_{est}"], xerr=xerr, yerr=df[f"S_{est}_err"],
                    marker="o", ms=5, lw=1.4, color=EST[est]["color"],
                    label=EST[est]["label"], capsize=2)
        for xi, yi, tag in zip(x, df[f"S_{est}"], df.rung):
            ax.annotate(SHORT[tag], (xi, yi), textcoords="offset points",
                        xytext=(4, -10), fontsize=6.5, color="0.3")
        ax.axhline(0, color="0.5", lw=0.7, ls=":")
        ax.set_xlim(-0.02, 1.14)      # headroom for the "filled" point label
        ax.set_xlabel(r"measured coupling  $f_{\rm wall}$   (fraction of $|\psi|^2$ in jellium)")
        ax.legend(fontsize=8, frameon=False, loc="upper left")
    axes[0].set_ylabel(style.axis_label("stopping_power", "$S$"))
    axes[0].set_title("(a) drift definition", fontsize=9)
    axes[1].set_title("(b) total-kinetic definition", fontsize=9)
    fig.suptitle(f"Stopping power vs interaction strength   (fit window $t={WIN[0]:.0f}-{WIN[1]:.0f}$ a.u.)",
                 fontsize=10)
    fig.tight_layout()
    style.save_presentation(fig, out / "F3_stopping_vs_coupling.png", dpi=dpi)


def f4_ratio_convergence(df, out, dpi):
    """How well each quantum definition reproduces the classical answer.

    One axis this time, because the CROSSING is the result: the two definitions
    bracket the classical value at weak coupling and agree with each other at
    strong coupling.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    x = df.f_wall_mean.to_numpy()
    xerr = np.vstack([x - df.f_wall_lo, df.f_wall_hi - x])
    for est in ("T1", "T2"):
        r = df[f"S_{est}"] / df.S_CL
        rerr = r * np.hypot(df[f"S_{est}_err"] / df[f"S_{est}"],
                            df.S_CL_err / df.S_CL)
        ax.errorbar(x, r, xerr=xerr, yerr=rerr, marker="o", ms=6, lw=1.6,
                    color=EST[est]["color"], capsize=2,
                    label=rf"$S({est})\,/\,S_{{\rm cl}}$")
    ax.axhline(1.0, color="0.4", lw=1.0, ls="--")
    ax.annotate("perfect agreement with classical", xy=(0.03, 1.02), fontsize=7.5,
                color="0.35")
    conv = df[df.f_wall_mean > 0.9]
    if len(conv):
        y = float(conv[["ratio_T1", "ratio_T2"]].to_numpy().mean())
        ax.axhline(y, color="0.55", lw=0.9, ls=":")
        ax.annotate(f"both converge to {y:.2f}", xy=(0.03, y + 0.02),
                    fontsize=7.5, color="0.35")
    # rung identity goes on a SECONDARY TOP axis, not as text near the bottom
    # spine — at r04/r00 the couplings are 0.93 and 0.995, so bottom labels
    # collided with each other and with the tick labels.
    top = ax.secondary_xaxis("top")
    top.set_xticks(list(x))
    top.set_xticklabels([SHORT[t] for t in df.rung], fontsize=7, color="0.3")
    top.set_xlabel("geometry (bore radius in units of $\\sigma_{\\rm WP}$)",
                   fontsize=8, color="0.3")
    ax.set_xlabel(r"measured coupling  $f_{\rm wall}$")
    ax.set_ylabel(r"$S_{\rm wavepacket}\,/\,S_{\rm classical}$")
    ax.set_ylim(-0.05, 1.15)
    ax.legend(fontsize=9, frameon=False, loc="upper right")
    ax.set_title("The two definitions disagree at weak coupling\nand converge at strong coupling",
                 fontsize=9.5)
    fig.tight_layout()
    style.save_presentation(fig, out / "F4_ratio_convergence.png", dpi=dpi)


def f5_variance_mechanism(data, df, out, dpi):
    """WHY they differ: T2 - T1 = var(p)/2m, the energy held in momentum spread.

    Left: the var term itself. Right: what fraction of the drift loss it accounts
    for — the single number that explains the convergence in F4.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    ax = axes[0]
    for tag, d in data.items():
        v = d["wp"].var_term_ev.to_numpy()
        ax.plot(d["s_wp"], v - v[0], lw=1.4, color=RUNG_COLOR[tag], label=NICE[tag])
    ax.axhline(0, color="0.5", lw=0.7, ls=":")
    ax.set_xlabel(style.axis_label("length", "path travelled"))
    ax.set_ylabel(r"$\Delta\,[\mathrm{var}(p)/2m]$  (eV)")
    ax.set_title("(a) energy diverted into momentum spread", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")

    ax = axes[1]
    frac = df.var_growth_ev / df.dE_T1_ev
    ax.plot(df.f_wall_mean, 100 * frac, marker="o", ms=6, lw=1.6, color="tab:purple")
    for xi, yi, tag in zip(df.f_wall_mean, 100 * frac, df.rung):
        ax.annotate(SHORT[tag], (xi, yi), textcoords="offset points",
                    xytext=(-8, 8), ha="right", fontsize=7, color="0.3")
    ax.set_xlim(0.02, 1.10)
    ax.set_xlabel(r"measured coupling  $f_{\rm wall}$")
    ax.set_ylabel(r"$\Delta\mathrm{var}(p)/2m$  as % of drift loss $\Delta T_1$")
    ax.set_title("(b) it dominates when coupling is weak", fontsize=9)
    ax.set_ylim(0, None)
    fig.suptitle(r"The gap between the definitions: $T_2-T_1=\mathrm{var}(p)/2m$",
                 fontsize=10)
    fig.tight_layout()
    style.save_presentation(fig, out / "F5_variance_mechanism.png", dpi=dpi)


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()

    style.apply_theme()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[results] loading rungs {RUNGS}")
    data, df = collect()
    if not len(df):
        print("[results] no rung present — nothing to draw")
        return 1

    df.to_csv(OUT / "results_summary.csv", index=False)
    f1_cumulative_loss(data, OUT, a.dpi)
    f2_local_stopping(data, OUT, a.dpi)
    f3_stopping_vs_coupling(df, OUT, a.dpi)
    f4_ratio_convergence(df, OUT, a.dpi)
    f5_variance_mechanism(data, df, OUT, a.dpi)

    (OUT / "window.json").write_text(json.dumps(
        {"t0_au": WIN[0], "t1_au": WIN[1], "v_floor": V_FLOOR,
         "rationale": "after wake build-up (pi/2 omega_p, omega_p=0.333 au); "
                      "before classical v drops below 0.85 v0 (t=20.6 at r00)"},
        indent=2))
    print(f"[results] wrote 5 figures + results_summary.csv -> {OUT}")
    print(df[["rung", "f_wall_mean", "S_T1", "S_T2", "S_CL",
              "ratio_T1", "ratio_T2", "r2_T2"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
