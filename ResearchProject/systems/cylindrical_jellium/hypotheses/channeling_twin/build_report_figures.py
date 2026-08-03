#!/usr/bin/env python3
"""Report-ready figures for the cylindrical-jellium channeling study.

Emits ONE PNG PER PANEL (never a pre-composed multi-panel image) so the panels
can be collated, reordered and resized in LaTeX later — the report-figures rule
"always produce individual plots, composition is a downstream LaTeX concern".

Three output sets, identical filenames across the two run sets so a reader can
diff them panel-for-panel:

    report_figures/setup/   ground-state density slices (xz, yz; linear + log)
    report_figures/twin/    the uncorrected LDA twin   (classical + wavepacket)
    report_figures/sic/     the SIC-PZ corrected wavepacket, SAME classical

The `sic` set reuses the twin's CLASSICAL half deliberately: the correction only
touches the wavepacket's Hamiltonian, so the classical reference is unchanged and
re-plotting it from the same CSVs guarantees the two sets share one reference.

HOUSE STANDARD APPLIED (report-figures skill / ADR 0004)
-------------------------------------------------------
- canonical theme `inqview.visualisation.style` — no ad-hoc rcParams anywhere;
- NO on-canvas titles (report figures carry their title in the LaTeX caption);
  no interpretive text on the canvas;
- time axes in **fs** (the canonical unit) with a secondary top axis in a.u., so
  the fit windows — which were chosen in a.u. — stay readable without the reader
  converting anything;
- **shared axis limits across the twin and sic sets**, computed in a first pass
  and applied in a second, so the two sets are visually comparable panel-for-panel
  (a per-figure autoscale would silently rescale the difference away);
- density / momentum maps ship a LINEAR and a LOG variant as separate files;
  signed difference maps use a symmetric SymLogNorm and a diverging cmap;
- colorbars carry a shared mathtext x10^n offset, <= 5 ticks, 2 s.f.;
- fixed canvas (`bbox_inches=None`) so the authored inch-width is the on-page
  width and fonts render at their true pt.

Usage
-----
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 build_report_figures.py
    ... [--only twin|sic|setup] [--dpi 600]
"""
from __future__ import annotations

import argparse
import importlib
import math
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import MaxNLocator, ScalarFormatter

HERE = Path(__file__).resolve().parent
SYSTEM = HERE.parents[1]                       # .../systems/cylindrical_jellium
SCRIPTS = SYSTEM / "scripts"
OUT_ROOT = HERE / "report_figures"

from inqview.visualisation import style           # noqa: E402
from inqview.visualisation.field_io import load_vti  # noqa: E402

AU_TO_FS = 0.024188843265857          # 1 atomic unit of time, in fs

# The user-chosen fit windows (a.u.), 2026-08-02. Kept verbatim from
# build_refined_notebook.py so the report figures and the notebook cannot drift.
WINDOWS = [
    ("T1  9-25",  "T1", (9.0, 25.0)),
    ("T2  21-30", "T2", (21.0, 30.0)),
    ("T2  5-20",  "T2", (5.0, 20.0)),
]
#: TWO colours only, across every stopping-power panel: one for the classical
#: projectile, one for the wavepacket. Windows are distinguished by LINE STYLE,
#: never by hue — so hue carries exactly one meaning (which projectile) and a
#: reader never has to ask whether a colour change means a different estimator or
#: a different run.
COLOR_CL = "tab:blue"
COLOR_WP = "tab:red"
WIN_STYLE = {"T1  9-25": "-", "T2  21-30": "--", "T2  5-20": ":"}

#: The defining formula for each estimator, rendered ON the panel next to its
#: number. A stopping power is meaningless without saying which derivative it is:
#: the classical number is a bath-side energy gain per unit path, the T1 number a
#: projectile-side drift loss, and they are only comparable because energy closes.
ESTIMATOR_RHS = {
    "T1": r"$-\,\mathrm{d}T_1/\mathrm{d}s$",
    "T2": r"$-\,\mathrm{d}T_2/\mathrm{d}s$",
    "CL": r"$+\,\mathrm{d}E_\mathrm{total}/\mathrm{d}s$",
}

ESTIMATOR_FORMULA = {
    "T1": r"$S=-\,\mathrm{d}T_1/\mathrm{d}s$",
    "T2": r"$S=-\,\mathrm{d}T_2/\mathrm{d}s$",
    "CL": r"$S=+\,\mathrm{d}E_\mathrm{total}/\mathrm{d}s$",
}

# The two run sets. `cl` is the SAME classical run in both — see module docstring.
RUN_SETS = {
    "twin": dict(wp_results=SCRIPTS / "channeling_twin/wp/results",
                 cl_results=SCRIPTS / "channeling_twin/classical/results",
                 wp_name="wp", cl_name="classical"),
    "sic":  dict(wp_results=SCRIPTS / "channeling_sic/wp/results",
                 cl_results=SCRIPTS / "channeling_twin/classical/results",
                 wp_name="wp_sic", cl_name="classical"),
}

GS_VTI = (SCRIPTS / "channeling_twin/classical/results/classical"
          / "raw/vti/density_gs_system/density_gs_system.vti")


# ===========================================================================
# plumbing
# ===========================================================================

def load_refined(wp_results: Path, cl_results: Path):
    """Import `refined` bound to a specific pair of result roots.

    The module reads CHAN_WP_RESULTS / CHAN_CL_RESULTS ONCE, at import, into
    module-level constants — so switching run sets requires purging and
    re-importing, not just re-setting the environment. (Same trap as the
    `_fresh_refined` helper in tests/test_refined.py.)
    """
    os.environ["CHAN_WP_RESULTS"] = str(wp_results)
    os.environ["CHAN_CL_RESULTS"] = str(cl_results)
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    for mod in ("refined", "channeling_stopping", "ks_stopping"):
        sys.modules.pop(mod, None)
    return importlib.import_module("refined")


class Ctx:
    """Everything one run set needs, computed once."""

    def __init__(self, tag: str, cfg: dict):
        self.tag = tag
        self.R = R = load_refined(cfg["wp_results"], cfg["cl_results"])
        self.CS = R.CS
        self.wp = R.wp_frame(cfg["wp_name"])
        self.cl = R.cl_frame(cfg["cl_name"])
        self.icl = R.interactions("classical", cfg["cl_name"])
        self.iwp = R.interactions("wp", cfg["wp_name"])
        self.imp = R.impulse_comparison(self.wp, self.cl)
        self.scl = R.combined_projectile_coupling(self.icl)
        self.swp = R.combined_projectile_coupling(self.iwp)

        md = R.momentum_slices(cfg["wp_name"])
        self.slices = R.nearest_slices(md, [0.0, 15.0, 30.0])

        self.wp_name = cfg["wp_name"]
        steps = R.available_wf_steps(cfg["wp_name"])
        self.wf_steps = [steps[0], steps[len(steps) // 2], steps[-1]]
        # finer transverse binning than the library default — see KPERP_BINS
        self.maps = {s: momentum_map_refined(self, s) for s in self.wf_steps}

        # local stopping power, centred +-1 a.u. stencil (same as the notebook)
        self.cl_path = self.cl.z_unwrapped.to_numpy()
        self.wp_path = self.wp.s_pintegral.to_numpy()
        self.S_cl = local_slope(-self.cl.d_e_total_ev.to_numpy(), self.cl_path,
                                self.cl.t.to_numpy())
        self.S_wp = local_slope(self.wp.T1_drift_ev.to_numpy(), self.wp_path,
                                self.wp.t.to_numpy())
        self.fits = self._fit()

    def _fit(self) -> pd.DataFrame:
        R, wp, cl = self.R, self.wp, self.cl
        ecol = {"T1": "T1_drift_ev", "T2": "T2_total_ev"}
        rows = []
        for label, est, (t0, t1) in WINDOWS:
            f_cl = R.fit_in_window(cl.z_unwrapped.to_numpy(), cl.ke_ev.to_numpy(),
                                   cl.t.to_numpy(), t0, t1)
            f = R.fit_in_window(wp.s_pintegral.to_numpy(), wp[ecol[est]].to_numpy(),
                                wp.t.to_numpy(), t0, t1)
            rows.append({"window": label, "estimator": est, "t0": t0, "t1": t1,
                         "S_wp": f["S"], "sigma_wp": f["sigma"], "r2_wp": f["r2"],
                         "S_cl": f_cl["S"], "sigma_cl": f_cl["sigma"],
                         "ratio": f["S"] / f_cl["S"] if f_cl["S"] else np.nan})
        return pd.DataFrame(rows)


def local_slope(energy_ev, path, t, half=1.0):
    """-dE/ds at each t; NaN where a full centred stencil is unavailable.

    Deliberately not np.gradient + boxcar: a `mode="same"` convolution fabricates
    the first and last half-window from a truncated kernel, and those are exactly
    the early-time points a window choice is most sensitive to.
    """
    energy_ev = np.asarray(energy_ev, float)
    path = np.asarray(path, float)
    t = np.asarray(t, float)
    out = np.full(t.size, np.nan)
    for i in range(t.size):
        if t[i] - half < t[0] or t[i] + half > t[-1]:
            continue
        lo = np.searchsorted(t, t[i] - half)
        hi = np.searchsorted(t, t[i] + half) - 1
        ds = path[hi] - path[lo]
        if ds != 0:
            out[i] = -(energy_ev[hi] - energy_ev[lo]) / ds
    return out


def t_axis(ax, *, secondary=True):
    """Canonical time axis: fs on the bottom, a.u. on top.

    The project unit standard is fs; the fit windows were chosen in a.u. Carrying
    both removes the conversion from the reader without putting an off-standard
    unit on the primary axis.
    """
    ax.set_xlabel(style.axis_label("time", "$t$"))
    if secondary:
        ax.tick_params(top=False)
        sec = ax.secondary_xaxis(
            "top", functions=(lambda v: v / AU_TO_FS, lambda v: v * AU_TO_FS))
        sec.set_xlabel("$t$ (a.u.)", labelpad=4)
    return ax


#: fixed inch geometry for every LINE panel — one identical data box across all
#: of them, so the panels align when collated. `style.figure_one_col()` cannot be
#: used directly: its rect leaves ~3.5 % headroom, which clips the secondary
#: a.u. axis (the style module documents this failure mode for titles/colorbars).
_LINE_W, _LINE_H = 3.5, 3.3
_LINE_MARGIN = dict(left=0.82, right=0.10, bottom=0.50, top=0.46)


def figure_line():
    """One-column line panel with headroom for the secondary a.u. axis."""
    m = _LINE_MARGIN
    ax_w = _LINE_W - m["left"] - m["right"]
    ax_h = _LINE_H - m["bottom"] - m["top"]
    fig = plt.figure(figsize=(_LINE_W, _LINE_H))
    ax = fig.add_axes((m["left"] / _LINE_W, m["bottom"] / _LINE_H,
                       ax_w / _LINE_W, ax_h / _LINE_H))
    return fig, ax


def figure_map(w_span: float, h_span: float, *, width_in: float = 3.5):
    """A 2-D map figure with EQUAL aspect and an explicit colorbar axes.

    A spatial map must not be stretched: `aspect="auto"` silently rescales x
    against z and makes a 40x60 Bohr cell look square. The axes box is therefore
    sized in INCHES from the data aspect, and the colorbar gets its own
    fixed-inch rectangle so its label can never spill off the canvas (the
    fixed-dimension idiom — `bbox_inches=None` at save time).
    """
    left_in, gap_in, cb_in, lab_in = 0.62, 0.10, 0.13, 0.62
    # top_in must clear the colorbar's shared "x10^n" offset text, which is drawn
    # ABOVE the bar and is clipped off-canvas by a bare fixed-rect save.
    bot_in, top_in = 0.52, 0.30
    ax_w = width_in - left_in - gap_in - cb_in - lab_in
    ax_h = ax_w * h_span / w_span
    fig_h = ax_h + bot_in + top_in
    fig = plt.figure(figsize=(width_in, fig_h))
    ax = fig.add_axes((left_in / width_in, bot_in / fig_h,
                       ax_w / width_in, ax_h / fig_h))
    cax = fig.add_axes(((left_in + ax_w + gap_in) / width_in, bot_in / fig_h,
                        cb_in / width_in, ax_h / fig_h))
    return fig, ax, cax


def symlog_ticks(vmax: float, linthresh: float, n_decades: int = 2) -> list:
    """Sparse +/- decade ticks for a SymLogNorm colorbar, plus zero.

    A SymLogNorm colorbar left to its default locator crowds the linear region:
    `0`, `1e-4` and `-1e-4` are drawn within a few points of each other and
    overprint into an unreadable smear. Explicit decade ticks from `vmax` down,
    floored at `linthresh`, keep the bar legible.
    """
    if not (vmax > 0 and linthresh > 0):
        return [0.0]
    top = math.floor(math.log10(vmax))
    lo = max(math.ceil(math.log10(linthresh)), top - n_decades + 1)
    decs = [10.0 ** e for e in range(int(lo), int(top) + 1)]
    return [-d for d in reversed(decs)] + [0.0] + decs


def sci_colorbar(fig, im, cax, label, *, sci=True, ticks=None):
    """Colorbar with a shared x10^n offset, <=5 ticks, no clipped decimals.

    `sci=False` for Log/SymLog norms, which keep their native 10^n decade ticks
    (a ScalarFormatter would destroy them). Pass `ticks` to override the locator.
    """
    cb = fig.colorbar(im, cax=cax, ticks=ticks)
    if sci:
        fmt = ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((0, 0))
        cb.locator = MaxNLocator(5)
        cb.formatter = fmt
        cb.update_ticks()
    cb.set_label(label)
    return cb


# ===========================================================================
# the panels — each draws ONE axes; the driver owns figure creation and saving
# ===========================================================================

def p_01a_trajectories(C, ax):
    cl, wp = C.cl, C.wp
    ax.plot(cl.t * AU_TO_FS, cl.z_unwrapped, lw=1.4, color=COLOR_CL,
            label="classical")
    ax.plot(wp.t * AU_TO_FS, wp.s_centroid, lw=1.2, ls="--", color=COLOR_WP,
            label=r"WP $s_\mathrm{centroid}$")
    ax.plot(wp.t * AU_TO_FS, wp.s_pintegral, lw=1.0, ls=":", color=COLOR_WP,
            label=r"WP $s_{\int\langle p\rangle}$")
    for zf in (-30, 30):
        ax.axhline(zf, color="0.5", lw=0.7, ls="-.")
    ax.set_ylim(-33, 33)
    ax.set_ylabel(style.axis_label("length", "$z$"))
    ax.legend(frameon=False, loc="lower right")
    t_axis(ax)


def p_01b_ehrenfest_residual(C, ax):
    wp = C.wp
    ax.plot(wp.t * AU_TO_FS, wp.ehrenfest_resid * 1e3, lw=1.2, color=COLOR_WP)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(r"$s_\mathrm{centroid}-s_{\int\langle p\rangle}$"
                  "\n" r"($10^{-3}$ Bohr)")
    t_axis(ax)


def p_01c_twin_separation(C, ax):
    cl, wp = C.cl, C.wp
    d = np.interp(cl.t, wp.t, wp.s_pintegral) - cl.z_unwrapped
    d2 = np.interp(cl.t, wp.t, wp.s_centroid) - cl.z_unwrapped
    ax.plot(cl.t * AU_TO_FS, d, lw=1.3, ls=":", color=COLOR_WP,
            label=r"$s_{\int\langle p\rangle}-z_\mathrm{cl}$")
    ax.plot(cl.t * AU_TO_FS, d2, lw=1.1, ls="--", color=COLOR_WP,
            label=r"$s_\mathrm{centroid}-z_\mathrm{cl}$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("length", "separation"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_02a_classical_energy_split(C, ax):
    cl = C.cl
    ax.plot(cl.t * AU_TO_FS, cl.d_e_total_ev, lw=1.3, color="tab:blue",
            label=r"$\Delta E_\mathrm{total}$ (bath)")
    ax.plot(cl.t * AU_TO_FS, cl.d_ke_ev, lw=1.3, color="tab:red",
            label=r"$\Delta(\frac{1}{2}mv^2)$")
    ax.plot(cl.t * AU_TO_FS, cl.closure_ev, lw=1.5, color="k", label="sum")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_02b_classical_closure_residual(C, ax):
    cl = C.cl
    ax.plot(cl.t * AU_TO_FS, cl.closure_ev * 1e6, lw=1.2, color="k")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(r"closure residual ($10^{-6}$ eV)")
    t_axis(ax)


def p_02c_classical_ke_absolute(C, ax):
    cl = C.cl
    ax.plot(cl.t * AU_TO_FS, cl.ke_ev, lw=1.3, color="tab:red")
    ax.axhline(cl.ke_ev.iloc[0], color="0.5", lw=0.7, ls="--")
    ax.set_ylabel(style.axis_label("energy", r"$\frac{1}{2}mv^2$"))
    t_axis(ax)


def p_03a_T1_T2(C, ax):
    wp = C.wp
    ax.plot(wp.t * AU_TO_FS, wp.T1_drift_ev, lw=1.3, color=COLOR_WP,
            label=r"$T_1=\langle p\rangle^2/2m$")
    ax.plot(wp.t * AU_TO_FS, wp.T2_total_ev, lw=1.3, ls="--", color=COLOR_WP,
            label=r"$T_2=T_1+\mathrm{var}(p)/2m$")
    ax.set_ylabel(style.axis_label("energy"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_03b_dT1_dT2(C, ax):
    wp = C.wp
    ax.plot(wp.t * AU_TO_FS, wp.d_T1_ev, lw=1.3, color=COLOR_WP,
            label=r"$\Delta T_1$")
    ax.plot(wp.t * AU_TO_FS, wp.d_T2_ev, lw=1.3, ls="--", color=COLOR_WP,
            label=r"$\Delta T_2$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_03c_var_term(C, ax):
    wp, R = C.wp, C.R
    ax.plot(wp.t * AU_TO_FS, wp.var_term_ev, lw=1.4, color="tab:red",
            label=r"$T_2-T_1$")
    ax.axhline(R.T2_MINUS_T1_FREE_EV, color="0.4", ls="--", lw=1.0,
               label=r"free ($1.28$ eV)")
    ax.set_ylabel(style.axis_label("energy", r"$\mathrm{var}(p)/2m$"))
    ax.legend(frameon=False)
    t_axis(ax)


def _var_components(C):
    wp, HA = C.wp, C.R.HA_TO_EV
    return (0.5 * wp.var_p3d * HA, 0.5 * wp.var_pz * HA, 0.5 * wp.var_perp * HA)


def p_04a_varp_total(C, ax):
    wp, R = C.wp, C.R
    tot, _, _ = _var_components(C)
    ax.plot(wp.t * AU_TO_FS, tot, lw=1.4, color="tab:red", label="3D")
    ax.axhline(R.T2_MINUS_T1_FREE_EV, color="0.4", ls="--", lw=1.0,
               label="free-evolution value")
    ax.set_ylabel(style.axis_label("energy", r"$\mathrm{var}(p)/2m$"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_04b_varp_split(C, ax):
    wp, R = C.wp, C.R
    _, vz, vperp = _var_components(C)
    HA = R.HA_TO_EV
    ax.plot(wp.t * AU_TO_FS, vz, lw=1.3, color="tab:blue",
            label=r"$\mathrm{var}(p_z)/2m$")
    ax.plot(wp.t * AU_TO_FS, vperp, lw=1.3, color="tab:cyan",
            label=r"$\mathrm{var}(p_\perp)/2m$")
    ax.axhline(0.5 * R.VAR_P_FREE * HA, color="tab:blue", ls="--", lw=0.9)
    ax.axhline(1.0 * R.VAR_P_FREE * HA, color="tab:cyan", ls="--", lw=0.9)
    ax.set_ylabel(style.axis_label("energy"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_04c_varp_growth(C, ax):
    wp = C.wp
    tot, vz, vperp = _var_components(C)
    ax.plot(wp.t * AU_TO_FS, tot - tot.iloc[0], lw=1.4, color="tab:red",
            label="total")
    ax.plot(wp.t * AU_TO_FS, vz - vz.iloc[0], lw=1.1, color="tab:blue", label=r"$z$")
    ax.plot(wp.t * AU_TO_FS, vperp - vperp.iloc[0], lw=1.1, color="tab:cyan",
            label=r"$\perp$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta[\mathrm{var}(p)/2m]$"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_05a_kinetic_vs_time(C, ax):
    cl, wp = C.cl, C.wp
    ax.plot(cl.t * AU_TO_FS, cl.d_ke_ev, lw=1.4, color=COLOR_CL,
            label="classical")
    ax.plot(wp.t * AU_TO_FS, wp.d_T1_ev, lw=1.3, color=COLOR_WP,
            label=r"WP $\Delta T_1$")
    ax.plot(wp.t * AU_TO_FS, wp.d_T2_ev, lw=1.3, ls="--", color=COLOR_WP,
            label=r"WP $\Delta T_2$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_05b_kinetic_vs_path(C, ax):
    cl, wp = C.cl, C.wp
    ax.plot(cl.z_unwrapped - cl.z_unwrapped.iloc[0], cl.d_ke_ev, lw=1.4,
            color=COLOR_CL, label="classical")
    s = wp.s_pintegral - wp.s_pintegral.iloc[0]
    ax.plot(s, wp.d_T1_ev, lw=1.3, color=COLOR_WP, label=r"WP $\Delta T_1$")
    ax.plot(s, wp.d_T2_ev, lw=1.3, ls="--", color=COLOR_WP,
            label=r"WP $\Delta T_2$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_xlabel(style.axis_label("length", "path travelled"))
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    ax.legend(frameon=False)


def _loss_curves(C):
    """Energy LOST by the projectile, all three estimators, sign-aligned.

    The three quantities are defined with opposite natural signs — the classical
    projectile's kinetic energy FALLS while the bath's total RISES — so plotting
    the raw deltas together (panels 05a/05b) puts the comparison below zero and
    makes the slopes read backwards against the stopping power, which is
    +dE_lost/ds. Negating all three puts them on one positive axis where the
    steeper curve is unambiguously the larger stopping power.
    """
    cl, wp = C.cl, C.wp
    return [
        (cl.t, cl.z_unwrapped - cl.z_unwrapped.iloc[0], -cl.d_ke_ev,
         COLOR_CL, "-", 1.5, r"classical $-\Delta(\frac{1}{2}mv^2)$"),
        (wp.t, wp.s_pintegral - wp.s_pintegral.iloc[0], -wp.d_T1_ev,
         COLOR_WP, "-", 1.4, r"WP $-\Delta T_1$"),
        (wp.t, wp.s_pintegral - wp.s_pintegral.iloc[0], -wp.d_T2_ev,
         COLOR_WP, "--", 1.4, r"WP $-\Delta T_2$"),
    ]


def p_05c_energy_loss_vs_time(C, ax):
    for t, _s, y, col, ls, lw, lab in _loss_curves(C):
        ax.plot(t * AU_TO_FS, y, lw=lw, ls=ls, color=col, label=lab)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel("energy lost by projectile (eV)")
    ax.legend(frameon=False, loc="upper left")
    t_axis(ax)


def p_05d_energy_loss_vs_path(C, ax):
    for _t, s, y, col, ls, lw, lab in _loss_curves(C):
        ax.plot(s, y, lw=lw, ls=ls, color=col, label=lab)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_xlabel(style.axis_label("length", "path travelled"))
    ax.set_ylabel("energy lost by projectile (eV)")
    ax.legend(frameon=False, loc="upper left")


def _slice_colors():
    return ["tab:blue", "tab:orange", "tab:red"]


def p_06a_momentum_dist_linear(C, ax):
    for (t_act, sl), c in zip(C.slices, _slice_colors()):
        ax.plot(sl.k, sl.n_wp, lw=1.2, color=c,
                label=f"$t={t_act*AU_TO_FS:.2f}$ fs")
    ax.axvline(C.CS.V0, color="0.4", ls="--", lw=0.9)
    ax.set_xlim(0, 2.5 * C.CS.V0)
    ax.set_xlabel(style.axis_label("momentum", "$k$"))
    ax.set_ylabel(r"$|\psi_\mathrm{WP}(k)|^2$")
    ax.legend(frameon=False)


def p_06b_momentum_dist_log(C, ax):
    for (t_act, sl), c in zip(C.slices, _slice_colors()):
        ax.semilogy(sl.k, np.maximum(sl.n_wp, 1e-16), lw=1.1, color=c,
                    label=f"$t={t_act*AU_TO_FS:.2f}$ fs")
    ax.axvline(C.CS.V0, color="0.4", ls="--", lw=0.9)
    ax.set_xlim(0, 2.5 * C.CS.V0)
    ax.set_xlabel(style.axis_label("momentum", "$k$"))
    ax.set_ylabel(r"$|\psi_\mathrm{WP}(k)|^2$")
    ax.legend(frameon=False)


def p_06c_varp_with_slices(C, ax):
    wp, R = C.wp, C.R
    ax.plot(wp.t * AU_TO_FS, wp.var_pz, lw=1.3, color="tab:red",
            label=r"$\mathrm{var}(p_z)$")
    ax.plot(wp.t * AU_TO_FS, wp.var_perp / 2.0, lw=1.1, color="tab:cyan",
            label=r"$\mathrm{var}(p_\perp)/2$")
    ax.axhline(R.VAR_P_FREE, color="0.4", ls="--", lw=1.0, label="free value")
    for (t_act, _), c in zip(C.slices, _slice_colors()):
        ax.axvline(t_act * AU_TO_FS, color=c, lw=0.9, ls=":")
    ax.set_ylabel(r"$\sigma_p^2$ (a.u.)")
    ax.legend(frameon=False)
    t_axis(ax)


def p_08a_kz_marginal(C, ax):
    for s, c in zip(C.wf_steps, _slice_colors()):
        kz, _, P = C.maps[s]
        ax.plot(kz, P.sum(axis=1), lw=1.2, color=c,
                label=f"$t={s*C.CS.DT*AU_TO_FS:.2f}$ fs")
    ax.axvline(C.CS.V0, color="0.4", ls="--", lw=0.9)
    ax.set_xlim(0.8, 3.0)
    ax.set_xlabel(style.axis_label("momentum", "$k_z$"))
    ax.set_ylabel(r"$P(k_z)$")
    ax.legend(frameon=False)


def p_08b_kperp_marginal(C, ax):
    for s, c in zip(C.wf_steps, _slice_colors()):
        _, kp, P = C.maps[s]
        ax.plot(kp, P.sum(axis=0), lw=1.2, color=c,
                label=f"$t={s*C.CS.DT*AU_TO_FS:.2f}$ fs")
    ax.set_xlim(0, 1.6)
    ax.set_xlabel(style.axis_label("momentum", r"$k_\perp$"))
    ax.set_ylabel(r"$P(k_\perp)$")
    ax.legend(frameon=False)


_PALETTE = {"e_ss": "tab:blue", "e_ps": "tab:red", "e_pp": "tab:purple",
            "e_sb": "tab:green", "e_pb": "tab:brown"}

#: SYMBOL-ONLY legend labels. `refined.TERM_LABEL` appends a gloss
#: ("$\\Delta E_{SB}$  bath$-$background") which is right for a notebook but ~4x
#: too wide for a 3.5 in report panel: a 2-column legend of those ran off both
#: edges and silently CLIPPED three of the five entries, leaving half the curves
#: unidentified. The gloss belongs in the LaTeX caption (report-figures rule 1:
#: no explanatory text on the canvas), so only the symbol is kept here.
TERM_SHORT = {t: rf"$\Delta E_{{{t.split('_')[1].upper()}}}$"
              for t in ("e_ss", "e_ps", "e_pp", "e_sb", "e_pb")}

#: Fractional headroom added to the TOP of the shared y-limit, per panel, so an
#: in-axes legend sits in reserved empty space instead of over the data. Applied
#: in `draw_set` AFTER the shared limits, which would otherwise overwrite it.
YHEADROOM = {
    "09a_interactions_classical": 0.30,
    "09b_interactions_wp": 0.30,
    "09c_interactions_differing": 0.32,
    "09d_interactions_both": 0.40,
    "09e_interactions_both_projectile": 0.38,
    "10a_ps_pb_separate": 0.32,
    "10b_ps_pb_combined": 0.22,
    "13a_T1_window_fits": 0.10,
    "13b_T2_window_fits": 0.10,
    "13d_classical_fit_vs_path": 0.12,
    "13e_classical_energy_vs_time": 0.18,
    "12c_local_stopping_power": 0.25,
}


def _legend(ax, **kw):
    """Legend on a minimal white bbox so it reads over any curve it must cross."""
    kw.setdefault("frameon", True)
    kw.setdefault("framealpha", 0.88)
    kw.setdefault("facecolor", "white")
    kw.setdefault("edgecolor", "none")
    kw.setdefault("borderpad", 0.25)
    kw.setdefault("handlelength", 1.4)
    kw.setdefault("columnspacing", 1.0)
    kw.setdefault("handletextpad", 0.5)
    return ax.legend(**kw)


def _interaction_panel(C, ax, df):
    n = 0
    for term in C.R.INTERACTION_TERMS:
        col = f"d_{term}_ev"
        if col not in df.columns:
            continue
        y = df[col].to_numpy()
        if np.allclose(y, 0.0):
            continue                       # identically zero: uniform background
        ax.plot(df.time_au * AU_TO_FS, y, lw=1.2, color=_PALETTE[term],
                label=TERM_SHORT[term])
        n += 1
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    _legend(ax, ncol=3 if n > 3 else n, loc="upper center", fontsize=8)
    t_axis(ax)


def p_09a_interactions_classical(C, ax):
    _interaction_panel(C, ax, C.icl)


def p_09b_interactions_wp(C, ax):
    _interaction_panel(C, ax, C.iwp)


def _interactions_both(C, ax, terms):
    """Both halves on one axes: hue = TERM, line style = representation.

    This is the REVERSE of the stopping-power panels, deliberately. There, hue
    meant "which projectile" because only two curves existed. Here there are five
    physically distinct pairwise terms, and collapsing them onto two colours
    would make E_SS indistinguishable from E_PS — the encoding has to follow the
    quantity being compared, so the representation moves to line style.

    Two terms are absent from the classical half BY CONSTRUCTION, and their
    absence is the physics: a rigid classical cloud has a fixed self-energy, so
    dE_PP = 0 exactly, and it does not move relative to a z-uniform background,
    so dE_PB = 0 exactly. Only the wavepacket has them.
    """
    from matplotlib.lines import Line2D
    seen = []
    for df, ls in ((C.icl, "-"), (C.iwp, "--")):
        for term in terms:
            col = f"d_{term}_ev"
            if col not in df.columns:
                continue
            y = df[col].to_numpy()
            if np.allclose(y, 0.0):
                continue                   # identically zero -> nothing to draw
            ax.plot(df.time_au * AU_TO_FS, y, lw=1.2, ls=ls, color=_PALETTE[term])
            if term not in seen:
                seen.append(term)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    handles = [Line2D([], [], color=_PALETTE[t], lw=1.4, label=TERM_SHORT[t])
               for t in seen]
    handles += [Line2D([], [], color="0.35", lw=1.4, ls="-", label="classical"),
                Line2D([], [], color="0.35", lw=1.4, ls="--", label="wavepacket")]
    _legend(ax, handles=handles, ncol=3, loc="upper center", fontsize=7.5)
    t_axis(ax)


def p_09d_interactions_both(C, ax):
    _interactions_both(C, ax, C.R.INTERACTION_TERMS)


def p_09e_interactions_both_projectile(C, ax):
    """The same comparison restricted to the PROJECTILE-coupling terms.

    On the full-term axes E_SS and E_SB span +-15 eV and compress everything
    else; the terms that actually differ between the twins (E_PS, E_PP, E_PB)
    live within +-3 eV. This panel is that window, so the classical-vs-WP
    difference is legible rather than a thickening of the zero line.
    """
    _interactions_both(C, ax, ("e_ps", "e_pp", "e_pb"))


def p_09c_interactions_differing(C, ax):
    icl, iwp = C.icl, C.iwp
    ax.plot(iwp.time_au * AU_TO_FS, iwp.d_e_pp_ev, lw=1.4, color="tab:purple",
            label=r"WP $\Delta E_{PP}$")
    ax.plot(iwp.time_au * AU_TO_FS, iwp.d_e_ps_ev, lw=1.1, ls="--", color="tab:red",
            label=r"WP $\Delta E_{PS}$")
    ax.plot(icl.time_au * AU_TO_FS, icl.d_e_ps_ev, lw=1.1, ls=":", color="tab:blue",
            label=r"cl. $\Delta E_{PS}$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    _legend(ax, ncol=2, loc="upper center", fontsize=8)
    t_axis(ax)


def p_10a_ps_pb_separate(C, ax):
    icl, iwp = C.icl, C.iwp
    ax.plot(icl.time_au * AU_TO_FS, icl.d_e_ps_ev, lw=1.1, ls=":", color=COLOR_CL,
            label=r"cl. $E_{PS}$")
    ax.plot(icl.time_au * AU_TO_FS, icl.d_e_pb_ev, lw=1.1, ls="-.", color=COLOR_CL,
            label=r"cl. $E_{PB}$")
    ax.plot(iwp.time_au * AU_TO_FS, iwp.d_e_ps_ev, lw=1.1, ls=":", color=COLOR_WP,
            label=r"WP $E_{PS}$")
    ax.plot(iwp.time_au * AU_TO_FS, iwp.d_e_pb_ev, lw=1.1, ls="-.", color=COLOR_WP,
            label=r"WP $E_{PB}$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    _legend(ax, ncol=2, loc="upper center", fontsize=8)
    t_axis(ax)


def p_10b_ps_pb_combined(C, ax):
    icl, iwp, scl, swp = C.icl, C.iwp, C.scl, C.swp
    n = min(len(scl), len(swp))
    ax.plot(icl.time_au * AU_TO_FS, scl, lw=1.4, color=COLOR_CL, label="classical")
    ax.plot(iwp.time_au * AU_TO_FS, swp, lw=1.4, color=COLOR_WP, label="wavepacket")
    ax.plot(icl.time_au.to_numpy()[:n] * AU_TO_FS, swp[:n] - scl[:n], lw=1.0,
            color="k", label="WP $-$ cl.")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta(E_{PS}+E_{PB})$"))
    _legend(ax, ncol=3, loc="upper center", fontsize=8)
    t_axis(ax)


def p_11a_cumulative_impulse(C, ax):
    imp = C.imp
    ax.plot(imp.t * AU_TO_FS, imp.dp_cl, lw=1.4, color=COLOR_CL, label="classical")
    ax.plot(imp.t * AU_TO_FS, imp.dp_wp, lw=1.4, color=COLOR_WP,
            label=r"WP $\Delta\langle p\rangle$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(r"$\Delta p$ (a.u.)")
    ax.legend(frameon=False)
    t_axis(ax)


def p_11b_impulse_ratio(C, ax):
    imp, wp, iwp = C.imp, C.wp, C.iwp
    epp = np.interp(wp.t, iwp.time_au, iwp.e_pp_ev)
    good = imp.t > 2.0
    ax.plot(imp.t[good] * AU_TO_FS, imp.impulse_ratio[good], lw=1.4, color="k",
            label=r"$\Delta p_\mathrm{WP}/\Delta p_\mathrm{cl}$")
    ax.plot(wp.t[good] * AU_TO_FS, (epp / epp[0])[good], lw=1.0, ls="--",
            color="tab:purple", label=r"$E_{PP}(t)/E_{PP}(0)$")
    ax.plot(wp.t[good] * AU_TO_FS, wp.f_bore[good], lw=1.0, ls=":", color="tab:red",
            label=r"$f_\mathrm{bore}$")
    ax.axhline(1.0, color="0.5", lw=0.7)
    ax.set_ylabel("ratio")
    ax.legend(frameon=False)
    t_axis(ax)


def p_11c_energy_gap(C, ax):
    cl, wp = C.cl, C.wp
    gap = np.interp(cl.t, wp.t, wp.d_T1_ev) - cl.d_ke_ev
    ax.plot(cl.t * AU_TO_FS, gap, lw=1.4, color="tab:brown")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label(
        "energy", r"$\Delta T_1-\Delta(\frac{1}{2}mv^2)$"))
    t_axis(ax)


def p_12a_fit_targets_time(C, ax):
    cl, wp = C.cl, C.wp
    ax.plot(cl.t * AU_TO_FS, cl.d_e_total_ev, lw=1.4, color=COLOR_CL,
            label=r"cl. $\Delta E_\mathrm{total}$")
    ax.plot(wp.t * AU_TO_FS, wp.d_T1_ev, lw=1.4, color=COLOR_WP,
            label=r"WP $\Delta T_1$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    ax.legend(frameon=False)
    t_axis(ax)


def p_12b_fit_targets_path(C, ax):
    cl, wp = C.cl, C.wp
    ax.plot(C.cl_path - C.cl_path[0], cl.d_e_total_ev, lw=1.4, color=COLOR_CL,
            label=r"cl. $\Delta E_\mathrm{total}$")
    ax.plot(C.wp_path - C.wp_path[0], -wp.d_T1_ev, lw=1.4, color=COLOR_WP,
            label=r"WP $-\Delta T_1$")
    ax.set_xlabel(style.axis_label("length", "path travelled"))
    ax.set_ylabel("energy lost by projectile (eV)")
    ax.legend(frameon=False)


def p_12c_local_stopping_power(C, ax):
    cl, wp, CS = C.cl, C.wp, C.CS
    ax.plot(cl.t * AU_TO_FS, C.S_cl, lw=1.3, color=COLOR_CL, label="classical")
    ax.plot(wp.t * AU_TO_FS, C.S_wp, lw=1.3, color=COLOR_WP, label="wavepacket")
    ax.axhline(0, color="0.5", lw=0.7)
    # HARD UPPER BOUND: the projectile re-approaches the disturbance it made at
    # launch, through the periodic image. The wake is 36 Bohr in a 60 Bohr box,
    # so this is not a small correction.
    d_img = np.abs(C.cl_path - (CS.LAUNCH_Z + CS.LZ))
    t_half = cl.t[d_img < CS.LAMBDA_P / 2.0]
    if len(t_half):
        ax.axvspan(float(t_half.iloc[0]) * AU_TO_FS, cl.t.iloc[-1] * AU_TO_FS,
                   color="tab:red", alpha=0.10, lw=0)
    ax.set_ylabel(style.axis_label("stopping_power", "local $S$"))
    _legend(ax, fontsize=8, loc="lower right")
    _formula_note(ax, "classical: " + ESTIMATOR_FORMULA["CL"] + "\n"
                      + "WP: " + ESTIMATOR_FORMULA["T1"], "upper left")
    t_axis(ax)


def _formula_note(ax, text, loc="lower left"):
    """Short formula/value note in axes coords, on a minimal white bbox."""
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    y, va = (0.04, "bottom") if "lower" in loc else (0.96, "top")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none",
                      alpha=0.88))


def _window_fit_panel(C, ax, est, ycol, ylabel):
    """WP fit target vs path, with the regressed line over each window."""
    wp, fits = C.wp, C.fits
    ax.plot(wp.s_pintegral, wp[ycol], lw=0.8, color="0.85")
    notes = [ESTIMATOR_FORMULA[est]]
    for label, e, (t0, t1) in WINDOWS:
        if e != est:
            continue
        m = (wp.t >= t0) & (wp.t <= t1)
        ax.plot(wp.s_pintegral[m], wp[ycol][m], lw=1.8, color=COLOR_WP,
                ls=WIN_STYLE[label], label=f"{t0:g}–{t1:g} a.u.")
        r = fits[fits.window == label].iloc[0]
        xs = wp.s_pintegral[m].to_numpy()
        ax.plot(xs, wp[ycol][m].to_numpy()[0] - r["S_wp"] * (xs - xs[0]),
                lw=0.9, ls="-", color="k", alpha=0.75)
        notes.append(rf"{t0:g}–{t1:g}: ${r['S_wp']:.4f}\pm{r['sigma_wp']:.4f}$")
    ax.set_xlabel(style.axis_label("length", "path"))
    ax.set_ylabel(ylabel)
    _legend(ax, fontsize=8, loc="upper right")
    _formula_note(ax, "\n".join(notes), "lower left")


def p_13a_T1_window_fits(C, ax):
    _window_fit_panel(C, ax, "T1", "T1_drift_ev",
                      style.axis_label("energy", "$T_1$"))


def p_13b_T2_window_fits(C, ax):
    _window_fit_panel(C, ax, "T2", "T2_total_ev",
                      style.axis_label("energy", "$T_2$"))


def p_13d_classical_fit_vs_path(C, ax):
    """The CLASSICAL stopping power, fitted the same way as the WP.

    The regression lives on the PATH axis because `S = dE/ds`; a straight line
    here has slope `S` directly. Fitted over the SAME three windows as the
    wavepacket, because both projectiles decelerate — a classical number from one
    window compared against a WP number from another compares two different
    velocities and means nothing.
    """
    cl, fits = C.cl, C.fits
    s = cl.z_unwrapped - cl.z_unwrapped.iloc[0]
    ax.plot(s, cl.d_e_total_ev, lw=0.8, color="0.85")
    notes = [ESTIMATOR_FORMULA["CL"]]
    for label, _est, (t0, t1) in WINDOWS:
        m = (cl.t >= t0) & (cl.t <= t1)
        ax.plot(s[m], cl.d_e_total_ev[m], lw=1.8, color=COLOR_CL,
                ls=WIN_STYLE[label], label=f"{t0:g}–{t1:g} a.u.")
        r = fits[fits.window == label].iloc[0]
        xs = s[m].to_numpy()
        ax.plot(xs, cl.d_e_total_ev[m].to_numpy()[0] + r["S_cl"] * (xs - xs[0]),
                lw=0.9, ls="-", color="k", alpha=0.75)
        notes.append(rf"{t0:g}–{t1:g}: ${r['S_cl']:.4f}\pm{r['sigma_cl']:.4f}$")
    ax.set_xlabel(style.axis_label("length", "path"))
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E_\mathrm{total}$"))
    _legend(ax, fontsize=8, loc="upper left")
    _formula_note(ax, "\n".join(notes), "lower right")


def p_13e_classical_energy_vs_time(C, ax):
    """Classical bath total energy vs TIME, with the fit windows marked.

    Companion to 13d: the slope that IS the stopping power lives on the path
    axis, so no straight line is drawn here — on a time axis the same fit is a
    curve, because the projectile decelerates and `s(t)` is not linear. Drawing a
    straight line here would misrepresent the fit, so the windows are shaded
    instead and the numbers quoted from the path regression.
    """
    cl = C.cl
    ax.plot(cl.t * AU_TO_FS, cl.d_e_total_ev, lw=1.5, color=COLOR_CL)
    for label, _est, (t0, t1) in WINDOWS:
        ax.axvspan(t0 * AU_TO_FS, t1 * AU_TO_FS, color=COLOR_CL, alpha=0.08, lw=0)
        ax.axvline(t0 * AU_TO_FS, color=COLOR_CL, lw=0.7,
                   ls=WIN_STYLE[label], alpha=0.7)
        ax.axvline(t1 * AU_TO_FS, color=COLOR_CL, lw=0.7,
                   ls=WIN_STYLE[label], alpha=0.7)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E_\mathrm{total}$"))
    _formula_note(ax, ESTIMATOR_FORMULA["CL"] + "\nshaded: fit windows",
                  "upper left")
    t_axis(ax)


def p_13c_stopping_bar(C, ax):
    """WP vs its matched classical reference, per window, with fit uncertainties.

    Two colours only. The estimator formula sits in the tick label, because the
    three rows are NOT the same quantity — two are T2 fits and one is T1 — and a
    bar chart that hides that invites reading them as one series.
    """
    sub = C.fits.reset_index(drop=True)
    y = np.arange(len(sub))
    ax.barh(y - 0.20, sub.S_wp, height=0.38, xerr=sub.sigma_wp, color=COLOR_WP,
            capsize=2.5, error_kw=dict(lw=0.9), label="wavepacket")
    ax.barh(y + 0.20, sub.S_cl, height=0.38, xerr=sub.sigma_cl, color=COLOR_CL,
            capsize=2.5, error_kw=dict(lw=0.9), label="classical")
    # the fit sigmas are ~1e-4, so 3 decimals would print every one as "0.000":
    # the uncertainty must be quoted at the digit it actually lives on.
    for i, r in sub.iterrows():
        ax.text(r.S_wp + r.sigma_wp + 0.004, i - 0.20,
                rf"${r.S_wp:.4f}\pm{r.sigma_wp:.4f}$", va="center", fontsize=7,
                color=COLOR_WP)
        ax.text(r.S_cl + r.sigma_cl + 0.004, i + 0.20,
                rf"${r.S_cl:.4f}\pm{r.sigma_cl:.4f}$", va="center", fontsize=7,
                color=COLOR_CL)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{ESTIMATOR_RHS[r.estimator]}\n{r.t0:g}–{r.t1:g} a.u."
                        for r in sub.itertuples()], fontsize=7)
    ax.set_xlim(0, float(max(sub.S_wp.max(), sub.S_cl.max())) * 1.75)
    ax.invert_yaxis()
    ax.set_xlabel(style.axis_label("stopping_power", "$S$"))
    # legend ABOVE the axes: this panel has no secondary time axis, so the top
    # margin is free, and inside the axes it collided with the bar annotations.
    _legend(ax, fontsize=8, ncol=2, loc="lower center",
            bbox_to_anchor=(0.5, 1.005), frameon=False)


LINE_PANELS = [
    ("01a_trajectories", p_01a_trajectories),
    ("01b_ehrenfest_residual", p_01b_ehrenfest_residual),
    ("01c_twin_separation", p_01c_twin_separation),
    ("02a_classical_energy_split", p_02a_classical_energy_split),
    ("02b_classical_closure_residual", p_02b_classical_closure_residual),
    ("02c_classical_ke_absolute", p_02c_classical_ke_absolute),
    ("03a_T1_T2", p_03a_T1_T2),
    ("03b_dT1_dT2", p_03b_dT1_dT2),
    ("03c_var_term", p_03c_var_term),
    ("04a_varp_total", p_04a_varp_total),
    ("04b_varp_split", p_04b_varp_split),
    ("04c_varp_growth", p_04c_varp_growth),
    ("05a_kinetic_vs_time", p_05a_kinetic_vs_time),
    ("05b_kinetic_vs_path", p_05b_kinetic_vs_path),
    ("05c_energy_loss_vs_time", p_05c_energy_loss_vs_time),
    ("05d_energy_loss_vs_path", p_05d_energy_loss_vs_path),
    ("06a_momentum_dist_linear", p_06a_momentum_dist_linear),
    ("06b_momentum_dist_log", p_06b_momentum_dist_log),
    ("06c_varp_with_slices", p_06c_varp_with_slices),
    ("08a_kz_marginal", p_08a_kz_marginal),
    ("08b_kperp_marginal", p_08b_kperp_marginal),
    ("09a_interactions_classical", p_09a_interactions_classical),
    ("09b_interactions_wp", p_09b_interactions_wp),
    ("09c_interactions_differing", p_09c_interactions_differing),
    ("09d_interactions_both", p_09d_interactions_both),
    ("09e_interactions_both_projectile", p_09e_interactions_both_projectile),
    ("10a_ps_pb_separate", p_10a_ps_pb_separate),
    ("10b_ps_pb_combined", p_10b_ps_pb_combined),
    ("11a_cumulative_impulse", p_11a_cumulative_impulse),
    ("11b_impulse_ratio", p_11b_impulse_ratio),
    ("11c_energy_gap", p_11c_energy_gap),
    ("12a_fit_targets_time", p_12a_fit_targets_time),
    ("12b_fit_targets_path", p_12b_fit_targets_path),
    ("12c_local_stopping_power", p_12c_local_stopping_power),
    ("13a_T1_window_fits", p_13a_T1_window_fits),
    ("13b_T2_window_fits", p_13b_T2_window_fits),
    ("13c_stopping_bar", p_13c_stopping_bar),
    ("13d_classical_fit_vs_path", p_13d_classical_fit_vs_path),
    ("13e_classical_energy_vs_time", p_13e_classical_energy_vs_time),
]


# ===========================================================================
# 2-D momentum maps  (their own driver: shared clim + linear/log variants)
# ===========================================================================

#: Transverse binning for the momentum maps. MEASURED, not guessed (2026-08-02).
#:
#: k_z needs no choice: it is the exact FFT grid, dk_z = 2 pi / L_z = 0.105
#: Bohr^-1, and that is a HARD limit of the 60 Bohr box — no analysis choice can
#: improve it. Only k_perp is binned, and the library default (one transverse
#: grid spacing, dk = 2 pi / 40 = 0.157) is 1.5x COARSER than the k_z axis, so
#: the map's resolution is transverse-limited.
#:
#: Counting the distinct |k_perp| = sqrt(kx^2+ky^2) lattice radii below 1.6
#: Bohr^-1 gives 45 of them from 325 grid points, so a finer binning is genuinely
#: supported: at dk_perp = 0.08 every bin is populated (minimum 1 point, ZERO
#: empty); at 0.05 four bins are already empty and the map combs. 0.08 — half a
#: grid spacing, double the default resolution — is therefore the finest honest
#: choice, and it is what is used here.
#:
#: KPERP_MAX is deliberately 2.0, ABOVE the 1.6 display limit, because
#: `kz_kperp_map` CLIPS all weight above `kperp_max` into the TOP bin. Cutting at
#: 1.6 would pile every high-k contribution into a fake spike at k_perp = 1.58,
#: inside the plotted range. Overflowing into [1.92, 2.0) puts that artefact
#: safely off-panel while keeping the bin edges aligned with the verified 0.08
#: grid. (Zero-padding the orbital to interpolate the k grid was REJECTED: 2.8 %
#: of the norm sits on a single box face at t = 30, so padding would impose a
#: discontinuity and manufacture ringing.)
KPERP_MAX = 2.0
KPERP_BINS = 25                       # 2.0 / 25 = 0.08 Bohr^-1
KPERP_DISPLAY = 1.6


def momentum_map_refined(C, step: int):
    """`C.maps[step]` recomputed at the finer transverse binning."""
    from inqview.visualisation.field_io import load_complex_vti, kz_kperp_map
    p = (Path(C.CS.WP_RESULTS) / C.wp_name / C.R.WF_SUBPATH
         / f"wavefunction_t{step:06d}.vti")
    return kz_kperp_map(load_complex_vti(p),
                        n_kperp_bins=KPERP_BINS, kperp_max=KPERP_MAX)


def momentum_map_clims(ctxs: dict) -> dict:
    """One clim for P and one for the difference, shared by BOTH run sets.

    Report-figures rule 7: any set of maps a reader compares must use ONE colour
    scale. Autoscaling per panel would make the twin and the SIC run look alike
    no matter how different they are.

    Restricted to the DISPLAYED window: the overflow bin at k_perp > 1.6 (see
    KPERP_MAX) and the far-k_z tails are never shown, so letting them set the
    colour scale would flatten the contrast on the structure that is.
    """
    pmax = 0.0
    dmax = 0.0
    for C in ctxs.values():
        kz, kp, P0 = C.maps[C.wf_steps[0]]
        vis = np.ix_((kz >= 0.8) & (kz <= 3.0), kp <= KPERP_DISPLAY)
        for s in C.wf_steps:
            P = C.maps[s][2]
            pmax = max(pmax, float(P[vis].max()))
            dmax = max(dmax, float(np.abs(P - P0)[vis].max()))
    return {"p_vmax": pmax, "d_vmax": dmax}


def draw_momentum_maps(C, clim, outdir: Path, dpi: int) -> list[str]:
    written = []
    P0 = C.maps[C.wf_steps[0]][2]
    pmax, dmax = clim["p_vmax"], clim["d_vmax"]
    pfloor = pmax / 1e4

    for s in C.wf_steps:
        kz, kp, P = C.maps[s]
        t_fs = s * C.CS.DT * AU_TO_FS
        tag = f"t{s * C.CS.DT:04.1f}au".replace(".", "p")

        kz_lo, kz_hi, kp_hi = 0.8, 3.0, KPERP_DISPLAY
        for scale in ("linear", "log"):
            fig, ax, cax = figure_map(kz_hi - kz_lo, kp_hi)
            if scale == "linear":
                im = ax.pcolormesh(kz, kp, P.T, shading="auto",
                                   cmap=style.cmap_for("sequential"),
                                   vmin=0, vmax=pmax)
            else:
                im = ax.pcolormesh(kz, kp, np.maximum(P.T, pfloor),
                                   shading="auto",
                                   cmap=style.cmap_for("sequential"),
                                   norm=LogNorm(vmin=pfloor, vmax=pmax))
            ax.axvline(C.CS.V0, color="w", ls="--", lw=0.9)
            ax.set_xlim(kz_lo, kz_hi)
            ax.set_ylim(0, kp_hi)
            ax.set_xlabel(style.axis_label("momentum", "$k_z$"))
            ax.set_ylabel(style.axis_label("momentum", r"$k_\perp$"))
            sci_colorbar(fig, im, cax, r"$P(k_z,k_\perp)$",
                         sci=(scale == "linear"))
            name = f"07a_P_map_{tag}_{scale}"
            fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches=None)
            plt.close(fig)
            written.append(name)

        if s == C.wf_steps[0]:
            continue                     # the difference from itself is zero
        D = P - P0
        fig, ax, cax = figure_map(kz_hi - kz_lo, kp_hi)
        lin = dmax / 100.0
        im = ax.pcolormesh(kz, kp, D.T, shading="auto",
                           cmap=style.cmap_for("diverging"),
                           norm=SymLogNorm(linthresh=lin, vmin=-dmax,
                                           vmax=dmax, base=10))
        ax.axvline(C.CS.V0, color="k", ls="--", lw=0.9)
        ax.set_xlim(kz_lo, kz_hi)
        ax.set_ylim(0, kp_hi)
        ax.set_xlabel(style.axis_label("momentum", "$k_z$"))
        ax.set_ylabel(style.axis_label("momentum", r"$k_\perp$"))
        sci_colorbar(fig, im, cax, r"$P(t)-P(0)$", sci=False,
                     ticks=symlog_ticks(dmax, lin))
        name = f"07b_P_diff_{tag}_symlog"
        fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches=None)
        plt.close(fig)
        written.append(name)
        _ = t_fs
    return written


# ===========================================================================
# ground-state density slices
# ===========================================================================

def draw_gs_density(outdir: Path, dpi: int) -> list[str]:
    """xz and yz slices of the ground-state bath density, linear + log.

    The VTI is loaded through the canonical `load_vti`, which returns PHYSICAL
    order and the cell-centred axes — never fftshifted
    (`.claude/rules/vti-coordinate-mapping.md`).
    """
    if not GS_VTI.is_file():
        print(f"  [skip] no ground-state VTI at {GS_VTI}")
        return []
    f = load_vti(GS_VTI)
    R_IN, R_OUT = 10.0, 14.0

    # Physical self-check: the annulus must actually be where the geometry says.
    # A centre<->edge index swap would put the maximum on the axis instead.
    mid_z = int(np.argmin(np.abs(f.z - 0.0)))
    prof = f.data[:, int(np.argmin(np.abs(f.y))), mid_z]
    r_at_max = abs(f.x[int(np.argmax(prof))])
    if not (R_IN - 1.0 <= r_at_max <= R_OUT + 1.0):
        raise AssertionError(
            f"ground-state density peaks at |x| = {r_at_max:.2f} Bohr, outside the "
            f"annulus [{R_IN}, {R_OUT}] — the VTI index->coordinate mapping is wrong")

    n = f.data
    vmax = float(np.percentile(n, 99.9))          # robust: kills single-cell spikes
    floor = vmax / 1e3
    written = []

    views = {
        "xz": (f.xz_slice(0.0), f.x, f.z, "$x$", "$z$"),
        "yz": (n[int(np.argmin(np.abs(f.x))), :, :].T, f.y, f.z, "$y$", "$z$"),
    }
    for view, (img, ha, va, hl, vl) in views.items():
        span_h = float(ha[-1] - ha[0])
        span_v = float(va[-1] - va[0])
        for scale in ("linear", "log"):
            fig, ax, cax = figure_map(span_h, span_v)
            extent = [ha[0], ha[-1], va[0], va[-1]]
            if scale == "linear":
                im = ax.imshow(img, origin="lower", extent=extent, aspect="equal",
                               cmap=style.cmap_for("sequential"),
                               vmin=0.0, vmax=vmax)
            else:
                im = ax.imshow(np.maximum(img, floor), origin="lower",
                               extent=extent, aspect="equal",
                               cmap=style.cmap_for("sequential"),
                               norm=LogNorm(vmin=floor, vmax=vmax))
            for r in (-R_OUT, -R_IN, R_IN, R_OUT):
                ax.axvline(r, color="w", lw=0.7, ls="--", alpha=0.8)
            ax.set_xlabel(style.axis_label("length", hl))
            ax.set_ylabel(style.axis_label("length", vl))
            sci_colorbar(fig, im, cax, r"$n$ (Bohr$^{-3}$)",
                         sci=(scale == "linear"))
            name = f"gs_density_{view}_{scale}"
            fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches=None)
            plt.close(fig)
            written.append(name)
    print(f"  ground state: peak |x| = {r_at_max:.2f} Bohr "
          f"(annulus {R_IN}-{R_OUT}), n_max(99.9%) = {vmax:.4e} Bohr^-3")
    return written


# ===========================================================================
# driver — two passes so the twin and sic sets share axis limits
# ===========================================================================

def collect_limits(ctxs: dict) -> dict:
    """Pass 1: draw every panel for every set, record the autoscaled limits."""
    lims: dict[str, dict] = {}
    for C in ctxs.values():
        for name, fn in LINE_PANELS:
            fig, ax = figure_line()
            try:
                fn(C, ax)
            except Exception as exc:                       # noqa: BLE001
                print(f"  [warn] {C.tag}/{name}: {exc}")
                plt.close(fig)
                continue
            xl, yl = ax.get_xlim(), ax.get_ylim()
            plt.close(fig)
            d = lims.setdefault(name, {"x": list(xl), "y": list(yl)})
            d["x"] = [min(d["x"][0], xl[0]), max(d["x"][1], xl[1])]
            d["y"] = [min(d["y"][0], yl[0]), max(d["y"][1], yl[1])]
    return lims


def _assert_legend_fits(fig, ax, where: str) -> None:
    """Fail loudly if the legend runs off the canvas or covers the y-axis label.

    This exists because a legend that overflows is NOT obviously broken in code:
    matplotlib clips it silently, so `09b_interactions_wp` shipped with three of
    its five entries cut off the left edge and nobody could tell which curve was
    which. A geometric check is the only thing that catches it without eyeballing
    every panel.
    """
    leg = ax.get_legend()
    if leg is None:
        return
    fig.canvas.draw()                      # legend extent is layout-dependent
    bb = leg.get_window_extent().transformed(fig.transFigure.inverted())
    if bb.x0 < -1e-3 or bb.x1 > 1 + 1e-3 or bb.y0 < -1e-3 or bb.y1 > 1 + 1e-3:
        raise AssertionError(
            f"{where}: legend extends outside the canvas "
            f"(x {bb.x0:.3f}..{bb.x1:.3f}, y {bb.y0:.3f}..{bb.y1:.3f}) — it will "
            f"be clipped. Shorten the labels or add YHEADROOM.")
    axbb = ax.get_position()
    if bb.x0 < axbb.x0 - 1e-3:
        raise AssertionError(
            f"{where}: legend starts left of the axes (x0 {bb.x0:.3f} < "
            f"{axbb.x0:.3f}) and will overlap the y-axis label.")


def draw_set(C, lims: dict, outdir: Path, dpi: int) -> list[str]:
    """Pass 2: redraw with the shared limits and save."""
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, fn in LINE_PANELS:
        fig, ax = figure_line()
        try:
            fn(C, ax)
        except Exception as exc:                           # noqa: BLE001
            print(f"  [FAIL] {C.tag}/{name}: {exc}")
            plt.close(fig)
            continue
        if name in lims:
            # a log y-axis must keep its own positive limits
            if ax.get_yscale() == "linear":
                lo, hi = lims[name]["y"]
                # reserve empty space at the top for an in-axes legend
                hi += YHEADROOM.get(name, 0.0) * (hi - lo)
                ax.set_ylim(lo, hi)
            ax.set_xlim(*lims[name]["x"])
        _assert_legend_fits(fig, ax, f"{C.tag}/{name}")
        fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches=None)
        plt.close(fig)
        written.append(name)
    return written


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["twin", "sic", "setup", "all"], default="all")
    ap.add_argument("--dpi", type=int, default=style.STYLE_CONFIG["save_dpi"])
    a = ap.parse_args()

    style.apply_theme()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, list[str]] = {}

    if a.only in ("setup", "all"):
        print("SETUP — ground-state density")
        d = OUT_ROOT / "setup"
        d.mkdir(parents=True, exist_ok=True)
        manifest["setup"] = draw_gs_density(d, a.dpi)

    tags = [t for t in ("twin", "sic") if a.only in ("all", t)]
    if tags:
        ctxs = {}
        for tag in tags:
            print(f"loading run set: {tag}")
            ctxs[tag] = Ctx(tag, RUN_SETS[tag])
            w = ctxs[tag].wp
            print(f"  {len(w)} steps, t = {w.t.iloc[0]:.2f}..{w.t.iloc[-1]:.2f} a.u.")

        print("pass 1: collecting shared axis limits across "
              f"{len(ctxs)} set(s)")
        lims = collect_limits(ctxs)
        clim = momentum_map_clims(ctxs)
        print(f"  shared momentum clim: P<={clim['p_vmax']:.3e}, "
              f"|dP|<={clim['d_vmax']:.3e}")

        for tag, C in ctxs.items():
            print(f"pass 2: writing {tag}")
            d = OUT_ROOT / tag
            names = draw_set(C, lims, d, a.dpi)
            names += draw_momentum_maps(C, clim, d, a.dpi)
            manifest[tag] = names
            print(f"  wrote {len(names)} figures to {d}")

        (OUT_ROOT / "shared_limits.json").write_text(
            json.dumps({"axis_limits": lims, "momentum_clim": clim}, indent=1))

    (OUT_ROOT / "manifest.json").write_text(json.dumps(manifest, indent=1))
    total = sum(len(v) for v in manifest.values())
    print(f"\n{total} figures under {OUT_ROOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
