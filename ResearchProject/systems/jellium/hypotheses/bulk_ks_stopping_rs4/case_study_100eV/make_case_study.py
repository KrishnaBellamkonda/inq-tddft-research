#!/usr/bin/env python3
"""Case study: 100 eV projectile in HIGH-DENSITY bulk jellium (r_s = 3.99 Bohr).

Builds the report-ready figure set comparing the CLASSICAL point projectile with
the KOHN-SHAM WAVEPACKET, for one matched twin pair, and writes the extracted
stopping powers with uncertainties to ``stopping_power.txt``.

WHAT IS PLOTTED, AND WHY EACH QUANTITY IS THE ONE IT IS
-------------------------------------------------------
The two halves are matched in every physical parameter (density, sigma, box,
grid, dt, functional) and differ ONLY in how the projectile is represented. That
makes the comparison meaningful, but it also means the two ledgers are NOT the
same object, so the analogue of "the projectile's kinetic energy" has to be
chosen deliberately:

  classical   T_cl = 1/2 m v_z^2                    (electron_track.csv)
  wavepacket  T_1  = <p^2>/2m                       total KS-orbital kinetic
              T_2  = <p>^2/2m                       DRIFT -- the classical analogue
              T_var= var(p)/2m = T_1 - T_2          momentum SPREAD, quantum-only

``T_2`` is the term to compare against the classical curve: it is the kinetic
energy of the centre-of-mass motion alone. ``T_1`` additionally contains the
packet's internal momentum spread, which at t=0 is the pure zero-point value
3/(4 sigma_psi^2) -- 5.10 eV at sigma_psi = 2 Bohr -- energy the packet has
before it has moved at all. Charging that to stopping would be wrong.
``T_var`` is a channel the classical projectile does not possess: it can absorb
drift energy without any of it reaching the bath.

POSITION. The classical z(t) is read directly. The wavepacket centroid is
obtained the way the user asked for and the way Ehrenfest's theorem licenses --
by integrating the mean momentum,

    z_c(t) = z_0 + (1/m) \\int_0^t <p_z>(t') dt'      (trapezoid, m = 1)

which is ``WPRun.s4``. This is preferred over the density centroid <z> because
it is immune to the periodic wrap at the cell face. The density centroid
(``WPRun.s3``, unwrapped in phase) is plotted alongside as an INDEPENDENT
cross-check: the two agree to ~0.08 Bohr over 68 Bohr of travel, which is the
numerical statement that Ehrenfest's theorem holds in this run.

STOPPING POWER. S = -dT/ds by OLS over the run's own transient-excluded window
[fit_t0, fit_t1], read from the run's wp_config.txt rather than retyped here.
The quoted uncertainty is the OLS slope standard error and a window-sensitivity
systematic (both window edges moved +/-3 a.u.) added in quadrature; the
systematic dominates and is the honest one, because where the transient ends is
a judgement call. See ``ks_stopping.fit_stopping``.

INTERACTION ENERGIES. The pairwise P/S/B decomposition (.claude/rules/
decomposed-interaction-energies.md). In BULK the background is uniform, so
phi_+ is identically zero and E_SB = E_PB = E_BB = 0; the physics is entirely in

    E_SS  bath-bath        (the wake the projectile drives)
    E_PP  projectile self-Hartree  (WP only in substance -- constant for a rigid
                                    classical cloud, decays as the packet spreads)
    E_PS  projectile-bath  (the interaction that does the stopping)

Style: the canonical project theme (ADR 0004, ``inqview.visualisation.style``),
600 dpi, fixed canvas, time axes in fs per the units standard with the a.u.
window stated in the annotation (the campaign quotes windows in a.u.).

Usage
-----
    venv/bin/python make_case_study.py                     # rs4 (high density)
    venv/bin/python make_case_study.py --family bulk_ks_stopping_rs4_sigma3
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
ENGINE = REPO / "ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping"
sys.path.insert(0, str(ENGINE))

import ks_stopping as K  # noqa: E402
from inqview.visualisation import style  # noqa: E402

HA_EV = K.HA_TO_EV
ATU_FS = 2.4188843265857e-2          # CODATA: 1 a.u. of time in femtoseconds

# One semantic colour per curve, held fixed across every figure so a reader can
# carry a colour between plots (the same discipline as the shared-clim rule).
C = {
    "classical": "#b2352c",
    "wp":        "#2f5d8f",
    "T1":        "#1a3f66",
    "T2":        "#4f97d0",
    "Tvar":      "#1f8a70",
    "e_ss":      "#2f5d8f",
    "e_pp":      "#b2352c",
    "e_ps":      "#1f8a70",
    # The three background terms are structurally zero in bulk; they get muted,
    # visually subordinate colours so they read as "present and flat" rather
    # than competing with the three terms that carry the physics.
    "e_sb":      "#8a6fb0",
    "e_pb":      "#c08a2e",
    "e_bb":      "#777777",
    "band":      "#999999",
    "fitline":   "#000000",
}


# ---------------------------------------------------------------------------
# Provenance — every constant is READ from the run, never retyped
# ---------------------------------------------------------------------------

def _parse_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    return out


class Meta:
    """Run parameters, parsed from the run's own artefacts."""

    def __init__(self, family: str, scripts: Path):
        self.family = family
        self.wp_dir = scripts / family / "wp"
        self.cl_dir = scripts / family / "classical"
        for d in (self.wp_dir, self.cl_dir):
            if not d.is_dir():
                raise FileNotFoundError(f"missing run half: {d}")

        rs_txt = _parse_kv(self.wp_dir / "results" / "run_summary.txt")
        cfg = _parse_kv(self.wp_dir / "results" / "raw" / "observables" / "wp_config.txt")

        cell = [float(x) for x in rs_txt["cell_bohr"].split("x")]
        self.lx, self.ly, self.lz = cell
        self.n_elec = int(rs_txt["n_electrons"])
        self.sigma_wp = float(rs_txt["wp_sigma_bohr"])
        self.energy_ev = float(rs_txt["wp_energy_ev"])
        self.k0 = float(rs_txt["wp_k0_bohr_inv"])
        self.dt = float(rs_txt["dt_au"])
        self.z0 = float(rs_txt["wp_center_bohr"].split()[2])
        self.fit_t0 = float(cfg["fit_t0_au"])
        self.fit_t1 = float(cfg["fit_t1_au"])

        volume = self.lx * self.ly * self.lz
        self.density = self.n_elec / volume
        self.r_s = (3.0 / (4.0 * np.pi * self.density)) ** (1.0 / 3.0)
        # Bulk plasma frequency of the interacting electron gas, omega_p = sqrt(4 pi n).
        self.omega_p_ev = np.sqrt(4.0 * np.pi * self.density) * HA_EV

    @property
    def subtitle(self) -> str:
        return (rf"$r_s$ = {self.r_s:.2f} Bohr, $\sigma$ = {self.sigma_wp:.0f} Bohr, "
                rf"{self.energy_ev:.0f} eV")


# ---------------------------------------------------------------------------
# Figure scaffolding — canonical theme, fixed canvas, room reserved for a title
# ---------------------------------------------------------------------------

# The theme's one-column rect is (0.180, 0.160, 0.785, 0.805), which runs to the
# very top of the canvas. These figures are standalone downloadables rather than
# slots in a LaTeX panel, so each carries a short identifying title; the rect is
# shortened to make room rather than falling back to bbox_inches="tight" (which
# would abandon the fixed on-page size the standard is built on).
#
# The reserved strip must clear the WHOLE two-line title, not just most of it.
# At 9 pt with linespacing 1.35 that is 2*9*1.35 = 24.3 pt of text plus the ~6 pt
# title pad = 0.42 in; on a 3.0 in canvas that is 0.14, and a further 0.03 keeps
# the ascenders off the boundary. Hence top of axes = 0.830. Verified by
# measuring the first inked row of every saved PNG -- an earlier 0.895 put the
# title flush against row 0 on 15 of 16 figures.
RECT = (0.180, 0.150, 0.785, 0.680)
RECT_WIDE = (0.098, 0.150, 0.888, 0.680)


def fig1(wide: bool = False):
    w = 6.8 if wide else 3.5
    fig = plt.figure(figsize=(w, 3.0))
    ax = fig.add_axes(RECT_WIDE if wide else RECT)
    return fig, ax


def fig2_stacked():
    """Two vertically stacked panels sharing a time axis (3.5 x 4.6 in)."""
    fig = plt.figure(figsize=(3.5, 4.6))
    # 0.42 in of title on a 4.6 in canvas is 0.091; top of axes = 0.870.
    # Left margin is WIDER than the one-column rect's 0.180: a stacked figure's
    # lower panel can carry signed decimal ticks ("-0.10"), which are far wider
    # than the "60" the 0.180 margin was sized for, and would push the y-label
    # off the canvas.
    top = fig.add_axes((0.215, 0.500, 0.750, 0.370))
    bot = fig.add_axes((0.215, 0.095, 0.750, 0.365), sharex=top)
    plt.setp(top.get_xticklabels(), visible=False)
    return fig, top, bot


def fig2_stacked_wide():
    """Two stacked panels at TWO-COLUMN width (6.8 x 4.2 in).

    The interaction panels carry SIX series each. At one-column width a
    six-entry legend is either half the panel tall (one column) or wider than
    the canvas (two columns) -- the margin check rejected both. Six series is a
    legitimate two-column figure; the extra width is what makes the legend fit
    beside the data instead of on top of it.
    """
    fig = plt.figure(figsize=(6.8, 4.2))
    top = fig.add_axes((0.115, 0.505, 0.860, 0.365))
    bot = fig.add_axes((0.115, 0.100, 0.860, 0.360), sharex=top)
    plt.setp(top.get_xticklabels(), visible=False)
    return fig, top, bot


def title(ax, main: str, meta: Meta) -> None:
    ax.set_title(f"{main}\n{meta.subtitle}", fontsize=9, linespacing=1.35)


def _shrink_overrunning_titles(fig) -> list[str]:
    """Reduce any title's font size until it fits inside the fixed canvas.

    A fixed-canvas figure (no bbox_inches="tight") silently CROPS a title that is
    wider than the figure -- the text simply runs off the edge and is lost at
    save time. Measuring the rendered extent and shrinking is what makes the
    fixed-size discipline safe for titles of unknown length, which matters
    because the subtitle carries run parameters that differ per family.
    """
    fig.canvas.draw()
    w_px = fig.get_size_inches()[0] * fig.dpi
    shrunk = []
    for ax in fig.axes:
        t = ax.title
        if not t.get_text():
            continue
        start = t.get_fontsize()
        for _ in range(20):
            bb = t.get_window_extent(fig.canvas.get_renderer())
            if bb.x0 >= 4.0 and bb.x1 <= w_px - 4.0:
                break
            t.set_fontsize(t.get_fontsize() - 0.25)
            fig.canvas.draw()
        if t.get_fontsize() < start:
            shrunk.append(f"{t.get_text().splitlines()[0]!r} "
                          f"{start:.2f}->{t.get_fontsize():.2f} pt")
    return shrunk


def legend(ax, loc: str, fontsize: float = 9, ncol: int | None = None):
    """Legend on a white patch, with real headroom made for it.

    Two separate problems, both visible in the first cut of these figures:
    the theme's frameless legend lets a curve run straight through the label
    text (unreadable), and an opaque patch placed over a curve HIDES data. The
    fix is to expand the axis limit on the legend's side so the legend sits in
    genuine empty space, and to back it with white as a second line of defence.
    """
    n = len(ax.get_legend_handles_labels()[0])
    # A six-entry legend is twice as tall as a three-entry one; two columns keep
    # it from eating half the panel.
    ncol = ncol if ncol is not None else (2 if n > 4 else 1)
    leg = ax.legend(loc=loc, fontsize=fontsize, frameon=True, facecolor="white",
                    edgecolor="none", framealpha=0.90, borderpad=0.3, ncol=ncol,
                    columnspacing=1.0, handlelength=1.6)
    if not (loc.startswith("upper") or loc.startswith("lower")):
        return leg

    # MEASURE the legend rather than assuming a fixed headroom. The first cut
    # used a flat 1.30x, which was sized for three entries and silently went
    # back to covering curves the moment the interaction panels grew to six.
    ax.figure.canvas.draw()
    frac = (leg.get_window_extent()
            .transformed(ax.transAxes.inverted()).height) + 0.05
    frac = min(frac, 0.55)                    # never compress data below ~45%
    lo, hi = ax.get_ylim()
    span = (hi - lo) / (1.0 - frac)
    if loc.startswith("upper"):
        ax.set_ylim(lo, lo + span)
    else:
        ax.set_ylim(hi - span, hi)
    return leg


def save(fig, path: Path) -> Path:
    for msg in _shrink_overrunning_titles(fig):
        print(f"    (title shrunk to fit: {msg})")
    fig.savefig(path, dpi=style.STYLE_CONFIG["save_dpi"], bbox_inches=None)
    plt.close(fig)
    print(f"  wrote {path.name}")
    return path


def shade_window(ax, meta: Meta, label: bool = True) -> None:
    """Mark the transient-excluded fit window on a TIME axis (in fs)."""
    ax.axvspan(meta.fit_t0 * ATU_FS, meta.fit_t1 * ATU_FS,
               color=C["band"], alpha=0.20, lw=0, zorder=0,
               label=(rf"fit window {meta.fit_t0:.1f}$-${meta.fit_t1:.1f} a.u."
                      if label else None))


def shade_window_s(ax, s: np.ndarray, t: np.ndarray, meta: Meta,
                   label: bool = True) -> None:
    """Mark the same window on a POSITION axis, via the trajectory s(t)."""
    m = (t >= meta.fit_t0) & (t <= meta.fit_t1)
    ax.axvspan(s[m][0], s[m][-1], color=C["band"], alpha=0.20, lw=0, zorder=0,
               label=(rf"fit window {meta.fit_t0:.1f}$-${meta.fit_t1:.1f} a.u."
                      if label else None))


def draw_fit(ax, fit: K.StoppingFit, colour: str = None) -> None:
    """Overlay the OLS line and its uncertainty band on a POSITION axis."""
    s, model = fit.s_fit, fit.T_model * HA_EV
    ax.plot(s, model, ls="--", lw=1.1, color=colour or C["fitline"], zorder=5)
    # Slope uncertainty pivots the line about the window centroid.
    delta = fit.uncertainty * np.abs(s - s.mean())
    ax.fill_between(s, model - delta, model + delta, color=colour or C["fitline"],
                    alpha=0.18, lw=0, zorder=4)


def annotate_S(ax, fit: K.StoppingFit, loc: str = "lower left") -> None:
    """Short quantitative annotation (permitted in clear whitespace, rule 2)."""
    txt = rf"$S$ = {fit.S_ev_per_bohr:.3f} $\pm$ {fit.uncertainty:.3f} eV/Bohr"
    xy = {"lower left": (0.04, 0.06), "lower right": (0.96, 0.06),
          "upper left": (0.04, 0.94), "upper right": (0.96, 0.94)}[loc]
    ha = "left" if "left" in loc else "right"
    va = "bottom" if "lower" in loc else "top"
    ax.text(*xy, txt, transform=ax.transAxes, fontsize=8, ha=ha, va=va,
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.82))


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def build(meta: Meta, out: Path) -> list[K.StoppingFit]:
    wp = K.load_wp_run(meta.wp_dir, meta.lz, meta.z0)
    cl = K.load_classical_run(meta.cl_dir, meta.lz)

    t_wp, t_cl = wp.t * ATU_FS, cl.t * ATU_FS
    T1, T2 = wp.T1 * HA_EV, wp.T2 * HA_EV
    Tvar = (wp.T1 - wp.T2) * HA_EV
    Tcl = cl.T * HA_EV
    zc, zc_dens, zcl = wp.s4, wp.s3, cl.z

    LT = style.axis_label("time")
    LZ = style.axis_label("length", "z")
    LE = style.axis_label("energy", "kinetic energy")

    L1 = r"$T_1 = \langle p^2 \rangle / 2m$"
    L2 = r"$T_2 = \langle p \rangle^2 / 2m$"
    LV = r"$T_\mathrm{var} = \mathrm{var}(p)/2m$"
    LCL = r"classical, $\frac{1}{2}mv_z^2$"

    # -- A. kinetic energy vs time -----------------------------------------
    fig, ax = fig1()
    ax.plot(t_cl, Tcl, color=C["classical"], label=LCL)
    ax.set_xlabel(LT); ax.set_ylabel(LE)
    title(ax, "Classical projectile: kinetic energy", meta)
    legend(ax, "upper right")
    save(fig, out / "01_classical_KE_vs_time.png")

    fig, ax = fig1()
    ax.plot(t_wp, T1, color=C["T1"], label=L1)
    ax.plot(t_wp, T2, color=C["T2"], ls="--", label=L2)
    ax.set_xlabel(LT); ax.set_ylabel(LE)
    title(ax, "Wavepacket: total vs drift kinetic energy", meta)
    legend(ax, "lower left")
    save(fig, out / "02_wp_T1_T2_vs_time.png")

    fig, ax = fig1()
    ax.plot(t_wp, Tvar, color=C["Tvar"], label=LV)
    ax.axhline(Tvar[0], color=C["band"], lw=0.8, ls=":")
    ax.text(0.97, 0.08, rf"$t=0$: $3/(4\sigma^2)$ = {Tvar[0]:.2f} eV",
            transform=ax.transAxes, fontsize=8, ha="right", va="bottom")
    ax.set_xlabel(LT); ax.set_ylabel(style.axis_label("energy", "spread energy"))
    title(ax, "Wavepacket: momentum-spread energy", meta)
    legend(ax, "upper left")
    save(fig, out / "03_wp_Tvar_vs_time.png")

    fig, top, bot = fig2_stacked()
    top.plot(t_wp, T1, color=C["T1"], label=L1)
    top.plot(t_wp, T2, color=C["T2"], ls="--", label=L2)
    top.plot(t_cl, Tcl, color=C["classical"], lw=1.0, alpha=0.85, label=LCL)
    top.set_ylabel(LE)
    legend(top, "lower left", 8)
    top.set_title("Kinetic-energy decomposition\n" + meta.subtitle,
                  fontsize=9, linespacing=1.35)
    bot.plot(t_wp, Tvar, color=C["Tvar"], label=LV)
    bot.set_xlabel(LT); bot.set_ylabel(style.axis_label("energy", "spread energy"))
    legend(bot, "upper left", 8)
    save(fig, out / "04_wp_kinetic_decomposition.png")

    # -- B. position vs time ------------------------------------------------
    fig, ax = fig1()
    ax.plot(t_cl, zcl, color=C["classical"], label="classical")
    ax.set_xlabel(LT); ax.set_ylabel(LZ)
    title(ax, "Classical projectile: position", meta)
    legend(ax, "upper left")
    save(fig, out / "05_classical_z_vs_time.png")

    # On one axes the two centroids are indistinguishable -- which IS the result,
    # but an invisible curve shows the reader nothing, and an inset placed over
    # the trajectory hides the very data it comments on. A second panel gives the
    # residual its own scale without occluding anything: the agreement becomes a
    # legible number of milli-Bohr instead of a claim in a caption.
    fig, top, bot = fig2_stacked()
    top.plot(t_wp, zc, color=C["wp"],
             label=r"$z_0 + \int \langle p_z \rangle\, dt\, / m$")
    top.plot(t_wp, zc_dens, color=C["Tvar"], ls=":", lw=1.4,
             label=r"density centroid $\langle z \rangle$")
    top.set_ylabel(LZ)
    top.set_title("Wavepacket centroid: Ehrenfest integration\n" + meta.subtitle,
                  fontsize=9, linespacing=1.35)
    legend(top, "upper left", 8)
    bot.plot(t_wp, wp.ehrenfest_residual, color=C["Tvar"],
             label=r"$\langle z \rangle - \int \langle p_z \rangle\, dt\, / m$")
    bot.axhline(0.0, color="k", lw=0.6, alpha=0.4)
    bot.set_xlabel(LT); bot.set_ylabel(style.axis_label("length", "residual"))
    legend(bot, "lower left", 8)
    save(fig, out / "06_wp_centroid_vs_time.png")

    fig, ax = fig1()
    ax.plot(t_cl, zcl, color=C["classical"], label="classical")
    ax.plot(t_wp, zc, color=C["wp"], label="wavepacket centroid")
    ax.set_xlabel(LT); ax.set_ylabel(LZ)
    title(ax, "Trajectories: classical vs wavepacket", meta)
    legend(ax, "upper left")
    save(fig, out / "07_position_vs_time_both.png")

    # -- C. kinetic energy vs position --------------------------------------
    fig, ax = fig1()
    ax.plot(zc, T1, color=C["T1"], label=L1)
    ax.plot(zc, T2, color=C["T2"], ls="--", label=L2)
    ax.set_xlabel(LZ); ax.set_ylabel(LE)
    title(ax, "Wavepacket: kinetic energy vs position", meta)
    legend(ax, "lower left")
    save(fig, out / "08_wp_T1_T2_vs_position.png")

    fig, ax = fig1()
    ax.plot(zcl, Tcl, color=C["classical"], label=LCL)
    ax.plot(zc, T1, color=C["T1"], label=L1)
    ax.plot(zc, T2, color=C["T2"], ls="--", label=L2)
    ax.set_xlabel(LZ); ax.set_ylabel(LE)
    title(ax, "Kinetic energy vs position: both projectiles", meta)
    legend(ax, "lower left", 8)
    save(fig, out / "09_T_vs_position_both.png")

    # -- D. the same plots, with the stopping-power fit overlaid ------------
    f_cl = K.fit_stopping(zcl, cl.T, cl.t, meta.fit_t0, meta.fit_t1,
                          "classical (1/2 m v^2)", v=cl.vz)
    f_T2 = K.fit_stopping(zc, wp.T2, wp.t, meta.fit_t0, meta.fit_t1,
                          "wavepacket drift T2", v=wp.pz)
    f_T1 = K.fit_stopping(zc, wp.T1, wp.t, meta.fit_t0, meta.fit_t1,
                          "wavepacket total T1", v=wp.pz)
    f_var = K.fit_stopping(zc, wp.T1 - wp.T2, wp.t, meta.fit_t0, meta.fit_t1,
                           "wavepacket spread T_var", v=wp.pz)

    fig, ax = fig1()
    ax.plot(zcl, Tcl, color=C["classical"], label=LCL)
    shade_window_s(ax, zcl, cl.t, meta)
    draw_fit(ax, f_cl)
    annotate_S(ax, f_cl, "lower left")
    ax.set_xlabel(LZ); ax.set_ylabel(LE)
    title(ax, "Stopping-power fit: classical", meta)
    legend(ax, "upper right", 8)
    save(fig, out / "10_fit_classical_T_vs_position.png")

    fig, ax = fig1()
    ax.plot(zc, T1, color=C["T1"], label=L1)
    ax.plot(zc, T2, color=C["T2"], ls="--", label=L2)
    shade_window_s(ax, zc, wp.t, meta)
    draw_fit(ax, f_T1, C["T1"])
    draw_fit(ax, f_T2, C["T2"])
    annotate_S(ax, f_T2, "lower left")
    ax.set_xlabel(LZ); ax.set_ylabel(LE)
    title(ax, r"Stopping-power fit: wavepacket", meta)
    legend(ax, "upper right", 8)
    save(fig, out / "11_fit_wp_T1_T2_vs_position.png")

    fig, ax = fig1()
    ax.plot(t_cl, Tcl, color=C["classical"], label=LCL)
    shade_window(ax, meta)
    m = (cl.t >= meta.fit_t0) & (cl.t <= meta.fit_t1)
    ax.plot(t_cl[m], f_cl.T_model * HA_EV, ls="--", lw=1.1, color=C["fitline"])
    annotate_S(ax, f_cl, "lower left")
    ax.set_xlabel(LT); ax.set_ylabel(LE)
    title(ax, "Stopping-power window: classical", meta)
    legend(ax, "upper right", 8)
    save(fig, out / "12_fit_classical_KE_vs_time.png")

    fig, ax = fig1()
    ax.plot(t_wp, T1, color=C["T1"], label=L1)
    ax.plot(t_wp, T2, color=C["T2"], ls="--", label=L2)
    shade_window(ax, meta)
    mw = (wp.t >= meta.fit_t0) & (wp.t <= meta.fit_t1)
    ax.plot(t_wp[mw], f_T1.T_model * HA_EV, ls="--", lw=1.1, color=C["fitline"])
    ax.plot(t_wp[mw], f_T2.T_model * HA_EV, ls="--", lw=1.1, color=C["fitline"])
    annotate_S(ax, f_T2, "lower left")
    ax.set_xlabel(LT); ax.set_ylabel(LE)
    title(ax, r"Stopping-power window: wavepacket", meta)
    legend(ax, "upper right", 8)
    save(fig, out / "13_fit_wp_kinetic_vs_time.png")

    # -- E. interaction energies -------------------------------------------
    ix_cl = K.load_interactions(meta.cl_dir, "classical")
    ix_wp = K.load_interactions(meta.wp_dir, "wp")

    # ALL SIX pairwise terms are plotted, not just the three that move.
    #
    # In BULK the background is uniform, so poisson(n_+) is pure G=0 -- which INQ
    # drops -- and phi_+ is identically zero, forcing E_SB = E_PB = E_BB = 0.
    # They are written as columns anyway so the schema matches the slab systems
    # (.claude/rules/decomposed-interaction-energies.md). Showing them flat at
    # zero on the SAME axis as the terms that do move is what distinguishes
    # "structurally zero" from "we forgot to compute them", and the magnified
    # lower panel turns that from an assertion into something the reader can
    # check: on a scale 50x finer than the E_SS excursion they still do not move.
    # Verified 2026-08-02: all three are bitwise 0.0 in both halves.
    # (column, colour, linestyle, linewidth, panel label, difference label).
    # The difference label is spelled out rather than derived from the panel
    # label by string surgery -- splitting on "$" leaves a bare unbalanced "$"
    # literal in the source, which is indistinguishable from a broken mathtext
    # span to the guard that checks them.
    ALL_TERMS = (
        ("e_ss", C["e_ss"], "-",  1.4, r"$E_\mathrm{SS}$ (bath$-$bath)",
         r"$\Delta E_\mathrm{SS}$"),
        ("e_pp", C["e_pp"], "--", 1.4, r"$E_\mathrm{PP}$ (projectile self)",
         r"$\Delta E_\mathrm{PP}$"),
        ("e_ps", C["e_ps"], "-.", 1.4, r"$E_\mathrm{PS}$ (projectile$-$bath)",
         r"$\Delta E_\mathrm{PS}$"),
        ("e_sb", C["e_sb"], "-",  0.9, r"$E_\mathrm{SB}$ (bath$-$background)",
         r"$\Delta E_\mathrm{SB}$"),
        ("e_pb", C["e_pb"], "--", 0.9, r"$E_\mathrm{PB}$ (projectile$-$bkg)",
         r"$\Delta E_\mathrm{PB}$"),
        ("e_bb", C["e_bb"], ":",  0.9, r"$E_\mathrm{BB}$ (background self)",
         r"$\Delta E_\mathrm{BB}$"),
    )
    BKG = ("e_sb", "e_pb", "e_bb")

    def _ix_panel(ix, name, fname):
        fig, top, bot = fig2_stacked_wide()
        t_fs = ix.t * ATU_FS
        for key, col, ls, lw, lab, _ in ALL_TERMS:
            top.plot(t_fs, getattr(ix, key) * HA_EV, color=col, ls=ls, lw=lw,
                     label=lab)
        top.axhline(0.0, color="k", lw=0.6, alpha=0.4)
        shade_window(top, meta, label=False)
        top.set_ylabel(style.axis_label("energy"))
        top.set_title(f"Interaction energies: {name}\n" + meta.subtitle,
                      fontsize=9, linespacing=1.35)
        legend(top, "upper left", 8, ncol=3)

        span = max(np.ptp(getattr(ix, k)) for k, *_rest in ALL_TERMS) * HA_EV
        worst = max(np.max(np.abs(getattr(ix, k))) for k in BKG) * HA_EV
        for key, col, ls, lw, lab, _ in ALL_TERMS:
            if key in BKG:
                bot.plot(t_fs, getattr(ix, key) * HA_EV, color=col, ls=ls,
                         lw=1.2, label=lab)
        bot.axhline(0.0, color="k", lw=0.6, alpha=0.4)
        shade_window(bot, meta, label=False)
        # Magnified 50x against the largest excursion in the panel above.
        bot.set_ylim(-0.02 * span, 0.02 * span)
        bot.set_xlabel(LT)
        bot.set_ylabel(style.axis_label("energy", "background terms"))
        bot.text(0.5, 0.86, rf"max $|E|$ = {worst:.1e} eV  (bulk: $\varphi_+ \equiv 0$)",
                 transform=bot.transAxes, fontsize=7, ha="center", va="top")
        legend(bot, "lower center", 8, ncol=3)
        return save(fig, out / fname)

    _ix_panel(ix_cl, "classical", "14_interactions_classical.png")
    _ix_panel(ix_wp, "wavepacket", "15_interactions_wp.png")

    # Difference: merge on step so the two cadences cannot silently misalign.
    keys = [k for k, *_rest in ALL_TERMS]
    d = pd.merge(
        pd.DataFrame({"step": ix_wp.step, "t": ix_wp.t,
                      **{k: getattr(ix_wp, k) for k in keys}}),
        pd.DataFrame({"step": ix_cl.step,
                      **{k: getattr(ix_cl, k) for k in keys}}),
        on="step", suffixes=("_wp", "_cl"))
    if d.empty:
        raise RuntimeError("no shared steps between the two halves' interactions.csv")

    fig, ax = fig1(wide=True)
    for key, col, ls, lw, _lab, dlab in ALL_TERMS:
        ax.plot(d["t"] * ATU_FS, (d[f"{key}_wp"] - d[f"{key}_cl"]) * HA_EV,
                color=col, ls=ls, lw=lw, label=dlab)
    ax.axhline(0.0, color="k", lw=0.6, alpha=0.4)
    shade_window(ax, meta, label=False)
    ax.set_xlabel(LT)
    ax.set_ylabel(style.axis_label("energy", "wavepacket $-$ classical"))
    title(ax, "Interaction-energy difference (all terms)", meta)
    legend(ax, "lower left", 8, ncol=3)
    save(fig, out / "16_interactions_difference.png")

    return [f_cl, f_T2, f_T1, f_var], (wp, cl, ix_wp, ix_cl, d)


# ---------------------------------------------------------------------------
# The results file
# ---------------------------------------------------------------------------

def write_results(meta: Meta, fits, data, out: Path) -> None:
    f_cl, f_T2, f_T1, f_var = fits
    wp, cl, ix_wp, ix_cl, d = data
    ratio = f_cl.S_ev_per_bohr / f_T2.S_ev_per_bohr
    # Ratio uncertainty by standard propagation of two independent fits.
    ratio_u = abs(ratio) * np.hypot(f_cl.uncertainty / f_cl.S_ev_per_bohr,
                                    f_T2.uncertainty / f_T2.S_ev_per_bohr)

    def dE(arr):
        return (arr[-1] - arr[0]) * HA_EV

    lines = [
        "ELECTRONIC STOPPING POWER -- 100 eV projectile in high-density bulk jellium",
        "=" * 78,
        "",
        f"family        : {meta.family}",
        f"r_s           : {meta.r_s:.3f} Bohr   (n = {meta.density:.6f} a.u.^-3, "
        f"{meta.n_elec} electrons in {meta.lx:.0f} x {meta.ly:.0f} x {meta.lz:.0f} Bohr)",
        f"omega_p       : {meta.omega_p_ev:.2f} eV",
        f"sigma_WP      : {meta.sigma_wp:.0f} Bohr   (classical UPF charge std "
        f"= sigma_WP/sqrt(2) = {meta.sigma_wp / np.sqrt(2):.3f} Bohr)",
        f"projectile KE : {meta.energy_ev:.0f} eV   (k0 = {meta.k0:.4f} Bohr^-1, "
        f"v0 = {meta.k0:.4f} a.u.)",
        f"dt            : {meta.dt} a.u.",
        f"fit window    : {meta.fit_t0:.2f} - {meta.fit_t1:.2f} a.u.  "
        f"({meta.fit_t0 * ATU_FS:.4f} - {meta.fit_t1 * ATU_FS:.4f} fs)",
        "",
        "METHOD",
        "-" * 78,
        "S = -dT/ds, ordinary least squares over the transient-excluded window",
        "above, with s the projectile path coordinate:",
        "  classical    s = z(t) from electron_track.csv",
        "  wavepacket   s = z_0 + integral of <p_z> dt / m  (Ehrenfest centroid)",
        "",
        "Quoted uncertainty = OLS slope standard error (stat) and a",
        "window-sensitivity systematic (syst) in quadrature. The systematic moves",
        "BOTH window edges independently by +/-3 a.u. and takes half the full",
        "spread of the resulting slopes; it prices the judgement call about where",
        "the transient ends, and it dominates the statistical error.",
        "",
        "RESULTS",
        "-" * 78,
    ]
    for f in fits:
        lines.append("  " + f.summary())
    lines += [
        "",
        f"  classical / wavepacket(T2) ratio = {ratio:.2f} +/- {ratio_u:.2f}",
        "",
        "READING THESE",
        "-" * 78,
        "T_2 is the comparable quantity: it is the drift kinetic energy of the",
        "packet's centre of mass, the direct analogue of 1/2 m v^2. T_1 also",
        "contains the internal momentum spread, which at t = 0 is the pure",
        f"zero-point value 3/(4 sigma^2) = {(wp.T1[0] - wp.T2[0]) * HA_EV:.3f} eV --",
        "energy the packet carries before it has moved. The T_var row is NOT a",
        "stopping power in the usual sense: it is the rate at which drift energy",
        "is converted into the packet's own momentum spread, a channel with no",
        "classical counterpart. A negative S there means the spread is GROWING.",
        "",
        "ENERGY BUDGET OVER THE WHOLE RUN (t = 0 to "
        f"{wp.t[-1]:.2f} a.u.)",
        "-" * 78,
        f"  classical   1/2 m v^2 : {cl.T[0] * HA_EV:8.3f} -> {cl.T[-1] * HA_EV:8.3f} eV"
        f"   (change {dE(cl.T):+8.3f} eV)",
        f"  wavepacket  T_1       : {wp.T1[0] * HA_EV:8.3f} -> {wp.T1[-1] * HA_EV:8.3f} eV"
        f"   (change {dE(wp.T1):+8.3f} eV)",
        f"  wavepacket  T_2       : {wp.T2[0] * HA_EV:8.3f} -> {wp.T2[-1] * HA_EV:8.3f} eV"
        f"   (change {dE(wp.T2):+8.3f} eV)",
        f"  wavepacket  T_var     : {(wp.T1[0] - wp.T2[0]) * HA_EV:8.3f} -> "
        f"{(wp.T1[-1] - wp.T2[-1]) * HA_EV:8.3f} eV   (change {dE(wp.T1 - wp.T2):+8.3f} eV)",
        "",
        "INTERACTION ENERGIES (change over the run, eV)",
        "-" * 78,
        f"  {'term':6s} {'classical':>12s} {'wavepacket':>12s} {'WP - cl':>12s}",
    ]
    for key, arr_cl, arr_wp in (("E_SS", ix_cl.e_ss, ix_wp.e_ss),
                                ("E_PP", ix_cl.e_pp, ix_wp.e_pp),
                                ("E_PS", ix_cl.e_ps, ix_wp.e_ps),
                                ("E_SB", ix_cl.e_sb, ix_wp.e_sb),
                                ("E_PB", ix_cl.e_pb, ix_wp.e_pb),
                                ("E_BB", ix_cl.e_bb, ix_wp.e_bb)):
        lines.append(f"  {key:6s} {dE(arr_cl):12.3f} {dE(arr_wp):12.3f} "
                     f"{dE(arr_wp) - dE(arr_cl):12.3f}")
    lines += [
        "",
        "  All SIX pairwise terms are listed. E_SB / E_PB / E_BB are bitwise zero",
        "  at every step in both halves -- verified, not assumed:",
        f"    classical  max|E_SB|,|E_PB|,|E_BB| = "
        f"{max(np.max(np.abs(a)) for a in (ix_cl.e_sb, ix_cl.e_pb, ix_cl.e_bb)):.1e} Ha",
        f"    wavepacket max|E_SB|,|E_PB|,|E_BB| = "
        f"{max(np.max(np.abs(a)) for a in (ix_wp.e_sb, ix_wp.e_pb, ix_wp.e_bb)):.1e} Ha",
    ]
    lines += [
        "",
        f"  E_PP(t=0) classical = {ix_cl.e_pp[0]:.6f} Ha, "
        f"wavepacket = {ix_wp.e_pp[0]:.6f} Ha",
        "  These agree because the classical Gaussian UPF is generated at",
        "  sigma_pot = sigma_WP/sqrt(2), so its charge cloud has the wavepacket's",
        "  t=0 density. That equality is a check on the sigma-matching convention",
        "  (.claude/rules/sigma-wp-convention.md), not an input to it.",
        "",
        "  Bulk jellium has a UNIFORM background, so phi_+ is identically zero and",
        "  E_SB = E_PB = E_BB = 0 by construction. Absolute E_PP carries the",
        "  charged-cell G=0 gauge; only the WP-minus-classical differences above",
        "  are gauge-clean.",
        "",
        "VALIDATION",
        "-" * 78,
        f"  Hartree closure vs INQ, classical : {np.nanmax(np.abs(ix_cl.closure)):.2e} Ha",
        f"  Hartree closure vs INQ, wavepacket: {np.nanmax(np.abs(ix_wp.closure)):.2e} Ha",
        f"  Ehrenfest residual max |<z> - int<p>dt| : "
        f"{np.max(np.abs(wp.ehrenfest_residual)):.4f} Bohr over "
        f"{abs(wp.s4[-1] - wp.s4[0]):.1f} Bohr of travel",
        f"  WP orbital norm range : {wp.norm.min():.9f} to {wp.norm.max():.9f}",
        f"  classical cloud clipping onset : "
        + ("never" if not np.isfinite(ix_cl.clip_time) else f"{ix_cl.clip_time:.2f} a.u.")
        + f"  (fit window ends {meta.fit_t1:.2f} a.u.)",
        f"  shared steps in the difference plot : {len(d)}",
        "",
        "PROVENANCE",
        "-" * 78,
        f"  wp run        : {meta.wp_dir}",
        f"  classical run : {meta.cl_dir}",
        f"  engine        : {ENGINE / 'ks_stopping.py'}",
        f"  builder       : {Path(__file__).resolve()}",
        "",
    ]
    path = out / "stopping_power.txt"
    path.write_text("\n".join(lines))
    print(f"  wrote {path.name}")
    print()
    for f in fits:
        print("   ", f.summary())
    print(f"    ratio classical/WP(T2) = {ratio:.2f} +/- {ratio_u:.2f}")


MIN_MARGIN_PX = 8


def verify_margins(out: Path) -> int:
    """Fail the build if any saved figure's ink runs into the canvas edge.

    A fixed-canvas figure crops silently: a title or y-label that overruns is
    simply gone from the PNG, and nothing in the build says so. This measures the
    ink bounding box of every figure actually written and reports the offenders,
    which is the only check that sees what a reader will see. It caught two
    distinct clipping bugs while this script was being written (a two-line title
    flush against row 0 on 15 of 16 figures, and a stacked panel's signed decimal
    ticks pushing the y-label off the left edge).
    """
    try:
        from PIL import Image
    except ImportError:
        print("  (margin check skipped: Pillow not available)")
        return 0

    bad = []
    for f in sorted(out.glob("*.png")):
        a = np.asarray(Image.open(f).convert("L"))
        rows = np.where(a.min(axis=1) < 200)[0]
        cols = np.where(a.min(axis=0) < 200)[0]
        if not len(rows) or not len(cols):
            bad.append((f.name, "blank")); continue
        m = min(rows[0], a.shape[0] - 1 - rows[-1],
                cols[0], a.shape[1] - 1 - cols[-1])
        if m < MIN_MARGIN_PX:
            bad.append((f.name, f"{m} px from an edge"))
    if bad:
        print(f"\n  FAIL: {len(bad)} figure(s) clipped at the canvas edge:")
        for name, why in bad:
            print(f"    {name}: {why}")
        return 1
    print(f"  margin check: all {len(list(out.glob('*.png')))} figures clear "
          f"(>= {MIN_MARGIN_PX} px)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--family", default="bulk_ks_stopping_rs4")
    ap.add_argument("--out", default=None,
                    help="output directory (default: this script's directory)")
    a = ap.parse_args()

    scripts = REPO / "ResearchProject/systems/jellium/scripts"
    meta = Meta(a.family, scripts)
    out = Path(a.out) if a.out else HERE
    out.mkdir(parents=True, exist_ok=True)

    style.apply_theme()
    print(f"case study: {meta.family}  ({meta.subtitle})")
    print(f"output    : {out}")
    fits, data = build(meta, out)
    write_results(meta, fits, data, out)
    return verify_margins(out)


if __name__ == "__main__":
    raise SystemExit(main())
