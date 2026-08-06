#!/usr/bin/env python3
"""Report-ready figures for the VACUUM self-interaction investigation.

    report_figures/sigma4/   the single sigma = 4 experiment (5 theory levels)
    report_figures/sweep/    how the error scales with sigma_WP (6 sigma x 4 levels)

One PNG per panel, never pre-composed — composition is a downstream LaTeX
concern. Same house standard as the channeling set: canonical
`inqview.visualisation.style`, no on-canvas titles, fixed canvas
(`bbox_inches=None`), legends geometrically checked so none can ship clipped.

WHAT THIS INVESTIGATION IS
--------------------------
A single electron has, EXACTLY, no self-interaction, so one electron alone in a
vacuum box must follow closed-form free dispersion. Propagating the SAME injected
Gaussian at several theory levels therefore measures the self-interaction error by
DIFFERENCE, against a reference that is exact:

    noninteracting   nothing                        <- the reference
    hartree          + Hartree self-interaction
    lda              + Hartree AND LDA xc self-interaction
    sic_h            lda MINUS the Hartree self-term only
    sic_pzrun        lda MINUS Hartree AND xc self-terms (full Perdew-Zunger)

Ratios are taken against the REFERENCE RUN rather than the analytic formula, so
grid and propagator error — identical between them by construction — cancels.

TIME AXIS. Unlike the channeling set these runs are stationary (k0 = 0) and the
sweep's duration scales as sigma^2, so a shared fs axis would compress five of the
six sweep panels into nothing. The sigma = 4 panels use fs (with the a.u. twin
axis, as the house standard requires); the SWEEP panels use the dimensionless
tau = t/sigma^2, which is the only axis on which the six runs are comparable at
all — stated on the axis label so it cannot be mistaken for a time in disguise.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from inqview.visualisation import style        # noqa: E402

import selfinteraction as SI                   # noqa: E402
import sigma_sweep as SW                       # noqa: E402

OUT_ROOT = HERE / "report_figures"
AU_TO_FS = 0.024188843265857

#: Hue = theory level, one meaning only. The reference is grey because it is the
#: baseline every other curve is measured against, not a competing result.
COLOR = {
    "noninteracting": "0.45",
    "hartree": "tab:orange",
    "lda": "tab:red",
    "sic_h": "tab:purple",
    "sic_pzrun": "tab:green",
}
LABEL = {
    "noninteracting": "non-interacting",
    "hartree": "Hartree",
    "lda": "Hartree + LDA xc",
    "sic_h": "SIC-H",
    "sic_pzrun": "SIC-PZ",
}
THREE = ("noninteracting", "hartree", "lda")
FIVE = ("noninteracting", "hartree", "lda", "sic_h", "sic_pzrun")

#: The closed-form results the whole investigation rests on, rendered on-panel.
F_SIGMA = r"$\sigma_\mathrm{dens}(t)=\sqrt{\sigma^2/2+t^2/2\sigma^2}$"
F_VARP = r"$\mathrm{var}(p_d)=1/2\sigma^2$ (constant)"
F_EPP = r"$E_{PP}=\left[\frac{1}{\sqrt{2\pi}}-\frac{\xi}{36}\right]\sigma^{-1}$"

_LINE_W, _LINE_H = 3.5, 3.3
_M = dict(left=0.82, right=0.10, bottom=0.50, top=0.46)


def figure_line():
    aw = _LINE_W - _M["left"] - _M["right"]
    ah = _LINE_H - _M["bottom"] - _M["top"]
    fig = plt.figure(figsize=(_LINE_W, _LINE_H))
    ax = fig.add_axes((_M["left"] / _LINE_W, _M["bottom"] / _LINE_H,
                       aw / _LINE_W, ah / _LINE_H))
    return fig, ax


def t_axis_fs(ax):
    ax.set_xlabel(style.axis_label("time", "$t$"))
    ax.tick_params(top=False)
    sec = ax.secondary_xaxis("top", functions=(lambda v: v / AU_TO_FS,
                                               lambda v: v * AU_TO_FS))
    sec.set_xlabel("$t$ (a.u.)", labelpad=4)


def legend(ax, **kw):
    kw.setdefault("frameon", True); kw.setdefault("framealpha", 0.88)
    kw.setdefault("facecolor", "white"); kw.setdefault("edgecolor", "none")
    kw.setdefault("borderpad", 0.25); kw.setdefault("handlelength", 1.5)
    kw.setdefault("fontsize", 8)
    return ax.legend(**kw)


def note(ax, text, loc="upper left"):
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    y, va = (0.04, "bottom") if "lower" in loc else (0.96, "top")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))


def check_legend(fig, ax, where):
    leg = ax.get_legend()
    if leg is None:
        return
    fig.canvas.draw()
    bb = leg.get_window_extent().transformed(fig.transFigure.inverted())
    if bb.x0 < -1e-3 or bb.x1 > 1 + 1e-3 or bb.y0 < -1e-3 or bb.y1 > 1 + 1e-3:
        raise AssertionError(f"{where}: legend clipped ({bb.x0:.3f}..{bb.x1:.3f})")


def save(fig, ax, outdir: Path, name: str, dpi: int, headroom: float = 0.0):
    if headroom and ax.get_yscale() == "linear":
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + headroom * (hi - lo))
    check_legend(fig, ax, name)
    fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches=None)
    plt.close(fig)
    return name


# ===========================================================================
# sigma = 4 — the single experiment
# ===========================================================================

def build_sigma4(outdir: Path, dpi: int) -> list[str]:
    outdir.mkdir(parents=True, exist_ok=True)
    runs = SI.load_all()
    ref = runs["noninteracting"]
    have = [t for t in FIVE if t in runs]
    w = []

    # --- 01 the measured Gaussian width, three levels (the headline)
    fig, ax = figure_line()
    ax.plot(ref.t * AU_TO_FS, ref.sigma_free(), lw=2.6, color="0.82", zorder=0,
            label="analytic free")
    for th in THREE:
        r = runs[th]
        ax.plot(r.t * AU_TO_FS, r.sigma_iso, lw=1.4, color=COLOR[th], label=LABEL[th])
    ax.set_ylabel(r"$\sigma_\mathrm{dens}$ (Bohr)")
    legend(ax, loc="upper left")
    note(ax, F_SIGMA, "lower right")
    t_axis_fs(ax)
    w.append(save(fig, ax, outdir, "01_width_three_levels", dpi, 0.06))

    # --- 01b absolute excess: the LDA curve is unreadable as a ratio near 1.08
    fig, ax = figure_line()
    for th in THREE:
        r = runs[th]
        n = min(len(r.t), len(ref.t))
        ax.plot(r.t[:n] * AU_TO_FS, r.sigma_iso[:n] - r.sigma_free()[:n], lw=1.4,
                color=COLOR[th], label=LABEL[th])
    ax.axhline(0.0, color="0.5", ls="--", lw=0.9)
    ax.set_ylabel(r"$\sigma_\mathrm{meas}-\sigma_\mathrm{free}$ (Bohr)")
    legend(ax, loc="upper left")
    t_axis_fs(ax)
    w.append(save(fig, ax, outdir, "01b_width_absolute_excess", dpi, 0.10))

    # --- 02 all five levels as a ratio: the SIC intervention test
    fig, ax = figure_line()
    for th in have:
        r = runs[th]
        n = min(len(r.t), len(ref.t))
        ax.plot(r.t[:n] * AU_TO_FS, r.sigma_iso[:n] / ref.sigma_iso[:n], lw=1.4,
                color=COLOR[th], label=LABEL[th])
    ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
    ax.set_ylabel(r"$\sigma\,/\,\sigma_\mathrm{non-int}$")
    legend(ax, loc="upper left", ncol=2)
    note(ax, "1.0 = no self-interaction", "lower right")
    t_axis_fs(ax)
    w.append(save(fig, ax, outdir, "02_excess_ratio_five_levels", dpi, 0.22))

    # --- 03 var(p): the SHARPEST gate — exactly conserved under free evolution
    fig, ax = figure_line()
    for th in have:
        r = runs[th]
        ax.plot(r.t * AU_TO_FS, r.var_p3d / (3.0 * SI.var_p_free(r.sigma_wp)),
                lw=1.4, color=COLOR[th], label=LABEL[th])
    ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
    ax.set_ylabel(r"$\mathrm{var}(p)\,/\,$ free value")
    legend(ax, loc="upper left", ncol=2)
    note(ax, F_VARP, "lower right")
    t_axis_fs(ax)
    w.append(save(fig, ax, outdir, "03_var_p_gate", dpi, 0.22))

    # --- 04 the energy: E_PP is a SOURCE, not a store
    fig, ax = figure_line()
    for th in THREE:
        r = runs[th]
        ax.plot(r.t * AU_TO_FS, r.e_pp_ev, lw=1.4, color=COLOR[th], label=LABEL[th])
    ax.set_ylabel(r"$E_{PP}$ (eV)")
    legend(ax, loc="upper right")
    note(ax, r"$E_{PP}=\frac{1}{2}\int n_\mathrm{wp}\phi_\mathrm{wp}$", "lower left")
    t_axis_fs(ax)
    w.append(save(fig, ax, outdir, "04_self_hartree_energy", dpi, 0.08))

    # --- 05 the closed identity: E_PP released == excess var(p)/2m gained
    fig, ax = figure_line()
    for th in ("hartree", "lda"):
        if th not in runs:
            continue
        e = SI.effect(runs[th], ref)
        ax.plot(e.t * AU_TO_FS, -e.d_e_pp_ev, lw=1.4, ls="--", color=COLOR[th],
                label=f"{LABEL[th]}: $-\\Delta E_{{PP}}$")
        ax.plot(e.t * AU_TO_FS, e.d_var_term_ev, lw=1.4, color=COLOR[th],
                label=f"{LABEL[th]}: excess $\\mathrm{{var}}(p)/2m$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_ylabel(style.axis_label("energy", r"$\Delta E$"))
    legend(ax, loc="upper left")
    note(ax, r"$k_0=0:\;\Delta[\mathrm{var}(p)/2m]=-\Delta E_{PP}$", "lower right")
    t_axis_fs(ax)
    w.append(save(fig, ax, outdir, "05_energy_identity", dpi, 0.30))
    return w


# ===========================================================================
# the sigma sweep
# ===========================================================================

def build_sweep(outdir: Path, dpi: int) -> list[str]:
    outdir.mkdir(parents=True, exist_ok=True)
    runs = SW.load_all()
    if not runs:
        print("  [skip] no sweep runs found")
        return []
    tab = SW.sigma_table(runs)
    gate = SW.protocol_gate(runs)
    w = []

    # --- 10 excess vs sigma: the trend, both interacting levels + the correction
    fig, ax = figure_line()
    ax.plot(tab.sigma, tab.excess_hartree - 1.0, "o-", lw=1.4, ms=4,
            color=COLOR["hartree"], label=LABEL["hartree"])
    ax.plot(tab.sigma, tab.excess_lda - 1.0, "s-", lw=1.4, ms=4,
            color=COLOR["lda"], label=LABEL["lda"])
    ax.plot(tab.sigma, tab.excess_sic_pzrun - 1.0, "^-", lw=1.4, ms=4,
            color=COLOR["sic_pzrun"], label=LABEL["sic_pzrun"])
    ax.axhline(0.0, color="0.4", ls="--", lw=1.0)
    ax.set_xlabel(r"$\sigma_\mathrm{WP}$ (Bohr)")
    ax.set_ylabel("excess spreading")
    legend(ax, loc="upper left")
    note(ax, r"at fixed $\tau=t/\sigma^2=1.875$", "lower right")
    w.append(save(fig, ax, outdir, "10_excess_vs_sigma", dpi, 0.10))

    # --- 11 the xc cancellation: NOT a constant of the functional
    fig, ax = figure_line()
    ax.plot(tab.sigma, 100.0 * tab.xc_cancellation, "o-", lw=1.5, ms=4.5,
            color=COLOR["lda"])
    ax.set_xlabel(r"$\sigma_\mathrm{WP}$ (Bohr)")
    ax.set_ylabel("xc cancellation (%)")
    note(ax, r"$1-\dfrac{\mathrm{excess_{lda}}-1}{\mathrm{excess_{hartree}}-1}$",
         "lower right")
    w.append(save(fig, ax, outdir, "11_xc_cancellation_vs_sigma", dpi, 0.12))

    # --- 12 the protocol gate: E_PP * sigma is a CONSTANT in the scaled box
    fig, ax = figure_line()
    ax.plot(gate.sigma, gate.epp0_x_sigma, "o", ms=5, color=COLOR["lda"],
            label="measured")
    pred = 1.0 / math.sqrt(2.0 * math.pi) - SW.XI_MADELUNG / 36.0
    ax.axhline(pred, color="0.4", ls="--", lw=1.1, label="analytic")
    ax.set_xlabel(r"$\sigma_\mathrm{WP}$ (Bohr)")
    ax.set_ylabel(r"$E_{PP}(0)\times\sigma$ (Ha$\cdot$Bohr)")
    ax.set_ylim(pred * 0.985, pred * 1.02)
    legend(ax, loc="upper right")
    note(ax, F_EPP, "lower left")
    w.append(save(fig, ax, outdir, "12_epp_scaling_gate", dpi))

    # --- 13 the trajectories themselves, on the only shared axis (tau)
    fig, ax = figure_line()
    for s in SW.SIGMAS:
        r_ref, r_lda = runs.get((s, "noninteracting")), runs.get((s, "lda"))
        if r_ref is None or r_lda is None:
            continue
        n = min(r_ref.tau.size, r_lda.tau.size)
        ax.plot(r_lda.tau[:n], r_lda.sigma_iso[:n] / r_ref.sigma_iso[:n],
                lw=1.3, label=rf"$\sigma={s:g}$")
    ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
    ax.set_xlabel(r"$\tau=t/\sigma^2$ (dimensionless)")
    ax.set_ylabel(r"$\sigma\,/\,\sigma_\mathrm{non-int}$")
    legend(ax, loc="upper left", ncol=2)
    w.append(save(fig, ax, outdir, "13_lda_excess_trajectories", dpi, 0.22))
    return w


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["sigma4", "sweep", "all"], default="all")
    ap.add_argument("--dpi", type=int, default=style.STYLE_CONFIG["save_dpi"])
    a = ap.parse_args()
    style.apply_theme()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    total = 0
    if a.only in ("sigma4", "all"):
        print("sigma = 4 experiment")
        n = build_sigma4(OUT_ROOT / "sigma4", a.dpi); total += len(n)
        print(f"  wrote {len(n)} figures")
    if a.only in ("sweep", "all"):
        print("sigma sweep")
        n = build_sweep(OUT_ROOT / "sweep", a.dpi); total += len(n)
        print(f"  wrote {len(n)} figures")
    print(f"\n{total} figures under {OUT_ROOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
