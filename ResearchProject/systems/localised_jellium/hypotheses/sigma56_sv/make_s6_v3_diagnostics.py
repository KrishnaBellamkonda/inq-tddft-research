#!/usr/bin/env python
"""sigma_WP = 6, v = 3.0 — the pairwise interaction ledger and the momentum deltas.

Two figures, both for the ONE twin pair `s6p0_v3p0` / `cl_s6p0_v3p0`:

  s6_v3_interactions.png   classical-vs-WP pairwise electrostatic ledger, built
                           to the same design as the cylindrical-jellium twin
                           figure (docs/reports/report2/.../cylindrical_jellium/
                           make_twin_interactions.py)
  s6_v3_momentum_delta.png the wavepacket's momentum deltas, longitudinal vs
                           transverse

--------------------------------------------------------------------------------
FIGURE 1 — THE LEDGER (design copied from the cylindrical twin figure)
--------------------------------------------------------------------------------
The decomposition (.claude/rules/decomposed-interaction-energies.md) splits the
electrostatic energy into three charge groups P (projectile), S (bath), B
(background):

    E_SS = 1/2 int n_S phi_S     bath-bath
    E_PP = 1/2 int n_P phi_P     projectile SELF-Hartree    <- quantum residual
    E_PS =     int n_S phi_P     projectile-bath            <- does the stopping
    E_SB = -   int n_S phi_+     bath-background
    E_PB = -   int n_P phi_+     projectile-background
    E_BB = 1/2 int n_+ phi_+     background self (constant) -- Delta identically 0

MAIN AXES  Delta E_SS and Delta E_SB, the two large bath terms. They move in near
mirror image: the bath polarises, gaining Hartree energy and losing background
attraction.

INSET  the projectile's own terms, also as Delta from t = 0:
  * E_PS + E_PB summed, because the individual E_PB carries the charged-cell
    G = 0 gauge and is not quotable alone; the sum cancels the monopole parts and
    the Delta cancels the remaining constant, so it is gauge-clean twice over.
    It is the net Coulomb energy the projectile feels from everything external.
  * E_PP kept SEPARATE — it is the projectile's self-Hartree, not an interaction
    with anything external, so folding it into the sum would corrupt the "what
    the projectile feels" reading. The twins differ QUALITATIVELY in it: the
    classical Gaussian is rigid, so its E_PP is constant and its Delta is a flat
    zero, while the packet's falls as it disperses (U ~ 1/a).

Twin parity check: every trace starts at the same value in BOTH halves, because
sigma_pot = sigma_WP/sqrt(2) makes the two projectiles source an identical
potential at injection.

--------------------------------------------------------------------------------
FIGURE 2 — THE MOMENTUM DELTAS
--------------------------------------------------------------------------------
All four curves are momenta in a.u., all as Delta from t = 0, so they share one
axis honestly:

    Delta <p_z>     the drift momentum actually transferred to the medium
    Delta sigma_px  } transverse momentum spread -- the packet being scattered
    Delta sigma_py  }   sideways
    Delta sigma_pz  longitudinal momentum spread

`sigma_pi = sqrt(sigma_pi2)` from wp_momentum_stats. The point the figure makes:
for a packet this wide the transverse channel is nearly inert. A Gaussian starts
at sigma_p = 1/(sqrt2 sigma) = 0.118 a.u. in every direction (verified against
the CSV at t = 0), and the transverse components barely move through the whole
transit while the longitudinal ones carry the interaction.

PLOT RANGE. Both figures stop at the 1 %-norm-loss time (t = 19.6 a.u. here,
`s56_stopping.sigma_r_window`). Past that the CAP is eating the packet and every
moment describes a decaying remnant, not the projectile. The slab transit —
t = (|z0| -+ 12.5)/v = 5.0 to 13.3 a.u. — is shaded so the reader can see that
the whole interaction happens inside the plotted range.

UNITS follow the report-2 convention (docs/reports/report2/drafts/draft1/CLAUDE.md):
atomic units for time, momentum and length; eV for energy.
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
# Every figure in this sweep is mirrored into the report draft (user, 2026-08-05).
# Saved to the report-2 house standard: 600 dpi, no on-canvas title, fixed axes
# rect (bbox_inches=None) -- see docs/reports/report2/drafts/draft1/CLAUDE.md §1.
REPORT_DIR = REPO / "docs/reports/report2/drafts/draft1/figures/jellium_slab"
REPORT_DPI = 600
sys.path.insert(0, str(HERE))
import s56_stopping as S                                        # noqa: E402
from inqview.visualisation import style                         # noqa: E402

style.apply_theme()

HA = 27.211386
SIGMA, V = 6.0, 3.0
LAUNCH_Z, SLAB_HALF = -27.5, 12.5
T_IN = (abs(LAUNCH_Z) - SLAB_HALF) / V           # 5.00 a.u.
T_OUT = (abs(LAUNCH_Z) + SLAB_HALF) / V          # 13.33 a.u.


def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:
    obs = S.run_dir(SIGMA, V, "wp") / "raw" / "observables"
    w = S._concat(obs, "interactions")
    mom = S._concat(obs, "wp_momentum_stats")
    rs = S._concat(obs, "wp_real_space_stats")
    c = S._concat(S.run_dir(SIGMA, V, "classical") / "raw" / "observables",
                  "interactions")
    _m, t_cut, _n, _ok = S.sigma_r_window(rs)
    return w, c, mom, t_cut


def _save(fig, *paths: Path) -> None:
    """Write one render to every path. bbox_inches is NEVER set: the report
    standard fixes the axes rect so on-page font sizes are predictable, and
    `tight` would silently change the figure's physical width."""
    for p in paths:
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=REPORT_DPI)
        print(f"wrote {p}")
    plt.close(fig)


def _d(series: pd.Series) -> np.ndarray:
    """Delta from t = 0, in eV."""
    a = series.to_numpy()
    return (a - a[0]) * HA


def fig_interactions(w: pd.DataFrame, c: pd.DataFrame, t_cut: float,
                     out: Path) -> None:
    w = w[w.time_au <= t_cut]
    c = c[c.time_au <= t_cut]

    fig = plt.figure(figsize=(7.0, 3.0))
    ax = fig.add_axes((0.63 / 7.0, 0.48 / 3.0,
                       1 - 0.63 / 7.0 - 0.1225 / 7.0,
                       1 - 0.48 / 3.0 - 0.105 / 3.0))
    col = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    c_ss, c_pp, c_net, c_sb = col[0], col[1], col[2], col[3]

    for d, ls in ((w, "-"), (c, "--")):
        ax.plot(d.time_au, _d(d.e_ss), ls, lw=1.2, color=c_ss)
        ax.plot(d.time_au, _d(d.e_sb), ls, lw=1.2, color=c_sb)
    ax.axvspan(T_IN, T_OUT, color="0.5", alpha=0.15, lw=0)
    ax.axhline(0.0, lw=0.6, color="0.6")
    ax.set_xlabel(r"$t$ (a.u.)")
    ax.set_ylabel(r"$\Delta E$ (eV)")
    ax.set_xlim(0, t_cut)

    # RESERVE the top half for the inset. Unlike the cylindrical twin (+-15 eV)
    # these bath terms reach +-100 eV and E_SB peaks in the upper RIGHT -- exactly
    # where an inset naturally goes, and the inset's opaque background then hides
    # the peak without any warning. Explicit headroom is the only fix that cannot
    # silently regress when the data changes.
    lo = min(_d(w.e_ss).min(), _d(w.e_sb).min(),
             _d(c.e_ss).min(), _d(c.e_sb).min())
    hi = max(_d(w.e_ss).max(), _d(w.e_sb).max(),
             _d(c.e_ss).max(), _d(c.e_sb).max())
    span = hi - lo
    ax.set_ylim(lo - 0.08 * span, hi + 1.10 * span)

    h = [plt.Line2D([], [], color=c_ss, lw=1.2, label=r"$E_{SS}$"),
         plt.Line2D([], [], color=c_sb, lw=1.2, label=r"$E_{SB}$"),
         plt.Line2D([], [], color="k", lw=1.2, ls="-", label="wavepacket"),
         plt.Line2D([], [], color="k", lw=1.2, ls="--", label="classical")]
    ax.legend(handles=h, fontsize=7, frameon=False, ncol=2, loc="lower left")

    axi = ax.inset_axes((0.30, 0.58, 0.67, 0.39))
    for d, ls in ((w, "-"), (c, "--")):
        axi.plot(d.time_au, _d(d.e_ps + d.e_pb), ls, lw=1.0, color=c_net)
        axi.plot(d.time_au, _d(d.e_pp), ls, lw=1.0, color=c_pp)
    axi.axvspan(T_IN, T_OUT, color="0.5", alpha=0.15, lw=0)
    axi.axhline(0.0, lw=0.5, color="0.6")
    axi.set_xlim(0, t_cut)
    axi.set_xlabel(r"$t$ (a.u.)", fontsize=7, labelpad=1)
    axi.set_ylabel(r"$\Delta E$ (eV)", fontsize=7, labelpad=1)
    axi.tick_params(labelsize=6)
    hi = [plt.Line2D([], [], color=c_net, lw=1.0, label=r"$E_{PS}+E_{PB}$"),
          plt.Line2D([], [], color=c_pp, lw=1.0, label=r"$E_{PP}$")]
    axi.legend(handles=hi, fontsize=6, frameon=False, ncol=2, loc="upper center")

    _save(fig, out, REPORT_DIR / "slab_s6_v3_interactions.png")


def fig_momentum(mom: pd.DataFrame, t_cut: float, out: Path) -> None:
    m = mom[mom.time_au <= t_cut]
    t = m.time_au.to_numpy()

    def dm(a):                     # Delta from t=0, a.u. of momentum
        a = np.asarray(a)
        return a - a[0]

    # CYLINDRICAL decomposition (user, 2026-08-05): the beam axis is z and the
    # slab is uniform in x-y, so x and y are physically equivalent -- the run
    # confirms it, sigma_px and sigma_py agree to every printed digit at every
    # step. Folding them into ONE transverse channel is therefore lossless and
    # removes a redundant curve:
    #     sigma_k_perp = sqrt(sigma_kx^2 + sigma_ky^2)      (2-D transverse)
    #     sigma_kz                                          (longitudinal)
    # <k_perp> itself is identically zero by symmetry, so the transverse channel
    # has no drift term -- only a spread. The longitudinal one has both.
    skperp = np.sqrt(m.sigma_px2.to_numpy() + m.sigma_py2.to_numpy())
    skz = np.sqrt(m.sigma_pz2.to_numpy())
    kz = m.pz_mean.to_numpy()

    fig, ax = style.figure_one_col()
    col = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax.plot(t, dm(kz), lw=1.3, color="k", label=r"$\Delta \langle k_z\rangle$")
    ax.plot(t, dm(skz), lw=1.2, color=col[0], label=r"$\Delta \sigma_{k_z}$")
    ax.plot(t, dm(skperp), lw=1.2, color=col[1],
            label=r"$\Delta \sigma_{k_\perp}$")
    ax.axvspan(T_IN, T_OUT, color="0.5", alpha=0.15, lw=0)
    ax.axhline(0.0, lw=0.6, color="0.6")
    ax.set_xlim(0, t_cut)
    ax.set_xlabel(r"$t$ (a.u.)")
    ax.set_ylabel(r"$\Delta k$ (a.u.)")
    ax.legend(fontsize=7, frameon=False, loc="best")
    _save(fig, out, REPORT_DIR / "slab_s6_v3_momentum_delta.png")

    # provenance the figure cannot carry
    s0 = 1.0 / (np.sqrt(2) * SIGMA)
    print(f"  t=0: sigma_kz = {skz[0]:.4f} (analytic 1/(sqrt2 sigma) = {s0:.4f}), "
          f"sigma_k_perp = {skperp[0]:.4f} (analytic sqrt2 x that = {np.sqrt(2)*s0:.4f})")
    print(f"  transit window t = {T_IN:.2f}-{T_OUT:.2f} a.u.;  plotted to {t_cut:.2f}")
    i = np.abs(t - T_OUT).argmin()
    for nm, a in (("<k_z>", kz), ("sigma_kz", skz), ("sigma_k_perp", skperp)):
        print(f"  {nm:<13s} Delta at slab exit = {a[i]-a[0]:+.5f}   "
              f"Delta at t_cut = {a[-1]-a[0]:+.5f} a.u.   "
              f"({100*(a[-1]-a[0])/a[0]:+.1f}% of its t=0 value)"
              if a[0] else "")


def main() -> int:
    w, c, mom, t_cut = load()
    fig_interactions(w, c, t_cut, HERE / "s6_v3_interactions.png")
    fig_momentum(mom, t_cut, HERE / "s6_v3_momentum_delta.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
