#!/usr/bin/env python
"""sigma_WP = 6, v = 3.0 (the case-study pair) — Delta E_total, classical vs WP.

    Delta E_total(t) = E_total(t) - E_GS        [eV]

The figure exists to show the PLATEAU: the deposit rises while the projectile is
inside the slab and then stops changing, which is what licenses quoting a single
number S = Delta E_total(t_final)/L_slab. A curve still drifting at t_final would
mean the estimator is being read off a transient.

THE NORM CORRECTION, AND WHY ONLY THE WP HALF GETS IT. INQ reports an orbital's
kinetic term as occ*<psi|T|psi>/<psi|psi> (inq/src/hamiltonian/energy.hpp:50-55).
As the CAP eats the packet its norm -> 0 while that ratio stays finite, so the
RAW E_total keeps a full particle's worth of kinetic energy for a packet that is
no longer there. Corrected:

    E_corr = E_raw - T_orb * (1 - norm)

A classical run has no WP orbital, so raw == corrected there and the classical
curve is untouched. Plotting the raw WP curve instead would show a large late
rise that is pure bookkeeping -- at v = 3.0 the raw deposit reaches 15 eV against
a corrected 9.5 eV. Both are drawn (raw as a faint dotted line) precisely so the
correction is visible rather than asserted.

CASE STUDY. sigma_WP = 6, v = 3.0 -- the pair already used for the interaction
ledger and the momentum map, so the three panels describe one physical event.
The shaded band is the slab transit, t = (|z0| -+ 12.5)/v = 5.0 to 13.3 a.u.

Outputs (see .claude memory `feedback-slab-figures-to-report2`):
    s6_v3_energy_delta.png                       standalone, beside this script
    <report>/jellium_slab/slab_s6_v3_energy_delta.png
    PANEL=1 -> <report>/jellium_slab/slab_panel/slab_s6_v3_energy_delta.png
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
REPORT_FIGS = REPO / "docs/reports/report2/drafts/draft1/figures"
REPORT_DIR = REPORT_FIGS / "jellium_slab"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPORT_FIGS))            # the width contract lives here
import s56_stopping as S                                        # noqa: E402
from _panel import panel_mode, slot_figure                      # noqa: E402
from inqview.visualisation import style                         # noqa: E402

style.apply_theme()

HA = 27.211386
SIGMA, V = 6.0, 3.0
LAUNCH_Z, SLAB_HALF = -27.5, 12.5
T_IN, T_OUT = (abs(LAUNCH_Z) - SLAB_HALF) / V, (abs(LAUNCH_Z) + SLAB_HALF) / V


def wp_curve(e_gs_ev: float):
    obs = S.run_dir(SIGMA, V, "wp") / "raw" / "observables"
    d = S._concat(obs, "observables")
    mom = S._concat(obs, "wp_momentum_stats")
    pos = S._concat(obs, "wp_real_space_stats")
    m = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))
    m = m[m.step.isin(d.step)]
    d = d[d.step.isin(m.step)]
    norm = (m["norm_check_r"] if "norm_check_r" in m else m["norm_check"]).to_numpy()
    t_orb = m["e_kin_ha"].to_numpy() * HA
    raw = d["energy_total"].to_numpy() * HA
    n = min(len(raw), len(norm))
    corr = raw[:n] - t_orb[:n] * (1.0 - norm[:n])
    return d["time_au"].to_numpy()[:n], raw[:n] - e_gs_ev, corr - e_gs_ev


def cl_curve(e_gs_ev: float):
    """(t, raw, corrected) for the classical half.

    THE CLASSICAL HALF NEEDS ITS OWN CORRECTION, and it is NOT the norm one.
    A classical projectile is an external potential that keeps flying (z reaches
    321 Bohr here), so E_PS -- its Coulomb interaction with the bath -- stays in
    `energy_total` for the whole run, decaying only as the 1/z monopole tail.
    The raw curve therefore never plateaus: measured drift over the last 10 % of
    this run is -1.04 eV, still falling at t_final. Subtracting the instantaneous
    E_PS leaves the energy actually deposited in the bath:

        Delta E_dep(t) = E_total(t) - E_GS - E_PS(t)

    which is the `S_deposit_eV_per_Bohr` column of s56_S_summary.csv and the
    estimator the report mandates (draft1/CLAUDE.md landmine 1). The WP half
    needs no such term: the CAP annihilates the packet, so E_PS(t_f) ~ 1e-5 eV.
    """
    obs = S.run_dir(SIGMA, V, "classical") / "raw" / "observables"
    d = S._concat(obs, "observables")
    ix = S._concat(obs, "interactions")
    m = d.merge(ix[["step", "e_ps"]], on="step")
    raw = m["energy_total"].to_numpy() * HA - e_gs_ev
    corr = raw - m["e_ps"].to_numpy() * HA
    return m["time_au"].to_numpy(), raw, corr


def main() -> int:
    e_gs_ev = S.e_gs_ha() * HA
    t_w, raw_w, corr_w = wp_curve(e_gs_ev)
    t_c, raw_c, corr_c = cl_curve(e_gs_ev)

    fig, ax = (slot_figure("half") if panel_mode() else style.figure_one_col())
    col = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    ax.plot(t_c, corr_c, lw=1.3, color=col[3], label="classical")
    ax.plot(t_w, corr_w, lw=1.3, color=col[0], label="wavepacket")
    ax.plot(t_c, raw_c, lw=0.8, ls=":", color=col[3], alpha=0.75,
            label="uncorrected")
    ax.plot(t_w, raw_w, lw=0.8, ls=":", color=col[0], alpha=0.75)
    ax.axvspan(T_IN, T_OUT, color="0.5", alpha=0.15, lw=0)
    ax.axhline(0.0, lw=0.6, color="0.6")
    ax.set_xlim(0, min(t_c.max(), t_w.max()))
    ax.set_xlabel(r"$t$ (a.u.)")
    ax.set_ylabel(r"$\Delta E_\mathrm{dep}$ (eV)")
    # The plateau is the POINT of this figure and it sits near 3-10 eV, while the
    # uncorrected transients peak at ~200 eV. Autoscaling on all four curves
    # squashes the plateau into the axis line, so the limit is set from the
    # CORRECTED pair and the dotted raw curves are allowed to run off the top.
    lo = min(corr_c.min(), corr_w.min())
    hi = max(corr_c.max(), corr_w.max())
    ax.set_ylim(lo - 0.10 * (hi - lo), hi + 0.40 * (hi - lo))
    ax.legend(fontsize=7, frameon=False, loc="upper right")

    out = [HERE / "s6_v3_energy_delta.png",
           REPORT_DIR / "slab_s6_v3_energy_delta.png"]
    if panel_mode():
        out = [REPORT_DIR / "slab_panel" / "slab_s6_v3_energy_delta.png"]
    for p in out:
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=600)
        print(f"wrote {p}")
    plt.close(fig)

    # plateau evidence, which the figure shows but does not quantify
    for nm, a in (("classical corr", corr_c), ("classical RAW ", raw_c),
                  ("wavepacket corr", corr_w), ("wavepacket RAW ", raw_w)):
        i90 = int(0.9 * (len(a) - 1))
        print(f"  {nm:<16s} final = {a[-1]:8.3f} eV   "
              f"drift over last 10% = {a[-1]-a[i90]:+8.4f} eV   "
              f"S = {a[-1]/25:6.3f} eV/Bohr")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
