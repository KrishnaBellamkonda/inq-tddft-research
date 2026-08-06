#!/usr/bin/env python
"""The four sigma_WP = 6 deliverables, as standalone files.

    s6_norm.png           total norm N(t) (and the WP's own norm)
    s6_energy_total.png   E_total(t) - E_GS, raw and corrected
    s6_interactions.png   the pairwise P/S/B interaction energies
    s6_gifs.md            index of the eight total-density GIFs

All eight sigma = 6 production runs (4 WP + 4 classical twins, v = 2.0-3.5).

CONVENTIONS THAT MATTER HERE
* S and the energy traces use E_PS-CORRECTED deposits where an S is quoted --
  see s56_stopping.measure(). The raw E_total curve is still plotted, because the
  gap between raw and corrected IS the point of the classical panel.
* WP norm columns are `norm_wp` / `norm_total`; the classical half writes
  `norm_slab` (the projectile is an external potential, never in the ledger).
  Not interchangeable -- handled explicitly rather than by a "first column
  matching 'norm'" guess, which is how an earlier version picked the wrong one.
* Interaction terms are from raw/observables/interactions.csv
  (.claude/rules/decomposed-interaction-energies.md). E_PP is identically 0 for a
  classical point charge and is the quantum residual for the WP.
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
import s56_stopping as S                                     # noqa: E402

SIGMA = 6.0
VS = (2.0, 2.5, 3.0, 3.5)
HA = 27.211386
COL = {2.0: "#4C72B0", 2.5: "#DD8452", 3.0: "#55A868", 3.5: "#C44E52"}


def obs_of(v: float, half: str) -> Path:
    return S.run_dir(SIGMA, v, half) / "raw" / "observables"


def norm_series(v: float, half: str):
    """(t, total-electron-count, wp-norm-or-None)."""
    ix = S._concat(obs_of(v, half), "interactions")
    t = ix.time_au.to_numpy()
    if half == "wp":
        # norm_total already includes the packet; norm_wp is the packet alone.
        return t, ix["norm_total"].to_numpy(), ix["norm_wp"].to_numpy()
    return t, ix["norm_slab"].to_numpy(), None


def main() -> int:
    # ---- 1. norm ---------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
    for v in VS:
        for half, ls in (("wp", "-"), ("classical", "--")):
            try:
                t, ntot, nwp = norm_series(v, half)
            except Exception as e:
                print(f"  skip norm {half} v={v}: {type(e).__name__}")
                continue
            ax[0].plot(t, ntot, ls=ls, color=COL[v], lw=1.3,
                       label=f"v={v} {'WP' if half=='wp' else 'cl'}")
            if nwp is not None:
                ax[1].semilogy(t, np.maximum(nwp, 1e-12), ls=ls, color=COL[v], lw=1.3,
                               label=f"v={v}")
    ax[0].set_ylabel("total electron count in cell")
    ax[0].set_title("(a) total norm — bath (+ packet)")
    ax[1].set_ylabel(r"$\|\psi_{\rm WP}\|^2$")
    ax[1].set_title("(b) wavepacket norm (WP half only)")
    for a in ax:
        a.set_xlabel("t (a.u.)")
        a.legend(fontsize="x-small", frameon=False, ncol=2)
    plt.tight_layout()
    fig.savefig(HERE / "s6_norm.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"wrote {HERE/'s6_norm.png'}")

    # ---- 2. total energy -------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
    for j, half in enumerate(("wp", "classical")):
        for v in VS:
            try:
                tr = S.energy_trace(SIGMA, v, half)
            except Exception as e:
                print(f"  skip energy {half} v={v}: {type(e).__name__}")
                continue
            ax[j].plot(tr.t, tr.dE_corr, color=COL[v], lw=1.4, label=f"v={v}")
            if half == "wp":
                ax[j].plot(tr.t, tr.dE_raw, color=COL[v], lw=0.8, alpha=0.35, ls=":")
            else:
                # The E_PS-corrected asymptote: what the curve WOULD settle to once
                # the projectile's monopole tail has decayed. Without this the
                # classical panel looks like it plateaus 4x too high.
                eps = S.e_ps_final(SIGMA, v, half)
                ax[j].axhline(tr.dE_corr.iloc[-1] - eps, color=COL[v], lw=0.8, ls=":")
        ti, to = S.transit_window(VS[0])
        ax[j].axvspan(ti, to, color="0.88", zorder=0)
        ax[j].set_xlabel("t (a.u.)")
        ax[j].set_ylabel(r"$E_{\rm total}(t)-E_{\rm GS}$ (eV)")
        ax[j].legend(fontsize="x-small", frameon=False, ncol=2)
    ax[0].set_title("(a) wavepacket — solid corrected, dotted raw")
    ax[1].set_title(r"(b) classical — dotted = $E_{PS}$-corrected asymptote")
    plt.tight_layout()
    fig.savefig(HERE / "s6_energy_total.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"wrote {HERE/'s6_energy_total.png'}")

    # ---- 3. interaction energies ----------------------------------------
    terms = [("e_ps", r"$E_{PS}$ projectile-bath"),
             ("e_pp", r"$E_{PP}$ projectile self"),
             ("e_ss", r"$E_{SS}$ bath-bath"),
             ("e_sb", r"$E_{SB}$ bath-background")]
    fig, axes = plt.subplots(2, 4, figsize=(16, 6.2), sharex=True)
    for row, half in enumerate(("wp", "classical")):
        for k, (c, lab) in enumerate(terms):
            a = axes[row, k]
            for v in VS:
                try:
                    ix = S._concat(obs_of(v, half), "interactions")
                except Exception:
                    continue
                if c not in ix:
                    continue
                a.plot(ix.time_au, (ix[c] - ix[c].iloc[0]) * HA,
                       color=COL[v], lw=1.2, label=f"v={v}")
            ti, to = S.transit_window(VS[0])
            a.axvspan(ti, to, color="0.88", zorder=0)
            a.set_title(f"{'WP' if half=='wp' else 'classical'} — {lab}",
                        fontsize="small")
            if row == 1:
                a.set_xlabel("t (a.u.)")
            if k == 0:
                a.set_ylabel(r"$\Delta E$ from $t=0$ (eV)")
    axes[0, 0].legend(fontsize="x-small", frameon=False, ncol=2)
    plt.tight_layout()
    fig.savefig(HERE / "s6_interactions.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"wrote {HERE/'s6_interactions.png'}")

    # ---- 4. GIF index ----------------------------------------------------
    lines = ["# sigma_WP = 6 — total-density GIFs\n",
             "Mid-y xz slice, physical order (never fftshifted).",
             "Slab faces dashed at z = +/-12.5; CAP inner edges at |z| = 40.\n",
             "| v | half | file |", "|---|---|---|"]
    for v in VS:
        for half in ("wp", "classical"):
            g = sorted((HERE / f"run_{half}_s6_v{v}_figs").glob("*_total_density.gif"))
            lines.append(f"| {v} | {half} | `{g[0].relative_to(HERE)}` |"
                         if g else f"| {v} | {half} | MISSING |")
    (HERE / "s6_gifs.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {HERE/'s6_gifs.md'}")

    # ---- provenance ------------------------------------------------------
    print("\nsigma = 6 summary (E_PS-corrected S):")
    print(f"{'v':>5} {'half':<10} {'S':>8} {'E_PS(tf) eV':>12} {'norm_final':>12}")
    for v in VS:
        for half in ("wp", "classical"):
            try:
                p = S.measure(SIGMA, v, half)
            except Exception:
                continue
            print(f"{v:>5.1f} {half:<10} {p.S_deposit_eV_per_Bohr:>8.3f} "
                  f"{p.e_ps_final_eV:>12.4f} {p.norm_final:>12.4e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
