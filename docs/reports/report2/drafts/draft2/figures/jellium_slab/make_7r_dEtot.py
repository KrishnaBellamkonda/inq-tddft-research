"""Draft-2 case study panel: ΔE_total vs time, WP + classical (fig 13).

sigma_WP=2, v=2.0, localised jellium r_s=4.18, WITH CAP.
Sources:
  WP        : hypotheses/wp_highdensity_sv/sweep_data/s2p0_v2p0
  classical : scripts/classical_highdensity_sv/dyn_direct_cap/results/
              s2p0_v2p0_cap/raw/observables

Projectile-background correction (classical only)
-------------------------------------------------
INQ's ``energy_total`` treats the classical projectile as an EXTERNAL potential
acting on the electrons, so its ledger carries the projectile-bath term
(``e_ps``, via ``energy_external``) but NOT the projectile-background term
``e_pb`` — the projectile is never a charge that INQ couples to the neutralising
background. Verified from the closure relations at t=0:
  classical : energy_external == e_sb + e_ps           (e_pb ABSENT)
  wavepacket: energy_external == e_sb + e_pb            (e_pb PRESENT)
In the WP run the projectile is a real KS orbital, so ``e_pb`` is already inside
``energy_total``. To compare like-for-like we therefore ADD ``e_pb`` back to the
classical total (from interactions.csv). Without it the classical ΔE_total
plateaus at an unphysical -78 eV (energy apparently leaving the system); with it
the plateau is a physical +12 eV of deposited energy. Only ``e_pb`` is time-
varying among the omitted terms (e_pp, e_bb are constant), so only it shifts ΔE.

Output: 7r_dEtot_vs_time_both_preview.png (600 DPI, bbox_inches=None)

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_7r_dEtot.py
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

REPO = Path(__file__).resolve().parents[7]
sys.path.insert(0, str(REPO / "inq-stack/python"))

from inqview.visualisation import style

HERE = Path(__file__).parent
LJ   = REPO / "ResearchProject/systems/localised_jellium"
HA   = 27.211386

WP_DIR = str(LJ / "hypotheses/wp_highdensity_sv/sweep_data/s2p0_v2p0")
CL_DIR = str(LJ / "scripts/classical_highdensity_sv/dyn_direct_cap/results"
                  "/s2p0_v2p0_cap/raw/observables")

WPC = style.WP         # tab:red  — wavepacket
CLC = style.CLASSICAL  # tab:blue — classical projectile


def cat(d: str, stem: str) -> pd.DataFrame | None:
    fs = sorted(glob.glob(f"{d}/{stem}.csv") + glob.glob(f"{d}/{stem}.from*.csv"))
    if not fs:
        return None
    return (pd.concat([pd.read_csv(f, comment="#") for f in fs])
            .drop_duplicates("step").sort_values("step").reset_index(drop=True))


def main() -> None:
    _wp_ob_corr = cat(WP_DIR, "observables_corrected")
    wp_ob = _wp_ob_corr if _wp_ob_corr is not None else cat(WP_DIR, "observables")
    if wp_ob is None:
        raise FileNotFoundError("Could not load WP observables.")

    wp_E = (wp_ob["energy_total_corrected"]
            if "energy_total_corrected" in wp_ob
            else wp_ob["energy_total"]).to_numpy()
    wp_t = wp_ob["time_au"].to_numpy()
    wp_dE = (wp_E - wp_E[0]) * HA

    n = len(wp_dE)
    wp_plat = float(np.mean(wp_dE[int(0.8 * n):]))

    # --- classical projectile: energy_total + e_pb (projectile-background) -----
    cl_ob = cat(CL_DIR, "observables")
    cl_it = cat(CL_DIR, "interactions")
    if cl_ob is None or cl_it is None:
        raise FileNotFoundError(f"Could not load classical observables/interactions from {CL_DIR}")
    cl = pd.merge(cl_ob[["step", "time_au", "energy_total"]],
                  cl_it[["step", "e_pb"]], on="step")
    cl_t = cl["time_au"].to_numpy()
    # INQ omits e_pb for a classical (external-potential) projectile; add it back
    # so the total matches the WP ledger (which already contains e_pb).
    cl_E = (cl["energy_total"] + cl["e_pb"]).to_numpy()
    cl_dE = (cl_E - cl_E[0]) * HA
    nc = len(cl_dE)
    cl_plat = float(np.mean(cl_dE[int(0.8 * nc):]))

    style.apply_theme()
    fig, ax = style.figure_one_col()

    ax.plot(wp_t, wp_dE, "-", color=WPC, lw=1.3, label="wavepacket")
    ax.axhline(wp_plat, color=WPC, lw=0.9, ls="--", zorder=1)

    ax.plot(cl_t, cl_dE, "-", color=CLC, lw=1.3, label="classical")
    ax.axhline(cl_plat, color=CLC, lw=0.9, ls="--", zorder=1)

    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\Delta E_\mathrm{total}$ (eV)")
    ax.legend(fontsize=7, frameon=False, loc="upper right")

    _tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.97, wp_plat + 2, f"${round(wp_plat):d}$ eV", transform=_tr,
            color=WPC, fontsize=7, va="bottom", ha="right")
    ax.text(0.97, cl_plat + 2, f"${round(cl_plat):d}$ eV", transform=_tr,
            color=CLC, fontsize=7, va="bottom", ha="right")

    out = HERE / "7r_dEtot_vs_time_both_preview.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {out}")
    print(f"  WP plateau:         {wp_plat:.1f} eV")
    print(f"  classical plateau:  {cl_plat:.1f} eV  (incl. e_pb; "
          f"raw energy_total would be {cl_plat - (cl['e_pb'].iloc[int(0.8*nc):].mean() - cl['e_pb'].iloc[0])*HA:.1f} eV)")


if __name__ == "__main__":
    main()
