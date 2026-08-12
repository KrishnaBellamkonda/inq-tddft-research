"""Draft-2 case study panel: ΔE_total vs time for WP + classical (fig 13).

Extracts the "7r" report panel from the make_case_study.py script (preview mode):
  sigma_WP=2, v=2.0, localised jellium r_s=4.18.

Canonical colour convention:
  WP        → tab:red   (WPC)
  Classical → tab:blue  (CLC)

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

# Canonical WP/classical colours (draft-2 standard)
WPC = "tab:red"    # wavepacket
CLC = "tab:blue"   # classical


def cat(d: str, stem: str) -> pd.DataFrame | None:
    fs = sorted(glob.glob(f"{d}/{stem}.csv") + glob.glob(f"{d}/{stem}.from*.csv"))
    if not fs:
        return None
    return (pd.concat([pd.read_csv(f, comment="#") for f in fs])
            .drop_duplicates("step").sort_values("step").reset_index(drop=True))


def main() -> None:
    _wp_ob_corr = cat(WP_DIR, "observables_corrected")
    wp_ob = _wp_ob_corr if _wp_ob_corr is not None else cat(WP_DIR, "observables")
    cl_ob = cat(CL_DIR, "observables")

    if wp_ob is None or cl_ob is None:
        raise FileNotFoundError(
            "Could not load observables. Ensure the sweep_data and classical run "
            "directories are accessible.")

    wp_E = (wp_ob["energy_total_corrected"]
            if "energy_total_corrected" in wp_ob
            else wp_ob["energy_total"]).to_numpy()
    wp_t = wp_ob["time_au"].to_numpy()
    wp_dE = (wp_E - wp_E[0]) * HA

    cl_t = cl_ob["time_au"].to_numpy()
    cl_dE = (cl_ob["energy_total"].to_numpy() - cl_ob["energy_total"].to_numpy()[0]) * HA

    n = len(wp_dE)
    wp_plat = float(np.mean(wp_dE[int(0.8 * n):]))
    cl_plat = float(np.mean(cl_dE[int(0.8 * len(cl_dE)):]))

    style.apply_theme()
    fig, ax = style.figure_one_col()

    ax.plot(wp_t, wp_dE, "-", color=WPC, lw=1.3, label="wavepacket")
    ax.plot(cl_t, cl_dE, "-", color=CLC, lw=1.3, label="classical")
    ax.axhline(wp_plat, color=WPC, lw=0.9, ls="--", zorder=1)
    ax.axhline(cl_plat, color=CLC, lw=0.9, ls="--", zorder=1)
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\Delta E_\mathrm{total}$ (eV)")
    ax.legend(fontsize=7, frameon=False, loc="upper right")

    _tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.97, wp_plat, f"${round(wp_plat):d}$ eV", transform=_tr,
            color=WPC, fontsize=7, va="bottom", ha="right")
    ax.text(0.97, cl_plat, f"${round(cl_plat):d}$ eV", transform=_tr,
            color=CLC, fontsize=7, va="top", ha="right")

    out = HERE / "7r_dEtot_vs_time_both_preview.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {out}")
    print(f"  WP plateau:  {wp_plat:.1f} eV")
    print(f"  Classical plateau: {cl_plat:.1f} eV")


if __name__ == "__main__":
    main()
