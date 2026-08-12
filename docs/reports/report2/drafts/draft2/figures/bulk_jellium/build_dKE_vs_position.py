"""Draft-2 ΔKE vs position figure (fig 6).

Three traces: classical + WP T1/T2 kinetic decomposition.
Canonical colour convention applied:
  Classical Δ(½mv²)   → tab:blue
  WP T1  Δ⟨p²⟩/2m    → tab:red   (total WP kinetic)
  WP T2  Δ⟨p⟩²/2m    → #762a83   (mean-momentum term, matches case_study T2)

Output: dKE_vs_position_rs5p7.png  (600 DPI, bbox_inches=None)

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 build_dKE_vs_position.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[7]   # tddft root
SYSTEM = REPO / "ResearchProject/systems/jellium"
KS_HYP = SYSTEM / "hypotheses/bulk_ks_stopping"

sys.path.insert(0, str(KS_HYP))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import ks_stopping as K
from inqview.visualisation import style

HA_TO_EV = 27.211386

BASE = SYSTEM / "scripts/bulk_ks_stopping"
LZ   = 80.0
Z0   = -32.0

T0_CL, T1_CL = 4.0, 18.9719
T0_T2, T1_T2 = 4.0,  9.37
T0_T1, T1_T1 = 9.37, 18.9719

# Canonical colours (draft-2 standard)
COLOR_CL = "tab:blue"   # classical
COLOR_T1 = "tab:red"    # WP total kinetic (T1 = ⟨p²⟩/2m)
COLOR_T2 = "#762a83"    # WP mean-momentum term (T2 = ⟨p⟩²/2m)

HERE = Path(__file__).parent


def main() -> None:
    cl = K.load_classical_run(BASE / "classical", box_length_z=LZ)
    wp = K.load_wp_run(BASE / "wp", box_length_z=LZ, z0=Z0)

    dT_cl = (cl.T  - cl.T[0])  * HA_TO_EV
    dT1   = (wp.T1 - wp.T1[0]) * HA_TO_EV
    dT2   = (wp.T2 - wp.T2[0]) * HA_TO_EV

    fit_cl = K.fit_classical(cl, T0_CL, T1_CL)
    fit_T1 = K.fit_stopping(wp.s3, wp.T1, wp.t, T0_T1, T1_T1,
                             label="S_T1 (<p^2>/2m)")
    fit_T2 = K.fit_stopping(wp.s3, wp.T2, wp.t, T0_T2, T1_T2,
                             label="S_T2 (<p>^2/2m)")

    z_cl0  = float(np.interp(T0_CL, cl.t, cl.z))
    z_cl1  = float(np.interp(T1_CL, cl.t, cl.z))
    z_t2_0 = float(np.interp(T0_T2, wp.t, wp.s3))
    z_t2_1 = float(np.interp(T1_T2, wp.t, wp.s3))
    z_t1_1 = float(np.interp(T1_T1, cl.t, cl.z))

    style.apply_theme()
    fig, ax = style.figure_one_col()

    ax.axvspan(z_t2_0, z_t2_1, color="#B8D4E8", alpha=0.55, zorder=0)
    ax.axvspan(z_t2_1, z_t1_1, color="#F7C5A8", alpha=0.45, zorder=0)

    ax.plot(cl.z,  dT_cl, lw=1.4, color=COLOR_CL, alpha=0.85, zorder=2,
            label=r"$\Delta(\frac{1}{2}mv^2)$ (classical)")
    ax.plot(wp.s3, dT1,   lw=1.4, color=COLOR_T1, alpha=0.85, zorder=2, ls="-",
            label=r"$\Delta\langle p^2\rangle/2m$")
    ax.plot(wp.s3, dT2,   lw=1.4, color=COLOR_T2, alpha=0.85, zorder=2, ls="--",
            label=r"$\Delta\langle p\rangle^2/2m$")

    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$z\ (\mathrm{Bohr})$")
    ax.set_ylabel(r"$\Delta KE\ (\mathrm{eV})$")
    ax.legend(fontsize=7, frameon=False, loc="lower left")

    out = HERE / "dKE_vs_position_rs5p7.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {out}")

    def _fmt_s(fit):
        s, u = fit.S_ev_per_bohr, fit.uncertainty
        if u > 0:
            import math
            mag = 10 ** math.floor(math.log10(abs(u)))
            u2 = round(u / mag) * mag
            s2 = round(s / mag) * mag
            dec = max(0, -int(math.floor(math.log10(mag))))
            return f"{s2:.{dec}f} ± {u2:.{dec}f}"
        return f"{s:.3f} ± ?"

    print(f"  Classical:  S = {_fmt_s(fit_cl)} eV/Bohr")
    print(f"  WP T1:      S = {_fmt_s(fit_T1)} eV/Bohr")
    print(f"  WP T2:      S = {_fmt_s(fit_T2)} eV/Bohr")


if __name__ == "__main__":
    main()
