"""ΔKE vs projectile position — figures-2 version (2.948 × 3.00 in).

Slot: 0.48\\textwidth = 2.948 in at 6.142 in textwidth.
Labels already in math mode; no title-case change needed.

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 build_dKE_vs_position.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
SYSTEM = REPO / "ResearchProject/systems/jellium"
KS_HYP = SYSTEM / "hypotheses/bulk_ks_stopping"

sys.path.insert(0, str(KS_HYP))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ks_stopping as K
from inqview.visualisation import style

style.apply_theme()

HERE = Path(__file__).parent
HA_TO_EV = 27.211386
BASE = SYSTEM / "scripts/bulk_ks_stopping"
LZ   = 80.0
Z0   = -32.0

T0_CL, T1_CL = 4.0, 18.9719
T0_T2, T1_T2 = 4.0,  9.37
T0_T1, T1_T1 = 9.37, 18.9719


def main() -> None:
    cl = K.load_classical_run(BASE / "classical", box_length_z=LZ)
    wp = K.load_wp_run(BASE / "wp", box_length_z=LZ, z0=Z0)

    dT_cl = (cl.T  - cl.T[0])  * HA_TO_EV
    dT1   = (wp.T1 - wp.T1[0]) * HA_TO_EV
    dT2   = (wp.T2 - wp.T2[0]) * HA_TO_EV

    z_t2_0 = float(np.interp(T0_T2, wp.t, wp.s3))
    z_t2_1 = float(np.interp(T1_T2, wp.t, wp.s3))
    z_t1_1 = float(np.interp(T1_T1, cl.t, cl.z))

    fig = plt.figure(figsize=(2.948, 3.0))
    ax  = fig.add_axes((0.20, 0.16, 0.77, 0.80))

    ax.axvspan(z_t2_0, z_t2_1, color="#B8D4E8", alpha=0.55, zorder=0)
    ax.axvspan(z_t2_1, z_t1_1, color="#F7C5A8", alpha=0.45, zorder=0)

    ax.plot(cl.z,  dT_cl, lw=1.4, color="C3", alpha=0.85, zorder=2,
            label=r"$\Delta(\frac{1}{2}mv^2)$")
    ax.plot(wp.s3, dT1,   lw=1.4, color="C0", alpha=0.85, zorder=2, ls="-",
            label=r"$\Delta\langle p^2\rangle/2m$")
    ax.plot(wp.s3, dT2,   lw=1.4, color="C0", alpha=0.85, zorder=2, ls="--",
            label=r"$\Delta\langle p\rangle^2/2m$")

    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$z\ (\mathrm{Bohr})$")
    ax.set_ylabel(r"$\Delta KE\ (\mathrm{eV})$")
    ax.legend(fontsize=7, frameon=False, loc="lower left")

    out = HERE / "dKE_vs_position_rs5p7.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
