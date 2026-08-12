"""Draft-2 momentum-transfer heatmap: sigma=6, v=3.0, (k_z, k_perp) plane (fig 14).

Runs in standalone (non-panel) mode: style.figure_one_col(with_colorbar=True),
3.5 × 3.0 in, 600 DPI, bbox_inches=None.

Output: slab_s6_v3_momentum_map.png

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_s6_v3_momentum_map.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[8]
HERE = Path(__file__).parent
SRC  = REPO / "ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"

sys.path.insert(0, str(SRC))
sys.path.insert(0, str(REPO / "docs/reports/report2/drafts/draft1/figures"))
sys.path.insert(0, str(REPO / "inq-stack/python"))

# Force non-panel mode before importing the source module
os.environ.pop("PANEL", None)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import s56_stopping as S
from inqview.visualisation import style
from inqview.visualisation.field_io import load_complex_vti, kz_kperp_map

style.apply_theme()

SIGMA, V = 6.0, 3.0
K0 = 3.0
STEP_EXIT = 348
DT = 0.04
VTI = S.run_dir(SIGMA, V, "wp") / "raw" / "vti" / "wavefunction_wp"


def load_map(step: int):
    f = load_complex_vti(VTI / f"wavefunction_t{step:06d}.vti")
    return kz_kperp_map(f)


def main() -> int:
    kz, kp, P0 = load_map(0)
    kz1, kp1, P1 = load_map(STEP_EXIT)
    assert np.allclose(kz, kz1) and np.allclose(kp, kp1), "grid changed between frames"

    mz = P0.sum(axis=1)
    kz_mean0 = float((kz * mz).sum() / mz.sum())
    assert abs(kz_mean0 - K0) < 1e-6, f"<k_z>(0) = {kz_mean0}, expected {K0}"
    mz1 = P1.sum(axis=1)
    kz_mean1 = float((kz * mz1).sum() / mz1.sum())

    dP = P1 - P0

    az, ap = np.abs(dP).sum(axis=1), np.abs(dP).sum(axis=0)
    iz = np.where(az > 1e-3 * az.max())[0]
    ip = np.where(ap > 1e-3 * ap.max())[0]

    fig, ax = style.figure_one_col(with_colorbar=True)
    d = dP * 1e3
    vmax = float(np.abs(d).max())
    m = ax.pcolormesh(kz, kp, d.T, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                      shading="nearest", rasterized=True)
    ax.axvline(K0, lw=0.7, color="0.35", ls="--")
    ax.text(K0, kp[ip[-1]] * 0.96, r"$k_0$", fontsize=7, color="0.35",
            ha="right", va="top")
    ax.set_xlim(kz[iz[0]], kz[iz[-1]])
    ax.set_ylim(0, kp[ip[-1]])
    ax.set_xlabel(r"$k_z$ (a.u.)")
    ax.set_ylabel(r"$k_\perp$ (a.u.)")
    cb = fig.colorbar(m, ax=ax, pad=0.02, fraction=0.06)
    cb.set_label(r"$\Delta |\psi(\mathbf{k})|^2$ ($\times 10^{-3}$)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    out = HERE / "slab_s6_v3_momentum_map.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {out}")
    print(f"  frames: t=0 and t={STEP_EXIT*DT:.2f} a.u. (step {STEP_EXIT})")
    print(f"  <k_z>: {kz_mean0:.7f} -> {kz_mean1:.7f}  (Delta={kz_mean1-kz_mean0:+.5f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
