#!/usr/bin/env python
"""sigma_WP = 6, v = 3.0 — momentum-transfer heatmap in the (k_z, k_perp) plane.

    Delta P = |psi(k)|^2 at the slab exit  -  |psi(k)|^2 at t = 0

One panel answering "where did the momentum go?": deceleration moves weight
LEFT along k_z, transverse scattering moves it UP the k_perp axis. Red = norm
gained at that momentum, blue = norm lost.

WHY (k_z, k_perp) AND NOT A RADIAL n(|k|). Binning by |k| alone folds the drift
direction into the same coordinate as the transverse spread, so a decelerating
packet and a sideways-heating one both read as "the peak moved left" — different
physics, same picture. This plane separates them. (The run's own
`momentum_distribution.csv` is the radial histogram and is unusable here for
exactly this reason — see docs/reports/report2/drafts/draft1/CLAUDE.md landmine
12c.)

MECHANICS, all delegated to the canonical helpers (landmine 13):
  * `load_complex_vti` — a plain `load_vti` silently returns ONLY the real part,
    which for a drifting packet is a cosine-modulated Gaussian: it looks like a
    plausible orbital and survives visual inspection.
  * `kz_kperp_map` — owns the `ifftshift`-before-`fftn` ordering (inqkit VTIs are
    in PHYSICAL order) and bins k_perp at one transverse grid spacing, the finest
    the data supports.
  * GATE: <k_z>(t=0) must equal k_0 = 3.0. Verified 3.0000000 at build time; the
    script asserts it, because every ordering error shows up here first.

THE JACOBIAN IS IN P, deliberately: P is a probability over (k_z, k_perp), i.e.
the shell sum, so an isotropic-in-plane Gaussian has a RAYLEIGH k_perp marginal
peaking at k_perp = sigma_p, not a Gaussian peaking at 0. That is the honest
"how much norm sits at this transverse momentum" and is what a difference map
must be built from. It is also why the t = 0 map peaks at k_perp ~ 0.27 rather
than 0.

DO NOT read moments off this map — along k_perp every point is assigned its bin
CENTRE and the Rayleigh tail falls steeply across a bin, so sum(k_perp^2 P) is
biased high by a few per cent. The quantitative deltas live in
`make_s6_v3_diagnostics.py`, which takes them from the exact CSV moments.

FRAME CHOICE. Wavefunctions are dumped every 116 steps (t = 4.64 a.u.), so the
nearest frame to the slab exit (t = 13.33) is step 348, t = 13.92 — 0.6 a.u.
past the far face. Stated in the caption rather than interpolated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
REPORT_DIR = REPO / "docs/reports/report2/drafts/draft1/figures/jellium_slab"
REPORT_DPI = 600
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "docs/reports/report2/drafts/draft1/figures"))
import s56_stopping as S                                        # noqa: E402
from _panel import panel_mode, SLOT_IN                          # noqa: E402
from inqview.visualisation import style                         # noqa: E402
from inqview.visualisation.field_io import (                    # noqa: E402
    load_complex_vti, kz_kperp_map)

style.apply_theme()

SIGMA, V = 6.0, 3.0
K0 = 3.0
STEP_EXIT = 348                 # nearest dump to the slab exit (t = 13.92 a.u.)
DT = 0.04
VTI = S.run_dir(SIGMA, V, "wp") / "raw" / "vti" / "wavefunction_wp"


def load_map(step: int):
    f = load_complex_vti(VTI / f"wavefunction_t{step:06d}.vti")
    return kz_kperp_map(f)


def main() -> int:
    kz, kp, P0 = load_map(0)
    kz1, kp1, P1 = load_map(STEP_EXIT)
    assert np.allclose(kz, kz1) and np.allclose(kp, kp1), "grid changed between frames"

    # GATE — every ordering/fftshift error surfaces as a wrong <k_z>(0).
    mz = P0.sum(axis=1)
    kz_mean0 = float((kz * mz).sum() / mz.sum())
    assert abs(kz_mean0 - K0) < 1e-6, f"<k_z>(0) = {kz_mean0}, expected {K0}"
    mz1 = P1.sum(axis=1)
    kz_mean1 = float((kz * mz1).sum() / mz1.sum())

    dP = P1 - P0

    # Crop to where the DIFFERENCE lives, not where P lives: P0's Rayleigh tail
    # runs to k_perp ~ 1.3 carrying no signal, and cropping on it leaves 60 % of
    # the panel empty.
    az, ap = np.abs(dP).sum(axis=1), np.abs(dP).sum(axis=0)
    iz = np.where(az > 1e-3 * az.max())[0]
    ip = np.where(ap > 1e-3 * ap.max())[0]

    # Scale by 1e3 rather than letting matplotlib emit an offset text: the house
    # standard forbids bbox_inches="tight", so a "x10^-2" above the colourbar is
    # CLIPPED by the fixed figure rect (it was, first build).
    if panel_mode():
        w, h = SLOT_IN["half"]
        fig = plt.figure(figsize=(w, h))
        left, bottom, width, height = style._ONE_COL_AXES_RECT
        ax = fig.add_axes((left, bottom, width * 0.86, height))
    else:
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

    outs = ([REPORT_DIR / "slab_panel" / "slab_s6_v3_momentum_map.png"]
            if panel_mode() else
            [HERE / "s6_v3_momentum_map.png",
             REPORT_DIR / "slab_s6_v3_momentum_map.png"])
    for p in outs:
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=REPORT_DPI)
        print(f"wrote {p}")
    plt.close(fig)

    print(f"  frames: t = 0 and t = {STEP_EXIT*DT:.2f} a.u. (step {STEP_EXIT})")
    print(f"  <k_z>: {kz_mean0:.7f} -> {kz_mean1:.7f}  (Delta = {kz_mean1-kz_mean0:+.5f})")
    print(f"  norm gained above k_perp = {kp[ip[0]]:.2f}: "
          f"{dP[:, ip].clip(min=0).sum():.4f}")
    print(f"  max |Delta P| = {vmax:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
