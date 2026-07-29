#!/usr/bin/env python3
"""build_gs_bore_comparison.py — how much the r_s=6 GS spreads into the bore vs r_s=2.

User question (2026-07-08): plot the difference in NORMALISED ground-state density
between the dense (r_s=2) and dilute (r_s=6) walls, to see how far the r_s=6 GS spills
into the hollow bore (where the on-axis projectile flies) and might contribute to
stopping.

Method.
- GS electron density = the **t=0 frame** of each classical run's `density_system`
  (the projectile is a pseudo-ion — it perturbs the potential, NOT the electron
  density — so frame-0 is the clean wall GS; verified frame-0 == the dedicated GS VTI
  for r_s=2). Loaded via canonical `inqview.load_vti` (PHYSICAL order, never fftshift).
- Cylindrically + axially averaged radial profile n(d), d=√(x²+y²), on a shared radial
  grid (both runs share the transverse 80×80 grid, so the difference is well-defined).
- NORMALISED by each wall's plateau density n_wall = ⟨n⟩ over 7<d<11 Bohr (≈ n₀), so
  the wall sits at 1.0 for both and the bore fraction is directly comparable.

Panels: (A) normalised profiles overlaid (semilogy) with the bore d<R_in shaded;
(B) the requested normalised DIFFERENCE Δ = (n/n_wall)_rs6 − (n/n_wall)_rs2;
(C) ABSOLUTE profiles (honesty: the projectile feels absolute density) — r_s=2's bore
density is still higher in absolute terms even though r_s=6 spreads more relative to
its bulk. Interpretation is left to the reader (provisional).

Writes gs_validation/gs_bore_spread_rs2_vs_rs6.png. Run:
    PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_gs_bore_comparison.py
"""
from __future__ import annotations
import glob
import math
import sys
from pathlib import Path

import numpy as np

STACK = "/local/data/public/skcb2/tddft/inq-stack/python"
if STACK not in sys.path:
    sys.path.insert(0, STACK)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from inqview import load_vti  # noqa: E402
from inqview.visualisation import style  # noqa: E402
try:
    style.apply()
except Exception:
    pass

SYS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium")
HERE = SYS / "hypotheses" / "annular_sv"
OUT = HERE / "gs_validation" / "gs_bore_spread_rs2_vs_rs6.png"
R_IN, R_OUT = 5.0, 13.0
RUNS = {6: "rs6_v0p15", 2: "rs2_v0p15"}     # frame-0 = clean wall GS
NBINS = 80


def gs_frame0(run):
    f = sorted(glob.glob(str(SYS / "annular_sv" / run / "results" / "**" /
                             "density_system" / "*.vti"), recursive=True))[0]
    vf = load_vti(f)
    return np.asarray(vf.data), np.asarray(vf.x), np.asarray(vf.y)


def radial_profile(n, x, y, nbins, dmax):
    X, Y = np.meshgrid(x, y, indexing="ij")
    D = np.sqrt(X ** 2 + Y ** 2)
    n_xy = n.mean(axis=2)                     # z-average (uniform tube)
    edges = np.linspace(0, dmax, nbins + 1)
    c = 0.5 * (edges[:-1] + edges[1:])
    prof = np.full(nbins, np.nan)
    for i in range(nbins):
        m = (D >= edges[i]) & (D < edges[i + 1])
        if m.any():
            prof[i] = n_xy[m].mean()
    return c, prof


def main():
    prof = {}
    for rs, run in RUNS.items():
        n, x, y = gs_frame0(run)
        dmax = float(min(x.max(), y.max()))
        c, p = radial_profile(n, x, y, NBINS, dmax)
        wall = np.nanmean(p[(c > 7) & (c < 11)])
        n0 = 3 / (4 * math.pi * rs ** 3)
        prof[rs] = dict(d=c, n=p, wall=wall, n0=n0, norm=p / wall)
        onaxis = np.nanmean(p[c < 1.5]) / wall
        print(f"rs={rs}: n0={n0:.5f} wall={wall:.5f} on-axis(d<1.5)/wall={onaxis:.3f}")

    d = prof[6]["d"]                          # shared grid
    diff = prof[6]["norm"] - prof[2]["norm"]  # (n/n_wall)_rs6 − (n/n_wall)_rs2
    bore = d < R_IN

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # (A) normalised profiles
    for rs, col in [(2, "C3"), (6, "C0")]:
        a.semilogy(prof[rs]["d"], prof[rs]["norm"], col, lw=1.8,
                   label=f"r_s={rs}  (wall n₀≈{prof[rs]['n0']:.4f})")
    for xx in (R_IN, R_OUT):
        a.axvline(xx, ls="--", lw=0.7, color="0.4")
    a.axvspan(0, R_IN, color="0.85", alpha=0.5)
    a.text(R_IN / 2, a.get_ylim()[1] * 0.4, "projectile", ha="center",
           fontsize=8, color="0.3")
    a.set_xlabel("radial distance d = √(x²+y²) (Bohr)")
    a.set_ylabel("normalised GS density  n(d)/n_wall")
    a.set_title("(A) Normalised radial GS density", fontsize=10)
    a.legend(fontsize=8); a.grid(alpha=.25, which="both")

    # (B) normalised difference — the requested plot
    b.plot(d, diff, "C4", lw=1.9)
    b.axhline(0, ls=":", color="0.5", lw=0.9)
    for xx in (R_IN, R_OUT):
        b.axvline(xx, ls="--", lw=0.7, color="0.4")
    b.axvspan(0, R_IN, color="0.85", alpha=0.5)
    b.fill_between(d, 0, diff, where=(diff > 0), color="C0", alpha=0.18)
    b.set_xlabel("radial distance d (Bohr)")
    b.set_ylabel(r"$\Delta$ = (n/n_wall)$_{r_s=6}$ − (n/n_wall)$_{r_s=2}$")
    b.set_title("(B) Normalised difference (rs6 − rs2)\n>0 ⇒ rs6 spreads more (rel. to bulk)",
                fontsize=10)
    b.grid(alpha=.25)

    # (C) absolute profiles (honesty)
    for rs, col in [(2, "C3"), (6, "C0")]:
        c.semilogy(prof[rs]["d"], prof[rs]["n"], col, lw=1.8, label=f"r_s={rs}")
    for xx in (R_IN, R_OUT):
        c.axvline(xx, ls="--", lw=0.7, color="0.4")
    c.axvspan(0, R_IN, color="0.85", alpha=0.5)
    c.set_xlabel("radial distance d (Bohr)")
    c.set_ylabel("absolute GS density  n(d) (a₀⁻³)")
    c.set_title("(C) Absolute density (what the projectile feels)", fontsize=10)
    c.legend(fontsize=8); c.grid(alpha=.25, which="both")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT, dpi=160)
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
