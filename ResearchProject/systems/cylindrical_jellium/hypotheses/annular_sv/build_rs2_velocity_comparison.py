#!/usr/bin/env python3
"""build_rs2_velocity_comparison.py — r_s=2 velocity-comparison overlay plots.

User (2026-07-08): for the r_s=2 tube, overlay the different-velocity runs as traces
in ONE plot, for (1) the total energy and (2) the total current density — to compare
across velocity and see whether the induced current correlates with projectile speed.

Velocities overlaid: v = {0.15, 0.30, 0.45} a.u. (0.30 added as the midpoint to sharpen
the current–velocity trend; the user named 0.15 and 0.45).

Method / what is plotted.
- Data: each run's `observables.csv` (time_au, energy_total, current_{x,y,z}).
- **Energy plot**: ΔE_total(t) = energy_total(t) − energy_total(0), one trace per
  velocity. Referenced to t=0 because all three runs share the IDENTICAL r_s=2 GS
  baseline E(0) = −213.25 Ha, and the absolute value (~−213) would hide the ~0.04–0.14
  Ha deposit on a raw axis. So ΔE = the energy the projectile deposits into the gas.
- **Current plot**: current_z(t), one trace per velocity. current_x, current_y are
  ~1000× smaller (on-axis symmetry), so the axial current_z IS effectively the total
  induced current density; peak |current_z| ∝ v is the flow→current correlation.

Writes two PNGs into rs2_velocity_figs/. Run:
    PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_rs2_velocity_comparison.py
"""
from __future__ import annotations
import glob
import sys
from pathlib import Path

STACK = "/local/data/public/skcb2/tddft/inq-stack/python"
if STACK not in sys.path:
    sys.path.insert(0, STACK)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import numpy as np  # noqa: E402
from inqview.visualisation import style  # noqa: E402
try:
    style.apply()
except Exception:
    pass

SYS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium")
OUTDIR = SYS / "hypotheses" / "annular_sv" / "rs2_velocity_figs"
VELS = [(0.15, "C0"), (0.30, "C1"), (0.45, "C3")]


def load(v):
    run = f"rs2_v{v:.2f}".replace(".", "p")
    o = glob.glob(str(SYS / "annular_sv" / run / "results" / "**" / "observables.csv"),
                  recursive=True)[0]
    return pd.read_csv(o).drop_duplicates("step").sort_values("time_au")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    data = {v: load(v) for v, _ in VELS}

    # (1) total energy: ΔE(t) per velocity
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for v, c in VELS:
        O = data[v]
        dE = O["energy_total"].to_numpy() - O["energy_total"].iloc[0]
        ax.plot(O["time_au"], dE, c, lw=1.6,
                label=f"v = {v:.2f}  (peak ΔE = {dE.max():.3f} Ha)")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\Delta E_\mathrm{total} = E(t) - E(0)$  (Ha)")
    ax.set_title("r_s = 2 tube: energy deposited vs time, by projectile velocity\n"
                 "(shared GS baseline E(0) = −213.25 Ha)", fontsize=10)
    ax.legend(title="projectile velocity"); ax.grid(alpha=.25)
    f1 = OUTDIR / "compare_energy_vs_velocity.png"
    fig.tight_layout(); fig.savefig(f1, dpi=160); plt.close(fig)

    # (2) total (axial) current density: current_z(t) per velocity
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    peaks = []
    for v, c in VELS:
        O = data[v]
        jz = O["current_z"].to_numpy()
        pk = np.abs(jz).max(); peaks.append((v, pk))
        ax.plot(O["time_au"], jz, c, lw=1.6,
                label=f"v = {v:.2f}  (peak |J_z| = {pk:.3f})")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"axial current density $J_z$ (a.u.)")
    ax.set_title("r_s = 2 tube: induced total (axial) current vs time, by velocity\n"
                 "(J_x, J_y ~1000× smaller — on-axis symmetry)", fontsize=10)
    ax.legend(title="projectile velocity"); ax.grid(alpha=.25)
    f2 = OUTDIR / "compare_current_vs_velocity.png"
    fig.tight_layout(); fig.savefig(f2, dpi=160); plt.close(fig)

    print("peak |J_z| vs v:", [(v, round(p, 4)) for v, p in peaks])
    print(f"wrote {f1}\nwrote {f2}")


if __name__ == "__main__":
    main()
