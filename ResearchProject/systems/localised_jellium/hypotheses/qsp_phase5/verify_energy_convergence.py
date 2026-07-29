#!/usr/bin/env python3
"""Verify that E_total(t) converges to a minimum for the valid phase-5 WP runs.

Criterion: under a CAP the only energy sink is absorption, so a healthy run has
E_total(t) decreasing monotonically and flattening to a plateau (its minimum),
i.e. late-time dE/dt -> 0 from below. A positive late slope = energy created
(impossible under a CAP) = the grid-aliasing pathology.
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HA = 27.211386
E_GS = -70.22568  # Ha, 90-box localised-jellium GS (handover)
ROOT = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts")

RUNS = [
    ("v1.3  (23 eV)",  ROOT/"qsp_phase5/wp/results/p5_wp_v1p3/raw/observables/observables.csv", True),
    ("v2.0  (54 eV)",  ROOT/"qsp_phase4/wp/results/p4_wp/raw/observables/observables.csv",       True),
    ("v3.0  (122 eV)", ROOT/"qsp_phase5/wp/results/p5_wp_v3p0/raw/observables/observables.csv",  True),
    ("v4.0  (218 eV)", ROOT/"qsp_phase5/wp/results/p5_wp_v4p0/raw/observables/observables.csv",  True),
    ("v5.0  (340 eV)", ROOT/"qsp_phase5/wp/results/p5_wp_v5p0/raw/observables/observables.csv",  False),  # aliased
    ("v6.0  (490 eV)", ROOT/"qsp_phase5/wp/results/p5_wp_v6p0/raw/observables/observables.csv",  False),  # aliased
]

def late_slope(t, E, frac=0.20):
    n = len(t); k = max(5, int(frac*n))
    m = np.polyfit(t[-k:], E[-k:], 1)[0]          # Ha/au
    return m*HA                                    # eV/au

print(f"{'run':14s} {'valid':5s} {'E0-EGS':>8s} {'Emin-EGS':>9s} {'Ef-EGS':>8s} "
      f"{'t_min/t_f':>9s} {'rise_after':>10s} {'late dE/dt':>10s}  verdict")
print("-"*108)

fig, ax = plt.subplots(figsize=(7.0,4.4))
for label, path, valid in RUNS:
    df = pd.read_csv(path)
    t = df["time_au"].to_numpy(); E = df["energy_total"].to_numpy()
    dep = (E - E_GS)*HA                            # deposited energy, eV
    imin = int(np.argmin(E)); tf = t[-1]
    rise_after = (E[imin:].max() - E[imin])*HA     # eV created after the minimum
    sl = late_slope(t, E)
    # verdict
    if rise_after > 1.0 and t[imin] < 0.9*tf:
        v = "FAIL: energy created after minimum (unphysical)"
    elif t[imin] >= 0.97*tf and abs(sl) < 0.05:
        v = "CONVERGED: flat plateau at minimum"
    elif t[imin] >= 0.95*tf and sl < 0:
        v = "OK: monotonic drain, still flattening (truncated -> upper bound)"
    else:
        v = "ambiguous"
    print(f"{label:14s} {str(valid):5s} {(E[0]-E_GS)*HA:8.1f} {(E[imin]-E_GS)*HA:9.1f} "
          f"{(E[-1]-E_GS)*HA:8.1f} {t[imin]/tf:9.3f} {rise_after:10.2f} {sl:10.3f}  {v}")
    ls = "-" if valid else "--"
    lw = 2.0 if valid else 1.4
    ax.plot(t, dep, ls, lw=lw, label=label)

ax.axhline(0, color="0.6", lw=0.8, ls=":")
ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$E_{\rm total}(t) - E_{\rm GS}$  (eV)")
ax.set_title("Energy convergence: deposited energy vs time")
ax.legend(fontsize=8, ncol=2)
out = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/figs")
out.mkdir(exist_ok=True)
fig.tight_layout(); fig.savefig(out/"energy_convergence_check.png", dpi=130)
print(f"\nplot -> {out/'energy_convergence_check.png'}")
