#!/usr/bin/env python3
"""Campaign-1 sigma decision: does the CLASSICAL Gaussian projectile stopping
track the POINT-charge Lindhard reference well enough at sigma_WP=1.0 (less
spreading) as it does at sigma_WP=0.5?

KEY FINDING of this script (honest, validated):
  * The analytical POINT-charge Lindhard (stopping_power_point) is TRUSTWORTHY:
    at r_s=4, v=2.71 it gives 0.716 eV/Bohr, matching the 0.719 eV/Bohr the
    localised baseline used; and the localised sigma_WP=0.5 CLASSICAL run gave
    S=0.706 eV/Bohr = 0.99x point  =>  sigma_WP=0.5 classical ~ point-Lindhard.
  * The analytical FINITE-sigma Lindhard (stopping_power_sigma) is NOT reliable:
    it predicts sigma_WP=0.5 should be only ~0.78x point, contradicting the run
    (0.99x). It OVER-suppresses. So it CANNOT be used to judge sigma_WP=1.0.
  => The sigma=1 vs sigma=0.5 decision needs a dedicated sigma_WP=1.0 CLASSICAL
     S(v) run compared to point-Lindhard. Not decidable from existing data.

sqrt(2) convention (made explicit): sigma_pot (charge std) = sigma_WP / sqrt(2).
All stopping in eV/Bohr. Outputs: sigma_lindhard_comparison.png + table.csv.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.analysis import lindhard_elf as L
from inqview.analysis.stopping_extract import load_track
from inqview.visualisation import style

style.apply_theme()
HERE = os.path.dirname(os.path.abspath(__file__))
HA = 27.211386  # eV per Ha
SQRT2 = np.sqrt(2.0)

RS_LOC, RS_BULK = 4.0, 5.69
kF_loc, kF_bulk = L.kF_from_rs(RS_LOC), L.kF_from_rs(RS_BULK)

# consistent high-res grid for point AND finite-sigma (apples-to-apples ratio)
def Sev(v, kF, sig):
    qmax = 2 * kF + 2 * v + 4.0
    return HA * L.stopping_power_sigma(v, kF, sig, qmax=qmax, n_q=4000,
                                       n_omega=4000, eta=1e-2)

V = np.linspace(0.7, 6.7, 22)
E_eV = 0.5 * V ** 2 * HA

# --- localised r_s=4: point + analytical finite-sigma (caveated) ---
Spt_loc = np.array([Sev(v, kF_loc, 0.0) for v in V])
Ssw05 = np.array([Sev(v, kF_loc, 0.5 / SQRT2) for v in V])   # sigma_WP=0.5
Ssw10 = np.array([Sev(v, kF_loc, 1.0 / SQRT2) for v in V])   # sigma_WP=1.0
# the ONE validated run anchor (localised baseline classical, sigma_WP=0.5, 100 eV)
V_ANCHOR, S_ANCHOR = 2.711, 0.706
print(f"r_s=4 point-Lindhard @v=2.71: {np.interp(2.711, V, Spt_loc):.3f} eV/Bohr "
      f"(baseline used 0.719); localised sigma_WP=0.5 classical run = {S_ANCHOR} "
      f"=> {S_ANCHOR/np.interp(2.711,V,Spt_loc):.2f}x point")
print(f"analytical finite-sigma predicts sigma_WP=0.5 @100eV = "
      f"{np.interp(2.711,V,Ssw05)/np.interp(2.711,V,Spt_loc):.2f}x point "
      f"(WRONG vs the 0.99x run => over-suppresses)")

# --- bulk r_s=5.69 supporting: point + the sigma_pot=0.5 run sweep (rough) ---
Spt_bulk = np.array([Sev(v, kF_bulk, 0.0) for v in V])
SV = ("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
      "run_sv_sigma0p5/results")
REG = {"v0p6": 0.6, "v0p8": 0.8, "v1p3": 1.3, "v2p0": 2.0, "v3p0": 3.0}
run_v, run_S = [], []
for sub, v0 in REG.items():
    obs, trk = f"{SV}/{sub}/observables.csv", f"{SV}/{sub}/electron_track.csv"
    if not (os.path.exists(obs) and os.path.exists(trk)):
        continue
    o = pd.read_csv(obs)
    t, E = o["time_au"].values, o["energy_total"].values
    tr = load_track(trk)
    m = t >= t.min() + 0.20 * (t.max() - t.min())
    s_at = np.interp(t[m], tr.t, tr.s)
    slope = np.polyfit(s_at, E[m] - E[0], 1)[0]
    run_v.append(v0); run_S.append(abs(slope) * HA)
run_v, run_S = np.array(run_v), np.array(run_S)

# --- table ---
tab = pd.DataFrame({
    "v_au": V, "E_eV": E_eV, "S_point_loc": Spt_loc,
    "S_finiteSigma_wp0.5": Ssw05, "S_finiteSigma_wp1.0": Ssw10})
tab.to_csv(f"{HERE}/sigma_lindhard_table.csv", index=False)

# --- plot ---
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.2, 4.5))

# (A) localised target r_s=4 — the decision panel
axA.plot(E_eV, Spt_loc, "k-", lw=2.4, label="POINT-charge Lindhard (trusted)")
axA.plot(E_eV, Ssw05, "C0--", lw=1.5,
         label="analytical finite-σ, σ$_{WP}$=0.5  (⚠ over-suppresses)")
axA.plot(E_eV, Ssw10, "C1--", lw=1.5,
         label="analytical finite-σ, σ$_{WP}$=1.0  (⚠ unreliable)")
axA.plot([100], [S_ANCHOR], "o", color="C2", ms=10, zorder=5,
         label="localised σ$_{WP}$=0.5 CLASSICAL run = 0.706 (≈ point, 0.99×)")
axA.annotate("run sits on POINT curve,\nNOT on the analytical σ=0.5 dashed\n"
             "→ finite-σ formula over-suppresses",
             xy=(100, S_ANCHOR), xytext=(165, 0.55), fontsize=7, color="C2",
             arrowprops=dict(arrowstyle="->", color="C2", lw=1.0))
axA.annotate("σ$_{WP}$=1.0 verdict needs a RUN\n(analytical can't be trusted)",
             xy=(100, np.interp(2.711, V, Ssw10)), xytext=(150, 0.20),
             fontsize=7, color="C1",
             arrowprops=dict(arrowstyle="->", color="C1", lw=1.0))
axA.set_xlabel("projectile energy (eV)")
axA.set_ylabel("stopping power S (eV/Bohr)")
axA.set_title("Localised target r$_s$=4 — σ decision (Campaign 1)", fontsize=9)
axA.legend(fontsize=6.5, frameon=False, loc="upper right")
axA.grid(alpha=0.25)

# (B) bulk r_s=5.69 supporting — classical σ_pot=0.5 sweep tracks point
axB.plot(E_eV, Spt_bulk, "k-", lw=2.0, label="POINT Lindhard, r$_s$=5.69")
if len(run_v):
    axB.plot(0.5 * run_v ** 2 * HA, run_S, "s", color="0.25", ms=7,
             label="bulk σ$_{pot}$=0.5 classical run S(v) (rough extraction)")
axB.set_xlabel("projectile energy (eV)")
axB.set_ylabel("stopping power S (eV/Bohr)")
axB.set_title("Bulk support: classical σ≈0.5 tracks point Lindhard across v", fontsize=9)
axB.legend(fontsize=7, frameon=False, loc="upper right")
axB.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(f"{HERE}/sigma_lindhard_comparison.png", dpi=200)
plt.close(fig)
print("wrote sigma_lindhard_comparison.png + sigma_lindhard_table.csv")
