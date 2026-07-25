#!/usr/bin/env python3
"""Q3 + Q8: stopping-power energy ledger from E_total (WP run), two baselines.

User-agreed method (2026-06-23): work with TOTAL energy only (the WP and bath cannot
be cleanly separated energetically). Deposited energy two ways:
  Formula 1 (subtract the 100 eV drift):  E_total(final) - E_total(0) + 100 eV
  Formula 2 (slab-GS baseline, exact):    E_total(final) - E_GS_slab
They are equal ONLY IF the WP's added energy is exactly 100 eV; their gap = the WP's
zero-point + self/interaction energy above the drift.  S = deposited / 25 Bohr.
Q8 cross-check: E_total(0)_WP - E_total(0)_classical (expected ~100 eV if the runs were
otherwise identical).

Sources: raw/observables/observables.csv (energy_total); E_GS from the slab GS run.
Outputs: qa_v_stopping_energy.png + printed ledger.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inqview.visualisation import style

style.apply_theme()
HA_EV = 27.211386
E_GS = -160.9920712992004        # slab ground-state total energy (no WP), Ha
E_DRIFT_HA = 0.5 * (2.7110633401 ** 2)   # 100 eV WP drift KE = 3.675 Ha
X_BOHR = 25.0
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"

ewp = pd.read_csv(f"{WP}/observables/observables.csv")
ecl = pd.read_csv(f"{CL}/observables/observables.csv")
E0w, Efw = ewp.energy_total.iloc[0], ewp.energy_total.iloc[-1]
E0c = ecl.energy_total.iloc[0]

# --- two-baseline deposited energy (WP run) ---
dep_f1 = (Efw - E0w) + E_DRIFT_HA          # subtract the 100 eV drift
dep_f2 = Efw - E_GS                          # slab-GS baseline (exact at full absorption)
gap = (E0w - E_GS) - E_DRIFT_HA              # WP zero-point + interaction above the drift
S1, S2 = dep_f1 / X_BOHR, dep_f2 / X_BOHR

# analytic WP zero-point KE (psi ~ exp(-r^2/2 sigma^2), sigma=0.5): 3/(4 sigma^2)
zpe = 3.0 / (4.0 * 0.5 ** 2)

print("=== Q3 — WP stopping via total energy (two baselines) ===")
print(f"E_total(0)_WP   = {E0w:.4f} Ha     E_total(final)_WP = {Efw:.4f} Ha")
print(f"E_GS_slab       = {E_GS:.4f} Ha     drift(100 eV)     = {E_DRIFT_HA:.4f} Ha")
print(f"Formula 1 (final-init+100eV): deposited = {dep_f1*HA_EV:7.1f} eV  -> S = {S1*HA_EV:.2f} eV/Bohr")
print(f"Formula 2 (final - E_GS)    : deposited = {dep_f2*HA_EV:7.1f} eV  -> S = {S2*HA_EV:.2f} eV/Bohr")
print(f"GAP (F2-F1) = WP energy above 100 eV drift = {gap*HA_EV:.1f} eV "
      f"(analytic zero-point KE alone = {zpe*HA_EV:.1f} eV)")
print("=== Q8 — initial total-energy comparison ===")
print(f"E_total(0)_WP - E_total(0)_classical = {(E0w-E0c)*HA_EV:.1f} eV  (NOT ~100 eV)")
print(f"E_total(0)_WP - E_GS                 = {(E0w-E_GS)*HA_EV:.1f} eV  (WP adds drift+zero-point)")
print(f"E_total(0)_cl - E_GS                 = {(E0c-E_GS)*HA_EV:.1f} eV  (ghost potential on unrelaxed GS)")

# --- residual-WP energy contamination of Formula 2 (Q7 test) ---
# At full absorption Ef = E_GS + deposited. At 62% absorption Ef ALSO holds the residual
# WP energy. So Formula 2 = deposited + E_residual_WP. Rough scale of the residual energy:
N_wp_final = 0.378
E_wp_per_norm = (E0w - E_GS)          # the WP's t=0 energy contribution per unit norm (Ha)
E_resid_est = N_wp_final * E_wp_per_norm
print(f"[Q7] Formula 2 = deposited + residual-WP energy. Residual WP norm {N_wp_final}, "
      f"rough residual energy ~{E_resid_est*HA_EV:.0f} eV (linear-in-norm estimate) -> "
      f"Formula 2 ({dep_f2*HA_EV:.0f} eV) is dominated by un-absorbed WP, NOT deposited.")

# --- figure: energy change from t=0 (+ the +100 eV 'slab-remaining' curve, Q9) ---
fig, ax = plt.subplots(figsize=(6.6, 4.4))
dEw = (ewp.energy_total - E0w) * HA_EV
dEc = (ecl.energy_total - E0c) * HA_EV
ax.plot(ewp.time_au, dEw, "C0-", lw=1.8,
        label=f"WP: E_total(t)-E_total(0)  (ends {(Efw-E0w)*HA_EV:.0f} eV — CAP removes WP)")
ax.plot(ewp.time_au, dEw + 100.0, "C0--", lw=1.6,
        label="WP: [E_total(t)-E_total(0)] + 100 eV  (drift credited -> Formula-1 'slab remaining')")
ax.plot(ecl.time_au, dEc, "C1-", lw=1.8,
        label=f"classical: E_total(t)-E_total(0)  (ends {(ecl.energy_total.iloc[-1]-E0c)*HA_EV:.0f} eV)")
ax.axhline(0, color="grey", lw=0.6)
ax.axhline(dep_f1 * HA_EV, color="C0", ls=":", lw=0.8)
ax.annotate(f"Formula-1 deposited {dep_f1*HA_EV:.0f} eV", xy=(1, dep_f1*HA_EV+3),
            fontsize=6.5, color="C0")
ax.set_xlabel("time (a.u.)"); ax.set_ylabel("energy change (eV)")
ax.set_title("Total-energy change vs time — WP (+100 eV drift-credited) vs classical", fontsize=9)
ax.legend(fontsize=6.5, frameon=False, loc="center left"); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/qa_v_stopping_energy.png", dpi=200); plt.close(fig)
print("wrote qa_v_stopping_energy.png")
