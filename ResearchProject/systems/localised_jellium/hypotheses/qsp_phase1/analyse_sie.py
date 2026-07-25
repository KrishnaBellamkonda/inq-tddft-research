#!/usr/bin/env python3
"""Phase-1 SIE diagnostic — compute SIE both ways + the zero-point cross-check.

Reads:
  GS:   scripts/qsp_phase1/gs/results/run_summary.txt   -> E_GS_slab (Ha)
  RT:   scripts/qsp_phase1/sie/results/sie_wp_far/raw/observables/
          observables.csv        -> energy_total(t=0)  (Ha)
          wp_momentum_stats.csv   -> e_kin_ha, pz_mean, px2/py2/pz2 (t=0)
          wp_real_space_stats.csv -> z_mean(0) (confirm launch ≈ −32)

SIE_a = E_tot(0) − (E_GS + 100 eV)        [user reference; = SIE + zero-point]
SIE_b = E_tot(0) − E_GS − KE_WP           [= SIE; KE_WP = ⟨p²⟩/2 measured]
cross-check: SIE_a − SIE_b = KE_WP − 100 eV  (≈ zero-point 3/(4σ²)=81.6 eV
             up to the drift discrepancy ½⟨pz⟩²−100).

Writes sie_results.csv in this directory + prints a summary.
"""
import os
import re
import numpy as np
import pandas as pd

HA = 27.211386
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
GS_SUMMARY = f"{ROOT}/scripts/qsp_phase1/gs/results/run_summary.txt"
RT = f"{ROOT}/scripts/qsp_phase1/sie/results/sie_wp_far/raw/observables"
HERE = os.path.dirname(os.path.abspath(__file__))
SIGMA = 0.5

# --- E_GS_slab ---
egs = None
with open(GS_SUMMARY) as f:
    for line in f:
        if line.startswith("ground_state_energy_ha"):
            egs = float(line.split("=")[1])
assert egs is not None, "E_GS not found"

# --- E_total(0) ---
obs = pd.read_csv(f"{RT}/observables.csv")
e0 = float(obs["energy_total"].iloc[0])

# --- KE_WP, drift, zero-point (t=0) ---
mom = pd.read_csv(f"{RT}/wp_momentum_stats.csv", comment="#")
ke = float(mom["e_kin_ha"].iloc[0])                 # ⟨p²⟩/2  (Ha)
pz = float(mom["pz_mean"].iloc[0])
drift = 0.5 * pz * pz                               # ½⟨pz⟩²  (Ha)
zp_measured = ke - drift                            # zero-point (Ha)

# --- launch confirm (norm from real-space stats = literal norm; momentum norm_check is scaled) ---
rsp = pd.read_csv(f"{RT}/wp_real_space_stats.csv", comment="#")
z0 = float(rsp["z_mean"].iloc[0])
norm0 = float(rsp["norm_check"].iloc[0]) if "norm_check" in rsp else np.nan

# --- SIE estimates ---
E_INJECT_EV = 100.0
SIE_a = (e0 - (egs + E_INJECT_EV / HA)) * HA        # user ref (eV)
SIE_b = (e0 - egs - ke) * HA                        # clean SIE (eV)
diff  = SIE_a - SIE_b                               # = KE_WP − 100 eV
zp_theory = 3.0 / (4.0 * SIGMA ** 2) * HA           # 3/(4σ²) in eV = 81.6

print(f"E_GS_slab        = {egs:.6f} Ha = {egs*HA:.2f} eV")
print(f"E_total(0)[WP far]= {e0:.6f} Ha = {e0*HA:.2f} eV")
print(f"launch z_mean(0)  = {z0:.3f} Bohr  (target −32)   WP norm0 = {norm0:.4f}")
print(f"KE_WP (⟨p²⟩/2)   = {ke:.6f} Ha = {ke*HA:.2f} eV")
print(f"  drift ½⟨pz⟩²    = {drift*HA:.2f} eV   (pz_mean={pz:.4f})")
print(f"  zero-point      = {zp_measured*HA:.2f} eV  (measured)   vs 3/4σ²={zp_theory:.2f} eV (theory)")
print("-" * 56)
print(f"SIE_a (E_GS+100eV ref) = {SIE_a:+.2f} eV   [= SIE + zero-point]")
print(f"SIE_b (clean SIE)      = {SIE_b:+.2f} eV   [= SIE]")
print(f"SIE_a − SIE_b          = {diff:+.2f} eV   (= KE_WP−100; expect ≈ zero-point {zp_measured*HA:.1f})")
print(f"=> THE SIE ≈ {SIE_b:.2f} eV  (cf. old r_s=4 value ≈ 4.5 eV)")

pd.DataFrame([{
    "E_GS_slab_Ha": egs, "E_total0_Ha": e0, "launch_z": z0, "wp_norm0": norm0,
    "KE_WP_eV": ke*HA, "drift_eV": drift*HA, "zero_point_meas_eV": zp_measured*HA,
    "zero_point_theory_eV": zp_theory, "SIE_a_eV": SIE_a, "SIE_b_eV": SIE_b,
    "SIE_a_minus_b_eV": diff,
}]).to_csv(f"{HERE}/sie_results.csv", index=False)
print(f"\nwrote {HERE}/sie_results.csv")
