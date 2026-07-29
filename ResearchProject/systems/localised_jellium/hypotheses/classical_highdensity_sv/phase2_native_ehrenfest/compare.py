#!/usr/bin/env python3
"""Phase 2 Test C — compare INQ native Ehrenfest (C1) vs perturbation (C2).

Overlays z(t) and vz(t), reports max|Δz|, max|Δvz| over the common trajectory,
and checks whether the electronic energy exchange matches. Writes a PNG and prints
a machine-readable summary consumed by native_ehrenfest_comparison.md.

Both CSVs share columns:
  step,time,z,vz,E_elec,E_total,E_kin,E_hartree,E_external,E_nonlocal,E_xc,E_ion
Row i (step i) records the projectile position AFTER i advances, so the two runs
are aligned step-for-step.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
try:
    from inqview.visualisation.style import apply_theme
    apply_theme()
except Exception as e:
    print(f"(style theme not applied: {e})")

HERE = Path(__file__).resolve().parent
SCRIPTS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
               "scripts/classical_highdensity_sv/phase2_native_ehrenfest")


def load(tag, csv):
    df = pd.read_csv(csv)
    print(f"{tag}: {len(df)} rows, z {df.z.iloc[0]:.4f}->{df.z.iloc[-1]:.4f}, "
          f"vz {df.vz.iloc[0]:.4f}->{df.vz.iloc[-1]:.4f}")
    return df


def main():
    nat_csv = SCRIPTS / "c1_native/results/native.csv"
    per_csv = SCRIPTS / "c2_pert/results/pert.csv"
    nat = load("NATIVE (C1)", nat_csv)
    per = load("PERT   (C2)", per_csv)

    # align on common step range
    n = min(len(nat), len(per))
    nat, per = nat.iloc[:n].reset_index(drop=True), per.iloc[:n].reset_index(drop=True)
    assert np.allclose(nat.step, per.step), "step misalignment"
    t = nat.time.values

    dz = np.abs(nat.z.values - per.z.values)
    dvz = np.abs(nat.vz.values - per.vz.values)
    max_dz, max_dvz = float(dz.max()), float(dvz.max())
    # energy exchange: electronic energy relative to its t=0 value
    dE_nat = nat.E_elec.values - nat.E_elec.values[0]
    dE_per = per.E_elec.values - per.E_elec.values[0]
    max_dE = float(np.abs(dE_nat - dE_per).max())
    HA = 27.211386

    # ghost-moved gate
    z_travel_nat = float(nat.z.values[-1] - nat.z.values[0])
    z_travel_per = float(per.z.values[-1] - per.z.values[0])

    # ---- plot ----
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0))
    ax[0].plot(t, nat.z, label="native (C1)", lw=1.8)
    ax[0].plot(t, per.z, "--", label="perturbation (C2)", lw=1.6)
    ax[0].set_xlabel("time (a.u.)"); ax[0].set_ylabel("z (Bohr)")
    ax[0].set_title(f"z(t)  max|Δz|={max_dz:.2e} Bohr"); ax[0].legend()

    ax[1].plot(t, nat.vz, label="native (C1)", lw=1.8)
    ax[1].plot(t, per.vz, "--", label="perturbation (C2)", lw=1.6)
    ax[1].set_xlabel("time (a.u.)"); ax[1].set_ylabel("vz (a.u.)")
    ax[1].set_title(f"vz(t)  max|Δvz|={max_dvz:.2e}"); ax[1].legend()

    ax[2].plot(t, dE_nat * HA, label="native ΔE_elec", lw=1.8)
    ax[2].plot(t, dE_per * HA, "--", label="pert ΔE_elec", lw=1.6)
    ax[2].set_xlabel("time (a.u.)"); ax[2].set_ylabel("ΔE_elec (eV)")
    ax[2].set_title(f"electronic energy exchange\nmax|ΔΔE|={max_dE*HA:.2e} eV"); ax[2].legend()

    fig.tight_layout()
    out = HERE / "native_ehrenfest_comparison.png"
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}")

    print("\n===== SUMMARY =====")
    print(f"n_common_steps   = {n}")
    print(f"max|dz|          = {max_dz:.6e} Bohr")
    print(f"max|dvz|         = {max_dvz:.6e} a.u.")
    print(f"max|ddE_elec|    = {max_dE:.6e} Ha = {max_dE*HA:.6e} eV")
    print(f"z_travel_native  = {z_travel_nat:.4f} Bohr (ghost moved: {abs(z_travel_nat)>1e-3})")
    print(f"z_travel_pert    = {z_travel_per:.4f} Bohr")
    print(f"final z: native={nat.z.values[-1]:.5f}  pert={per.z.values[-1]:.5f}")
    print(f"final vz: native={nat.vz.values[-1]:.5f}  pert={per.vz.values[-1]:.5f}")
    # relative scale
    zspan = float(np.ptp(nat.z.values))
    vspan = float(np.ptp(nat.vz.values)) if np.ptp(nat.vz.values) > 0 else 1.0
    print(f"max|dz|/z_span   = {max_dz/max(zspan,1e-30):.3e}")
    print(f"max|dvz|/|v0|    = {max_dvz/abs(nat.vz.values[0]):.3e}")


if __name__ == "__main__":
    main()
