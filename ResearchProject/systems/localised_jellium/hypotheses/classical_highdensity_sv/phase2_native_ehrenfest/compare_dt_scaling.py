#!/usr/bin/env python3
"""dt-scaling check for Phase 2 Test C.

If the native-vs-perturbation gap is dominated by the O(dt) intra-step ordering,
halving dt should shrink max|Δz| and max|Δvz|. Compares the (native, pert) gap at
dt=0.02 vs dt=0.01 on the SAME physical trajectory (interpolated onto a common
time grid so the two step counts are comparable).
"""
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPTS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
               "scripts/classical_highdensity_sv/phase2_native_ehrenfest")


def gap(nat_csv, per_csv):
    nat = pd.read_csv(nat_csv); per = pd.read_csv(per_csv)
    # common time grid (both share identical time stamps step*dt)
    n = min(len(nat), len(per))
    nat, per = nat.iloc[:n], per.iloc[:n]
    assert np.allclose(nat.time.values, per.time.values)
    dz = np.abs(nat.z.values - per.z.values)
    dvz = np.abs(nat.vz.values - per.vz.values)
    return float(dz.max()), float(dvz.max()), len(nat)


def main():
    c1 = SCRIPTS / "c1_native"
    c2 = SCRIPTS / "c2_pert"
    dz2, dvz2, n2 = gap(c1 / "results/native_dt0p02.csv", c2 / "results/pert_dt0p02.csv")
    dz1, dvz1, n1 = gap(c1 / "results_dthalf/native.csv", c2 / "results_dthalf/pert.csv")

    print(f"dt=0.02 ({n2} steps): max|dz|={dz2:.6e}  max|dvz|={dvz2:.6e}")
    print(f"dt=0.01 ({n1} steps): max|dz|={dz1:.6e}  max|dvz|={dvz1:.6e}")
    print(f"ratio dz(0.02)/dz(0.01)   = {dz2/dz1:.3f}   (≈2 ⇒ O(dt))")
    print(f"ratio dvz(0.02)/dvz(0.01) = {dvz2/dvz1:.3f}   (≈2 ⇒ O(dt))")


if __name__ == "__main__":
    main()
