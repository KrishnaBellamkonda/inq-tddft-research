#!/usr/bin/env python3
"""Phase 2 PROPER deposited-energy ledger (replaces the invalid centroid proxy).

WP deposited energy (per-step CSVs only, no field data needed):
    dE_dep = E_wp(0) - [ INT e_kin(t) * (-d norm_wp) + e_kin(end) * norm_wp(end) ]
i.e. initial packet KE minus KE carried away — through the CAP (each absorbed
norm parcel carries the instantaneous norm-weighted mean KE e_kin_ha) and by
the surviving fraction. CAVEAT (recorded): the CAP absorbs the fast transmitted
front first, so e_kin(t) OVERestimates the KE of late (slow) absorbed flux —
dE_dep is therefore a LOWER bound on deposition; the front-bias is worst where
reflection is strong (E25). Cross-checks printed: total CAP removal from
energy_total (E_tot(0)-E_tot(end)) and bath norm loss (secondary emission).

Classical: dE = KE(-8) - KE(+8) from electron_track (image-force-cancelling
window), as in build_phase2_ratio.py.
"""
from pathlib import Path
import csv

import numpy as np

HERE = Path(__file__).parent
SYS = HERE.parent.parent
HA_EV = 27.211386245988
ENERGIES = [25, 100, 300]
W_HALF = 8.0


def wp_ledger(ev):
    obs = SYS / f"scripts/twodef_sv/wp/results/E{ev}/raw/observables"
    mom = np.genfromtxt(obs / "wp_momentum_stats.csv", delimiter=",", names=True, skip_header=1)
    ix = np.genfromtxt(obs / "interactions.csv", delimiter=",", names=True)
    n = min(len(mom["e_kin_ha"]), len(ix["norm_wp"]))
    ekin, nw = mom["e_kin_ha"][:n], ix["norm_wp"][:n]
    obs_csv = np.genfromtxt(obs / "observables.csv", delimiter=",", names=True)
    etot = obs_csv["energy_total"][: n]

    # energy carried out through the CAP by WP flux: sum e_kin * (-d norm_wp)
    dn = -np.diff(nw)
    carried_cap = float(np.sum(0.5 * (ekin[1:] + ekin[:-1]) * dn))
    carried_surv = float(ekin[-1] * nw[-1])
    e0 = float(ekin[0])
    de_dep = (e0 - carried_cap - carried_surv) * HA_EV

    cap_total = float((etot[0] - etot[-1]) * HA_EV)      # ALL CAP removal (WP + secondaries)
    bath_norm_loss = float((ix["norm_total"][0] - nw[0]) - (ix["norm_total"][n-1] - nw[n-1]))
    return de_dep, e0 * HA_EV, carried_cap * HA_EV, carried_surv * HA_EV, cap_total, bath_norm_loss, float(nw[n-1])


def classical_de(ev):
    t = np.genfromtxt(SYS / f"scripts/twodef_sv/classical/results/E{ev}/raw/observables/electron_track.csv",
                      delimiter=",", names=True)
    ke_in = np.interp(-W_HALF, t["z"], t["ke_ion_ha"])
    ke_out = np.interp(+W_HALF, t["z"], t["ke_ion_ha"])
    return (ke_in - ke_out) * HA_EV


def main():
    rows = []
    print(f"{'E(eV)':>6} {'dE_cl':>7} {'dE_wp':>7} {'ratio':>6} | {'E_wp0':>7} {'out_cap':>8} {'out_surv':>8} {'capTOT':>7} {'d_bathN':>8} {'nw_end':>7}")
    for ev in ENERGIES:
        de_c = classical_de(ev)
        de_w, e0, ccap, csurv, captot, dbn, nwend = wp_ledger(ev)
        r = de_w / de_c if de_c else float("nan")
        rows.append((ev, de_c, de_w, r, e0, ccap, csurv, captot, dbn, nwend))
        print(f"{ev:>6} {de_c:>7.2g} {de_w:>7.2g} {r:>6.2g} | {e0:>7.3g} {ccap:>8.3g} {csurv:>8.3g} {captot:>7.3g} {dbn:>8.3g} {nwend:>7.3g}")

    with open(HERE / "phase2_ledger.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_eV", "dE_classical_eV", "dE_wp_deposited_eV", "ratio",
                    "E_wp0_eV", "carried_out_cap_eV", "carried_surviving_eV",
                    "cap_total_removed_eV", "bath_norm_loss", "norm_wp_end"])
        w.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from inqview.visualisation import style
        style.apply()
    except Exception:
        pass
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.4))
    evs = [r[0] for r in rows]
    ax1.plot(evs, [r[1] for r in rows], "o-", label="classical (±8 Bohr window)")
    ax1.plot(evs, [r[2] for r in rows], "s-", label="WP (flux ledger, lower bound)")
    ax1.set_xlabel("E (eV)"); ax1.set_ylabel(r"$\Delta E$ deposited (eV)")
    ax1.set_xscale("log"); ax1.legend(frameon=False, fontsize=8)
    ax2.plot(evs, [r[3] for r in rows], "d-", color="k")
    ax2.axhline(1.0, ls=":", lw=1); ax2.axhline(2.2, ls="--", lw=1, label="bulk jellium 2.2")
    ax2.set_xlabel("E (eV)"); ax2.set_ylabel("WP / classical")
    ax2.set_xscale("log"); ax2.legend(frameon=False, fontsize=8)
    fig.suptitle("Phase 2 twins — deposited-energy ledger (bilayer graphene)", fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "phase2_ledger.png", dpi=180)
    print(f"\nwrote {HERE/'phase2_ledger.csv'} and phase2_ledger.png")


if __name__ == "__main__":
    main()
