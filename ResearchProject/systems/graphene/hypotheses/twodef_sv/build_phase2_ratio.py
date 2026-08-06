#!/usr/bin/env python3
"""Phase 2 twin-sanity ratio: WP vs classical energy deposit, bilayer graphene.

Plan: docs/plans/real-material-stopping-comparison.md (Phase 2 gate: the ratio
is MEASURED and documented, not asserted; hypothesis WP/classical ~ 1, bulk
jellium gave ~2.2).

Method (matched z-window, image-force-cancelling):
  ΔE = KE(z=-W_HALF) - KE(z=+W_HALF) for each projectile, W_HALF = 8 Bohr —
  symmetric around the bilayer (layers at z=±3.165), so the conservative
  image/core attraction gained on approach is repaid on exit and the window
  difference isolates DEPOSITED energy.
  * classical: KE from electron_track.csv (ke_ion_ha at first crossing of ∓W).
  * WP: KE proxy = <p_z>²/2 at the centroid's crossings (wp_real_space_stats
    z_mean + wp_momentum_stats pz). CAVEAT (qsp5 lesson, recorded): after CAP
    contact the norm-weighted moments are biased; both crossings happen BEFORE
    significant absorption (checked via norm column), so the bias is bounded
    by the printed norm drop across the window.
Output: phase2_ratio.csv + phase2_ratio.png + printed 2-s.f. table.
"""
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SYS = HERE.parent.parent
W_HALF = 8.0
ENERGIES = [25, 100, 300]
HA_EV = 27.211386245988


def crossing_ke_classical(track_csv, w=W_HALF):
    rows = np.genfromtxt(track_csv, delimiter=",", names=True)
    z, ke, norm_note = rows["z"], rows["ke_ion_ha"], None
    ke_in = np.interp(-w, z, ke)   # z monotonic pre-stall for these runs
    ke_out = np.interp(+w, z, ke)
    return ke_in, ke_out, norm_note


def crossing_ke_wp(obs_dir, w=W_HALF):
    # stats CSVs carry a '# wp_state_index=...' comment line before the header
    pos = np.genfromtxt(obs_dir / "wp_real_space_stats.csv", delimiter=",",
                        names=True, skip_header=1)
    mom = np.genfromtxt(obs_dir / "wp_momentum_stats.csv", delimiter=",",
                        names=True, skip_header=1)
    zc = pos["z_mean_circ"]
    n = min(len(zc), len(mom["pz_mean"]))
    zc, pz = zc[:n], mom["pz_mean"][:n]
    norm = pos["norm_check"][:n]
    # first crossings while centroid is monotone through the slab
    def at(zq):
        idx = np.argmax(zc >= zq)
        if idx == 0 and zc[0] < zq:
            return None, None
        return 0.5 * pz[idx] ** 2, norm[idx]
    ke_in, n_in = at(-w)
    ke_out, n_out = at(+w)
    return ke_in, ke_out, (n_in, n_out)


def main():
    out_rows = []
    for ev in ENERGIES:
        cl_track = SYS / f"scripts/twodef_sv/classical/results/E{ev}/raw/observables/electron_track.csv"
        wp_obs = SYS / f"scripts/twodef_sv/wp/results/E{ev}/raw/observables"
        ke_in_c, ke_out_c, _ = crossing_ke_classical(cl_track)
        de_c = (ke_in_c - ke_out_c) * HA_EV
        ke_in_w, ke_out_w, norms = crossing_ke_wp(wp_obs)
        if ke_out_w is None:
            de_w, ratio, normnote = float("nan"), float("nan"), "WP centroid never crossed +W (backscatter/absorbed)"
        else:
            de_w = (ke_in_w - ke_out_w) * HA_EV
            ratio = de_w / de_c if de_c else float("nan")
            normnote = f"norm {norms[0]:.3f}->{norms[1]:.3f} across window"
        out_rows.append((ev, de_c, de_w, ratio, normnote))

    with open(HERE / "phase2_ratio.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_eV", "dE_classical_eV", "dE_wp_eV", "ratio_wp_over_classical", "note"])
        w.writerows(out_rows)

    print(f"{'E (eV)':>7} {'dE_cl (eV)':>11} {'dE_wp (eV)':>11} {'ratio':>7}  note")
    for ev, dc, dw, r, note in out_rows:
        print(f"{ev:>7} {dc:>11.2g} {dw:>11.2g} {r:>7.2g}  {note}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from inqview.visualisation import style
        style.apply()
    except Exception:
        pass
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.4))
    evs = [r[0] for r in out_rows]
    ax1.plot(evs, [r[1] for r in out_rows], "o-", label="classical")
    ax1.plot(evs, [r[2] for r in out_rows], "s-", label="wavepacket")
    ax1.set_xlabel("E (eV)"); ax1.set_ylabel(r"$\Delta E$ per bilayer transit (eV)")
    ax1.set_xscale("log"); ax1.legend(frameon=False)
    ax2.plot(evs, [r[3] for r in out_rows], "d-", color="k")
    ax2.axhline(1.0, ls=":", lw=1)
    ax2.axhline(2.2, ls="--", lw=1, label="bulk jellium (2.2)")
    ax2.set_xlabel("E (eV)"); ax2.set_ylabel("WP / classical")
    ax2.set_xscale("log"); ax2.legend(frameon=False)
    fig.suptitle("Phase 2 twin sanity — bilayer graphene, matched ±8 Bohr window", fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "phase2_ratio.png", dpi=180)
    print(f"\nwrote {HERE/'phase2_ratio.csv'} and phase2_ratio.png")


if __name__ == "__main__":
    main()
