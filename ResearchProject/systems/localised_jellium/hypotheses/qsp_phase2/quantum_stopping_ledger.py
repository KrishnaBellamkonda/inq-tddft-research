#!/usr/bin/env python3
"""Quantum (wavepacket) stopping-power energy ledger — qsp_phase2 (26-6-26 meeting).

Computes the QUANTUM electronic stopping power at E=100 eV from the p2_wp / p2_classical
twin runs (tau=40 a.u., 50x50x70 box, sigma_WP=0.5) via the user's RETAINED-ENERGY
(bath) definition:

    S_WP = [E_total(t_f) - E_jellium(0)] / L_z ,   L_z = 25 Bohr (slab thickness)
    E_jellium(0) = E_total(0) - <T_WP> - E_SIE     (strip WP kinetic + self-interaction)
    E_total(t_f) = E_jellium(t_f)   (valid once the CAP has absorbed the WP remnants)

and the curiosity check E_jellium(0) vs E_GS (bare-slab ground state).

Classical slab stopping (p2_classical) is the Ehrenfest ion DELTA-KE, reported TWO ways:
  * slab-centre KE minimum (first transit) = the user's "lowest energy point" estimate;
  * equal-potential-face loss (|z|=12.5) = the conservative-well-corrected value.
The classical E_total energy method does NOT apply (lowest E_total is t=0; the ion-bath
interaction swamps the deposited energy; the ion wraps and re-enters by t_f).

Sources: campaign docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md;
handover docs/handovers/localised-jellium.md (P1.1 E_GS, P1.3 SIE, retained-energy def).

Numbers presented to 2 s.f. (3 s.f. for near-equal differences) per
.claude/rules/number-rounding.md. Full precision kept in the returned dict.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/local/data/public/skcb2/tddft")
sys.path.insert(0, str(ROOT / "inq-stack/python"))
from inqview.visualisation import style  # noqa: E402

HA_TO_EV = 27.211386245988
L_Z_BOHR = 25.0          # slab thickness |z|<12.5 -> 25 Bohr (the stopping path length)
SLAB_FACE = 12.5         # |z| of the slab faces (same geometry for both runs)
_SCR = ROOT / "ResearchProject/systems/localised_jellium/scripts"

# Per-run configuration. Both are sigma_WP=0.5, E=100 eV (k0=2.711) WP/classical
# twins with a slab |z|<12.5; they differ in slab density and run length.
RUNS = {
    # p2: the headline r_s=5.67 (N=82) twin, tau=40 a.u.
    "p2": dict(
        wp_obs=_SCR / "qsp_phase2/wp/results/p2_wp/raw/observables",
        cl_obs=_SCR / "qsp_phase2/classical/results/p2_classical/raw/observables",
        E_GS_HA=-45.75885,    # shared_gs/slab_n82_L50x50x70 (handover P1.1)
        E_SIE_EV=4.40,        # locked SIE_b (handover P1.3)
        n_bath=82, r_s=5.667, tau_au=40, label="p2 (N=82, r_s=5.67, tau=40)",
    ),
    # p5: the older r_s=4 (N=234) twin, tau=18 a.u. ("the 14-20 a.u. run")
    "p5": dict(
        wp_obs=_SCR / "fullsuite_wp/results/p5_wp/raw/observables",
        cl_obs=_SCR / "fullsuite_classical/results/p5_classical/raw/observables",
        E_GS_HA=-160.99207,   # static-run log, shared_gs/slab_n234_L50
        E_SIE_EV=4.50,        # r_s=4 SIE (handover; gives E_jellium(0)-E_GS~+0.5 eV)
        n_bath=234, r_s=3.995, tau_au=18, label="p5 (N=234, r_s=4.0, tau=18)",
    ),
    # p3: the BIG-BOX production twin (50x50x90, tau=100 a.u.) — same r_s=5.67/N=82
    # density as p2, built for the energy method (longer time, full absorption).
    "p3": dict(
        wp_obs=_SCR / "qsp_phase3/wp/results/p3_wp/raw/observables",
        cl_obs=_SCR / "qsp_phase3/classical/results/p3_classical/raw/observables",
        E_GS_HA=-70.22568216820937,  # shared_gs/slab_n82_L50x50x90 run_summary
        E_SIE_EV=4.40,               # same r_s=5.67 slab as p2
        n_bath=82, r_s=5.667, tau_au=100,
        label="p3 (N=82, r_s=5.67, tau=100, big box 50x50x90)",
    ),
}


def cfg(run="p2"):
    return RUNS[run]


def _read_obs(obs_dir):
    """observables.csv -> DataFrame (step,time_au,energy_total,...)."""
    return pd.read_csv(obs_dir / "observables.csv")


def _read_wp_ekin0(obs_dir):
    """<T_WP>(0): the run-measured WP kinetic energy at step 0 (e_kin_ha column).
    The file's first line is a '# wp_state_index=.. write_every=..' comment."""
    df = pd.read_csv(obs_dir / "wp_momentum_stats.csv", comment="#")
    row0 = df.loc[df["step"] == 0].iloc[0]
    return float(row0["e_kin_ha"])


def compute_wp_ledger(run="p2"):
    """The quantum energy ledger. Returns a dict of Ha + eV values (full precision).

    Headline S uses the FULL-LEDGER, full-absorption assumption (user decision
    2026-06-26): S = (E_total(t_f) - E_GS)/L_z, reported as an UPPER BOUND because
    the sigma=0.5 packet's 82 eV zero-point energy inflates the bath gain and the
    WP is not fully absorbed. The drift-credit alternative is also returned for the
    transparency table but is NOT the headline (it gives an impossible negative for
    the more-absorbed run)."""
    c = cfg(run)
    WP_OBS = c["wp_obs"]
    E_GS_HA, E_SIE_EV = c["E_GS_HA"], c["E_SIE_EV"]
    obs = _read_obs(WP_OBS)
    t = obs["time_au"].to_numpy()
    Etot = obs["energy_total"].to_numpy()
    E0 = float(Etot[0])
    Ef = float(Etot[-1])
    tf = float(t[-1])
    T_wp = _read_wp_ekin0(WP_OBS)
    E_sie_ha = E_SIE_EV / HA_TO_EV
    T_drift_ev = 100.0                        # 1/2 k0^2 for k0=2.711 (E=100 eV)
    T_zp_ev = T_wp * HA_TO_EV - T_drift_ev    # zero-point = full <T_WP> - drift

    E_jell0 = E0 - T_wp - E_sie_ha            # reconstructed initial bath energy
    # FULL LEDGER (headline): bath gain referenced to the bare-slab GS.
    dE_ha = Ef - E_GS_HA                       # = E_total(t_f) - E_GS
    S_wp = dE_ha * HA_TO_EV / L_Z_BOHR         # eV/Bohr (UPPER BOUND)
    # drift-credit alternative (zero-point-persists assumption) for the table only:
    dE_driftcredit_ev = (Ef - E0) * HA_TO_EV + T_drift_ev + E_SIE_EV
    S_driftcredit = dE_driftcredit_ev / L_Z_BOHR

    # --- convergence gate diagnostics ---
    # absorbed norm from electron_number.csv (N_total(0) - N_total(t)); WP not fully
    # absorbed at tau=40 -> S is an UPPER BOUND.
    en = pd.read_csv(WP_OBS / "electron_number.csv")
    ncol = [c for c in en.columns if c not in ("step", "time_au")][-1]
    N = en[ncol].to_numpy()
    N_time = en["time_au"].to_numpy()
    n_absorbed = float(N[0] - N[-1])
    wp_norm_remaining = 1.0 - n_absorbed     # ~1 WP electron launched
    # late E_total slope (last 25% of time) — plateau test
    m = t >= (0.75 * tf)
    slope_ha_per_au = float(np.polyfit(t[m], Etot[m], 1)[0])
    slope_ev_per_au = slope_ha_per_au * HA_TO_EV
    gate_met = (wp_norm_remaining < 0.02) and (abs(slope_ev_per_au) < 0.05)

    return dict(
        run=run, label=c["label"], r_s=c["r_s"],
        E_total_0_ha=E0, E_total_f_ha=Ef, t_f_au=tf,
        T_wp_ha=T_wp, T_drift_ev=T_drift_ev, T_zp_ev=T_zp_ev,
        E_sie_ev=E_SIE_EV, E_sie_ha=E_sie_ha,
        E_jellium_0_ha=E_jell0, E_GS_ha=E_GS_HA,
        E_jellium0_minus_GS_ha=E_jell0 - E_GS_HA,
        E_jellium0_minus_GS_ev=(E_jell0 - E_GS_HA) * HA_TO_EV,
        dE_ha=dE_ha, dE_ev=dE_ha * HA_TO_EV,
        S_wp_ev_per_bohr=S_wp, S_is_upper_bound=not gate_met,
        dE_driftcredit_ev=dE_driftcredit_ev, S_driftcredit_ev_per_bohr=S_driftcredit,
        n_absorbed=n_absorbed, wp_norm_remaining=wp_norm_remaining,
        slope_ev_per_au=slope_ev_per_au, gate_met=gate_met,
        time_au=t, E_total_ha=Etot, N_total=N, N_time_au=N_time,
    )


def compute_classical_slab(run="p2"):
    """Classical Ehrenfest-ion stopping through the slab, two ways.
    electron_track.csv columns: step,time_au,x,y,z,vx,vy,vz,ke_ion_ha."""
    CL_OBS = cfg(run)["cl_obs"]
    trk = pd.read_csv(CL_OBS / "electron_track.csv")
    trk = trk.drop_duplicates(subset="step", keep="last").reset_index(drop=True)
    t = trk["time_au"].to_numpy()
    z = trk["z"].to_numpy()
    ke_ev = trk["ke_ion_ha"].to_numpy() * HA_TO_EV
    z_launch = float(z[0])
    ke_launch = float(ke_ev[0])

    # FIRST transit only: from launch until the ion first exits the slab (z > +12.5).
    exit_idx = np.argmax(z > SLAB_FACE) if np.any(z > SLAB_FACE) else len(z) - 1
    first = slice(0, exit_idx + 1)
    zc = z[first]
    kec = ke_ev[first]
    i_min = int(np.argmin(kec))           # lowest KE during first transit = slab centre
    z_center = float(zc[i_min])
    ke_center = float(kec[i_min])
    t_center = float(t[first][i_min])
    # user's "lowest energy point" estimate: launch -> lowest-KE point
    S_center = (ke_launch - ke_center) / (z_center - z_launch)

    # did the ion cleanly TRAVERSE the slab (cross the far face +12.5)?
    traversed = bool(np.any(z >= SLAB_FACE))

    # equal-potential-face method: KE at z=-12.5 (entry) vs z=+12.5 (exit).
    # Only meaningful for a clean traversal.
    def ke_at_z(z_target):
        idx = int(np.argmax(z >= z_target))
        return float(ke_ev[idx]), float(t[idx])
    ke_entry, t_entry = ke_at_z(-SLAB_FACE)
    if traversed:
        ke_exit, t_exit = ke_at_z(+SLAB_FACE)
        S_face = (ke_entry - ke_exit) / (2 * SLAB_FACE)
    else:
        ke_exit, t_exit, S_face = float("nan"), float("nan"), float("nan")

    # Ehrenfest energy balance (robust, geometry-free): the electronic system gains
    # exactly the ion KE it loses. dKE_ion over the whole run = energy into electrons.
    ke_final = float(ke_ev[-1])
    dKE_ion_ev = ke_launch - ke_final
    # classical bath gain via E_total directly (ion not absorbed; contaminated by the
    # ion-bath interaction if the ion is still near the slab at t_f).
    clo = pd.read_csv(CL_OBS / "observables.csv")
    E0_cl = float(clo["energy_total"].iloc[0])
    Ef_cl = float(clo["energy_total"].iloc[-1])
    dE_total_cl_ev = (Ef_cl - E0_cl) * HA_TO_EV

    return dict(
        z_launch=z_launch, ke_launch_ev=ke_launch, ke_final_ev=ke_final,
        z_max=float(z.max()), z_final=float(z[-1]), traversed=traversed,
        z_center=z_center, ke_center_ev=ke_center, t_center_au=t_center,
        S_center_ev_per_bohr=S_center,
        ke_entry_ev=ke_entry, ke_exit_ev=ke_exit, t_entry_au=t_entry, t_exit_au=t_exit,
        S_face_ev_per_bohr=S_face,
        dKE_ion_ev=dKE_ion_ev, S_dKE_ev_per_bohr=dKE_ion_ev / L_Z_BOHR,
        E_total_0_cl_ha=E0_cl, E_total_f_cl_ha=Ef_cl,
        dE_total_cl_ev=dE_total_cl_ev, S_etotal_cl_ev_per_bohr=dE_total_cl_ev / L_Z_BOHR,
        time_au=t, z_bohr=z, ke_ev=ke_ev,
    )


# velocity / energy of the 100 eV projectile point (m_e=1: E = 1/2 v^2)
E_POINT_EV = 100.0
V_POINT_AU = float(np.sqrt(2.0 * E_POINT_EV / HA_TO_EV))   # = 2.711 a.u.


def _fmt(x, sf=2):
    """Round to sf significant figures for display."""
    if x == 0:
        return "0"
    from math import floor, log10
    d = sf - int(floor(log10(abs(x)))) - 1
    return f"{round(x, d):g}"


def _wp_ekin_series(run="p2"):
    """<T_WP>(t) per-norm kinetic series (e_kin_ha) for the residual diagnostic."""
    df = pd.read_csv(cfg(run)["wp_obs"] / "wp_momentum_stats.csv", comment="#")
    return df["time_au"].to_numpy(), df["e_kin_ha"].to_numpy()


def make_figures(figdir, run="p2"):
    """Two diagnostic figures: (1) WP convergence/absorption; (2) classical KE(z)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style.apply_theme()
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    wp = compute_wp_ledger(run)
    cl = compute_classical_slab(run)
    tk, ek = _wp_ekin_series(run)

    # ---- FIG 1: WP convergence (why S is an upper bound) ----
    fig, axs = plt.subplots(3, 1, figsize=(6.2, 7.0), sharex=True)
    t, E = wp["time_au"], wp["E_total_ha"] * HA_TO_EV
    axs[0].plot(t, E, "-", color="C0", lw=1.6)
    axs[0].axhline(wp["E_jellium_0_ha"] * HA_TO_EV, ls="--", color="k", lw=1.0,
                   label=r"$E_{\rm jellium}(0)$")
    axs[0].axhline(wp["E_GS_ha"] * HA_TO_EV, ls=":", color="0.5", lw=1.0, label=r"$E_{\rm GS}$")
    axs[0].set_ylabel(r"$E_{\rm total}$ (eV)")
    axs[0].legend(fontsize=7, frameon=False)
    axs[0].set_title(rf"WP run is NOT converged at $\tau={wp['t_f_au']:.0f}$ a.u. "
                     r"$\Rightarrow$ $S$ is an upper bound")
    axs[1].plot(wp["N_time_au"], wp["N_total"], "-", color="C1", lw=1.6)
    axs[1].set_ylabel(r"$N_{\rm total}$")
    axs[1].annotate(f"{wp['n_absorbed']:.2f} e absorbed\n(WP norm rem. "
                    f"{wp['wp_norm_remaining']:.3f}; gate <0.02)",
                    xy=(0.5, 0.4), xycoords="axes fraction", fontsize=7.5)
    axs[2].plot(tk, ek * HA_TO_EV, "-", color="C3", lw=1.6)
    axs[2].axhline(wp["dE_ev"], ls="--", color="k", lw=1.0,
                   label=rf"$\Delta E$ = {wp['dE_ev']:.0f} eV")
    axs[2].set_ylabel(r"$\langle T_{\rm WP}\rangle$ (eV)")
    axs[2].set_xlabel("time (a.u.)")
    axs[2].legend(fontsize=7, frameon=False)
    axs[2].annotate(r"residual WP KE($t_f$) $\approx \Delta E$:" "\n"
                    "the 'deposit' is mostly the\nun-absorbed packet, not the bath",
                    xy=(0.30, 0.55), xycoords="axes fraction", fontsize=7.5)
    fig.tight_layout()
    p1 = figdir / "wp_convergence.png"
    fig.savefig(p1, dpi=140, bbox_inches="tight")
    plt.close(fig)

    # ---- FIG 2: classical KE(z) ----
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    z, ke = cl["z_bohr"], cl["ke_ev"]
    ax.plot(z, ke, "-", color="0.6", lw=1.0, alpha=0.6, label="full trajectory")
    for zf in (-SLAB_FACE, SLAB_FACE):
        ax.axvline(zf, ls="--", color="C2", lw=1.0)
    if cl["traversed"]:
        m = (z >= cl["z_launch"]) & (np.arange(len(z)) <= np.argmax(z > SLAB_FACE))
        ax.plot(z[m], ke[m], "-", color="C0", lw=1.8, label="first transit")
        ax.plot([-SLAB_FACE, SLAB_FACE], [cl["ke_entry_ev"], cl["ke_exit_ev"]], "o",
                color="k", ms=6, label=f"equal-pot. faces: S={cl['S_face_ev_per_bohr']:.2f}")
        ax.plot(cl["z_center"], cl["ke_center_ev"], "v", color="C3", ms=9,
                label=f"slab-centre min: S={cl['S_center_ev_per_bohr']:.1f}")
        ax.set_title("Classical projectile: conservative well vs true (face) stopping")
    else:
        ax.plot(cl["z_max"], 0, "X", color="C3", ms=11,
                label=f"STALLED at z={cl['z_max']:.1f} (KE$\\to$0), reversed")
        ax.set_title(f"Classical projectile ANOMALY: did NOT traverse "
                     f"(ΔKE_ion={cl['dKE_ion_ev']:.0f} eV — trapped/reflected)")
    ax.set_xlabel("ion position  z  (Bohr)")
    ax.set_ylabel(r"$KE_{\rm ion}$ (eV)")
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    p2 = figdir / "classical_ke_z.png"
    fig.savefig(p2, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return {"wp_convergence": str(p1), "classical_ke_z": str(p2)}


if __name__ == "__main__":
    import sys as _sys
    run = _sys.argv[1] if len(_sys.argv) > 1 else "p2"
    wp = compute_wp_ledger(run)
    cl = compute_classical_slab(run)
    print(f"\n=== QUANTUM (WP) ENERGY LEDGER — {wp['label']} ===")
    print(f"  E_total(0)      = {_fmt(wp['E_total_0_ha'],4)} Ha  ({_fmt(wp['E_total_0_ha']*HA_TO_EV,3)} eV)")
    print(f"  <T_WP>(0)       = {_fmt(wp['T_wp_ha'],3)} Ha  ({_fmt(wp['T_wp_ha']*HA_TO_EV,3)} eV "
          f"= drift {_fmt(wp['T_drift_ev'],3)} + zero-point {_fmt(wp['T_zp_ev'],2)})")
    print(f"  E_SIE           = {_fmt(wp['E_sie_ev'],2)} eV")
    print(f"  E_jellium(0)    = {_fmt(wp['E_jellium_0_ha'],5)} Ha")
    print(f"  E_GS (bare)     = {_fmt(wp['E_GS_ha'],5)} Ha")
    print(f"  E_jellium(0)-E_GS = {_fmt(wp['E_jellium0_minus_GS_ev'],2)} eV   (consistency check)")
    print(f"  E_total(t_f)    = {_fmt(wp['E_total_f_ha'],5)} Ha")
    print(f"  dE = E_total(t_f)-E_GS = {_fmt(wp['dE_ev'],2)} eV   [FULL LEDGER]")
    print(f"  S_WP (full ledger) = {_fmt(wp['S_wp_ev_per_bohr'],2)} eV/Bohr  (UPPER BOUND)")
    print(f"  [alt] drift-credit dE={_fmt(wp['dE_driftcredit_ev'],2)} eV -> "
          f"S={_fmt(wp['S_driftcredit_ev_per_bohr'],2)} eV/Bohr (zero-point-persists; not headline)")
    print(f"  gate: WP norm remaining={_fmt(wp['wp_norm_remaining'],2)} (need <0.02), "
          f"late slope={_fmt(wp['slope_ev_per_au'],2)} eV/a.u. -> "
          f"{'MET' if wp['gate_met'] else 'NOT MET (upper bound)'}")
    print(f"\n=== CLASSICAL SLAB ({run}) ===")
    print(f"  KE launch={_fmt(cl['ke_launch_ev'],3)} eV @ z={_fmt(cl['z_launch'],3)}; "
          f"final KE={_fmt(cl['ke_final_ev'],2)} eV @ z={_fmt(cl['z_final'],3)}")
    print(f"  traversed slab (z>+12.5)? {cl['traversed']}  (z_max={_fmt(cl['z_max'],3)})")
    if cl["traversed"]:
        print(f"  faces: KE(-12.5)={_fmt(cl['ke_entry_ev'],3)}  KE(+12.5)={_fmt(cl['ke_exit_ev'],3)} "
              f"-> S_face={_fmt(cl['S_face_ev_per_bohr'],2)} eV/Bohr")
    else:
        print("  faces: N/A — ion did NOT traverse (trapped/reflected) -> face method invalid")
    print(f"  dKE_ion (launch-final)={_fmt(cl['dKE_ion_ev'],2)} eV -> S_dKE={_fmt(cl['S_dKE_ev_per_bohr'],2)}")
    print(f"  E_total method: dE={_fmt(cl['dE_total_cl_ev'],2)} eV -> S={_fmt(cl['S_etotal_cl_ev_per_bohr'],2)}")
