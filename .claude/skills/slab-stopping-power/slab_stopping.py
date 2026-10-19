"""slab-stopping-power — E_absorbed / L_z for a localised-jellium SLAB run.

Self-contained kernels for the `slab-stopping-power` skill. Point them at a run's
observables directory and get the energy absorbed by the slab (eV) and the
stopping power S (eV/Bohr) for either the classical or the wavepacket half, with
the gates that make the number trustworthy.

DEFINITION
    S = E_absorbed / L_z            (L_z = slab thickness, default 25 Bohr)

E_absorbed = energy retained in the jellium slab once the projectile is gone,
measured against the PROJECTILE-FREE ground state E_GS -- NOT against E_total(t=0).
E_total(0) carries the projectile's launch energy (classical: the static E_PS(0);
WP: the wavepacket formation energy), which E_GS does not; only E_GS isolates the
deposit. See SKILL.md.

    classical :  E_absorbed = [ E_total_INQ(t_f) + E_PB(t_f) ] - E_GS
                 INQ's energy_total omits E_PB (the projectile is an external
                 potential, not a density), so it must be added from the ledger.

    wavepacket:  E_absorbed = [ E_total_INQ(t_f) - T1*(1-norm) ]_{t_f} - E_GS
                 E_PB is already inside energy_total (the WP is a real orbital).
                 The norm term removes INQ's per-particle-mean inflation as the CAP
                 drains the packet. Valid ONLY if the WP is fully absorbed (norm~0).

All energies are read in Hartree from the CSVs and returned in eV.
"""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd

HA_TO_EV = 27.211386
L_Z_DEFAULT = 25.0


def concat_obs(obs_dir: str | Path, base: str) -> pd.DataFrame:
    """`base.csv` + every resume segment `base.fromNNNN.csv`, in step order.

    Concatenating is mandatory: reading only `base.csv` silently truncates a
    RESUMED run at its first kill point, and the deposit lives at t_final.
    """
    obs_dir = Path(obs_dir)
    files = sorted(glob.glob(str(obs_dir / f"{base}.csv"))
                   + glob.glob(str(obs_dir / f"{base}.from*.csv")))
    if not files:
        raise FileNotFoundError(f"no {base}.csv under {obs_dir}")
    df = pd.concat([pd.read_csv(f, comment="#") for f in files])
    return df.drop_duplicates("step").sort_values("step").reset_index(drop=True)


def extensive_kinetic_total(obs_dir: str | Path,
                            *, identity_tol_ha: float = 5e-3) -> pd.DataFrame | None:
    """Per-step norm-corrected total energy for ANY slab run (classical OR wp),
    applying the ALL-ORBITAL extensive-kinetic fix. Returns None if the run did
    not emit `orbital_kinetic_stats.csv` (legacy runs -> caller falls back).

    THE ALL-ORBITAL CORRECTION (supersedes the single-WP-orbital fix). INQ reports
    the kinetic energy per-particle-normalized, Σ_i occ_i·T_i/⟨ψ_i|ψ_i⟩
    (energy.hpp:55). Under a CAP EVERY draining orbital -- the wavepacket AND the
    bath orbitals it excites into the absorber -- has its kinetic pinned at the
    surviving remnant's mean and inflates energy_total. `OrbitalKineticStats`
    (inqkit) logs the EXTENSIVE kinetic kin_bare_total = Σ_i occ_i·T_i (no norm
    division, all orbitals). The corrected total swaps INQ's normalized kinetic for
    the extensive one:

        E_total_corr = energy_total_INQ − energy_kinetic_INQ + kin_bare_total

    This is the single-WP-orbital formula E − T1·(1−norm) generalised to the whole
    KS set, and it is applied identically to the classical and the WP half.

    GATE: kin_normdiv_total (OrbitalKineticStats' own reconstruction of INQ's
    reduction) must equal INQ's reported energy_kinetic -- proves the extensive
    kinetic is the same operator INQ used, differenced only by the norm division.
    """
    obs_dir = Path(obs_dir)
    if not (obs_dir / "orbital_kinetic_stats.csv").exists() \
       and not list(obs_dir.glob("orbital_kinetic_stats.from*.csv")):
        return None
    o = concat_obs(obs_dir, "observables")
    k = concat_obs(obs_dir, "orbital_kinetic_stats")
    m = pd.merge(o, k, on=["step", "time_au"])
    resid = float((m["kin_normdiv_total_ha"] - m["energy_kinetic"]).abs().max())
    if resid > identity_tol_ha:
        raise ValueError(
            f"extensive-kinetic identity failed: max|kin_normdiv − energy_kinetic| "
            f"= {resid:.2e} Ha > {identity_tol_ha:.0e}. OrbitalKineticStats is not "
            f"reproducing INQ's kinetic reduction — do not trust the correction.")
    m["e_total_corr_ha"] = (m["energy_total"] - m["energy_kinetic"]
                            + m["kin_bare_total_ha"])
    return m[["step", "time_au", "energy_total", "energy_kinetic",
              "kin_bare_total_ha", "kin_normdiv_total_ha", "norm_total",
              "e_total_corr_ha"]]


def e_absorbed_classical(obs_dir: str | Path, e_gs_ha: float,
                         *, closure_tol_ha: float = 1e-6) -> float:
    """E_absorbed (eV) for a classical projectile, with E_PB restored.

    Requires `observables.csv` and `interactions.csv` (the pairwise ledger).
    GATE: INQ's energy_external must equal E_SB + E_PS (i.e. E_PB genuinely
    absent) -- if it does not, this is not a plain classical direct-potential run
    and the +E_PB reconstruction does not apply.
    """
    o = concat_obs(obs_dir, "observables")
    ix = concat_obs(obs_dir, "interactions")
    m = pd.merge(o, ix, on=["step", "time_au"])
    resid = float((m["energy_external"] - (m["e_sb"] + m["e_ps"])).abs().max())
    if resid > closure_tol_ha:
        raise ValueError(
            f"classical closure failed: max|E_ext-(E_SB+E_PS)| = {resid:.2e} Ha "
            f"> {closure_tol_ha:.0e}. Either E_PB is already in energy_total "
            f"(is this a WP run?) or the ledger convention differs.")
    # ALL-ORBITAL norm correction if the run emitted orbital_kinetic_stats.csv;
    # else the raw INQ total (legacy runs: bath norm loss is negligible anyway).
    corr = extensive_kinetic_total(obs_dir)
    e_tot_f = (float(corr["e_total_corr_ha"].iloc[-1]) if corr is not None
               else float(m["energy_total"].iloc[-1]))
    e_pb_f = float(m["e_pb"].iloc[-1])
    return (e_tot_f + e_pb_f - e_gs_ha) * HA_TO_EV


def e_absorbed_wp(obs_dir: str | Path, e_gs_ha: float,
                  *, norm_max: float = 1e-3) -> tuple[float, float]:
    """(E_absorbed eV, norm_final) for a wavepacket projectile.

    PREFERRED: `observables.csv` + `orbital_kinetic_stats.csv` -> the ALL-ORBITAL
    extensive-kinetic correction (see extensive_kinetic_total). FALLBACK (legacy
    sigma56 / wp_highdensity_sv runs): `wp_momentum_stats.csv` (e_kin_ha) +
    `wp_real_space_stats.csv` -> the single-WP-orbital correction E − T1·(1−norm).
    Either way E_PB needs no handling -- it is inside energy_total already.

    GATE: the WP must be fully absorbed (norm_final <= norm_max). If it is not,
    E_total(t_f) still holds the WP's formation energy and E_absorbed is a
    mid-transit excitation, not the deposit -- EXTEND the run (LJ_RESUME=1).
    """
    # PREFERRED: all-orbital extensive-kinetic correction (orbital_kinetic_stats.csv).
    # Gate WP absorption on the WP orbital norm (norm_wp in interactions.csv), NOT
    # on norm_total (which is ~N — all electrons — and never falls to 0).
    ext = extensive_kinetic_total(obs_dir)
    if ext is not None:
        ix = concat_obs(obs_dir, "interactions")
        norm_f = (float(ix["norm_wp"].iloc[-1]) if "norm_wp" in ix
                  else float(ext["norm_total"].iloc[-1]))
        if norm_f > norm_max:
            raise ValueError(
                f"WP not absorbed: norm_wp_final = {norm_f:.2e} > {norm_max:.0e}. "
                f"E_total(t_f) still carries WP formation energy; extend the run.")
        return (float(ext["e_total_corr_ha"].iloc[-1]) - e_gs_ha) * HA_TO_EV, norm_f

    # LEGACY fallback: single-WP-orbital correction (wp_momentum_stats.csv), used by
    # the sigma56 / wp_highdensity_sv sweeps that predate OrbitalKineticStats.
    o = concat_obs(obs_dir, "observables")
    mom = concat_obs(obs_dir, "wp_momentum_stats")
    rs = concat_obs(obs_dir, "wp_real_space_stats")
    m = pd.merge(mom, rs, on=["step", "time_au"], suffixes=("_p", "_r"))
    m = m[m.step.isin(o.step)]
    norm = (m["norm_check_r"] if "norm_check_r" in m else m["norm_check"]).to_numpy()
    T1 = m["e_kin_ha"].to_numpy()                 # Hartree
    e = o["energy_total"].to_numpy()              # Hartree
    n = min(len(e), len(norm))
    corr = e[:n] - T1[:n] * (1.0 - norm[:n])
    norm_f = float(norm[n - 1])
    if norm_f > norm_max:
        raise ValueError(
            f"WP not absorbed: norm_final = {norm_f:.2e} > {norm_max:.0e}. "
            f"E_total(t_f) still carries WP formation energy; extend the run.")
    return (corr[-1] - e_gs_ha) * HA_TO_EV, norm_f


def stopping_power(e_absorbed_ev: float, L_z: float = L_Z_DEFAULT) -> float:
    """S (eV/Bohr) = E_absorbed / L_z."""
    return e_absorbed_ev / L_z


def measure(obs_dir: str | Path, half: str, e_gs_ha: float,
            L_z: float = L_Z_DEFAULT) -> dict:
    """One run -> {E_absorbed_eV, S_eV_per_Bohr, norm_final?}. `half` in
    {'classical','wp'}."""
    if half == "classical":
        e_abs = e_absorbed_classical(obs_dir, e_gs_ha)
        out = {"half": half, "E_absorbed_eV": e_abs}
    elif half == "wp":
        e_abs, norm_f = e_absorbed_wp(obs_dir, e_gs_ha)
        out = {"half": half, "E_absorbed_eV": e_abs, "norm_final": norm_f}
    else:
        raise ValueError(f"half must be 'classical' or 'wp', got {half!r}")
    out["S_eV_per_Bohr"] = stopping_power(e_abs, L_z)
    out["L_z"] = L_z
    return out
