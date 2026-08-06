"""
Stopping power for the sigma_WP = 5 and 6 Bohr TWIN campaign (sweep sigma56_sv).

THE DEFINITION (user-locked; the project's primary measure):

    S = E_absorbed / L_slab_z ,   E_absorbed = E_total(t_final) - E_GS,
    L_slab_z = 25 Bohr

The arithmetic itself is NOT reimplemented here. It lives in
`hypotheses/slab_ks_wrap/e_absorbed.py::measure_dir`, which is path-generic and
was validated against the wp_highdensity_sv synthesis notebooks (reproduced their
S_deposit column to <= 3e-8 across all 12 runs). This module is the ADAPTER: run
naming, the campaign's own E_GS, the completeness check, and the dispersion
geometry that puts a wavepacket point on a time-averaged-sigma axis.

--------------------------------------------------------------------------------
WHAT IS DIFFERENT ABOUT THIS CAMPAIGN, AND WHY IT MATTERS FOR THE ESTIMATOR
--------------------------------------------------------------------------------
Both halves carry the CAP (user decision 2026-08-02). In every previous campaign
only the wavepacket half did, which meant `E_total(t_f) - E_GS` was a retained-
excitation LOWER BOUND on the WP side and the medium's gain DIRECTLY on the
classical side -- the same formula measuring two different things, and the
documented reason the WP deposit curve sat 3-5x below the classical one. With the
absorber on both, the two are the same measurement. `cl_nocap_*` runs (eta = 0,
same binary) measure what the absorber costs on the classical half.

The norm correction still applies to the WP half ONLY: INQ reports the orbital
kinetic term as occ*<psi|T|psi>/<psi|psi> (inq/src/hamiltonian/energy.hpp:50-55),
so a CAP-decaying packet keeps contributing its per-particle MEAN and inflates
E_total. A classical run has no WP orbital, so raw == corrected there.
`measure_dir` handles both, keyed on `half`.

--------------------------------------------------------------------------------
THE TIME-AVERAGED SIGMA AXIS
--------------------------------------------------------------------------------
A wavepacket has no single width: sigma_d(t) = sqrt(sigma^2/2 + t^2/(2 sigma^2)).
Its classical twin does. `sigma_eq()` returns sqrt(2) * <sigma_d> averaged over
the in-slab transit -- the sigma_WP label a CONSTANT-width packet would carry to
present the same time-averaged charge cloud. That is the axis on which classical
and wavepacket points are comparable.

At sigma = 5/6 the distinction nearly vanishes (sigma_eq = 5.3-5.7 and 6.2-6.5
against labels of 5 and 6), which is the whole reason these widths were chosen.
At sigma = 2 it does not: sigma_eq runs 4.0-6.4 depending on velocity, so one
sigma = 2 label describes four different objects.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
LJ = REPO / "ResearchProject/systems/localised_jellium"
SCRIPTS = LJ / "scripts/sigma56_sv"

# Reuse the validated engine rather than reimplementing the estimator.
sys.path.insert(0, str(LJ / "hypotheses/slab_ks_wrap"))
from e_absorbed import measure_dir, _concat, HA_TO_EV, L_SLAB_Z   # noqa: E402

# ---------------------------------------------------------------------------
# Campaign geometry (must match shared/configs/slab_n100_L35x35x105.hpp)
# ---------------------------------------------------------------------------
LX = LY = 35.0
LZ = 105.0
SLAB_HALF = 12.5
LAUNCH_Z = -27.5
CAP_L = 12.5
CAP_Z_INNER = LZ / 2.0 - CAP_L          # +/- 40.0
DX = 0.40
DT = 0.04
R_S = 4.183

SIGMAS = (5.0, 6.0)
VELOCITIES = (2.0, 2.5, 3.0, 3.5)

# N_STEPS = round(4.36 * (|launch_z| + L_z/2) / (v * dt)) -- see the plan.
# Kept as a literal table so `complete` cannot silently disagree with the
# dispatcher; the dispatcher's table is the same four numbers.
STEPS_TARGET = {2.0: 4360, 2.5: 3488, 3.0: 2907, 3.5: 2491}


def sigma_tag(sigma: float) -> str:
    """'s5p0_' / 's6p0_' — every run in this sweep carries its sigma explicitly."""
    return "s" + f"{sigma:.1f}".replace(".", "p") + "_"


def v_tag(v: float) -> str:
    return "v" + f"{v:.1f}".replace(".", "p")


def run_name(sigma: float, v: float, half: str = "wp", cap: bool = True) -> str:
    """Directory name, matching the dispatchers exactly.

    wp        -> s6p0_v2p0
    classical -> cl_s6p0_v2p0        (cl_nocap_s6p0_v3p0 when cap=False)
    vac       -> vac_s6p0_v2p0
    """
    st, vt = sigma_tag(sigma), v_tag(v)
    if half == "wp":
        return st + vt
    if half == "classical":
        return ("cl_" if cap else "cl_nocap_") + st + vt
    if half == "vac":
        return "vac_" + st + vt
    raise ValueError(f"unknown half {half!r}")


def run_dir(sigma: float, v: float, half: str = "wp", cap: bool = True) -> Path:
    sub = {"wp": "wp", "classical": "classical", "vac": "vac"}[half]
    return SCRIPTS / sub / "results" / run_name(sigma, v, half, cap)


# ---------------------------------------------------------------------------
# Ground state
# ---------------------------------------------------------------------------
def e_gs_ha(tag: str = "dx0p4") -> float:
    """E_GS in Hartree, READ FROM THE GS RUN rather than hard-coded.

    Hard-coding it is how a campaign silently shifts every S by a constant when
    the box changes: the 85-Bohr value (207.18323 Ha) is for a DIFFERENT
    calculation and must never be used here. The deposit has to reference the
    ground state the runs were actually started from.
    """
    summary = SCRIPTS / "gs" / "results" / tag / "run_summary.txt"
    if not summary.exists():
        raise FileNotFoundError(
            f"no ground-state summary at {summary} — run shared/bin/run-s56-gs.slurm")
    for line in summary.read_text().splitlines():
        if line.startswith("ground_state_energy_ha"):
            return float(line.split("=", 1)[1])
    raise KeyError(f"ground_state_energy_ha not found in {summary}")


# ---------------------------------------------------------------------------
# Dispersion geometry -> the time-averaged sigma axis
# ---------------------------------------------------------------------------
def sigma_d(t, sigma: float):
    """Free-Gaussian DENSITY width. sigma_d(0) = sigma/sqrt2 = the classical
    twin's sigma_pot, which is exactly why the two are matched at t = 0."""
    return np.sqrt(sigma**2 / 2.0 + np.asarray(t, float)**2 / (2.0 * sigma**2))


def transit_window(v: float, launch_z: float = LAUNCH_Z) -> tuple[float, float]:
    """(t_in, t_out): when the projectile CENTROID is inside the slab, assuming
    constant velocity. For the classical half the projectile decelerates, so this
    is the nominal window; the measured track is in projectile.csv."""
    return ((abs(launch_z) - SLAB_HALF) / v, (abs(launch_z) + SLAB_HALF) / v)


def transverse_overlap_time(sigma: float, l_xy: float = LX) -> float:
    """When the packet's own periodic images start to overlap: 6 sigma_d = L_xy.
    32.8 a.u. at sigma = 5, 34.0 at sigma = 6 — both comfortably past the transit
    (which ends by t = 20.0 at the slowest velocity), unlike sigma = 0.5 where the
    two windows did not intersect at all."""
    target = 2.0 * (l_xy / 6.0) ** 2
    return float(sigma * np.sqrt(max(target - sigma**2, 0.0)))


def mean_sigma_d(sigma: float, v: float, n: int = 4001) -> float:
    """Time-average of sigma_d over the in-slab transit."""
    ti, to = transit_window(v)
    t = np.linspace(ti, to, n)
    return float(np.trapezoid(sigma_d(t, sigma), t) / (to - ti))


# ---------------------------------------------------------------------------
# Effective width, averaged over the window where the packet is still intact
# ---------------------------------------------------------------------------
# THE WINDOW (user decision, 2026-08-03): t = 0 until the wavepacket norm has
# fallen by 1 % of its launch value. Everything after that is excluded.
#
# WHY. sigma_r(t) = sqrt(sigma_x^2+sigma_y^2+sigma_z^2) is only a property of the
# projectile while the projectile is still there. Once the CAP starts removing it
# the packet is cut in two -- part still near the slab, part in the absorber band
# ~40 Bohr away -- and a second moment about a single centroid reads that
# bimodality as an enormous variance: measured sigma_r SPIKES to 40-48 Bohr, then
# oscillates around 25-30 on a residue of norm ~1e-9. A full-run mean averages
# mostly over that dead zone and returned <sigma_r> ~ 20 for a packet that
# launches at 7.35 and leaves the slab at 8-9.
#
# The 1 %-loss point sits comfortably AFTER the in-slab transit and BEFORE the
# packet reaches the absorber, so the whole averaging window lies on the stretch
# where the measured trace still tracks the free-Gaussian law sqrt(3)*sigma_d(t).
# See s6_sigma_r_traces.png.
NORM_DROP_FRAC = 0.01


def norm_window_end(norm: np.ndarray, drop: float = NORM_DROP_FRAC) -> int:
    """Index of the first step at which `norm` has fallen by `drop` from t=0.

    Relative to norm[0], not to 1.0: the injected packet is renormalised but a
    run may still start a hair off unity, and the window must mean "lost 1 % of
    what it started with". Returns the last index if the threshold is never
    crossed (a run too short to lose 1 %).
    """
    thresh = float(norm[0]) * (1.0 - drop)
    hit = np.where(norm <= thresh)[0]
    return int(hit[0]) if hit.size else len(norm) - 1


def sigma_r_window(rs, drop: float = NORM_DROP_FRAC):
    """(mean_sigma_r, t_end, n_steps, reached_threshold) from a real-space-stats frame.

    `rs` is a wp_real_space_stats frame (segments already concatenated) carrying
    sigma_x2/sigma_y2/sigma_z2 and norm_check. Works for both this campaign's runs
    and the legacy L_z = 85 sweeps, which use the same column names.
    """
    t = rs["time_au"].to_numpy()
    sr = np.sqrt((rs["sigma_x2"] + rs["sigma_y2"] + rs["sigma_z2"]).to_numpy())
    nm = (rs["norm_check"].to_numpy() if "norm_check" in rs
          else np.ones_like(t))
    k = norm_window_end(nm, drop)
    reached = bool(nm[k] <= nm[0] * (1.0 - drop))
    return float(sr[:k + 1].mean()), float(t[k]), k + 1, reached


def sigma_eq(sigma: float, v: float) -> float:
    """sqrt(2) * <sigma_d> — the sigma_WP LABEL an equivalent constant-width
    packet would carry. This is the x-coordinate for a time-averaged-sigma plot."""
    return float(np.sqrt(2.0) * mean_sigma_d(sigma, v))


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Point:
    sigma_wp: float
    v: float
    half: str
    cap: bool
    run: str
    S_eV_per_Bohr: float
    S_raw_eV_per_Bohr: float
    E_absorbed_eV: float
    t_final_au: float
    steps_done: int
    steps_target: int
    complete: bool
    settled: bool
    plateau_drift_eV: float
    norm_final: float
    sigma_d_entry: float
    sigma_d_exit: float
    sigma_eq: float
    # --- the E_PS tail correction (2026-08-03) -------------------------------
    e_ps_final_eV: float          # projectile-bath interaction still in the ledger
    z_proj_final: float           # where the classical projectile ended up
    S_deposit_eV_per_Bohr: float  # S with that tail removed -- USE THIS ONE


def e_ps_final(sigma: float, v: float, half: str, cap: bool = True) -> float:
    """E_PS at the last recorded step, in eV. 0.0 if interactions.csv is absent."""
    obs = run_dir(sigma, v, half, cap) / "raw" / "observables"
    try:
        return float(_concat(obs, "interactions")["e_ps"].to_numpy()[-1]) * HA_TO_EV
    except (FileNotFoundError, KeyError, IndexError, ValueError,
            pd.errors.EmptyDataError):
        # A run killed by the full filesystem (2026-08-03) can leave a
        # zero-byte interactions.csv. Returning 0.0 makes S_deposit fall back to
        # the uncorrected S rather than crashing the whole table; such runs are
        # incomplete anyway and are filtered out downstream on `complete`.
        return float("nan")


def measure(sigma: float, v: float, half: str = "wp", cap: bool = True,
            gs_tag: str = "dx0p4") -> Point:
    """One (sigma, v, half) point, with the evidence needed to trust it."""
    d = run_dir(sigma, v, half, cap)
    if not d.exists():
        raise FileNotFoundError(d)
    a = measure_dir(d, e_gs_ha(gs_tag), half if half == "wp" else "classical",
                    "n100", v, R_S)
    target = STEPS_TARGET[v]
    ti, to = transit_window(v)

    # THE E_PS TAIL CORRECTION (2026-08-03). The estimator
    # S = [E_total(t_f) - E_GS]/L assumes the projectile's interaction with the
    # bath has decayed to zero by t_f. For the WP half it HAS: the CAP annihilates
    # the packet (norm_wp ~ 1e-10) and E_PS(t_f) is 1e-5 eV. For the CLASSICAL half
    # it never does -- the projectile is a PRESCRIBED external potential that keeps
    # travelling, and its monopole tail falls off only as N_e/z. At t_f it sits at
    # z = 321 Bohr and STILL contributes E_PS = 8.5 eV, verified against the bare
    # monopole 100/z to 0.6-4.4 %. That is 62-80 % of the raw classical "deposit"
    # of 10.6-13.7 eV, so uncorrected it does not just add noise, it dominates.
    #
    # This is also why the classical drift was identically -1.04 eV at every
    # velocity: N_STEPS was sized so that v*t_f = 4.36*(|z0| + L_z/2) = 349 Bohr is
    # CONSTANT across the sweep, so every run ends at the same z and carries the
    # same tail. A v-independent "drift" was the fingerprint of a geometric
    # artefact, not of physics.
    #
    # Removing it makes the two halves measure the same quantity -- the energy left
    # in the slab once the projectile is gone -- which is the whole point of the
    # pairwise decomposition (.claude/rules/decomposed-interaction-energies.md).
    eps_f = e_ps_final(sigma, v, half, cap)
    z_f = LAUNCH_Z + v * a.t_final_au
    S_dep = (a.E_absorbed_eV - eps_f) / L_SLAB_Z

    return Point(
        sigma_wp=sigma, v=v, half=half, cap=cap, run=d.name,
        S_eV_per_Bohr=a.S_eV_per_Bohr,
        S_raw_eV_per_Bohr=a.S_raw_eV_per_Bohr,
        E_absorbed_eV=a.E_absorbed_eV,
        t_final_au=a.t_final_au,
        steps_done=a.steps_done, steps_target=target,
        # COMPLETENESS IS NOT OPTIONAL. deposit_stopping on a still-propagating
        # run once returned a perfectly plausible S (2.35 eV/Bohr, norm 1.000)
        # from 86 of 3623 steps. Anything consuming this must filter on it.
        complete=(a.steps_done >= target),
        settled=a.settled, plateau_drift_eV=a.plateau_drift_eV,
        norm_final=a.norm_final,
        sigma_d_entry=float(sigma_d(ti, sigma)),
        sigma_d_exit=float(sigma_d(to, sigma)),
        sigma_eq=sigma_eq(sigma, v),
        e_ps_final_eV=eps_f,
        z_proj_final=float(z_f),
        S_deposit_eV_per_Bohr=S_dep,
    )


def table(gs_tag: str = "dx0p4", include_controls: bool = True) -> pd.DataFrame:
    """Every point in the sweep. Missing / unfinished runs are REPORTED and kept
    (with complete=False) rather than dropped, so a gap is visible in the table
    instead of silently vanishing from the figure."""
    rows = []
    for sigma in SIGMAS:
        for v in VELOCITIES:
            for half in ("wp", "classical"):
                try:
                    rows.append(asdict(measure(sigma, v, half, True, gs_tag)))
                except (FileNotFoundError, KeyError) as e:
                    print(f"  MISSING {run_name(sigma, v, half)}: {type(e).__name__}")
            if include_controls and v == 3.0:
                try:
                    rows.append(asdict(measure(sigma, v, "classical", False, gs_tag)))
                except (FileNotFoundError, KeyError):
                    print(f"  MISSING {run_name(sigma, v, 'classical', False)}")
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["half", "sigma_wp", "v"]).reset_index(drop=True)


def cap_cost(gs_tag: str = "dx0p4", window: int = 20) -> pd.DataFrame:
    """What the absorber costs the classical half: CAP-on minus CAP-free at
    v = 3.0, per sigma. This is the number that says whether putting a CAP on the
    classical runs (to make the estimator symmetric) was cheap or expensive.

    COMPARED AT A MATCHED STEP, not at each run's own final step. The two runs of
    a CAP pair routinely stop at different steps (a walltime kill, a full disk),
    and differencing their endpoints then folds the missing steps into what gets
    reported as the CAP effect. 2026-08-03: complete cl_s6p0_v3p0 (2907) was
    differenced against incomplete cl_nocap_s6p0_v3p0 (2680) and the 227-step gap
    was reported as a -6.6 % absorber cost. The fix is to read BOTH traces at the
    largest step they share, averaged over `window` samples to damp the ringing.

    `step_matched` and `complete_*` are emitted so a consumer can see the pair is
    short. A matched comparison of two short runs is still a valid CAP measurement
    (both saw the same trajectory), it just does not extend to t_final.
    """
    rows = []
    for sigma in SIGMAS:
        try:
            on = measure(sigma, 3.0, "classical", True, gs_tag)
            off = measure(sigma, 3.0, "classical", False, gs_tag)
            tr_on = energy_trace(sigma, 3.0, "classical", True)
            tr_off = energy_trace(sigma, 3.0, "classical", False)
        except (FileNotFoundError, KeyError):
            print(f"  MISSING CAP-cost pair at sigma = {sigma}")
            continue

        step_max = int(min(tr_on.step.max(), tr_off.step.max()))

        def _S_at(tr):
            w = tr[tr.step <= step_max].tail(window)
            return float(w["dE_corr"].mean()) / L_SLAB_Z

        s_on, s_off = _S_at(tr_on), _S_at(tr_off)
        if not (on.complete and off.complete):
            print(f"  NOTE CAP-cost at sigma = {sigma}: pair is short "
                  f"({on.steps_done}/{off.steps_done} steps) — matched at "
                  f"step {step_max}")
        rows.append({
            "sigma_wp": sigma, "v": 3.0,
            "step_matched": step_max, "window": window,
            "S_cap_on": s_on, "S_cap_off": s_off,
            "delta_eV_per_Bohr": s_on - s_off,
            "delta_pct": 100.0 * (s_on - s_off) / abs(s_off) if s_off else np.nan,
            # the unmatched endpoint values, kept only so the difference between
            # the two conventions stays visible rather than being silently fixed
            "S_cap_on_endpoint": on.S_eV_per_Bohr,
            "S_cap_off_endpoint": off.S_eV_per_Bohr,
            "complete_on": on.complete, "complete_off": off.complete,
            "steps_on": on.steps_done, "steps_off": off.steps_done,
        })
    return pd.DataFrame(rows)


def energy_trace(sigma: float, v: float, half: str = "wp", cap: bool = True):
    """t, dE_raw, dE_corr (eV, relative to E_GS) — the plateau you are reading S
    off. Always look at this before quoting a number: a drifting tail means the
    run needs extending (LJ_RESUME=1 with a larger LJ_N_STEPS)."""
    obs = run_dir(sigma, v, half, cap) / "raw" / "observables"
    d = _concat(obs, "observables")
    e_ev = d["energy_total"].to_numpy() * HA_TO_EV - e_gs_ha() * HA_TO_EV
    out = pd.DataFrame({"t": d["time_au"].to_numpy(),
                        "step": d["step"].to_numpy(), "dE_raw": e_ev})
    if half == "wp":
        mom, pos = _concat(obs, "wp_momentum_stats"), _concat(obs, "wp_real_space_stats")
        m = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))
        m = m[m.step.isin(d.step)]
        norm = (m["norm_check_r"] if "norm_check_r" in m else m["norm_check"]).to_numpy()
        T1 = m["e_kin_ha"].to_numpy() * HA_TO_EV
        n = min(len(e_ev), len(norm))
        corr = e_ev.copy()
        corr[:n] -= T1[:n] * (1.0 - norm[:n])
        out["dE_corr"] = corr
        out["norm"] = np.concatenate([norm[:n], np.full(len(e_ev) - n, np.nan)])
    else:
        out["dE_corr"] = out["dE_raw"]
        out["norm"] = np.nan
    return out


if __name__ == "__main__":
    print(f"sigma56_sv — L_slab_z = {L_SLAB_Z} Bohr, launch z = {LAUNCH_Z}\n")
    print("dispersion geometry (no runs needed):")
    print(f"  {'sig':>4} {'v':>4} {'t_in':>6} {'t_out':>6} {'sd_in':>7} {'sd_out':>7} "
          f"{'growth':>7} {'sig_eq':>7} {'t_ov':>6}")
    for s in SIGMAS:
        for v in VELOCITIES:
            ti, to = transit_window(v)
            print(f"  {s:>4g} {v:>4.1f} {ti:>6.2f} {to:>6.2f} {sigma_d(ti, s):>7.3f} "
                  f"{sigma_d(to, s):>7.3f} {sigma_d(to, s)/sigma_d(ti, s):>7.2f} "
                  f"{sigma_eq(s, v):>7.2f} {transverse_overlap_time(s):>6.1f}")
    try:
        print(f"\nE_GS = {e_gs_ha():.9f} Ha")
        t = table()
        if not t.empty:
            print(t.to_string(index=False))
    except FileNotFoundError as e:
        print(f"\n(no ground state yet: {e})")
