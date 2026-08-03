"""Refined per-half diagnostics for the annular-tube channeling twin.

Plan: docs/plans/cylindrical-channeling-ks-stopping.md, section 8.

WHY THIS EXISTS ALONGSIDE ``channeling_stopping.py``
----------------------------------------------------
``channeling_stopping.py`` answers "was the aim met?" with a fit window the
ANALYSIS derives (first breach of f_bore >= 0.95). This module answers a
different, earlier question: *what do the raw diagnostics look like, so a human
can choose the window themselves, per half*. It therefore does NO window
selection and NO verdict. It loads, converts units, and forms the handful of
derived quantities the diagnostics are built from — each one an exact algebraic
identity that ``tests/test_refined.py`` pins.

THE LABEL SWAP — READ THIS BEFORE USING T1/T2
---------------------------------------------
This module uses the USER'S convention (2026-08-02), which is the REVERSE of
``ks_stopping.py``'s:

    THIS MODULE                              ks_stopping.py
    T1 = |<p>|^2 / 2m        drift only   =  T2
    T2 = |<p>|^2/2m + var(p)/2m           =  T1   ( = <p^2>/2m = e_kin_ha )

They are the same two physical quantities under swapped names. The swap is not
cosmetic: the study's conclusion is that the DRIFT channel is the trustworthy
stopping estimator, so reading a T1 curve under the wrong convention inverts the
result. Every frame this module returns names them ``T1_drift_ev`` and
``T2_total_ev`` — never a bare ``T1`` — so a mis-read has to be deliberate.

    T2 - T1 = var(p)/2m = (1/2) sum_d sigma_pd^2

is exact by construction of ``sigma_pd^2 = <p_d^2> - <p_d>^2``, not an
approximation, and is asserted in the tests.

THE TWO POSITION DEFINITIONS
----------------------------
For the wavepacket the "position" is not unique and the two candidates fail
differently, which is why both are carried:

  s_centroid   circular (Resta-phase) centroid of |psi|^2, unwrapped. Periodic-
               exact — it is meaningful even while the packet straddles a cell
               face, which this one does at t=0 (launched 2 Bohr from the -z
               face). It is a property of the DENSITY, so a packet that splits
               into a transmitted and a reflected lobe reports their weighted
               mean, which is a physically real but no longer particle-like
               position.
  s_pintegral  z0 + integral <p_z> dt. A property of the MOMENTUM. Under exact
               Ehrenfest dynamics with m = 1 these must agree; where they part
               company is exactly where "the wavepacket has a trajectory" stops
               being true, so their difference is a diagnostic, not an error.

For the classical half the position is unambiguous and both reduce to it.

UNITS
-----
Everything returned with an ``_ev`` suffix is in eV; positions are Bohr; times
are atomic units. The electron mass is 1 a.u., so p = v and KE = p^2/2 exactly —
no mass factor appears anywhere below, and that is correct, not an omission.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import channeling_stopping as CS  # noqa: E402
import ks_stopping as K  # noqa: E402  (re-exported by channeling_stopping's path setup)

HA_TO_EV = CS.HA_TO_EV

# The free-evolution reference values, from the launched Gaussian. A packet that
# never interacts holds these forever; departure from them IS the interaction.
VAR_P_FREE = CS.VAR_P_FREE          # 0.03125 per Cartesian direction, a.u.
VAR_P_FREE_3D = 3.0 * VAR_P_FREE    # 0.09375 summed over x, y, z
T2_MINUS_T1_FREE_EV = 0.5 * VAR_P_FREE_3D * HA_TO_EV   # 1.2755 eV


# ---------------------------------------------------------------------------
# Wavepacket half
# ---------------------------------------------------------------------------

def wp_frame(name: str = "wp") -> pd.DataFrame:
    """Per-step wavepacket diagnostics, one tidy frame, user T1/T2 convention.

    Columns
    -------
    step, t                     step index, time (a.u.)
    pz, px, py                  <p_d>(t), a.u.  (m = 1, so these are velocities)
    var_pz, var_perp, var_p3d   sigma_pd^2 summed as indicated, a.u.
    T1_drift_ev                 |<p>|^2 / 2m                              [eV]
    T2_total_ev                 |<p>|^2/2m + var(p)/2m = <p^2>/2m         [eV]
    var_term_ev                 T2 - T1 = var(p)/2m                       [eV]
    d_T1_ev, d_T2_ev            the same, minus their t=0 values          [eV]
    s_centroid, s_pintegral     the two position definitions              [Bohr]
    ehrenfest_resid             s_centroid - s_pintegral                  [Bohr]
    f_bore, f_wall, r_mean, sigma_r    the channeling premise, measured
    sigma_z_circ                longitudinal packet width (circular)      [Bohr]
    e_total_ev, d_e_total_ev    INQ's conserved total energy              [eV]
    norm                        orbital norm, the propagation gate
    """
    run = CS.load_wp(name)
    b = run.base
    obs = CS._obs_dir(run.run_dir)

    mom = K._concat_segments(obs, "wp_momentum_stats")
    pos = K._concat_segments(obs, "wp_real_space_stats")
    df = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))

    px = df["px_mean"].to_numpy()
    py = df["py_mean"].to_numpy()
    pz = df["pz_mean"].to_numpy()
    var_px = df["sigma_px2"].to_numpy()
    var_py = df["sigma_py2"].to_numpy()
    var_pz = df["sigma_pz2"].to_numpy()

    # |<p>|^2 uses all three components so that T2 - T1 = var/2m is EXACT.
    # <p_x> and <p_y> are ~1e-12 here (the tube is axisymmetric and the packet is
    # launched on axis), so numerically T1 is the z drift — but writing it as the
    # full vector is what makes the identity hold by construction rather than to
    # within the transverse residual.
    drift2 = px**2 + py**2 + pz**2
    var_p3d = var_px + var_py + var_pz

    T1 = 0.5 * drift2                          # Ha
    T2 = 0.5 * (drift2 + var_p3d)              # Ha; == e_kin_ha identically

    out = pd.DataFrame({
        "step": df["step"].to_numpy(),
        "t": df["time_au"].to_numpy(),
        "px": px, "py": py, "pz": pz,
        "var_px": var_px, "var_py": var_py, "var_pz": var_pz,
        "var_perp": var_px + var_py,
        "var_p3d": var_p3d,
        "T1_drift_ev": T1 * HA_TO_EV,
        "T2_total_ev": T2 * HA_TO_EV,
        "var_term_ev": (T2 - T1) * HA_TO_EV,
        "e_kin_inq_ev": df["e_kin_ha"].to_numpy() * HA_TO_EV,   # cross-check col
        "s_centroid": b.s3,
        "s_pintegral": b.s4,
        "sigma_z_circ": df["sigma_z_circ"].to_numpy(),
        "norm": b.norm,
    })
    out["d_T1_ev"] = out["T1_drift_ev"] - out["T1_drift_ev"].iloc[0]
    out["d_T2_ev"] = out["T2_total_ev"] - out["T2_total_ev"].iloc[0]
    out["ehrenfest_resid"] = out["s_centroid"] - out["s_pintegral"]

    # The channeling premise, on the same time axis. radial_occupancy is written
    # at the same cadence, but merge on step rather than assume it.
    rad = K._concat_segments(obs, "wp_radial_occupancy")
    out = out.merge(
        rad[["step", "f_bore", "f_wall", "r_mean", "sigma_r"]], on="step", how="left")

    ener = K._concat_segments(obs, "observables")[["step", "energy_total"]]
    out = out.merge(ener, on="step", how="left")
    out["e_total_ev"] = out["energy_total"] * HA_TO_EV
    out["d_e_total_ev"] = out["e_total_ev"] - out["e_total_ev"].iloc[0]
    return out.drop(columns=["energy_total"])


def momentum_slices(name: str = "wp") -> pd.DataFrame:
    """Long-format ``momentum_distribution.csv``: step, t, k, n_wp, n_total.

    ``n_wp`` is |psi_WP(k)|^2 for the projectile orbital alone; ``n_total`` is the
    same for the whole electron system. Written on a coarser cadence than the
    scalars (every 15 steps here), so ``nearest_slices`` snaps a requested time to
    an available one rather than interpolating a distribution.

    DO NOT LOAD THIS WITH ``K._concat_segments``. That helper ends with
    ``drop_duplicates(subset="step", keep="last")``, which is right for a
    scalar-per-step observable (a resume boundary legitimately writes one step
    twice) and CATASTROPHIC here: this file carries 128 rows per step, one per
    k-bin, so de-duplicating on ``step`` alone keeps ONE bin per time and
    silently discards the other 127.

    The symptom is not an exception. It is a plausible-looking table in which
    every distribution has collapsed to the single highest-k bin and integrates
    to zero — observed 2026-08-02 while building refined_analysis.ipynb. The
    resume-segment de-duplication is still done, but on ``(step, k)``.
    """
    obs = CS._obs_dir(CS.WP_RESULTS / name)
    files = sorted(Path(obs).glob("momentum_distribution*.csv"))
    if not files:
        raise FileNotFoundError(f"no momentum_distribution*.csv under {obs}")
    df = pd.concat([pd.read_csv(f, comment="#") for f in files], ignore_index=True)
    df = df.rename(columns={"time_au": "t", "k_bohr_inv": "k"})
    df = (df.sort_values(["step", "k"])
            .drop_duplicates(subset=["step", "k"], keep="last")
            .reset_index(drop=True))

    # Guard the invariant rather than trust it: every time slice must carry the
    # same number of k-bins. A ragged frame here means the de-duplication key is
    # wrong again, and it must fail loudly rather than plot a collapsed curve.
    per_step = df.groupby("step").size()
    if per_step.nunique() != 1:
        raise ValueError(
            f"ragged momentum_distribution under {obs}: bins per step range "
            f"{per_step.min()}-{per_step.max()}. Expected one k-grid per step.")
    if per_step.iloc[0] < 2:
        raise ValueError(
            f"momentum_distribution under {obs} has {per_step.iloc[0]} bin(s) per "
            f"step — the k axis has been collapsed (see this docstring).")
    return df


def nearest_slices(md: pd.DataFrame, times) -> list[tuple[float, pd.DataFrame]]:
    """Pick the available momentum-distribution slices nearest to ``times``.

    Returns [(actual_time, slice_frame_sorted_by_k)]. Snapping rather than
    interpolating is deliberate: a linear blend of two distributions at different
    times is not a distribution the run ever had, and a broadening feature would
    be smeared by the blend itself.
    """
    avail = np.sort(md["t"].unique())
    out = []
    for want in times:
        t_act = float(avail[int(np.argmin(np.abs(avail - float(want))))])
        sl = md[md["t"] == t_act].sort_values("k")
        out.append((t_act, sl.reset_index(drop=True)))
    return out


# ---------------------------------------------------------------------------
# Classical half
# ---------------------------------------------------------------------------

def cl_frame(name: str = "classical") -> pd.DataFrame:
    """Per-step classical diagnostics plus the energy-budget closure.

    THE CLOSURE IS THE POINT. The projectile is an external moving-charge
    perturbation, not an INQ ion, so ``energy_ion`` and ``energy_ion_kinetic``
    are identically zero and ``energy_total`` is the BATH energy alone. The
    projectile's kinetic energy is tracked separately by the Ehrenfest
    integrator. Energy is conserved overall, so

        d(E_total_bath) + d(KE_projectile) = 0

    to integrator accuracy, and ``closure_ev`` is that residual. It is the one
    number that certifies the classical stopping power is a real energy transfer
    and not a bookkeeping artefact of the moving perturbation. Measured
    2026-08-02 on the production run: max |closure| = 2.2e-5 eV over 1501 steps.

    Columns
    -------
    step, t
    z, z_unwrapped, vz          projectile trajectory (Bohr, a.u.)
    x, y                        transverse drift — the channeling check
    force_z                     Hellmann-Feynman drag force (Ha/Bohr)
    ke_ev, d_ke_ev              1/2 m v^2 of the projectile                 [eV]
    e_total_ev, d_e_total_ev    INQ total (bath) energy                     [eV]
    closure_ev                  d_e_total_ev + d_ke_ev; must be ~0          [eV]
    """
    run_dir = CS.CL_RESULTS / name
    obs = CS._obs_dir(run_dir)

    proj = K._concat_segments(obs, "projectile")
    ener = K._concat_segments(obs, "observables")[["step", "energy_total"]]
    df = proj.merge(ener, on="step", how="left")

    z_raw = df["proj_z"].to_numpy()

    # DO NOT use the run's own ``proj_z_unwrapped`` column. Measured 2026-08-02 on
    # the production run: proj_z_unwrapped[i] == proj_z[i+1] EXACTLY for all 1501
    # rows, and proj_z_unwrapped[0] = -27.9617 while the launch point is -28.0. It
    # is the POST-update position written against the PRE-update step index, so
    # pairing it with the same row's velocity or kinetic energy mixes two steps.
    # The tell is velocity-Verlet consistency against the recorded proj_vz:
    # proj_z closes to 1.2e-9, proj_z_unwrapped only to 2.1e-6.
    #
    # Consequence for results already published: the offset is a CONSTANT 0.038
    # Bohr (= v dt), so it moves a linear fit's intercept and not its slope — no
    # reported S is affected. It is corrected here because the two columns stop
    # being interchangeable the moment a run actually wraps (this one has
    # n_wraps == 0 throughout, which is why the bug stayed invisible).
    z_unwrapped = K.unwrap_periodic(z_raw, CS.LZ)
    z_unwrapped = z_unwrapped - z_unwrapped[0] + z_raw[0]

    out = pd.DataFrame({
        "step": df["step"].to_numpy(),
        "t": df["time_au"].to_numpy(),
        "z": z_raw,
        "z_unwrapped": z_unwrapped,
        "vz": df["proj_vz"].to_numpy(),
        "x": df["proj_x"].to_numpy() if "proj_x" in df else np.zeros(len(df)),
        "y": df["proj_y"].to_numpy() if "proj_y" in df else np.zeros(len(df)),
        "force_z": df["force_z"].to_numpy() if "force_z" in df else np.nan,
        # Carried only so the one-step lag documented above stays visible to
        # anyone comparing against an older analysis; never used for a fit.
        "z_unwrapped_raw": (df["proj_z_unwrapped"].to_numpy()
                            if "proj_z_unwrapped" in df.columns else np.nan),
        "ke_ev": df["energy_proj_ke"].to_numpy() * HA_TO_EV,
        "e_total_ev": df["energy_total"].to_numpy() * HA_TO_EV,
    })
    out["d_ke_ev"] = out["ke_ev"] - out["ke_ev"].iloc[0]
    out["d_e_total_ev"] = out["e_total_ev"] - out["e_total_ev"].iloc[0]
    out["closure_ev"] = out["d_e_total_ev"] + out["d_ke_ev"]
    return out


# ---------------------------------------------------------------------------
# Interaction energies — both halves, same schema
# ---------------------------------------------------------------------------

INTERACTION_TERMS = ("e_ss", "e_ps", "e_pp", "e_sb", "e_pb")

TERM_LABEL = {
    "e_ss": r"$\Delta E_{SS}$  bath$-$bath",
    "e_ps": r"$\Delta E_{PS}$  projectile$-$bath",
    "e_pp": r"$\Delta E_{PP}$  projectile self-Hartree",
    "e_sb": r"$\Delta E_{SB}$  bath$-$background",
    "e_pb": r"$\Delta E_{PB}$  projectile$-$background",
}


def interactions(half: str, name: str | None = None) -> pd.DataFrame:
    """``interactions.csv`` with ``d_<term>_ev`` delta-from-t0 columns.

    Deltas, not absolutes: E_SB / E_PB / E_BB carry the charged-cell G=0 gauge
    (.claude/rules/decomposed-interaction-energies.md), so their absolute values
    are not comparable between the two representations while their CHANGES are.
    E_SS / E_PS / E_PP are gauge-clean but are shown as deltas too, so every
    curve on the figure starts at zero and can be read against the others.
    """
    if half not in ("wp", "classical"):
        raise ValueError(f"half must be 'wp' or 'classical', got {half!r}")
    root = CS.WP_RESULTS if half == "wp" else CS.CL_RESULTS
    return CS.load_interactions(root / (name or half))


# ---------------------------------------------------------------------------
# Window-scoped stopping fit — used only once the USER has chosen a window
# ---------------------------------------------------------------------------

def fit_in_window(x: np.ndarray, y_ev: np.ndarray, t: np.ndarray,
                  t0: float, t1: float) -> dict:
    """OLS slope of ``-dy/dx`` restricted to t in [t0, t1]. Returns eV/Bohr.

    Thin wrapper so every window the user tries is fitted by ONE code path,
    whichever half and whichever pair of (energy, path) definitions it is
    applied to. ``S`` is the stopping power: energy LOST per unit path, hence
    the sign flip on the slope.
    """
    m = (t >= t0) & (t <= t1)
    n = int(m.sum())
    if n < 3:
        return {"S": float("nan"), "sigma": float("nan"), "r2": float("nan"),
                "n": n, "t0": t0, "t1": t1}
    xs, ys = np.asarray(x)[m], np.asarray(y_ev)[m]
    A = np.vstack([xs, np.ones_like(xs)]).T
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    slope, icpt = coef
    resid = ys - (slope * xs + icpt)
    dof = max(n - 2, 1)
    s2 = float(resid @ resid) / dof
    sxx = float(((xs - xs.mean()) ** 2).sum())
    sigma = float(np.sqrt(s2 / sxx)) if sxx > 0 else float("nan")
    ss_tot = float(((ys - ys.mean()) ** 2).sum())
    r2 = 1.0 - float(resid @ resid) / ss_tot if ss_tot > 0 else float("nan")
    return {"S": float(-slope), "sigma": sigma, "r2": r2,
            "n": n, "t0": float(t0), "t1": float(t1),
            "x0": float(xs[0]), "x1": float(xs[-1])}


# ---------------------------------------------------------------------------
# 2-D momentum maps (k_z, k_perp) — added 2026-08-02
# ---------------------------------------------------------------------------

WF_SUBPATH = "raw/vti/wavefunction_wp"


def momentum_map(name: str = "wp", step: int = 0):
    """(k_z, k_perp, P) for the WP orbital at one saved wavefunction step.

    Reconstructed from the COMPLEX orbital dump, not from
    ``momentum_distribution.csv``. That CSV carries only the radial ``n(|k|)``,
    which folds the drift direction into the same coordinate as the transverse
    spread — so deceleration (weight moving down in k_z) and transverse heating
    (weight moving up in k_perp) are indistinguishable in it. Separating them is
    the entire point of this map.

    Verified 2026-08-02: the t=0 dump round-trips to the recorded moments
    exactly — <k_z> = 1.917011, var(k_z) = 0.031250, var(k_perp) = 0.062500,
    matching ``wp_momentum_stats.csv`` to all printed digits.
    """
    from inqview.visualisation.field_io import load_complex_vti, kz_kperp_map
    p = Path(CS.WP_RESULTS) / name / WF_SUBPATH / f"wavefunction_t{step:06d}.vti"
    if not p.is_file():
        raise FileNotFoundError(f"no wavefunction dump at {p}")
    return kz_kperp_map(load_complex_vti(p))


def available_wf_steps(name: str = "wp") -> list[int]:
    d = Path(CS.WP_RESULTS) / name / WF_SUBPATH
    if not d.is_dir():
        return []
    out = []
    for f in sorted(d.glob("wavefunction_t*.vti")):
        try:
            out.append(int(f.stem.split("_t")[1]))
        except (IndexError, ValueError):
            continue
    return out


def kz_asymmetry(kz: np.ndarray, P: np.ndarray) -> dict:
    """Asymmetry of the longitudinal momentum distribution.

    Answers: *is the impulse the same for every momentum channel?* A rigid
    translation (every channel decelerated equally) preserves the SHAPE — same
    width, same skew. Any change in those means the impulse was
    momentum-dependent.

    ``frac_above`` is computed by INTERPOLATING THE CDF at the mean, not by
    summing the grid points above it. That matters here: the packet is carried
    by only ~8 resolved k_z points, so a hard ``kz > mean`` comparison is
    dominated by where those few points happen to fall — it returns 0.454 for a
    packet that is EXACTLY symmetric (measured at t=0, where the true answer is
    0.500 and the skewness is -0.0000). The interpolated version returns 0.4987.
    """
    kz = np.asarray(kz, dtype=float)
    p = np.asarray(P, dtype=float)
    if p.ndim == 2:
        p = p.sum(axis=1)
    tot = p.sum()
    if not tot > 0:
        raise ValueError("kz_asymmetry: zero weight")
    p = p / tot

    mean = float((kz * p).sum())
    var = float(((kz - mean) ** 2 * p).sum())
    std = float(np.sqrt(var))
    skew = float(((kz - mean) ** 3 * p).sum() / std**3) if std > 0 else float("nan")

    dk = float(kz[1] - kz[0])
    edges = np.concatenate([[kz[0] - dk / 2], kz + dk / 2])
    cdf = np.concatenate([[0.0], np.cumsum(p)])
    frac_below = float(np.interp(mean, edges, cdf))
    median = float(np.interp(0.5, cdf, edges))

    return {"mean_kz": mean, "sigma_kz": std, "skewness": skew,
            "frac_above_mean": 1.0 - frac_below, "frac_below_mean": frac_below,
            "median_minus_mean": median - mean}


# ---------------------------------------------------------------------------
# Twin comparison: impulse, drag, and the combined projectile coupling
# ---------------------------------------------------------------------------

def impulse_comparison(wp: pd.DataFrame, cl: pd.DataFrame) -> pd.DataFrame:
    """Cumulative impulse of each half and their ratio, per step.

    THE IMPULSE IS THE PRIMITIVE QUANTITY, not the energy. Both halves start at
    the same p, and T1 = p^2/2m EXACTLY for the drift channel, so the entire
    T1-vs-classical energy gap is algebraically a gap in delta-p. Reporting the
    ratio directly removes the p^2 nonlinearity that makes the energy curves
    hard to read.

    No smoothing: this is a cumulative quantity, so it is already integrated and
    a boxcar would only add edge bias.
    """
    t = wp["t"].to_numpy()
    dp_wp = wp["pz"].to_numpy() - wp["pz"].to_numpy()[0]
    dp_cl_full = cl["vz"].to_numpy() - cl["vz"].to_numpy()[0]
    dp_cl = np.interp(t, cl["t"].to_numpy(), dp_cl_full)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(np.abs(dp_cl) > 0, dp_wp / dp_cl, np.nan)
    return pd.DataFrame({"t": t, "dp_wp": dp_wp, "dp_cl": dp_cl,
                         "impulse_ratio": ratio})


def combined_projectile_coupling(i_half: pd.DataFrame) -> np.ndarray:
    """delta(E_PS + E_PB) in eV — the projectile's coupling to EVERYTHING else.

    Why the SUM is the meaningful quantity and the two terms separately are not:

        E_PS + E_PB = integral n_P (phi_S - phi_plus) = integral n_P phi_(S+B)

    i.e. the projectile's interaction with the NET charge density of the system
    (bath electrons plus neutralising background). The system is neutral, so
    phi_(S+B) is the screened, short-ranged potential — a well-defined physical
    field. Split apart, phi_S and phi_plus are each the potential of a charged
    subsystem and carry the charged-cell G=0 gauge individually.

    This is why the two halves look wildly different term-by-term and agree in
    the sum (measured 2026-08-02: individual terms differ by up to 2.7 eV, the
    sum by at most 0.24 eV).
    """
    return (i_half["d_e_ps_ev"] + i_half["d_e_pb_ev"]).to_numpy()
