"""Stopping-power extraction for the annular-tube CHANNELING TWIN.

Plan: docs/plans/cylindrical-channeling-ks-stopping.md

Deterministic engine: it does the arithmetic and nothing else, so every notebook
and every table gets bit-identical numbers from one implementation. The four
KS-orbital stopping definitions are NOT re-implemented here — they live in
``ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/ks_stopping.py``
and are imported. Applying the BULK definitions unchanged to a new geometry is
the point of the study; re-deriving them here would defeat it.

WHAT THIS MODULE ADDS OVER THE BULK ENGINE
------------------------------------------

1. LAYOUT. These runs write to
   ``<scripts>/channeling_twin/{wp,classical}/results/<name>/raw/observables/``.

2. NO IN-MEDIUM PATH CORRECTION, AND THAT IS A RESULT NOT AN OMISSION.
   The slab study needed s5 = integral f v dt because 25 of its 85 Bohr were
   vacuum, so a centroid-path fit averaged drag over medium AND vacuum. The tube
   is UNIFORM along z: the medium fills every z the projectile visits, so the
   in-medium path IS the path and -dT/ds3 / -dT/ds4 are already forces. The
   correction the slab needed does not exist here.

3. THE CHANNELING WINDOW, MEASURED. What the tube CAN silently violate is its
   own premise — that the packet stays in the hollow bore. f_bore(t) from
   ``wp_radial_occupancy.csv`` (inqkit::observables::radial_occupancy) is that
   premise as a curve, and ``channeling_window()`` derives the fit window from
   it rather than from the free-dispersion formula. The formula is reported
   alongside as a cross-check, never used as the authority: the packet stops
   being Gaussian as soon as it scatters.

4. THE MECHANISM DIAGNOSTIC. In bulk, -dT1/ds was not a stopping power because
   var(p) GREW through interaction with the bath (+6.8 eV) while the drift term
   stayed flat. var(p) is conserved under free evolution, so that growth was
   interaction. ``var_p_freeze()`` measures the same quantity here. The study's
   claim is not merely "S_WP = S_classical" but "S_WP = S_classical BECAUSE
   channeling froze var(p)", and only this function can support the second half.

5. THE SAME-WINDOW COMPARISON. Both twins are light projectiles and both
   decelerate (.claude/rules/light-projectile-stopping.md), so the honest
   comparison fits BOTH over the SAME time window. ``compare()`` does that, and
   also reports the classical initial-drag fit (v >= 0.85 v0) that the rule
   prescribes for a standalone classical number, so the two conventions can be
   seen not to disagree.

WHAT IS DELIBERATELY NOT DONE
-----------------------------
No deposit-based estimator (Method A of the stopping-power-extraction skill) for
the WP half. There is no CAP, so ``energy_total`` is a CONSERVED quantity rather
than a fit target — for the wavepacket its fit target is identically zero because
the projectile is inside the system. Energy conservation is the correctness GATE;
the -dT/ds slopes are the measurement.

CARRIED-OVER CONCLUSION from docs/handovers/bulk-jellium-ks-stopping.md, to be
re-tested here and not re-derived: T1 - T2 grew at a fixed rate in bulk that did
NOT scale with bath density, most plausibly self-interaction error (the packet is
an occupied KS orbital whose own charge enters the Hartree potential — the E_PP
column of interactions.csv). If that reading is right, channeling should NOT
remove it, because SIE is a property of the orbital and not of its environment.
Whether E_PP stays flat while T1-T2 does or does not is therefore a discriminating
test, and ``load_interactions()`` exists to run it.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
_KS = REPO / "ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping"
if str(_KS) not in sys.path:
    sys.path.insert(0, str(_KS))

import ks_stopping as K  # noqa: E402

HA_TO_EV = K.HA_TO_EV

# ---------------------------------------------------------------------------
# Geometry and physics. Mirrored from
# ResearchProject/systems/cylindrical_jellium/shared/configs/channeling_tube_rs3.hpp.
# run_summary.txt is the authority if these ever disagree; summary_of() reads it.
# ---------------------------------------------------------------------------
LX = LY = 40.0
LZ = 60.0
R_IN = 10.0
R_OUT = 14.0
EDGE_W = 0.5
DX = 0.5
DT = 0.02
N_ELEC = 160
SIGMA_WP = 4.0
LAUNCH_Z = -28.0
V0 = 1.91701127                 # 50 eV
N_STEPS = 1500

V_ANNULUS = math.pi * (R_OUT**2 - R_IN**2) * LZ      # 18095.5737 Bohr^3
N0 = N_ELEC / V_ANNULUS                              # 8.84194e-3
RS = (3.0 / (4.0 * math.pi * N0)) ** (1.0 / 3.0)     # 3.000000
OMEGA_P = math.sqrt(4.0 * math.pi * N0)              # 0.333333 a.u. = 9.07 eV
T_PLASMON = 2.0 * math.pi / OMEGA_P                  # 18.85 a.u.
V_FERMI = (3.0 * math.pi**2 * N0) ** (1.0 / 3.0)     # 0.639719
LAMBDA_P = 2.0 * math.pi * V0 / OMEGA_P              # 36.14 Bohr

SIGMA_POT = SIGMA_WP / math.sqrt(2.0)                # 2.82843 = WP density std
SIGMA_P = 1.0 / (math.sqrt(2.0) * SIGMA_WP)          # 0.176777 momentum std
VAR_P_FREE = SIGMA_P**2                              # 0.03125 = 1/(2 sigma^2)
LOCALISATION_EV = 3.0 / (4.0 * SIGMA_WP**2) * HA_TO_EV   # 1.2755 eV = (T1-T2)(0)

SCRIPTS = REPO / "ResearchProject/systems/cylindrical_jellium/scripts/channeling_twin"
# Overridable so the notebook builder and the test suite can be smoke-run against
# synthetic runs without touching the real results tree.
WP_RESULTS = Path(os.environ.get("CHAN_WP_RESULTS", SCRIPTS / "wp" / "results"))
CL_RESULTS = Path(os.environ.get("CHAN_CL_RESULTS", SCRIPTS / "classical" / "results"))

# The channeling window is DERIVED from the measured f_bore(t); this is the
# threshold that defines "still channeling". 0.95 = at most 5 % of the packet has
# reached the wall. Chosen so the window ends before the free-dispersion estimate
# below (a stricter, measured criterion), not after it.
F_BORE_MIN = 0.95
# Fraction of the channeling window dropped as launch transient. The WP is
# orthogonalised against the occupied manifold at t=0 and the bath needs ~one
# screening time to respond, so the first slice is not a steady drag.
TRANSIENT_FRAC = 0.15


def sigma_d(t):
    """Free-dispersion density std of the launched Gaussian, Bohr."""
    return np.sqrt(SIGMA_WP**2 / 2.0 + np.asarray(t, dtype=float) ** 2 / (2.0 * SIGMA_WP**2))


def t_when_sigma_d(target: float) -> float:
    """Time at which sigma_d(t) reaches `target` (Bohr). Velocity-independent."""
    inner = target**2 - SIGMA_WP**2 / 2.0
    return float(np.sqrt(2.0 * SIGMA_WP**2 * inner)) if inner > 0 else 0.0


# 2 sigma_d reaches the bore wall: t = 23.32 a.u. (step 1166 of 1500).
T_2SIGMA_AT_WALL = t_when_sigma_d(R_IN / 2.0)
# Transverse periodic images overlap when 6 sigma_d = L_xy: t = 34.15 a.u., i.e.
# AFTER the run ends. Unlike the slab study, this run never drags an array.
T_TRANSVERSE_OVERLAP = t_when_sigma_d(LX / 6.0)


def summary_of(run_dir: Path) -> dict:
    """run_summary.txt as a flat dict of strings (the run's own provenance)."""
    out: dict[str, str] = {}
    p = Path(run_dir) / "run_summary.txt"
    if not p.is_file():
        return out
    for line in p.read_text().splitlines():
        if "=" not in line:
            continue
        k, _, v = line.partition("=")
        out[k.strip()] = v.strip()
    return out


def is_complete(run_dir: Path) -> bool:
    return summary_of(run_dir).get("run_completed", "false").lower() == "true"


# ---------------------------------------------------------------------------
# Wavepacket half
# ---------------------------------------------------------------------------

@dataclass
class ChannelWPRun:
    """A channeling wavepacket run: the bulk WPRun plus the radial channel."""
    name: str
    run_dir: Path
    base: K.WPRun               # t, T1, T2, pz, s3, s4, norm, sigma_z
    f_bore: np.ndarray          # fraction of |psi|^2 inside r_perp < R_in
    f_wall: np.ndarray          # fraction inside the jellium wall
    r_mean: np.ndarray          # <r_perp>(t), Bohr
    sigma_r: np.ndarray         # transverse spread of the packet, Bohr
    var_pz: np.ndarray          # sigma_pz^2(t) — the mechanism diagnostic
    e_total: np.ndarray         # Ha, conserved (no CAP)
    complete: bool
    steps_done: int
    steps_target: int

    @property
    def t(self) -> np.ndarray:
        return self.base.t

    @property
    def energy_drift_ev(self) -> float:
        """E_total(end) - E_total(0), eV. Must be ~0: no CAP => H is Hermitian."""
        if self.e_total.size < 2:
            return float("nan")
        return float((self.e_total[-1] - self.e_total[0]) * HA_TO_EV)

    @property
    def norm_drift(self) -> float:
        return float(abs(self.base.norm[-1] - self.base.norm[0]))

    @property
    def localisation_ev(self) -> np.ndarray:
        """(T1 - T2)(t) in eV. Constant at 1.2755 eV under free evolution."""
        return (self.base.T1 - self.base.T2) * HA_TO_EV

    def channeling_window(self, f_min: float = F_BORE_MIN,
                          transient_frac: float = TRANSIENT_FRAC) -> tuple[float, float]:
        """[t0, t1] over which the packet is MEASURED to be channeling.

        t1 = the last time f_bore has not yet fallen below f_min (i.e. the first
        breach ends the window); t0 = transient_frac of it. Falls back to the whole
        run if f_bore never breaches, which is the ideal outcome.
        """
        t = self.t
        breach = np.flatnonzero(self.f_bore < f_min)
        t1 = float(t[breach[0]]) if breach.size else float(t[-1])
        return (transient_frac * t1, t1)

    def full_window(self, transient_frac: float = 0.05) -> tuple[float, float]:
        """Whole run minus a short launch transient — the sensitivity comparison."""
        return (transient_frac * float(self.t[-1]), float(self.t[-1]))


def _obs_dir(run_dir: Path) -> Path:
    obs = Path(run_dir) / "raw" / "observables"
    if not obs.is_dir():
        raise FileNotFoundError(f"no observables under {run_dir}")
    return obs


def load_wp(name: str = "wp") -> ChannelWPRun:
    run_dir = WP_RESULTS / name
    obs = _obs_dir(run_dir)

    # K.load_wp_run hard-codes an extra "results" path level this layout does not
    # have, so the same steps are done inline. The DEFINITIONS still come from the
    # bulk module (unwrap_periodic, WPRun, fit_stopping) — only the lookup differs.
    mom = K._concat_segments(obs, "wp_momentum_stats")
    pos = K._concat_segments(obs, "wp_real_space_stats")
    df = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))

    t = df["time_au"].to_numpy()
    T1 = df["e_kin_ha"].to_numpy()
    px, py, pz = (df[c].to_numpy() for c in ("px_mean", "py_mean", "pz_mean"))
    T2 = 0.5 * (px**2 + py**2 + pz**2)

    if "z_mean_circ" not in df.columns:
        raise KeyError(
            "z_mean_circ missing from wp_real_space_stats.csv. The packet is "
            "launched 2 Bohr from the -z face and straddles it at t=0, so the "
            "naive z_mean is meaningless for this run — the circular centroid is "
            "not optional here."
        )
    s3 = K.unwrap_periodic(df["z_mean_circ"].to_numpy(), LZ)
    s3 = s3 - s3[0] + df["z_mean_circ"].to_numpy()[0]
    s4 = LAUNCH_Z + np.concatenate(
        [[0.0], np.cumsum(0.5 * (pz[1:] + pz[:-1]) * np.diff(t))])

    base = K.WPRun(
        run_dir=run_dir, box_length_z=LZ, t=t, step=df["step"].to_numpy(),
        T1=T1, T2=T2, pz=pz, s3=s3, s3_naive=df["z_mean"].to_numpy(), s4=s4,
        norm=df["norm_check_r"].to_numpy() if "norm_check_r" in df
             else df["norm_check"].to_numpy(),
        sigma_z=df["sigma_z_circ"].to_numpy(),
        parseval=df["norm_check_p"].to_numpy() if "norm_check_p" in df
                 else np.full(len(t), np.nan),
    )

    # --- the radial (channeling) channel -----------------------------------
    occ = K._concat_segments(obs, "wp_radial_occupancy")
    occ = occ.drop_duplicates(subset="step").set_index("step")
    steps = df["step"].to_numpy()

    def _align(col: str) -> np.ndarray:
        a = occ[col].reindex(steps).to_numpy()
        if np.isnan(a).any():          # written every step; a gap means truncation
            a = pd.Series(a).interpolate(limit_direction="both").to_numpy()
        return a

    f_bore = _align("f_bore")
    f_wall = _align("f_wall")
    r_mean = _align("r_mean")
    sig_r = _align("sigma_r")

    # --- energy conservation ------------------------------------------------
    try:
        en = K._concat_segments(obs, "observables")
        e_total = en["energy_total"].to_numpy()
    except (FileNotFoundError, KeyError):
        e_total = np.array([])

    steps_done = int(steps[-1])
    target = int(float(summary_of(run_dir).get("n_steps", N_STEPS)))
    return ChannelWPRun(
        name=name, run_dir=run_dir, base=base,
        f_bore=f_bore, f_wall=f_wall, r_mean=r_mean, sigma_r=sig_r,
        var_pz=df["sigma_pz2"].to_numpy() if "sigma_pz2" in df else np.full(len(t), np.nan),
        e_total=e_total, complete=steps_done >= target - 1,
        steps_done=steps_done, steps_target=target,
    )


# ---------------------------------------------------------------------------
# Classical half
# ---------------------------------------------------------------------------

@dataclass
class ChannelClassicalRun:
    name: str
    run_dir: Path
    base: K.ClassicalRun        # t, z (unwrapped), vz, T
    x: np.ndarray               # Bohr — must stay ~0 (channeling stability)
    y: np.ndarray
    force_z: np.ndarray         # Ha/Bohr, the instantaneous drag
    force_x: np.ndarray         # Ha/Bohr, must be ~0 by tube symmetry
    u_proj_bg: np.ndarray       # Ha
    e_total: np.ndarray         # Ha (electronic only)
    complete: bool
    steps_done: int
    steps_target: int

    @property
    def t(self) -> np.ndarray:
        return self.base.t

    @property
    def v_fraction(self) -> np.ndarray:
        return self.base.v_fraction

    @property
    def off_axis_max(self) -> float:
        """Largest transverse excursion, Bohr. Tube symmetry says ~0."""
        return float(np.max(np.hypot(self.x, self.y)))

    @property
    def conserved(self) -> np.ndarray:
        """E_electronic + KE_proj + U_proj_bg, Ha. Flat if force and potential agree."""
        if self.e_total.size != self.base.T.size:
            return np.array([])
        return self.e_total + self.base.T + self.u_proj_bg

    @property
    def conserved_drift_ev(self) -> float:
        c = self.conserved
        return float((c[-1] - c[0]) * HA_TO_EV) if c.size >= 2 else float("nan")


def load_classical(name: str = "classical") -> ChannelClassicalRun:
    run_dir = CL_RESULTS / name
    obs = _obs_dir(run_dir)

    trk = K._concat_segments(obs, "electron_track")
    z = K.unwrap_periodic(trk["z"].to_numpy(), LZ)
    z = z - z[0] + trk["z"].to_numpy()[0]
    base = K.ClassicalRun(
        run_dir=run_dir, t=trk["time_au"].to_numpy(), step=trk["step"].to_numpy(),
        z=z, vz=trk["vz"].to_numpy(), T=trk["ke_ion_ha"].to_numpy(),
    )

    # projectile.csv carries what electron_track.csv does not: U_proj_bg, the
    # force components, and proj_z_unwrapped (which cross-checks the unwrap above).
    try:
        pj = K._concat_segments(obs, "projectile")
        pj = pj.drop_duplicates(subset="step").set_index("step")
        steps = trk["step"].to_numpy()

        def _col(c, default=np.nan):
            if c not in pj.columns:
                return np.full(len(steps), default)
            a = pj[c].reindex(steps).to_numpy()
            return pd.Series(a).interpolate(limit_direction="both").to_numpy() \
                if np.isnan(a).any() else a

        u_bg = _col("energy_proj_bg_ideal", 0.0)
        fz, fx = _col("force_z", 0.0), _col("force_x", 0.0)
        # proj_z_unwrapped is written by the binary and needs no heuristic; prefer
        # it over the unwrap above when present.
        if "proj_z_unwrapped" in pj.columns:
            zu = _col("proj_z_unwrapped")
            if np.isfinite(zu).all():
                base.z = zu
    except (FileNotFoundError, KeyError):
        u_bg = np.zeros(len(base.t))
        fz = fx = np.full(len(base.t), np.nan)

    try:
        en = K._concat_segments(obs, "observables")
        e_total = en["energy_total"].to_numpy()
    except (FileNotFoundError, KeyError):
        e_total = np.array([])

    steps_done = int(trk["step"].to_numpy()[-1])
    target = int(float(summary_of(run_dir).get("n_steps", N_STEPS)))
    return ChannelClassicalRun(
        name=name, run_dir=run_dir, base=base,
        x=trk["x"].to_numpy(), y=trk["y"].to_numpy(),
        force_z=fz, force_x=fx, u_proj_bg=u_bg, e_total=e_total,
        complete=steps_done >= target - 1, steps_done=steps_done, steps_target=target,
    )


# ---------------------------------------------------------------------------
# Pairwise interaction ledger (both halves share a 12-column schema)
# ---------------------------------------------------------------------------

def load_interactions(run_dir: Path) -> pd.DataFrame:
    """interactions.csv with the closure residuals appended.

    The residual columns are the point: E_PP is the projectile SELF-Hartree, which
    exists only for the wavepacket (a classical external charge has no
    self-interaction in the KS Hamiltonian), and is the leading candidate for any
    residual classical/WP discrepancy in S.
    """
    df = K._concat_segments(_obs_dir(run_dir), "interactions")
    for c in ("e_ss", "e_pp", "e_ps", "e_sb", "e_pb", "e_bb"):
        if c in df:
            df[c + "_ev"] = df[c] * HA_TO_EV
    # Deltas from t=0: the absolute E_SB/E_PB/E_BB carry the charged-cell G=0
    # gauge and are NOT comparable across representations; their CHANGES are.
    for c in ("e_ss", "e_pp", "e_ps", "e_sb", "e_pb"):
        if c in df:
            df["d_" + c + "_ev"] = (df[c] - df[c].iloc[0]) * HA_TO_EV
    return df


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

@dataclass
class FreezeResult:
    """Is var(p) frozen? The MECHANISM behind a successful channeling result."""
    var_start: float
    var_end: float
    var_free: float
    growth_pct: float
    localisation_start_ev: float
    localisation_end_ev: float
    localisation_drift_ev: float
    frozen: bool
    window: tuple[float, float]

    def summary(self) -> str:
        return (f"var(p_z): {self.var_start:.5f} -> {self.var_end:.5f} "
                f"({self.growth_pct:+.1f} %; free value {self.var_free:.5f}) | "
                f"(T1-T2): {self.localisation_start_ev:.3f} -> "
                f"{self.localisation_end_ev:.3f} eV "
                f"(drift {self.localisation_drift_ev:+.3f} eV) | "
                f"{'FROZEN' if self.frozen else 'NOT frozen'}")


def var_p_freeze(run: ChannelWPRun, window: tuple[float, float] | None = None,
                 tol_pct: float = 10.0) -> FreezeResult:
    """Growth of var(p_z) over `window` (default: the channeling window).

    var(p) is CONSERVED under free evolution, so any growth is interaction with
    the bath. In bulk it grew and contaminated T1; the channeling claim is that it
    does not grow here. `tol_pct` is what counts as frozen.
    """
    w = window or run.channeling_window()
    m = (run.t >= w[0]) & (run.t <= w[1])
    if m.sum() < 2:
        m = np.ones_like(run.t, dtype=bool)
    v = run.var_pz[m]
    loc = run.localisation_ev[m]
    growth = 100.0 * (v[-1] - v[0]) / v[0] if np.isfinite(v[0]) and v[0] != 0 else float("nan")
    return FreezeResult(
        var_start=float(v[0]), var_end=float(v[-1]), var_free=VAR_P_FREE,
        growth_pct=float(growth),
        localisation_start_ev=float(loc[0]), localisation_end_ev=float(loc[-1]),
        localisation_drift_ev=float(loc[-1] - loc[0]),
        frozen=bool(np.isfinite(growth) and abs(growth) <= tol_pct),
        window=w,
    )


@dataclass
class ChannelingResult:
    """Did the packet stay in the bore? The PREMISE of the study."""
    f_bore_start: float
    f_bore_end: float
    f_bore_min: float
    f_wall_end: float
    r_mean_start: float
    r_mean_end: float
    t_breach: float             # first time f_bore < F_BORE_MIN (inf if never)
    t_2sigma_at_wall: float     # the free-dispersion estimate, for comparison
    channeling: bool

    def summary(self) -> str:
        br = "never" if not np.isfinite(self.t_breach) else f"t = {self.t_breach:.1f} a.u."
        return (f"f_bore: {self.f_bore_start:.3f} -> {self.f_bore_end:.3f} "
                f"(min {self.f_bore_min:.3f}); f_wall(end) = {self.f_wall_end:.3f}; "
                f"<r_perp>: {self.r_mean_start:.2f} -> {self.r_mean_end:.2f} Bohr; "
                f"breach {br} (free-dispersion estimate "
                f"{self.t_2sigma_at_wall:.1f}) | "
                f"{'CHANNELING' if self.channeling else 'left the bore'}")


def channeling_check(run: ChannelWPRun, f_min: float = F_BORE_MIN) -> ChannelingResult:
    breach = np.flatnonzero(run.f_bore < f_min)
    t_breach = float(run.t[breach[0]]) if breach.size else float("inf")
    return ChannelingResult(
        f_bore_start=float(run.f_bore[0]), f_bore_end=float(run.f_bore[-1]),
        f_bore_min=float(np.min(run.f_bore)), f_wall_end=float(run.f_wall[-1]),
        r_mean_start=float(run.r_mean[0]), r_mean_end=float(run.r_mean[-1]),
        t_breach=t_breach, t_2sigma_at_wall=T_2SIGMA_AT_WALL,
        channeling=bool(not breach.size or t_breach > 0.5 * float(run.t[-1])),
    )


# ---------------------------------------------------------------------------
# The comparison — the deliverable
# ---------------------------------------------------------------------------

@dataclass
class Comparison:
    window: tuple[float, float]
    wp_fits: dict[str, K.StoppingFit]        # S_13, S_14, S_23, S_24
    cl_same_window: K.StoppingFit            # classical over the SAME window
    cl_initial_drag: K.StoppingFit           # classical, v >= 0.85 v0 (the rule)
    wp_full_window: dict[str, K.StoppingFit] # sensitivity: the whole run
    freeze: FreezeResult
    channel: ChannelingResult
    wp_energy_drift_ev: float
    cl_conserved_drift_ev: float
    cl_off_axis_max: float
    verdict: str = ""
    aim_met: bool = False
    agreement_pct: float = float("nan")

    def table(self) -> pd.DataFrame:
        rows = []
        for key, f in self.wp_fits.items():
            rows.append(dict(estimator=key, half="wp", S_ev_per_bohr=f.S_ev_per_bohr,
                             uncertainty=f.uncertainty, r2=f.r2, n=f.n_points,
                             t0=f.t_window[0], t1=f.t_window[1], mean_v=f.mean_v))
        for key, f in (("S_cl_same_window", self.cl_same_window),
                       ("S_cl_initial_drag", self.cl_initial_drag)):
            rows.append(dict(estimator=key, half="classical", S_ev_per_bohr=f.S_ev_per_bohr,
                             uncertainty=f.uncertainty, r2=f.r2, n=f.n_points,
                             t0=f.t_window[0], t1=f.t_window[1], mean_v=f.mean_v))
        return pd.DataFrame(rows)


# S_24 = -d(<p>^2/2m)/d(integral <p> dt) is the DRIFT-vs-DRIFT combination: both
# sides built from <p_z> alone, so it is the definition that is a stopping power
# whether or not var(p) is frozen. The others are reported and are the
# cross-checks; this is the headline.
PRIMARY_ESTIMATOR = "S_24"


def compare(wp: ChannelWPRun, cl: ChannelClassicalRun,
            f_min: float = F_BORE_MIN,
            agreement_tol_pct: float = 20.0) -> Comparison:
    """Fit both halves over the measured channeling window and judge the aim.

    `agreement_tol_pct` is the band within which the two S values count as
    reproducing each other. 20 % is deliberately generous: the deliverable is a
    demonstration that a KS-orbital definition lands ON the classical curve rather
    than a factor ~2 away (the bulk result), not a percent-level metrology claim.
    """
    w = wp.channeling_window(f_min=f_min)

    wp_fits = K.fit_all_wp(wp.base, w[0], w[1])
    wp_full = K.fit_all_wp(wp.base, *wp.full_window())
    cl_same = K.fit_stopping(cl.base.z, cl.base.T, cl.base.t, w[0], w[1],
                             label="S_cl (same window as the WP fit)", v=cl.base.vz)
    cl_drag = K.fit_classical_early(cl.base)

    freeze = var_p_freeze(wp, window=w)
    channel = channeling_check(wp, f_min=f_min)

    s_wp = wp_fits[PRIMARY_ESTIMATOR].S_ev_per_bohr
    s_cl = cl_same.S_ev_per_bohr
    agree = 100.0 * abs(s_wp - s_cl) / abs(s_cl) if s_cl not in (0.0,) and np.isfinite(s_cl) else float("nan")

    ok_s = bool(np.isfinite(agree) and agree <= agreement_tol_pct)
    aim = ok_s and freeze.frozen and channel.channeling

    bits = [
        f"{PRIMARY_ESTIMATOR} = {s_wp:.2f} eV/Bohr vs classical {s_cl:.2f} eV/Bohr "
        f"over t = {w[0]:.1f}-{w[1]:.1f} a.u.  (differ by {agree:.0f} %"
        f"{'' if ok_s else f', outside the {agreement_tol_pct:.0f} % band'})",
        channel.summary(),
        freeze.summary(),
    ]
    if aim:
        head = ("AIM MET: the channeling wavepacket reproduces the classical stopping "
                "power, and it does so for the stated reason — the packet stayed in "
                "the bore and var(p) stayed frozen.")
    elif ok_s and not (freeze.frozen and channel.channeling):
        head = ("AIM PARTLY MET: the two stopping powers agree, but the mechanism "
                "check failed — the agreement is not yet explained by channeling, "
                "so treat it as unexplained rather than confirmed.")
    elif not ok_s and freeze.frozen and channel.channeling:
        head = ("AIM NOT MET: the packet did channel cleanly and var(p) stayed "
                "frozen, yet the two stopping powers still differ. That is the "
                "interesting failure — it points at something other than "
                "interaction-driven momentum spreading (start with E_PP, the WP "
                "self-Hartree, in interactions.csv).")
    else:
        head = ("AIM NOT MET: neither the stopping powers nor the channeling "
                "premise held. Check f_bore(t) first — if the packet left the bore, "
                "this run did not test the hypothesis at all.")

    return Comparison(
        window=w, wp_fits=wp_fits, cl_same_window=cl_same, cl_initial_drag=cl_drag,
        wp_full_window=wp_full, freeze=freeze, channel=channel,
        wp_energy_drift_ev=wp.energy_drift_ev,
        cl_conserved_drift_ev=cl.conserved_drift_ev,
        cl_off_axis_max=cl.off_axis_max,
        verdict=head + "\n  " + "\n  ".join(bits), aim_met=aim, agreement_pct=agree,
    )


def load_pair(wp_name: str = "wp", cl_name: str = "classical"):
    return load_wp(wp_name), load_classical(cl_name)


if __name__ == "__main__":
    wp, cl = load_pair()
    c = compare(wp, cl)
    print(c.table().to_string(index=False))
    print()
    print(c.verdict)
