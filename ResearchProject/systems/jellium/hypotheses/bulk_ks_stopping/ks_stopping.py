"""Stopping-power extraction for the bulk-jellium KS-orbital twin pair.

Plan: docs/plans/bulk-jellium-ks-stopping.md

This is the deterministic engine: it does the arithmetic and nothing else. The
notebooks narrate; this module computes, so both run notebooks and any later
cross-run comparison get bit-identical numbers from one implementation.

THE FOUR DEFINITIONS
--------------------
Kinetic energy of the projectile (wavepacket run), per step, from
``wp_momentum_stats.csv``:

    T1 = <p^2>/2m          -> column ``e_kin_ha``       (full orbital KE)
    T2 = <p>^2/2m          -> 0.5*(px^2+py^2+pz^2)      (drift only)

T1 - T2 = (3/2) * sum_d sigma_pd^2 is the localisation + scattering energy;
at t=0 it equals 3/(8 sigma^2) = 2.551 eV for sigma_WP = 2 Bohr. Its CHANGE
along the trajectory is the momentum-broadening contribution to apparent
stopping, which is the physics contrast the twin pair is built to expose.

Position of the projectile, per step:

    s3 = density centroid of the WP KS orbital, from ``wp_real_space_stats.csv``
    s4 = z0 + integral of <p_z> dt   (cumulative trapezoid of ``pz_mean``)

s3 USES ``z_mean_circ``, NOT ``z_mean``. The naive integral of z|psi|^2 is
discontinuous across a periodic face and slides smoothly to a WRONG value while
the packet straddles it; the circular (phase) estimator is exact in a periodic
cell. Both are loaded so the notebook can show where they part company.

EHRENFEST IDENTITY. The wavepacket run has no ions and no absorbing potential,
so the KS Hamiltonian is purely local and d<z>/dt = <p_z>/m holds EXACTLY.
s3 and s4 must therefore agree to numerical precision: ``ehrenfest_residual``
is a validation diagnostic, not a physical result. A growing residual means the
WP orbital is leaking norm into the bath (or the packet has wrapped and the
naive centroid was used by mistake).

Stopping power is then, for each (i, j) in {1,2} x {3,4},

    S_ij = -d(T_i)/d(s_j)

fitted by ordinary least squares over the analysis window and reported in
eV/Bohr. The classical twin gives the single unambiguous reference
S_cl = -d(1/2 m v^2)/dz from ``electron_track.csv``.

UNCERTAINTY. Three numbers are reported and they mean different things:
  * ``stderr``      - OLS standard error of the slope: how well a straight line
                      is determined by the points, nothing more.
  * ``window_syst`` - half-spread of the slope as the fit window edges are moved
                      +/- 3 a.u. This is the honest dominant uncertainty; it is
                      usually much larger than stderr and is what should be
                      quoted.
  * ``r2``          - if this is far below 1, T(s) is not linear over the window
                      and a single S does not describe it. Report the curvature
                      rather than the slope in that case.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

HA_TO_EV = 27.211386


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _concat_segments(run_dir: Path, stem: str) -> pd.DataFrame:
    """Concatenate ``<stem>.csv`` and any ``<stem>.fromNNN.csv`` resume segments.

    Segment CSVs come from the resume path (.claude/rules/final-timestep-
    checkpoint.md). They are ordered by step and de-duplicated on the boundary
    step, which appears in both the segment that ended there and the one that
    restarted from it.
    """
    files = sorted(run_dir.glob(f"{stem}*.csv"))
    if not files:
        raise FileNotFoundError(f"no {stem}*.csv under {run_dir}")
    frames = [pd.read_csv(f, comment="#") for f in files]
    df = pd.concat(frames, ignore_index=True)
    df = df.sort_values("step").drop_duplicates(subset="step", keep="last")
    return df.reset_index(drop=True)


def unwrap_periodic(z: np.ndarray, box_length: float) -> np.ndarray:
    """Unwrap a periodic coordinate into a continuous trajectory.

    The circular centroid is reported in (-L/2, +L/2]; a projectile crossing the
    face reappears at the other end as a jump of one box length. Unwrapping in
    PHASE (not in z) is what makes this exact: convert to theta = 2 pi z / L,
    np.unwrap, convert back.
    """
    theta = 2.0 * np.pi * np.asarray(z, dtype=float) / box_length
    return np.unwrap(theta) * box_length / (2.0 * np.pi)


@dataclass
class WPRun:
    """Per-step kinematics of the wavepacket projectile."""
    run_dir: Path
    box_length_z: float
    t: np.ndarray            # a.u.
    step: np.ndarray
    T1: np.ndarray           # Ha, <p^2>/2m
    T2: np.ndarray           # Ha, <p>^2/2m
    pz: np.ndarray           # Bohr^-1
    s3: np.ndarray           # Bohr, circular centroid, unwrapped
    s3_naive: np.ndarray     # Bohr, raw integral of z|psi|^2 (for comparison)
    s4: np.ndarray           # Bohr, integral of <p_z> dt
    norm: np.ndarray         # WP orbital norm from REAL space (should stay ~1)
    sigma_z: np.ndarray      # Bohr, circular spread (periodic-safe width)
    parseval: np.ndarray     # momentum-space norm_check: a large FFT constant,
                             # NOT 1 — see load_wp_run.

    @property
    def ehrenfest_residual(self) -> np.ndarray:
        """s3 - s4. Must be ~0; see the module docstring."""
        return self.s3 - self.s4

    @property
    def localisation_energy(self) -> np.ndarray:
        """T1 - T2 in eV: the momentum-width (localisation + scattering) energy."""
        return (self.T1 - self.T2) * HA_TO_EV


def load_wp_run(run_dir: str | Path, box_length_z: float, z0: float) -> WPRun:
    run_dir = Path(run_dir)
    obs = run_dir / "results" / "raw" / "observables"
    mom = _concat_segments(obs, "wp_momentum_stats")
    pos = _concat_segments(obs, "wp_real_space_stats")

    # The two files are written on the same cadence, but merge on step rather
    # than assume alignment — a resume can leave them off by a boundary row.
    df = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))

    t = df["time_au"].to_numpy()
    T1 = df["e_kin_ha"].to_numpy()
    px, py, pz = (df[c].to_numpy() for c in ("px_mean", "py_mean", "pz_mean"))
    T2 = 0.5 * (px**2 + py**2 + pz**2)

    if "z_mean_circ" not in df.columns:
        raise KeyError(
            "z_mean_circ missing from wp_real_space_stats.csv — this run predates "
            "the periodic-aware centroid (added 2026-07-30). The naive z_mean "
            "cannot be used near a cell face; re-run or restrict the analysis."
        )
    s3 = unwrap_periodic(df["z_mean_circ"].to_numpy(), box_length_z)
    # Pin the unwrapped branch to the launch position: np.unwrap fixes the
    # increments, not the absolute offset.
    s3 = s3 - s3[0] + df["z_mean_circ"].to_numpy()[0]

    # s4: integrate <p_z> dt. m = 1 (electron) so v_z = <p_z> exactly.
    s4 = z0 + np.concatenate([[0.0], np.cumsum(0.5 * (pz[1:] + pz[:-1]) * np.diff(t))])

    return WPRun(
        run_dir=run_dir, box_length_z=box_length_z,
        t=t, step=df["step"].to_numpy(), T1=T1, T2=T2, pz=pz,
        s3=s3, s3_naive=df["z_mean"].to_numpy(), s4=s4,
        # BOTH files carry a `norm_check`, and they mean DIFFERENT things. The
        # real-space one (suffix _r) is the physical orbital norm and must be 1.
        # The momentum-space one (_p) is an unnormalised Parseval sum over the
        # FFT grid — here ~4.9e7 — and gating it on 1 would report a healthy run
        # as catastrophically broken. Take the real-space column; keep the other
        # only as a constancy diagnostic.
        norm=df["norm_check_r"].to_numpy() if "norm_check_r" in df
             else df["norm_check"].to_numpy(),
        sigma_z=df["sigma_z_circ"].to_numpy(),
        parseval=df["norm_check_p"].to_numpy() if "norm_check_p" in df
                 else np.full(len(t), np.nan),
    )


@dataclass
class ClassicalRun:
    """Per-step kinematics of the classical projectile."""
    run_dir: Path
    t: np.ndarray
    step: np.ndarray
    z: np.ndarray            # Bohr
    vz: np.ndarray           # Bohr/atu
    T: np.ndarray            # Ha, 1/2 m v^2

    @property
    def v_fraction(self) -> np.ndarray:
        """v_z / v_z(0) — the deceleration the light-projectile rule warns about."""
        return self.vz / self.vz[0]


def load_classical_run(run_dir: str | Path, box_length_z: float) -> ClassicalRun:
    run_dir = Path(run_dir)
    obs = run_dir / "results" / "raw" / "observables"
    df = _concat_segments(obs, "electron_track")
    z = unwrap_periodic(df["z"].to_numpy(), box_length_z)
    z = z - z[0] + df["z"].to_numpy()[0]
    return ClassicalRun(
        run_dir=run_dir,
        t=df["time_au"].to_numpy(), step=df["step"].to_numpy(),
        z=z, vz=df["vz"].to_numpy(), T=df["ke_ion_ha"].to_numpy(),
    )


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

@dataclass
class StoppingFit:
    label: str
    S_ev_per_bohr: float
    stderr: float                 # OLS slope standard error, eV/Bohr
    window_syst: float = 0.0      # window-sensitivity systematic, eV/Bohr
    r2: float = 0.0
    n_points: int = 0
    t_window: tuple[float, float] = (0.0, 0.0)
    s_range: tuple[float, float] = (0.0, 0.0)
    mean_v: float = 0.0
    # Raw fit arrays, kept so the notebook can draw the fit and its residuals.
    s_fit: np.ndarray = field(default_factory=lambda: np.array([]), repr=False)
    T_fit: np.ndarray = field(default_factory=lambda: np.array([]), repr=False)
    T_model: np.ndarray = field(default_factory=lambda: np.array([]), repr=False)

    @property
    def uncertainty(self) -> float:
        """The number to quote: stderr and window systematic in quadrature."""
        return math.hypot(self.stderr, self.window_syst)

    def summary(self) -> str:
        return (f"{self.label}: S = {self.S_ev_per_bohr:.2f} +/- "
                f"{self.uncertainty:.2f} eV/Bohr "
                f"(stat {self.stderr:.2f}, syst {self.window_syst:.2f}; "
                f"r2 = {self.r2:.3f}, n = {self.n_points}, "
                f"t = {self.t_window[0]:.1f}-{self.t_window[1]:.1f} a.u.)")


def _ols(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, np.ndarray]:
    """Slope, its standard error, r^2 and the model values. No dependencies."""
    n = len(x)
    if n < 3:
        return float("nan"), float("nan"), float("nan"), np.full(n, np.nan)
    xm, ym = x.mean(), y.mean()
    sxx = np.sum((x - xm) ** 2)
    if sxx <= 0:
        return float("nan"), float("nan"), float("nan"), np.full(n, np.nan)
    slope = np.sum((x - xm) * (y - ym)) / sxx
    intercept = ym - slope * xm
    model = intercept + slope * x
    resid = y - model
    ss_res = np.sum(resid ** 2)
    ss_tot = np.sum((y - ym) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    stderr = math.sqrt(ss_res / (n - 2) / sxx) if n > 2 else float("nan")
    return slope, stderr, r2, model


def fit_stopping(s: np.ndarray, T: np.ndarray, t: np.ndarray,
                 t0: float, t1: float, label: str,
                 v: np.ndarray | None = None,
                 window_scan: float = 3.0) -> StoppingFit:
    """S = -dT/ds by OLS over t in [t0, t1], with a window-sensitivity systematic.

    ``T`` in Hartree, ``s`` in Bohr; the returned S is in eV/Bohr.

    The systematic moves BOTH window edges independently by +/- ``window_scan``
    a.u. and takes half the full spread of the resulting slopes. That is the
    dominant uncertainty in practice: where the transient ends and where
    interference begins are judgement calls, and this prices them.
    """
    m = (t >= t0) & (t <= t1)
    s_w, T_w = s[m], T[m]
    slope, stderr, r2, model = _ols(s_w, T_w)
    S = -slope * HA_TO_EV

    slopes = []
    for d0 in (-window_scan, 0.0, window_scan):
        for d1 in (-window_scan, 0.0, window_scan):
            a, b = t0 + d0, t1 + d1
            if b - a < 5.0:
                continue
            mm = (t >= a) & (t <= b)
            if mm.sum() < 10:
                continue
            sl, _, _, _ = _ols(s[mm], T[mm])
            if np.isfinite(sl):
                slopes.append(-sl * HA_TO_EV)
    syst = 0.5 * (max(slopes) - min(slopes)) if len(slopes) > 1 else 0.0

    return StoppingFit(
        label=label, S_ev_per_bohr=S, stderr=abs(stderr * HA_TO_EV),
        window_syst=syst, r2=r2, n_points=int(m.sum()),
        t_window=(t0, t1),
        s_range=(float(s_w[0]), float(s_w[-1])) if m.sum() else (0.0, 0.0),
        mean_v=float(np.mean(v[m])) if v is not None and m.sum() else 0.0,
        s_fit=s_w, T_fit=T_w, T_model=model,
    )


def fit_all_wp(run: WPRun, t0: float, t1: float) -> dict[str, StoppingFit]:
    """The four (KE definition) x (position definition) combinations."""
    out = {}
    for ti, (Tname, T) in enumerate((("T1", run.T1), ("T2", run.T2)), start=1):
        for si, (sname, s) in enumerate((("s3", run.s3), ("s4", run.s4)), start=3):
            key = f"S_{ti}{si}"
            desc = {("T1", "s3"): "<p^2>/2m vs centroid",
                    ("T1", "s4"): "<p^2>/2m vs integral<p>dt",
                    ("T2", "s3"): "<p>^2/2m vs centroid",
                    ("T2", "s4"): "<p>^2/2m vs integral<p>dt"}[(Tname, sname)]
            out[key] = fit_stopping(s, T, run.t, t0, t1,
                                    label=f"{key}  ({desc})", v=run.pz)
    return out


def fit_classical(run: ClassicalRun, t0: float, t1: float) -> StoppingFit:
    return fit_stopping(run.z, run.T, run.t, t0, t1,
                        label="S_cl  (classical 1/2 m v^2 vs z)", v=run.vz)


def fit_classical_early(run: ClassicalRun, v_frac: float = 0.85) -> StoppingFit:
    """Initial drag over the near-constant-velocity window.

    Required by .claude/rules/light-projectile-stopping.md: a light Ehrenfest
    projectile decelerates, so a full-trajectory regression averages S over every
    velocity between v0 and the final one instead of giving S AT v0. This
    restricts to v >= v_frac * v0, widening the threshold if too few points
    survive.
    """
    for frac in (v_frac, 0.70, 0.50):
        m = run.v_fraction >= frac
        if m.sum() >= 30:
            t_sel = run.t[m]
            return fit_stopping(run.z, run.T, run.t, t_sel.min(), t_sel.max(),
                                label=f"S_cl,initial (v >= {frac:.2f} v0)",
                                v=run.vz, window_scan=1.0)
    return fit_stopping(run.z, run.T, run.t, run.t[0], run.t[-1],
                        label="S_cl,initial (FALLBACK: full range)", v=run.vz)


# ---------------------------------------------------------------------------
# Energy ledger
# ---------------------------------------------------------------------------

def load_energies(run_dir: str | Path) -> pd.DataFrame:
    """observables.csv with every energy component, plus derived columns."""
    obs = Path(run_dir) / "results" / "raw" / "observables"
    df = _concat_segments(obs, "observables")
    ecols = [c for c in df.columns if c.startswith("energy_")]
    for c in ecols:
        df[c + "_ev"] = df[c] * HA_TO_EV
    if "energy_total" in df:
        df["delta_e_total_ev"] = (df["energy_total"] - df["energy_total"].iloc[0]) * HA_TO_EV
    return df


def conservation_check(df: pd.DataFrame,
                       projectile_ke_loss_ev: float | None = None) -> dict[str, float]:
    """Energy bookkeeping — the numerical figure of merit (Tier B3).

    There is no absorbing potential anywhere in this study, so nothing removes
    energy from the simulation and the bookkeeping must close. What "close"
    means differs between the two halves of the twin pair:

    * **Wavepacket run** — the projectile IS an electron orbital, so it is inside
      ``energy_total``. The system is closed and ``energy_total`` must be
      CONSTANT. Its drift is the pure numerical error, and every stopping number
      inherits it.

    * **Classical run** — VERIFIED 2026-07-30: INQ leaves ``energy_ion_kinetic``
      at zero (the Ehrenfest ion propagator does not feed it back into
      ``energy::total()``). ``energy_total`` is therefore the ELECTRONIC energy
      only, and it is *supposed* to rise — by exactly the kinetic energy the
      projectile gives up. Reading its increase as "drift" would be wrong.
      The meaningful test is the CLOSURE

          dE_electronic  ==  -dT_projectile

      which for this run holds to 0.049 eV out of a 22.58 eV transfer (0.22 %).
      Pass ``projectile_ke_loss_ev`` to have that computed.
    """
    if "energy_total" not in df:
        return {}
    e = df["energy_total"].to_numpy() * HA_TO_EV
    out = {
        "e_total_initial_ev": float(e[0]),
        "e_total_final_ev": float(e[-1]),
        "drift_ev": float(e[-1] - e[0]),
        "max_excursion_ev": float(np.max(np.abs(e - e[0]))),
    }
    if projectile_ke_loss_ev is not None:
        gain = float(e[-1] - e[0])
        out["projectile_ke_loss_ev"] = float(projectile_ke_loss_ev)
        out["electronic_gain_ev"] = gain
        out["closure_mismatch_ev"] = abs(gain - projectile_ke_loss_ev)
        out["closure_mismatch_pct"] = (
            100.0 * abs(gain - projectile_ke_loss_ev) / abs(projectile_ke_loss_ev)
            if projectile_ke_loss_ev else float("nan"))
    return out


# ---------------------------------------------------------------------------
# Pairwise interaction energies (.claude/rules/decomposed-interaction-energies.md)
# ---------------------------------------------------------------------------

# norm_proj threshold below which the classical projectile's Gaussian cloud is
# considered clipped by the +z box face. MEASURED, not guessed (sigma=3
# r_s=5.702, 2026-08-01): with norm_proj == 1 exactly, E_PP is bit-exactly
# constant; at >= 1-1e-9 it holds to 1.8e-11 Ha; relaxing to 1-1e-6 admits the
# clipping shoulder. The first few rows sit a few 1e-9 under 1.0 purely from
# discretising the Gaussian on the grid, which is why only the CONTIGUOUS TAIL
# counts as clipping. Identical logic to scripts/verify_interactions_closure.py.
CLIP_TOL = 1e-9


@dataclass
class Interactions:
    """Per-step P/S/B pairwise electrostatic decomposition of one run half.

    Terms (Ha), with P = projectile, S = system/bath electrons, B = background:

        e_ss = 1/2 int n_S phi_S     bath-bath
        e_pp = 1/2 int n_P phi_P     projectile SELF-Hartree  <- the quantum residual
        e_ps =     int n_S phi_P     projectile-bath          <- what stops it

    BULK ONLY: the background is uniform, so poisson(n+) is pure G=0, which INQ
    drops. phi+ is IDENTICALLY ZERO and e_sb = e_pb = e_bb = 0. They are carried
    as columns anyway so the schema matches the slab systems.
    """
    run_dir: Path
    half: str                    # "wp" | "classical"
    t: np.ndarray                # a.u.
    step: np.ndarray
    e_ss: np.ndarray             # Ha
    e_pp: np.ndarray             # Ha
    e_ps: np.ndarray             # Ha
    e_sb: np.ndarray             # Ha, 0 in bulk
    e_pb: np.ndarray             # Ha, 0 in bulk
    e_bb: np.ndarray             # Ha, 0 in bulk
    norm: np.ndarray             # norm_wp (wp) | norm_proj (classical)
    proj_z: np.ndarray | None    # classical only, Bohr (raw, may wrap)
    closure: np.ndarray          # Ha, residual against INQ's own Hartree energy

    @property
    def clip_index(self) -> int | None:
        """First index of the CONTIGUOUS TAIL where the projectile cloud is clipped.

        ``None`` for a wavepacket half (no rigid cloud to clip) and for a
        classical half that never reaches the face.
        """
        if self.half != "classical":
            return None
        clean = self.norm >= 1.0 - CLIP_TOL
        tail = len(clean)
        while tail > 0 and not clean[tail - 1]:
            tail -= 1
        return None if tail >= len(clean) else int(tail)

    @property
    def clip_time(self) -> float:
        """Time (a.u.) at which clipping starts; ``inf`` if it never does.

        This is a HARD UPPER BOUND on any fit window using this run — it is not
        in any config header, because it is a property of the trajectory meeting
        the box, not of the physics being modelled.
        """
        i = self.clip_index
        return float("inf") if i is None else float(self.t[i])

    def in_window(self, t0: float, t1: float) -> np.ndarray:
        """Boolean mask for [t0, t1], truncated at the clipping onset."""
        return (self.t >= t0) & (self.t <= min(t1, self.clip_time))


def load_interactions(run_dir: str | Path, half: str) -> Interactions:
    """Load ``interactions.csv`` (+ resume segments) for one half of a pair.

    The closure residual is computed against whichever INQ scalar that half
    exposes, because the two representations put the projectile in DIFFERENT
    ledger terms:

        classical : e_ss                  == energy_hartree   (col e_hartree_inq)
        wp        : e_ss + e_ps + e_pp    == energy_hartree   (col e_hartree_check)

    For the WP half ``e_ss + e_ps + e_pp == e_hartree_check`` is an algebraic
    identity of ``compute_coulomb_wp`` (e_ss is *defined* as e_hartree_check -
    cross + e_pp), so it is exact by construction and proves nothing on its own.
    The informative comparison is against INQ's own ``energy_hartree`` from
    observables.csv, which is what this function reports.
    """
    if half not in ("wp", "classical"):
        raise ValueError(f"half must be 'wp' or 'classical', got {half!r}")
    run_dir = Path(run_dir)
    obs = run_dir / "results" / "raw" / "observables"
    df = _concat_segments(obs, "interactions")

    for c in ("e_ss", "e_pp", "e_ps"):
        if c not in df.columns:
            raise KeyError(
                f"{c} missing from interactions.csv under {obs} — this run predates "
                "the pairwise decomposition (.claude/rules/decomposed-interaction-"
                "energies.md). Re-run with the wired run.cpp to gain it."
            )

    # Cross-check against INQ's own Hartree energy. observables.csv is written on
    # a different cadence, so merge on step and compare only the shared rows.
    inq_col = "e_hartree_check" if half == "wp" else "e_hartree_inq"
    ours = df[inq_col].to_numpy() if inq_col in df.columns else np.full(len(df), np.nan)
    try:
        en = _concat_segments(obs, "observables")
        m = df[["step"]].merge(en[["step", "energy_hartree"]], on="step", how="left")
        closure = ours - m["energy_hartree"].to_numpy()
    except (FileNotFoundError, KeyError):
        closure = np.full(len(df), np.nan)

    zero = np.zeros(len(df))
    norm_col = "norm_wp" if half == "wp" else "norm_proj"
    return Interactions(
        run_dir=run_dir, half=half,
        t=df["time_au"].to_numpy(), step=df["step"].to_numpy(),
        e_ss=df["e_ss"].to_numpy(), e_pp=df["e_pp"].to_numpy(),
        e_ps=df["e_ps"].to_numpy(),
        e_sb=df["e_sb"].to_numpy() if "e_sb" in df else zero,
        e_pb=df["e_pb"].to_numpy() if "e_pb" in df else zero,
        e_bb=df["e_bb"].to_numpy() if "e_bb" in df else zero,
        norm=df[norm_col].to_numpy() if norm_col in df else np.full(len(df), np.nan),
        proj_z=df["proj_z"].to_numpy() if "proj_z" in df else None,
        closure=closure,
    )


# ---------------------------------------------------------------------------
# Phase space — the classical/WP comparison on one footing
# ---------------------------------------------------------------------------

def local_stopping(s: np.ndarray, T: np.ndarray, half_width: int = 12
                   ) -> np.ndarray:
    """S(s) = -dT/ds as a CENTRED ROLLING OLS SLOPE, in eV/Bohr.

    A raw finite difference of T(s) is unusable here: the per-step energy change
    is ~1e-5 Ha against a step in s of ~0.1 Bohr, so point-to-point noise swamps
    the trend. A rolling least-squares slope over (2*half_width+1) samples is the
    smallest honest estimator — it is exactly the fit ``fit_stopping`` performs,
    evaluated locally instead of over one wide window.

    Edges are filled with the nearest interior value rather than extrapolated, so
    the first/last ``half_width`` points repeat and MUST NOT be read as structure.
    """
    s = np.asarray(s, dtype=float)
    T = np.asarray(T, dtype=float)
    n = len(s)
    out = np.full(n, np.nan)
    if n < 2 * half_width + 1:
        return out
    for i in range(half_width, n - half_width):
        sl = slice(i - half_width, i + half_width + 1)
        x, y = s[sl], T[sl]
        xc = x - x.mean()
        denom = float(xc @ xc)
        if denom <= 0:
            continue
        out[i] = -float(xc @ (y - y.mean())) / denom * HA_TO_EV
    out[:half_width] = out[half_width]
    out[n - half_width:] = out[n - half_width - 1]
    return out


@dataclass
class PairPhase:
    """Classical and wavepacket projectile kinematics on a common footing.

    WHAT IS COMPARABLE TO WHAT. The classical projectile has one kinetic energy,
    T_cl = 1/2 m v^2. The wavepacket has two, and only one of them is its
    counterpart:

        T2 = <p>^2/2m   DRIFT kinetic energy   -> the classical analogue
        T1 = <p^2>/2m   TOTAL orbital KE       -> drift + internal
        T1 - T2         INTERNAL (momentum-width) energy, NO classical counterpart

    Comparing T_cl against T1 would charge the wavepacket for its own zero-point
    spread and overstate its energy by ~2.5 eV at sigma_WP = 2 Bohr before it has
    moved at all. The phase-space comparison therefore uses (z, v) and (z, T2),
    and reports T1 - T2 separately as the quantum-only channel.

    Velocities are directly comparable without conversion: the projectile is an
    electron, m = 1 in atomic units, so v = <p_z> exactly.
    """
    family: str
    wp: WPRun
    cl: ClassicalRun
    ix_wp: Interactions | None = None
    ix_cl: Interactions | None = None
    v0: float = field(init=False)

    def __post_init__(self) -> None:
        # Both halves are launched at the same k0, so either is a valid v0; take
        # the classical one, which is exact by construction (the WP centroid
        # momentum can differ in the last digits from grid discretisation).
        self.v0 = float(self.cl.vz[0])

    def divergence(self, frac: float = 0.05) -> dict[str, float]:
        """Where the two velocity histories part company.

        Defined as the first time |v_wp - v_cl| exceeds ``frac`` of v0, on the
        COMMON time grid. Returns the time, both positions, and both velocities
        there. All-NaN if they never diverge by that much.

        This is a descriptive marker for the phase portrait, not a fitted
        physical quantity — the threshold is a reading aid and is stated on the
        figure so it cannot be mistaken for a result.
        """
        t = self.cl.t
        v_wp = np.interp(t, self.wp.t, self.wp.pz)
        d = np.abs(v_wp - self.cl.vz)
        hit = np.nonzero(d > frac * abs(self.v0))[0]
        if len(hit) == 0:
            return dict(t=float("nan"), z_cl=float("nan"), z_wp=float("nan"),
                        v_cl=float("nan"), v_wp=float("nan"), frac=frac)
        i = int(hit[0])
        return dict(
            t=float(t[i]), z_cl=float(self.cl.z[i]),
            z_wp=float(np.interp(t[i], self.wp.t, self.wp.s3)),
            v_cl=float(self.cl.vz[i]), v_wp=float(v_wp[i]), frac=frac,
        )

    def epp_on_z(self, half: str) -> tuple[np.ndarray, np.ndarray]:
        """(z, E_PP in eV) for one half, with E_PP mapped onto the trajectory.

        interactions.csv and the trajectory files are written on different
        cadences, so E_PP is interpolated onto the trajectory's own time grid
        rather than index-matched.
        """
        ix = self.ix_wp if half == "wp" else self.ix_cl
        if ix is None:
            return np.array([]), np.array([])
        if half == "wp":
            z, t = self.wp.s3, self.wp.t
        else:
            z, t = self.cl.z, self.cl.t
        return z, np.interp(t, ix.t, ix.e_pp) * HA_TO_EV


def load_pair(scripts_dir: str | Path, family: str, box_length_z: float,
              z0: float) -> PairPhase:
    """Load both halves of one twin pair, with interactions where available."""
    scripts_dir = Path(scripts_dir)
    wp = load_wp_run(scripts_dir / family / "wp", box_length_z, z0)
    cl = load_classical_run(scripts_dir / family / "classical", box_length_z)

    def _try(half: str) -> Interactions | None:
        try:
            return load_interactions(scripts_dir / family / half, half)
        except (FileNotFoundError, KeyError):
            return None

    return PairPhase(family=family, wp=wp, cl=cl,
                     ix_wp=_try("wp"), ix_cl=_try("classical"))
