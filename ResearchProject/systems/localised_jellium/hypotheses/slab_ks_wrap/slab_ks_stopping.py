"""Stopping-power extraction for the CAP-free wrap-around slab KS study.

Plan: docs/plans/slab-ks-orbital-stopping-wrap.md

Deterministic engine: it does the arithmetic and nothing else, so every notebook
and every cross-run table gets bit-identical numbers from one implementation.
The four KS-orbital definitions themselves are NOT re-implemented here — they
live in ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/ks_stopping.py
(``ks_stopping``) and are imported, which is the whole point of the study: the
BULK definitions are being applied unchanged to a SLAB.

WHAT THIS MODULE ADDS OVER THE BULK ENGINE
------------------------------------------

1. LAYOUT. These runs write to
   ``<scripts>/slab_ks_wrap/{wp,classical}/results/<name>/raw/observables/``.

2. s5, THE IN-SLAB PATH — the reason the study exists.

   Stopping power is a force: energy lost per unit path INSIDE the medium. The
   slab is 25 of the 85 Bohr the packet traverses, and at sigma_WP = 2 the packet
   is wider than the slab from t ~ 35 a.u. onward
   (sigma_d(t) = sqrt(sigma^2/2 + t^2/(2 sigma^2)) = sqrt(2 + t^2/8)), so it is
   inside and outside at the same time. Fitting -dT/ds against the centroid path
   then averages the drag over slab AND vacuum and under-reports it.

       f(t)        = in-slab fraction of |psi|^2, MEASURED on the grid every step
                     (wp_slab_occupancy.csv, from inqkit::observables::slab_occupancy)
       s5(t)       = integral of f(t') <p_z>(t')/m dt'

   With dT/dt = -F v f and ds5/dt = f v, we get -dT/ds5 = F exactly, in the
   localised limit (f -> 1, where s5 reduces to the ordinary path) AND in the
   delocalised one (f -> 25/85 = 0.294, where it applies the filling factor).

3. TWO WINDOWS, reported side by side (user decision 2026-07-31).

   WINDOW A — first pass. t in [transient, min(t_wrap, t_spread)] with
   t_wrap = (L_z/2 - z0)/v the first +z face crossing and t_spread = 34.64 a.u.
   the moment sigma_d reaches the slab half-thickness. This is the localised
   regime and is what the classical single-pass benchmark is comparable with.

   WINDOW B — whole run, fitted against s5. This is the long, defensible window
   the study was set up to obtain; it is only meaningful against s5, never
   against s3/s4.

4. THE WRAPPED CLASSICAL TWIN. projectile.csv here carries proj_z (wrapped) and
   proj_z_unwrapped (continuous), so no unwrapping heuristic is needed; the
   module cross-checks them against each other.

WHAT IS DELIBERATELY NOT DONE
-----------------------------
No CAP correction: there is no CAP. That also means ``energy_total`` is a
conserved quantity rather than a fit target, so the deposit-based estimator
(Method A of the stopping-power-extraction skill) is UNDEFINED for the WP half —
its fit target is identically zero because the wavepacket is inside the system.
The four -dT/ds slopes are the measurement; energy conservation is the gate.

CARRIED-OVER CONCLUSION from docs/handovers/bulk-jellium-ks-stopping.md, to be
re-tested and not re-derived: T1 - T2 grows at a fixed ~0.043 eV/Bohr that does
NOT scale with bath density, most plausibly self-interaction error (the packet is
an occupied KS orbital whose own charge enters the Hartree potential). So S2 is
the defensible stopping power and S1 must be reported as "S2 minus a spreading
term", never as an independent measurement.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
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
# Geometry. Mirrored from shared/configs/slab_n{100,40}_L35x35x85.hpp;
# run_summary.txt is the authority if these ever disagree.
# ---------------------------------------------------------------------------
LX = LY = 35.0
LZ = 85.0
SLAB_HALF = 12.5
SLAB_THICKNESS = 2.0 * SLAB_HALF
FILLING_FACTOR = SLAB_THICKNESS / LZ          # 0.2941
DX = 0.40
DT = 0.04
SIGMA_WP = 2.0
LAUNCH_Z = -24.0
V_SLAB = LX * LY * SLAB_THICKNESS              # 30625 Bohr^3

SCRIPTS = REPO / "ResearchProject/systems/localised_jellium/scripts/slab_ks_wrap"
WP_RESULTS = SCRIPTS / "wp" / "results"
CL_RESULTS = SCRIPTS / "classical" / "results"

DENSITIES = (100, 40)
VELOCITIES = (2.0, 2.5, 3.0, 3.5)
N_STEPS = {2.0: 4529, 2.5: 3623, 3.0: 3020, 3.5: 2588}

# Moment at which the packet's density std reaches the slab HALF-thickness, i.e.
# when it stops being small compared with the medium. sigma_d(t) = sqrt(2+t^2/8)
# = 12.5 at t = sqrt(8*(12.5^2 - 2)) = 34.64 a.u. Velocity-independent: dispersion
# is set by sigma, not by v.
T_SPREAD_LIMIT = float(np.sqrt(8.0 * (SLAB_HALF**2 - SIGMA_WP**2 / 2.0)))

# Transverse periodic images overlap when 6 sigma_d = L_xy: t = 16.0 a.u. Past
# this the object being dragged is a periodic ARRAY of packets. Reported in every
# notebook; not correctable without changing L_xy, which would change r_s.
T_TRANSVERSE_OVERLAP = float(np.sqrt(8.0 * ((LX / 6.0) ** 2 - SIGMA_WP**2 / 2.0)))


def n0_for(n_elec: int) -> float:
    return n_elec / V_SLAB


def rs_for(n_elec: int) -> float:
    return float((3.0 / (4.0 * np.pi * n0_for(n_elec))) ** (1.0 / 3.0))


def plasma_period_for(n_elec: int) -> float:
    """2 pi / omega_p in a.u. (omega_p = sqrt(4 pi n) for a HEG)."""
    return float(2.0 * np.pi / np.sqrt(4.0 * np.pi * n0_for(n_elec)))


def sigma_d(t: np.ndarray | float) -> np.ndarray | float:
    """Free-dispersion density std of the launched Gaussian, Bohr."""
    return np.sqrt(SIGMA_WP**2 / 2.0 + np.asarray(t) ** 2 / (2.0 * SIGMA_WP**2))


def run_name(n_elec: int, v: float) -> str:
    return f"n{n_elec}_v{str(v).replace('.', 'p')}"


# ---------------------------------------------------------------------------
# Wavepacket runs
# ---------------------------------------------------------------------------

@dataclass
class SlabWPRun:
    """A wavepacket run: the bulk WPRun plus the in-slab occupancy channel."""
    name: str
    n_elec: int
    v0: float
    base: K.WPRun               # t, T1, T2, pz, s3, s4, norm, sigma_z
    f_in_slab: np.ndarray       # in-slab fraction of |psi|^2, per step
    s5: np.ndarray              # Bohr, in-slab path = integral f <p_z>/m dt
    e_total: np.ndarray         # Ha, from observables.csv (conserved: no CAP)
    complete: bool
    steps_done: int
    steps_target: int

    @property
    def t(self) -> np.ndarray:
        return self.base.t

    @property
    def rs(self) -> float:
        return rs_for(self.n_elec)

    @property
    def energy_drift_ev(self) -> float:
        """E_total(end) - E_total(0) in eV. Must be ~0: no CAP => H Hermitian."""
        if self.e_total.size < 2:
            return float("nan")
        return float((self.e_total[-1] - self.e_total[0]) * HA_TO_EV)

    @property
    def norm_drift(self) -> float:
        return float(abs(self.base.norm[-1] - self.base.norm[0]))

    def window_a(self, transient_frac: float = 0.2) -> tuple[float, float]:
        """First-pass window: launch transient dropped, ends at the earlier of the
        first face crossing and the moment the packet outgrows the slab."""
        t_wrap = (LZ / 2.0 - LAUNCH_Z) / self.v0
        t1 = min(t_wrap, T_SPREAD_LIMIT)
        return (transient_frac * t1, t1)

    def window_b(self, transient_frac: float = 0.05) -> tuple[float, float]:
        """Whole run, minus a short launch transient."""
        t_end = float(self.t[-1])
        return (transient_frac * t_end, t_end)


def _load_wp(name: str, n_elec: int, v0: float) -> SlabWPRun:
    run_dir = WP_RESULTS / name
    obs = run_dir / "raw" / "observables"
    if not obs.is_dir():
        raise FileNotFoundError(f"no observables under {run_dir}")

    # K.load_wp_run hard-codes an extra "results" path level that this layout does
    # not have, so the same steps are done inline here. The DEFINITIONS still come
    # from the bulk module (unwrap_periodic, WPRun, fit_stopping) — only the file
    # lookup differs.
    mom = K._concat_segments(obs, "wp_momentum_stats")
    pos = K._concat_segments(obs, "wp_real_space_stats")
    df = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))

    t = df["time_au"].to_numpy()
    T1 = df["e_kin_ha"].to_numpy()
    px, py, pz = (df[c].to_numpy() for c in ("px_mean", "py_mean", "pz_mean"))
    T2 = 0.5 * (px**2 + py**2 + pz**2)

    if "z_mean_circ" not in df.columns:
        raise KeyError("z_mean_circ missing — this run predates the circular centroid.")
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

    # --- in-slab occupancy and the in-slab path ---------------------------
    occ = K._concat_segments(obs, "wp_slab_occupancy")
    occ = occ.drop_duplicates(subset="step").set_index("step")
    f = occ["f_in_slab"].reindex(df["step"].to_numpy()).to_numpy()
    if np.isnan(f).any():
        # Written every step by the binary; a gap means a truncated segment.
        f = pd.Series(f).interpolate(limit_direction="both").to_numpy()
    # s5 = integral f v dt, v = <p_z>/m with m = 1. Trapezoid, same grid as s4.
    integrand = f * pz
    s5 = np.concatenate(
        [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(t))])

    # --- energy conservation ----------------------------------------------
    try:
        en = K._concat_segments(obs, "observables")
        e_total = en["energy_total"].to_numpy()
    except (FileNotFoundError, KeyError):
        e_total = np.array([])

    steps_done = int(df["step"].to_numpy()[-1])
    steps_target = N_STEPS.get(v0, steps_done)
    return SlabWPRun(
        name=name, n_elec=n_elec, v0=v0, base=base, f_in_slab=f, s5=s5,
        e_total=e_total, complete=steps_done >= steps_target - 1,
        steps_done=steps_done, steps_target=steps_target,
    )


def load_wp(n_elec: int, v: float) -> SlabWPRun:
    return _load_wp(run_name(n_elec, v), n_elec, v)


# ---------------------------------------------------------------------------
# Wrapped classical twins
# ---------------------------------------------------------------------------

@dataclass
class SlabClassicalRun:
    name: str
    n_elec: int
    v0: float
    t: np.ndarray
    step: np.ndarray
    z_wrapped: np.ndarray
    z: np.ndarray            # continuous path, straight from the binary
    vz: np.ndarray
    T: np.ndarray            # Ha, 1/2 m v^2
    n_wraps: np.ndarray

    @property
    def v_fraction(self) -> np.ndarray:
        return self.vz / self.vz[0]

    @property
    def stopped_at_bohr(self) -> float:
        """Path length at which |v| first falls below 10 % of v0 — a mass-1
        electron at v = 2 is EXPECTED to stop after ~2 slab crossings."""
        below = np.where(np.abs(self.v_fraction) < 0.1)[0]
        return float(self.z[below[0]] - self.z[0]) if below.size else float("nan")

    def in_slab_mask(self) -> np.ndarray:
        return np.abs(self.z_wrapped) <= SLAB_HALF

    @property
    def s_in_slab(self) -> np.ndarray:
        """The classical analogue of s5: path accumulated while inside the slab.
        For a point-like projectile the occupancy is 0 or 1, so this is just the
        trajectory clipped to the slab."""
        inside = self.in_slab_mask().astype(float)
        integrand = inside * self.vz
        return np.concatenate(
            [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(self.t))])


def load_classical(n_elec: int, v: float) -> SlabClassicalRun:
    name = run_name(n_elec, v)
    obs = CL_RESULTS / name / "raw" / "observables"
    df = K._concat_segments(obs, "projectile")
    if "proj_z_unwrapped" not in df.columns:
        raise KeyError(
            "proj_z_unwrapped missing from projectile.csv — this run predates the "
            "wrapped classical twin. The wrapped proj_z cannot be differenced."
        )
    return SlabClassicalRun(
        name=name, n_elec=n_elec, v0=v,
        t=df["time_au"].to_numpy(), step=df["step"].to_numpy(),
        z_wrapped=df["proj_z"].to_numpy(), z=df["proj_z_unwrapped"].to_numpy(),
        vz=df["proj_vz"].to_numpy(), T=df["energy_proj_ke"].to_numpy(),
        n_wraps=df["n_wraps"].to_numpy(),
    )


# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------

def fit_wp(run: SlabWPRun) -> dict[str, K.StoppingFit]:
    """Every (T, s) combination in both windows.

    Window A keys are S13/S14/S23/S24 (the bulk study's four), plus S15/S25
    against the in-slab path. Window B keys are S15_B/S25_B — the whole-run fit,
    which is ONLY meaningful against s5.
    """
    out: dict[str, K.StoppingFit] = {}
    tA0, tA1 = run.window_a()
    tB0, tB1 = run.window_b()

    for i, T in ((1, run.base.T1), (2, run.base.T2)):
        for j, s in ((3, run.base.s3), (4, run.base.s4), (5, run.s5)):
            out[f"S{i}{j}"] = K.fit_stopping(
                s, T, run.t, tA0, tA1, f"S{i}{j} (window A, first pass)",
                v=run.base.pz)
        out[f"S{i}5_B"] = K.fit_stopping(
            run.s5, T, run.t, tB0, tB1, f"S{i}5 (window B, whole run, in-slab path)",
            v=run.base.pz)
    return out


def fit_classical_windows(run: SlabClassicalRun,
                          transient_frac: float = 0.2) -> dict[str, K.StoppingFit]:
    """Initial drag (window A) and the whole-run in-slab fit (window B).

    A light free-Ehrenfest projectile decelerates by design
    (.claude/rules/light-projectile-stopping.md), so the whole-run number is an
    average over the swept velocity range, not S at v0. Window A is the S(v0)
    to quote; window B is the like-for-like comparator for the WP's window B.
    """
    out: dict[str, K.StoppingFit] = {}
    t_wrap = (LZ / 2.0 - LAUNCH_Z) / run.v0
    t1 = min(t_wrap, T_SPREAD_LIMIT)
    out["S_A"] = K.fit_stopping(run.z, run.T, run.t, transient_frac * t1, t1,
                                "S_classical (window A, first pass)", v=run.vz)
    # Window B stops where the projectile does: fitting a stopped projectile adds
    # points with zero path and zero energy loss and biases the slope to zero.
    moving = np.where(np.abs(run.v_fraction) >= 0.3)[0]
    tB1 = float(run.t[moving[-1]]) if moving.size else float(run.t[-1])
    out["S_B"] = K.fit_stopping(run.s_in_slab, run.T, run.t,
                                0.05 * tB1, tB1,
                                "S_classical (window B, in-slab path, while moving)",
                                v=run.vz)
    return out


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def summarise() -> pd.DataFrame:
    """One row per (density, velocity): both halves, both windows, plus gates."""
    rows = []
    for n_elec in DENSITIES:
        for v in VELOCITIES:
            row: dict[str, object] = {
                "n_elec": n_elec, "r_s": round(rs_for(n_elec), 3),
                "n0": n0_for(n_elec), "v": v,
                "T_plasmon_au": round(plasma_period_for(n_elec), 2),
                "n_plasma_periods": round(N_STEPS[v] * DT / plasma_period_for(n_elec), 2),
            }
            try:
                wp = load_wp(n_elec, v)
                fits = fit_wp(wp)
                row.update(
                    wp_complete=wp.complete,
                    wp_steps=f"{wp.steps_done}/{wp.steps_target}",
                    energy_drift_ev=wp.energy_drift_ev,
                    norm_drift=wp.norm_drift,
                    ehrenfest_resid_bohr=float(np.max(np.abs(wp.base.ehrenfest_residual))),
                    f_final=float(wp.f_in_slab[-1]),
                    f_geometric=FILLING_FACTOR,
                    **{k: f.S_ev_per_bohr for k, f in fits.items()},
                    **{f"{k}_r2": f.r2 for k, f in fits.items()},
                )
            except (FileNotFoundError, KeyError) as e:
                row["wp_error"] = str(e)[:80]
            try:
                cl = load_classical(n_elec, v)
                cfits = fit_classical_windows(cl)
                row.update(
                    cl_S_A=cfits["S_A"].S_ev_per_bohr,
                    cl_S_B=cfits["S_B"].S_ev_per_bohr,
                    cl_S_A_r2=cfits["S_A"].r2,
                    cl_v_final_frac=float(cl.v_fraction[-1]),
                    cl_n_wraps=int(cl.n_wraps[-1]),
                    cl_stopped_after_bohr=cl.stopped_at_bohr,
                )
            except (FileNotFoundError, KeyError) as e:
                row["cl_error"] = str(e)[:80]
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    df = summarise()
    print(df.to_string(index=False))
    out = HERE / "S_summary.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")
