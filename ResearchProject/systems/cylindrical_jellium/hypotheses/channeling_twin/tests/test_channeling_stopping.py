#!/usr/bin/env python3
"""Known-case tests for channeling_stopping.py — the analysis engine of the
annular-tube channeling twin (docs/plans/cylindrical-channeling-ks-stopping.md).

WHY THESE ARE EXACT AND NOT REGRESSION BASELINES
------------------------------------------------
A constant stopping power S has a closed-form trajectory, so the synthetic runs
below are not "some plausible numbers" — they are the analytic solution, and the
engine must return S back to machine precision.

    dT/ds = -S,   T = p^2 / 2m,   ds = (p/m) dt,   m = 1
    =>  p dp/dt = -S p  =>  dp/dt = -S        (CONSTANT deceleration)
    =>  p(t) = p0 - S t
        s(t) = s0 + p0 t - S t^2 / 2
        T(t) = p(t)^2 / 2 = T0 - S (s - s0)

p(t) is LINEAR, so the trapezoid rule the engine uses to rebuild s4 = integral
<p_z> dt is exact on it, and the OLS fit of T against s recovers the input S with
no discretisation error at all. Any deviation is a bug in the engine, not noise.

WHAT IS COVERED
  * classical S and the four WP definitions all recover the input S
  * s3 (circular centroid, unwrapped) reproduces s4 — this exercises the periodic
    unwrap on a trajectory that DOES cross the cell face, which the production run
    can do on a resume
  * the channeling window is derived from the MEASURED f_bore, and shortens when
    the packet leaves the bore
  * var_p_freeze distinguishes a frozen var(p) from a growing one
  * compare() returns AIM MET only when all three conditions hold, and returns the
    specific "clean channeling but S still differs" diagnosis when they do not
  * resume segments (observables.fromNNN.csv) are concatenated, not duplicated

Run:  <repo>/venv/bin/python3 -m pytest <this file> -q
      or standalone:  <repo>/venv/bin/python3 <this file>
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import channeling_stopping as CS  # noqa: E402

HA = CS.HA_TO_EV
S_TRUE_EV = 0.20                      # eV/Bohr — the input stopping power
S_TRUE_HA = S_TRUE_EV / HA            # 0.0073498 Ha/Bohr
P0 = CS.V0                            # 1.917 a.u.
Z0 = CS.LAUNCH_Z                      # -28
N = CS.N_STEPS                        # 1500 steps -> 1501 rows
DT = CS.DT


# ---------------------------------------------------------------------------
# Synthetic run construction
# ---------------------------------------------------------------------------

def analytic_trajectory(n_steps: int = N):
    """(step, t, p, s, T) for constant stopping power S_TRUE."""
    step = np.arange(n_steps + 1)
    t = step * DT
    p = P0 - S_TRUE_HA * t
    s = Z0 + P0 * t - 0.5 * S_TRUE_HA * t**2
    T = 0.5 * p**2
    return step, t, p, s, T


def wrap_into_cell(z: np.ndarray) -> np.ndarray:
    """Fold into [-L/2, L/2), the node convention of the circular centroid."""
    return (np.asarray(z) + CS.LZ / 2.0) % CS.LZ - CS.LZ / 2.0


def write_wp_run(root: Path, name: str = "wp", *,
                 f_bore_breach_t: float | None = None,
                 var_growth_pct: float = 0.0,
                 n_steps: int = N,
                 split_at: int | None = None) -> Path:
    """A synthetic wavepacket run directory.

    f_bore_breach_t : time at which f_bore falls below F_BORE_MIN (None = never)
    var_growth_pct  : linear growth of var(p_z) across the run, per cent
    split_at        : if given, write observables as two resume segments so the
                      segment-concatenation path is exercised
    """
    step, t, p, s, T2 = analytic_trajectory(n_steps)
    var = CS.VAR_P_FREE * (1.0 + var_growth_pct / 100.0 * t / t[-1])
    # T1 = <p^2>/2 = (<p_z>^2 + var_x + var_y + var_z)/2 ; the transverse variances
    # are held at the free value so T1 - T2 tracks var(p_z) alone.
    T1 = T2 + 0.5 * (2.0 * CS.VAR_P_FREE + var)

    obs = root / name / "raw" / "observables"
    obs.mkdir(parents=True, exist_ok=True)

    pd.DataFrame({
        "step": step, "time_au": t,
        "px_mean": 0.0, "py_mean": 0.0, "pz_mean": p,
        "px2_mean": 0.0, "py2_mean": 0.0, "pz2_mean": p**2 + var,
        "sigma_px2": CS.VAR_P_FREE, "sigma_py2": CS.VAR_P_FREE, "sigma_pz2": var,
        "e_kin_ha": T1, "norm_check": 4.9e7,
    }).to_csv(obs / "wp_momentum_stats.csv", index=False)

    pd.DataFrame({
        "step": step, "time_au": t,
        "x_mean": 0.0, "y_mean": 0.0, "z_mean": wrap_into_cell(s),
        "x2_mean": 0.0, "y2_mean": 0.0, "z2_mean": 0.0,
        "sigma_x2": CS.SIGMA_POT**2, "sigma_y2": CS.SIGMA_POT**2, "sigma_z2": CS.SIGMA_POT**2,
        "norm_check": 1.0,
        "x_mean_circ": 0.0, "y_mean_circ": 0.0, "z_mean_circ": wrap_into_cell(s),
        "R_x": 1.0, "R_y": 1.0, "R_z": 0.9,
        "sigma_x_circ": CS.SIGMA_POT, "sigma_y_circ": CS.SIGMA_POT,
        "sigma_z_circ": CS.sigma_d(t),
    }).to_csv(obs / "wp_real_space_stats.csv", index=False)

    f_bore = np.full_like(t, 0.998)
    if f_bore_breach_t is not None:
        # Linear decay to 0.80 from the breach time onwards.
        m = t >= f_bore_breach_t
        if m.any():
            frac = (t[m] - f_bore_breach_t) / max(t[-1] - f_bore_breach_t, 1e-12)
            f_bore[m] = 0.949 - 0.15 * frac      # first point already below 0.95
    pd.DataFrame({
        "step": step, "time_au": t,
        "f_bore": f_bore, "f_wall": 1.0 - f_bore, "f_outside": 0.0,
        "r_mean": CS.SIGMA_POT * math.sqrt(math.pi / 2.0) * np.ones_like(t),
        "r2_mean": 2.0 * CS.SIGMA_POT**2 * np.ones_like(t),
        "sigma_r": CS.SIGMA_POT * math.sqrt(2.0 - math.pi / 2.0) * np.ones_like(t),
        "norm_total": 1.0,
    }).to_csv(obs / "wp_radial_occupancy.csv", index=False)

    energies = pd.DataFrame({
        "step": step, "time_au": t, "energy_total": -45.0,
        "energy_kinetic": 10.0, "energy_hartree": 5.0, "energy_xc": -3.0,
        "energy_external": -20.0, "energy_nonlocal": 0.0, "energy_ion": 0.0,
    })
    if split_at is None:
        energies.to_csv(obs / "observables.csv", index=False)
    else:
        energies.iloc[: split_at + 1].to_csv(obs / "observables.csv", index=False)
        # The boundary step appears in BOTH files, as a real resume does.
        energies.iloc[split_at:].to_csv(obs / f"observables.from{split_at}.csv", index=False)

    pd.DataFrame({
        "step": step, "time_au": t,
        "e_ss": 1.0, "e_pp": 0.09, "e_ps": np.linspace(0.5, 0.2, len(t)),
        "e_sb": -2.0, "e_pb": -0.3, "e_bb": 0.7,
        "e_hartree_check": 1.19, "e_external_check": -2.3,
        "norm_proj": 1.0, "norm_electrons": 161.0,
    }).to_csv(obs / "interactions.csv", index=False)

    (root / name / "run_summary.txt").write_text(
        f"run_completed = true\nn_steps = {n_steps}\nrepresentation = wavepacket\n"
        f"projectile = wavepacket_orbital\n")
    return root / name


def write_classical_run(root: Path, name: str = "classical", *,
                        s_scale: float = 1.0, n_steps: int = N) -> Path:
    """A synthetic classical run whose stopping power is s_scale x S_TRUE."""
    step = np.arange(n_steps + 1)
    t = step * DT
    S = S_TRUE_HA * s_scale
    v = P0 - S * t
    z = Z0 + P0 * t - 0.5 * S * t**2
    T = 0.5 * v**2

    obs = root / name / "raw" / "observables"
    obs.mkdir(parents=True, exist_ok=True)

    pd.DataFrame({
        "step": step, "time_au": t,
        "x": 0.0, "y": 0.0, "z": wrap_into_cell(z),
        "vx": 0.0, "vy": 0.0, "vz": v, "ke_ion_ha": T,
    }).to_csv(obs / "electron_track.csv", index=False)

    pd.DataFrame({
        "step": step, "time_au": t,
        "proj_z": wrap_into_cell(z), "proj_vz": v,
        "energy_proj_ke": T, "energy_proj_bg_ideal": 0.0,
        "proj_z_unwrapped": z, "n_wraps": 0,
        "proj_x": 0.0, "proj_y": 0.0,
        "force_x": 0.0, "force_y": 0.0, "force_z": -S,
    }).to_csv(obs / "projectile.csv", index=False)

    # energy_total chosen so E_electronic + KE_proj + U_proj_bg is exactly flat.
    pd.DataFrame({
        "step": step, "time_au": t, "energy_total": -45.0 - T,
        "energy_kinetic": 10.0, "energy_hartree": 5.0, "energy_xc": -3.0,
        "energy_external": -20.0, "energy_nonlocal": 0.0, "energy_ion": 0.0,
    }).to_csv(obs / "observables.csv", index=False)

    pd.DataFrame({
        "step": step, "time_au": t,
        "e_ss": 1.0, "e_pp": 0.0, "e_ps": np.linspace(0.5, 0.2, len(t)),
        "e_sb": -2.0, "e_pb": -0.3, "e_bb": 0.7,
        "e_hartree_check": 1.0, "e_external_check": -2.3,
        "norm_proj": 1.0, "norm_electrons": 160.0,
    }).to_csv(obs / "interactions.csv", index=False)

    (root / name / "run_summary.txt").write_text(
        f"run_completed = true\nn_steps = {n_steps}\nrepresentation = perturbation\n"
        f"projectile = gaussian_charge_perturbation\n")
    return root / name


def _point(tmp_path, **wp_kw):
    """Build a twin pair under tmp_path and point the module at it."""
    CS.WP_RESULTS = tmp_path
    CS.CL_RESULTS = tmp_path
    write_wp_run(tmp_path, **wp_kw)
    return tmp_path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_classical_fit_recovers_the_input_stopping_power(tmp_path):
    CS.WP_RESULTS = CS.CL_RESULTS = tmp_path
    write_classical_run(tmp_path)
    cl = CS.load_classical()

    fit = CS.K.fit_stopping(cl.base.z, cl.base.T, cl.base.t, 5.0, 25.0, "test", v=cl.base.vz)
    assert fit.S_ev_per_bohr == pytest_approx(S_TRUE_EV, rel=1e-9)
    assert fit.r2 == pytest_approx(1.0, abs=1e-12)
    # The trajectory is exactly on-axis, which is what the tube's symmetry demands.
    assert cl.off_axis_max == 0.0
    # E_electronic + KE_proj + U_proj_bg was constructed flat.
    assert abs(cl.conserved_drift_ev) < 1e-9


def test_wp_all_four_definitions_recover_the_input(tmp_path):
    _point(tmp_path)
    wp = CS.load_wp()

    # s3 (circular centroid, unwrapped) must reproduce s4 (integral <p> dt). The
    # synthetic trajectory crosses the +z face, so this is a real unwrap test.
    assert np.max(np.abs(wp.base.s3 - wp.base.s4)) < 1e-8, "s3 and s4 disagree"

    fits = CS.K.fit_all_wp(wp.base, 5.0, 25.0)
    for key, f in fits.items():
        # T1 differs from T2 by a CONSTANT here (var frozen), so all four slopes
        # are identical and equal to the input.
        assert f.S_ev_per_bohr == pytest_approx(S_TRUE_EV, rel=1e-8), key


def test_channeling_window_uses_the_measured_f_bore(tmp_path):
    # Clean run: f_bore never breaches, so the window runs to the end.
    _point(tmp_path)
    wp = CS.load_wp()
    t0, t1 = wp.channeling_window()
    assert t1 == pytest_approx(wp.t[-1])
    assert t0 == pytest_approx(CS.TRANSIENT_FRAC * t1)
    assert CS.channeling_check(wp).channeling is True

    # Breached run: the window must END AT THE BREACH, not at the end of the run,
    # and not at the free-dispersion estimate either.
    _point(tmp_path, f_bore_breach_t=10.0)
    wp2 = CS.load_wp()
    _, t1b = wp2.channeling_window()
    assert t1b == pytest_approx(10.0, abs=DT)
    assert t1b < CS.T_2SIGMA_AT_WALL, "the measured window must be able to beat the formula"
    ch = CS.channeling_check(wp2)
    assert ch.t_breach == pytest_approx(10.0, abs=DT)
    assert ch.channeling is False        # breached in the first third of the run


def test_var_p_freeze_separates_frozen_from_growing(tmp_path):
    _point(tmp_path)
    frozen = CS.var_p_freeze(CS.load_wp())
    assert frozen.frozen is True
    assert abs(frozen.growth_pct) < 1e-9
    # (T1 - T2) is the localisation energy and must sit at its closed-form value.
    assert frozen.localisation_start_ev == pytest_approx(CS.LOCALISATION_EV, rel=1e-9)
    assert abs(frozen.localisation_drift_ev) < 1e-9

    _point(tmp_path, var_growth_pct=50.0)
    growing = CS.var_p_freeze(CS.load_wp())
    assert growing.frozen is False
    # The window is a sub-interval of the run, so the growth seen over it is a
    # fraction of the 50 % imposed across the whole run — but must be well outside
    # the 10 % tolerance and correctly signed.
    assert growing.growth_pct > 10.0
    assert growing.localisation_drift_ev > 0.0


def test_compare_reports_aim_met_only_when_all_three_hold(tmp_path):
    CS.WP_RESULTS = CS.CL_RESULTS = tmp_path
    write_classical_run(tmp_path)

    # (a) matched S, frozen var, clean channel -> AIM MET
    write_wp_run(tmp_path)
    c = CS.compare(CS.load_wp(), CS.load_classical())
    assert c.aim_met is True
    assert c.agreement_pct < 1e-6
    assert "AIM MET" in c.verdict
    assert c.wp_fits[CS.PRIMARY_ESTIMATOR].S_ev_per_bohr == pytest_approx(S_TRUE_EV, rel=1e-8)
    assert c.cl_same_window.S_ev_per_bohr == pytest_approx(S_TRUE_EV, rel=1e-8)
    assert set(c.table()["estimator"]) >= {"S_13", "S_14", "S_23", "S_24",
                                           "S_cl_same_window", "S_cl_initial_drag"}

    # (b) clean channel + frozen var but a classical twin stopping 2.5x harder —
    #     the "interesting failure" branch, which must NOT be reported as met.
    write_classical_run(tmp_path, name="classical_hard", s_scale=2.5)
    c2 = CS.compare(CS.load_wp(), CS.load_classical("classical_hard"))
    assert c2.aim_met is False
    assert c2.agreement_pct > 20.0
    assert "AIM NOT MET" in c2.verdict
    assert "self-Hartree" in c2.verdict     # points at E_PP, the right next step

    # (c) matched S but the packet left the bore -> agreement is unexplained
    write_wp_run(tmp_path, name="wp_breach", f_bore_breach_t=6.0)
    c3 = CS.compare(CS.load_wp("wp_breach"), CS.load_classical())
    assert c3.aim_met is False
    assert "AIM PARTLY MET" in c3.verdict


def test_resume_segments_are_concatenated_without_duplication(tmp_path):
    _point(tmp_path, split_at=900)
    wp = CS.load_wp()
    # 1501 unique steps across the two segments, with the boundary step counted once.
    assert wp.e_total.size == N + 1
    assert wp.steps_done == N
    assert wp.complete is True


def test_load_interactions_adds_gauge_safe_deltas(tmp_path):
    _point(tmp_path)
    df = CS.load_interactions(CS.WP_RESULTS / "wp")
    for c in ("e_pp_ev", "d_e_ps_ev", "d_e_pp_ev"):
        assert c in df.columns
    assert df["e_pp_ev"].iloc[0] == pytest_approx(0.09 * HA)
    # Deltas start at zero by construction — that is what makes them gauge-safe.
    assert df["d_e_ps_ev"].iloc[0] == 0.0
    assert df["d_e_ps_ev"].iloc[-1] == pytest_approx((0.2 - 0.5) * HA)


def test_geometry_constants_match_the_locked_design():
    """The module's mirrored constants must still describe r_s = 3."""
    assert CS.RS == pytest_approx(3.0, abs=1e-5)
    assert CS.N0 == pytest_approx(160.0 / (math.pi * (14.0**2 - 10.0**2) * 60.0))
    assert CS.OMEGA_P * HA == pytest_approx(9.07, abs=0.01)
    assert CS.V0 / CS.V_FERMI == pytest_approx(3.0, abs=0.01)
    assert CS.SIGMA_POT == pytest_approx(4.0 / math.sqrt(2.0))
    assert CS.LOCALISATION_EV == pytest_approx(1.2755, abs=1e-3)
    # 2 sigma_d reaches the bore wall inside the run; images overlap after it ends.
    assert 0.0 < CS.T_2SIGMA_AT_WALL < N * DT
    assert CS.T_TRANSVERSE_OVERLAP > N * DT


# ---------------------------------------------------------------------------
# pytest-free fallback so the file runs standalone on a node without pytest
# ---------------------------------------------------------------------------

try:
    from pytest import approx as pytest_approx
except ImportError:                                        # pragma: no cover
    class _Approx:
        def __init__(self, v, rel=None, abs=None):
            self.v, self.rel, self.abs = v, rel, abs
        def __eq__(self, other):
            tol = self.abs if self.abs is not None else 0.0
            if self.rel is not None:
                tol = max(tol, self.rel * builtins_abs(self.v))
            if tol == 0.0:
                tol = 1e-9
            return builtins_abs(other - self.v) <= tol
    builtins_abs = abs
    def pytest_approx(v, rel=None, abs=None):              # noqa: A002
        return _Approx(v, rel=rel, abs=abs)


if __name__ == "__main__":
    import tempfile, traceback
    fails = 0
    tests = [(n, o) for n, o in sorted(globals().items())
             if n.startswith("test_") and callable(o)]
    for name, fn in tests:
        with tempfile.TemporaryDirectory() as td:
            try:
                if fn.__code__.co_argcount:
                    fn(Path(td))
                else:
                    fn()
                print(f"  PASS  {name}")
            except Exception:
                fails += 1
                print(f"  FAIL  {name}")
                traceback.print_exc()
    print(f"\n{len(tests) - fails}/{len(tests)} passed")
    raise SystemExit(1 if fails else 0)
