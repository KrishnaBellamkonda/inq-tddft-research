"""Known-case tests for the in-slab-path stopping estimator (s5).

Plan: docs/plans/slab-ks-orbital-stopping-wrap.md

WHY THESE EXIST. s5 = integral f(t) <p_z>/m dt is the one piece of NEW analysis
logic in this study; every stopping power reported by Window B is -dT/ds5. If s5
is wrong, every number is wrong by the same factor and nothing downstream would
notice, because a wrong-by-a-constant stopping power still looks perfectly
plausible. So it is pinned here against synthetic runs whose answer is known by
construction.

The runs below are built by hand rather than loaded from disk: a real run cannot
tell you what its own stopping power "should" be, so testing against one would be
circular.

    python3 -m pytest ResearchProject/systems/localised_jellium/hypotheses/slab_ks_wrap/tests/
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import slab_ks_stopping as S  # noqa: E402


def _synthetic(v0: float, t_end: float, f_of_t, S_true_ev_per_bohr: float,
               n_elec: int = 100, dt: float = S.DT):
    """A run at constant velocity whose energy falls at exactly S_true per Bohr
    of IN-SLAB path, with a prescribed occupancy f(t).

        dT/dt = -(S_true/HA_TO_EV) * f(t) * v      [Ha per a.u.]

    so -dT/ds5 == S_true by construction, whatever f does — which is the property
    under test. -dT/ds4 is then S_true * <f>, i.e. WRONG unless f == 1.
    """
    t = np.arange(0.0, t_end + 0.5 * dt, dt)
    f = np.asarray(f_of_t(t), dtype=float)
    pz = np.full_like(t, v0)

    S_ha = S_true_ev_per_bohr / S.HA_TO_EV
    dTdt = -S_ha * f * v0
    T = np.concatenate([[0.5 * v0**2],
                        0.5 * v0**2 + np.cumsum(0.5 * (dTdt[1:] + dTdt[:-1]) * np.diff(t))])

    s4 = S.LAUNCH_Z + v0 * t
    integ = f * pz
    s5 = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(t))])

    base = S.K.WPRun(
        run_dir=Path("/synthetic"), box_length_z=S.LZ, t=t,
        step=np.arange(len(t)), T1=T, T2=T, pz=pz,
        s3=s4.copy(), s3_naive=s4.copy(), s4=s4,
        norm=np.ones_like(t), sigma_z=S.sigma_d(t), parseval=np.full_like(t, np.nan),
    )
    return S.SlabWPRun(
        name="synthetic", n_elec=n_elec, v0=v0, base=base, f_in_slab=f, s5=s5,
        e_total=np.zeros_like(t), complete=True,
        steps_done=len(t) - 1, steps_target=len(t) - 1,
    )


# ---------------------------------------------------------------------------

def test_fully_inside_slab_s5_equals_the_ordinary_path():
    """f == 1: the packet never leaves the medium, so the in-slab path IS the
    path and s5 must reduce to the ordinary trajectory."""
    run = _synthetic(v0=2.0, t_end=20.0, f_of_t=lambda t: np.ones_like(t),
                     S_true_ev_per_bohr=1.0)
    # s5 is measured from 0, s4 from the launch point: compare the displacements.
    assert np.allclose(run.s5, run.base.s4 - S.LAUNCH_Z, atol=1e-9)


def test_delocalised_packet_s5_recovers_the_true_stopping_power():
    """THE DECISIVE CASE. A packet spread uniformly over the box sits in the slab
    only 25/85 of the time, so a centroid-path fit under-reports the force by
    exactly that filling factor. s5 must undo it."""
    ff = S.FILLING_FACTOR
    S_true = 0.8
    run = _synthetic(v0=3.0, t_end=120.0,
                     f_of_t=lambda t: np.full_like(t, ff),
                     S_true_ev_per_bohr=S_true)

    fits = S.fit_wp(run)
    # Window B, in-slab path: the true force.
    assert fits["S25_B"].S_ev_per_bohr == pytest.approx(S_true, rel=1e-6)
    assert fits["S25_B"].r2 == pytest.approx(1.0, abs=1e-9)

    # Window A against the centroid path: diluted by the filling factor. This is
    # the error the estimator exists to prevent, so it is asserted explicitly.
    assert fits["S24"].S_ev_per_bohr == pytest.approx(S_true * ff, rel=1e-6)
    assert fits["S24"].S_ev_per_bohr < 0.35 * S_true


def test_top_hat_occupancy_recovers_S_while_centroid_path_does_not():
    """A localised packet crossing a slab: f is a top hat in time. -dT/ds5 must
    still give the true force, while -dT/ds4 averages in the vacuum flight."""
    v0, S_true = 2.0, 1.2

    def f_of_t(t):
        z = S.LAUNCH_Z + v0 * t                       # centroid, first pass only
        return (np.abs(z) <= S.SLAB_HALF).astype(float)

    run = _synthetic(v0=v0, t_end=(S.LZ / 2.0 - S.LAUNCH_Z) / v0, f_of_t=f_of_t,
                     S_true_ev_per_bohr=S_true)
    fits = S.fit_wp(run)

    assert fits["S25"].S_ev_per_bohr == pytest.approx(S_true, rel=1e-3)
    # The in-slab path over one full crossing is the slab thickness.
    assert run.s5[-1] == pytest.approx(S.SLAB_THICKNESS, rel=2e-3)
    # The centroid-path fit is diluted and must NOT be quoted as S.
    assert fits["S24"].S_ev_per_bohr < 0.75 * S_true


def test_s5_is_monotonic_and_never_exceeds_the_total_path():
    """f is a fraction in [0,1], so the in-slab path can only grow, and can never
    outrun the distance actually travelled."""
    rng = np.random.default_rng(0)
    run = _synthetic(v0=2.5, t_end=100.0,
                     f_of_t=lambda t: 0.5 * (1.0 + np.sin(t / 3.0)) * rng.uniform(0.5, 1.0, t.size),
                     S_true_ev_per_bohr=0.5)
    assert np.all(np.diff(run.s5) >= -1e-12)
    total_path = run.v0 * (run.t - run.t[0])
    assert np.all(run.s5 <= total_path + 1e-9)


def test_geometry_constants_match_the_configuration():
    """Guards against a silent edit to the geometry block: these are the numbers
    the plan, the run.cpp defaults and the SLURM tables all assume."""
    assert S.FILLING_FACTOR == pytest.approx(25.0 / 85.0)
    assert S.rs_for(100) == pytest.approx(4.181, abs=1e-3)
    assert S.rs_for(40) == pytest.approx(5.675, abs=1e-3)
    assert S.plasma_period_for(100) == pytest.approx(31.02, abs=0.02)
    assert S.plasma_period_for(40) == pytest.approx(49.04, abs=0.02)
    # sigma_d(t) = sqrt(2 + t^2/8) reaches the slab half-thickness at 35.1 a.u.
    assert S.sigma_d(S.T_SPREAD_LIMIT) == pytest.approx(S.SLAB_HALF, rel=1e-9)
    # ... and the transverse images overlap (6 sigma_d = L_xy) at 16.0 a.u.
    assert 6.0 * S.sigma_d(S.T_TRANSVERSE_OVERLAP) == pytest.approx(S.LX, rel=1e-9)
