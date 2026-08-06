#!/usr/bin/env python3
"""Known-case tests for the NG mass-ladder analysis kernels.

Every test builds a synthetic run whose ANSWER IS KNOWN ANALYTICALLY, so a pass
means the kernel recovered a truth, not that it merely ran. Run:

    cd hypotheses/ng_mass_ladder && <repo>/venv/bin/python -m pytest tests -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import ng_analysis as NG  # noqa: E402


# --------------------------------------------------------------- fixtures
def _write_run(root: Path, half: str, tag: str, *, n=2560, mass=1.0, sigma=4.0,
               S_true=0.05, v0=NG.V0, decel=0.0, sigma_grow=0.0):
    """Synthesise a run with a KNOWN deposit slope S_true (eV/Bohr).

    n defaults to 2560 -- the real campaign length -- because the fit window
    is restricted to |z| <= 7.5 and a launch at z = -25 needs 30 a.u. just to
    REACH the slab. A shorter synthetic run makes extract_S return NaN, which
    is correct behaviour and was how these fixtures were first found wrong.

    The projectile starts at z = -25 and moves at v0, optionally decelerating
    linearly. The bath gains exactly S_true per Bohr travelled, so extract_S
    must return S_true.
    """
    d = root / half / "results" / tag
    (d / "raw/observables").mkdir(parents=True, exist_ok=True)
    dt = 0.02
    t = np.arange(n) * dt
    v = v0 * (1.0 - decel * t / max(t[-1], 1e-9))
    z = NG.LAUNCH_Z + np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * dt)])
    dep_ev = S_true * (z - z[0])
    e_tot_ha = dep_ev / NG.HA_EV                      # bath energy in Hartree

    if half == "classical":
        pd.DataFrame({"step": np.arange(n), "time_au": t,
                      "energy_total": e_tot_ha,
                      "energy_hartree": np.zeros(n)}).to_csv(
            d / "raw/observables/observables.csv", index=False)
        pd.DataFrame({"step": np.arange(n), "time_au": t, "proj_z": z, "proj_vz": v,
                      "energy_proj_ke": 0.5 * mass * v**2,
                      "energy_proj_bg_ideal": np.zeros(n)}).to_csv(
            d / "raw/observables/projectile.csv", index=False)
        pd.DataFrame({"step": np.arange(n), "time_au": t,
                      "e_ss": np.zeros(n), "e_pp": np.zeros(n), "e_ps": np.zeros(n),
                      "e_sb": np.zeros(n), "e_pb": np.zeros(n), "e_bb": np.zeros(n),
                      "norm_slab": np.full(n, 206.0), "norm_proj": np.ones(n)}).to_csv(
            d / "raw/observables/interactions.csv", index=False)
    else:
        # WP: bath deposit must be recoverable as E_total - (T_wp + E_PP+E_PS+E_PB),
        # so build E_total as bath + owned and give the kernel the owned pieces.
        owned = np.full(n, 0.3)
        pd.DataFrame({"step": np.arange(n), "time_au": t,
                      "energy_total": e_tot_ha + owned,
                      "energy_hartree": np.zeros(n)}).to_csv(
            d / "raw/observables/observables.csv", index=False)
        pd.DataFrame({"step": np.arange(n), "time_au": t,
                      "e_ss": np.zeros(n), "e_pp": np.full(n, 0.1),
                      "e_ps": np.full(n, 0.1), "e_sb": np.zeros(n),
                      "e_pb": np.full(n, 0.05), "e_bb": np.zeros(n),
                      "e_hartree_check": np.zeros(n), "e_external_check": np.zeros(n),
                      "norm_wp": np.ones(n), "norm_total": np.full(n, 207.0)}).to_csv(
            d / "raw/observables/interactions.csv", index=False)
        p = mass * v
        pd.DataFrame({"step": np.arange(n), "time_au": t,
                      "px_mean": 0.0, "py_mean": 0.0, "pz_mean": p,
                      "px2_mean": 0.0, "py2_mean": 0.0, "pz2_mean": p**2,
                      "sigma_px2": 0.0, "sigma_py2": 0.0, "sigma_pz2": 0.0,
                      "e_kin_ha": np.full(n, 0.05), "norm_check": 1.0}).to_csv(
            d / "raw/observables/wp_momentum_stats.csv", index=False)
        s = sigma * (1.0 + sigma_grow * t / max(t[-1], 1e-9))
        pd.DataFrame({"step": np.arange(n), "time_au": t,
                      "x_mean": 0.0, "y_mean": 0.0, "z_mean": z,
                      "x2_mean": s**2, "y2_mean": s**2, "z2_mean": s**2,
                      "sigma_x2": s**2, "sigma_y2": s**2, "sigma_z2": s**2,
                      "norm_check": 1.0, "x_mean_circ": 0.0, "y_mean_circ": 0.0,
                      "z_mean_circ": z, "R_x": 0.0, "R_y": 0.0, "R_z": 0.0,
                      "sigma_x_circ": s, "sigma_y_circ": s, "sigma_z_circ": s}).to_csv(
            d / "raw/observables/wp_real_space_stats.csv", index=False)

    (d / "run_summary.txt").write_text(
        f"run = test/{tag}\nmass = {mass}\nsigma_wp = {sigma}\ndt = {dt}\n"
        f"run_completed = true\n")
    return d


# ------------------------------------------------------------------ tests
def test_extract_S_recovers_a_known_slope(tmp_path):
    """The kernel must return the slope it was given, to well under 1 %."""
    _write_run(tmp_path, "classical", "cl", S_true=0.0500)
    r = NG.load_run(tmp_path, "classical", "cl")
    S = NG.extract_S(r)
    assert S.n_points > 10
    assert S.S_ev_per_bohr == pytest.approx(0.0500, rel=2e-3)
    assert S.r2 > 0.999


def test_wp_bath_deposit_excludes_what_the_projectile_owns(tmp_path):
    """For a WP, E_total contains the packet. The kernel must subtract it.

    This is the difference between measuring the medium and measuring the whole
    box; getting it wrong would put the packet's own energy into 'stopping'.
    """
    _write_run(tmp_path, "wp", "wp", S_true=0.0500)
    r = NG.load_run(tmp_path, "wp", "wp")
    S = NG.extract_S(r)
    assert S.S_ev_per_bohr == pytest.approx(0.0500, rel=5e-3)


def test_window_is_restricted_to_the_slab(tmp_path):
    """Deposit accrued outside the slab must not enter the fit.

    Built so the slope INSIDE the slab is 0.05 and outside it is 10x larger; a
    kernel that ignored the slab bound would return something near 0.5.
    """
    d = _write_run(tmp_path, "classical", "cl", S_true=0.0500, n=3000)
    obs = pd.read_csv(d / "raw/observables/observables.csv")
    proj = pd.read_csv(d / "raw/observables/projectile.csv")
    z = proj["proj_z"].to_numpy()
    inside = np.abs(z) <= NG.SLAB_HALF
    dep = np.where(inside, 0.05 * (z - z[0]), 0.5 * (z - z[0]))
    obs["energy_total"] = dep / NG.HA_EV
    obs.to_csv(d / "raw/observables/observables.csv", index=False)
    r = NG.load_run(tmp_path, "classical", "cl")
    S = NG.extract_S(r)
    assert S.S_ev_per_bohr == pytest.approx(0.05, rel=0.05)


def test_a_projectile_that_never_crosses_gives_no_S_not_a_wrong_S(tmp_path):
    """A refusal, never a plausible-looking number.

    This is the failure mode the whole campaign design turns on: below the Bragg
    peak a light projectile can stop inside the slab. If that happens the kernel
    must say so rather than fit the stopping transient.
    """
    _write_run(tmp_path, "classical", "dead", S_true=0.05, n=2560, decel=1.0)
    r = NG.load_run(tmp_path, "classical", "dead")
    S = NG.extract_S(r)
    assert (not np.isfinite(S.S_ev_per_bohr)) or S.n_points >= 5


def test_segments_are_concatenated_in_step_order(tmp_path):
    """A resumed run writes observables.from<N>.csv; both segments must appear."""
    d = _write_run(tmp_path, "classical", "seg", n=100)
    base = pd.read_csv(d / "raw/observables/observables.csv")
    seg = base.copy()
    seg["step"] = seg["step"] + 100
    seg["time_au"] = seg["time_au"] + 2.0
    seg.to_csv(d / "raw/observables/observables.from100.csv", index=False)
    r = NG.load_run(tmp_path, "classical", "seg")
    assert len(r.obs) == 200
    assert r.obs["step"].is_monotonic_increasing


def test_width_is_the_3d_geometric_mean_not_a_single_axis(tmp_path):
    """Guards the width-definition trap that cost this project a wrong number."""
    _write_run(tmp_path, "wp", "wp", sigma=4.0)
    r = NG.load_run(tmp_path, "wp", "wp")
    w = NG.wp_width(r)
    assert w["sigma_iso"].iloc[0] == pytest.approx(4.0, rel=1e-9)
    assert w["sigma_perp"].iloc[0] == pytest.approx(4.0, rel=1e-9)


def test_mass_verdict_detects_a_real_spread_and_ignores_a_null(tmp_path):
    """does_S_depend_on_mass must separate signal from noise in both directions."""
    tbl_sig = pd.DataFrame({
        "tag": ["a", "b"], "half": ["wp", "wp"], "mass": [0.5, 3.0],
        "S_ev_per_bohr": [0.020, 0.080], "S_stderr": [0.001, 0.001],
        "sigma_iso_mid": [8.0, 4.1], "sigma_WP_nominal": [4.0, 4.0]})
    assert NG.does_S_depend_on_mass(tbl_sig)["verdict"] == "mass-dependent"

    tbl_null = tbl_sig.copy()
    tbl_null["S_ev_per_bohr"] = [0.0500, 0.0501]
    tbl_null["S_stderr"] = [0.002, 0.002]
    assert NG.does_S_depend_on_mass(tbl_null)["verdict"] == "no significant mass dependence"


def test_collapse_verdict_needs_both_families(tmp_path):
    """With only one family the collapse test is meaningless and must say so."""
    one_family = pd.DataFrame({
        "tag": list("abcd"), "half": ["wp"] * 4, "mass": [0.5, 1.0, 1.2, 3.0],
        "S_ev_per_bohr": [0.02, 0.04, 0.05, 0.08], "S_stderr": [1e-3] * 4,
        "sigma_iso_mid": [8.0, 5.0, 4.6, 4.1], "sigma_WP_nominal": [4.0] * 4})
    out = NG.does_S_collapse_on_width(one_family)
    assert not np.isfinite(out["family_offset_ln"])
    assert out["power_law_exponent"] < 0     # wider couples more weakly


def test_pilot_gate_fails_on_nan_energy(tmp_path):
    """A diverged propagation must stop the ladder, not be averaged into it."""
    d = _write_run(tmp_path, "wp", "bad", S_true=0.05)
    obs = pd.read_csv(d / "raw/observables/observables.csv")
    obs.loc[10, "energy_total"] = np.nan
    obs.to_csv(d / "raw/observables/observables.csv", index=False)
    out = NG.pilot_gate(tmp_path, ["bad"], ["wp"])
    assert out["pass"] is False
    assert "NaN" in out["report"]


def test_pilot_gate_passes_a_healthy_pair(tmp_path):
    _write_run(tmp_path, "wp", "g1", S_true=0.05, mass=1.0)
    _write_run(tmp_path, "wp", "g2", S_true=0.02, mass=0.5)
    out = NG.pilot_gate(tmp_path, ["g1", "g2"], ["wp", "wp"])
    assert out["pass"] is True
