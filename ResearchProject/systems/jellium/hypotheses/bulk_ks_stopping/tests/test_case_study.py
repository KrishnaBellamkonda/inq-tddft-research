"""Tests for the 100 eV high-density case study (``case_study_100eV/make_case_study.py``).

Three things are checked, and they fail for different reasons:

1. PROVENANCE PARSING. ``Meta`` reads every physical constant out of the run's own
   artefacts instead of hardcoding them, so a parsing slip would silently retitle
   and mislabel all sixteen figures. Tested against synthetic artefacts with
   PLANTED values whose r_s and omega_p are known analytically.

2. THE KINEMATIC KERNELS. The Ehrenfest centroid is an integral, and T_2 is a
   contraction of three means; both are cheap to get subtly wrong (trapezoid vs
   rectangle, <p>^2 vs <p^2>) and neither error would look wrong on a plot.
   Tested on synthetic runs whose exact answer is known in closed form.

3. THE LAYOUT GUARDS. A fixed-canvas figure CROPS overrunning text silently --
   the title is simply absent from the PNG and the build reports success. Both
   guards therefore carry a NEGATIVE self-test: each is shown to fire on the real
   defect it exists to catch, because a guard that never fires is worse than no
   guard (it reads as evidence).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

FAMILY = "bulk_ks_stopping_rs4"

# NOTE ON PATHS: hypotheses/<family>/tests is a SYMLINK to
# hypotheses/bulk_ks_stopping/tests, so Path(__file__).resolve() always lands in
# the canonical suite regardless of which family directory the test was invoked
# through. Anchor on the hypotheses root and name the family explicitly rather
# than walking up from __file__ and assuming the parent is the right family.
HERE = Path(__file__).resolve().parent          # .../hypotheses/bulk_ks_stopping/tests
HYPOTHESES = HERE.parents[1]                    # .../systems/jellium/hypotheses
SCRIPTS = HERE.parents[2] / "scripts"           # .../systems/jellium/scripts
CASE_DIR = HYPOTHESES / FAMILY / "case_study_100eV"
BUILDER = CASE_DIR / "make_case_study.py"


def _load_builder():
    spec = importlib.util.spec_from_file_location("_case_study_builder", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


pytest.importorskip("matplotlib")
M = _load_builder()


# ---------------------------------------------------------------------------
# 1. Provenance parsing
# ---------------------------------------------------------------------------

def _synthetic_run(tmp_path: Path, *, n_elec: int, lz: float = 80.0,
                   lxy: float = 40.0, sigma: float = 2.0) -> Path:
    """A run tree carrying only what ``Meta`` reads."""
    for half in ("wp", "classical"):
        (tmp_path / "fam" / half / "results" / "raw" / "observables").mkdir(
            parents=True, exist_ok=True)
    res = tmp_path / "fam" / "wp" / "results"
    (res / "run_summary.txt").write_text(
        "RUN SUMMARY\n===========\n\n"
        "3. System configuration\n-----------------------\n"
        f"cell_bohr       = {lxy:g} x {lxy:g} x {lz:g}\n"
        f"n_electrons     = {n_elec}\n"
        f"wp_sigma_bohr   = {sigma:g}\n"
        "wp_center_bohr  = 0 0 -32\n"
        "wp_k0_bohr_inv  = 2.71106334010243\n"
        "wp_energy_ev    = 100\n"
        "dt_au           = 0.04\n")
    (res / "raw" / "observables" / "wp_config.txt").write_text(
        "fit_t0_au      = 4\nfit_t1_au      = 18.42703328144701\n")
    return tmp_path


def test_meta_parses_the_runs_own_artefacts(tmp_path):
    root = _synthetic_run(tmp_path, n_elec=482)
    meta = M.Meta("fam", root)
    assert (meta.lx, meta.ly, meta.lz) == (40.0, 40.0, 80.0)
    assert meta.n_elec == 482
    assert meta.sigma_wp == 2.0
    assert meta.z0 == -32.0
    assert meta.dt == 0.04
    assert meta.fit_t0 == 4.0
    assert meta.fit_t1 == pytest.approx(18.42703328144701)


def test_r_s_matches_the_analytic_value_for_a_planted_density(tmp_path):
    # Choose N so that r_s is exactly 2 Bohr: n = 3/(4 pi r_s^3), N = n * V.
    lxy, lz, r_s_target = 40.0, 80.0, 2.0
    volume = lxy * lxy * lz
    n_exact = 3.0 / (4.0 * np.pi * r_s_target ** 3)
    n_elec = int(round(n_exact * volume))
    root = _synthetic_run(tmp_path, n_elec=n_elec, lz=lz, lxy=lxy)
    meta = M.Meta("fam", root)
    # Rounding N to an integer perturbs r_s at the 1e-5 level, no more.
    assert meta.r_s == pytest.approx(r_s_target, rel=1e-4)


def test_omega_p_is_sqrt_4_pi_n_in_ev(tmp_path):
    root = _synthetic_run(tmp_path, n_elec=482)
    meta = M.Meta("fam", root)
    expected = np.sqrt(4.0 * np.pi * meta.density) * M.HA_EV
    assert meta.omega_p_ev == pytest.approx(expected, rel=1e-12)


def test_meta_refuses_a_missing_half(tmp_path):
    root = _synthetic_run(tmp_path, n_elec=482)
    (root / "fam" / "classical").rename(root / "fam" / "classical_moved")
    with pytest.raises(FileNotFoundError):
        M.Meta("fam", root)


# ---------------------------------------------------------------------------
# 2. Kinematic kernels — closed-form cases
# ---------------------------------------------------------------------------

def _synthetic_wp(tmp_path: Path, t, px, py, pz, z_of_t) -> Path:
    """Write a minimal pair of WP stats CSVs with EXACT planted values.

    Written via pandas rather than f-string formatting: under numpy 2,
    ``repr(np.float64(0.3))`` is ``"np.float64(0.3)"``, which pandas reads back
    as a string column and turns a numeric test into a TypeError. ``to_csv``
    round-trips full float precision, which these exactness assertions need.
    """
    import pandas as pd

    obs = tmp_path / "wp" / "results" / "raw" / "observables"
    obs.mkdir(parents=True, exist_ok=True)
    step = np.arange(len(t))
    pd.DataFrame({
        "step": step, "time_au": t,
        "px_mean": px, "py_mean": py, "pz_mean": pz,
        # <p^2>/2m with zero momentum spread planted, so T_1 == T_2 exactly.
        "e_kin_ha": 0.5 * (np.asarray(px) ** 2 + np.asarray(py) ** 2
                           + np.asarray(pz) ** 2),
        "norm_check": 1.0,
    }).to_csv(obs / "wp_momentum_stats.csv", index=False)
    pd.DataFrame({
        "step": step, "time_au": t,
        "z_mean": z_of_t, "z_mean_circ": z_of_t,
        "sigma_z_circ": 2.0, "norm_check": 1.0,
    }).to_csv(obs / "wp_real_space_stats.csv", index=False)
    return tmp_path / "wp"


def test_ehrenfest_centroid_is_exact_for_constant_momentum(tmp_path):
    """z_0 + integral <p_z> dt with constant p_z must be exactly z_0 + p_z t.

    The trapezoid rule is exact on a constant integrand, so any deviation here
    is an indexing or dt bug, not discretisation.
    """
    t = np.arange(0.0, 10.0001, 0.04)
    p0, z0 = 2.71106334010243, -32.0
    pz = np.full(len(t), p0)
    run = _synthetic_wp(tmp_path, t, np.zeros(len(t)), np.zeros(len(t)), pz,
                        z0 + p0 * t)
    wp = M.K.load_wp_run(run, box_length_z=80.0, z0=z0)
    np.testing.assert_allclose(wp.s4, z0 + p0 * t, rtol=0, atol=1e-12)


def test_ehrenfest_centroid_is_exact_for_linear_momentum(tmp_path):
    """A linearly decelerating p_z integrates to a parabola, again exactly."""
    t = np.arange(0.0, 10.0001, 0.04)
    p0, a, z0 = 2.7, -0.05, -32.0
    pz = p0 + a * t
    run = _synthetic_wp(tmp_path, t, np.zeros(len(t)), np.zeros(len(t)), pz,
                        np.zeros(len(t)))
    wp = M.K.load_wp_run(run, box_length_z=80.0, z0=z0)
    np.testing.assert_allclose(wp.s4, z0 + p0 * t + 0.5 * a * t ** 2,
                               rtol=0, atol=1e-10)


def test_T2_is_the_square_of_the_mean_not_the_mean_of_the_square(tmp_path):
    """T_2 = <p>^2/2m must use the MEAN momentum in all three components."""
    t = np.array([0.0, 0.04, 0.08])
    px = np.array([0.3, 0.3, 0.3]); py = np.array([-0.4, -0.4, -0.4])
    pz = np.array([2.0, 2.0, 2.0])
    run = _synthetic_wp(tmp_path, t, px, py, pz, np.zeros(3))
    wp = M.K.load_wp_run(run, box_length_z=80.0, z0=0.0)
    np.testing.assert_allclose(wp.T2, 0.5 * (0.09 + 0.16 + 4.0), rtol=0, atol=1e-14)


def test_spread_energy_is_the_difference_of_the_two(tmp_path):
    """With zero spread planted, T_1 - T_2 must vanish identically."""
    t = np.arange(0.0, 1.0001, 0.04)
    pz = np.full(len(t), 2.5)
    run = _synthetic_wp(tmp_path, t, np.zeros(len(t)), np.zeros(len(t)), pz,
                        np.zeros(len(t)))
    wp = M.K.load_wp_run(run, box_length_z=80.0, z0=0.0)
    np.testing.assert_allclose(wp.T1 - wp.T2, 0.0, atol=1e-14)


# ---------------------------------------------------------------------------
# 3. Layout guards, each with a negative self-test
# ---------------------------------------------------------------------------

def _png(tmp_path: Path, name: str, *, touch_edge: bool):
    from PIL import Image
    a = np.full((200, 300), 255, dtype=np.uint8)
    if touch_edge:
        a[0, 40:260] = 0            # ink on the very first row
    else:
        a[60:140, 60:240] = 0       # ink well inside
    Image.fromarray(a).save(tmp_path / name)


def test_margin_check_passes_a_clean_figure(tmp_path):
    _png(tmp_path, "clean.png", touch_edge=False)
    assert M.verify_margins(tmp_path) == 0


def test_margin_check_FIRES_on_a_figure_clipped_at_the_edge(tmp_path):
    """Negative self-test: the real defect must be detected.

    This is the exact failure the guard was written for -- a two-line title
    flush against row 0, which crops at save time and leaves no trace in the
    build log.
    """
    _png(tmp_path, "clipped.png", touch_edge=True)
    assert M.verify_margins(tmp_path) == 1


def test_margin_check_reports_a_blank_figure(tmp_path):
    from PIL import Image
    Image.fromarray(np.full((100, 100), 255, dtype=np.uint8)).save(
        tmp_path / "blank.png")
    assert M.verify_margins(tmp_path) == 1


def test_overrunning_title_is_shrunk_to_fit():
    import matplotlib
    matplotlib.use("Agg")
    fig, ax = M.fig1()
    ax.set_title("W" * 200, fontsize=9)
    shrunk = M._shrink_overrunning_titles(fig)
    assert shrunk, "a 200-character title must be reported as shrunk"
    assert ax.title.get_fontsize() < 9.0
    matplotlib.pyplot.close(fig)


def test_a_short_title_is_left_alone():
    """Positive control: the guard must not fiddle with a title that fits."""
    import matplotlib
    matplotlib.use("Agg")
    fig, ax = M.fig1()
    ax.set_title("short", fontsize=9)
    assert M._shrink_overrunning_titles(fig) == []
    assert ax.title.get_fontsize() == 9.0
    matplotlib.pyplot.close(fig)


def test_every_math_span_in_the_builder_renders():
    """Mathtext is valid TeX-shaped but NOT TeX: ``\\frac12`` parses as Python and
    reads as TeX yet raises in matplotlib. Only the real parser can judge."""
    import ast
    import re
    from matplotlib.mathtext import MathTextParser

    tree = ast.parse(BUILDER.read_text(), filename=str(BUILDER))
    values = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            values.append(node.value)
        elif isinstance(node, ast.JoinedStr):
            # RUNTIME value: an f-string's placeholders are not in the source.
            values.append("".join(
                v.value if isinstance(v, ast.Constant) and isinstance(v.value, str)
                else "1.0" for v in node.values))

    parser = MathTextParser("agg")
    spans, failures = 0, []
    for v in values:
        if "$" not in v:
            continue
        assert v.count("$") % 2 == 0, f"unbalanced $ in {v!r}"
        for m in re.finditer(r"\$([^$]+)\$", v):
            spans += 1
            try:
                parser.parse(f"${m.group(1)}$", dpi=72, prop=None)
            except Exception as exc:
                failures.append((m.group(0), f"{type(exc).__name__}: {exc}"))
    assert spans > 20, f"expected the builder to carry math spans, found {spans}"
    assert not failures, f"mathtext failures: {failures}"


def test_the_mathtext_guard_is_not_vacuous():
    """Negative self-test for the check above."""
    from matplotlib.mathtext import MathTextParser
    with pytest.raises(Exception):
        MathTextParser("agg").parse(r"$\frac12 mv^2$", dpi=72, prop=None)


def test_no_control_characters_in_builder_string_literals():
    r"""``\a``/``\b``/``\f``/``\v`` from ``\approx``/``\beta``/``\frac``/``\varphi``
    in a non-raw string compile fine and render as garbage."""
    import ast
    tree = ast.parse(BUILDER.read_text(), filename=str(BUILDER))
    bad = {chr(c) for c in range(32)} - {"\t", "\n"}
    offenders = [n.value for n in ast.walk(tree)
                 if isinstance(n, ast.Constant) and isinstance(n.value, str)
                 and bad & set(n.value)]
    assert not offenders, f"control characters in: {offenders!r}"


# ---------------------------------------------------------------------------
# 4. Integration against the real run (skipped if it is not on disk)
# ---------------------------------------------------------------------------

FAMILY = "bulk_ks_stopping_rs4"
_have_run = (SCRIPTS / FAMILY / "wp" / "results" / "raw" / "observables"
             / "wp_momentum_stats.csv").exists()
needs_run = pytest.mark.skipif(not _have_run, reason=f"{FAMILY} run not on disk")


@needs_run
def test_the_three_stopping_powers_are_additive():
    """S(T_1) = S(T_2) + S(T_var), because OLS slope is linear in the ordinate.

    T_1 = T_2 + T_var identically, so the three fitted slopes over the SAME
    window must add. A failure means the fits were run over different windows or
    against different path coordinates -- exactly the kind of mismatch that
    produces a plausible but wrong S.
    """
    meta = M.Meta(FAMILY, SCRIPTS)
    wp = M.K.load_wp_run(meta.wp_dir, meta.lz, meta.z0)
    kw = dict(t=wp.t, t0=meta.fit_t0, t1=meta.fit_t1)
    s_T1 = M.K.fit_stopping(wp.s4, wp.T1, label="T1", **kw).S_ev_per_bohr
    s_T2 = M.K.fit_stopping(wp.s4, wp.T2, label="T2", **kw).S_ev_per_bohr
    s_var = M.K.fit_stopping(wp.s4, wp.T1 - wp.T2, label="var", **kw).S_ev_per_bohr
    assert s_T1 == pytest.approx(s_T2 + s_var, rel=1e-10)


@needs_run
def test_the_packet_starts_at_its_zero_point_momentum_spread():
    """T_1 - T_2 at t=0 must equal 3/(4 sigma_psi^2) for a minimum-uncertainty
    Gaussian of psi-width sigma_psi (.claude/rules/sigma-wp-convention.md)."""
    meta = M.Meta(FAMILY, SCRIPTS)
    wp = M.K.load_wp_run(meta.wp_dir, meta.lz, meta.z0)
    expected_ha = 3.0 / (4.0 * meta.sigma_wp ** 2)
    assert (wp.T1[0] - wp.T2[0]) == pytest.approx(expected_ha, rel=1e-6)


@needs_run
def test_drift_kinetic_energy_starts_at_the_nominal_beam_energy():
    """T_2(0) must be the run's stated projectile energy -- NOT T_1(0), which is
    higher by the zero-point spread. Getting this backwards is the single
    easiest way to mislabel the whole comparison."""
    meta = M.Meta(FAMILY, SCRIPTS)
    wp = M.K.load_wp_run(meta.wp_dir, meta.lz, meta.z0)
    assert wp.T2[0] * M.HA_EV == pytest.approx(meta.energy_ev, rel=1e-6)
    assert wp.T1[0] * M.HA_EV > meta.energy_ev


@needs_run
def test_classical_and_wavepacket_share_the_t0_projectile_self_energy():
    """E_PP(0) must agree between the halves: the classical UPF is generated at
    sigma_pot = sigma_WP/sqrt(2) precisely so its cloud matches the packet's t=0
    density. This is a check ON the sigma convention, not an input to it."""
    meta = M.Meta(FAMILY, SCRIPTS)
    cl = M.K.load_interactions(meta.cl_dir, "classical")
    wp = M.K.load_interactions(meta.wp_dir, "wp")
    assert cl.e_pp[0] == pytest.approx(wp.e_pp[0], rel=1e-4)


@needs_run
@pytest.mark.parametrize("half", ["classical", "wp"])
def test_the_three_background_terms_are_exactly_zero_in_bulk(half):
    """Figures 14-16 plot all SIX pairwise terms and assert on the canvas that
    E_SB, E_PB and E_BB vanish. Bulk has a uniform background, so poisson(n_+)
    is pure G=0 -- which INQ drops -- and phi_+ is identically zero.

    Asserted BITWISE, not to a tolerance: these are structural zeros, so any
    non-zero value at all means the background was not uniform (or the columns
    were mis-wired), and the annotation on those figures would be false.
    """
    meta = M.Meta(FAMILY, SCRIPTS)
    ix = M.K.load_interactions(meta.wp_dir if half == "wp" else meta.cl_dir, half)
    for name in ("e_sb", "e_pb", "e_bb"):
        arr = getattr(ix, name)
        assert np.all(arr == 0.0), f"{half}/{name} max|E| = {np.max(np.abs(arr)):.3e}"


@needs_run
def test_classical_projectile_self_energy_is_constant():
    """A rigid Gaussian UPF cannot spread, so E_PP must not move. If it does, the
    classical cloud is being rebuilt at the wrong width or the wrong position."""
    meta = M.Meta(FAMILY, SCRIPTS)
    cl = M.K.load_interactions(meta.cl_dir, "classical")
    assert np.ptp(cl.e_pp) < 1e-9
