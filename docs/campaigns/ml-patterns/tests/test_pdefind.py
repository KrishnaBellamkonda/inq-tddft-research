"""Known-case code-tests for the PDE-FIND kernel (T9 gate, ADR 0012).

Recover three KNOWN governing PDEs from synthetic data:
  1. advection      u_t  = -c u_x          (1st order)
  2. diffusion      u_t  =  nu u_xx         (1st order)
  3. wave/plasma    u_tt =  c^2 u_xx - w^2 u (2nd order)
plus checks that forward-integration (Wall 2) and bootstrap stability (Wall 3)
behave. Run: venv/bin/python3 -m pytest docs/campaigns/ml-patterns/tests/test_pdefind.py -q
or execute directly (has a __main__ runner).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from kernels import pdefind as PF


def _coeff(model, name):
    return dict(zip(model.names, model.coeffs)).get(name, 0.0)


def _spurious_max(model, keep):
    """Largest |coeff| among terms NOT in `keep`."""
    return max((abs(c) for c, n in zip(model.coeffs, model.names) if n not in keep),
              default=0.0)


# ---------------------------------------------------------------------------
def make_advection(c=1.3, Nx=200, T=80, L=20.0, dt=0.05, noise=0.0, seed=1):
    x = np.linspace(-L / 2, L / 2, Nx)
    dx = x[1] - x[0]
    t = np.arange(T) * dt
    # translating Gaussian bump: u = exp(-((x-c t)^2)/(2 s^2))
    s = 1.5
    u = np.exp(-((x[None, :] - c * t[:, None]) ** 2) / (2 * s * s))
    if noise:
        u = u + noise * np.random.default_rng(seed).standard_normal(u.shape)
    return u, dx, dt


def make_diffusion(nu=0.4, Nx=200, T=80, L=20.0, dt=0.02, noise=0.0, seed=2):
    x = np.linspace(-L / 2, L / 2, Nx)
    dx = x[1] - x[0]
    t = np.arange(T) * dt
    # sum of decaying Fourier modes on [-L/2,L/2] with period L
    u = np.zeros((T, Nx))
    for k in [1, 2, 3]:
        kk = 2 * np.pi * k / L
        u += np.exp(-nu * kk * kk * t)[:, None] * np.sin(kk * x)[None, :]
    if noise:
        u = u + noise * np.random.default_rng(seed).standard_normal(u.shape)
    return u, dx, dt


def make_wave(c=1.0, w=2.0, Nx=200, T=120, L=20.0, dt=0.02, noise=0.0, seed=3):
    x = np.linspace(-L / 2, L / 2, Nx)
    dx = x[1] - x[0]
    t = np.arange(T) * dt
    # standing modes: dispersion Omega_k^2 = c^2 k^2 + w^2
    u = np.zeros((T, Nx))
    for k in [1, 2, 3]:
        kk = 2 * np.pi * k / L
        Om = np.sqrt(c * c * kk * kk + w * w)
        u += np.cos(Om * t)[:, None] * np.sin(kk * x)[None, :]
    if noise:
        u = u + noise * np.random.default_rng(seed).standard_normal(u.shape)
    return u, dx, dt


# ---------------------------------------------------------------------------
def test_advection():
    u, dx, dt = make_advection(c=1.3)
    m = PF.discover_pde_1d(u, dx, dt, order=1, threshold=0.05, x_margin=6, t_margin=3)
    c = _coeff(m, "u_x")
    assert abs(c - (-1.3)) < 0.15, f"u_x coeff {c} != -1.3\n{m.pretty()}"
    assert _spurious_max(m, {"u_x"}) < 0.15, f"spurious term present: {m.pretty()}"
    print("advection OK:", m.pretty())


def test_diffusion():
    u, dx, dt = make_diffusion(nu=0.4)
    m = PF.discover_pde_1d(u, dx, dt, order=1, threshold=0.03, x_margin=6, t_margin=3)
    c = _coeff(m, "u_xx")
    assert abs(c - 0.4) < 0.08, f"u_xx coeff {c} != 0.4\n{m.pretty()}"
    assert _spurious_max(m, {"u_xx"}) < 0.1, f"spurious term present: {m.pretty()}"
    print("diffusion OK:", m.pretty())


def test_wave_plasma():
    u, dx, dt = make_wave(c=1.0, w=2.0)
    m = PF.discover_pde_1d(u, dx, dt, order=2, threshold=0.05, x_margin=6, t_margin=3)
    cxx = _coeff(m, "u_xx")
    cu = _coeff(m, "u")
    assert abs(cxx - 1.0) < 0.2, f"u_xx coeff {cxx} != 1.0\n{m.pretty()}"
    assert abs(cu - (-4.0)) < 0.5, f"u coeff {cu} != -4.0 (w^2)\n{m.pretty()}"
    assert _spurious_max(m, {"u_xx", "u"}) < 0.3, f"spurious: {m.pretty()}"
    print("wave/plasma OK:", m.pretty())


def test_forward_prediction():
    """Wall 2: the discovered advection PDE should forward-predict held-out time."""
    u, dx, dt = make_advection(c=1.3, T=100)
    m = PF.discover_pde_1d(u, dx, dt, order=1, threshold=0.05, x_margin=6, t_margin=3)
    rel, pred = PF.forward_score(m, u, fit_frac=0.5)
    assert rel < 0.2, f"forward-prediction rel L2 {rel} too high"
    print(f"forward-prediction OK: rel L2 = {rel:.3f}")


def test_bootstrap_stability():
    """Wall 3: the true term is active in ~all resamples; spurious terms are not."""
    u, dx, dt = make_diffusion(nu=0.4)
    bs = PF.bootstrap_stability(u, dx, dt, order=1, n_boot=12, threshold=0.03,
                                x_margin=6, t_margin=3)
    fa = bs["frac_active"]
    assert fa.get("u_xx", 0) > 0.8, f"u_xx unstable: {fa.get('u_xx')}"
    print(f"bootstrap OK: u_xx active {fa.get('u_xx'):.2f} of resamples")


def test_noise_robustness():
    """Diffusion recovered under modest noise with time-smoothing."""
    u, dx, dt = make_diffusion(nu=0.4, noise=0.01)
    m = PF.discover_pde_1d(u, dx, dt, order=1, threshold=0.05, smooth_t=1.0,
                           smooth_x=1.0, x_margin=8, t_margin=4)
    c = _coeff(m, "u_xx")
    assert abs(c - 0.4) < 0.15, f"noisy u_xx coeff {c} != 0.4\n{m.pretty()}"
    print("noise-robust OK:", m.pretty())


ALL = [test_advection, test_diffusion, test_wave_plasma, test_forward_prediction,
       test_bootstrap_stability, test_noise_robustness]

if __name__ == "__main__":
    fails = 0
    for t in ALL:
        try:
            t()
        except AssertionError as e:
            fails += 1
            print(f"FAIL {t.__name__}: {e}")
    print(f"\n{len(ALL) - fails}/{len(ALL)} passed")
    sys.exit(1 if fails else 0)
