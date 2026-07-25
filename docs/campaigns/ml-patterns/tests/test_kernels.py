"""Known-case code-tests for ml-patterns kernels (T1 pre-gate).

Run: venv/bin/python3 docs/campaigns/ml-patterns/tests/test_kernels.py
Each test embeds an INDEPENDENTLY known answer (synthetic ground truth), per the
code-test skill: write -> known-case-test -> confirm. No circular reuse of the
kernel to generate its own expected value.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from kernels import pod as P, dmd as D, formfactor as FF  # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")


# ---------------------------------------------------------------- POD
def test_pod_rank2():
    """A field built from exactly 2 spatial modes -> POD recovers 2 modes,
    ~100% energy in the first two, and the recovered 2D subspace spans the truth."""
    rng = np.random.default_rng(1)
    n, m = 500, 80
    phi1 = rng.standard_normal(n); phi1 /= np.linalg.norm(phi1)
    phi2 = rng.standard_normal(n)
    phi2 -= (phi2 @ phi1) * phi1; phi2 /= np.linalg.norm(phi2)
    t = np.linspace(0, 4 * np.pi, m)
    a1, a2 = 3.0 * np.cos(t), 1.0 * np.sin(0.7 * t)
    X = np.outer(phi1, a1) + np.outer(phi2, a2)
    res = P.pod(X, rank=10)
    e2 = res.cumulative_energy[1]
    # subspace recovery: projection of phi1,phi2 onto first 2 POD modes ~ 1
    Umodes = res.modes[:, :2]
    proj1 = np.linalg.norm(Umodes.T @ phi1)
    proj2 = np.linalg.norm(Umodes.T @ phi2)
    check("POD: first 2 modes capture ~all energy", e2 > 0.999, f"E2={e2:.5f}")
    check("POD: mode-3 energy negligible", res.energy_fraction[2] < 1e-6,
          f"f3={res.energy_fraction[2]:.2e}")
    check("POD: truth subspace recovered", proj1 > 0.999 and proj2 > 0.999,
          f"p1={proj1:.4f} p2={proj2:.4f}")


def test_pod_randomized_matches():
    rng = np.random.default_rng(2)
    X = rng.standard_normal((300, 60)).astype(np.float32)
    a = P.pod(X, rank=5, randomized=False)
    b = P.pod(X, rank=5, randomized=True)
    # singular values agree
    rel = np.abs(a.singular_values[:5] - b.singular_values[:5]) / a.singular_values[:5]
    check("POD: randomized ~ deterministic SVD", np.max(rel) < 0.05,
          f"max rel={np.max(rel):.3f}")


# ---------------------------------------------------------------- DMD
def test_dmd_damped_sinusoid():
    """Snapshots x(t)=exp(-g t)[cos(w t) p1 + sin(w t) p2]; DMD must recover w and -g."""
    rng = np.random.default_rng(3)
    n = 400
    p1 = rng.standard_normal(n); p2 = rng.standard_normal(n)
    w_true, g_true = 1.3, 0.07
    dt = 0.05
    T = 200
    t = np.arange(T) * dt
    X = np.empty((n, T))
    for k in range(T):
        X[:, k] = np.exp(-g_true * t[k]) * (np.cos(w_true * t[k]) * p1
                                            + np.sin(w_true * t[k]) * p2)
    res = D.dmd(X, dt=dt, rank=2)
    i, w_rec, g_rec, amp = res.dominant()
    check("DMD: angular frequency recovered", abs(w_rec - w_true) / w_true < 0.01,
          f"w={w_rec:.4f} vs {w_true}")
    check("DMD: decay rate recovered", abs(g_rec - (-g_true)) < 0.01,
          f"g={g_rec:.4f} vs {-g_true}")


def test_dmd_window():
    rng = np.random.default_rng(4)
    n = 200
    p = rng.standard_normal(n)
    dt = 0.1
    T = 100
    t = np.arange(T) * dt
    w_true = 2.0
    X = np.outer(p, np.cos(w_true * t)) + np.outer(rng.standard_normal(n) * 0.0,
                                                   np.sin(w_true * t))
    # add a second pattern in quadrature for a clean complex pair
    p2 = rng.standard_normal(n)
    X = X + np.outer(p2, np.sin(w_true * t))
    res = D.dmd(X, dt=dt, rank=2, window=(10, 90))
    _, w_rec, _, _ = res.dominant()
    check("DMD: windowed frequency recovered", abs(w_rec - w_true) / w_true < 0.02,
          f"w={w_rec:.4f} vs {w_true}")


# ---------------------------------------------------------------- form factor
def test_FWP_analytic():
    q = np.array([0.0, 0.5, 1.0, 2.0])
    s = 0.7071067811865475
    expect = np.exp(-0.5 * (q * s) ** 2)
    got = FF.F_WP(q, s)
    check("F_WP: matches exp(-q^2 sigma^2/2)", np.allclose(got, expect),
          f"max err={np.max(np.abs(got-expect)):.2e}")


def test_radial_spectrum_gaussian():
    """FFT magnitude of a real-space Gaussian of std s is Gaussian exp(-q^2 s^2/2);
    radial_power_spectrum must recover the width s within a few %."""
    N, dx = 64, 0.5
    s = 1.5
    ax = (np.arange(N) - N // 2) * dx
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    r2 = X ** 2 + Y ** 2 + Z ** 2
    field = np.exp(-r2 / (2 * s ** 2))
    q, amp = FF.radial_power_spectrum(field, dx, nbins=40)
    amp = amp / amp[0]
    sel = (q > 0) & (amp > 1e-3)
    # fit ln amp = -q^2 s_rec^2/2  -> slope on q^2
    coef = np.polyfit(q[sel] ** 2, np.log(amp[sel]), 1)
    s_rec = np.sqrt(-2 * coef[0])
    check("radial spectrum: Gaussian width recovered", abs(s_rec - s) / s < 0.05,
          f"s_rec={s_rec:.3f} vs {s}")


def test_q_ratio_gaussians():
    """R(q) of two Gaussian fields (s1,s2) follows exp(-q^2(s1^2-s2^2)/2)."""
    N, dx = 64, 0.5
    s1, s2 = 1.8, 1.0
    ax = (np.arange(N) - N // 2) * dx
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    r2 = X ** 2 + Y ** 2 + Z ** 2
    f1 = np.exp(-r2 / (2 * s1 ** 2))
    f2 = np.exp(-r2 / (2 * s2 ** 2))
    q, R = FF.q_ratio(f1, f2, dx, nbins=40)
    R = R / R[0]
    pred = np.exp(-0.5 * q ** 2 * (s1 ** 2 - s2 ** 2))
    pred = pred / pred[0]
    sel = (q > 0) & (q < 1.5)
    rel = np.abs(R[sel] - pred[sel]) / np.maximum(pred[sel], 1e-3)
    check("q_ratio: Gaussian ratio follows exp(-q^2 ds^2/2)", np.median(rel) < 0.1,
          f"median rel={np.median(rel):.3f}")


def test_foncv_upf():
    upf = ("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
           "shared/pseudopotentials/electron-ONCV-1.2.upf")
    if not os.path.isfile(upf):
        check("F_ONCV: UPF present", False, "missing UPF"); return
    q = np.linspace(0, 6, 200)
    qmax, F = FF.foncv_unity_range(upf, q, z=1.0, tol=0.05)
    # F at very low q ~ 1; rolls off (decreases) at high q
    qpos = q > 0
    check("F_ONCV: ~1 at low q", abs(F[qpos][0] - 1.0) < 0.05, f"F(qmin)={F[qpos][0]:.3f}")
    check("F_ONCV: unity-range positive", qmax > 0.3, f"qmax(|F-1|<5%)={qmax:.2f} 1/Bohr")
    print(f"     F_ONCV unity (5%) q-range: 0 < q <= {qmax:.2f} 1/Bohr")


if __name__ == "__main__":
    for fn in [test_pod_rank2, test_pod_randomized_matches,
               test_dmd_damped_sinusoid, test_dmd_window,
               test_FWP_analytic, test_radial_spectrum_gaussian,
               test_q_ratio_gaussians, test_foncv_upf]:
        try:
            fn()
        except Exception as e:
            FAIL.append(fn.__name__)
            print(f"[ERROR] {fn.__name__}: {e!r}")
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    sys.exit(1 if FAIL else 0)
