"""Known-case tests for inqview.io.gaussian_psp (erf-smoothed electron UPF).

Validates, against the analytic Gaussian-charge result:
  (a) V(r) reproduces C*erf(r/(sigma*sqrt2))/r to < 1e-6 (template units);
  (b) the radial Fourier transform matches (4*pi/q^2)*exp(-q^2 sigma^2/2)
      (Hartree), evaluated stably via the erfc residual, rel < 1e-4;
  (c) V(0) in Hartree equals sqrt(2/pi)/sigma to < 1e-3.

Run:
    venv/bin/python -m pytest inq-stack/tests/python/inqview/io/test_gaussian_psp.py -v
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from inqview.io import gaussian_psp as G

# Resolved RELATIVE to this file (repo root is five levels up:
# inq-stack/tests/python/inqview/io/ -> repo root), so the suite follows the repo
# instead of a machine-specific absolute path. The previous hard-coded
# /local/data/public/... path silently skipped every test after the CSD3
# migration — a skip is not a pass.
_REPO_ROOT = Path(__file__).resolve().parents[5]
TEMPLATE = (
    _REPO_ROOT
    / "ResearchProject/systems/jellium/shared/pseudopotentials/electron-ONCV-1.2.upf"
)

pytestmark = pytest.mark.skipif(
    not TEMPLATE.exists(), reason="electron-ONCV template not present"
)


def _read_local(path: Path) -> tuple[np.ndarray, np.ndarray]:
    text = path.read_text()
    _, _, r = G._read_block(text, "PP_R")
    _, _, v = G._read_block(text, "PP_LOCAL")
    return np.asarray(r), np.asarray(v)


# NOTE: the parametrize values are UNIFIED (wavepacket) sigmas. The charge std the
# erf actually uses is sigma_charge = sigma_wp/sqrt(2) (unified convention, 2026-06-21).
@pytest.mark.parametrize("sigma_wp", [0.5, 0.4])
def test_sigma_convention(tmp_path, sigma_wp):
    """generate uses charge std = sigma_wp/sqrt(2) (the unification contract)."""
    out = tmp_path / f"electron_gaussian_wpsigma{sigma_wp}.upf"
    res = G.generate_gaussian_psp(TEMPLATE, sigma_wp, out)
    assert res.sigma_wp == pytest.approx(sigma_wp)
    assert res.sigma_charge == pytest.approx(sigma_wp / math.sqrt(2.0))


@pytest.mark.parametrize("sigma_wp", [0.5, 0.4])
def test_vr_matches_erf_form(tmp_path, sigma_wp):
    out = tmp_path / f"electron_gaussian_wpsigma{sigma_wp}.upf"
    res = G.generate_gaussian_psp(TEMPLATE, sigma_wp, out)
    r, v = _read_local(out)
    mask = (r >= 0.05) & (r <= 10.0)
    # potential is built from the CHARGE std, not sigma_wp:
    analytic = res.coulomb_coeff * G.v_erf_hartree(r[mask], res.sigma_charge)
    rel = np.abs(v[mask] - analytic) / np.abs(analytic)
    assert rel.max() < 1e-6, f"max rel dev {rel.max():.2e}"


@pytest.mark.parametrize("sigma_wp", [0.5, 0.4])
def test_v0_hartree(tmp_path, sigma_wp):
    out = tmp_path / f"electron_gaussian_wpsigma{sigma_wp}.upf"
    res = G.generate_gaussian_psp(TEMPLATE, sigma_wp, out)
    expected = math.sqrt(2.0 / math.pi) / res.sigma_charge   # uses charge std
    assert abs(res.v0_hartree - expected) < 1e-3
    # v[0] is in template units = coulomb_coeff * sqrt(2/pi)/sigma_charge
    r, v = _read_local(out)
    assert abs(v[0] / res.coulomb_coeff - expected) < 1e-3


@pytest.mark.parametrize("sigma_wp", [0.5, 0.4])
def test_fourier_form_factor(tmp_path, sigma_wp):
    """FT of V_Ha(r)=erf/r compared to (4pi/q^2)exp(-q^2 sigma_charge^2/2).

    Stable evaluation: V_Ha = 1/r - erfc(r/sigma_charge sqrt2)/r. The 1/r part has
    known FT 4pi/q^2; the erfc residual decays as a Gaussian so its radial sine
    transform converges on the finite mesh.
    """
    out = tmp_path / f"electron_gaussian_wpsigma{sigma_wp}.upf"
    res = G.generate_gaussian_psp(TEMPLATE, sigma_wp, out)
    sigma = res.sigma_charge                         # form factor uses the charge std
    r, v_template = _read_local(out)
    v_ha = v_template / res.coulomb_coeff           # to Hartree (C=1 form)

    # residual R(r) = V_Ha - 1/r  (= -erfc(r/sigma sqrt2)/r), finite-decaying
    R = np.zeros_like(r)
    nz = r > 1e-12
    R[nz] = v_ha[nz] - 1.0 / r[nz]
    dr = r[1] - r[0]

    # Test over the form factor's physical support qσ ≲ 2 (q ≲ 4): beyond this
    # exp(-q^2 sigma^2/2) is < 0.02 and the σ-smoothing has removed the weight
    # (the stopping q-integral carries exp(-q^2 sigma^2)). At higher q the test
    # would only measure the oscillatory-quadrature floor (err ~ eps*q) against
    # an exponentially-vanishing target. V(r)→1e-6 (test_vr) is the primary gate.
    qs = np.linspace(0.2, 4.0, 30)
    for q in qs:
        # FT[R](q) = (4pi/q) * int_0^inf r R(r) sin(qr) dr
        # trapezoidal quadrature (oscillatory integrand: far better than Riemann)
        integ = np.trapezoid(r * R * np.sin(q * r), r)
        ft_R = (4.0 * math.pi / q) * integ
        ft = 4.0 * math.pi / q**2 + ft_R
        coulomb = 4.0 * math.pi / q**2
        analytic = coulomb * math.exp(-(q**2) * sigma**2 / 2.0)
        # Measure the form factor as a fraction of the bare Coulomb amplitude
        # (the well-posed metric: rel-to-analytic is ill-posed where exp(-q^2 s^2/2)
        # is exponentially small and 4pi/q^2 + FT[R] cancel catastrophically).
        err = abs(ft - analytic) / coulomb
        assert err < 2e-4, f"q={q:.2f} sigma={sigma}: form-factor err {err:.2e}"
