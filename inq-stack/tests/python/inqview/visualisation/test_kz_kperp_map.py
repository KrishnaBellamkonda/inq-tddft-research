"""Known-case tests for inqview.visualisation.field_io.kz_kperp_map.

The expectations are ANALYTIC. For a drifting Gaussian orbital

    psi(r) = exp(-|r - b|^2 / (2 sigma^2)) exp(i k0 z)

the momentum density factorises exactly:

    longitudinal   k_z  ~  N(k0, sigma_p^2)          sigma_p = 1/(sqrt2 sigma)
    transverse     k_x, k_y each ~ N(0, sigma_p^2)

so the SHELL-SUMMED transverse marginal (the Jacobian is included by
construction — see the function's docstring) is a RAYLEIGH distribution:

    P(k_perp) = (k_perp / sigma_p^2) exp(-k_perp^2 / (2 sigma_p^2))

with mode at sigma_p, mean sigma_p sqrt(pi/2) and <k_perp^2> = 2 sigma_p^2.

Getting the Rayleigh right is the whole point of the test: a naive reading
expects the transverse marginal to peak at k_perp = 0 (it is a Gaussian in
k_x and k_y, after all), and it does not. A function that dropped the shell
Jacobian would peak at zero and would still look perfectly plausible.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "inq-stack" / "python"))

from inqview.visualisation.field_io import (  # noqa: E402
    ComplexVtiField, kz_kperp_map, kz_marginal)

# Deliberately anisotropic (Lz != Lxy) so an x/z axis swap cannot pass.
# Grid chosen so the TRANSVERSE distribution is actually resolved:
# sigma_p / dk_perp = 3.0. At the obvious first choice (NX=32, sigma=3) that
# ratio is 0.60 -- the whole transverse distribution is narrower than one bin
# and <k_perp^2> comes out 92 % high. That is a property of the grid, not of
# the function, and it is why this fixture is not the "natural-looking" one.
NX = NY = 64
NZ = 80
DX = 0.5
SIGMA = 1.2
K0 = 1.4
SIGMA_P = 1.0 / (np.sqrt(2.0) * SIGMA)


def _gaussian_field(k0z: float = K0, center=(0.0, 0.0, 0.0)) -> ComplexVtiField:
    """A drifting Gaussian in PHYSICAL order, as a VTI loader would return it."""
    x = (np.arange(NX) - NX // 2) * DX
    y = (np.arange(NY) - NY // 2) * DX
    z = (np.arange(NZ) - NZ // 2) * DX
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    dx_, dy_, dz_ = X - center[0], Y - center[1], Z - center[2]
    psi = np.exp(-(dx_**2 + dy_**2 + dz_**2) / (2.0 * SIGMA**2)) * np.exp(1j * k0z * dz_)
    psi = psi / np.sqrt((np.abs(psi) ** 2).sum())
    return ComplexVtiField(data=psi, x=x, y=y, z=z,
                           origin=(x[0], y[0], z[0]), spacing=(DX, DX, DX))


def test_map_is_normalised_and_shaped_correctly():
    kz, kperp, P = kz_kperp_map(_gaussian_field())
    assert P.shape == (kz.size, kperp.size)
    assert kz.size == NZ
    assert P.sum() == pytest.approx(1.0, rel=1e-12)
    assert (P >= 0).all()
    assert np.all(np.diff(kz) > 0), "k_z must be returned sorted ascending"


def test_longitudinal_marginal_recovers_the_drift():
    """Summing out k_perp must give N(k0, sigma_p^2)."""
    kz, _, P = kz_kperp_map(_gaussian_field())
    pz = P.sum(axis=1)
    mean = float((kz * pz).sum())
    var = float(((kz - mean) ** 2 * pz).sum())
    assert mean == pytest.approx(K0, abs=1e-6)
    # EXACT, not approximate: k_z is the native FFT axis and is never binned.
    assert var == pytest.approx(SIGMA_P**2, rel=1e-6)


def test_transverse_marginal_is_rayleigh_not_gaussian():
    """<k_perp^2> = 2 sigma_p^2, and the mode sits at sigma_p, NOT at zero.

    This is what pins the shell Jacobian. Dropping it would put the peak at
    k_perp = 0 and halve <k_perp^2>.
    """
    _, kperp, P = kz_kperp_map(_gaussian_field())
    pp = P.sum(axis=0)

    # The SHAPE is the assertion that matters and it is unambiguous: a Rayleigh
    # marginal peaks away from the origin, a Gaussian one would peak at it.
    mode = kperp[int(np.argmax(pp))]
    assert pp[0] < pp[int(np.argmax(pp))], "a Rayleigh marginal does not peak at 0"
    assert mode == pytest.approx(SIGMA_P, abs=1.5 * (kperp[1] - kperp[0]))

    # <k_perp^2> is only approximate BY CONSTRUCTION -- every point in a bin is
    # assigned the bin centre, and a Rayleigh tail falls steeply across a bin, so
    # the estimate is biased HIGH by a few per cent (+6.3 % here). Verified
    # separately that the raw unbinned moment is exact to 0.00 % on this grid, so
    # this is bin-centre assignment and nothing else. Asserted loosely, and with
    # the sign of the bias pinned, so a regression that lost the Jacobian
    # (which would HALVE it) still fails.
    m2 = float((kperp**2 * pp).sum())
    assert m2 == pytest.approx(2.0 * SIGMA_P**2, rel=0.12)
    assert m2 > 2.0 * SIGMA_P**2, "bin-centre bias is positive; a low value means lost Jacobian"


def test_agrees_with_the_existing_1d_kz_marginal():
    """Cross-check against kz_marginal, which is independently tested."""
    field = _gaussian_field()
    kz_a, p_a = kz_marginal(field)
    kz_b, _, P = kz_kperp_map(field)
    dk = kz_b[1] - kz_b[0]
    assert np.allclose(kz_a, kz_b)
    assert np.allclose(p_a, P.sum(axis=1) / dk, rtol=1e-10, atol=1e-12)


def test_drift_moves_weight_along_kz_only():
    """A slower packet shifts in k_z and does NOT change the transverse spread.

    This is the discrimination the function exists for: pure deceleration must
    be visible on the k_z axis and invisible on the k_perp axis.
    """
    _, kperp, P_fast = kz_kperp_map(_gaussian_field(k0z=K0))
    kz, _, P_slow = kz_kperp_map(_gaussian_field(k0z=0.5 * K0))

    mean_fast = float((kz * P_fast.sum(axis=1)).sum())
    mean_slow = float((kz * P_slow.sum(axis=1)).sum())
    assert mean_slow == pytest.approx(0.5 * K0, abs=1e-6)
    assert mean_fast - mean_slow == pytest.approx(0.5 * K0, abs=1e-6)

    m2 = lambda P: float((kperp**2 * P.sum(axis=0)).sum())
    assert m2(P_fast) == pytest.approx(m2(P_slow), rel=1e-9)


def test_binning_default_is_one_transverse_grid_spacing():
    """Finer bins than the grid spacing manufacture a spiky comb; don't."""
    _, kperp, _ = kz_kperp_map(_gaussian_field())
    dk_xy = 2.0 * np.pi / (NX * DX)
    assert (kperp[1] - kperp[0]) == pytest.approx(dk_xy, rel=0.02)


def test_rejects_an_empty_field():
    f = _gaussian_field()
    with pytest.raises(ValueError, match="zero total momentum weight"):
        kz_kperp_map(ComplexVtiField(data=np.zeros_like(f.data), x=f.x, y=f.y,
                                     z=f.z, origin=f.origin, spacing=f.spacing))
