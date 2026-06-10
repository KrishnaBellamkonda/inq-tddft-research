"""Tests for Python centre-of-density (``inqview.analysis.center_of_density``).

Pure numpy on synthetic densities (IV-M02). Asserts the node-convention centroid
against an analytic value, documents the E04 dx/2 offset, and checks the
WP/total/bath three-way comparison with bath = total − wp.
"""
from __future__ import annotations

import numpy as np
import pytest

import inqview.analysis.center_of_density as cod

pytestmark = pytest.mark.analysis

ORIGIN = (-5.0, -5.0, -5.0)
SPACING = (0.5, 0.5, 0.5)
NX = NY = NZ = 20


def _gaussian(cx, cy, cz, sigma=1.0):
    """A Gaussian density centred at (cx,cy,cz) Bohr on the node grid."""
    x = ORIGIN[0] + np.arange(NX) * SPACING[0]
    y = ORIGIN[1] + np.arange(NY) * SPACING[1]
    z = ORIGIN[2] + np.arange(NZ) * SPACING[2]
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    return np.exp(-((X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2) / (2 * sigma ** 2))


def test_single_point_centroid_is_its_node_coordinate():
    rho = np.zeros((NX, NY, NZ))
    rho[8, 12, 4] = 1.0
    c = cod.center_of_density(rho, ORIGIN, SPACING)
    assert (c.x, c.y, c.z) == pytest.approx(
        (ORIGIN[0] + 8 * SPACING[0], ORIGIN[1] + 12 * SPACING[1],
         ORIGIN[2] + 4 * SPACING[2]))


def test_gaussian_centroid_recovers_centre():
    # tol 1e-2 absorbs the small centroid pull from asymmetric tail truncation
    # at the box edge (a fixture artifact; the COD math is exact — see the
    # single-point and half-cell tests).
    c = cod.center_of_density(_gaussian(1.0, -1.0, 0.5), ORIGIN, SPACING)
    assert (c.x, c.y, c.z) == pytest.approx((1.0, -1.0, 0.5), abs=1e-2)
    assert c.total_weight > 0.0


def test_half_cell_convention_is_offset_by_exactly_dx_over_2():
    """Documents E04: the (i+½)·dx convention shifts COD by +dx/2 per axis."""
    rho = _gaussian(0.0, 0.0, 0.0)
    node = cod.center_of_density(rho, ORIGIN, SPACING, half_cell=False)
    half = cod.center_of_density(rho, ORIGIN, SPACING, half_cell=True)
    assert half.x - node.x == pytest.approx(0.5 * SPACING[0], abs=1e-9)
    assert half.y - node.y == pytest.approx(0.5 * SPACING[1], abs=1e-9)
    assert half.z - node.z == pytest.approx(0.5 * SPACING[2], abs=1e-9)


def test_bath_is_total_minus_wp():
    wp = _gaussian(2.0, 0.0, 0.0, sigma=0.8)
    bath = _gaussian(-2.0, 0.0, 0.0, sigma=1.0)   # contained well inside the box
    total = wp + bath
    cmp = cod.compare(total, wp, ORIGIN, SPACING)
    assert cmp.wp.x == pytest.approx(2.0, abs=2e-2)
    assert cmp.bath.x == pytest.approx(-2.0, abs=2e-2)
    # total centroid lies strictly between the two
    assert -2.0 < cmp.total.x < 2.0


def test_cod_offset_vs_inqkit_recovers_dx_over_2():
    """From-run cross-check: inqkit (half-cell) COD − python (node) COD averaged
    over a series equals (dx,dy,dz)/2 (documents E04). Pure, no VTK/real run."""
    from inqview.pipeline.cod import cod_offset_vs_inqkit

    spacing = np.array([0.4, 0.5, 0.6])
    rng = np.arange(5)
    node = np.stack([rng * 0.1, rng * 0.2, -rng * 0.05], axis=1)   # arbitrary node COD series
    inqkit = node + 0.5 * spacing[None, :]                         # half-cell shift
    off = cod_offset_vs_inqkit(node, inqkit)
    assert np.allclose(off, 0.5 * spacing, atol=1e-12)


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
