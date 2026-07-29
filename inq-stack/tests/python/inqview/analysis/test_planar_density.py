"""Known-case tests for inqview.analysis.planar_density (Δn(z,t))."""
from __future__ import annotations

import numpy as np
import pytest

from inqview.analysis.planar_density import (
    planar_delta_map,
    planar_profile,
)


def test_planar_profile_sums_transverse_plane():
    """planar_profile keeps the z axis and sums x,y; cell_area scales it."""
    nx, ny, nz = 3, 4, 5
    cube = np.ones((nx, ny, nz))
    prof = planar_profile(cube, axis=2)
    assert prof.shape == (nz,)
    assert np.allclose(prof, nx * ny)                  # 12 per z-plane
    # cell_area turns the sum into an integral
    prof2 = planar_profile(cube, axis=2, cell_area=0.25)
    assert np.allclose(prof2, nx * ny * 0.25)


def test_planar_profile_picks_a_known_z_slab():
    """A density localized on one z-slice shows up only there."""
    cube = np.zeros((2, 2, 4))
    cube[:, :, 2] = 3.0                                 # all mass at z-index 2
    prof = planar_profile(cube, axis=2)
    assert np.allclose(prof, [0, 0, 12, 0])            # 2*2*3 = 12


def test_delta_map_reference_column_is_zero_and_increment_tracked():
    nz, nt = 4, 3
    z = np.linspace(-6, 6, nz)
    t = np.array([0.0, 0.1, 0.2])
    base = np.ones((2, 2, nz))
    f0 = base.copy()
    f1 = base.copy(); f1[:, :, 1] += 0.5               # +0.5 at z-index 1
    f2 = base.copy(); f2[:, :, 1] += 1.0
    dmap = planar_delta_map([f0, f1, f2], t, z, axis=2)
    assert dmap.dn.shape == (nz, nt)
    assert np.allclose(dmap.dn[:, 0], 0.0)             # t0 reference
    # z-index 1 has 4 cells * 0.5 = 2.0 at t1, 4.0 at t2
    assert np.isclose(dmap.dn[1, 1], 2.0)
    assert np.isclose(dmap.dn[1, 2], 4.0)
    assert np.allclose(dmap.dn[0], 0.0)               # untouched slab stays 0


def test_delta_map_extent_matches_axes():
    z = np.array([-3.0, 0.0, 3.0]); t = np.array([0.0, 1.0])
    cubes = [np.zeros((1, 1, 3)), np.ones((1, 1, 3))]
    dmap = planar_delta_map(cubes, t, z)
    assert dmap.extent == (0.0, 1.0, -3.0, 3.0)


def test_bad_shapes_raise():
    with pytest.raises(ValueError):
        planar_profile(np.ones((2, 2)), axis=1)        # not 3-D
    with pytest.raises(ValueError):
        planar_delta_map([], [], [])                    # no cubes
