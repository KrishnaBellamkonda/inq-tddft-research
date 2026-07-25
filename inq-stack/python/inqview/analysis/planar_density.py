"""planar_density.py — planar-integrated density profiles and Δn(z, t).

The paper's Fig. 1 trace: the *planar-integrated* electron-density change as a
function of the propagation coordinate ``z`` and time ``t``,

    Δn(z, t) = ∫∫ [ n(x, y, z, t) − n(x, y, z, t₀) ] dx dy
             = (dx·dy) · ∑_{x,y} [ n(·, t) − n(·, t₀) ] .

This kernel operates on **already-loaded** density cubes (numpy), so it stays in
the deps-clean ``analysis`` layer; VTI loading (which needs VTK) is the caller's
job in the ``visualisation`` / pipeline layer.

Sign convention: Δn > 0 = electron accumulation (pile-up), Δn < 0 = depletion.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["PlanarDeltaN", "planar_profile", "planar_delta_map"]


def planar_profile(cube: np.ndarray, axis: int = 2,
                   cell_area: float = 1.0) -> np.ndarray:
    """Planar-integrate a 3-D density cube onto ``axis``.

    cube      : (nx, ny, nz) density.
    axis      : the kept axis (default 2 = z); the other two are summed.
    cell_area : dx·dy [bohr²] so the sum becomes an integral (default 1 = raw
                sum). Pass dx*dy to get ∫∫ n dx dy.

    Returns a 1-D profile of length ``cube.shape[axis]``.
    """
    cube = np.asarray(cube, dtype=np.float64)
    if cube.ndim != 3:
        raise ValueError(f"cube must be 3-D; got shape {cube.shape}")
    sum_axes = tuple(a for a in range(3) if a != axis)
    return cube.sum(axis=sum_axes) * cell_area


@dataclass(frozen=True)
class PlanarDeltaN:
    """Δn(z, t) map plus its axes.

    dn    : (nz, nt) array — Δn integrated in the transverse plane, per z, per
            frame (column t is frame t minus frame 0).
    z     : (nz,) coordinate along the kept axis [bohr].
    t     : (nt,) time of each frame [a.u.].
    """

    dn: np.ndarray
    z: np.ndarray
    t: np.ndarray

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """(t_min, t_max, z_min, z_max) for imshow(origin='lower')."""
        return (float(self.t[0]), float(self.t[-1]),
                float(self.z[0]), float(self.z[-1]))


def planar_delta_map(cubes, times, z, axis: int = 2,
                     cell_area: float = 1.0) -> PlanarDeltaN:
    """Build Δn(z, t) from a time-ordered sequence of density cubes.

    cubes     : iterable of (nx, ny, nz) cubes, frame 0 is the reference t₀.
    times     : matching frame times [a.u.].
    z         : coordinate axis [bohr] for the kept ``axis`` (length nz).
    axis      : kept axis (default 2 = z).
    cell_area : dx·dy [bohr²] (default 1 = raw sum).

    Returns PlanarDeltaN with ``dn`` shape (nz, nt).
    """
    profiles = [planar_profile(c, axis=axis, cell_area=cell_area) for c in cubes]
    if not profiles:
        raise ValueError("no cubes provided")
    P = np.stack(profiles, axis=1)        # (nz, nt)
    dn = P - P[:, [0]]                     # subtract the t0 reference column
    t = np.asarray(times, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    if dn.shape != (z.size, t.size):
        raise ValueError(f"shape mismatch: dn {dn.shape} vs (z={z.size}, t={t.size})")
    return PlanarDeltaN(dn=dn, z=z, t=t)
