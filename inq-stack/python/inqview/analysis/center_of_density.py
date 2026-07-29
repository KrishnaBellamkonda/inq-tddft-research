"""Centre-of-density in Python, correct node convention (IV-M02).

Computes ⟨r⟩ = ∫ r ρ / ∫ ρ from a saved density grid using INQ's NODE
convention (r_i = origin + i·dx), NOT the half-cell (i+½)·dx convention that the
inqkit C++ ``center_of_density`` still uses (the deferred E04 bug). Because
inqview is post-processing, it can compute COD correctly without waiting on the
production-side fix — and the difference between the two conventions is exactly
dx/2, which a test pins down to document E04.

Provides the WP / total / bath three-way comparison (bath = total − wp, the
canonical run-independent bath density). Pure numpy.

Density arrays use the inqkit x-slowest / z-fastest layout: shape (nx, ny, nz).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

Triple = "tuple[float, float, float]"


@dataclass(frozen=True)
class COD:
    x: float
    y: float
    z: float
    total_weight: float


@dataclass(frozen=True)
class CODComparison:
    wp: COD
    total: COD
    bath: COD


def center_of_density(rho, origin, spacing, *, half_cell: bool = False) -> COD:
    """⟨r⟩ of a density grid (Bohr). ``half_cell=True`` reproduces the E04
    (i+½)·dx convention — used only to demonstrate the dx/2 offset."""
    rho = np.asarray(rho, dtype=float)
    nx, ny, nz = rho.shape
    off = 0.5 if half_cell else 0.0
    x = origin[0] + (np.arange(nx) + off) * spacing[0]
    y = origin[1] + (np.arange(ny) + off) * spacing[1]
    z = origin[2] + (np.arange(nz) + off) * spacing[2]
    w = float(rho.sum())
    if w <= 0.0:
        return COD(0.0, 0.0, 0.0, 0.0)
    cx = float((rho.sum(axis=(1, 2)) * x).sum() / w)
    cy = float((rho.sum(axis=(0, 2)) * y).sum() / w)
    cz = float((rho.sum(axis=(0, 1)) * z).sum() / w)
    return COD(cx, cy, cz, w)


def compare(total_rho, wp_rho, origin, spacing) -> CODComparison:
    """WP / total / bath centres of density. bath = total − wp (canonical)."""
    total_rho = np.asarray(total_rho, dtype=float)
    wp_rho = np.asarray(wp_rho, dtype=float)
    bath_rho = total_rho - wp_rho
    return CODComparison(
        wp=center_of_density(wp_rho, origin, spacing),
        total=center_of_density(total_rho, origin, spacing),
        bath=center_of_density(bath_rho, origin, spacing),
    )
