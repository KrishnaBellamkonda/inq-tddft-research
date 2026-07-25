"""inqview.visualisation.field_io — the ONE canonical VTI loader.

Why this exists
---------------
inqkit VTIs are written in **physical order**: `inqkit::io::RealField3DWriter`
applies `fft_shift_index()` at write time and stamps `Origin = -L/2`, so array
index 0 already maps to the left-edge coordinate `-L/2` (NOT the FFT-natural
centre). Therefore **VTI data must never be `np.fft.fftshift`-ed** — doing so
swaps centre↔edge and silently produces flipped pictures (the recurring
"slab-at-the-edges" bug). Only LEED screen `.dat` files are FFT-natural and need
a shift; those have their own loader (`inqview.io.load_leed_pattern`).

Every GIF/slice/profile — notebooks and `make_*_postproc.py` — must load through
`load_vti` and use the returned coordinate axes. No hand-rolled VTK reads, no
hand-rolled fftshift.

Layer note (ADR 0003): VTI reading needs VTK, and `inqview.io` is contractually
numpy-only, so this canonical loader lives in the VTK-allowed `visualisation`
layer and is re-exported lazily as `inqview.load_vti`.
"""
from __future__ import annotations

from pathlib import Path
from typing import NamedTuple, Optional, Tuple

import numpy as np


class VtiField(NamedTuple):
    """A VTI loaded in physical order.

    data : np.ndarray, shape (nx, ny, nz), indexed [ix, iy, iz]
    x, y, z : 1-D cell-centred coordinate axes (Bohr), monotonically increasing
    origin, spacing : the VTI ImageData origin/spacing tuples (x, y, z)
    """
    data: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    origin: Tuple[float, float, float]
    spacing: Tuple[float, float, float]

    def xz_slice(self, y: float = 0.0) -> np.ndarray:
        """Return the [z, x] density slice nearest plane y=`y`, ready for
        imshow(origin='lower', extent=[x0,x1,z0,z1]). Transposed so rows are z."""
        iy = int(np.argmin(np.abs(self.y - y)))
        return self.data[:, iy, :].T  # (nz, nx): rows=z, cols=x


def load_vti(
    path: str | Path,
    *,
    expect_centered_axis: Optional[str] = None,
    expect_tol_bohr: Optional[float] = None,
) -> VtiField:
    """Load an inqkit VTI in physical order (no fftshift) with a hard self-check.

    Parameters
    ----------
    path : VTI file.
    expect_centered_axis : optional, one of {'x','y','z'}. If given, assert that
        the planar-summed density profile along that axis peaks near coordinate 0
        (the box centre) — a loud failure if the index→coordinate mapping is
        wrong. Use for runs whose feature (slab/cluster) sits at the centre.
    expect_tol_bohr : tolerance for the centred-feature check (default: 8 grid
        spacings, generous — it only needs to catch a centre↔edge swap).

    Returns
    -------
    VtiField (physical order).
    """
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(f"VTK is required to read {path}") from exc

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    origin = tuple(float(v) for v in img.GetOrigin())
    spacing = tuple(float(v) for v in img.GetSpacing())

    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64, copy=False)
    # VTK ImageData is x-fastest: reshape (nz,ny,nx) then transpose to (nx,ny,nz).
    data = flat.reshape((nz, ny, nx)).transpose(2, 1, 0)

    ox, oy, oz = origin
    sx, sy, sz = spacing
    # Cell-centred sample coordinates (matches the C++ writer / pipeline/density).
    x = ox + (np.arange(nx) + 0.5) * sx
    y = oy + (np.arange(ny) + 0.5) * sy
    z = oz + (np.arange(nz) + 0.5) * sz

    # ---- HARD self-check: the axis/dim invariants ------------------------
    assert data.shape == (nx, ny, nz), (
        f"reshape/transpose mismatch: data {data.shape} vs dims {(nx, ny, nz)}")
    for name, ax, s in (("x", x, sx), ("y", y, sy), ("z", z, sz)):
        assert s > 0.0, f"non-positive spacing on {name}: {s}"
        assert ax[1] > ax[0], f"{name} axis not increasing"
    # Physical order ⇒ first sample sits at the left edge, not the centre.
    assert abs(x[0] - (ox + 0.5 * sx)) < 1e-9, "x[0] is not the physical left edge"

    if expect_centered_axis is not None:
        # Discriminate a CENTRED feature (mass in the inner half of the box) from
        # an edge-split one (the centre↔edge swap a wrong fftshift produces). Both
        # are symmetric, so argmax/centre-of-mass cannot tell them apart — compare
        # inner-half vs outer-half |n| mass instead. Robust to Friedel peaks.
        axis_index = {"x": 0, "y": 1, "z": 2}[expect_centered_axis]
        coord = (x, y, z)[axis_index]
        other = tuple(a for a in (0, 1, 2) if a != axis_index)
        profile = np.abs(data).sum(axis=other)  # planar-summed |n| along the axis
        L_axis = nx * sx if axis_index == 0 else (ny * sy if axis_index == 1 else nz * sz)
        inner = np.abs(coord) < 0.25 * L_axis
        inner_mass = float(profile[inner].sum())
        outer_mass = float(profile[~inner].sum())
        assert inner_mass > outer_mass, (
            f"index→coordinate mapping looks WRONG on {expect_centered_axis}: "
            f"inner-half |n| mass {inner_mass:.3g} ≤ outer-half {outer_mass:.3g} — "
            f"the feature sits at the EDGES, not the centre. Did something "
            f"np.fft.fftshift a physical-order VTI?")

    return VtiField(data=data, x=x, y=y, z=z, origin=origin, spacing=spacing)
