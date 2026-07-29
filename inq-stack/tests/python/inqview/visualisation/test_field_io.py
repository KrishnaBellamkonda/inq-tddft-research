"""Tests for the canonical VTI loader `inqview.visualisation.load_vti`.

Known-case checks (validation-gates rule): synthesise VTIs with VTK whose
geometry is known, then assert the loader returns PHYSICAL order (no fftshift),
correct cell-centred axes, and that the `expect_centered_axis` self-check passes
for a centred feature and FIRES for an edge-split one (the centre↔edge swap a
wrong fftshift produces).
"""
from __future__ import annotations

import numpy as np
import pytest

vtk = pytest.importorskip("vtk")
from vtk.util.numpy_support import numpy_to_vtk  # noqa: E402

from inqview.visualisation import load_vti  # noqa: E402


def _write_vti(path, data_xyz, origin, spacing):
    """Write data[ix,iy,iz] to a .vti with the given origin/spacing (physical)."""
    nx, ny, nz = data_xyz.shape
    img = vtk.vtkImageData()
    img.SetDimensions(nx, ny, nz)
    img.SetOrigin(*origin)
    img.SetSpacing(*spacing)
    # VTK is x-fastest: flatten as [iz,iy,ix] then ravel.
    flat = np.ascontiguousarray(data_xyz.transpose(2, 1, 0).ravel(), dtype=np.float64)
    arr = numpy_to_vtk(flat, deep=True)
    arr.SetName("density")
    img.GetPointData().SetScalars(arr)
    w = vtk.vtkXMLImageDataWriter()
    w.SetFileName(str(path))
    w.SetInputData(img)
    w.Write()


def _grid(L=50.0, dx=0.5):
    n = int(round(L / dx))
    origin = (-L / 2, -L / 2, -L / 2)
    spacing = (dx, dx, dx)
    z = origin[2] + (np.arange(n) + 0.5) * dx
    return n, origin, spacing, z


def test_load_vti_physical_order_and_axes(tmp_path):
    n, origin, spacing, z = _grid()
    data = np.zeros((n, n, n))
    # centred slab: nonzero where |z| < 12.5
    data[:, :, np.abs(z) < 12.5] = 1.0
    p = tmp_path / "slab.vti"
    _write_vti(p, data, origin, spacing)

    f = load_vti(p, expect_centered_axis="z")  # must NOT raise
    assert f.data.shape == (n, n, n)
    # physical order: first sample is the left edge, axes increasing & cell-centred
    assert abs(f.x[0] - (origin[0] + 0.5 * spacing[0])) < 1e-9
    assert f.z[1] > f.z[0]
    np.testing.assert_allclose(f.z, z)
    # round-trips the data unchanged (no fftshift applied)
    np.testing.assert_allclose(f.data, data)


def test_load_vti_does_not_fftshift(tmp_path):
    """A feature at a known asymmetric coordinate must come back at THAT
    coordinate — i.e. the loader applies no shift."""
    n, origin, spacing, z = _grid()
    data = np.zeros((n, n, n))
    iz = int(np.argmin(np.abs(z - 8.0)))  # bright plane at z≈+8
    data[:, :, iz] = 1.0
    p = tmp_path / "plane.vti"
    _write_vti(p, data, origin, spacing)

    f = load_vti(p)
    prof = f.data.sum(axis=(0, 1))
    assert abs(f.z[int(prof.argmax())] - 8.0) < spacing[2]  # stays at +8, not flipped


def test_centered_guard_fires_on_edge_split(tmp_path):
    """An edge-split feature (what a wrong fftshift makes) must trip the guard."""
    n, origin, spacing, z = _grid()
    data = np.zeros((n, n, n))
    data[:, :, np.abs(z) > 0.375 * 50.0] = 1.0  # mass in the outer eighths
    p = tmp_path / "edges.vti"
    _write_vti(p, data, origin, spacing)

    with pytest.raises(AssertionError, match="EDGES"):
        load_vti(p, expect_centered_axis="z")


if __name__ == "__main__":
    import tempfile, pathlib
    d = pathlib.Path(tempfile.mkdtemp())
    test_load_vti_physical_order_and_axes(d)
    test_load_vti_does_not_fftshift(d)
    test_centered_guard_fires_on_edge_split(d)
    print("load_vti tests passed.")
