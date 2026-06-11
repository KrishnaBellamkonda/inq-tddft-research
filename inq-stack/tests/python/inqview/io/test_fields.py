"""Unit tests for the field dataclasses (``inqview.fields``).

Pure (numpy only). Checks derived geometry and the post-init validation guards
that reject malformed arrays — these are the structures every loader returns,
so their invariants protect the whole io layer.

Current path ``inqview.fields`` → moves to ``inqview.io`` in the restructure.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.io.fields import ComplexField3D, FieldMeta, RealField3D

pytestmark = pytest.mark.io


def _meta(field_type="real_field_3d", dtype="float64"):
    return FieldMeta(
        field_type=field_type, dtype=dtype, nx=2, ny=3, nz=4,
        origin_bohr=(0.0, 0.0, 0.0), spacing_bohr=(0.5, 0.25, 2.0),
        layout="x_fastest",
    )


def test_meta_derived_geometry():
    m = _meta()
    assert m.shape == (2, 3, 4)
    assert m.num_points == 24
    assert m.voxel_volume_bohr3 == pytest.approx(0.5 * 0.25 * 2.0)
    assert m.numpy_dtype == np.dtype("float64")
    assert m.expected_real_bytes == 24 * 8        # float64 = 8 bytes
    assert m.is_real and not m.is_complex


def test_meta_bad_dtype_raises():
    with pytest.raises(ValueError):
        _meta(dtype="not_a_dtype").numpy_dtype


def test_realfield_accepts_matching_array():
    m = _meta()
    rf = RealField3D(meta=m, array=np.zeros(m.shape, dtype=np.float64))
    assert rf.shape == (2, 3, 4)
    assert rf.min == 0.0 and rf.max == 0.0 and rf.mean == 0.0


def test_realfield_rejects_shape_mismatch():
    m = _meta()
    with pytest.raises(ValueError):
        RealField3D(meta=m, array=np.zeros((2, 2, 2), dtype=np.float64))


def test_realfield_rejects_non_float():
    m = _meta()
    with pytest.raises(ValueError):
        RealField3D(meta=m, array=np.zeros(m.shape, dtype=np.int32))


def test_complexfield_magnitude_phase_array():
    m = _meta(field_type="complex_field_3d")
    real = np.full(m.shape, 3.0)
    imag = np.full(m.shape, 4.0)
    cf = ComplexField3D(meta=m, real=real, imag=imag)
    assert np.allclose(cf.magnitude, 5.0)             # sqrt(3²+4²)
    assert np.allclose(cf.phase, np.arctan2(4.0, 3.0))
    assert np.allclose(cf.array, 3.0 + 4.0j)


def test_complexfield_rejects_shape_mismatch():
    m = _meta(field_type="complex_field_3d")
    with pytest.raises(ValueError):
        ComplexField3D(meta=m, real=np.zeros(m.shape), imag=np.zeros((1, 1, 1)))


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
