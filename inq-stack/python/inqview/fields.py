from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np


@dataclass(frozen=True)
class FieldMeta:
    """
    Metadata describing a 3D field written by the C++ side.

    Notes
    -----
    - `dtype` is the scalar dtype used in the raw file(s).
    - For complex fields, the real and imaginary parts are stored separately.
    """

    field_type: str
    dtype: str
    nx: int
    ny: int
    nz: int

    origin_bohr: tuple[float, float, float]
    spacing_bohr: tuple[float, float, float]
    layout: str

    field_name: str | None = None
    units: str | None = None
    time_au: float | None = None

    # Optional file references from metadata
    value_file: str | None = None
    real_file: str | None = None
    imag_file: str | None = None

    extra: Mapping[str, str] = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, int, int]:
        return (self.nx, self.ny, self.nz)

    @property
    def num_points(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def numpy_dtype(self) -> np.dtype:
        try:
            return np.dtype(self.dtype)
        except TypeError as exc:
            raise ValueError(f"Unsupported dtype in metadata: {self.dtype!r}") from exc

    @property
    def voxel_volume_bohr3(self) -> float:
        dx, dy, dz = self.spacing_bohr
        return dx * dy * dz

    @property
    def expected_real_bytes(self) -> int:
        return self.num_points * self.numpy_dtype.itemsize

    @property
    def is_real(self) -> bool:
        return self.field_type == "real_field_3d"

    @property
    def is_complex(self) -> bool:
        return self.field_type == "complex_field_3d"


@dataclass
class RealField3D:
    """
    A real-valued 3D scalar field.
    """

    meta: FieldMeta
    array: np.ndarray

    def __post_init__(self) -> None:
        if self.array.ndim != 3:
            raise ValueError(f"RealField3D expects a 3D array, got ndim={self.array.ndim}")
        if self.array.shape != self.meta.shape:
            raise ValueError(
                f"Array shape {self.array.shape} does not match metadata shape {self.meta.shape}"
            )
        if not np.issubdtype(self.array.dtype, np.floating):
            raise ValueError(
                f"RealField3D expects a floating-point array, got dtype={self.array.dtype}"
            )

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.array.shape

    @property
    def dtype(self) -> np.dtype:
        return self.array.dtype

    @property
    def min(self) -> float:
        return float(np.min(self.array))

    @property
    def max(self) -> float:
        return float(np.max(self.array))

    @property
    def mean(self) -> float:
        return float(np.mean(self.array))


@dataclass
class ComplexField3D:
    """
    A complex-valued 3D field stored as separate real and imaginary arrays.
    """

    meta: FieldMeta
    real: np.ndarray
    imag: np.ndarray

    def __post_init__(self) -> None:
        if self.real.ndim != 3 or self.imag.ndim != 3:
            raise ValueError("ComplexField3D expects 3D arrays for both real and imag parts")
        if self.real.shape != self.meta.shape or self.imag.shape != self.meta.shape:
            raise ValueError(
                "ComplexField3D array shapes do not match metadata shape "
                f"{self.meta.shape}: real={self.real.shape}, imag={self.imag.shape}"
            )
        if not np.issubdtype(self.real.dtype, np.floating):
            raise ValueError(f"Real part must be floating-point, got {self.real.dtype}")
        if not np.issubdtype(self.imag.dtype, np.floating):
            raise ValueError(f"Imag part must be floating-point, got {self.imag.dtype}")

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.real.shape

    @property
    def dtype(self) -> np.dtype:
        return np.result_type(self.real.dtype, self.imag.dtype, np.complex64)

    @property
    def array(self) -> np.ndarray:
        return self.real + 1j * self.imag

    @property
    def magnitude(self) -> np.ndarray:
        return np.sqrt(self.real * self.real + self.imag * self.imag)

    @property
    def phase(self) -> np.ndarray:
        return np.arctan2(self.imag, self.real)
