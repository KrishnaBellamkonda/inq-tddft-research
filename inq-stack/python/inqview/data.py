# This file handles the proper loading of the data
# from different file formats in which it is written.
#

from dataclasses import dataclass
from pathlib import Path
import numpy as np


class DataError(RuntimeError):
    pass


@dataclass
class ComplexField3D:
    meta_path: Path
    real_path: Path
    imag_path: Path
    field_name: str
    shape: tuple[int, int, int]
    spacing_bohr: tuple[float, float, float]
    kpoint: tuple[float, float, float]
    orbital_index: int
    spin_index: int
    values: np.ndarray  # complex128 array of shape (nx, ny, nz)


@dataclass
class RealField3D:
    meta_path: Path
    value_path: Path
    field_name: str
    shape: tuple[int, int, int]
    spacing_bohr: tuple[float, float, float]
    values: np.ndarray  # float64 array of shape (nx, ny, nz)


# General function for reading metadata of any type
def _read_key_value_meta(path: str | Path) -> dict[str, str]:
    path = Path(path)
    if not path.exists():
        raise DataError(f"Metadata file does not exist: {path}")

    meta: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        # Ignoring empty lines or lines with comments
        if not line or line.startswith("#"):
            continue
        # Filtering non assignment lines
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        meta[key.strip()] = value.strip()

    return meta


def load_complex_field(meta_path: str | Path) -> ComplexField3D:
    meta_path = Path(meta_path)
    meta = _read_key_value_meta(meta_path)

    # List of required properties
    required = [
        "type",
        "dtype",
        "field_name",
        "nx",
        "ny",
        "nz",
        "spacing_bohr",
        "layout",
        "real_file",
        "imag_file",
    ]

    # Find if any missing keys
    missing = [key for key in required if key not in meta]
    if missing:
        raise DataError(f"Missing required metadata keys in {meta_path}: {missing}")
    # File compatibility checks
    if meta["type"] != "complex_field_3d":
        raise DataError(f"Unsupported field type in {meta_path}: {meta['type']}")

    if meta["dtype"] != "float64":
        raise DataError(f"Unsupported dtype in {meta_path}: {meta['dtype']}")

    if meta["layout"] != "x_slowest_z_fastest":
        raise DataError(f"Unsupported layout in {meta_path}: {meta['layout']}")

    # Extracting data from meta, real and imag files

    # Grid dimensions
    nx = int(meta["nx"])
    ny = int(meta["ny"])
    nz = int(meta["nz"])

    # Grid spacing
    spacing = tuple(float(x) for x in meta["spacing_bohr"].split())
    if len(spacing) != 3:
        raise DataError(f"spacing_bohr must have 3 entries in {meta_path}")

    real_path = meta_path.parent / meta["real_file"]
    imag_path = meta_path.parent / meta["imag_file"]

    if not real_path.exists():
        raise DataError(f"Real raw file does not exist: {real_path}")
    if not imag_path.exists():
        raise DataError(f"Imag raw file does not exist: {imag_path}")

    # The real and imag files contain C++ arrays in
    # binary format. These are read directly into
    # nd arrays using the following function
    real = np.fromfile(real_path, dtype=np.float64)
    imag = np.fromfile(imag_path, dtype=np.float64)

    # Sanity checks for the size of the array
    expected_size = nx * ny * nz
    if real.size != expected_size:
        raise DataError(
            f"Real raw file size mismatch: expected {expected_size}, got {real.size}"
        )
    if imag.size != expected_size:
        raise DataError(
            f"Imag raw file size mismatch: expected {expected_size}, got {imag.size}"
        )

    real = real.reshape((nx, ny, nz))
    imag = imag.reshape((nx, ny, nz))
    values = real + 1j * imag

    # k point is obtained from the meta data
    # with a default value of 0 0 0 (Gamma)
    kpoint = tuple(float(x) for x in meta.get("kpoint", "0 0 0").split())
    if len(kpoint) != 3:
        # kpoint = (0.0, 0.0, 0.0)
        raise DataError(f"kpoint object has ${len(kpoint)} entries instead of 3.")

    orbital_index = int(meta.get("orbital_index", -1))
    spin_index = int(meta.get("spin_index", 0))

    return ComplexField3D(
        meta_path=meta_path,
        real_path=real_path,
        imag_path=imag_path,
        field_name=meta["field_name"],
        shape=(nx, ny, nz),
        spacing_bohr=(spacing[0], spacing[1], spacing[2]),
        kpoint=(kpoint[0], kpoint[1], kpoint[2]),
        orbital_index=orbital_index,
        spin_index=spin_index,
        values=values,
    )


def load_real_field(meta_path: str | Path) -> RealField3D:
    meta_path = Path(meta_path)
    meta = _read_key_value_meta(meta_path)

    required = [
        "type",
        "dtype",
        "field_name",
        "nx",
        "ny",
        "nz",
        "spacing_bohr",
        "layout",
        "value_file",
    ]
    missing = [key for key in required if key not in meta]
    if missing:
        raise DataError(f"Missing required metadata keys in {meta_path}: {missing}")

    if meta["type"] != "real_field_3d":
        raise DataError(f"Unsupported field type in {meta_path}: {meta['type']}")

    if meta["dtype"] != "float64":
        raise DataError(f"Unsupported dtype in {meta_path}: {meta['dtype']}")

    if meta["layout"] != "x_slowest_z_fastest":
        raise DataError(f"Unsupported layout in {meta_path}: {meta['layout']}")

    nx = int(meta["nx"])
    ny = int(meta["ny"])
    nz = int(meta["nz"])

    spacing = tuple(float(x) for x in meta["spacing_bohr"].split())
    if len(spacing) != 3:
        raise DataError("spacing_bohr does not have 3 entries")

    value_path = meta_path.parent / meta["value_file"]
    if not value_path.exists():
        raise DataError(f"Raw field file does not exist: {value_path}")

    values = np.fromfile(value_path, dtype=np.float64)

    expected_size = nx * ny * nz
    if values.size != expected_size:
        raise DataError(
            f"Raw field size mismatch: expected {expected_size}, got {values.size}"
        )

    values = values.reshape((nx, ny, nz))

    return RealField3D(
        meta_path=meta_path,
        value_path=value_path,
        field_name=meta["field_name"],
        shape=(nx, ny, nz),
        spacing_bohr=(spacing[0], spacing[1], spacing[2]),
        values=values,
    )
