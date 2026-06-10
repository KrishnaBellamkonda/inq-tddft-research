# Handles the parsing of the field data values from the outputed files
# and stores them in the dataclasses defined in fields.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from .fields import ComplexField3D, FieldMeta, RealField3D


class DataError(RuntimeError):
    """
    Raised when inqview cannot load or validate simulation data.
    """


@dataclass
class FieldSeries:
    """
    Lightweight descriptor for a directory containing field metadata files.
    """

    root: Path
    files: list[Path]
    field_name: str
    kind: str = "raw_meta"

    def __len__(self) -> int:
        return len(self.files)

    def __iter__(self):
        return iter(self.files)


class SimulationData:
    """
    Minimal simulation-directory helper for field discovery/loading.
    """

    def __init__(self, root: str | Path):
        self.root = Path(root).expanduser().resolve()
        if not self.root.exists():
            raise DataError(f"Simulation root does not exist: {self.root}")

    def require(self, relative_path: str | Path) -> Path:
        path = self.root / relative_path
        if not path.exists():
            raise DataError(f"Required path does not exist: {path}")
        return path

    def field_series(self, relative_dir: str | Path) -> FieldSeries:
        directory = self.require(relative_dir)
        if not directory.is_dir():
            raise DataError(f"Field series path is not a directory: {directory}")

        meta_files = sorted(directory.glob("*.meta.txt"), key=_natural_sort_key)
        if not meta_files:
            raise DataError(f"No '*.meta.txt' files found in: {directory}")

        first_meta = load_meta(meta_files[0])
        field_name = first_meta.field_name or meta_files[0].stem.replace(".meta", "")

        return FieldSeries(
            root=directory,
            files=meta_files,
            field_name=field_name,
            kind="raw_meta",
        )

    def load_real_field(
        self,
        raw_relative_path: str | Path | None = None,
        meta_relative_path: str | Path | None = None,
    ) -> RealField3D:
        raw_path = self.require(raw_relative_path) if raw_relative_path is not None else None
        meta_path = self.require(meta_relative_path) if meta_relative_path is not None else None
        return load_real_field(raw_path=raw_path, meta_path=meta_path)

    def load_complex_field(
        self,
        real_relative_path: str | Path | None = None,
        imag_relative_path: str | Path | None = None,
        meta_relative_path: str | Path | None = None,
    ) -> ComplexField3D:
        real_path = self.require(real_relative_path) if real_relative_path is not None else None
        imag_path = self.require(imag_relative_path) if imag_relative_path is not None else None
        meta_path = self.require(meta_relative_path) if meta_relative_path is not None else None
        return load_complex_field(
            real_raw_path=real_path,
            imag_raw_path=imag_path,
            meta_path=meta_path,
        )


def infer_meta_path(raw_path: str | Path) -> Path:
    """
    Infer metadata sidecar path from a raw data filename.

    Examples
    --------
    density_total.raw -> density_total.meta.txt
    orbital_0000_real.raw -> orbital_0000.meta.txt
    orbital_0000_imag.raw -> orbital_0000.meta.txt
    """
    raw_path = Path(raw_path)

    if raw_path.suffix != ".raw":
        raise DataError(f"Expected a '.raw' file, got: {raw_path}")

    stem = raw_path.stem
    if stem.endswith("_real"):
        stem = stem[:-5]
    elif stem.endswith("_imag"):
        stem = stem[:-5]

    return raw_path.with_name(f"{stem}.meta.txt")


def load_meta(meta_path: str | Path) -> FieldMeta:
    """
    Load a field metadata sidecar.

    Supported keys
    --------------
    Required:
      - type
      - dtype
      - nx
      - ny
      - nz
      - spacing_bohr
      - layout

    Optional but strongly preferred:
      - origin_bohr

    Optional:
      - field_name
      - units
      - time_au
      - value_file
      - real_file
      - imag_file
    """
    meta_path = Path(meta_path)

    if not meta_path.exists():
        raise DataError(f"Metadata file does not exist: {meta_path}")

    data = _parse_meta_file(meta_path)

    required = ["type", "dtype", "nx", "ny", "nz", "spacing_bohr", "layout"]
    missing = [key for key in required if key not in data]
    if missing:
        raise DataError(f"Missing required metadata keys in {meta_path}: {missing}")

    nx = _parse_int(data["nx"], key="nx", meta_path=meta_path)
    ny = _parse_int(data["ny"], key="ny", meta_path=meta_path)
    nz = _parse_int(data["nz"], key="nz", meta_path=meta_path)

    if "origin_bohr" in data:
        origin_bohr = _parse_float_triplet(
            data["origin_bohr"], key="origin_bohr", meta_path=meta_path
        )
    else:
        origin_bohr = (0.0, 0.0, 0.0)

    spacing_bohr = _parse_float_triplet(
        data["spacing_bohr"], key="spacing_bohr", meta_path=meta_path
    )

    time_au = None
    if "time_au" in data:
        time_au = _parse_float(data["time_au"], key="time_au", meta_path=meta_path)

    known_keys = {
        "type",
        "dtype",
        "nx",
        "ny",
        "nz",
        "origin_bohr",
        "spacing_bohr",
        "layout",
        "time_au",
        "field_name",
        "units",
        "value_file",
        "real_file",
        "imag_file",
    }
    extra = {k: v for k, v in data.items() if k not in known_keys}

    return FieldMeta(
        field_type=data["type"].strip(),
        dtype=data["dtype"].strip(),
        nx=nx,
        ny=ny,
        nz=nz,
        origin_bohr=origin_bohr,
        spacing_bohr=spacing_bohr,
        layout=data["layout"].strip(),
        field_name=data.get("field_name", None),
        units=data.get("units", None),
        time_au=time_au,
        value_file=data.get("value_file", None),
        real_file=data.get("real_file", None),
        imag_file=data.get("imag_file", None),
        extra=extra,
    )


def load_real_field(
    raw_path: str | Path | None = None,
    meta_path: str | Path | None = None,
) -> RealField3D:
    """
    Load a real-valued 3D field.

    Usage patterns
    --------------
    1. load_real_field(raw_path="density_total.raw")
       -> infers metadata path

    2. load_real_field(meta_path="density_total.meta.txt")
       -> uses value_file from metadata if present

    3. load_real_field(raw_path="...", meta_path="...")
       -> explicit override
    """
    if raw_path is None and meta_path is None:
        raise DataError("load_real_field requires at least one of raw_path or meta_path")

    raw_path = Path(raw_path) if raw_path is not None else None
    meta_path = Path(meta_path) if meta_path is not None else None

    if meta_path is None and raw_path is not None:
        meta_path = infer_meta_path(raw_path)

    assert meta_path is not None
    meta = load_meta(meta_path)

    if not meta.is_real:
        raise DataError(
            f"Metadata type says {meta.field_type!r}, but load_real_field expects "
            f"'real_field_3d'. Metadata file: {meta_path}"
        )

    if raw_path is None:
        raw_path = _resolve_real_value_path_from_meta(meta_path, meta)

    flat = _read_flat_real_array(raw_path, meta)
    array = _reshape_flat_array(flat, meta)

    return RealField3D(meta=meta, array=array)


def load_complex_field(
    real_raw_path: str | Path | None = None,
    imag_raw_path: str | Path | None = None,
    meta_path: str | Path | None = None,
) -> ComplexField3D:
    """
    Load a complex-valued 3D field from separate real and imaginary raw files.

    Usage patterns
    --------------
    1. load_complex_field(real_raw_path="orbital_0000_real.raw")
       -> infers imag path and metadata path

    2. load_complex_field(meta_path="orbital_0000.meta.txt")
       -> uses real_file and imag_file from metadata if present

    3. load_complex_field(real_raw_path="...", imag_raw_path="...", meta_path="...")
       -> explicit override
    """
    if real_raw_path is None and meta_path is None:
        raise DataError("load_complex_field requires at least one of real_raw_path or meta_path")

    real_raw_path = Path(real_raw_path) if real_raw_path is not None else None
    imag_raw_path = Path(imag_raw_path) if imag_raw_path is not None else None
    meta_path = Path(meta_path) if meta_path is not None else None

    if meta_path is None and real_raw_path is not None:
        meta_path = infer_meta_path(real_raw_path)

    assert meta_path is not None
    meta = load_meta(meta_path)

    if not meta.is_complex:
        raise DataError(
            f"Metadata type says {meta.field_type!r}, but load_complex_field expects "
            f"'complex_field_3d'. Metadata file: {meta_path}"
        )

    if real_raw_path is None:
        real_raw_path = _resolve_complex_part_path_from_meta(meta_path, meta, which="real")

    if imag_raw_path is None:
        if meta.imag_file:
            imag_raw_path = _resolve_complex_part_path_from_meta(meta_path, meta, which="imag")
        else:
            imag_raw_path = _infer_imag_path(real_raw_path)

    real_flat = _read_flat_real_array(real_raw_path, meta)
    imag_flat = _read_flat_real_array(imag_raw_path, meta)

    real = _reshape_flat_array(real_flat, meta)
    imag = _reshape_flat_array(imag_flat, meta)

    return ComplexField3D(meta=meta, real=real, imag=imag)


def _parse_meta_file(meta_path: Path) -> dict[str, str]:
    result: dict[str, str] = {}

    with meta_path.open("r", encoding="utf-8") as handle:
        for lineno, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line:
                continue
            if line.startswith("#"):
                continue

            if "=" not in line:
                raise DataError(
                    f"Invalid metadata line in {meta_path} at line {lineno}: {raw_line.rstrip()!r}"
                )

            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            if not key:
                raise DataError(
                    f"Empty metadata key in {meta_path} at line {lineno}: {raw_line.rstrip()!r}"
                )

            result[key] = value

    return result


def _parse_int(text: str, *, key: str, meta_path: Path) -> int:
    try:
        return int(text)
    except ValueError as exc:
        raise DataError(
            f"Could not parse integer for key {key!r} in {meta_path}: value={text!r}"
        ) from exc


def _parse_float(text: str, *, key: str, meta_path: Path) -> float:
    try:
        return float(text)
    except ValueError as exc:
        raise DataError(
            f"Could not parse float for key {key!r} in {meta_path}: value={text!r}"
        ) from exc


def _parse_float_triplet(text: str, *, key: str, meta_path: Path) -> tuple[float, float, float]:
    parts = text.split()
    if len(parts) != 3:
        raise DataError(
            f"Expected three values for key {key!r} in {meta_path}, got {len(parts)}: {text!r}"
        )

    try:
        values = tuple(float(x) for x in parts)
    except ValueError as exc:
        raise DataError(
            f"Could not parse float triplet for key {key!r} in {meta_path}: {text!r}"
        ) from exc

    return values  # type: ignore[return-value]


def _read_flat_real_array(raw_path: Path, meta: FieldMeta) -> np.ndarray:
    if not raw_path.exists():
        raise DataError(f"Raw data file does not exist: {raw_path}")

    actual_bytes = raw_path.stat().st_size
    expected_bytes = meta.expected_real_bytes

    if actual_bytes != expected_bytes:
        raise DataError(
            "Raw file size does not match metadata.\n"
            f"  file: {raw_path}\n"
            f"  expected bytes: {expected_bytes}\n"
            f"  actual bytes:   {actual_bytes}\n"
            f"  shape:          {meta.shape}\n"
            f"  dtype:          {meta.dtype}"
        )

    flat = np.fromfile(raw_path, dtype=meta.numpy_dtype)

    if flat.size != meta.num_points:
        raise DataError(
            f"Unexpected number of values in {raw_path}: got {flat.size}, "
            f"expected {meta.num_points}"
        )

    return flat


def _reshape_flat_array(flat: np.ndarray, meta: FieldMeta) -> np.ndarray:
    """
    Reshape flat data into (nx, ny, nz).

    Supported layouts
    -----------------
    x_slowest_z_fastest:
        flat index = ix * (ny * nz) + iy * nz + iz
        So C-order reshape to (nx, ny, nz) is correct.

    x_fastest_z_slowest:
        flat index = iz * (ny * nx) + iy * nx + ix
        So reshape as (nz, ny, nx) then transpose.
    """
    layout = meta.layout.strip()

    if layout in {"x_slowest_z_fastest", "z_fastest"}:
        return flat.reshape(meta.shape, order="C")

    if layout in {"x_fastest_z_slowest", "x_fastest"}:
        return flat.reshape((meta.nz, meta.ny, meta.nx), order="C").transpose(2, 1, 0)

    raise DataError(
        f"Unsupported field layout {layout!r}. "
        "Supported layouts are: 'x_slowest_z_fastest', 'z_fastest', "
        "'x_fastest_z_slowest', 'x_fastest'."
    )


def _resolve_real_value_path_from_meta(meta_path: Path, meta: FieldMeta) -> Path:
    if meta.value_file:
        return (meta_path.parent / meta.value_file).resolve()

    fallback = meta_path.with_suffix("")
    fallback = fallback.with_name(fallback.name.replace(".meta", "") + ".raw")
    return fallback.resolve()


def _resolve_complex_part_path_from_meta(meta_path: Path, meta: FieldMeta, *, which: str) -> Path:
    if which == "real":
        if meta.real_file:
            return (meta_path.parent / meta.real_file).resolve()
        fallback = meta_path.with_name(meta_path.name.replace(".meta.txt", "_real.raw"))
        return fallback.resolve()

    if which == "imag":
        if meta.imag_file:
            return (meta_path.parent / meta.imag_file).resolve()
        fallback = meta_path.with_name(meta_path.name.replace(".meta.txt", "_imag.raw"))
        return fallback.resolve()

    raise ValueError(f"Invalid complex part selector: {which!r}")


def _infer_imag_path(real_raw_path: Path) -> Path:
    name = real_raw_path.name
    if name.endswith("_real.raw"):
        return real_raw_path.with_name(name.replace("_real.raw", "_imag.raw"))
    raise DataError(
        "Could not infer imaginary raw path from real raw path. "
        f"Expected a filename ending in '_real.raw', got: {real_raw_path}"
    )


def _natural_sort_key(path: Path) -> list[object]:
    parts = re.split(r"(\d+)", path.name)
    key: list[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part)
    return key
