"""LEGACY — Python VTI writer (deprecated; ADR 0003 cleanup).

VTI emission is now done in C++ (`inqkit::io::VTIImageDataWriter`), so this
module is redundant for new runs. It is **retained, not cut**, because
`visualisation.paraview` (`convert_real_series_to_vti`) and the `defaults` movie
wrappers still call it to convert any raw+meta field series for ParaView
rendering. Do not build new features on it; prefer the C++ VTI output. Slated for
removal once the ParaView path is confirmed unused.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from xml.sax.saxutils import escape

import numpy as np

from ..io.data import FieldSeries, load_meta, load_real_field
from ..io.fields import RealField3D


@dataclass
class VTISeriesResult:
    root: Path
    files: list[Path]
    array_name: str
    data_min: float
    data_max: float


def write_vti(
    field: RealField3D,
    output_path: str | Path,
    array_name: str | None = None,
) -> Path:
    """
    Write a real 3D field to a VTK XML ImageData (.vti) file.

    Notes
    -----
    - This first implementation writes ASCII VTI for simplicity and robustness.
    - The field is written as POINT_DATA.
    - VTK ImageData expects x to vary fastest in point order, so we transpose
      from our in-memory (nx, ny, nz) layout before flattening.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    meta = field.meta
    array_name = array_name or meta.field_name or "density"

    nx, ny, nz = meta.shape
    ox, oy, oz = meta.origin_bohr
    dx, dy, dz = meta.spacing_bohr

    extent = f"0 {nx - 1} 0 {ny - 1} 0 {nz - 1}"
    vtk_type = _vtk_xml_scalar_type(field.array.dtype)

    values = _vtk_point_data_flatten(field.array)
    values_text = _format_ascii_values(values)

    xml = f"""<?xml version="1.0"?>
<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
  <ImageData WholeExtent="{extent}" Origin="{ox:.16g} {oy:.16g} {oz:.16g}" Spacing="{dx:.16g} {dy:.16g} {dz:.16g}">
    <Piece Extent="{extent}">
      <PointData Scalars="{escape(array_name)}">
        <DataArray type="{vtk_type}" Name="{escape(array_name)}" format="ascii">
{values_text}
        </DataArray>
      </PointData>
      <CellData>
      </CellData>
    </Piece>
  </ImageData>
</VTKFile>
"""

    output_path.write_text(xml, encoding="utf-8")
    return output_path


def convert_real_meta_to_vti(
    meta_path: str | Path,
    output_path: str | Path | None = None,
    array_name: str | None = None,
) -> Path:
    """
    Load one real-field metadata/raw pair and write one .vti file.
    """
    meta_path = Path(meta_path)
    field = load_real_field(meta_path=meta_path)

    if output_path is None:
        stem = _meta_stem(meta_path)
        output_path = meta_path.with_name(f"{stem}.vti")

    return write_vti(field=field, output_path=output_path, array_name=array_name)


def convert_real_series_to_vti(
    series: FieldSeries | Iterable[str | Path],
    output_dir: str | Path,
    array_name: str | None = None,
) -> VTISeriesResult:
    """
    Convert a series of real-field metadata files into a VTI series.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    meta_paths = _series_meta_paths(series)
    if not meta_paths:
        raise RuntimeError("No metadata files were provided for VTI conversion.")

    out_files: list[Path] = []
    global_min = np.inf
    global_max = -np.inf
    final_array_name: str | None = array_name

    for meta_path in meta_paths:
        field = load_real_field(meta_path=meta_path)

        current_array_name = array_name or field.meta.field_name or "density"
        if final_array_name is None:
            final_array_name = current_array_name

        stem = _meta_stem(meta_path)
        out_path = output_dir / f"{stem}.vti"
        write_vti(field=field, output_path=out_path, array_name=current_array_name)
        out_files.append(out_path)

        global_min = min(global_min, float(np.min(field.array)))
        global_max = max(global_max, float(np.max(field.array)))

    assert final_array_name is not None

    return VTISeriesResult(
        root=output_dir,
        files=out_files,
        array_name=final_array_name,
        data_min=float(global_min),
        data_max=float(global_max),
    )


def _series_meta_paths(series: FieldSeries | Iterable[str | Path]) -> list[Path]:
    if isinstance(series, FieldSeries):
        return [Path(p) for p in series.files]
    return [Path(p) for p in series]


def _meta_stem(meta_path: Path) -> str:
    name = meta_path.name
    if name.endswith(".meta.txt"):
        return name[:-9]
    if name.endswith(".txt"):
        return meta_path.stem
    return meta_path.name


def _vtk_xml_scalar_type(dtype: np.dtype) -> str:
    dtype = np.dtype(dtype)

    if dtype == np.dtype(np.float64):
        return "Float64"
    if dtype == np.dtype(np.float32):
        return "Float32"
    if dtype == np.dtype(np.int64):
        return "Int64"
    if dtype == np.dtype(np.int32):
        return "Int32"
    if dtype == np.dtype(np.uint64):
        return "UInt64"
    if dtype == np.dtype(np.uint32):
        return "UInt32"
    if dtype == np.dtype(np.uint16):
        return "UInt16"
    if dtype == np.dtype(np.uint8):
        return "UInt8"

    raise TypeError(f"Unsupported dtype for VTI writing: {dtype}")


def _vtk_point_data_flatten(array_xyz: np.ndarray) -> np.ndarray:
    """
    Convert an array shaped (nx, ny, nz) into VTK point ordering.

    Our in-memory convention:
      array[ix, iy, iz]

    VTK ImageData point order is x-fastest, then y, then z.

    So we transpose to (z, y, x) and flatten in C-order.
    """
    if array_xyz.ndim != 3:
        raise ValueError(f"Expected a 3D array, got shape={array_xyz.shape}")
    return np.asarray(array_xyz).transpose(2, 1, 0).ravel(order="C")


def _format_ascii_values(values: np.ndarray, values_per_line: int = 8) -> str:
    parts: list[str] = []
    n = values.size
    for start in range(0, n, values_per_line):
        chunk = values[start : start + values_per_line]
        parts.append("          " + " ".join(f"{float(x):.16g}" for x in chunk))
    return "\n".join(parts)
