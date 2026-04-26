"""
Verify that the four .vti files produced by smoketest.cpp round-trip
correctly: dimensions, origin, spacing, and per-point values.

Reads:
  out_real_ascii.vti       (4x4x4, value(ix,iy,iz) = 100*ix + 10*iy + iz)
  out_real_binary.vti      (same field, binary format)
  out_complex_ascii.vti    (3x3x3, two arrays psi_real / psi_imag)
  out_complex_binary.vti   (same field, binary format)

Uses xml.etree + numpy + base64 directly so it does not depend on a system
VTK install. The point ordering check is the load-bearing test: VTK
ImageData stores PointData with x fastest, so point index N = ix + nx*(iy + ny*iz).
The expected real field at point N decodes to f(ix,iy,iz) = 100*ix + 10*iy + iz.
"""

from __future__ import annotations

import base64
import struct
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np


def parse_vti(path: Path) -> dict:
    tree = ET.parse(path)
    root = tree.getroot()
    if root.tag != "VTKFile":
        raise SystemExit(f"{path}: root tag is {root.tag}, not VTKFile")
    image_data = root.find("ImageData")
    extent = image_data.attrib["WholeExtent"].split()
    e = [int(x) for x in extent]
    nx = e[1] - e[0] + 1
    ny = e[3] - e[2] + 1
    nz = e[5] - e[4] + 1
    origin = tuple(float(x) for x in image_data.attrib["Origin"].split())
    spacing = tuple(float(x) for x in image_data.attrib["Spacing"].split())

    arrays = {}
    point_data = image_data.find("Piece").find("PointData")
    for da in point_data.findall("DataArray"):
        name = da.attrib["Name"]
        fmt = da.attrib["format"]
        text = da.text or ""
        if fmt == "ascii":
            values = np.array([float(t) for t in text.split()], dtype=np.float64)
        elif fmt == "binary":
            blob = base64.b64decode("".join(text.split()))
            (n_bytes,) = struct.unpack("<Q", blob[:8])
            if n_bytes != len(blob) - 8:
                raise SystemExit(
                    f"{path}/{name}: header says {n_bytes} bytes but blob has "
                    f"{len(blob) - 8}"
                )
            values = np.frombuffer(blob[8:], dtype="<f8").copy()
        else:
            raise SystemExit(f"{path}/{name}: unknown format {fmt}")

        if values.size != nx * ny * nz:
            raise SystemExit(
                f"{path}/{name}: got {values.size} values, expected "
                f"{nx*ny*nz}"
            )
        arrays[name] = values
    return {
        "nx": nx, "ny": ny, "nz": nz,
        "origin": origin, "spacing": spacing,
        "arrays": arrays,
    }


def expected_real_field_xfastest(nx: int, ny: int, nz: int) -> np.ndarray:
    """f(ix,iy,iz) = 100*ix + 10*iy + iz, x-fastest stream order."""
    out = np.empty(nx * ny * nz, dtype=np.float64)
    k = 0
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                out[k] = 100.0 * ix + 10.0 * iy + iz
                k += 1
    return out


def expected_complex_xfastest(nx: int, ny: int, nz: int):
    re = np.empty(nx * ny * nz, dtype=np.float64)
    im = np.empty_like(re)
    k = 0
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                re[k] = ix + iy + iz
                im[k] = ix - iy + iz
                k += 1
    return re, im


def assert_close(label: str, a: np.ndarray, b: np.ndarray, atol: float = 1e-12):
    if a.shape != b.shape:
        raise SystemExit(f"{label}: shape mismatch {a.shape} vs {b.shape}")
    diff = np.max(np.abs(a - b))
    if diff > atol:
        raise SystemExit(f"{label}: max abs diff {diff} exceeds tol {atol}")
    print(f"  {label}: max abs diff = {diff:.3e}")


def main() -> int:
    here = Path(__file__).resolve().parent

    print("Real field (4x4x4) — ASCII:")
    a = parse_vti(here / "out_real_ascii.vti")
    assert (a["nx"], a["ny"], a["nz"]) == (4, 4, 4)
    assert a["origin"] == (-1.5, -1.5, -1.5)
    assert a["spacing"] == (1.0, 1.0, 1.0)
    expected = expected_real_field_xfastest(4, 4, 4)
    assert_close("ascii values", a["arrays"]["scalar"], expected)

    print("Real field (4x4x4) — BINARY:")
    b = parse_vti(here / "out_real_binary.vti")
    assert (b["nx"], b["ny"], b["nz"]) == (4, 4, 4)
    assert_close("binary values", b["arrays"]["scalar"], expected)
    assert_close("ascii vs binary", a["arrays"]["scalar"], b["arrays"]["scalar"])

    print("Complex field (3x3x3) — ASCII:")
    c = parse_vti(here / "out_complex_ascii.vti")
    assert (c["nx"], c["ny"], c["nz"]) == (3, 3, 3)
    re_ref, im_ref = expected_complex_xfastest(3, 3, 3)
    assert_close("ascii psi_real", c["arrays"]["psi_real"], re_ref)
    assert_close("ascii psi_imag", c["arrays"]["psi_imag"], im_ref)

    print("Complex field (3x3x3) — BINARY:")
    d = parse_vti(here / "out_complex_binary.vti")
    assert_close("binary psi_real", d["arrays"]["psi_real"], re_ref)
    assert_close("binary psi_imag", d["arrays"]["psi_imag"], im_ref)
    assert_close("ascii vs binary psi_real",
                 c["arrays"]["psi_real"], d["arrays"]["psi_real"])
    assert_close("ascii vs binary psi_imag",
                 c["arrays"]["psi_imag"], d["arrays"]["psi_imag"])

    print("\nALL OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
