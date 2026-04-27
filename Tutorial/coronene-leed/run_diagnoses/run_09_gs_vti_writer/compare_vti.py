"""
Diff every .vti emitted by run_09 (C++ native VTI writer) against the
matching .vti from run_08 (inqview Python conversion). Both runs use the
same SCF input and writer contract, so every (origin, spacing, dimensions,
per-point value) should match exactly.

The only legitimate source of difference is floating-point reproducibility
of the GS SCF — same code, same hardware, same compiler flags, so we
expect bit-identical output. A tolerance of 1e-12 absolute is generous.
"""

from __future__ import annotations

import base64
import struct
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

REF_DIR = Path("/local/data/public/skcb2/tddft/Tutorial/coronene-leed/"
               "run_diagnoses/run_08_gs_only_wp_check/results")
OUR_DIR = Path(__file__).resolve().parent / "results"


def read_vti(path: Path) -> dict:
    tree = ET.parse(path)
    root = tree.getroot()
    img = root.find("ImageData")
    e = [int(x) for x in img.attrib["WholeExtent"].split()]
    nx = e[1] - e[0] + 1
    ny = e[3] - e[2] + 1
    nz = e[5] - e[4] + 1
    origin = tuple(float(x) for x in img.attrib["Origin"].split())
    spacing = tuple(float(x) for x in img.attrib["Spacing"].split())
    arrays = {}
    for da in img.find("Piece").find("PointData").findall("DataArray"):
        name = da.attrib["Name"]
        fmt = da.attrib["format"]
        text = da.text or ""
        if fmt == "ascii":
            values = np.fromstring(text, sep=" ", dtype=np.float64)
        elif fmt == "binary":
            blob = base64.b64decode("".join(text.split()))
            (n_bytes,) = struct.unpack("<Q", blob[:8])
            values = np.frombuffer(blob[8:8 + n_bytes], dtype="<f8").copy()
        else:
            raise SystemExit(f"{path}: unknown DataArray format {fmt}")
        arrays[name] = values
    return {
        "shape": (nx, ny, nz),
        "origin": origin,
        "spacing": spacing,
        "arrays": arrays,
    }


def compare_one(label: str, our_path: Path, ref_path: Path,
                tol: float = 1e-12) -> bool:
    if not ref_path.exists():
        print(f"  SKIP {label}: no reference at {ref_path}")
        return True
    a = read_vti(our_path)
    b = read_vti(ref_path)
    ok = True
    if a["shape"] != b["shape"]:
        print(f"  FAIL {label}: shape {a['shape']} vs {b['shape']}")
        ok = False
    if not np.allclose(a["origin"], b["origin"], atol=tol):
        print(f"  FAIL {label}: origin {a['origin']} vs {b['origin']}")
        ok = False
    if not np.allclose(a["spacing"], b["spacing"], atol=tol):
        print(f"  FAIL {label}: spacing {a['spacing']} vs {b['spacing']}")
        ok = False
    common = set(a["arrays"]).intersection(b["arrays"])
    if not common:
        print(f"  FAIL {label}: no shared array names "
              f"({list(a['arrays'])} vs {list(b['arrays'])})")
        return False
    for name in sorted(common):
        diff = float(np.max(np.abs(a["arrays"][name] - b["arrays"][name])))
        status = "OK  " if diff <= tol else "FAIL"
        if diff > tol:
            ok = False
        print(f"  {status} {label} [{name}]  max|Δ| = {diff:.3e}")
    return ok


def main() -> int:
    if not OUR_DIR.exists():
        print(f"FATAL: {OUR_DIR} does not exist; run run_09 first.",
              file=sys.stderr)
        return 2

    our_vtis = sorted(OUR_DIR.rglob("*.vti"))
    if not our_vtis:
        print("FATAL: no .vti files under results/", file=sys.stderr)
        return 2

    print(f"Comparing {len(our_vtis)} .vti files against {REF_DIR}\n")
    all_ok = True
    for our_path in our_vtis:
        rel = our_path.relative_to(OUR_DIR)
        ref_path = REF_DIR / rel
        ok = compare_one(str(rel), our_path, ref_path)
        all_ok = all_ok and ok

    print("\n" + ("ALL OK" if all_ok else "SOME FILES DIFFER"))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
