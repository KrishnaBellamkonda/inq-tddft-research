from __future__ import annotations

"""
Diagnostic-run analysis for run_06_centred_writer_check.

This run is designed to test the hypothesis that inqkit's density writer
publishes an FFT-natural array (cell centre at array index 0) but tags it
with origin = -L/2, so the resulting .vti is rendered with the molecule at
the corners.

Outputs produced here:
  results/vti/                          standard (buggy) VTI conversion
                                        (molecule expected at corners)
  results/vti_fftshifted/               experimental: same .raw but with an
                                        np.fft.fftshift applied along x,y,z
                                        (molecule expected to render at the
                                        metadata centre, i.e. correctly)
  results/grid_diagnostics_summary.txt  textual summary of where the density
                                        peak sits in the raw INQ array (so
                                        the user can see at a glance whether
                                        FFT-natural ordering is in use)

If the buggy .vti renders the molecule at the corners and the fftshifted one
renders it at the centre, hypothesis H1 is confirmed and we know the fix
location: inqkit::fields::density::total/orbital must shift the array
before storing it.
"""

import re
import sys
from pathlib import Path

import numpy as np

RUN_DIR = Path(__file__).resolve().parent
RESULTS_DIR = RUN_DIR / "results"
VTI_ROOT = RESULTS_DIR / "vti"
VTI_FFTSHIFT_ROOT = RESULTS_DIR / "vti_fftshifted"

# run_diagnoses/run_06_*/ -> parents[3] is the repo root
REPO_ROOT = RUN_DIR.parents[3] if len(RUN_DIR.parents) >= 4 else RUN_DIR
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))

try:
    import inqview
    from inqview.data import load_real_field
    from inqview.fields import RealField3D, FieldMeta
    from inqview.vti import write_vti
except Exception as exc:
    raise SystemExit(
        "Could not import inqview. Make sure inq-stack/python is in the path "
        "or that inqview is installed in your environment.\n"
        f"Original import error: {exc}"
    )


META_SUFFIX = ".meta.txt"


def _natural_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return [int(p) if p.isdigit() else p for p in parts]


def _paired_raw_exists(meta_path: Path) -> bool:
    name = meta_path.name
    if not name.endswith(META_SUFFIX):
        return False
    raw_name = name[: -len(META_SUFFIX)] + ".raw"
    return (meta_path.parent / raw_name).exists()


def _find_series_roots(results_dir: Path) -> list[Path]:
    roots: set[Path] = set()
    for meta_path in results_dir.rglob(f"*{META_SUFFIX}"):
        if not _paired_raw_exists(meta_path):
            continue
        if VTI_ROOT in meta_path.parents:
            continue
        if VTI_FFTSHIFT_ROOT in meta_path.parents:
            continue
        roots.add(meta_path.parent)
    return sorted(roots)


def _meta_stem(meta_path: Path) -> str:
    name = meta_path.name
    if name.endswith(".meta.txt"):
        return name[:-9]
    return meta_path.stem


def _array_name_for(meta_path: Path) -> str:
    parent = meta_path.parent.name
    return "density"


def _convert_one(meta_path: Path, out_dir: Path, fftshift: bool) -> Path:
    field = load_real_field(meta_path=meta_path)
    arr = np.asarray(field.array)
    if fftshift:
        # Shift so that array index 0 corresponds to position -L/2 in the
        # rendered .vti, matching the metadata Origin field. This undoes
        # INQ's FFT-natural layout where index 0 is the cell centre.
        arr = np.fft.fftshift(arr, axes=(0, 1, 2))

    shifted_field = RealField3D(meta=field.meta, array=arr)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{_meta_stem(meta_path)}.vti"
    write_vti(field=shifted_field, output_path=out_path, array_name=_array_name_for(meta_path))
    return out_path


def _convert_series_pair(series_dir: Path) -> tuple[int, int]:
    metas = sorted(
        (p for p in series_dir.glob(f"*{META_SUFFIX}") if _paired_raw_exists(p)),
        key=_natural_key,
    )
    if not metas:
        return 0, 0

    rel = series_dir.relative_to(RESULTS_DIR)
    out_buggy = VTI_ROOT / rel
    out_fixed = VTI_FFTSHIFT_ROOT / rel

    n_buggy = 0
    n_fixed = 0
    for meta_path in metas:
        _convert_one(meta_path, out_buggy, fftshift=False)
        n_buggy += 1
        _convert_one(meta_path, out_fixed, fftshift=True)
        n_fixed += 1
    return n_buggy, n_fixed


def _summary_for_density(series_dir: Path) -> str | None:
    metas = sorted(
        (p for p in series_dir.glob(f"*{META_SUFFIX}") if _paired_raw_exists(p)),
        key=_natural_key,
    )
    if not metas:
        return None
    field = load_real_field(meta_path=metas[0])
    arr = np.asarray(field.array)
    nx, ny, nz = arr.shape
    ix, iy, iz = np.unravel_index(int(np.argmax(arr)), arr.shape)
    centre_val = float(arr[0, 0, 0])
    corner_val = float(arr[nx // 2, ny // 2, nz // 2])

    sym = lambda i, n: i if i < (n + 1) // 2 else i - n
    peak_x_fft = sym(int(ix), nx)
    peak_y_fft = sym(int(iy), ny)
    peak_z_fft = sym(int(iz), nz)

    return (
        f"  series       : {series_dir.relative_to(RESULTS_DIR)}\n"
        f"  shape        : ({nx}, {ny}, {nz})\n"
        f"  max array idx: ({ix}, {iy}, {iz})\n"
        f"  max value    : {float(arr.max()):.6e}\n"
        f"  max sym pos  : ({peak_x_fft}, {peak_y_fft}, {peak_z_fft})  (FFT-natural)\n"
        f"  rho[0,0,0]              = {centre_val:.6e}  (FFT-natural: cell centre)\n"
        f"  rho[nx/2,ny/2,nz/2]     = {corner_val:.6e}  (FFT-natural: -L/2 corner)\n"
    )


def main() -> int:
    if not RESULTS_DIR.exists():
        print(f"No results/ directory found at {RESULTS_DIR}")
        return 1

    series_roots = _find_series_roots(RESULTS_DIR)
    if not series_roots:
        print("No .raw/.meta.txt field series found under results/")
        return 1

    print("=== run_06 analysis: dual VTI conversion (buggy + fftshifted) ===")
    print(f"Run dir   : {RUN_DIR}")
    print(f"Results   : {RESULTS_DIR}")
    print(f"VTI buggy : {VTI_ROOT}")
    print(f"VTI fixed : {VTI_FFTSHIFT_ROOT}")
    print(f"Series    : {len(series_roots)}")

    converted = 0
    for series_dir in series_roots:
        rel = series_dir.relative_to(RESULTS_DIR)
        try:
            n_buggy, n_fixed = _convert_series_pair(series_dir)
            print(f"  OK  {rel}  |  buggy: {n_buggy}  |  fftshifted: {n_fixed}")
            converted += 1
        except Exception as exc:
            print(f"  FAIL {rel}  ->  {exc}")

    # Diagnostic summary for the total density only.
    total_density_dir = RESULTS_DIR / "density"
    summary_path = RESULTS_DIR / "grid_diagnostics_summary.txt"
    if total_density_dir.exists():
        summary = _summary_for_density(total_density_dir)
        if summary:
            summary_path.write_text(
                "# run_06 grid diagnostics summary (total density)\n"
                "#\n"
                "# If 'max sym pos' is near (0,0,0) and rho[0,0,0] >> rho[nx/2,ny/2,nz/2],\n"
                "# the INQ array is FFT-natural-ordered (cell centre at index 0). The\n"
                "# 'buggy' .vti will then render the molecule at the corners; the\n"
                "# 'fftshifted' .vti should render it at the metadata centre.\n\n"
                + summary
            )
            print(f"\nWrote summary: {summary_path.relative_to(RUN_DIR)}")

    print(f"Done. Converted {converted} / {len(series_roots)} series.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
