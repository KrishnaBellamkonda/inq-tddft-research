from __future__ import annotations

"""
Post-fix analysis for run_06_centred_writer_check.

After the inqkit::fields::density fft-shift fix
(inq-stack/include/inqkit/fields/density.hpp), the .raw arrays already
contain the density in a contiguous left-to-right physical layout, so the
.vti conversion needs no further shifting. Convert every .raw + .meta.txt
field series under results/ into VTI files mirroring the source directory
structure under results/vti/.

For diagnosis, also writes results/grid_diagnostics_summary.txt with the
location of the density peak in the raw array. After the fix, the peak
must be near (nx/2, ny/2, nz/2) (= the metadata cell centre).
"""

import re
import sys
from pathlib import Path

import numpy as np

RUN_DIR = Path(__file__).resolve().parent
RESULTS_DIR = RUN_DIR / "results"
VTI_ROOT = RESULTS_DIR / "vti"

REPO_ROOT = RUN_DIR.parents[3] if len(RUN_DIR.parents) >= 4 else RUN_DIR
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))

try:
    import inqview
    from inqview.data import load_real_field
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
        roots.add(meta_path.parent)
    return sorted(roots)


def _collect_meta_series(series_dir: Path) -> list[Path]:
    metas = [p for p in series_dir.glob(f"*{META_SUFFIX}") if _paired_raw_exists(p)]
    return sorted(metas, key=_natural_key)


def _target_vti_dir(series_dir: Path) -> Path:
    rel = series_dir.relative_to(RESULTS_DIR)
    return VTI_ROOT / rel


def convert_series(series_dir: Path) -> tuple[Path, int]:
    meta_paths = _collect_meta_series(series_dir)
    if not meta_paths:
        raise RuntimeError(f"No .raw/.meta.txt pairs found in {series_dir}")
    out_dir = _target_vti_dir(series_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result = inqview.convert_real_series_to_vti(
        meta_paths, out_dir, array_name="density"
    )
    return out_dir, len(result.files)


def _summary_for_density(series_dir: Path) -> str | None:
    metas = _collect_meta_series(series_dir)
    if not metas:
        return None
    field = load_real_field(meta_path=metas[0])
    arr = np.asarray(field.array)
    nx, ny, nz = arr.shape
    ix, iy, iz = np.unravel_index(int(np.argmax(arr)), arr.shape)
    return (
        f"  series        : {series_dir.relative_to(RESULTS_DIR)}\n"
        f"  shape         : ({nx}, {ny}, {nz})\n"
        f"  metadata centre (nx/2, ny/2, nz/2) = ({nx//2}, {ny//2}, {nz//2})\n"
        f"  max array idx : ({ix}, {iy}, {iz})\n"
        f"  max value     : {float(arr.max()):.6e}\n"
        f"  rho[0,0,0]                 = {float(arr[0,0,0]):.6e}  "
        f"(should be ~0 after fix)\n"
        f"  rho[nx/2,ny/2,nz/2]        = {float(arr[nx//2, ny//2, nz//2]):.6e}  "
        f"(should be substantial after fix)\n"
    )


def main() -> int:
    if not RESULTS_DIR.exists():
        print(f"No results/ directory found at {RESULTS_DIR}")
        return 1

    series_roots = _find_series_roots(RESULTS_DIR)
    if not series_roots:
        print("No .raw/.meta.txt field series found under results/")
        return 1

    print("=== run_06 analysis (post-fix): VTI conversion ===")
    print(f"Run dir : {RUN_DIR}")
    print(f"Results : {RESULTS_DIR}")
    print(f"VTI out : {VTI_ROOT}")
    print(f"Series  : {len(series_roots)}")

    converted = 0
    for series_dir in series_roots:
        rel = series_dir.relative_to(RESULTS_DIR)
        try:
            out_dir, nfiles = convert_series(series_dir)
            print(f"  OK  {rel}  ->  {out_dir.relative_to(RUN_DIR)}  |  {nfiles} file(s)")
            converted += 1
        except Exception as exc:
            print(f"  FAIL {rel}  ->  {exc}")

    total_density_dir = RESULTS_DIR / "density"
    summary_path = RESULTS_DIR / "grid_diagnostics_summary.txt"
    if total_density_dir.exists():
        summary = _summary_for_density(total_density_dir)
        if summary:
            summary_path.write_text(
                "# run_06 grid diagnostics summary (post-fix, total density)\n"
                "#\n"
                "# After the writer fix, the density peak in the .raw array must\n"
                "# sit near (nx/2, ny/2, nz/2) (the metadata centre), and rho[0,0,0]\n"
                "# (cell -L/2 corner) must be ~0.\n\n" + summary
            )
            print(f"\nWrote summary: {summary_path.relative_to(RUN_DIR)}")

    print(f"Done. Converted {converted} / {len(series_roots)} series.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
