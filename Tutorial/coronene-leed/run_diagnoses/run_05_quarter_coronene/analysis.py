from __future__ import annotations

"""
Convert every .raw + .meta.txt field series under results/ into VTI files,
mirroring the source directory structure under results/vti/.

Run from the simulation directory (the directory containing results/).
"""

import re
import sys
from pathlib import Path
from typing import Iterable

RUN_DIR = Path(__file__).resolve().parent
RESULTS_DIR = RUN_DIR / "results"
VTI_ROOT = RESULTS_DIR / "vti"

# run_diagnoses/run_XX/ -> parents[3] is the repo root
REPO_ROOT = RUN_DIR.parents[3] if len(RUN_DIR.parents) >= 4 else RUN_DIR
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))

try:
    import inqview
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


def _parse_meta(meta_path: Path) -> dict[str, str]:
    meta: dict[str, str] = {}
    for raw_line in meta_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        meta[key.strip()] = value.strip()
    return meta


def _paired_raw_exists(meta_path: Path) -> bool:
    name = meta_path.name
    if not name.endswith(META_SUFFIX):
        return False
    raw_name = name[: -len(META_SUFFIX)] + ".raw"
    return (meta_path.parent / raw_name).exists()


def _infer_array_name(meta_paths: list[Path]) -> str:
    if not meta_paths:
        return "density"
    meta = _parse_meta(meta_paths[0])
    field_name = meta.get("field_name", "").strip()
    if field_name:
        return field_name
    parent_name = meta_paths[0].parent.name
    if parent_name.startswith("orbital_"):
        return "density"
    return parent_name or "density"


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
    metas = [
        p for p in series_dir.glob(f"*{META_SUFFIX}")
        if _paired_raw_exists(p)
    ]
    return sorted(metas, key=_natural_key)


def _target_vti_dir(series_dir: Path) -> Path:
    rel = series_dir.relative_to(RESULTS_DIR)
    return VTI_ROOT / rel


def convert_series(series_dir: Path) -> tuple[Path, int, str]:
    meta_paths = _collect_meta_series(series_dir)
    if not meta_paths:
        raise RuntimeError(f"No .raw/.meta.txt pairs found in {series_dir}")

    out_dir = _target_vti_dir(series_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    array_name = _infer_array_name(meta_paths)
    result = inqview.convert_real_series_to_vti(
        meta_paths,
        out_dir,
        array_name=array_name,
    )
    return out_dir, len(result.files), array_name


def main() -> int:
    if not RESULTS_DIR.exists():
        print(f"No results/ directory found at {RESULTS_DIR}")
        return 1

    series_roots = _find_series_roots(RESULTS_DIR)
    if not series_roots:
        print("No .raw/.meta.txt field series found under results/")
        return 1

    print("=== Convert all .raw series to VTI ===")
    print(f"Run dir : {RUN_DIR}")
    print(f"Results : {RESULTS_DIR}")
    print(f"VTI out : {VTI_ROOT}")
    print(f"Series found: {len(series_roots)}")

    converted = 0
    for series_dir in series_roots:
        rel = series_dir.relative_to(RESULTS_DIR)
        try:
            out_dir, nfiles, array_name = convert_series(series_dir)
            print(f"  OK  {rel}  ->  {out_dir.relative_to(RUN_DIR)}  |  {nfiles} file(s)  |  array='{array_name}'")
            converted += 1
        except Exception as exc:
            print(f"  FAIL {rel}  ->  {exc}")

    print(f"Done. Converted {converted} / {len(series_roots)} series.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
