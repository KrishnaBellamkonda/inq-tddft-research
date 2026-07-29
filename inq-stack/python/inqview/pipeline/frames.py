"""Density frame-series loader (data-coupled glue for the analysis kernels).

Reads a run's ``results/raw/vti/<category>/density_t*.vti`` series into a numpy
stack so the deps-clean analysis kernels (`plasmon_spectrum.spectrum_3d_binned`,
`center_of_density.compare`, `wp_integrity` ipr) can run on real runs. Lives in
``pipeline`` (not ``analysis``/``io``) because VTI reading needs VTK — kept out
of the deps-clean layers. VTK is imported lazily via the existing
``density._load_vti_array``, so importing this module pulls no VTK until called.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np

_STEP_RE = re.compile(r"_t(\d+)\.vti$")


def _step_of(path: Path) -> int:
    m = _STEP_RE.search(path.name)
    return int(m.group(1)) if m else -1


def load_density_series(
    results_dir,
    category: str = "density_total",
    *,
    max_frames: Optional[int] = None,
):
    """Load a VTI density series → ``(cubes, steps, meta)``.

    ``cubes`` is (n_t, nx, ny, nz); ``steps`` is the per-frame integer step
    parsed from the filenames (``_t<step>``), ordered ascending (time_au =
    step·dt — the caller multiplies by the run dt); ``meta`` is the last frame's
    VTI metadata. Raises ``FileNotFoundError`` if the category has no frames.
    """
    from . import _common
    from .density import _load_vti_array

    cat_dir = Path(results_dir) / "raw" / "vti" / category
    files = _common.list_vti_series(cat_dir, category)
    if not files:
        raise FileNotFoundError(f"no VTI frames under {cat_dir}")
    if max_frames is not None:
        files = files[:max_frames]

    cubes, steps, meta = [], [], None
    for f in files:
        cube, meta = _load_vti_array(f)
        cubes.append(np.asarray(cube, dtype=float))
        steps.append(_step_of(Path(f)))
    return np.stack(cubes), np.asarray(steps, dtype=int), meta


def load_dn_series(results_dir, category: str = "density_total", **kw):
    """Load δn(r,t) = n(r,t) − n(r,0) for a category → ``(dn, steps, meta)``."""
    cubes, steps, meta = load_density_series(results_dir, category, **kw)
    return cubes - cubes[0][None, ...], steps, meta
