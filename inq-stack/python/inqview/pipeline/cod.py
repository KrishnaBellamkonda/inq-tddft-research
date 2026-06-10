"""From-run centre-of-density glue (IV-M02).

Loads a run's density VTI series and computes the COD time-series with the
correct NODE convention (`analysis.center_of_density`), plus a cross-check
helper against the inqkit-written `observables.csv` cod_* columns (which use the
E04 half-cell convention). VTK loading lives here in `pipeline`, not in the
deps-clean `analysis` kernel.
"""
from __future__ import annotations

import numpy as np


def cod_series_from_run(
    results_dir,
    category: str = "density_total",
    *,
    origin=None,
    spacing=None,
    half_cell: bool = False,
):
    """Load the density series and return ``(steps, list[COD])`` (node convention).

    ``origin``/``spacing`` default to the VTI metadata. Set ``half_cell=True`` to
    reproduce the inqkit (i+½)·dx convention.
    """
    from ..analysis.center_of_density import center_of_density
    from .frames import load_density_series

    cubes, steps, meta = load_density_series(results_dir, category)
    o = origin if origin is not None else _meta_get(meta, "origin")
    s = spacing if spacing is not None else _meta_get(meta, "spacing")
    cods = [center_of_density(c, o, s, half_cell=half_cell) for c in cubes]
    return steps, cods


def cod_offset_vs_inqkit(node_xyz, inqkit_xyz) -> np.ndarray:
    """Mean per-axis ``inqkit_csv_COD − python_node_COD`` over the series.

    Both arrays are (n, 3) in Bohr. The result documents the **E04** bug: it
    should equal ``(dx, dy, dz)/2`` because the inqkit CSV uses the half-cell
    convention while the python COD uses the node convention.
    """
    return (np.asarray(inqkit_xyz, dtype=float)
            - np.asarray(node_xyz, dtype=float)).mean(axis=0)


def _meta_get(meta, key):
    if meta is None or key not in meta:
        raise KeyError(
            f"VTI metadata has no {key!r}; pass origin/spacing explicitly "
            "(e.g. from run_summary)."
        )
    return meta[key]
