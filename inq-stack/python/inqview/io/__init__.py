"""inqview.io — loaders and field/format dataclasses (ADR 0003).

The lowest layer: numpy only, NO matplotlib/VTK/scipy. Reads INQ/inqkit run
artefacts (field binaries + metadata, observables, LEED screens) into plain
dataclasses that the analysis and visualisation layers consume.
"""
from __future__ import annotations

from .fields import ComplexField3D, FieldMeta, RealField3D
from .data import (
    DataError,
    FieldSeries,
    SimulationData,
    infer_meta_path,
    load_complex_field,
    load_meta,
    load_real_field,
)
from .leed import LeedPattern, load_leed_pattern

__all__ = [
    "ComplexField3D",
    "FieldMeta",
    "RealField3D",
    "DataError",
    "FieldSeries",
    "SimulationData",
    "infer_meta_path",
    "load_complex_field",
    "load_meta",
    "load_real_field",
    "LeedPattern",
    "load_leed_pattern",
]
