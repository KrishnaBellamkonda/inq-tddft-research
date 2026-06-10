"""inqview.analysis — numeric post-processing kernels (ADR 0003).

Imports only numpy / scipy / pandas; returns plain frozen dataclasses. NEVER
imports matplotlib or VTK (the deps-clean invariant — a headless node computes
observables without plotting deps). Renderers that consume these results live in
``inqview.visualisation``.
"""
from __future__ import annotations

# NB: do NOT re-export the bare `center_of_density` function name here — it would
# shadow the submodule of the same name. Use compare_cod, or import the function
# from inqview.analysis.center_of_density directly.
from .center_of_density import COD, CODComparison, compare as compare_cod
from .energy_components import EnergyComponents, compute as compute_energy_components
from .wp_integrity import (
    WPIntegrity,
    assemble_from_run as assemble_wp_integrity,
    ipr,
    kl_series,
    momentum_kl,
)
from .plasmon_spectrum import (
    PlasmonSpectrum,
    extract_axial_nq,
    spectrum_3d_binned,
    spectrum_from_nq,
)
from .fourier import FourierResult, FourierTransform, WindowSpec

__all__ = [
    "COD",
    "CODComparison",
    "compare_cod",
    "EnergyComponents",
    "compute_energy_components",
    "WPIntegrity",
    "assemble_wp_integrity",
    "ipr",
    "kl_series",
    "momentum_kl",
    "PlasmonSpectrum",
    "extract_axial_nq",
    "spectrum_3d_binned",
    "spectrum_from_nq",
    "FourierResult",
    "FourierTransform",
    "WindowSpec",
]
