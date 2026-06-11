"""inqview — post-processing & visualisation for INQ/inqkit TDDFT runs.

Public API is layered (ADR 0003):

- ``inqview.io``            loaders + field/format dataclasses (numpy only)
- ``inqview.analysis``      numeric kernels → frozen dataclasses (numpy/scipy/pandas)
- ``inqview.visualisation`` all rendering (matplotlib, VTK/paraview, GIF) + theme
- ``inqview.pipeline``      thin phase orchestration

The names below are kept importable from the top level (``from inqview import
RealField3D``) for backward compatibility, but they are resolved **lazily** via
PEP 562 ``__getattr__``: importing ``inqview`` — or any deps-clean subpackage
such as ``inqview.analysis`` — pulls in NO matplotlib and NO VTK until a
plotting/IO name is actually accessed. This is what makes the ADR-0003
deps-clean invariant hold (see ``tests/test_deps_clean.py``).
"""
from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

# Public name -> dotted submodule (relative to this package) that defines it.
# Grouped by layer so the eventual io/visualisation re-homing is a one-line edit.
_LAZY_EXPORTS = {
    # NOTE: legacy `config` (Theme/DEFAULT_THEME/PlotDefaults/RenderDefaults) and
    # `defaults` (one-call movie wrappers) are NO LONGER part of the public API
    # (ADR 0004 — the canonical theme is `inqview.visualisation.style`). They
    # remain as internal modules (`plots` uses `config` internally) but are not
    # re-exported from the top level.
    # io: fields
    "FieldMeta": "io.fields",
    "RealField3D": "io.fields",
    "ComplexField3D": "io.fields",
    # io: data loaders
    "DataError": "io.data",
    "FieldSeries": "io.data",
    "SimulationData": "io.data",
    "infer_meta_path": "io.data",
    "load_meta": "io.data",
    "load_real_field": "io.data",
    "load_complex_field": "io.data",
    # io: LEED screen loader (renamed from the colliding `screens` phase name)
    "LeedPattern": "io.leed",
    "load_leed_pattern": "io.leed",
    # analysis: fourier kernel (moved out of the flat layout — ADR 0003)
    "FourierResult": "analysis.fourier",
    "FourierTransform": "analysis.fourier",
    "WindowSpec": "analysis.fourier",
    # visualisation: VTI writer (legacy Python writer; verify-then-cut)
    "write_vti": "visualisation.vti",
    "convert_real_meta_to_vti": "visualisation.vti",
    "convert_real_series_to_vti": "visualisation.vti",
    "VTISeriesResult": "visualisation.vti",
    # visualisation: paraview pipeline
    "VolumeRenderSpec": "visualisation.paraview",
    "AnimationSpec": "visualisation.paraview",
    "ParaViewPipeline": "visualisation.paraview",
    # visualisation: matplotlib plots
    "load_observables": "visualisation.plots",
    "plot_energy_vs_time": "visualisation.plots",
    "plot_total_energy_vs_time": "visualisation.plots",
    "plot_all_energy_components_vs_time": "visualisation.plots",
    "plot_current_vs_time": "visualisation.plots",
    "plot_dipole_vs_time": "visualisation.plots",
    "plot_observables_summary": "visualisation.plots",
    "plot_spectrum": "visualisation.plots",
    "plot_spectrum_summary": "visualisation.plots",
    "plot_density_slice": "visualisation.plots",
    "plot_leed_pattern": "visualisation.plots",
    # pipeline: default render-spec wrappers (used by run analysis.py scripts —
    # kept public; only the legacy config/theme names were dropped)
    "default_density_movie": "defaults",
    "default_wavepacket_movie": "defaults",
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str):  # PEP 562
    target = _LAZY_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = importlib.import_module(f"{__name__}.{target}")
    value = getattr(module, name)
    globals()[name] = value  # cache so subsequent access skips __getattr__
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY_EXPORTS))


if TYPE_CHECKING:  # let type-checkers/IDEs still see the names eagerly
    from .io.fields import ComplexField3D, FieldMeta, RealField3D
    from .io.data import (
        DataError,
        FieldSeries,
        SimulationData,
        infer_meta_path,
        load_complex_field,
        load_meta,
        load_real_field,
    )
    from .io.leed import LeedPattern, load_leed_pattern
    from .analysis.fourier import FourierResult, FourierTransform, WindowSpec
    from .visualisation.vti import (
        VTISeriesResult,
        convert_real_meta_to_vti,
        convert_real_series_to_vti,
        write_vti,
    )
    from .visualisation.paraview import AnimationSpec, ParaViewPipeline, VolumeRenderSpec
    from .visualisation.plots import (
        load_observables,
        plot_all_energy_components_vs_time,
        plot_current_vs_time,
        plot_density_slice,
        plot_dipole_vs_time,
        plot_energy_vs_time,
        plot_leed_pattern,
        plot_observables_summary,
        plot_spectrum,
        plot_spectrum_summary,
        plot_total_energy_vs_time,
    )
    from .defaults import default_density_movie, default_wavepacket_movie
