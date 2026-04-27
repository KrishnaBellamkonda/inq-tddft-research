"""inqview.postprocess — generalisable post-processing pipeline.

Consumes a ``results/`` directory laid out per
``docs/results_folder_structure_spec.md`` and produces every figure / GIF /
summary artefact the spec requires. The pipeline is split into phases so a
caller (e.g. the coronene-specific ``coronene_postprocess.py`` CLI) can
choose to run all of them or a subset.

Phases (in dependency order):
    summary      — verify and (re)build ``results/run_summary.txt``
    gs           — ground-state plots (orbital gallery, density)
    observables  — time-domain + FFT plots from ``raw/observables/``
    density      — 2D slice GIFs for total / system / wp densities
    screens      — total / instantaneous / time-windowed LEED + coord checks
    overlap      — WP-overlap-with-GS-orbitals bar-chart GIF
    orbitals     — RT orbital gallery (only if RT orbital VTIs were emitted)
    paraview     — 3D volume renders via pvbatch (slowest; opt-out)

Public entry point: :func:`run`.
"""

from __future__ import annotations

from .pipeline import PHASES, PipelineResult, run

__all__ = ["PHASES", "PipelineResult", "run"]
