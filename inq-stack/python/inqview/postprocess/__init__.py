"""DEPRECATED alias — this package was renamed to ``inqview.pipeline`` (ADR-0003).

Existing per-run ``analyse.py`` scripts import ``inqview.postprocess.*`` (e.g.
``from inqview.postprocess import pipeline`` / ``density_fourier``). This shim
forwards those imports to the new ``inqview.pipeline`` location so the run
scripts keep working unchanged. Prefer ``inqview.pipeline`` in new code.

Mechanism: ``__path__`` is pointed at the ``inqview/pipeline/`` source directory,
so ``import inqview.postprocess.<submodule>`` loads the renamed module from its
new home, and the package-level entry points are re-exported below.
"""
from __future__ import annotations

import warnings

from inqview import pipeline as _pipeline

# Submodule imports (inqview.postprocess.<name>) resolve to the pipeline/ dir.
__path__ = list(_pipeline.__path__)

# Package-level public entry points (inqview.postprocess.run, etc.).
from inqview.pipeline import PHASES, PipelineResult, run  # noqa: E402,F401

__all__ = ["PHASES", "PipelineResult", "run"]

warnings.warn(
    "inqview.postprocess was renamed to inqview.pipeline (ADR-0003); "
    "update imports to inqview.pipeline.",
    DeprecationWarning,
    stacklevel=2,
)
