"""Wake-rendering helpers (the viz half of the wake phase, ADR 0003 split).

The bath-density math + VTI loading stay in ``inqview.pipeline.wake`` (VTK-
coupled). The shared colour-scale helper is a pure plotting concern and lives
here. Encodes the shared-colorbar rule (memory ``feedback_shared_colorbar_rule``):
directly-compared panels use ONE symmetric (vmin,vmax) about zero, with optional
percentile clipping to suppress lone spikes.
"""
from __future__ import annotations

import numpy as np


def shared_clim(*arrays, symmetric: bool = True, pct: float = 100.0):
    """One (vmin,vmax) over ALL arrays — for directly-compared panels.

    ``symmetric=True`` → (−m, m) about zero (signed Δn). ``pct<100`` clips to a
    percentile to suppress lone spikes (e.g. WP self-spike near the boundary).
    """
    m = 0.0
    for a in arrays:
        a = np.asarray(a)
        v = np.percentile(np.abs(a), pct) if pct < 100 else np.abs(a).max()
        m = max(m, float(v))
    return (-m, m) if symmetric else (0.0, m)
