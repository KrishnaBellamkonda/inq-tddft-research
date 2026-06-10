"""Unit tests for the shared colour-scale helper
(``inqview.pipeline.wake.shared_clim``).

Pure (numpy). Encodes the shared-colorbar rule: panels compared directly use
ONE symmetric (vmin,vmax) about zero; percentile clipping suppresses lone
spikes. (The VTI-reading bath math needs a run fixture and is covered by the
integration tier, not here.)

The phase now lives at ``inqview.pipeline.wake`` (relocate step of ADR-0003);
when it is split, ``shared_clim`` (a viz helper) follows to ``visualisation``
and the bath math to ``analysis`` — only the import line will move.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.pipeline.wake import shared_clim

pytestmark = pytest.mark.analysis


def test_symmetric_about_zero_uses_global_max_abs():
    """One symmetric scale over ALL arrays = (−m, m), m = global max |value|."""
    vmin, vmax = shared_clim(np.array([1.0, -3.0]), np.array([2.0, 0.5]))
    assert (vmin, vmax) == pytest.approx((-3.0, 3.0))


def test_asymmetric_mode_starts_at_zero():
    vmin, vmax = shared_clim(np.array([0.0, 2.0, 5.0]), symmetric=False)
    assert (vmin, vmax) == pytest.approx((0.0, 5.0))


def test_percentile_clips_lone_spike():
    """pct<100 clips to a percentile so a single outlier doesn't set the scale."""
    a = np.concatenate([np.ones(99), [1000.0]])      # one huge spike
    _, vmax = shared_clim(a, symmetric=False, pct=95.0)
    assert vmax < 1000.0                              # spike suppressed
    assert vmax == pytest.approx(np.percentile(np.abs(a), 95.0))


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
