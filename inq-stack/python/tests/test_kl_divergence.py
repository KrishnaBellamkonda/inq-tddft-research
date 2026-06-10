"""Analytic tests for the KL-divergence kernel helpers
(``inqview.pipeline.kl_divergence._normalise`` / ``._kl``).

Pure-numpy, machine-independent. KL has exact closed-form values for known
distributions, so we assert against them directly. This underpins IV-M05
(``WPIntegrity.kl_mom`` reuses the momentum ``KL(P_t‖P_0)``).

The phase now lives at ``inqview.pipeline.kl_divergence`` (relocate step of
ADR-0003); when split, these pure helpers move to ``inqview.analysis`` — only
the import line changes.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.pipeline.kl_divergence import _kl, _normalise

pytestmark = pytest.mark.analysis


def test_normalise_sums_to_one():
    p = _normalise(np.array([1.0, 1.0, 2.0]))
    assert p.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(p, [0.25, 0.25, 0.5])


def test_normalise_zero_input_returns_zeros():
    """Degenerate all-zero histogram → zeros (no divide-by-zero)."""
    p = _normalise(np.zeros(4))
    assert np.all(p == 0.0)


def test_kl_self_is_zero():
    """KL(P‖P) = 0 for any distribution (Gibbs equality)."""
    p = _normalise(np.array([0.1, 0.4, 0.3, 0.2]))
    assert _kl(p, p) == pytest.approx(0.0, abs=1e-12)


def test_kl_known_value_nats():
    """KL([0.5,0.5] ‖ [0.25,0.75]) = 0.5 ln2 + 0.5 ln(2/3) ≈ 0.14384 nats."""
    p = np.array([0.5, 0.5])
    q = np.array([0.25, 0.75])
    expected = 0.5 * np.log(0.5 / 0.25) + 0.5 * np.log(0.5 / 0.75)
    assert _kl(p, q) == pytest.approx(expected, rel=1e-10)
    assert expected == pytest.approx(0.143841, abs=1e-5)


def test_kl_nonnegative():
    """KL(P‖Q) ≥ 0 for any P, Q (Gibbs' inequality)."""
    rng = np.random.default_rng(0)              # seeded → machine-independent
    for _ in range(20):
        p = _normalise(rng.random(8))
        q = _normalise(rng.random(8))
        assert _kl(p, q) >= -1e-12


def test_kl_is_asymmetric():
    """KL is not a metric: KL(P‖Q) ≠ KL(Q‖P) in general."""
    p = _normalise(np.array([0.7, 0.2, 0.1]))
    q = _normalise(np.array([0.2, 0.3, 0.5]))
    assert _kl(p, q) != pytest.approx(_kl(q, p), rel=1e-3)


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
