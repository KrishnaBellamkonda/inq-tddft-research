"""t=0 identity + conservation tests for GS-projected occupations (IV-M09).

n_i^GS(t) = Σ_j f_j(0)·|⟨ψ_i^GS|ψ_j(t)⟩|². At t=0 the evolved states ARE the GS
states, so the squared-overlap matrix is the identity and n_i^GS(0) = f_i — the
core sanity check for the projected-occupation phase. Pure numpy; tests the
compute core `_project_full` with analytically known matrices.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.pipeline.gs_projected_occupations import _project_full

pytestmark = pytest.mark.analysis


def test_t0_identity_returns_initial_occupations():
    """O_sq = I (t=0) ⇒ n_i^GS(0) = f_i exactly."""
    f = np.array([1.0, 1.0, 0.5, 0.0, 0.0])
    n0 = _project_full(np.eye(f.size), f)
    assert np.allclose(n0, f, atol=1e-12)


def test_complete_basis_conserves_total_occupation():
    """Column-stochastic O_sq (complete basis Σ_i O_ij = 1) conserves Σ f."""
    f = np.array([1.0, 1.0, 0.0, 0.0])
    O = np.array([[0.7, 0.3, 0.0, 0.0],
                  [0.3, 0.7, 0.0, 0.0],
                  [0.0, 0.0, 0.6, 0.4],
                  [0.0, 0.0, 0.4, 0.6]])
    n = _project_full(O, f)
    assert np.isclose(n.sum(), f.sum(), atol=1e-12)


def test_swap_moves_occupation_to_excited_state():
    """A swap overlap moves occupation from the occupied to the empty state."""
    f = np.array([1.0, 0.0])
    swap = np.array([[0.0, 1.0], [1.0, 0.0]])
    n = _project_full(swap, f)
    assert np.allclose(n, [0.0, 1.0], atol=1e-12)
