"""Pure KL-divergence helpers (analysis layer, ADR 0003 split).

The numeric core of the ``kl_divergence`` phase: normalise a histogram and the
forward KL ``Σ_k p_k log(p_k/q_k)`` over bins with ``q_k>0``. Pure numpy, no
matplotlib — so headless analysis can compute KL series. The phase's plotting +
run() stay in ``inqview.pipeline.kl_divergence``. (Mirrors the KL reimplemented
inline in ``wp_integrity`` to keep that module deps-clean.)
"""
from __future__ import annotations

import numpy as np

EPS = 1e-300


def _normalise(p: np.ndarray) -> np.ndarray:
    """Return p / sum(p) (safe; returns zeros if sum is zero)."""
    s = float(np.asarray(p, dtype=float).sum())
    if s <= 0:
        return np.zeros_like(np.asarray(p, dtype=float))
    return np.asarray(p, dtype=float) / s


def _kl(p: np.ndarray, q: np.ndarray) -> float:
    """Σ_k p_k log(p_k / q_k), over bins with q_k > 0 only (p, q pre-normalised)."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    mask = q > 0
    if not mask.any():
        return float("nan")
    pp = np.where(mask, p, 0.0) + EPS
    qq = np.where(mask, q, 1.0)
    contrib = np.where(mask, pp * np.log(pp / qq), 0.0)
    return float(contrib.sum())
