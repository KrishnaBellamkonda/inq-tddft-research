"""Proper Orthogonal Decomposition (POD / PCA / Karhunen-Loeve) via truncated SVD.

Campaign-local kernel (ml-patterns). PRE-GATED in T1 (code-test + formula-validation).

Formula (Lumley 1967; review arXiv:2111.04829 "Modal analysis"; Brunton & Kutz,
"Data-Driven Science and Engineering", Ch. 1):
  Snapshot matrix  X in R^{n x m}, n = #spatial points (features), m = #snapshots.
  Economy SVD       X = U S V^T,  U in R^{n x r}, S diag, V in R^{m x r}.
  POD spatial modes = columns of U (orthonormal spatial fields).
  Mode energy       E_i = S_i^2,  energy fraction f_i = S_i^2 / sum_j S_j^2.
  Temporal coeffs   a = S V^T  (row i = time evolution of mode i).
  Rank-k reconstruction  X_k = U[:,:k] S[:k] V[:,:k]^T minimises ||X - X_k||_F
  over all rank-k matrices (Eckart-Young).

We optionally mean-centre columns (PCA convention) before the SVD; for the campaign
the GS subtraction already removes the dominant offset, so centring is optional.

Implementation uses scipy/numpy SVD with optional randomized truncation
(Halko, Martinsson, Tropp, SIAM Rev. 53, 217 (2011)) for memory-bounded large n.
float32 throughout (campaign memory budget).
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass


@dataclass(frozen=True)
class PODResult:
    modes: np.ndarray          # (n, r) spatial modes (columns), orthonormal
    singular_values: np.ndarray  # (r,)
    coeffs: np.ndarray         # (r, m) temporal coefficients a = S V^T
    energy_fraction: np.ndarray  # (r,) S_i^2 / sum S^2
    mean: np.ndarray | None    # (n,) column mean if centred else None

    @property
    def cumulative_energy(self) -> np.ndarray:
        return np.cumsum(self.energy_fraction)

    def n_modes_for(self, frac: float) -> int:
        """Smallest #modes whose cumulative energy >= frac."""
        return int(np.searchsorted(self.cumulative_energy, frac) + 1)


def pod(X: np.ndarray, rank: int | None = None, center: bool = False,
        randomized: bool = False, n_oversamples: int = 10,
        random_state: int = 0) -> PODResult:
    """POD of a snapshot matrix.

    Parameters
    ----------
    X : (n_features, n_snapshots) array. Columns are snapshots (flattened fields).
    rank : truncation rank r (None -> min(n, m)).
    center : subtract column mean (PCA convention) before SVD.
    randomized : use randomized SVD truncation (memory-bounded for large n).
    """
    X = np.asarray(X, dtype=np.float32)
    if X.ndim != 2:
        raise ValueError("X must be 2D (n_features, n_snapshots)")
    n, m = X.shape
    full_rank = min(n, m)
    r = full_rank if rank is None else min(rank, full_rank)

    mean = None
    Xc = X
    if center:
        mean = X.mean(axis=1, keepdims=True)
        Xc = X - mean

    if randomized and r < full_rank:
        U, S, Vt = _randomized_svd(Xc, r, n_oversamples, random_state)
    else:
        # economy SVD then truncate
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        U, S, Vt = U[:, :r], S[:r], Vt[:r, :]

    total = float(np.sum(S.astype(np.float64) ** 2))
    energy = (S.astype(np.float64) ** 2) / total if total > 0 else np.zeros_like(S, dtype=np.float64)
    coeffs = (S[:, None] * Vt).astype(np.float32)
    return PODResult(
        modes=U.astype(np.float32),
        singular_values=S.astype(np.float32),
        coeffs=coeffs,
        energy_fraction=energy.astype(np.float64),
        mean=(mean.ravel() if mean is not None else None),
    )


def _randomized_svd(X, r, p, seed, n_iter=4):
    """Randomized SVD with subspace (power) iterations (Halko et al. 2011, Alg 4.4).

    Power iterations sharpen the captured subspace when the spectrum decays slowly.
    """
    rng = np.random.default_rng(seed)
    n, m = X.shape
    ell = min(r + p, m)
    Omega = rng.standard_normal((m, ell)).astype(np.float32)
    Y = X @ Omega
    Q, _ = np.linalg.qr(Y)
    for _ in range(n_iter):
        Z, _ = np.linalg.qr(X.T @ Q)
        Q, _ = np.linalg.qr(X @ Z)
    B = Q.T @ X
    Ub, S, Vt = np.linalg.svd(B, full_matrices=False)
    U = Q @ Ub
    return U[:, :r], S[:r], Vt[:r, :]
