"""Dynamic Mode Decomposition (exact DMD) — windowed for non-stationary data.

Campaign-local kernel (ml-patterns). PRE-GATED in T1 (code-test + formula-validation).

Formula (Schmid, J. Fluid Mech. 656, 5 (2010); Tu, Rowley, Luchtenburg, Brunton,
Kutz, "On Dynamic Mode Decomposition", J. Comput. Dyn. 1, 391 (2014) — exact DMD;
Brunton & Kutz, "Data-Driven Science and Engineering" Ch. 7):

  Snapshots x_0 .. x_m sampled at uniform interval dt.
  X  = [x_0 ... x_{m-1}]  (n x m)
  X' = [x_1 ... x_m]      (n x m)         (shifted by one dt)
  Seek best-fit linear operator A with X' ~ A X.
  Reduced SVD  X = U S V*  truncated to rank r.
  Atilde = U* X' V S^{-1}                  (r x r)
  Eigendecomposition  Atilde W = W Lambda.
  DMD (exact) modes  Phi = X' V S^{-1} W   (columns; each tied to one eigenvalue).
  Discrete eigenvalues lambda_i -> continuous  omega_i = ln(lambda_i) / dt.
    angular frequency  Re/Im:  oscillation frequency  f_i = Im(omega_i)/(2 pi),
    angular freq       w_i     = Im(omega_i),
    growth/decay rate  g_i     = Re(omega_i)  (g<0 => decaying).
  Mode amplitudes b solve  Phi b = x_0  (least squares).

Windowing: for a decelerating light projectile the dynamics are non-stationary
(project rule light-projectile-stopping). DMD is applied over a user-supplied
window of snapshots (the early near-constant-velocity stretch). Nyquist guard:
caller must ensure dt < pi/omega_p so the plasmon mode is resolved.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass


@dataclass(frozen=True)
class DMDResult:
    eigenvalues: np.ndarray    # (r,) discrete lambda_i (complex)
    omega: np.ndarray          # (r,) continuous ln(lambda)/dt (complex)
    modes: np.ndarray          # (n, r) DMD spatial modes (complex)
    amplitudes: np.ndarray     # (r,) b from Phi b = x_0
    dt: float

    @property
    def frequency(self) -> np.ndarray:
        """Oscillation frequency f = Im(omega)/(2 pi)  [1/time]."""
        return np.imag(self.omega) / (2.0 * np.pi)

    @property
    def angular_frequency(self) -> np.ndarray:
        """Angular frequency Im(omega)  [rad/time]."""
        return np.imag(self.omega)

    @property
    def growth_rate(self) -> np.ndarray:
        """Re(omega): <0 decaying, >0 growing."""
        return np.real(self.omega)

    def dominant(self, positive_freq_only: bool = True):
        """Index of the largest-amplitude mode with |Im(omega)|>0.

        Returns (index, angular_frequency, growth_rate, amplitude).
        Eigenvalues come in conjugate pairs; restricting to positive angular
        frequency picks one representative per physical oscillation.
        """
        amp = np.abs(self.amplitudes)
        w = self.angular_frequency
        mask = np.abs(w) > 1e-9
        if positive_freq_only:
            mask &= (w > 0)
        if not np.any(mask):
            mask = np.abs(w) > 1e-9
        idx_candidates = np.where(mask)[0]
        if idx_candidates.size == 0:
            i = int(np.argmax(amp))
        else:
            i = int(idx_candidates[np.argmax(amp[idx_candidates])])
        return i, float(abs(w[i])), float(self.growth_rate[i]), float(amp[i])


def dmd(snapshots: np.ndarray, dt: float, rank: int | None = None,
        window: tuple[int, int] | None = None) -> DMDResult:
    """Exact DMD of a uniformly-sampled snapshot sequence.

    Parameters
    ----------
    snapshots : (n_features, n_times) real or complex array (columns = times).
    dt : sampling interval.
    rank : SVD truncation rank (None -> full).
    window : (start, stop) slice of time columns to use (early constant-v window).
    """
    Xall = np.asarray(snapshots)
    if Xall.ndim != 2:
        raise ValueError("snapshots must be 2D (n_features, n_times)")
    if window is not None:
        Xall = Xall[:, window[0]:window[1]]
    Xall = Xall.astype(np.complex128 if np.iscomplexobj(Xall) else np.float64)
    if Xall.shape[1] < 3:
        raise ValueError("need >= 3 snapshots for DMD")

    X = Xall[:, :-1]
    Xp = Xall[:, 1:]
    n, m = X.shape
    full = min(n, m)
    r = full if rank is None else min(rank, full)

    U, S, Vh = np.linalg.svd(X, full_matrices=False)
    U, S, Vh = U[:, :r], S[:r], Vh[:r, :]
    V = Vh.conj().T
    Sinv = np.diag(1.0 / S)

    Atilde = U.conj().T @ Xp @ V @ Sinv
    lam, W = np.linalg.eig(Atilde)
    # exact DMD modes
    Phi = Xp @ V @ Sinv @ W

    omega = np.log(lam.astype(np.complex128)) / dt
    # amplitudes b: Phi b = x_0
    x0 = Xall[:, 0]
    b, *_ = np.linalg.lstsq(Phi, x0, rcond=None)

    return DMDResult(eigenvalues=lam, omega=omega, modes=Phi, amplitudes=b, dt=float(dt))
