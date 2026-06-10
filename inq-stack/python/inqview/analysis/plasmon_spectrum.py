"""Plasmon spectral locator (IV-M01 / IV-M04).

Renamed from the old "loss function": the quantity |n_q(ω)|²/q² is a plasmon
**PEAK-LOCATOR** (peak positions ω_p(q) are reliable), NOT the quantitative
energy-loss function −Im[1/ε] (validated in
``docs/validation/loss-function-formula-validation.md``). It is the induced-
density power spectrum, q-weighted.

Two bugs from the old density_fourier are fixed here:
- IV-E01: the time-FFT uses the COMPLEX n_q(t) (``np.fft.fft``), not the real
  part only — so ±ω are not folded together.
- IV-E02: the returned quantity is |n_q(ω)|²/q² (power, q-weighted), not |n_q|.

Pure numpy. Both **axial** (q∥z, ``F[0,0,m]``) and **3d_binned** (all reciprocal
modes shell-averaged by |q|) modes are implemented. The VTI frame loader that
feeds either lives in ``inqview.pipeline`` (VTK is not deps-clean for analysis).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

HA_TO_EV = 27.211386245988


@dataclass(frozen=True)
class PlasmonSpectrum:
    """Induced-density power spectrum, q-weighted (peak-locator)."""

    q: np.ndarray            # |q| per mode (Bohr⁻¹)
    omega: np.ndarray        # angular frequency axis (Ha)
    power: np.ndarray        # |n_q(ω)|²              shape (n_omega, n_q)
    loss: np.ndarray         # |n_q(ω)|² / q²         shape (n_omega, n_q)
    peak_omega: np.ndarray   # |ω| at the spectral peak per mode (Ha)
    mode: str = "axial"

    @property
    def peak_omega_ev(self) -> np.ndarray:
        return self.peak_omega * HA_TO_EV


def extract_axial_nq(dn_series: np.ndarray, Lz: float, m_max: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """Axial Fourier modes n_q_m(t) = FFT3D(δn)[0,0,m] for m=1..m_max.

    ``dn_series`` is (n_t, nx, ny, nz) of δn(r,t)=n(r,t)−n(r,0) in the inqkit
    **x-slowest / z-fastest** layout, so axis −1 is z and ``F[0,0,m]`` is the
    longitudinal (q∥z) mode. Returns ``(n_q_t[n_t, m_max] complex,
    q_vals[m_max])`` with q_m = 2π m / Lz.
    """
    dn = np.asarray(dn_series)
    nt = dn.shape[0]
    n_q = np.empty((nt, m_max), dtype=complex)
    for i in range(nt):
        F = np.fft.fftn(dn[i])              # 3D spatial FFT
        for m in range(1, m_max + 1):
            n_q[i, m - 1] = F[0, 0, m]      # axial: q = (0,0,m)
    q_vals = np.array([2.0 * np.pi * m / Lz for m in range(1, m_max + 1)])
    return n_q, q_vals


def spectrum_from_nq(
    n_q_t: np.ndarray,
    dt: float,
    q_vals: np.ndarray,
    *,
    window: bool = True,
    t_start_index: int = 0,
    mode: str = "axial",
) -> PlasmonSpectrum:
    """Time-FFT each complex n_q(t) → |n_q(ω)|² and the q-weighted locator.

    ``n_q_t`` is (n_t, n_q) complex. A Hann window and an optional transient cut
    are applied. The peak ``|ω|`` per mode locates the plasmon at that q.
    """
    sig = np.asarray(n_q_t, dtype=complex)[t_start_index:]
    nt = sig.shape[0]
    if nt < 2:
        raise ValueError("need ≥2 time samples for the spectrum")
    win = np.hanning(nt)[:, None] if window else 1.0
    S = np.fft.fft(sig * win, axis=0)              # COMPLEX FFT (fixes IV-E01)
    omega = 2.0 * np.pi * np.fft.fftfreq(nt, d=dt)  # angular freq (Ha)

    power = np.abs(S) ** 2                          # |n_q(ω)|²  (fixes IV-E02)
    q = np.asarray(q_vals, dtype=float)
    loss = power / (q[None, :] ** 2)               # |n_q|² / q²
    peak = np.array([abs(omega[int(np.argmax(power[:, j]))]) for j in range(q.size)])
    return PlasmonSpectrum(q=q, omega=omega, power=power, loss=loss,
                           peak_omega=peak, mode=mode)


def spectrum_3d_binned(
    dn_series: np.ndarray,
    dt: float,
    spacing: float,
    *,
    n_bins: int = 32,
    q_max: float | None = None,
    window: bool = True,
    t_start_index: int = 0,
) -> PlasmonSpectrum:
    """Isotropic plasmon spectrum: shell-average |n_q(ω)|² over ALL reciprocal modes.

    ``dn_series`` is (n_t, nx, ny, nz) of δn(r,t). Each frame is 3D-FFT'd, the
    stack is time-FFT'd, the per-mode power |n_q(ω)|² is binned into ``n_bins``
    spherical |q| shells (cubic-grid spacing assumed isotropic). The DC mode
    (|q|=0) is excluded. Returns a :class:`PlasmonSpectrum` whose ``q`` are the
    shell centres and whose ``power``/``loss`` are shell-averaged.

    Memory note: holds the full (n_omega, nx, ny, nz) spectrum — use a coarse /
    strided density series for large runs (e.g. ``density_delta_coarse``).
    """
    dn = np.asarray(dn_series, dtype=float)
    if dn.ndim != 4:
        raise ValueError(f"dn_series must be (n_t,nx,ny,nz); got {dn.shape}")
    sig = dn[t_start_index:]
    nt, nx, ny, nz = sig.shape
    if nt < 2:
        raise ValueError("need ≥2 time samples for the spectrum")

    F = np.fft.fftn(sig, axes=(1, 2, 3))                      # spatial FFT/frame
    win = np.hanning(nt)[:, None, None, None] if window else 1.0
    S = np.fft.fft(F * win, axis=0)                           # time FFT (complex)
    omega = 2.0 * np.pi * np.fft.fftfreq(nt, d=dt)
    powmode = (np.abs(S) ** 2).reshape(nt, -1)               # (n_omega, n_modes)

    qx = 2.0 * np.pi * np.fft.fftfreq(nx, d=spacing)
    qy = 2.0 * np.pi * np.fft.fftfreq(ny, d=spacing)
    qz = 2.0 * np.pi * np.fft.fftfreq(nz, d=spacing)
    qmag = np.sqrt(qx[:, None, None] ** 2 + qy[None, :, None] ** 2
                   + qz[None, None, :] ** 2).ravel()

    pos = qmag > 0
    if q_max is None:
        q_max = float(qmag[pos].max())
    edges = np.linspace(0.0, q_max, n_bins + 1)
    q_centers = 0.5 * (edges[:-1] + edges[1:])
    shell = np.digitize(qmag, edges) - 1                     # shell per mode

    power_shell = np.zeros((nt, n_bins))
    for b in range(n_bins):
        mask = pos & (shell == b)
        if mask.any():
            power_shell[:, b] = powmode[:, mask].mean(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        loss = np.where(q_centers[None, :] > 0,
                        power_shell / q_centers[None, :] ** 2, 0.0)
    peak = np.array([abs(omega[int(np.argmax(power_shell[:, b]))])
                     for b in range(n_bins)])
    return PlasmonSpectrum(q=q_centers, omega=omega, power=power_shell,
                           loss=loss, peak_omega=peak, mode="3d_binned")
