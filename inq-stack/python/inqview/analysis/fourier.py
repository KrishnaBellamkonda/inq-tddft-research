"""
fourier.py — FFT-based post-processing for inqview observables.

Provides FourierTransform (class-based) with tunable windowing, detrending,
and convenience methods for energy, current, and dipole columns.
"""

# TODO: Need to investigate carefully if the FT windowing, detrending, and convenience
# methods are skewing the findings somehow. Especially for the QuantumKickExtension 
# run. 

# TODO: Also, for what kind of tasks is fourier.py being used now?


from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.signal import detrend as scipy_detrend
from scipy.signal import get_window

_VALID_WINDOWS = ("boxcar", "hann", "hamming", "blackman", "tukey", "kaiser", "flattop")


@dataclass
class WindowSpec:
    """Window function selection with tunable parameters.

    name   : one of boxcar | hann | hamming | blackman | tukey | kaiser | flattop
    alpha  : Tukey taper fraction in [0, 1]  (used only when name='tukey')
    beta   : Kaiser shape parameter >= 0     (used only when name='kaiser')
    """
    name: str = "hann"
    alpha: float = 0.5
    beta: float = 14.0

    def __post_init__(self) -> None:
        if self.name not in _VALID_WINDOWS:
            raise ValueError(
                f"Unknown window '{self.name}'. Valid choices: {_VALID_WINDOWS}"
            )

    def build(self, n: int) -> np.ndarray:
        """Return a length-n window array."""
        if self.name == "tukey":
            return get_window(("tukey", self.alpha), n)
        if self.name == "kaiser":
            return get_window(("kaiser", self.beta), n)
        return get_window(self.name, n)


@dataclass
class FourierResult:
    """Result of a single Fourier transform.

    frequency_au : positive-frequency axis in atomic units (1/time_au = Ha/hbar)
    amplitude    : |FFT| one-sided, normalised by N
    power        : amplitude**2 (power spectrum)
    column       : source observable column name
    dt_au        : timestep used
    window       : WindowSpec that was applied
    """
    frequency_au: np.ndarray
    amplitude: np.ndarray
    power: np.ndarray
    column: str
    dt_au: float
    window: WindowSpec
    t_start_au: float = 0.0       # §13.6 transient cutoff applied (0 = none)
    subtract: str = "mean"        # baseline removed before windowing (IV-M12)


class FourierTransform:
    """Windowed FFT of TDDFT time-series observables.

    Parameters
    ----------
    window : WindowSpec
        Window function and its parameters. Default: Hann window.
    detrend : bool | None
        Legacy switch, retained for back-compat. None (default) → the canonical
        'mean' baseline is used (user verdict 2026-06-25). True → 'detrend'
        (linear trend removed); False → 'none'. Prefer the explicit ``subtract=``.
    zero_pad : int
        Multiplier on signal length applied via zero-padding before the FFT.
        Larger values give a smoother frequency-axis interpolation (no extra
        spectral information, just denser sampling). Default 4 for spectra
        comparable to QBall's smooth-noise output. Set 1 to disable.
    smooth_sigma_bins : float
        Optional Gaussian smoothing kernel width (in frequency bins) applied
        to the magnitude spectrum AFTER the FFT. Default 0 (no smoothing).
        Mild values (0.5–1.0) clean up high-frequency Gibbs ripple at the
        cost of slight resolution loss; pair with a Hann or Kaiser window.
    t_start_au : float
        Transient-region cutoff: discard all samples with ``time_au -
        time_au[0] < t_start_au`` before windowing + FFT (per
        observables_reference §13.6). Default 0 (no skip). Use this to
        exclude the WP-injection shake-up so the spectrum reflects the
        steady-state response. The cutoff is recorded on the returned
        ``FourierResult.t_start_au`` field for downstream captions / CSV
        headers.
    """

    _SUBTRACT = ("initial", "mean", "detrend", "none")

    def __init__(
        self,
        window: WindowSpec | None = None,
        detrend: bool | None = None,
        zero_pad: int = 4,
        smooth_sigma_bins: float = 0.0,
        t_start_au: float = 0.0,
        subtract: str | None = None,
    ) -> None:
        self.window = window if window is not None else WindowSpec("hann")
        # Baseline removal (IV-M12). `subtract` supersedes the legacy `detrend`
        # bool. CANONICAL DEFAULT = 'mean' (user verdict 2026-06-25, overriding
        # the fft-drift-removal dossier's per-observable split): a single uniform
        # pipeline — mean removal → Hann → 4x pad → coherent-gain rfft — for EVERY
        # observable (dipole, current, energy). `detrend` is retained as an opt-in
        # comparison (the panel overlays it dashed). Back-compat: an explicit
        # `subtract=` wins; an explicit `detrend=True/False` maps to
        # 'detrend'/'none'; only when BOTH are unset do we fall to 'mean'.
        self.subtract = (
            subtract if subtract is not None
            else "detrend" if detrend is True
            else "none" if detrend is False
            else "mean"
        )
        if self.subtract not in self._SUBTRACT:
            raise ValueError(
                f"unknown subtract={self.subtract!r}; valid: {self._SUBTRACT}")
        self.detrend = self.subtract == "detrend"
        self.zero_pad = max(1, int(zero_pad))
        self.smooth_sigma_bins = float(smooth_sigma_bins)
        self.t_start_au = float(t_start_au)

    def _apply_subtract(self, values: np.ndarray) -> np.ndarray:
        """Remove the chosen baseline before windowing (IV-M12).

        'initial' enforces the induced-response IC Δs(0)=0 (Yabana-Bertsch);
        'mean' removes the DC level; 'detrend' removes a linear trend (best
        when a genuine drift exists, e.g. energy); 'none' leaves the signal —
        with which a constant offset hijacks the ω≈0 bin (todo.txt #2).
        """
        if self.subtract == "none":
            return values.copy()
        if self.subtract == "initial":
            return values - values[0]
        if self.subtract == "mean":
            return values - float(np.mean(values))
        return scipy_detrend(values, type="linear")        # 'detrend'

    # ------------------------------------------------------------------
    # Core
    # ------------------------------------------------------------------

    def transform(
        self,
        time_au: np.ndarray,
        values: np.ndarray,
        column: str = "",
    ) -> FourierResult:
        """Transform a time-series array into the frequency domain.

        Parameters
        ----------
        time_au : 1-D array of sample times in atomic units (must be uniform).
        values  : 1-D array of real-valued observable samples.
        column  : label stored in the returned FourierResult.
        """
        time_au = np.asarray(time_au, dtype=float)
        values = np.asarray(values, dtype=float)

        if time_au.ndim != 1 or values.ndim != 1:
            raise ValueError("time_au and values must be 1-D arrays.")
        if len(time_au) != len(values):
            raise ValueError("time_au and values must have the same length.")
        if len(time_au) < 2:
            raise ValueError("Need at least 2 samples for FFT.")

        # §13.6 transient skip: drop samples with t < t_start_au (relative
        # to the first time point) before windowing / FFT.
        if self.t_start_au > 0.0:
            t_cut = float(time_au[0]) + self.t_start_au
            mask = time_au >= t_cut
            if mask.sum() < 4:
                raise ValueError(
                    f"FourierTransform: t_start_au={self.t_start_au} a.u. "
                    f"leaves <4 samples ({mask.sum()})")
            time_au = time_au[mask]
            values = values[mask]

        n = len(time_au)
        dt = float(time_au[1] - time_au[0])

        sig = self._apply_subtract(values)
        win = self.window.build(n)
        sig_win = sig * win

        # Zero-pad to the requested length (better frequency-axis interpolation)
        n_pad = n * self.zero_pad
        if self.zero_pad > 1:
            sig_padded = np.zeros(n_pad, dtype=float)
            sig_padded[:n] = sig_win
        else:
            sig_padded = sig_win

        raw = np.fft.rfft(sig_padded)
        freq = np.fft.rfftfreq(n_pad, d=dt)

        # One-sided amplitude with WINDOW COHERENT-GAIN normalisation (IV-E03):
        # divide by Σwin (the coherent gain), NOT by n, so a unit-amplitude
        # tone returns ~1.0 for ANY window (Harris 1978). Backward-compatible:
        # boxcar has Σwin = n. Interior bins ×2 for the one-sided spectrum.
        # (Zero-padding only interpolates the axis; Σwin uses the length-n
        # window, so the calibration is independent of zero_pad.)
        amplitude = np.abs(raw) / win.sum()
        amplitude[1:-1] *= 2.0

        # Optional Gaussian smoothing in the frequency-bin domain.
        if self.smooth_sigma_bins > 0.0:
            from scipy.ndimage import gaussian_filter1d
            amplitude = gaussian_filter1d(
                amplitude, sigma=self.smooth_sigma_bins, mode="nearest")

        return FourierResult(
            frequency_au=freq,
            amplitude=amplitude,
            power=amplitude ** 2,
            column=column,
            dt_au=dt,
            window=self.window,
            t_start_au=self.t_start_au,
            subtract=self.subtract,
        )

    # ------------------------------------------------------------------
    # Convenience methods
    # ------------------------------------------------------------------

    def transform_column(
        self,
        observables: pd.DataFrame,
        column: str,
    ) -> FourierResult:
        """Transform a named column from a loaded observables DataFrame."""
        if "time_au" not in observables.columns:
            raise ValueError("DataFrame must contain a 'time_au' column.")
        if column not in observables.columns:
            raise ValueError(
                f"Column '{column}' not found. Available: {list(observables.columns)}"
            )
        return self.transform(
            observables["time_au"].values, observables[column].values, column=column
        )

    def transform_energy(self, observables: pd.DataFrame) -> FourierResult:
        """Transform the total energy (energy_total, or hartree+xc sum as fallback)."""
        if "energy_total" in observables.columns:
            return self.transform_column(observables, "energy_total")
        if "energy_hartree" in observables.columns and "energy_xc" in observables.columns:
            combined = observables["energy_hartree"].values + observables["energy_xc"].values
            return self.transform(
                observables["time_au"].values, combined, column="energy_hartree+xc"
            )
        raise ValueError(
            "No energy column found (need energy_total, or energy_hartree+energy_xc)."
        )

    def transform_current(
        self,
        observables: pd.DataFrame,
        component: str = "x",
    ) -> FourierResult:
        """Transform a current component ('x', 'y', or 'z')."""
        return self.transform_column(observables, f"current_{component}")

    def transform_dipole(
        self,
        observables: pd.DataFrame,
        component: str = "x",
    ) -> FourierResult:
        """Transform a dipole component ('x', 'y', or 'z')."""
        return self.transform_column(observables, f"dipole_{component}")
