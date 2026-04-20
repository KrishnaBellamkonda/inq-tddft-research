"""
fourier.py — FFT-based post-processing for inqview observables.

Provides FourierTransform (class-based) with tunable windowing, detrending,
and convenience methods for energy, current, and dipole columns.
"""

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


class FourierTransform:
    """Windowed FFT of TDDFT time-series observables.

    Parameters
    ----------
    window : WindowSpec
        Window function and its parameters. Default: Hann window.
    detrend : bool
        If True, subtract a linear trend before transforming (removes DC drift).
    """

    def __init__(
        self,
        window: WindowSpec | None = None,
        detrend: bool = True,
    ) -> None:
        self.window = window if window is not None else WindowSpec("hann")
        self.detrend = detrend

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

        n = len(time_au)
        dt = float(time_au[1] - time_au[0])

        sig = scipy_detrend(values, type="linear") if self.detrend else values.copy()
        win = self.window.build(n)
        sig_win = sig * win

        raw = np.fft.rfft(sig_win)
        freq = np.fft.rfftfreq(n, d=dt)

        # One-sided normalisation: divide by N; double non-DC and non-Nyquist bins
        amplitude = np.abs(raw) / n
        amplitude[1:-1] *= 2.0

        return FourierResult(
            frequency_au=freq,
            amplitude=amplitude,
            power=amplitude ** 2,
            column=column,
            dt_au=dt,
            window=self.window,
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
