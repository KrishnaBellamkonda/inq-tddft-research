"""Synthetic analytic signals with known spectra (test ground truth).

Importable helper shared by ``conftest.py`` and the analysis-kernel tests.
Pure numpy; expected spectra are closed-form so tests assert against them
directly (anti-circularity, IV-M10).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Tone:
    """A pure real tone sampled on a uniform time grid, with an on-bin freq.

    ``f0`` is placed exactly on an FFT bin (``f0 = n_cycles / (n * dt)``) so the
    one-sided amplitude lands in a single bin with no scalloping loss — the
    cleanest case for asserting recovered amplitude == ``A`` and peak == ``f0``.
    """

    t: np.ndarray          # sample times (a.u.)
    y: np.ndarray          # signal samples
    A: float               # tone amplitude
    f0: float              # tone frequency (cycles per a.u.) — exactly on a bin
    dt: float              # sample spacing (a.u.)
    n: int                 # number of samples


def make_tone(
    *,
    A: float = 2.5,
    n_cycles: int = 20,
    n: int = 256,
    dt: float = 0.5,
    offset: float = 0.0,
    drift: float = 0.0,
) -> Tone:
    """Build ``A*cos(2π f0 t) + offset + drift*t`` with ``f0`` exactly on a bin.

    offset : constant baseline (DC) added — used to probe subtraction/detrend.
    drift  : linear-trend slope added — used to probe detrend vs initial/mean.
    """
    t = np.arange(n) * dt
    f0 = n_cycles / (n * dt)                # on-bin by construction
    y = A * np.cos(2.0 * np.pi * f0 * t) + offset + drift * t
    return Tone(t=t, y=y, A=A, f0=f0, dt=dt, n=n)
