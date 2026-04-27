"""
screens.py — load and represent LEED/transmission screen patterns.

Reads .dat files written by inqkit::screens::LeedPatternAccumulator::save().

File format:
  # label=LABEL z=Z_BOHR total_time=T_AU n_accum=N
  # nx=NX ny=NY dx=DX dy=DY origin_x=OX origin_y=OY
  v00 v01 ... v0(NX-1)
  v10 ...
  ...
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

import numpy as np

PathLike = Union[str, Path]


@dataclass
class LeedPattern:
    """Accumulated 2D density pattern from a PlaneScreen.

    data      : shape (ny, nx) — ∑_t ρ(x,y,z,t)·dt  [bohr⁻³·a.u.]
    z_bohr    : screen z-position (bohr)
    label     : screen label from the C++ run
    total_time_au : total integration time (a.u.)
    n_accum   : number of accumulation calls
    nx, ny    : grid dimensions
    dx_bohr, dy_bohr : real-space grid spacing (bohr)
    origin_x_bohr, origin_y_bohr : grid origin (bohr)
    """

    data: np.ndarray
    z_bohr: float
    label: str
    total_time_au: float
    n_accum: int
    nx: int
    ny: int
    dx_bohr: float
    dy_bohr: float
    origin_x_bohr: float
    origin_y_bohr: float

    @property
    def x_axis(self) -> np.ndarray:
        """Physical x-coordinates of grid columns (bohr)."""
        return self.origin_x_bohr + np.arange(self.nx) * self.dx_bohr

    @property
    def y_axis(self) -> np.ndarray:
        """Physical y-coordinates of grid rows (bohr)."""
        return self.origin_y_bohr + np.arange(self.ny) * self.dy_bohr

    @property
    def extent_bohr(self) -> tuple[float, float, float, float]:
        """(x_min, x_max, y_min, y_max) for matplotlib imshow extent."""
        return (
            float(self.x_axis[0]),
            float(self.x_axis[-1] + self.dx_bohr),
            float(self.y_axis[0]),
            float(self.y_axis[-1] + self.dy_bohr),
        )

    def inverse_fft(self, method: str = "patterson",
                    hann: bool = True) -> np.ndarray:
        """Inverse-FFT this LEED screen back to a real-space density estimate.

        Thin wrapper around inqview.postprocess._ifft.reconstruct_real_space.
        See that module's docstring for the available methods
        ('patterson' — Patterson autocorrelation, default; 'amp_only' —
        phase-less amplitude reconstruction).
        """
        from .postprocess._ifft import reconstruct_real_space
        return reconstruct_real_space(self, method=method, hann=hann)


def _parse_header_kv(line: str) -> dict[str, str]:
    """Parse 'key=value key=value ...' from a header comment line."""
    kv: dict[str, str] = {}
    for token in line.lstrip("# ").split():
        if "=" in token:
            k, v = token.split("=", 1)
            kv[k] = v
    return kv


def load_leed_pattern(path: PathLike) -> LeedPattern:
    """Load a .dat file written by LeedPatternAccumulator.save().

    Parameters
    ----------
    path : path to the .dat file.

    Returns
    -------
    LeedPattern
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Screen pattern file not found: {path}")

    with path.open() as fh:
        lines = fh.readlines()

    if len(lines) < 2:
        raise ValueError(f"Screen pattern file too short: {path}")

    # Header line 1: label, z, total_time (accumulated) or t (snapshot), n_accum
    h1 = _parse_header_kv(lines[0])
    label      = h1.get("label", "")
    z_bohr     = float(h1.get("z", "0"))
    total_time = float(h1.get("total_time", h1.get("t", "0")))
    n_accum    = int(h1.get("n_accum", "1"))

    # Header line 2 is optional: snapshot files write only 1 header line.
    # Detect by checking whether line[1] starts with '#'.
    if lines[1].startswith("#"):
        h2 = _parse_header_kv(lines[1])
        nx       = int(h2.get("nx", "0"))
        ny       = int(h2.get("ny", "0"))
        dx_bohr  = float(h2.get("dx", "1"))
        dy_bohr  = float(h2.get("dy", "1"))
        origin_x = float(h2.get("origin_x", "0"))
        origin_y = float(h2.get("origin_y", "0"))
        data_start = 2
    else:
        nx = ny = 0
        dx_bohr = dy_bohr = 1.0
        origin_x = origin_y = 0.0
        data_start = 1

    # Data rows
    data_lines = [l for l in lines[data_start:] if l.strip() and not l.startswith("#")]
    rows = []
    for dl in data_lines:
        rows.append([float(v) for v in dl.split()])

    data = np.array(rows, dtype=np.float64)

    # Infer nx/ny from data if not in header (snapshot format)
    if nx == 0 and data.ndim == 2:
        ny, nx = data.shape
        dx_bohr = dy_bohr = 1.0

    if data.ndim != 2 or data.shape[0] != ny or data.shape[1] != nx:
        raise ValueError(
            f"Unexpected data shape {data.shape} for nx={nx}, ny={ny} in {path}"
        )

    # ── FFT-shift to centre the diffraction peak ──────────────────────────
    # LeedPatternAccumulator (inq-stack/include/inqkit/screens/...) writes the
    # screen in INQ's FFT-natural order: array index (0, 0) = physical origin
    # x = 0, y = 0; then positive coordinates first, then wrapped negative.
    # Without np.fft.fftshift the diffraction peak lands at a corner with its
    # 4-way symmetric tails distributed to the other three corners (the
    # "four-corner-split" failure mode the spec §17.6 warns about).
    #
    # Reference correct loader (the proven path):
    # ResearchProject/systems/coronene/run_propagate_paper_replica/analysis.py
    # `_load_screen_centred`.
    data = np.fft.fftshift(data)

    # The C++ writer emits origin_x = origin_y = 0 (assuming raw-index
    # convention). After fftshift, array index (0, 0) corresponds to physical
    # (-Lx/2, -Ly/2). Override the origin so LeedPattern.extent_bohr spans
    # [-Lx/2, +Lx/2, -Ly/2, +Ly/2] automatically.
    origin_x = -0.5 * nx * dx_bohr
    origin_y = -0.5 * ny * dy_bohr

    return LeedPattern(
        data=data,
        z_bohr=z_bohr,
        label=label,
        total_time_au=total_time,
        n_accum=n_accum,
        nx=nx,
        ny=ny,
        dx_bohr=dx_bohr,
        dy_bohr=dy_bohr,
        origin_x_bohr=origin_x,
        origin_y_bohr=origin_y,
    )
