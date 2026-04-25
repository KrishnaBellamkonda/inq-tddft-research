"""
Loader and plotter for KS overlap matrix snapshots written by
inqkit::observables::OrbitalOverlapMatrix.

Schema (per docs/observables_reference.md):
    O_ij(t) = | dV * sum_r conj(psi_i^GS(r)) * psi_j(r, t) |^2
    rows i in [0, n_ref): GS reference orbitals stored at t=0
    cols j in [0, n_evolved): evolved orbitals at t (col n_ref-1+1 is WP slot)

The C++ driver passes ``n_ref = wp_idx`` (the index of the slot the WP was
injected into), which means *n_ref includes any unoccupied extra-state
ground orbitals*. For visualisation, callers typically want to filter to
the *meaningful* submatrix: rows = occupied orbitals only, columns =
occupied + WP slot.

Helpers:
    iter_overlap_series(overlap_dir)
    pick_meaningful_columns(n_occupied, wp_idx, n_evolved)
    plot_overlap_column_gif(...)
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import imageio.v3 as iio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


@dataclass
class OverlapSnapshot:
    step: int
    time_au: float
    matrix: np.ndarray  # shape (n_ref, n_evolved)
    file: Path


def load_overlap_csv(path: Path) -> np.ndarray | None:
    """
    Load one ``overlap_NNNNNN.csv`` written by the C++ writer. Returns the
    matrix as a 2D ``np.ndarray`` of shape ``(n_ref, n_evolved)`` or None
    if the file cannot be parsed.
    """
    if not path.exists():
        return None
    rows: list[list[float]] = []
    n_evolved = None
    with path.open() as fh:
        for raw in fh:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                vals = [float(v) for v in stripped.split(",")]
            except ValueError:
                # malformed row -> abort, signal failure
                return None
            if n_evolved is None:
                n_evolved = len(vals)
            elif len(vals) != n_evolved:
                # ragged row -> abort, signal failure
                return None
            rows.append(vals)
    if not rows:
        return None
    return np.array(rows, dtype=float)


def iter_overlap_series(overlap_dir: Path) -> Iterator[OverlapSnapshot]:
    """Iterate snapshots in time order using ``index.csv``."""
    index_csv = overlap_dir / "index.csv"
    if not index_csv.exists():
        return
    expected_shape: tuple[int, int] | None = None
    with index_csv.open() as fh:
        for row in csv.DictReader(fh):
            try:
                step = int(row["step"])
                t = float(row["time_au"])
            except (KeyError, ValueError):
                continue
            file_name = row.get("file")
            if not file_name:
                continue
            mat_path = overlap_dir / file_name
            mat = load_overlap_csv(mat_path)
            if mat is None:
                continue
            if expected_shape is None:
                expected_shape = mat.shape
            elif mat.shape != expected_shape:
                # silently dropping inconsistent shapes is what bit the
                # previous GIF code; here we surface them via stderr.
                import sys
                print(
                    f"WARN: overlap snapshot {mat_path.name} has shape {mat.shape}, "
                    f"expected {expected_shape}; skipping",
                    file=sys.stderr,
                )
                continue
            yield OverlapSnapshot(step=step, time_au=t, matrix=mat, file=mat_path)


def pick_meaningful_columns(
    n_occupied: int, wp_idx: int, n_evolved: int
) -> list[tuple[int, int]]:
    """
    Decide which columns of the raw overlap matrix are physically
    meaningful for visualisation.

    Returns a list of ``(j_global, j_label)`` tuples where:
        j_global is the index into the matrix column dimension
        j_label is the human-friendly orbital index used in filenames.

    Selected columns are:
        the n_occupied occupied orbitals (j_global in [0, n_occupied)),
        plus the WP slot (j_global = wp_idx).

    Raises ValueError if wp_idx is outside the column range.
    """
    if wp_idx >= n_evolved or wp_idx < 0:
        raise ValueError(
            f"wp_idx={wp_idx} is out of range for n_evolved={n_evolved}"
        )
    out: list[tuple[int, int]] = []
    for j in range(min(n_occupied, n_evolved)):
        out.append((j, j))
    if 0 <= wp_idx < n_evolved and wp_idx >= n_occupied:
        out.append((wp_idx, wp_idx))
    return out


def plot_overlap_column_gif(
    out_path: Path,
    n_ref_rows: int,
    col_data: list[np.ndarray],
    times: list[float],
    title_prefix: str = "Overlap",
    fps: int = 12,
) -> Path:
    """
    Build a GIF of bar charts: one frame per timestep, x = GS orbital
    index i in [0, n_ref_rows), y = ``|O_ij|^2`` for the requested column j.
    """
    if not col_data:
        return out_path
    if any(c.shape != (n_ref_rows,) for c in col_data):
        raise ValueError(
            f"col_data has inconsistent row count; expected each entry of shape "
            f"({n_ref_rows},) but got {[c.shape for c in col_data]}"
        )
    ymax = max(float(c.max()) for c in col_data) * 1.05 + 1e-9
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    for t, col in zip(times, col_data):
        fig, ax = plt.subplots(figsize=(6, 3), dpi=80)
        ax.bar(np.arange(n_ref_rows), col, color="steelblue", edgecolor="none")
        ax.set_xlim(-0.5, n_ref_rows - 0.5)
        ax.set_ylim(0, ymax)
        ax.set_xlabel("GS orbital index i")
        ax.set_ylabel(r"$|O_{ij}|^2$")
        ax.set_title(f"{title_prefix}  t={t:.2f} a.u.")
        fig.tight_layout()
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        frames.append(frame)
        plt.close(fig)
    iio.imwrite(out_path, frames, fps=fps, loop=0)
    return out_path
