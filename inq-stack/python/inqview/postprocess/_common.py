"""Shared helpers for the postprocess phases."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def sigfigs(x: float, n: int = 3) -> str:
    """Round x to n significant figures and format without trailing junk."""
    if not np.isfinite(x):
        return str(x)
    if x == 0:
        return "0"
    from math import floor, log10

    d = n - int(floor(log10(abs(x)))) - 1
    return f"{round(x, d):g}"


def fs_from_au(t_au: float) -> float:
    return t_au / 41.341374575751


def title(run_name: str, what: str, *, step: int | None = None,
          total_steps: int | None = None, time_au: float | None = None) -> str:
    """Standard plot title: ``run_name: what [, step k/N, t = X.XX fs]``."""
    parts = [f"{run_name}: {what}"]
    if step is not None and total_steps is not None:
        parts.append(f"step {step:d}/{total_steps:d}")
    if time_au is not None:
        parts.append(f"t = {sigfigs(fs_from_au(time_au))} fs")
    return ", ".join(parts)


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def need_rebuild(out: Path, rebuild: bool) -> bool:
    return rebuild or not out.exists()


def list_vti_series(directory: Path, prefix: str | None = None) -> list[Path]:
    """Return sorted list of ``<directory>/*_tNNNNNN.vti`` files.

    The ``prefix`` argument is kept for backwards compat; the C++ writer uses
    the field_name (e.g. ``density``), not the directory name, so trying to
    match ``{cat}_t*.vti`` would miss every file. We glob ``*_t*.vti`` and
    require the ``_tNNNNNN`` suffix via the trailing ``_t*.vti`` pattern.
    """
    if not directory.exists():
        return []
    return list(sorted(directory.glob("*_t*.vti")))


def list_screen_files(directory: Path) -> list[Path]:
    if not directory.exists():
        return []
    return sorted(directory.glob("screen_*.dat"))
