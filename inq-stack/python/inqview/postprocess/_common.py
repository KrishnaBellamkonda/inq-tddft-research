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
          total_steps: int | None = None, time_au: float | None = None,
          multiline: bool = True) -> str:
    """Plot title.

    Two modes:

    * ``multiline=True`` (default for animations): line 1 = ``run_name: what``;
      line 2 = ``step k/N, t = X.XX fs`` (only those tokens that are non-None
      appear). Matches the visualisation rule (TODO 1a).
    * ``multiline=False``: legacy single-line ``run_name: what, step ..., t = ...``.
    """
    head = f"{run_name}: {what}"
    sub_parts: list[str] = []
    if step is not None and total_steps is not None:
        sub_parts.append(f"step {step:d}/{total_steps:d}")
    if time_au is not None:
        sub_parts.append(f"t = {sigfigs(fs_from_au(time_au))} fs")
    sub = ", ".join(sub_parts)
    if not sub:
        return head
    return f"{head}\n{sub}" if multiline else f"{head}, {sub}"


def write_animation(out_stem: Path, png_paths: list[Path], *,
                    fps: int = 8, also_mp4: bool = True) -> dict[str, Path]:
    """Render an animation from ``png_paths`` to both GIF and MP4.

    Returns the dict ``{"gif": Path, "mp4": Path | None}``. MP4 is silently
    skipped if ``ffmpeg`` is not on PATH or imageio's libx264 plugin is
    missing.

    Each PNG file becomes one frame at the given fps. ``out_stem`` is the
    path *without* extension (the helper appends ``.gif`` and ``.mp4``).
    """
    import imageio.v2 as imageio
    import shutil
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)

    gif_path = out_stem.with_suffix(".gif")
    with imageio.get_writer(gif_path, mode="I", fps=fps, loop=0) as wr:
        for p in png_paths:
            wr.append_data(imageio.imread(p))

    mp4_path: Path | None = None
    if also_mp4 and shutil.which("ffmpeg"):
        try:
            mp4_path = out_stem.with_suffix(".mp4")
            with imageio.get_writer(
                mp4_path, fps=fps, codec="libx264", quality=8,
                pixelformat="yuv420p", macro_block_size=1
            ) as wr:
                for p in png_paths:
                    wr.append_data(imageio.imread(p))
        except Exception:
            mp4_path = None

    return {"gif": gif_path, "mp4": mp4_path}


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
