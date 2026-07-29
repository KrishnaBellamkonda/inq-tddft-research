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


# ──────────────────────────────────────────────────────────────────────────
# Post-IFW shading (campaign-universal helper, plan rule §4)
#
# The interference-free window (IFW) is the part of a WP trajectory during
# which the WP's far-face Gaussian tail has not yet hit the box boundary;
# beyond t_IFW the response of the bath is contaminated by self-image
# scattering and every quantitative claim should explicitly *exclude* the
# post-IFW segment. Plots that include both regions shade the post-IFW one
# in light grey so the reader can see the unshaded region is the "trusted"
# fit window.
#
# Geometry, from ``shared/configs/boundary_rule.hpp`` (both standard and
# relaxed rules agree on ifw_end_z; only stop_z differs trivially):
#   ifw_end_z = +L/2 - 3 sigma          (Gaussian-3sigma-at-far-face)
#   stop_z    = +L/2 -   sigma          (centroid stops 1 sigma from face)
#   t_IFW     = (ifw_end_z - launch_z) / |v|
#   t_total   = (stop_z    - launch_z) / |v|
# where |v| = |k_0| in atomic units (m_e = 1).
# ──────────────────────────────────────────────────────────────────────────


def post_ifw_window_au(*, launch_z_bohr: float, l_bohr: float,
                        sigma_bohr: float, v_au: float) -> tuple[float, float]:
    """Compute (t_IFW, t_total) in atomic time units.

    Inputs are all in atomic units. ``v_au`` is the projectile speed along
    the launch axis (i.e. |k_0| in Bohr^-1 — same number as the velocity in
    Bohr / a.u. because m_e = 1).

    Raises ``ValueError`` if ``v_au <= 0`` or the box geometry leaves no
    IFW region.
    """
    if v_au <= 0:
        raise ValueError(f"post_ifw_window_au: v_au must be > 0, got {v_au}")
    ifw_end_z = 0.5 * l_bohr - 3.0 * sigma_bohr
    stop_z    = 0.5 * l_bohr -       sigma_bohr
    t_ifw   = (ifw_end_z - launch_z_bohr) / v_au
    t_total = (stop_z    - launch_z_bohr) / v_au
    if not (t_total > t_ifw >= 0):
        # E.g. launch_z too close to the far face (relaxed-rule σ=8 case).
        # Caller decides whether to skip or shade everything.
        pass
    return t_ifw, t_total


def post_ifw_window_from_summary(results_dir: Path,
                                  ) -> tuple[float, float] | None:
    """Read ``run_summary.txt`` and derive (t_IFW, t_total) in a.u.

    Returns ``None`` if any required field is missing. Looks for:

    * ``wp_center_bohr  = x y z``    → launch_z = z component
    * ``wp_sigma_bohr   = <float>``
    * ``wp_k0_bohr_inv  = kx ky kz`` → v_au = sqrt(kx^2 + ky^2 + kz^2)
    * ``cell_bohr       = L^3 (...)`` → L
    """
    import re

    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
    text = rs.read_text()

    def _floats(key: str) -> list[float] | None:
        m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(.*)$",
                      text, flags=re.MULTILINE)
        if not m:
            return None
        # Strip parens-suffix (e.g. "50^3 (cubic, periodic)" reduces to "50^3")
        tail = m.group(1).split("(")[0].strip()
        # Numeric tokens: handles "0 0 -10", "5", or "50^3" (split on '^').
        tokens = re.split(r"\s+|\^", tail)
        try:
            return [float(t) for t in tokens if t]
        except ValueError:
            return None

    center = _floats("wp_center_bohr")
    sigma  = _floats("wp_sigma_bohr")
    k0     = _floats("wp_k0_bohr_inv")
    cell   = _floats("cell_bohr")

    if not (center and sigma and k0 and cell):
        return None
    if len(center) < 3 or len(k0) < 3 or not cell:
        return None

    import math
    launch_z = center[2]
    sigma_b  = sigma[0]
    v_au     = math.sqrt(k0[0]**2 + k0[1]**2 + k0[2]**2)
    l_bohr   = cell[0]

    try:
        return post_ifw_window_au(launch_z_bohr=launch_z, l_bohr=l_bohr,
                                  sigma_bohr=sigma_b, v_au=v_au)
    except ValueError:
        return None


def post_ifw_shade(ax, t_ifw_au: float, t_total_au: float, *,
                    color: str = "0.85", alpha: float = 0.5,
                    label: str | None = "post-IFW") -> None:
    """Shade the post-IFW slab of a matplotlib time-axis (x-axis = time_au).

    Idempotent re: the axis state — just calls ``ax.axvspan(...)``. Pass
    ``label=None`` to suppress the legend entry.
    """
    if not (t_total_au > t_ifw_au):
        return
    ax.axvspan(t_ifw_au, t_total_au, color=color, alpha=alpha,
               label=label, zorder=0, lw=0)


def ifw_highlight(ax, t_ifw_au: float, *,
                   color: str = "#fff4cc", alpha: float = 0.55,
                   label: str | None = "IFW (focus region)") -> None:
    """Highlight the interference-free window itself (t in [0, t_IFW]).

    Use this when the IFW is the region of interest and the post-IFW
    region should be visible-but-deemphasised (no shade). The default
    soft-yellow tint makes the IFW pop without obscuring the data.

    Pass ``label=None`` to suppress the legend entry.
    """
    if not (t_ifw_au > 0):
        return
    ax.axvspan(0.0, t_ifw_au, color=color, alpha=alpha,
               label=label, zorder=0, lw=0)
