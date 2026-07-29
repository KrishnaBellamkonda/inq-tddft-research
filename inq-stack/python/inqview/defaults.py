"""
defaults.py — one-call pipelines for common inqview workflows.

Wraps the full FieldSeries → VTI → ParaView → GIF pipeline into a
single function call with sensible defaults.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from .io.data import FieldSeries
    from .visualisation.paraview import AtomSpec, VolumeRenderSpec

PathLike = Union[str, Path]


def default_density_movie(
    series: "FieldSeries",
    output_dir: PathLike,
    pv_executable: PathLike | None = None,
    fps: int = 12,
    render: "VolumeRenderSpec | None" = None,
    image_size: tuple[int, int] = (1600, 1200),
    atoms: "AtomSpec | None" = None,
    frame_stride: int = 1,
) -> dict[str, Path]:
    """Convert a real-density FieldSeries to a GIF via ParaView.

    Pipeline: FieldSeries -> VTI series -> ParaView PNG frames -> GIF.

    Parameters
    ----------
    series       : FieldSeries pointing to the density frames directory.
    output_dir   : root output directory; subdirs vti/, frames/, gif/ are created.
    pv_executable: path to pvbatch/pvpython. Uses PARAVIEW_EXEC env var if None.
    fps          : frames per second for the output GIF.
    render       : VolumeRenderSpec override. None uses default colour/opacity settings.
    image_size   : (width, height) in pixels.
    atoms        : AtomSpec for optional CPK-sphere overlay.
    frame_stride : render every Nth frame (1 = all frames).

    Returns
    -------
    dict with keys 'gif', 'frames_dir', 'vti_dir'.
    """
    import os
    from .visualisation.vti import convert_real_series_to_vti
    from .visualisation.paraview import AnimationSpec, ParaViewPipeline, VolumeRenderSpec

    output_dir = Path(output_dir)
    vti_dir    = output_dir / "vti"
    frames_dir = output_dir / "frames"
    gif_path   = output_dir / "density.gif"

    vti_result = convert_real_series_to_vti(series, vti_dir)

    if render is None:
        render = VolumeRenderSpec(
            array_name="density",
            scalar_range=(float(vti_result.data_min), float(vti_result.data_max)),
        )
    elif render.scalar_range is None:
        from dataclasses import replace
        render = replace(
            render,
            scalar_range=(float(vti_result.data_min), float(vti_result.data_max)),
        )

    pv_exec = pv_executable or os.environ.get("PARAVIEW_EXEC", "pvbatch")
    pv = ParaViewPipeline(pv_executable=str(pv_exec))

    animation = AnimationSpec(
        output_frames_dir=frames_dir,
        image_size=image_size,
        frame_stride=frame_stride,
    )

    pv.render_vti_series(
        vti_files=vti_result.files,
        render=render,
        animation=animation,
        atoms=atoms,
    )

    pv.build_gif(frames_dir=frames_dir, output_path=gif_path, fps=fps)

    return {"gif": gif_path, "frames_dir": frames_dir, "vti_dir": vti_dir}


def default_wavepacket_movie(
    series: "FieldSeries",
    output_dir: PathLike,
    pv_executable: PathLike | None = None,
    fps: int = 12,
    render: "VolumeRenderSpec | None" = None,
    image_size: tuple[int, int] = (1600, 1200),
    atoms: "AtomSpec | None" = None,
    frame_stride: int = 1,
) -> dict[str, Path]:
    """Convert a wavepacket orbital-density FieldSeries to a GIF via ParaView.

    Identical to default_density_movie() -- provided as a separate entry point
    so callers can distinguish density and WP movies by name.
    """
    return default_density_movie(
        series=series,
        output_dir=output_dir,
        pv_executable=pv_executable,
        fps=fps,
        render=render,
        image_size=image_size,
        atoms=atoms,
        frame_stride=frame_stride,
    )
