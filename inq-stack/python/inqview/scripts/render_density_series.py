#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from inqview.io.data import SimulationData
from inqview.visualisation.paraview import AnimationSpec, ParaViewPipeline, VolumeRenderSpec


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a real-field raw/meta series to VTI and render PNG frames with ParaView."
    )

    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path("."),
        help="Root directory of the simulation results.",
    )
    parser.add_argument(
        "--series-dir",
        type=str,
        required=True,
        help="Directory containing the raw/meta density series, relative to results-root "
             "(example: results/density or results/orbital_density).",
    )
    parser.add_argument(
        "--vti-dir",
        type=Path,
        required=True,
        help="Directory to write the VTI series into.",
    )
    parser.add_argument(
        "--frames-dir",
        type=Path,
        required=True,
        help="Directory to write PNG frames into.",
    )
    parser.add_argument(
        "--array-name",
        type=str,
        default="density",
        help="Array name to store in the VTI and render in ParaView.",
    )
    parser.add_argument(
        "--pv-exe",
        type=Path,
        default=None,
        help="Path to pvbatch/pvpython, or to a ParaView bin directory.",
    )
    parser.add_argument(
        "--frame-stride",
        type=int,
        default=1,
        help="Render every Nth VTI frame.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1600,
        help="Output image width in pixels.",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=1200,
        help="Output image height in pixels.",
    )
    parser.add_argument(
        "--gif",
        type=Path,
        default=None,
        help="Optional GIF output path.",
    )
    parser.add_argument(
        "--mp4",
        type=Path,
        default=None,
        help="Optional MP4 output path.",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=12,
        help="FPS for GIF/MP4 output.",
    )
    parser.add_argument(
        "--opengl-backend",
        type=str,
        default=None,
        choices=["GLX", "EGL", "OSMesa"],
        help="Optional ParaView OpenGL backend.",
    )

    args = parser.parse_args()

    sim = SimulationData(args.results_root)
    series = sim.field_series(args.series_dir)

    pv = ParaViewPipeline(
        pv_executable=args.pv_exe,
        prefer_pvbatch=True,
        force_offscreen=True,
        opengl_backend=args.opengl_backend,
    )

    render = VolumeRenderSpec(array_name=args.array_name)
    animation = AnimationSpec(
        output_frames_dir=args.frames_dir,
        image_size=(args.width, args.height),
        frame_stride=args.frame_stride,
        filename_prefix="frame",
    )

    frames = pv.render_density_from_meta_series(
        series=series,
        vti_output_dir=args.vti_dir,
        render=render,
        animation=animation,
    )

    print(f"Wrote {len(frames)} PNG frames to {args.frames_dir}")

    if args.gif is not None:
        gif_path = pv.build_gif(args.frames_dir, args.gif, fps=args.fps)
        print(f"Wrote GIF: {gif_path}")

    if args.mp4 is not None:
        mp4_path = pv.build_mp4(args.frames_dir, args.mp4, fps=args.fps)
        print(f"Wrote MP4: {mp4_path}")


if __name__ == "__main__":
    main()
