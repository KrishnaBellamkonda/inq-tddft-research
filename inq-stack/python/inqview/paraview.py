from __future__ import annotations

from dataclasses import dataclass, field, replace
import json
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
from typing import Iterable

from .data import FieldSeries
from .vti import convert_real_series_to_vti

_ANG_TO_BOHR = 1.8897259886

# CPK colours (RGB 0–1), standard chemistry convention
_CPK_COLORS: dict[str, tuple[float, float, float]] = {
    "H":  (1.000, 1.000, 1.000),
    "C":  (0.565, 0.565, 0.565),
    "N":  (0.188, 0.314, 0.973),
    "O":  (1.000, 0.051, 0.051),
    "F":  (0.565, 0.878, 0.314),
    "S":  (1.000, 1.000, 0.188),
    "Cl": (0.122, 0.941, 0.122),
    "Br": (0.647, 0.161, 0.161),
    "I":  (0.580, 0.000, 0.580),
    "P":  (1.000, 0.502, 0.000),
    "Si": (0.941, 0.784, 0.627),
    "Na": (0.671, 0.361, 0.949),
    "Mg": (0.541, 1.000, 0.000),
    "Ca": (0.239, 1.000, 0.000),
}
_DEFAULT_CPK_COLOR: tuple[float, float, float] = (0.667, 0.667, 0.667)

# Van der Waals radii in Å (used to set sphere size in ParaView)
_VDW_RADII_ANG: dict[str, float] = {
    "H":  1.20,
    "C":  1.70,
    "N":  1.55,
    "O":  1.52,
    "F":  1.47,
    "S":  1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I":  1.98,
    "P":  1.80,
    "Si": 2.10,
    "Na": 2.27,
    "Mg": 1.73,
    "Ca": 2.31,
}
_DEFAULT_VDW_ANG: float = 1.70


@dataclass
class AtomSpec:
    """
    Atom positions and species for rendering as CPK-coloured spheres in ParaView.

    Positions must be in bohr, in the same Cartesian frame as the VTI grid
    (i.e. origin at the (0,0,0) corner of the simulation box).
    """
    positions: list[list[float]]    # shape (N, 3), bohr
    symbols:   list[str]            # element symbols, length N
    radius_scale:   float = 0.4    # fraction of VDW radius used for sphere size
    opacity:        float = 1.0
    specular:       float = 0.3
    specular_power: float = 20.0


@dataclass
class VolumeRenderSpec:
    array_name: str
    scalar_range: tuple[float, float] | None = None
    color_preset: str | None = "Cividis (matplotlib)"
    show_scalar_bar: bool = False

    # Camera / view
    camera_azimuth_deg: float = 35.0
    camera_elevation_deg: float = 25.0

    # Background
    background_rgb: tuple[float, float, float] = (1.0, 1.0, 1.0)

    # Simple opacity control points as (scalar, opacity)
    # If None, a default ramp is constructed from scalar_range.
    opacity_points: list[tuple[float, float]] | None = None


@dataclass
class AnimationSpec:
    output_frames_dir: Path
    image_size: tuple[int, int] = (1600, 1200)
    frame_stride: int = 1
    filename_prefix: str = "frame"


class ParaViewPipeline:
    """
    Minimal ParaView batch-render wrapper for VTI series.
    """

    def __init__(
        self,
        pv_executable: str | Path | None = None,
        prefer_pvbatch: bool = True,
        force_offscreen: bool = True,
        opengl_backend: str | None = None,
    ):
        self.prefer_pvbatch = prefer_pvbatch
        self.force_offscreen = force_offscreen
        self.opengl_backend = opengl_backend
        self.pv_executable = self._find_paraview_executable(pv_executable)

    def render_density_from_meta_series(
        self,
        series: FieldSeries | Iterable[str | Path],
        vti_output_dir: str | Path,
        render: VolumeRenderSpec,
        animation: AnimationSpec,
        atoms: AtomSpec | None = None,
    ) -> list[Path]:
        """
        Convenience pipeline:
          raw/meta -> VTI series -> ParaView PNG frames

        If *atoms* is provided, CPK-coloured spheres are rendered at the
        given positions alongside the volume density.
        """
        vti_result = convert_real_series_to_vti(
            series=series,
            output_dir=vti_output_dir,
            array_name=render.array_name,
        )

        effective_render = render
        if effective_render.scalar_range is None:
            effective_render = replace(
                effective_render,
                scalar_range=(vti_result.data_min, vti_result.data_max),
            )

        return self.render_vti_series(
            vti_files=vti_result.files,
            render=effective_render,
            animation=animation,
            atoms=atoms,
        )

    def render_vti_series(
        self,
        vti_files: Iterable[str | Path],
        render: VolumeRenderSpec,
        animation: AnimationSpec,
        atoms: AtomSpec | None = None,
    ) -> list[Path]:
        """
        Render a VTI series to a directory of PNG frames using ParaView batch Python.
        """
        files = [Path(p).resolve() for p in vti_files]
        if not files:
            raise RuntimeError("No VTI files were provided for rendering.")

        animation.output_frames_dir.mkdir(parents=True, exist_ok=True)

        scalar_range = render.scalar_range
        if scalar_range is None:
            raise RuntimeError(
                "render_vti_series requires render.scalar_range to be set. "
                "Use render_density_from_meta_series(...) for automatic range derivation."
            )

        opacity_points = render.opacity_points
        if opacity_points is None:
            opacity_points = _default_opacity_points(*scalar_range)

        config: dict = {
            "files": [str(p) for p in files],
            "array_name": render.array_name,
            "scalar_range": [float(scalar_range[0]), float(scalar_range[1])],
            "color_preset": render.color_preset,
            "show_scalar_bar": bool(render.show_scalar_bar),
            "camera_azimuth_deg": float(render.camera_azimuth_deg),
            "camera_elevation_deg": float(render.camera_elevation_deg),
            "background_rgb": [float(x) for x in render.background_rgb],
            "opacity_points": [[float(x), float(y)] for x, y in opacity_points],
            "output_frames_dir": str(animation.output_frames_dir.resolve()),
            "image_size": [int(animation.image_size[0]), int(animation.image_size[1])],
            "frame_stride": int(animation.frame_stride),
            "filename_prefix": animation.filename_prefix,
        }

        if atoms is not None:
            config["atoms"] = {
                "positions": [[float(x) for x in pos] for pos in atoms.positions],
                "symbols": list(atoms.symbols),
                "radius_scale": float(atoms.radius_scale),
                "opacity": float(atoms.opacity),
                "specular": float(atoms.specular),
                "specular_power": float(atoms.specular_power),
                "cpk_colors": {k: list(v) for k, v in _CPK_COLORS.items()},
                "vdw_radii_ang": dict(_VDW_RADII_ANG),
                "default_cpk_color": list(_DEFAULT_CPK_COLOR),
                "default_vdw_ang": _DEFAULT_VDW_ANG,
                "ang_to_bohr": _ANG_TO_BOHR,
            }

        with tempfile.TemporaryDirectory(prefix="inqview_paraview_") as tmpdir_str:
            tmpdir = Path(tmpdir_str)
            cfg_path = tmpdir / "config.json"
            script_path = tmpdir / "render_vti_series.py"

            cfg_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
            script_path.write_text(_paraview_batch_script(), encoding="utf-8")

            cmd = [str(self.pv_executable)]
            if self.force_offscreen:
                cmd.append("--force-offscreen-rendering")
            if self.opengl_backend:
                cmd.extend(["--opengl-window-backend", self.opengl_backend])

            cmd.extend([str(script_path), str(cfg_path)])

            subprocess.run(cmd, check=True)

        return sorted(animation.output_frames_dir.glob(f"{animation.filename_prefix}_*.png"))

    def build_gif(
        self,
        frames_dir: str | Path,
        output_path: str | Path,
        fps: int = 12,
    ) -> Path:
        """
        Build a GIF from rendered frames.

        Requires imageio + Pillow support in the current Python environment.
        """
        try:
            import imageio.v3 as iio
        except ImportError as exc:
            raise RuntimeError(
                "build_gif requires imageio. Install it with: python -m pip install imageio"
            ) from exc

        frames_dir = Path(frames_dir)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        frame_paths = sorted(frames_dir.glob("*.png"))
        if not frame_paths:
            raise RuntimeError(f"No PNG frames found in {frames_dir}")

        images = [iio.imread(path) for path in frame_paths]
        duration_ms = int(round(1000 / fps))
        iio.imwrite(output_path, images, duration=duration_ms, loop=0)

        return output_path

    def build_mp4(
        self,
        frames_dir: str | Path,
        output_path: str | Path,
        fps: int = 12,
    ) -> Path:
        """
        Build an MP4 from rendered frames.

        Requires imageio and an ffmpeg backend available to imageio.
        """
        try:
            import imageio.v3 as iio
        except ImportError as exc:
            raise RuntimeError(
                "build_mp4 requires imageio. Install it with: "
                "python -m pip install imageio imageio-ffmpeg"
            ) from exc

        frames_dir = Path(frames_dir)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        frame_paths = sorted(frames_dir.glob("*.png"))
        if not frame_paths:
            raise RuntimeError(f"No PNG frames found in {frames_dir}")

        images = [iio.imread(path) for path in frame_paths]
        iio.imwrite(output_path, images, fps=fps)

        return output_path

    def _find_paraview_executable(self, explicit: str | Path | None) -> Path:
        if explicit is not None:
            p = Path(explicit).expanduser().resolve()
            if p.is_dir():
                candidate = self._pick_from_bin_dir(p)
                if candidate is None:
                    raise RuntimeError(f"No pvbatch/pvpython found in directory: {p}")
                return candidate
            if p.exists():
                return p
            raise RuntimeError(f"ParaView executable path does not exist: {p}")

        env_candidates = []
        for key in ("PARAVIEW_PVBATCH", "PARAVIEW_PVPYTHON", "PARAVIEW_BIN_DIR"):
            value = os.environ.get(key)
            if value:
                env_candidates.append(Path(value).expanduser())

        for path in env_candidates:
            if path.is_dir():
                candidate = self._pick_from_bin_dir(path)
                if candidate is not None:
                    return candidate
            elif path.exists():
                return path.resolve()

        names = ["pvbatch", "pvpython"] if self.prefer_pvbatch else ["pvpython", "pvbatch"]
        for name in names:
            found = shutil.which(name)
            if found:
                return Path(found).resolve()

        raise RuntimeError(
            "Could not find a ParaView batch executable. Provide pv_executable=..., "
            "or set PARAVIEW_PVBATCH / PARAVIEW_PVPYTHON / PARAVIEW_BIN_DIR."
        )

    def _pick_from_bin_dir(self, directory: Path) -> Path | None:
        candidates = ["pvbatch", "pvpython"] if self.prefer_pvbatch else ["pvpython", "pvbatch"]
        for name in candidates:
            p = directory / name
            if p.exists():
                return p.resolve()
        return None


def _default_opacity_points(vmin: float, vmax: float) -> list[tuple[float, float]]:
    if vmax <= vmin:
        return [(vmin, 0.0), (vmin + 1e-12, 0.2)]

    span = vmax - vmin
    return [
        (vmin, 0.0),
        (vmin + 0.03 * span, 0.0),
        (vmin + 0.15 * span, 0.03),
        (vmin + 0.50 * span, 0.12),
        (vmax, 0.35),
    ]


def _paraview_batch_script() -> str:
    return r'''
from __future__ import annotations

import json
from pathlib import Path
import sys

from paraview.simple import *  # noqa: F401,F403


def main() -> None:
    cfg_path = Path(sys.argv[1])
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))

    files = [Path(p) for p in cfg["files"]]
    array_name = cfg["array_name"]
    scalar_min, scalar_max = cfg["scalar_range"]
    color_preset = cfg["color_preset"]
    show_scalar_bar = cfg["show_scalar_bar"]
    azimuth_deg = cfg["camera_azimuth_deg"]
    elevation_deg = cfg["camera_elevation_deg"]
    background_rgb = cfg["background_rgb"]
    opacity_points = cfg["opacity_points"]
    output_frames_dir = Path(cfg["output_frames_dir"])
    image_size = cfg["image_size"]
    frame_stride = int(cfg["frame_stride"])
    filename_prefix = cfg["filename_prefix"]

    output_frames_dir.mkdir(parents=True, exist_ok=True)

    render_view = GetActiveViewOrCreate("RenderView")
    render_view.ViewSize = image_size
    render_view.Background = background_rgb
    render_view.OrientationAxesVisibility = 0

    # Create atom spheres once; they persist across all frames.
    if "atoms" in cfg:
        atom_data = cfg["atoms"]
        ang_to_bohr = float(atom_data["ang_to_bohr"])
        radius_scale = float(atom_data["radius_scale"])
        for sym, pos in zip(atom_data["symbols"], atom_data["positions"]):
            color = atom_data["cpk_colors"].get(sym, atom_data["default_cpk_color"])
            r_ang = float(atom_data["vdw_radii_ang"].get(sym, atom_data["default_vdw_ang"]))
            r_bohr = r_ang * ang_to_bohr * radius_scale
            sphere = Sphere()
            sphere.Center = [float(c) for c in pos]
            sphere.Radius = float(r_bohr)
            sphere.ThetaResolution = 24
            sphere.PhiResolution = 24
            d = Show(sphere, render_view)
            d.ColorArrayName = ['POINTS', '']   # solid color, no scalar array
            d.DiffuseColor = [float(c) for c in color]
            d.Specular = float(atom_data["specular"])
            d.SpecularPower = float(atom_data["specular_power"])
            d.Opacity = float(atom_data["opacity"])

    did_adjust_camera = False
    frame_counter = 0

    for _, vti_path in enumerate(files[::frame_stride]):
        reader = XMLImageDataReader(FileName=[str(vti_path)])
        try:
            reader.PointArrayStatus = [array_name]
        except Exception:
            pass

        display = Show(reader, render_view)
        display.Representation = "Volume"

        ColorBy(display, ("POINTS", array_name))

        lut = GetColorTransferFunction(array_name)
        pwf = GetOpacityTransferFunction(array_name)

        try:
            if color_preset:
                lut.ApplyPreset(color_preset, True)
        except Exception:
            pass

        lut.RescaleTransferFunction(scalar_min, scalar_max)
        pwf.RescaleTransferFunction(scalar_min, scalar_max)

        flat_points = []
        for x, y in opacity_points:
            flat_points.extend([float(x), float(y), 0.5, 0.0])
        pwf.Points = flat_points

        try:
            display.SetScalarBarVisibility(render_view, show_scalar_bar)
        except Exception:
            pass

        Render(render_view)

        if not did_adjust_camera:
            try:
                render_view.ResetCamera()
                cam = GetActiveCamera()
                cam.Azimuth(float(azimuth_deg))
                cam.Elevation(float(elevation_deg))
                render_view.ResetCameraClippingRange()
            except Exception:
                pass
            did_adjust_camera = True
            Render(render_view)

        out_path = output_frames_dir / f"{filename_prefix}_{frame_counter:06d}.png"
        SaveScreenshot(str(out_path), render_view, ImageResolution=image_size)

        Hide(reader, render_view)
        Delete(reader)

        frame_counter += 1


if __name__ == "__main__":
    main()
'''
