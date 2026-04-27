"""Phase: ``paraview_3d`` — overlay (system + WP) volume rendering, two cameras.

Reads ``raw/vti/density_rt_system/*.vti`` and ``raw/vti/density_rt_wp/*.vti``,
renders both volumes in the same ParaView scene with:

* log-scale colour LUT per volume,
* opacity tied to the density value (piecewise linear with three control points),
* two camera setups:
  - **view_headon**: camera at the −L_z/2 face looking toward +z; the WP
    starts deep in the scene and approaches the eye over time.
  - **view_3q**: yaw +45° / elevate +25° from view_headon, slight dolly out.

Outputs (per run, under ``results/analysis/density/paraview_3d/``):

* ``volume_overlay_view_headon.{gif,mp4}``
* ``volume_overlay_view_3q.{gif,mp4}``

The phase is gated on ``pvbatch`` being on PATH (or available at the
canonical ParaView 6.1 install location). It is opt-in (skipped unless
``--with-paraview``).
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from . import _common
from . import pipeline as _pipeline


_PVBATCH_CANDIDATES = [
    "/local/data/public/skcb2/tddft/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch",
    shutil.which("pvbatch") or "",
]


def _find_pvbatch() -> str | None:
    for p in _PVBATCH_CANDIDATES:
        if p and Path(p).is_file():
            return p
    return None


_TIME_RE = re.compile(r"_t(\d{6})\.vti$")


def _step_from_filename(p: Path) -> int:
    m = _TIME_RE.search(p.name)
    return int(m.group(1)) if m else -1


def _pvbatch_script() -> str:
    """The pvbatch driver. Reads its config from sys.argv[1] (a JSON file)."""
    return r"""
import json, sys
from paraview.simple import *

cfg = json.loads(open(sys.argv[1]).read())

# Helper: build a Volume rep for one VTI series.
def make_volume(files, array_name, color_preset, opacity_points, scalar_range):
    src = XMLImageDataReader(FileName=files)
    src.PointArrayStatus = [array_name]
    rep = Show(src, GetActiveViewOrCreate('RenderView'))
    rep.Representation = 'Volume'
    # Note: ParaView 6.1's UniformGridRepresentation no longer exposes
    # SelectMapper as a Python attribute; the default Smart mapper is used.
    ColorBy(rep, ('POINTS', array_name))
    lut = GetColorTransferFunction(array_name)
    lut.RescaleTransferFunction(*scalar_range)
    try:
        lut.UseLogScale = 1
    except Exception:
        pass
    try:
        lut.ApplyPreset(color_preset, True)
    except Exception:
        pass
    pwf = GetOpacityTransferFunction(array_name)
    # opacity_points is a list of [value, opacity] pairs
    flat = []
    for v, a in opacity_points:
        flat.extend([float(v), float(a), 0.5, 0.0])
    pwf.Points = flat
    return src, rep

view = GetActiveViewOrCreate('RenderView')
view.Background = [0.0, 0.0, 0.0]
view.OrientationAxesVisibility = 1
view.UseColorPaletteForBackground = 0
view.ViewSize = cfg['image_size']

system_src, system_rep = make_volume(
    cfg['system_files'], cfg['array_name'],
    cfg['system_color_preset'],
    cfg['system_opacity_points'],
    cfg['system_scalar_range'])
wp_src, wp_rep = make_volume(
    cfg['wp_files'], cfg['array_name'],
    cfg['wp_color_preset'],
    cfg['wp_opacity_points'],
    cfg['wp_scalar_range'])

# Time-aware animation across the union of times from both series.
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
n_frames = scene.NumberOfFrames

# Frame the bounds.
ResetCamera(view)

# Camera position: explicit so the two views are reproducible regardless
# of ParaView's azimuth/elevation accumulation order.
import math
def set_camera(view, azimuth_deg, elevation_deg, distance_factor):
    cam = view.GetActiveCamera()
    # Start from the head-on view: camera at -Z far, looking toward +Z.
    # ParaView's default ResetCamera looks down -Z from +Z; we invert it.
    fp = list(cam.GetFocalPoint())
    pos = list(cam.GetPosition())
    # Distance from focal point to current position
    dx = pos[0]-fp[0]; dy = pos[1]-fp[1]; dz = pos[2]-fp[2]
    d = math.sqrt(dx*dx + dy*dy + dz*dz) * distance_factor
    # Place camera at -Z relative to focal point.
    cam.SetPosition(fp[0], fp[1], fp[2] - d)
    cam.SetFocalPoint(fp[0], fp[1], fp[2])
    cam.SetViewUp(0.0, 1.0, 0.0)
    # Now apply yaw + elevation about the up vector.
    cam.Azimuth(azimuth_deg)
    cam.Elevation(elevation_deg)
    view.Update()

set_camera(view, cfg['camera_azimuth_deg'], cfg['camera_elevation_deg'],
           cfg['distance_factor'])

# Render each timestep.
times = scene.TimeKeeper.TimestepValues
out_dir = cfg['frames_dir']
prefix = cfg['filename_prefix']
import os; os.makedirs(out_dir, exist_ok=True)
n_emitted = 0
for i, t in enumerate(times):
    scene.AnimationTime = float(t)
    scene.UpdateAnimationUsingDataTimeSteps()
    Render(view)
    fname = os.path.join(out_dir, prefix + ('_%04d.png' % i))
    SaveScreenshot(fname, view, ImageResolution=cfg['image_size'])
    n_emitted += 1

print('rendered_frames=%d' % n_emitted)
"""


def _opacity_points(p50: float, p99: float) -> list[list[float]]:
    """Three-point piecewise linear opacity ramp (from the plan §2.7)."""
    return [
        [0.0, 0.0],
        [float(p50), 0.05],
        [float(p99), 0.6],
    ]


def _scalar_range(files: list[Path]) -> tuple[float, float, float]:
    """Return (vmin, p50, p99) for a series of VTIs."""
    import numpy as np
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    # Sample the centre frame to estimate the range cheaply.
    centre = files[len(files) // 2]
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(centre))
    reader.Update()
    arr = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(0))
    vmin = float(arr.min()); p50 = float(np.percentile(arr, 50))
    p99 = float(np.percentile(arr, 99))
    if p99 <= vmin:
        p99 = vmin + 1.0
    return vmin, p50, p99


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    pvbatch = _find_pvbatch()
    if pvbatch is None:
        _pipeline.skip("pvbatch not found; install ParaView 6.1 or set PATH")

    raw_vti = results_dir / "raw" / "vti"
    sys_dir = raw_vti / "density_rt_system"
    wp_dir = raw_vti / "density_rt_wp"
    sys_files = _common.list_vti_series(sys_dir, "density_rt_system")
    wp_files = _common.list_vti_series(wp_dir, "density_rt_wp")
    if not sys_files or not wp_files:
        _pipeline.skip(
            f"missing VTI series: system={len(sys_files)}, wp={len(wp_files)}")

    out_dir = _common.ensure_dir(
        results_dir / "analysis" / "density" / "paraview_3d")

    # Estimate per-series scalar range + opacity control points.
    sys_min, sys_p50, sys_p99 = _scalar_range(sys_files)
    wp_min, wp_p50, wp_p99 = _scalar_range(wp_files)

    # The two camera setups — see plan §2.7.
    views = [
        ("view_headon", dict(camera_azimuth_deg=0.0,
                             camera_elevation_deg=0.0,
                             distance_factor=1.0)),
        ("view_3q",     dict(camera_azimuth_deg=45.0,
                             camera_elevation_deg=25.0,
                             distance_factor=1.15)),
    ]

    notes: dict = {"out_dir": str(out_dir), "pvbatch": pvbatch}

    for view_name, cam in views:
        out_stem = out_dir / f"volume_overlay_{view_name}"
        gif_path = out_stem.with_suffix(".gif")
        if not _common.need_rebuild(gif_path, rebuild):
            notes[view_name] = str(gif_path) + " (cached)"
            continue
        with tempfile.TemporaryDirectory(prefix="inq_pv3d_") as tmp:
            tmp_path = Path(tmp)
            frames_dir = tmp_path / "frames"
            cfg = {
                "system_files": [str(p) for p in sys_files],
                "wp_files":     [str(p) for p in wp_files],
                "array_name":   "density",
                "system_color_preset": "Blue to Red Rainbow",  # Blues-style
                "wp_color_preset":     "Orange",
                "system_opacity_points": _opacity_points(sys_p50, sys_p99),
                "wp_opacity_points":     _opacity_points(wp_p50, wp_p99),
                "system_scalar_range":   [max(sys_min, 1e-12), sys_p99],
                "wp_scalar_range":       [max(wp_min, 1e-12),  wp_p99],
                "image_size":            [800, 800],
                "frames_dir":            str(frames_dir),
                "filename_prefix":       view_name,
                **cam,
            }
            cfg_path = tmp_path / "config.json"
            cfg_path.write_text(json.dumps(cfg, indent=2))
            script_path = tmp_path / "render.py"
            script_path.write_text(_pvbatch_script())

            cmd = [pvbatch, "--force-offscreen-rendering",
                   str(script_path), str(cfg_path)]
            try:
                subprocess.run(cmd, check=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               timeout=3600)
            except subprocess.CalledProcessError as e:
                err = (e.stderr or b"").decode("utf-8", errors="replace")[:2000]
                _pipeline.skip(f"{view_name} pvbatch failed: {err}")

            png_paths = sorted(frames_dir.glob(f"{view_name}_*.png"))
            if not png_paths:
                _pipeline.skip(f"{view_name}: pvbatch produced 0 frames")
            outs = _common.write_animation(out_stem, png_paths, fps=8)
            notes[view_name] = {
                "gif": str(outs["gif"]),
                "mp4": str(outs["mp4"]) if outs["mp4"] else None,
                "n_frames": len(png_paths),
            }

    return notes
