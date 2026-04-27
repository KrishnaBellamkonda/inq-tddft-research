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
    """The pvbatch driver. Reads its config from sys.argv[1] (a JSON file).

    Each series gets a *unique* registered array name (``density_system``,
    ``density_wp``) so ParaView builds a separate LUT/PWF per series even
    though both VTI files contain a scalar field called ``density``. This
    avoids the all-blue failure mode where both volumes shared the system's
    LUT and the WP became invisible.
    """
    return r"""
import json, sys, os
from paraview.simple import *

cfg = json.loads(open(sys.argv[1]).read())


def _make_lut_pwf(array_name, scalar_range, rgb_points, opacity_points):
    '''Build a colour LUT and opacity PWF with explicit RGB control points.

    rgb_points : list of [normalised_position_in_0_1, R, G, B] (RGB in [0,1])
    opacity_points : list of [density_value, opacity, midpoint, sharpness]
    '''
    lo, hi = float(scalar_range[0]), float(scalar_range[1])
    if hi <= lo:
        hi = lo + 1.0
    span = hi - lo
    # Linear LUT (no log) with explicit RGB points placed at
    # absolute density values lo + frac * span.
    flat = []
    for frac, r, g, b in rgb_points:
        flat.extend([lo + float(frac) * span,
                     float(r), float(g), float(b)])
    lut = GetColorTransferFunction(array_name)
    lut.RGBPoints = flat
    lut.ColorSpace = 'RGB'
    try:
        lut.UseLogScale = 0
    except Exception:
        pass
    pwf = GetOpacityTransferFunction(array_name)
    pflat = []
    for v, a in opacity_points:
        pflat.extend([float(v), float(a), 0.5, 0.0])
    pwf.Points = pflat
    return lut, pwf


def make_volume(files, registered_name, source_array,
                rgb_points, opacity_points, scalar_range):
    src = XMLImageDataReader(FileName=files)
    src.PointArrayStatus = [source_array]
    # Rename the active scalar so each volume has its own LUT.
    calc = Calculator(Input=src)
    calc.ResultArrayName = registered_name
    calc.Function = source_array
    rep = Show(calc, GetActiveViewOrCreate('RenderView'))
    rep.Representation = 'Volume'
    ColorBy(rep, ('POINTS', registered_name))
    _make_lut_pwf(registered_name, scalar_range, rgb_points, opacity_points)
    return src, calc, rep


view = GetActiveViewOrCreate('RenderView')
view.Background = [1.0, 1.0, 1.0]    # white background — black hides the
                                     # axes grid badly
view.OrientationAxesVisibility = 1
view.UseColorPaletteForBackground = 0
view.ViewSize = cfg['image_size']

# ---- Axes grid (data axes with bohr ticks) ----------------------------
ag = view.AxesGrid
ag.Visibility = 1
ag.AxesToLabel = 63   # all six faces
ag.XTitle = 'x (bohr)'
ag.YTitle = 'y (bohr)'
ag.ZTitle = 'z (bohr)'
try:
    ag.GridColor = [0.4, 0.4, 0.4]
    ag.XLabelColor = [0.0, 0.0, 0.0]
    ag.YLabelColor = [0.0, 0.0, 0.0]
    ag.ZLabelColor = [0.0, 0.0, 0.0]
    ag.XTitleColor = [0.0, 0.0, 0.0]
    ag.YTitleColor = [0.0, 0.0, 0.0]
    ag.ZTitleColor = [0.0, 0.0, 0.0]
except Exception:
    pass

# Build the two volumes.
system_src, system_calc, system_rep = make_volume(
    cfg['system_files'], 'density_system',
    cfg['array_name'],
    cfg['system_rgb_points'], cfg['system_opacity_points'],
    cfg['system_scalar_range'])
wp_src, wp_calc, wp_rep = make_volume(
    cfg['wp_files'], 'density_wp',
    cfg['array_name'],
    cfg['wp_rgb_points'], cfg['wp_opacity_points'],
    cfg['wp_scalar_range'])

# Time-aware animation across the union of times from both series.
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

# Frame the bounds before camera placement.
ResetCamera(view)

# Camera placement: explicit so the views are reproducible.
import math
def set_camera(view, azimuth_deg, elevation_deg, distance_factor):
    cam = view.GetActiveCamera()
    fp = list(cam.GetFocalPoint())
    pos = list(cam.GetPosition())
    dx = pos[0]-fp[0]; dy = pos[1]-fp[1]; dz = pos[2]-fp[2]
    d = math.sqrt(dx*dx + dy*dy + dz*dz) * distance_factor
    # Place camera at -Z relative to focal point so we look toward +Z.
    cam.SetPosition(fp[0], fp[1], fp[2] - d)
    cam.SetFocalPoint(fp[0], fp[1], fp[2])
    cam.SetViewUp(0.0, 1.0, 0.0)
    cam.Azimuth(azimuth_deg)
    cam.Elevation(elevation_deg)
    view.Update()

set_camera(view, cfg['camera_azimuth_deg'], cfg['camera_elevation_deg'],
           cfg['distance_factor'])

# Render every timestep.
times = scene.TimeKeeper.TimestepValues
out_dir = cfg['frames_dir']
prefix = cfg['filename_prefix']
os.makedirs(out_dir, exist_ok=True)
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


# Explicit RGB control points for the two volume series. Each list is
# [normalised_position_in_0_1, R, G, B] with R/G/B in [0, 1]. The pvbatch
# script places these at absolute density values lo + frac * (hi - lo).
#
# System ramp: white at the low end (faint molecular tail) → mid-blue at
# the median → deep blue at the peak. Picks contrast against a white
# background.
_SYSTEM_RGB_POINTS = [
    [0.00, 0.95, 0.95, 1.00],
    [0.50, 0.40, 0.55, 0.85],
    [1.00, 0.05, 0.15, 0.55],
]
# WP ramp: pale yellow → orange → deep red-orange. Visually distinct from
# the blue system without overlap.
_WP_RGB_POINTS = [
    [0.00, 1.00, 0.95, 0.75],
    [0.50, 1.00, 0.55, 0.10],
    [1.00, 0.65, 0.10, 0.00],
]


def _scalar_range(files: list[Path]) -> tuple[float, float, float, float]:
    """Return (p10, p50, p90, p99) across a 3-frame sample of the series.

    Sampling the start, middle, and end frames captures the dynamic range
    of the WP density (which is highly localised early but spreads later).
    The lower-bound for the LUT is p10 (not vmin) because vmin ≈ 0 plus
    log-or-linear-near-zero density makes the bulk of the volume render
    at the same colour bin.
    """
    import numpy as np
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    n = len(files)
    sample_idx = sorted({0, n // 2, n - 1})
    samples: list[np.ndarray] = []
    for i in sample_idx:
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(str(files[i]))
        reader.Update()
        samples.append(vtk_to_numpy(
            reader.GetOutput().GetPointData().GetArray(0)))
    arr = np.concatenate(samples)
    p10 = float(np.percentile(arr, 10))
    p50 = float(np.percentile(arr, 50))
    p90 = float(np.percentile(arr, 90))
    p99 = float(np.percentile(arr, 99))
    if p99 <= p10:
        p99 = p10 + 1.0
    return p10, p50, p90, p99


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

    # Estimate per-series scalar range. Use p10 (not min) as the lower
    # bound so the bulk near-zero volume doesn't dominate the LUT.
    sys_p10, sys_p50, sys_p90, sys_p99 = _scalar_range(sys_files)
    wp_p10,  wp_p50,  wp_p90,  wp_p99  = _scalar_range(wp_files)

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
                "system_rgb_points":     _SYSTEM_RGB_POINTS,
                "wp_rgb_points":         _WP_RGB_POINTS,
                "system_opacity_points": _opacity_points(sys_p50, sys_p99),
                "wp_opacity_points":     _opacity_points(wp_p50, wp_p99),
                # Linear (not log) range, clamped to [p10, p99] per series so
                # the bulk near-zero density doesn't crush the LUT.
                "system_scalar_range":   [sys_p10, sys_p99],
                "wp_scalar_range":       [wp_p10,  wp_p99],
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
