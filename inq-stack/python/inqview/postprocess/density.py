"""Phase: ``density`` — 2D slice GIFs of total / system / WP densities.

For each density category (``density_rt_total``, ``density_rt_system``,
``density_rt_wp``) and each slice plane (xy, xz, yz), produces one GIF
under ``results/analysis/density/`` with a fixed colour scale across frames
per ``docs/visualisation-instructions-v1.md``.

The frames are loaded from ``results/raw/vti/<category>/`` (the C++ runs
write VTI directly). VTI loading uses the same vtkXMLImageDataReader that
``inqview.vti`` already wraps — see :func:`_load_vti_array`.

If a category is missing, that GIF is skipped silently (logged to the
phase result).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from . import _common
from . import pipeline as _pipeline

_PLANES = {
    "xy": (2, "x (bohr)", "y (bohr)"),
    "xz": (1, "x (bohr)", "z (bohr)"),
    "yz": (0, "y (bohr)", "z (bohr)"),
}

_CATEGORIES = ("density_rt_total", "density_rt_system", "density_rt_wp")


def _load_vti_array(path: Path) -> tuple[np.ndarray, dict]:
    """Return (3D float array shape (nx,ny,nz), meta dict).

    Reads via vtkXMLImageDataReader. Falls back gracefully if VTK isn't
    importable; the phase will be skipped in that case by the caller.
    """
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(f"VTK is required to read {path}") from exc

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    pd = img.GetPointData()
    arr_vtk = pd.GetArray(0)
    flat = vtk_to_numpy(arr_vtk).astype(np.float64, copy=False)
    # VTK ImageData uses x-fastest; reshape to (nz,ny,nx) then transpose.
    cube = flat.reshape((nz, ny, nx)).transpose(2, 1, 0)
    meta = {
        "nx": nx, "ny": ny, "nz": nz,
        "origin": tuple(img.GetOrigin()),
        "spacing": tuple(img.GetSpacing()),
    }
    return cube, meta


_TIME_RE = re.compile(r"_t(\d{6})\.vti$")


def _step_from_filename(p: Path) -> int:
    m = _TIME_RE.search(p.name)
    return int(m.group(1)) if m else -1


def _global_vmin_vmax(arrs: list[np.ndarray], percentile: float | None) -> tuple[float, float]:
    if percentile is None:
        vmin = min(float(a.min()) for a in arrs)
        vmax = max(float(a.max()) for a in arrs)
    else:
        # Use percentile clipping for visibility on peaked WP density.
        flat = np.concatenate([a.ravel() for a in arrs])
        lo = float(np.percentile(flat, 100 - percentile))
        hi = float(np.percentile(flat, percentile))
        vmin, vmax = lo, hi
    if vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def run(results_dir: Path, *, run_name: str, rebuild: bool,
        percentile: float | None = 99.0, dt_au: float = 0.020,
        write_every: int = 10, **_) -> dict:
    out_dir = _common.ensure_dir(results_dir / "analysis" / "density")
    raw_vti = results_dir / "raw" / "vti"

    notes: dict = {"out_dir": str(out_dir)}

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        _pipeline.skip(f"missing matplotlib: {exc}")

    for cat in _CATEGORIES:
        cat_dir = raw_vti / cat
        files = _common.list_vti_series(cat_dir, cat)
        if not files:
            notes[cat] = "no VTI files"
            continue

        # Load all frames once (memory: ~Nx*Ny*Nz*8 bytes per frame).
        frames: list[tuple[int, float, np.ndarray, dict]] = []
        for f in files:
            step = _step_from_filename(f)
            time_au = step * dt_au if step > 0 else 0.0
            cube, meta = _load_vti_array(f)
            frames.append((step, time_au, cube, meta))

        # Per-plane: emit linear + log animation pairs (each pair as gif+mp4).
        for plane, (axis, xlab, ylab) in _PLANES.items():
            slices = [f[2].take(f[2].shape[axis] // 2, axis=axis)
                      for f in frames]
            vmin_lin, vmax_lin = _global_vmin_vmax(slices, percentile)
            # log1p domain: floor at 0; same vmax derived from log1p(percentile).
            vmin_log = 0.0
            vmax_log = float(np.log1p(vmax_lin))

            for scale_label, scale_kwargs, vmin, vmax, transform, cbar_label in [
                ("",     dict(), vmin_lin, vmax_lin, lambda a: a,
                 r"density (bohr$^{-3}$)"),
                ("_log", dict(), vmin_log, vmax_log, np.log1p,
                 r"log$_{10}$(1 + density)"),
            ]:
                out_stem = out_dir / f"{cat[len('density_rt_'):]}_{plane}{scale_label}"
                gif_path = out_stem.with_suffix(".gif")
                if not _common.need_rebuild(gif_path, rebuild):
                    notes[f"{cat}_{plane}{scale_label}"] = str(gif_path) + " (cached)"
                    continue

                tmp_dir = _common.ensure_dir(
                    out_dir / f".__tmp_{cat}_{plane}{scale_label}")
                png_paths: list[Path] = []
                for i, (step, t_au, cube, meta) in enumerate(frames):
                    slc = cube.take(cube.shape[axis] // 2, axis=axis)
                    slc_t = transform(slc)
                    fig, ax = plt.subplots(figsize=(5, 5), dpi=120)
                    ox, oy, oz = meta["origin"]
                    dx, dy, dz = meta["spacing"]
                    if plane == "xy":
                        extent = [ox, ox + meta["nx"] * dx,
                                  oy, oy + meta["ny"] * dy]
                    elif plane == "xz":
                        extent = [ox, ox + meta["nx"] * dx,
                                  oz, oz + meta["nz"] * dz]
                    else:  # yz
                        extent = [oy, oy + meta["ny"] * dy,
                                  oz, oz + meta["nz"] * dz]
                    im = ax.imshow(slc_t.T, origin="lower", extent=extent,
                                   cmap="viridis", aspect="equal",
                                   vmin=vmin, vmax=vmax)
                    plt.colorbar(im, ax=ax, label=cbar_label)
                    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
                    ax.set_title(_common.title(
                        run_name, f"{cat} {plane}{scale_label or ''}",
                        step=step, total_steps=frames[-1][0],
                        time_au=t_au,
                    ))
                    fig.tight_layout()
                    p = tmp_dir / f"frame_{i:04d}.png"
                    fig.savefig(p)
                    plt.close(fig)
                    png_paths.append(p)

                outs = _common.write_animation(out_stem, png_paths, fps=8)
                for p in png_paths:
                    p.unlink(missing_ok=True)
                tmp_dir.rmdir()
                notes[f"{cat}_{plane}{scale_label}"] = {
                    "gif": str(outs["gif"]),
                    "mp4": str(outs["mp4"]) if outs["mp4"] else None,
                }

    return notes
