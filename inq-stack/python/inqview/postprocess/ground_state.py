"""Phase: ``gs`` — ground-state plots from ``raw/ground_state/``.

Produces under ``results/analysis/ground_state/``:

* ``density_gs_system_xy.png`` — midplane slice of the GS system density.
* ``gs_orbital_gallery.png``    — small grid of the per-orbital densities,
  if VTI orbitals are available.

If the only GS data is the VTIs in ``raw/vti/density_gs_system/`` (and not
the older ``.raw + .meta.txt`` pair), this phase is currently a soft-skip
because ``inqview.fields.RealField3D`` is loaded from raw+meta. Slice
plotting from VTI directly is a future enhancement; for now we record the
skip reason and move on.
"""

from __future__ import annotations

from pathlib import Path

from . import _common
from . import pipeline as _pipeline


def _load_vti_cube(path: Path):
    """Read a binary VTI into (cube[nx,ny,nz], origin, spacing)."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    import numpy as np
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64,
                                                                copy=False)
    cube = flat.reshape((nz, ny, nx)).transpose(2, 1, 0)
    return cube, tuple(img.GetOrigin()), tuple(img.GetSpacing())


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    out_dir = _common.ensure_dir(results_dir / "analysis" / "ground_state")
    raw_gs = results_dir / "raw" / "ground_state"
    raw_vti_gs = results_dir / "raw" / "vti" / "density_gs_system"
    raw_vti_orb = results_dir / "raw" / "vti" / "density_gs_orbitals"

    has_raw_density = (raw_gs / "density_system").exists() and any(
        (raw_gs / "density_system").glob("*.meta.txt")
    )
    has_vti = raw_vti_gs.exists() and any(raw_vti_gs.glob("*.vti"))
    has_vti_orbitals = raw_vti_orb.exists() and any(raw_vti_orb.glob("*.vti"))

    if not (has_raw_density or has_vti or has_vti_orbitals):
        _pipeline.skip(
            "no GS density data found "
            f"(neither {raw_gs}/density_system, {raw_vti_gs}, nor {raw_vti_orb})"
        )

    notes = {"out_dir": str(out_dir), "has_raw_density": has_raw_density,
             "has_vti": has_vti, "has_vti_orbitals": has_vti_orbitals}

    import matplotlib.pyplot as plt

    # ---- GS system density xy slice ------------------------------------
    out_xy = out_dir / "density_gs_system_xy.png"
    if has_raw_density:
        from .. import load_real_field, plot_density_slice
        meta = next((raw_gs / "density_system").glob("*.meta.txt"))
        field = load_real_field(meta)
        if _common.need_rebuild(out_xy, rebuild):
            plot_density_slice(field, out_xy, axis=2)
        notes["density_xy"] = str(out_xy)
    elif has_vti:
        if _common.need_rebuild(out_xy, rebuild):
            vti = next(raw_vti_gs.glob("*.vti"))
            cube, origin, spacing = _load_vti_cube(vti)
            slc = cube[:, :, cube.shape[2] // 2]
            ox, oy, oz = origin
            dx, dy, dz = spacing
            extent = [ox, ox + cube.shape[0] * dx,
                      oy, oy + cube.shape[1] * dy]
            fig, ax = plt.subplots(figsize=(5, 5), dpi=120)
            im = ax.imshow(slc.T, origin="lower", extent=extent,
                           cmap="viridis", aspect="equal")
            plt.colorbar(im, ax=ax, label="density (bohr^-3)")
            ax.set_xlabel("x (bohr)"); ax.set_ylabel("y (bohr)")
            ax.set_title(_common.title(run_name, "GS system density (xy mid-plane)"))
            fig.tight_layout()
            fig.savefig(out_xy); plt.close(fig)
        notes["density_xy"] = str(out_xy)

    # ---- GS orbital gallery (raw or VTI) -------------------------------
    out_gallery = out_dir / "gs_orbital_gallery.png"
    raw_orb = raw_gs / "density_gs_orbitals"

    meta_files: list[Path] = []
    if raw_orb.exists():
        meta_files = sorted(raw_orb.glob("orbital_*.meta.txt"))

    if meta_files:
        from .. import load_real_field
        n = len(meta_files)
        cols = 6; rows = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols,
                                 figsize=(2 * cols, 2 * rows), dpi=120)
        axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
        for i, meta in enumerate(meta_files):
            f = load_real_field(meta)
            arr = f.array
            slc = arr[:, :, arr.shape[2] // 2]
            ax = axes[i]
            ax.imshow(slc, origin="lower", aspect="equal", cmap="viridis")
            ax.set_title(_common.title(run_name, f"orb {i:02d}"),
                         fontsize="x-small")
            ax.set_xticks([]); ax.set_yticks([])
        for j in range(n, rows * cols):
            axes[j].axis("off")
        fig.tight_layout()
        fig.savefig(out_gallery); plt.close(fig)
        notes["orbital_gallery"] = str(out_gallery)
    elif has_vti_orbitals and _common.need_rebuild(out_gallery, rebuild):
        vti_files = sorted(raw_vti_orb.glob("orbital_*.vti"))
        n = len(vti_files)
        cols = 6; rows = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols,
                                 figsize=(2 * cols, 2 * rows), dpi=120)
        axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
        for i, vti in enumerate(vti_files):
            cube, _, _ = _load_vti_cube(vti)
            slc = cube[:, :, cube.shape[2] // 2]
            ax = axes[i]
            ax.imshow(slc.T, origin="lower", aspect="equal", cmap="viridis")
            ax.set_title(_common.title(run_name, f"orb {i:02d}"),
                         fontsize="x-small")
            ax.set_xticks([]); ax.set_yticks([])
        for j in range(n, rows * cols):
            axes[j].axis("off")
        fig.tight_layout()
        fig.savefig(out_gallery); plt.close(fig)
        notes["orbital_gallery"] = str(out_gallery)

    return notes
