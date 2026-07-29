"""
inqview.pipeline.orbitals_per_kpoint
========================================

Reads the band-major per-(band, k) VTI dump produced by
`inqkit_smoke::dump_orbitals_per_kpoint` (helpers at
`Tutorial/_inqkit_tests/_orbital_dump_helpers.hpp`) and produces:

  band_NNN/re_psi_grid.png       — mid-cell-slice grid of Re psi across all k
  band_NNN/im_psi_grid.png       — same for Im psi
  band_NNN/density_grid.png      — same for |psi|^2
  bands_summary.png              — evalue vs k-point index for each chosen band
  paraview_recipe.md             — short note on driving the headline 3D view

Wired into inqview.pipeline.pipeline as phase `orbitals_per_kpoint`,
gated on the existence of an `orbitals_per_kpoint/orbital_index.csv` file
under `<run_dir>/results/analysis/ground_state/`.

Public API:
    run(run_dir)            — run the whole phase on one run dir.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

# TODO: Can a band structure plot be made using this script?


# ----- VTI reader (vendored to avoid a hard inqview-wide dependency) -------

def _read_vti(path: Path) -> tuple[np.ndarray, tuple[float, float, float], tuple[float, float, float]]:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    dx, dy, dz = img.GetSpacing()
    ox, oy, oz = img.GetOrigin()
    arr = vtk_to_numpy(img.GetPointData().GetArray(0))
    field = arr.reshape((nz, ny, nx)).transpose(2, 1, 0)  # (nx,ny,nz)
    return field, (dx, dy, dz), (ox, oy, oz)


def _slice_extent(origin: tuple, spacing: tuple, shape: tuple,
                  axis_a: int, axis_b: int) -> list[float]:
    return [origin[axis_a],
            origin[axis_a] + shape[axis_a] * spacing[axis_a],
            origin[axis_b],
            origin[axis_b] + shape[axis_b] * spacing[axis_b]]


# ----- Public phase entry point ---------------------------------------------

def run(run_dir: str | Path) -> int:
    """Run the orbitals-per-kpoint visualisation phase on a run dir.

    Returns 0 on success, non-zero on failure (no figure produced).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    run_dir = Path(run_dir).resolve()
    root = run_dir / "results" / "analysis" / "ground_state" / "orbitals_per_kpoint"
    csv = root / "orbital_index.csv"
    if not csv.exists():
        print(f"  [orbitals_per_kpoint] skipped: {csv} not found")
        return 1

    df = pd.read_csv(csv)
    print(f"  [orbitals_per_kpoint] {len(df)} rows from {csv}")

    bands = sorted(df["band"].unique().tolist())
    n_kpts = df["kpoint_index"].nunique()
    print(f"  [orbitals_per_kpoint] {len(bands)} bands × {n_kpts} k-points")

    # ----- per-band figures: re, im, density grids ------------------------
    for band in bands:
        band_df = df[df["band"] == band].sort_values("kpoint_index")
        band_dir = root / f"band_{band:03d}"
        if not band_dir.is_dir():
            print(f"    band {band}: directory missing, skipping")
            continue

        # Read all VTIs first, on a common scale per component.
        re_fields, im_fields, den_fields = [], [], []
        meta = None
        for ki in band_df["kpoint_index"].astype(int).tolist():
            re_p  = band_dir / f"re_psi_t{ki:06d}.vti"
            im_p  = band_dir / f"im_psi_t{ki:06d}.vti"
            den_p = band_dir / f"density_t{ki:06d}.vti"
            if not (re_p.exists() and im_p.exists() and den_p.exists()):
                print(f"    band {band} k {ki}: missing VTI, skipping band")
                break
            f_re,  spc, org = _read_vti(re_p)
            f_im,  _,   _   = _read_vti(im_p)
            f_den, _,   _   = _read_vti(den_p)
            re_fields.append(f_re)
            im_fields.append(f_im)
            den_fields.append(f_den)
            if meta is None:
                meta = (spc, org, f_re.shape)
        else:
            _make_band_grid(band, band_df, re_fields,  meta, "re_psi",  band_dir / "re_psi_grid.png",  cmap="RdBu_r", symmetric=True)
            _make_band_grid(band, band_df, im_fields,  meta, "im_psi",  band_dir / "im_psi_grid.png",  cmap="RdBu_r", symmetric=True)
            _make_band_grid(band, band_df, den_fields, meta, "density", band_dir / "density_grid.png", cmap="viridis", symmetric=False)
            continue
        # Loop ended via `break`.
        continue

    # ----- bands_summary.png ----------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for band in bands:
        sub = df[df["band"] == band].sort_values("kpoint_index")
        ax.plot(sub["kpoint_index"], sub["evalue_ev"], "o-", lw=1.0, label=f"band {band}")
    ax.set_xlabel("k-point index")
    ax.set_ylabel("eigenvalue (eV)")
    ax.set_title("Eigenvalues vs k-point index — selected bands")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    summary_png = root / "bands_summary.png"
    plt.savefig(summary_png, dpi=150)
    plt.close()
    print(f"    wrote {summary_png}")

    # ----- paraview_recipe.md --------------------------------------------
    recipe = root / "paraview_recipe.md"
    recipe.write_text(_PARAVIEW_RECIPE)
    print(f"    wrote {recipe}")

    return 0


_PARAVIEW_RECIPE = """\
# ParaView recipe for per-(band, k) orbitals

Each `band_NNN/` directory holds three VTI series, one per component:

* `re_psi_t<KKK>.vti`   — Re psi_{n,k}(r)
* `im_psi_t<KKK>.vti`   — Im psi_{n,k}(r)
* `density_t<KKK>.vti`  — |psi_{n,k}(r)|^2

The "step" axis (KKK) is the kpoint index, not real time. Animating in
ParaView sweeps the orbital across the BZ at fixed band — this is the
headline pedagogical view.

## Headline view: Bloch phase across the BZ for a single band

1. File → Open → `band_001/re_psi_t..vti` (use the file series; ParaView
   will collapse the names automatically).
2. Apply.
3. Volume render with a diverging colourmap (e.g. *cool to warm*).
4. Use the time-step slider to step through the k-points. Re psi should
   show clear sinusoidal modulation at the lattice scale for k-points
   far from Gamma; the modulation period scales as 1/|k|.
5. Add a second View, open `im_psi_t..vti`, link the time controls.
   Im psi should be shifted by pi/2 of the same wavelength.
6. Add a third View with `density_t..vti`. |psi|^2 will be much more
   uniform across k — the Bloch phase is the carrier of k.

## Side-by-side comparison: same k-point across multiple bands

1. Open `band_001/re_psi_t000000.vti`, `band_012/re_psi_t000000.vti`, etc.
2. Use linked views with the same colour table.

The 2D grid PNGs (`re_psi_grid.png`, `im_psi_grid.png`,
`density_grid.png`) in this directory are pre-rendered summaries useful
for reports — open them with any image viewer.
"""


def _make_band_grid(band: int,
                    band_df: pd.DataFrame,
                    fields: list[np.ndarray],
                    meta,
                    label: str,
                    out_png: Path,
                    cmap: str,
                    symmetric: bool):
    """One PNG per (band, component): N_k panels of mid-cell slices."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    spc, org, shape = meta
    nk = len(fields)
    if nk == 0:
        return

    # Fixed colour limits across panels.
    if symmetric:
        vmax = max(float(np.abs(f).max()) for f in fields)
        vmin = -vmax
    else:
        vmax = max(float(f.max()) for f in fields)
        vmin = 0.0

    n_cols = min(nk, 4)
    n_rows = int(np.ceil(nk / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.2 * n_cols, 3.2 * n_rows),
                             squeeze=False)
    for ax in axes.flat:
        ax.axis("off")

    # Slice along z through the plane of maximum sum_i |f[:, :, i]|, summed
    # across all k-fields. This finds where the orbital actually has support
    # — robust whether the atom is at cell-centre or cell-edge.
    z_support = sum(np.abs(f).sum(axis=(0, 1)) for f in fields)
    nz_mid = int(np.argmax(z_support))
    extent = _slice_extent(org, spc, shape, axis_a=0, axis_b=1)

    for i, f in enumerate(fields):
        r, c = divmod(i, n_cols)
        ax = axes[r, c]
        slc = f[:, :, nz_mid].T
        im = ax.imshow(slc, origin="lower", extent=extent,
                       cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
        ki = int(band_df.iloc[i]["kpoint_index"])
        kx = band_df.iloc[i]["kx"]
        ky = band_df.iloc[i]["ky"]
        kz = band_df.iloc[i]["kz"]
        ax.set_title(f"k {ki}: ({kx:.2f}, {ky:.2f}, {kz:.2f})", fontsize=8)
        ax.axis("on")
        ax.set_xticks([])
        ax.set_yticks([])

    z_at = org[2] + nz_mid * spc[2]
    fig.suptitle(
        f"band {band}: {label}  "
        f"(xy slice at z={z_at:.2f} Bohr — chosen to maximise |field|, "
        f"fixed colour scale)",
        fontsize=11)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7,
                 label=label)
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    wrote {out_png}")


# ----- pipeline.py wiring ---------------------------------------------------

PHASE_NAME = "orbitals_per_kpoint"


def is_applicable(run_dir: str | Path) -> bool:
    return (Path(run_dir) / "results" / "analysis" / "ground_state"
            / "orbitals_per_kpoint" / "orbital_index.csv").exists()


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("run_dir", type=Path)
    args = p.parse_args()
    raise SystemExit(run(args.run_dir))


if __name__ == "__main__":
    main()
