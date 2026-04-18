"""plot_orbitals_3d.py — 3D isosurface rendering of jellium KS orbitals.

Requires: pyvista, vtk, numpy  (quantum-wave-packet pyenv)
Run with:
    pyenv activate quantum-wave-packet   # (NOT active when running inq-run)
    python plot_orbitals_3d.py

For each shell (|n|² = 0, 1, 2, 3, 4) this script:
  1. Reads the _real.txt orbital file header to extract k = (kx, ky, kz) and L.
  2. Reconstructs Re[ψ_k(r)] = cos(k·r)/√Ω analytically on a 3D grid (no large
     binary data needed from C++; plane-wave formula is exact).
  3. Renders two semi-transparent isosurfaces at ±isovalue, coloured red/blue.
  4. Saves one PNG per shell to results/orbital_3d_n2_<M>.png and a tiled
     summary figure to results/orbitals_3d_summary.png.

Notes:
  - |n|²=0 has k=(0,0,0) so Re[ψ]=const — no isosurface; shown as volume slice.
  - Runs fully headless (pv.OFF_SCREEN=True); no display required.
  - Grid resolution is controlled by N_RENDER_PTS (default 60³ — fast; raise to
    120 for publication-quality renders).
"""

import glob
import os
import re
import sys

import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True   # headless rendering — no display needed

# ── Parameters ───────────────────────────────────────────────────────────────

L          = 40.0          # cell side (bohr); must match run.cpp
N_RENDER   = 80            # grid points per axis for 3D render (80³ ≈ fast)
ISO_FRAC   = 0.40          # isosurface at ±ISO_FRAC × max(|Re[ψ]|)
OPACITY    = 0.65          # surface opacity (0=transparent, 1=opaque)
OUT_DIR    = "results"

# Colours
C_POS = (0.85, 0.20, 0.15)   # red — positive lobe
C_NEG = (0.15, 0.35, 0.80)   # blue — negative lobe
C_BG  = (0.97, 0.97, 0.97)   # light grey background

# ── Helpers ──────────────────────────────────────────────────────────────────

def parse_kvec_from_header(fname):
    """Extract (kx, ky, kz) and L from the _real.txt file header."""
    k = None
    with open(fname) as f:
        for line in f:
            if not line.startswith('#'):
                break
            # Line looks like: # k = (0.157080, 0.000000, 0.000000) bohr⁻¹
            m = re.search(r'k\s*=\s*\(([\d.eE+\-]+),\s*([\d.eE+\-]+),\s*([\d.eE+\-]+)\)', line)
            if m:
                k = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    return k


def build_3d_field(kx, ky, kz, L, N):
    """Return Re[ψ_k] on an N×N×N grid covering [0, L)³."""
    Omega = L ** 3
    norm  = 1.0 / np.sqrt(Omega)
    x = np.linspace(0, L, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')   # shape (N,N,N)
    kr = kx * X + ky * Y + kz * Z
    return norm * np.cos(kr)   # Re[ψ_k]


def render_orbital(oi, n2, kx, ky, kz, out_path):
    """Render one orbital's 3D isosurface and save to out_path (PNG)."""
    field = build_3d_field(kx, ky, kz, L, N_RENDER)
    fmax  = np.abs(field).max()

    # Build PyVista uniform grid
    grid = pv.ImageData()
    grid.dimensions = np.array(field.shape) + 1      # cell corners
    grid.spacing    = (L / N_RENDER,) * 3
    grid.origin     = (0.0, 0.0, 0.0)
    grid.cell_data['Re_psi'] = field.ravel(order='F')   # Fortran order for VTK

    pl = pv.Plotter(off_screen=True, window_size=(900, 800))
    pl.set_background(C_BG)
    pl.add_axes(line_width=3, color='black')

    # Bounding-box wire frame for spatial reference
    pl.add_mesh(pv.Box(bounds=(0, L, 0, L, 0, L)),
                style='wireframe', color='grey', opacity=0.3, line_width=1)

    if n2 == 0:
        # k=(0,0,0) → Re[ψ] is constant — show a mid-plane slice instead
        slice_mesh = grid.slice(normal='z')
        pl.add_mesh(slice_mesh, scalars='Re_psi', cmap='RdBu_r',
                    show_scalar_bar=True, scalar_bar_args={'title': 'Re[ψ]'})
        pl.add_text(f'|n|²=0  k=(0,0,0)  Re[ψ]=const  (z-slice shown)',
                    position='upper_left', font_size=11, color='black')
    else:
        iso_val = ISO_FRAC * fmax

        # Convert cell data → point data for smooth isosurface
        pgrid  = grid.cell_data_to_point_data()
        surf_p = pgrid.contour([+iso_val], scalars='Re_psi')
        surf_n = pgrid.contour([-iso_val], scalars='Re_psi')

        if surf_p.n_points > 0:
            pl.add_mesh(surf_p, color=C_POS, opacity=OPACITY, smooth_shading=True)
        if surf_n.n_points > 0:
            pl.add_mesh(surf_n, color=C_NEG, opacity=OPACITY, smooth_shading=True)

        # Semi-transparent mid-plane slice for context
        sl = grid.slice(normal='z', origin=(0, 0, L/2))
        pl.add_mesh(sl, scalars='Re_psi', cmap='RdBu_r', opacity=0.25,
                    show_scalar_bar=False)

        kstr = f'({kx/(2*np.pi/L):.0f},{ky/(2*np.pi/L):.0f},{kz/(2*np.pi/L):.0f})k₀'
        pl.add_text(
            f'|n|²={n2}   k={kstr}   iso=±{iso_val:.4f} bohr⁻³/²',
            position='upper_left', font_size=11, color='black'
        )

    # Isometric-style camera
    pl.camera_position = 'iso'
    pl.camera.azimuth  = 30
    pl.camera.elevation = 20
    pl.camera.zoom(0.85)

    pl.screenshot(out_path)
    pl.close()
    print(f'  Saved {out_path}')


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    re_files = sorted(glob.glob('results/orbitals/orbital_*_real.txt'))
    if not re_files:
        print('ERROR: No orbital _real.txt files found. Run inq-run in this directory first.')
        sys.exit(1)

    png_paths = []
    for rf in re_files:
        # Parse index and n2 from filename  orbital_N_n2_M_real.txt
        m = re.search(r'orbital_(\d+)_n2_(\d+)_real\.txt', rf)
        if not m:
            continue
        oi, n2 = int(m.group(1)), int(m.group(2))

        k = parse_kvec_from_header(rf)
        if k is None:
            print(f'  [skip] Could not parse k-vector from {rf}')
            continue
        kx, ky, kz = k

        out_png = f'{OUT_DIR}/orbital_3d_n2_{n2}.png'
        print(f'Rendering orbital {oi}  |n|²={n2}  k=({kx:.4f}, {ky:.4f}, {kz:.4f}) ...')
        render_orbital(oi, n2, kx, ky, kz, out_png)
        png_paths.append((n2, out_png))

    # Tile into a summary figure using matplotlib (no extra deps)
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    n = len(png_paths)
    if n == 0:
        print('No images rendered.')
        return

    fig, axes = plt.subplots(1, n, figsize=(4*n, 4.5))
    if n == 1:
        axes = [axes]

    for ax, (n2, path) in zip(axes, png_paths):
        img = mpimg.imread(path)
        ax.imshow(img)
        ax.set_title(f'$|n|^2={n2}$', fontsize=12)
        ax.axis('off')

    fig.suptitle(
        f'Jellium KS orbitals — Re[$\\psi_k$] isosurfaces\n'
        f'$L={L}\\,a_0$,  iso level = {ISO_FRAC:.0%}$\\times$max,  '
        f'grid {N_RENDER}³,  red=positive, blue=negative',
        fontsize=11, y=1.01
    )
    plt.tight_layout()
    summary_path = f'{OUT_DIR}/orbitals_3d_summary.png'
    fig.savefig(summary_path, bbox_inches='tight', dpi=120)
    plt.close(fig)
    print(f'\nSaved tiled summary: {summary_path}')


if __name__ == '__main__':
    main()
