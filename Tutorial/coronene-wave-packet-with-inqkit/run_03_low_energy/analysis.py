"""
run_03_low_energy: analysis.py
Validates inqkit WavePacket injection and visualises GS + WP orbital densities.

Run with:
  /local/data/public/skcb2/pyenv/versions/3.10.19/envs/quantum-wave-packet/bin/python3 analysis.py

Expected:
  N_elec (GS density) ≈ 108.0  (coronene 108 valence electrons)
  wp_norm             ≈  1.0   (normalised Gaussian after orthogonalisation)
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.data import SimulationData, load_real_field
from inqview.vti import convert_real_meta_to_vti
from inqview.paraview import ParaViewPipeline, VolumeRenderSpec, AnimationSpec
from inqview.fields import RealField3D

RESULTS   = Path(__file__).parent / "results"
VISDIR    = RESULTS / "visualisation"
VTIDIR    = VISDIR / "vti"
FRAMESDIR = VISDIR / "frames"
SLICEDIR  = VISDIR / "slices"
PV_EXE    = Path("/local/data/public/skcb2/tddft/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch")

for d in (VTIDIR, FRAMESDIR, SLICEDIR):
    d.mkdir(parents=True, exist_ok=True)

# ── 1. Load simulation data ──────────────────────────────────────────────────
sim = SimulationData(RESULTS)

density_series  = sim.field_series("density")
density_meta    = density_series.files[0]
density         = load_real_field(meta_path=density_meta)

wp_series       = sim.field_series("wp_density")
wp_meta         = wp_series.files[0]
wp_density      = load_real_field(meta_path=wp_meta)

# ── 2. Validation ────────────────────────────────────────────────────────────
dV      = density.meta.voxel_volume_bohr3
n_elec  = density.array.sum() * dV
wp_norm = wp_density.array.sum() * dV

print(f"[GS density]   N_electrons = {n_elec:.3f}   (expect ≈ 108.0)")
print(f"[WP density]   norm        = {wp_norm:.4f}   (expect ≈ 1.0)")

PASS = "✓"
FAIL = "✗"
print("\n=== Validation ===")
print(f"  {PASS if abs(n_elec - 108.0) < 1.0 else FAIL}  N_electrons = {n_elec:.3f}  (expect 108 ± 1)")
print(f"  {PASS if abs(wp_norm - 1.0) < 0.05 else FAIL}  WP norm     = {wp_norm:.4f}  (expect 1.0 ± 0.05)")

# ── 3. Read wp_params.txt ────────────────────────────────────────────────────
params_path = RESULTS / "wp_params.txt"
if params_path.exists():
    print("\n=== WP parameters ===")
    for line in params_path.read_text().splitlines():
        print(f"  {line}")

# ── 4. Matplotlib slice plots ────────────────────────────────────────────────
def save_slices(arr, label, filename):
    nx, ny, nz = arr.shape
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(f"run_03_low_energy — {label}")

    im0 = axes[0].imshow(arr[:, ny // 2, :].T, origin="lower", cmap="hot", aspect="auto")
    axes[0].set_title("x-z slice (y=mid)")
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(arr[nx // 2, :, :].T, origin="lower", cmap="hot", aspect="auto")
    axes[1].set_title("y-z slice (x=mid)")
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    im2 = axes[2].imshow(arr[:, :, nz // 2].T, origin="lower", cmap="hot", aspect="auto")
    axes[2].set_title("x-y slice (z=mid, coronene plane)")
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    fig.tight_layout()
    out = SLICEDIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out}")


print("\n=== Matplotlib slices ===")
save_slices(density.array,    "GS total density (bohr⁻³)", "density_slices.png")
save_slices(wp_density.array, "WP orbital density |ψ_wp|²", "wp_density_slices.png")

# ── 5. VTI conversion ────────────────────────────────────────────────────────
print("\n=== VTI conversion ===")
vti_density = convert_real_meta_to_vti(
    density_meta,
    output_path=VTIDIR / "density_total.vti",
    array_name="total_density",
)
print(f"  Saved: {vti_density}")

vti_wp = convert_real_meta_to_vti(
    wp_meta,
    output_path=VTIDIR / "wp_density.vti",
    array_name="wp_density",
)
print(f"  Saved: {vti_wp}")

# ── 6. ParaView renders ──────────────────────────────────────────────────────
print("\n=== ParaView renders ===")
try:
    pv = ParaViewPipeline(pv_executable=PV_EXE)

    pv.render_density_from_meta_series(
        density_series,
        vti_output_dir=VTIDIR / "density_series",
        render=VolumeRenderSpec(array_name="total_density"),
        animation=AnimationSpec(
            output_frames_dir=FRAMESDIR / "density",
            image_size=(1600, 1200)),
    )
    print(f"  Density frames: {FRAMESDIR / 'density'}")

    pv.render_density_from_meta_series(
        wp_series,
        vti_output_dir=VTIDIR / "wp_density_series",
        render=VolumeRenderSpec(array_name="wp_density"),
        animation=AnimationSpec(
            output_frames_dir=FRAMESDIR / "wp_density",
            image_size=(1600, 1200)),
    )
    print(f"  WP density frames: {FRAMESDIR / 'wp_density'}")

except Exception as exc:
    print(f"  ParaView render failed: {exc}", file=sys.stderr)

# ── 7. Output summary ────────────────────────────────────────────────────────
print("\n=== Output files ===")
for p in sorted(VISDIR.rglob("*")):
    if p.is_file():
        print(f"  {p}")
