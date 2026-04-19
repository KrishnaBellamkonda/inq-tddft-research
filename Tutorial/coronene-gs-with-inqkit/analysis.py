"""
coronene-gs-with-inqkit: analysis.py
Validates the inqkit write → inqview read → visualise pipeline for coronene C24H12.

Run with: /local/data/public/skcb2/pyenv/versions/3.10.19/envs/quantum-wave-packet/bin/python3 analysis.py
Expected norms: N_elec ≈ 108, HOMO density norm ≈ 1.0, HOMO |ψ|² norm ≈ 1.0
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.data import SimulationData, load_real_field, load_complex_field
from inqview.vti import convert_real_meta_to_vti, write_vti
from inqview.paraview import ParaViewPipeline, VolumeRenderSpec, AnimationSpec
from inqview.fields import RealField3D

RESULTS   = Path(__file__).parent / "results"
VISDIR    = RESULTS / "visualisation"
VTIDIR    = VISDIR / "vti"
FRAMESDIR = VISDIR / "frames"
PV_EXE    = Path("/local/data/public/skcb2/tddft/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch")

VISDIR.mkdir(parents=True, exist_ok=True)
VTIDIR.mkdir(parents=True, exist_ok=True)

# ── 1. SimulationData handle ─────────────────────────────────────────────────
sim = SimulationData(RESULTS)

# ── 2. Total electron density ────────────────────────────────────────────────
density_series   = sim.field_series("density")
density_meta     = density_series.files[0]
density          = load_real_field(meta_path=density_meta)
dV               = density.meta.voxel_volume_bohr3
n_elec           = density.array.sum() * dV
print(f"[density]   grid shape = {density.array.shape}")
print(f"[density]   N_electrons = {n_elec:.3f}   (expect ≈ 108.0)")

# ── 3. HOMO orbital density (|ψ_53|²) ───────────────────────────────────────
orb_density_series = sim.field_series("orbital_density")
orb_density_meta   = orb_density_series.files[0]
rho_homo           = load_real_field(meta_path=orb_density_meta)
homo_density_norm  = rho_homo.array.sum() * dV
print(f"[HOMO ρ]    density norm = {homo_density_norm:.4f}   (expect ≈ 1.0)")

# ── 4. HOMO complex wavefunction ─────────────────────────────────────────────
orbital_series   = sim.field_series("orbitals")
orbital_meta     = orbital_series.files[0]
psi_homo         = load_complex_field(meta_path=orbital_meta)
psi_norm         = (psi_homo.real**2 + psi_homo.imag**2).sum() * dV
print(f"[HOMO ψ]    |ψ|² norm = {psi_norm:.4f}   (expect ≈ 1.0)")

# ── 5. Validation summary ────────────────────────────────────────────────────
PASS = "\u2713"
FAIL = "\u2717"
print("\n=== Validation summary ===")
print(f"  {PASS if abs(n_elec - 108.0) < 1.0 else FAIL}  N_electrons = {n_elec:.3f}  (expect 108 ± 1)")
print(f"  {PASS if abs(homo_density_norm - 1.0) < 0.05 else FAIL}  HOMO density norm = {homo_density_norm:.4f}  (expect 1.0 ± 0.05)")
print(f"  {PASS if abs(psi_norm - 1.0) < 0.05 else FAIL}  HOMO |ψ|² norm = {psi_norm:.4f}  (expect 1.0 ± 0.05)")

# ── 6. Matplotlib slice plots ────────────────────────────────────────────────
def _save_slices(arr, label, filename):
    nx, ny, nz = arr.shape
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(f"Coronene — {label}")

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
    out = VISDIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out}")


print("\n=== Matplotlib slices ===")
_save_slices(density.array, "total density (bohr⁻³)", "density_slices.png")
_save_slices(rho_homo.array, "HOMO orbital density", "orbital_density_slices.png")
_save_slices(psi_homo.real, "Re(HOMO ψ)", "orbital_psi_re_slices.png")

# ── 7. VTI conversion ────────────────────────────────────────────────────────
print("\n=== VTI conversion ===")
vti_density = convert_real_meta_to_vti(
    density_meta,
    output_path=VTIDIR / "density_total.vti",
    array_name="total_density",
)
print(f"  Saved: {vti_density}")

vti_orb_density = convert_real_meta_to_vti(
    orb_density_meta,
    output_path=VTIDIR / "orbital_density.vti",
    array_name="homo_density",
)
print(f"  Saved: {vti_orb_density}")

# wavefunction magnitude |ψ|²
mag2 = (psi_homo.real**2 + psi_homo.imag**2).astype(np.float64)
mag_field = RealField3D(meta=psi_homo.meta, array=mag2)
vti_psi_mag = write_vti(mag_field, VTIDIR / "orbital_magnitude.vti", array_name="homo_magnitude")
print(f"  Saved: {vti_psi_mag}")

# ── 8. ParaView renders ──────────────────────────────────────────────────────
print("\n=== ParaView renders ===")
try:
    pv = ParaViewPipeline(pv_executable=PV_EXE)

    density_frames = FRAMESDIR / "density"
    pv.render_density_from_meta_series(
        density_series,
        vti_output_dir=VTIDIR / "density_series",
        render=VolumeRenderSpec(array_name="total_density"),
        animation=AnimationSpec(output_frames_dir=density_frames, image_size=(1600, 1200)),
    )
    print(f"  Density frames: {density_frames}")

    orb_frames = FRAMESDIR / "orbital_density"
    pv.render_density_from_meta_series(
        orb_density_series,
        vti_output_dir=VTIDIR / "orbital_density_series",
        render=VolumeRenderSpec(array_name="homo_density"),
        animation=AnimationSpec(output_frames_dir=orb_frames, image_size=(1600, 1200)),
    )
    print(f"  Orbital density frames: {orb_frames}")

    print("  ParaView renders complete.")
except Exception as exc:
    print(f"  ParaView render failed: {exc}", file=sys.stderr)

# ── 9. Output summary ────────────────────────────────────────────────────────
print("\n=== Output files ===")
for p in sorted(VISDIR.rglob("*")):
    if p.is_file():
        print(f"  {p}")
