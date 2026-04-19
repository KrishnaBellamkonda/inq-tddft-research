"""
N2 ground-state analysis.

Loads the density and orbital files written by run.cpp, converts them to VTI,
renders volume images via ParaView, and saves diagnostic plots.

Run with:
    python analysis.py

Outputs go to results/visualisation/.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from inqview.data import SimulationData, load_real_field, load_complex_field
from inqview.vti import convert_real_meta_to_vti
from inqview.paraview import AnimationSpec, AtomSpec, ParaViewPipeline, VolumeRenderSpec

# ── Paths ──────────────────────────────────────────────────────────────────────
HERE    = Path(__file__).parent.resolve()
RESULTS = HERE / "results"
VIS_DIR = RESULTS / "visualisation"
VTI_DIR = VIS_DIR / "vti"
PV_EXE  = "/local/data/public/skcb2/tddft/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch"

VTI_DIR.mkdir(parents=True, exist_ok=True)

# ── N atom positions in simulation-box frame ───────────────────────────────────
# Cell: L = 30 bohr, origin at (0,0,0).
# run.cpp inserts atoms at (0,0,+half_bond) and (0,0,-half_bond) where
# half_bond = 1.10 Å / 2 = 0.55 Å. INQ wraps the negative-z atom to
# z = L - 0.55*1.8897 bohr ≈ 28.96 bohr (the known finite-cell wrap issue).
_ANG_TO_BOHR = 1.8897259886
_HALF_BOND_BOHR = 0.55 * _ANG_TO_BOHR   # ≈ 1.039 bohr
_L_BOHR = 30.0
_N_ATOMS = AtomSpec(
    positions=[
        [0.0, 0.0, _HALF_BOND_BOHR],          # positive-z atom
        [0.0, 0.0, _L_BOHR - _HALF_BOND_BOHR], # wrapped atom (near opposite face)
    ],
    symbols=["N", "N"],
    radius_scale=0.35,
    opacity=1.0,
)

pv = ParaViewPipeline(pv_executable=PV_EXE)

# ── Helper: standard volume render spec for a real density ─────────────────────
def _render_spec(array_name: str, vmax: float) -> VolumeRenderSpec:
    return VolumeRenderSpec(
        array_name=array_name,
        scalar_range=(0.0, vmax * 0.4),
        color_preset="Cividis (matplotlib)",
        show_scalar_bar=False,
        camera_azimuth_deg=30.0,
        camera_elevation_deg=20.0,
        opacity_points=[
            (0.000,         0.00),
            (0.005 * vmax,  0.00),
            (0.030 * vmax,  0.15),
            (0.100 * vmax,  0.40),
            (0.250 * vmax,  0.70),
            (0.400 * vmax,  0.90),
        ],
    )

# ── 1. Total density ───────────────────────────────────────────────────────────
print("── Total density ─────────────────────────────────────────────────────────")
density_meta = RESULTS / "density" / "density_total.meta.txt"
density_vti  = VTI_DIR / "total_density" / "density_total.vti"

vti_path = convert_real_meta_to_vti(
    meta_path=density_meta,
    output_path=density_vti,
    array_name="total_density",
)
print(f"  VTI written: {vti_path}")

rho = load_real_field(meta_path=density_meta)
dV  = rho.meta.voxel_volume_bohr3
n_e = float(rho.array.sum()) * dV
print(f"  N_electrons = {n_e:.4f}  (expected 10)")
vmax_rho = float(rho.array.max())
print(f"  Density max = {vmax_rho:.4f} e/bohr³")

frames_dir = VIS_DIR / "frames" / "total_density"
pv.render_vti_series(
    vti_files=[vti_path],
    render=_render_spec("total_density", vmax_rho),
    animation=AnimationSpec(
        output_frames_dir=frames_dir,
        image_size=(800, 600),
        filename_prefix="frame",
    ),
    atoms=_N_ATOMS,
)
frames = sorted(frames_dir.glob("frame_*.png"))
print(f"  Rendered {len(frames)} frame(s) → {frames_dir}")

# ── 2. Orbital density |ψ₀|² ──────────────────────────────────────────────────
print("\n── Orbital density (KS orbital 0) ────────────────────────────────────────")
orb_dens_meta = RESULTS / "orbital_density" / "orbital_0000_density.meta.txt"
orb_dens_vti  = VTI_DIR / "orbital_density" / "orbital_0000_density.vti"

vti_path_od = convert_real_meta_to_vti(
    meta_path=orb_dens_meta,
    output_path=orb_dens_vti,
    array_name="orbital_density",
)
print(f"  VTI written: {vti_path_od}")

rho_orb = load_real_field(meta_path=orb_dens_meta)
vmax_od = float(rho_orb.array.max())
print(f"  Orbital density max = {vmax_od:.4f} e/bohr³")

frames_dir_od = VIS_DIR / "frames" / "orbital_density"
pv.render_vti_series(
    vti_files=[vti_path_od],
    render=_render_spec("orbital_density", vmax_od),
    animation=AnimationSpec(
        output_frames_dir=frames_dir_od,
        image_size=(800, 600),
        filename_prefix="frame",
    ),
    atoms=_N_ATOMS,
)
frames_od = sorted(frames_dir_od.glob("frame_*.png"))
print(f"  Rendered {len(frames_od)} frame(s) → {frames_dir_od}")

# ── 3. Complex KS orbital ψ₀ — diagnostic plots ────────────────────────────────
print("\n── Complex orbital ψ₀ — diagnostic plots ─────────────────────────────────")
orb_meta = RESULTS / "orbitals" / "orbital_0000.meta.txt"
orbital  = load_complex_field(meta_path=orb_meta)

psi     = orbital.values                      # complex (nx, ny, nz)
density = np.abs(psi) ** 2
print(f"  Shape       : {orbital.shape}")
print(f"  Spacing     : {orbital.spacing_bohr} bohr")
print(f"  Orbital idx : {orbital.orbital_index}   spin: {orbital.spin_index}")
print(f"  Norm (sum*dV): {density.sum() * np.prod(orbital.spacing_bohr):.6f}")

meta   = rho.meta
oz, dz_sp, nz = meta.origin_bohr[2], meta.spacing_bohr[2], meta.nz
dx_sp, dy_sp  = meta.spacing_bohr[0], meta.spacing_bohr[1]
z_grid = oz + np.arange(nz) * dz_sp

# z-profile of real, imag, and |ψ|²
fig, axes = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
for label, arr, ax in [
    ("Re ψ₀",   psi.real,  axes[0]),
    ("Im ψ₀",   psi.imag,  axes[1]),
    ("|ψ₀|²",  density,   axes[2]),
]:
    profile = arr.sum(axis=(0, 1)) * dx_sp * dy_sp
    ax.plot(z_grid, profile, lw=1.2)
    ax.axhline(0, color="gray", lw=0.6, ls="--")
    ax.set_ylabel(label)
    ax.axvline(_HALF_BOND_BOHR,          color="k", lw=0.8, ls=":", label="N atom")
    ax.axvline(_L_BOHR - _HALF_BOND_BOHR, color="k", lw=0.8, ls=":")
axes[-1].set_xlabel("z (bohr)")
axes[0].set_title("N₂: KS orbital 0 — z-profiles")
axes[0].legend(fontsize=7)
fig.tight_layout()
out = VIS_DIR / "orbital_z_profile.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"  Saved: {out}")

# x-z density slice at y = ny//2
ny = orbital.shape[1]
nx = orbital.shape[0]
ox = meta.origin_bohr[0]
x_grid = ox + np.arange(nx) * dx_sp
density_xz = density[:, ny // 2, :]

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for ax, (label, arr) in zip(axes, [
    ("Re ψ₀", psi.real[:, ny // 2, :]),
    ("Im ψ₀", psi.imag[:, ny // 2, :]),
    ("|ψ₀|²", density_xz),
]):
    im = ax.pcolormesh(z_grid, x_grid, arr, cmap="coolwarm" if "Re" in label or "Im" in label else "inferno", rasterized=True)
    ax.set_title(label)
    ax.set_xlabel("z (bohr)")
    ax.set_ylabel("x (bohr)")
    ax.axvline(_HALF_BOND_BOHR,          color="w" if "|" in label else "k", lw=0.7, ls=":")
    ax.axvline(_L_BOHR - _HALF_BOND_BOHR, color="w" if "|" in label else "k", lw=0.7, ls=":")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
fig.suptitle("N₂: KS orbital 0 — x-z slice (y = centre)", fontsize=10)
fig.tight_layout()
out = VIS_DIR / "orbital_xz_slice.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"  Saved: {out}")

# ── Summary ────────────────────────────────────────────────────────────────────
print("\n── Summary ───────────────────────────────────────────────────────────────")
print(f"  N_electrons (total density) : {n_e:.4f}")
print(f"  Total density VTI           : {density_vti}")
print(f"  Orbital density VTI         : {orb_dens_vti}")
print(f"  ParaView frames             : {VIS_DIR}/frames/")
print(f"  Diagnostic plots            : {VIS_DIR}/")
