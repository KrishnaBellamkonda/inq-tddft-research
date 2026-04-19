"""
N2 real-time density analysis.

Loads the 62-frame density series written by run.cpp, validates electron count
at each frame, plots density z-profiles and 2D slices at selected times, converts
to a VTI series, renders a ParaView animation with CPK nitrogen spheres, and
assembles a GIF.

Run with:
    pyenv activate quantum-wave-packet
    python3 analysis.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from inqview.data import SimulationData, load_real_field, load_meta
from inqview.vti import convert_real_series_to_vti
from inqview.paraview import (
    AnimationSpec,
    AtomSpec,
    ParaViewPipeline,
    VolumeRenderSpec,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
HERE         = Path(__file__).parent.resolve()
RESULTS      = HERE / "results"
DENSITY_DIR  = RESULTS / "real_time" / "density"
VIS_DIR      = RESULTS / "visualisation"
VTI_DIR      = VIS_DIR / "vti"
FRAMES_DIR   = VIS_DIR / "frames"
PV_EXE       = "/local/data/public/skcb2/tddft/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch"

VIS_DIR.mkdir(parents=True, exist_ok=True)

# ── N atom positions in physical frame (symmetric basis, origin at -L/2) ──────
# Cell: L=20 bohr, atoms centred at (0,0,±half_bond_bohr)
# half_bond = 0.575 Å × 1.8897259886 bohr/Å = 1.0866 bohr
_HALF_BOND_BOHR = 0.575 * 1.8897259886
_N_ATOMS = AtomSpec(
    positions=[[0.0, 0.0, -_HALF_BOND_BOHR],
               [0.0, 0.0, +_HALF_BOND_BOHR]],
    symbols=["N", "N"],
    radius_scale=0.35,
    opacity=1.0,
)

# ── 1. Load field series ───────────────────────────────────────────────────────
print("Loading density series …")
sim = SimulationData(RESULTS)
series = sim.field_series("real_time/density")
print(f"  Found {len(series)} metadata files")

# ── 2. Validate N_electrons across all frames ──────────────────────────────────
print("\nValidating electron count …")
times:   list[float] = []
n_elecs: list[float] = []

for meta_path in series.files:
    field = load_real_field(meta_path=meta_path)
    dV    = field.meta.voxel_volume_bohr3
    n_e   = float(field.array.sum()) * dV
    t     = field.meta.time_au or 0.0
    times.append(t)
    n_elecs.append(n_e)

n_arr = np.array(n_elecs)
print(f"  N_electrons: mean={n_arr.mean():.4f}  min={n_arr.min():.4f}  max={n_arr.max():.4f}")
if not np.all(np.abs(n_arr - 10.0) < 0.5):
    print("  WARNING: electron count drifts significantly from 10!")
else:
    print("  OK — all frames within 0.5 of 10 electrons")

# Check time_au is monotonically increasing
times_arr = np.array(times)
if np.all(np.diff(times_arr) > 0):
    print(f"  time_au: monotonically increasing  ({times_arr[0]:.1f} → {times_arr[-1]:.1f} au)")
else:
    print("  WARNING: time_au is not monotonically increasing!")

# ── 3. Compute dipole moment dz(t) = integral z * rho(r,t) dV ─────────────────
# This is the key observable for a z-kick: it oscillates at the N2 electronic
# excitation frequencies. dz(t) - dz(0) reveals the response.
print("\nComputing z-dipole moment …")
dipoles: list[float] = []
meta0 = load_meta(series.files[0])
oz, dz_sp = meta0.origin_bohr[2], meta0.spacing_bohr[2]
nz_pts    = meta0.nz
z_grid    = oz + np.arange(nz_pts) * dz_sp   # physical z coords [bohr]
dx, dy    = meta0.spacing_bohr[0], meta0.spacing_bohr[1]
dV_full   = dx * dy * dz_sp

for meta_path in series.files:
    field = load_real_field(meta_path=meta_path)
    # dz = integral z * rho dV: contract z_grid over the z-axis
    linear_rho = field.array.sum(axis=(0, 1)) * dx * dy   # [e/bohr]
    dz_val     = float((linear_rho * z_grid).sum() * dz_sp)
    dipoles.append(dz_val)

dip_arr  = np.array(dipoles)
dip_ref  = dip_arr[0]
ddip     = dip_arr - dip_ref   # induced dipole (zero at t=0)
print(f"  dz(0) = {dip_ref:.6f} bohr·e")
print(f"  Δdz range: [{ddip.min():.4e}, {ddip.max():.4e}] bohr·e")

fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(times_arr, ddip, lw=1.2, color="steelblue")
ax.axhline(0, color="gray", lw=0.8, ls="--")
ax.set_xlabel("time (au)")
ax.set_ylabel("Δdz (bohr·e)")
ax.set_title("N₂: induced z-dipole moment after kick")
fig.tight_layout()
out = VIS_DIR / "dipole_z_vs_time.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"  Saved: {out}")

# ── 4. Plot N_electrons vs time ────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(times_arr, n_arr, lw=1.2)
ax.axhline(10.0, color="gray", lw=0.8, ls="--", label="expected")
ax.set_xlabel("time (au)")
ax.set_ylabel("N_electrons")
ax.set_title("N₂: electron count conservation")
ax.legend(fontsize=8)
fig.tight_layout()
out = VIS_DIR / "n_electrons_vs_time.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"\n  Saved: {out}")

# ── 4. Z-profile of density (summed over x,y) at several times ─────────────────
# Sample 5 evenly-spaced frames
n_frames   = len(series.files)
sample_idx = np.linspace(0, n_frames - 1, 5, dtype=int)

fig, ax = plt.subplots(figsize=(7, 4))
first_field = load_real_field(meta_path=series.files[sample_idx[0]])
meta0 = first_field.meta
oz    = meta0.origin_bohr[2]
dz    = meta0.spacing_bohr[2]
nz    = meta0.nz
z_coords = oz + np.arange(nz) * dz

for idx in sample_idx:
    field = load_real_field(meta_path=series.files[idx])
    t     = field.meta.time_au or 0.0
    # Sum over x and y, scale by dx*dy to get linear density [e/bohr]
    dx, dy = field.meta.spacing_bohr[0], field.meta.spacing_bohr[1]
    profile = field.array.sum(axis=(0, 1)) * dx * dy
    ax.plot(z_coords, profile, label=f"t={t:.0f} au")

ax.axvline(-_HALF_BOND_BOHR, color="k", ls=":", lw=0.8, label="N atoms")
ax.axvline(+_HALF_BOND_BOHR, color="k", ls=":", lw=0.8)
ax.set_xlabel("z (bohr)")
ax.set_ylabel("linear density (e/bohr)")
ax.set_title("N₂: z-profile of electron density")
ax.legend(fontsize=7)
fig.tight_layout()
out = VIS_DIR / "density_z_profile.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"  Saved: {out}")

# ── 5. 2D x-z density slice at several times ──────────────────────────────────
fig, axes = plt.subplots(1, 5, figsize=(15, 3), sharey=True)
for ax, idx in zip(axes, sample_idx):
    field = load_real_field(meta_path=series.files[idx])
    t     = field.meta.time_au or 0.0
    ny    = field.meta.ny
    # x-z slice at y = ny//2
    sl = field.array[:, ny // 2, :]
    nx = field.meta.nx
    ox = field.meta.origin_bohr[0]
    dx = field.meta.spacing_bohr[0]
    x_coords = ox + np.arange(nx) * dx
    ax.pcolormesh(z_coords, x_coords, sl, cmap="inferno", rasterized=True)
    ax.set_title(f"t={t:.0f} au", fontsize=9)
    ax.set_xlabel("z (bohr)", fontsize=8)
    ax.axvline(-_HALF_BOND_BOHR, color="w", lw=0.5, ls=":")
    ax.axvline(+_HALF_BOND_BOHR, color="w", lw=0.5, ls=":")
axes[0].set_ylabel("x (bohr)", fontsize=8)
fig.suptitle("N₂: x-z density slice (y = centre)", fontsize=10)
fig.tight_layout()
out = VIS_DIR / "density_xz_slices.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"  Saved: {out}")

# ── 6. VTI series ──────────────────────────────────────────────────────────────
print("\nConverting to VTI series …")
vti_result = convert_real_series_to_vti(
    series=series,
    output_dir=VTI_DIR,
    array_name="total_density",
)
print(f"  Written {len(vti_result.files)} VTI files to {VTI_DIR}")
print(f"  Scalar range: [{vti_result.data_min:.4f}, {vti_result.data_max:.4f}]")

# ── 7. ParaView render ─────────────────────────────────────────────────────────
print("\nLaunching ParaView render …")
pv = ParaViewPipeline(pv_executable=PV_EXE)

_rho_max = vti_result.data_max
render_spec = VolumeRenderSpec(
    array_name="total_density",
    # Cap at 40% of max so outer cloud (bonding/lone-pair) occupies the upper opacity range.
    scalar_range=(0.0, _rho_max * 0.4),
    color_preset="Cividis (matplotlib)",
    show_scalar_bar=False,
    camera_azimuth_deg=30.0,
    camera_elevation_deg=20.0,
    # Aggressive opacity ramp — VisRTX needs high values to show anything.
    opacity_points=[
        (0.000,            0.00),
        (0.005 * _rho_max, 0.00),
        (0.030 * _rho_max, 0.15),
        (0.100 * _rho_max, 0.40),
        (0.250 * _rho_max, 0.70),
        (0.400 * _rho_max, 0.90),
    ],
)
anim_spec = AnimationSpec(
    output_frames_dir=FRAMES_DIR,
    image_size=(800, 600),
    frame_stride=1,
    filename_prefix="frame",
)

pv.render_density_from_meta_series(
    series=series,
    vti_output_dir=VTI_DIR,
    render=render_spec,
    animation=anim_spec,
    atoms=_N_ATOMS,
)
print(f"  Frames written to {FRAMES_DIR}")

# ── 8. Build GIF ───────────────────────────────────────────────────────────────
print("\nBuilding GIF …")
gif_path = VIS_DIR / "n2_density_rt.gif"
pv.build_gif(frames_dir=FRAMES_DIR, output_path=gif_path, fps=10)
print(f"  GIF: {gif_path}")

# ── Summary ────────────────────────────────────────────────────────────────────
print("\n── Summary ───────────────────────────────────────────────────────────────")
print(f"  Frames loaded       : {n_frames}")
print(f"  Time span           : {times_arr[0]:.1f} – {times_arr[-1]:.1f} au")
print(f"  N_electrons (mean)  : {n_arr.mean():.4f}")
print(f"  Figures             : {VIS_DIR}")
print(f"  VTI series          : {VTI_DIR}  ({len(vti_result.files)} files)")
print(f"  ParaView frames     : {FRAMES_DIR}")
print(f"  GIF                 : {gif_path}")
