"""
make_density_gifs.py -- Generate xy-plane density GIFs from run_002 snapshots.

GIFs produced:
  density_animations/total_density_z_flake.gif  -- Total electron density at z=z_flake
  density_animations/delta_density_z_flake.gif  -- Dn(t) = n(t)-n(0) at z=z_flake

z_flake = Lz/2 = 29.952 bohr = 15.85 Ang  (coronene plane, centred geometry)
z_mid   = (z_flake + z_obs)/2 = 35.95 bohr = 19.03 Ang  -- NOT saved in run_002.

Notes:
- WP orbital slice files (wp_orbital/) contain stale GPU-buffer values and are
  NOT used. The total density snapshots are physically correct (confirmed by energy).
- Delta-density Dn = n(t)-n(0) isolates the WP-induced perturbation above the
  static coronene background.
- Cell: 18.4 x 18.4 Ang, grid 100x100, 21 snapshots from t=0 to t=10.0 a.u.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import imageio.v2 as imageio
import io

# ---- Paths ------------------------------------------------------------------
RUN_DIR = Path(__file__).parent
RESULTS = RUN_DIR / "results"
SNAP_DIR = RESULTS / "density_snapshots"
OUT_DIR  = RESULTS / "density_animations"
OUT_DIR.mkdir(exist_ok=True)

# ---- Physical constants -----------------------------------------------------
LX_ANG   = 18.4
LY_ANG   = 18.4
BOHR_TO_ANG = 0.529177210903
AU_TO_FS    = 0.024188843

Z_FLAKE_BOHR = 29.952157
Z_OBS_BOHR   = 41.952157
Z_FLAKE_ANG  = Z_FLAKE_BOHR * BOHR_TO_ANG   # 15.85 Ang
Z_MID_ANG    = (Z_FLAKE_BOHR + Z_OBS_BOHR) / 2 * BOHR_TO_ANG  # 19.03 Ang

# WP designed to arrive at coronene at t1 = 3.18 a.u.
T1_AU = 3.183

# ---- Helpers ----------------------------------------------------------------
def read_snapshot(fpath):
    """Return (t_au, z_bohr, data_2d [NX x NY])."""
    with open(fpath) as f:
        header = f.readline().strip()
        t_au   = float(header.split("t=")[1].split()[0])
        z_bohr = float(header.split("z=")[1])
        data   = np.loadtxt(f)
    return t_au, z_bohr, data


def make_frame(data, t_au, title, label, cmap, vmin, vmax, lx, ly):
    """Render one heatmap frame as PNG bytes."""
    fig, ax = plt.subplots(figsize=(5, 5), dpi=120)
    im = ax.imshow(
        data.T, origin="lower",
        extent=[0, lx, 0, ly],
        cmap=cmap, vmin=vmin, vmax=vmax,
        aspect="equal", interpolation="bilinear",
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label, fontsize=9)
    ax.set_xlabel("x (Ang)", fontsize=10)
    ax.set_ylabel("y (Ang)", fontsize=10)
    t_as = t_au * AU_TO_FS * 1000
    ax.set_title(f"{title}\nt = {t_au:.3f} a.u.  ({t_as:.2f} as)", fontsize=9)
    ax.tick_params(labelsize=8)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.read()


# ---- Load all frames --------------------------------------------------------
print("Reading total density snapshots...")
snap_files  = sorted(SNAP_DIR.glob("snapshot_t*.txt"))
snap_frames = [read_snapshot(f) for f in snap_files]

t_vals    = np.array([f[0] for f in snap_frames])
densities = np.array([f[2] for f in snap_frames])   # (N_frames, NX, NY)
NF, NX, NY = densities.shape
print(f"Grid: {NX} x {NY}, z = {Z_FLAKE_ANG:.2f} Ang")
print(f"Frames: {NF}, t = {t_vals[0]:.3f} - {t_vals[-1]:.3f} a.u.")

# Reference: static ground state (t=0, WP far from coronene plane)
n0 = densities[0]

# Difference density
delta = densities - n0[np.newaxis, :, :]

# Colour scales
n_p995        = np.percentile(densities, 99.5)
delta_abs_max = np.max(np.abs(delta)) * 1.05

print(f"Total density 99.5th pctile: {n_p995:.4e}")
print(f"Delta density abs max:       {delta_abs_max:.4e}")


# ---- GIF 1: Total electron density ------------------------------------------
print("\nGenerating total density GIF (z_flake)...")
imgs1 = []
for i in range(NF):
    t_au, _, data = snap_frames[i]
    frame = make_frame(
        data, t_au,
        title=f"Total density  z = {Z_FLAKE_ANG:.2f} Ang (coronene plane)",
        label="n(x,y) (bohr^-3)",
        cmap="inferno",
        vmin=0, vmax=n_p995,
        lx=LX_ANG, ly=LY_ANG,
    )
    imgs1.append(imageio.imread(io.BytesIO(frame)))

gif1 = OUT_DIR / "total_density_z_flake.gif"
imageio.mimsave(gif1, imgs1, duration=0.35, loop=0)
print(f"  Saved: {gif1}")


# ---- GIF 2: Difference density Dn = n(t) - n(0) ----------------------------
print("\nGenerating difference density GIF (z_flake)...")
imgs2 = []
for i in range(NF):
    t_au = t_vals[i]
    frame = make_frame(
        delta[i], t_au,
        title=f"Dn = n(t)-n(0)  z = {Z_FLAKE_ANG:.2f} Ang (coronene plane)",
        label="Dn(x,y) (bohr^-3)",
        cmap="RdBu_r",
        vmin=-delta_abs_max, vmax=delta_abs_max,
        lx=LX_ANG, ly=LY_ANG,
    )
    imgs2.append(imageio.imread(io.BytesIO(frame)))

gif2 = OUT_DIR / "delta_density_z_flake.gif"
imageio.mimsave(gif2, imgs2, duration=0.35, loop=0)
print(f"  Saved: {gif2}")


# ---- Static frame: side-by-side at t1 (WP at coronene) ---------------------
idx_t1 = int(np.argmin(np.abs(t_vals - T1_AU)))
t_t1   = t_vals[idx_t1]
t_t1_as = t_t1 * AU_TO_FS * 1000

fig, axes = plt.subplots(1, 2, figsize=(11, 5), dpi=130)
extent = [0, LX_ANG, 0, LY_ANG]

im0 = axes[0].imshow(densities[idx_t1].T, origin="lower", extent=extent,
                      cmap="inferno", vmin=0, vmax=n_p995, aspect="equal")
fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04).set_label(
    "n(x,y) (bohr^-3)", fontsize=9)
axes[0].set_title(
    f"Total density  (t1: WP at coronene)\n"
    f"t = {t_t1:.3f} a.u. ({t_t1_as:.1f} as)", fontsize=9)
axes[0].set_xlabel("x (Ang)"); axes[0].set_ylabel("y (Ang)")

d_max = np.max(np.abs(delta[idx_t1]))
im1 = axes[1].imshow(delta[idx_t1].T, origin="lower", extent=extent,
                      cmap="RdBu_r", vmin=-d_max, vmax=d_max, aspect="equal")
fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04).set_label(
    "Dn(x,y) (bohr^-3)", fontsize=9)
axes[1].set_title(
    f"Dn = n(t)-n(0)  (WP-induced change)\n"
    f"t = {t_t1:.3f} a.u. ({t_t1_as:.1f} as)", fontsize=9)
axes[1].set_xlabel("x (Ang)"); axes[1].set_ylabel("y (Ang)")

fig.suptitle(
    f"Coronene plane z = {Z_FLAKE_ANG:.2f} Ang  --  run_002  200 eV  d=1.4 Ang",
    fontsize=11)
fig.tight_layout()
static_path = OUT_DIR / "density_at_t1.png"
fig.savefig(static_path, bbox_inches="tight")
plt.close(fig)
print(f"\n  Static summary (t1 frame): {static_path}")

print(f"\nNote: z_mid = {Z_MID_ANG:.2f} Ang was not saved in run_002 (only z_flake).")
print("Done.")
