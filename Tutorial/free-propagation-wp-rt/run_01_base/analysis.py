"""
analysis.py — free-propagation run_01_base post-processing.

Validates:
  1. N_elec = 1.0 per frame
  2. z-profile heatmap; σ(t) vs analytic formula
  3. Observables summary plot
  4. Screen pattern plots (3 screens)
  5. Density GIF via default_density_movie()
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys

# Allow running from repo root or this directory
REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))

import inqview

# ── Paths ─────────────────────────────────────────────────────────────────────
RUN_DIR  = Path(__file__).parent
RESULTS  = RUN_DIR / "results"
VISDIR   = RESULTS / "visualisation"
VISDIR.mkdir(parents=True, exist_ok=True)

# ── Physical constants and run params ─────────────────────────────────────────
LZ_BOHR      = 89.856
DT_AU        = 0.02
N_STEPS      = 10000
WRITE_EVERY  = 100
WP_SIGMA_ANG = 0.53
ANG_TO_BOHR  = 1.8897259886
SIGMA0       = WP_SIGMA_ANG * ANG_TO_BOHR   # initial width in bohr

PV_EXE = REPO_ROOT / "ParaView-6.1.0-MPI-Linux-Python3.12-x86_64" / "bin" / "pvbatch"

# ── 1. Load density series ────────────────────────────────────────────────────
print("Loading density series...")
series = inqview.FieldSeries(RESULTS / "density_rt")
n_frames = len(series.frames)
print(f"  {n_frames} frames loaded")

# ── 2. Validate N_elec per frame ──────────────────────────────────────────────
print("Checking N_elec per frame...")
nelec_vals = []
for i, frame in enumerate(series.frames):
    field = inqview.load_real_field(frame.data_path)
    dx, dy, dz = field.meta.spacing_bohr
    nelec = float(field.array.sum() * dx * dy * dz)
    nelec_vals.append(nelec)
nelec_arr = np.array(nelec_vals)
print(f"  N_elec mean={nelec_arr.mean():.4f} min={nelec_arr.min():.4f} max={nelec_arr.max():.4f}")
assert np.all(np.abs(nelec_arr - 1.0) < 0.005), "N_elec deviated > 0.005 from 1.0"
print("  PASS: N_elec = 1.0 ± 0.005")

# ── 3. z-profile and σ(t) extraction ─────────────────────────────────────────
print("Computing z-profiles and sigma(t)...")
times = [f.time_au for f in series.frames]
sigma_t = []

for i, frame in enumerate(series.frames):
    field = inqview.load_real_field(frame.data_path)
    dx, dy, dz = field.meta.spacing_bohr
    # Integrate over x,y → ρ(z)
    rho_z = field.array.sum(axis=(0, 1)) * dx * dy   # shape (nz,)
    nz = rho_z.shape[0]
    z_coords = field.meta.origin_bohr[2] + np.arange(nz) * dz

    # Second moment: σ² = ∫z²ρ(z)dz − (∫zρ(z)dz)²
    norm = float(rho_z.sum() * dz)
    if norm < 1e-12:
        sigma_t.append(np.nan)
        continue
    z_mean  = float((z_coords * rho_z).sum() * dz) / norm
    z2_mean = float((z_coords**2 * rho_z).sum() * dz) / norm
    sigma_t.append(np.sqrt(max(z2_mean - z_mean**2, 0.0)))

times_arr = np.array(times)
sigma_arr = np.array(sigma_t)

# Analytic prediction: σ(t) = σ0 √(1 + t²/σ0⁴)
sigma_analytic = SIGMA0 * np.sqrt(1.0 + (times_arr / SIGMA0**2)**2)

# z-profile heatmap
print("  Plotting z-profile heatmap...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: σ(t) comparison
axes[0].plot(times_arr, sigma_arr, label="TDDFT σ(t)")
axes[0].plot(times_arr, sigma_analytic, "--", label="Analytic (free particle)")
axes[0].set_xlabel("Time (a.u.)")
axes[0].set_ylabel("σ (bohr)")
axes[0].set_title("WP width σ(t) — run_01_base")
axes[0].legend()

# Right: density vs z vs t (selected frames)
stride = max(1, len(series.frames) // 20)
rho_z_all = []
for i, frame in enumerate(series.frames[::stride]):
    field = inqview.load_real_field(frame.data_path)
    dx, dy, dz = field.meta.spacing_bohr
    rho_z = field.array.sum(axis=(0, 1)) * dx * dy
    rho_z_all.append(rho_z)
rho_z_mat = np.array(rho_z_all)   # (n_sampled, nz)

field0 = inqview.load_real_field(series.frames[0].data_path)
dz0 = field0.meta.spacing_bohr[2]
nz  = rho_z_mat.shape[1]
z_ax = field0.meta.origin_bohr[2] + np.arange(nz) * dz0
t_ax = times_arr[::stride]

im = axes[1].imshow(
    rho_z_mat.T, origin="lower", aspect="auto",
    extent=[t_ax[0], t_ax[-1], z_ax[0], z_ax[-1]],
    cmap="viridis",
)
plt.colorbar(im, ax=axes[1], label="ρ(z) (bohr⁻¹)")
axes[1].set_xlabel("Time (a.u.)")
axes[1].set_ylabel("z (bohr)")
axes[1].set_title("Density z-profile vs time")

fig.tight_layout()
fig.savefig(VISDIR / "sigma_t_and_zprofile.png", dpi=150)
plt.close(fig)
print(f"  Saved {VISDIR / 'sigma_t_and_zprofile.png'}")

# ── 4. Observables summary ────────────────────────────────────────────────────
obs_path = RESULTS / "observables.csv"
if obs_path.exists():
    print("Plotting observables...")
    inqview.plot_observables_summary(obs_path, VISDIR / "observables_summary.png")
    print(f"  Saved {VISDIR / 'observables_summary.png'}")

# ── 5. Screen pattern plots ───────────────────────────────────────────────────
screens_dir = RESULTS / "screens"
screen_out  = VISDIR / "screens"
screen_out.mkdir(parents=True, exist_ok=True)
for dat_file in sorted(screens_dir.glob("*.dat")):
    pattern = inqview.load_leed_pattern(dat_file)
    out_png = screen_out / (dat_file.stem + ".png")
    inqview.plot_leed_pattern(pattern, out_png)
    print(f"  Screen plot: {out_png.name}")

# ── 6. Density GIF ────────────────────────────────────────────────────────────
if PV_EXE.exists():
    print("Generating density GIF...")
    paths = inqview.default_density_movie(series, VISDIR / "density", pv_executable=PV_EXE)
    print(f"  GIF: {paths['gif']}")
else:
    print(f"  Skipping GIF: pvbatch not found at {PV_EXE}")

print("\nAnalysis complete.")
