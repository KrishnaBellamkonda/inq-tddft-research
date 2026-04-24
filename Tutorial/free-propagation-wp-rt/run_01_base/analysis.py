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

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))

import inqview

RUN_DIR  = Path(__file__).parent
RESULTS  = RUN_DIR / "results"
VISDIR   = RESULTS / "visualisation"
VISDIR.mkdir(parents=True, exist_ok=True)

WP_SIGMA_ANG = 0.53
SIGMA0       = WP_SIGMA_ANG * 1.8897259886   # bohr
PV_EXE = REPO_ROOT / "ParaView-6.1.0-MPI-Linux-Python3.12-x86_64" / "bin" / "pvbatch"

# ── 1. Load density series ────────────────────────────────────────────────────
print("Loading density series...")
series = inqview.SimulationData(RUN_DIR).field_series("results/density_rt")
n_frames = len(series.files)
print(f"  {n_frames} frames loaded")

# ── 2. Validate N_elec per frame ──────────────────────────────────────────────
print("Checking N_elec per frame...")
nelec_vals = []
for meta_path in series.files:
    field = inqview.load_real_field(meta_path=meta_path)
    dx, dy, dz = field.meta.spacing_bohr
    nelec_vals.append(float(field.array.sum() * dx * dy * dz))
nelec_arr = np.array(nelec_vals)
print(f"  N_elec mean={nelec_arr.mean():.4f} min={nelec_arr.min():.4f} max={nelec_arr.max():.4f}")
assert np.all(np.abs(nelec_arr - 1.0) < 0.005), "N_elec deviated > 0.005 from 1.0"
print("  PASS: N_elec = 1.0 ± 0.005")

# ── 3. z-profile and σ(t) extraction ─────────────────────────────────────────
print("Computing z-profiles and sigma(t)...")
times, sigma_t = [], []

for meta_path in series.files:
    field = inqview.load_real_field(meta_path=meta_path)
    dx, dy, dz = field.meta.spacing_bohr
    rho_z = field.array.sum(axis=(0, 1)) * dx * dy
    nz = rho_z.shape[0]
    z_coords = field.meta.origin_bohr[2] + np.arange(nz) * dz

    norm = float(rho_z.sum() * dz)
    if norm > 1e-12:
        z_mean  = float((z_coords * rho_z).sum() * dz) / norm
        z2_mean = float((z_coords**2 * rho_z).sum() * dz) / norm
        sigma_t.append(np.sqrt(max(z2_mean - z_mean**2, 0.0)))
    else:
        sigma_t.append(np.nan)
    times.append(field.meta.time_au if field.meta.time_au is not None else np.nan)

times_arr  = np.array(times)
sigma_arr  = np.array(sigma_t)
sigma_analytic = SIGMA0 * np.sqrt(1.0 + (times_arr / SIGMA0**2)**2)

# z-profile heatmap
print("  Plotting sigma(t) and z-profile...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(times_arr, sigma_arr, label="TDDFT sigma(t)")
axes[0].plot(times_arr, sigma_analytic, "--", label="Analytic (free particle)")
axes[0].set_xlabel("Time (a.u.)")
axes[0].set_ylabel("sigma (bohr)")
axes[0].set_title("WP width sigma(t) — run_01_base")
axes[0].legend()

stride = max(1, n_frames // 20)
rho_z_all = []
for meta_path in series.files[::stride]:
    field = inqview.load_real_field(meta_path=meta_path)
    dx, dy, dz = field.meta.spacing_bohr
    rho_z_all.append(field.array.sum(axis=(0, 1)) * dx * dy)
rho_z_mat = np.array(rho_z_all)

field0 = inqview.load_real_field(meta_path=series.files[0])
dz0  = field0.meta.spacing_bohr[2]
nz   = rho_z_mat.shape[1]
z_ax = field0.meta.origin_bohr[2] + np.arange(nz) * dz0
t_ax = times_arr[::stride]

im = axes[1].imshow(
    rho_z_mat.T, origin="lower", aspect="auto",
    extent=[t_ax[0], t_ax[-1], z_ax[0], z_ax[-1]], cmap="viridis",
)
plt.colorbar(im, ax=axes[1], label="rho(z) (bohr^-1)")
axes[1].set_xlabel("Time (a.u.)")
axes[1].set_ylabel("z (bohr)")
axes[1].set_title("Density z-profile vs time")

fig.tight_layout()
fig.savefig(VISDIR / "sigma_t_and_zprofile.png", dpi=150)
plt.close(fig)
print(f"  Saved sigma_t_and_zprofile.png")

# ── 4. Observables summary ────────────────────────────────────────────────────
obs_path = RESULTS / "observables.csv"
if obs_path.exists():
    print("Plotting observables...")
    inqview.plot_observables_summary(obs_path, VISDIR / "observables_summary.png")
    print(f"  Saved observables_summary.png")

# ── 5. Screen pattern plots ───────────────────────────────────────────────────
screens_dir = RESULTS / "screens"
screen_out  = VISDIR / "screens"
screen_out.mkdir(parents=True, exist_ok=True)
for dat_file in sorted(screens_dir.glob("*.dat")):
    pattern = inqview.load_leed_pattern(dat_file)
    inqview.plot_leed_pattern(pattern, screen_out / (dat_file.stem + ".png"))
    print(f"  Screen plot: {dat_file.stem}")

# ── 6. Density GIF ────────────────────────────────────────────────────────────
if PV_EXE.exists():
    print("Generating density GIF...")
    paths = inqview.default_density_movie(series, VISDIR / "density", pv_executable=PV_EXE)
    print(f"  GIF: {paths['gif']}")
else:
    print(f"  Skipping GIF: pvbatch not found at {PV_EXE}")

print("\nAnalysis complete.")
