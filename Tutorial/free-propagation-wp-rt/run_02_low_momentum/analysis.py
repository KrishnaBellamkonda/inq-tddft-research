"""analysis.py — run_02_low_momentum post-processing (sigma=0.53A, 50 eV)."""
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

SIGMA0   = 0.53 * 1.8897259886   # bohr
PV_EXE   = REPO_ROOT / "ParaView-6.1.0-MPI-Linux-Python3.12-x86_64" / "bin" / "pvbatch"

print("Loading density series...")
series = inqview.SimulationData(RUN_DIR).field_series("results/density_rt")
print(f"  {len(series.files)} frames loaded")

# Validate N_elec
nelec_vals = []
for meta_path in series.files:
    f = inqview.load_real_field(meta_path=meta_path)
    dx, dy, dz = f.meta.spacing_bohr
    nelec_vals.append(float(f.array.sum() * dx * dy * dz))
nelec_arr = np.array(nelec_vals)
print(f"  N_elec mean={nelec_arr.mean():.4f} min={nelec_arr.min():.4f}")
assert np.all(np.abs(nelec_arr - 1.0) < 0.005), "N_elec check failed"
print("  PASS: N_elec = 1.0")

# sigma(t) from z-profile
times, sigma_t = [], []
for meta_path in series.files:
    f = inqview.load_real_field(meta_path=meta_path)
    dx, dy, dz = f.meta.spacing_bohr
    rho_z = f.array.sum(axis=(0, 1)) * dx * dy
    nz = rho_z.shape[0]
    z  = f.meta.origin_bohr[2] + np.arange(nz) * dz
    norm = float(rho_z.sum() * dz)
    if norm > 1e-12:
        z_mean = float((z * rho_z).sum() * dz) / norm
        z2_mean = float((z**2 * rho_z).sum() * dz) / norm
        sigma_t.append(np.sqrt(max(z2_mean - z_mean**2, 0.0)))
    else:
        sigma_t.append(np.nan)
    times.append(inqview.load_real_field(meta_path=meta_path).meta.time_au)

t = np.array(times)
s = np.array(sigma_t)
s_an = SIGMA0 * np.sqrt(1.0 + (t / SIGMA0**2)**2)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t, s, label="TDDFT sigma(t)")
ax.plot(t, s_an, "--", label="Analytic (free particle)")
ax.set_xlabel("Time (a.u.)"); ax.set_ylabel("sigma (bohr)")
ax.set_title("run_02_low_momentum: WP width sigma(t)")
ax.legend()
fig.tight_layout()
fig.savefig(VISDIR / "sigma_t.png", dpi=150)
plt.close(fig)

# Observables
obs_path = RESULTS / "observables.csv"
if obs_path.exists():
    inqview.plot_observables_summary(obs_path, VISDIR / "observables_summary.png")
    print("  Observables summary saved")

# Screen plots
screens_dir = RESULTS / "screens"
(VISDIR / "screens").mkdir(exist_ok=True)
for dat in sorted(screens_dir.glob("*.dat")):
    pattern = inqview.load_leed_pattern(dat)
    inqview.plot_leed_pattern(pattern, VISDIR / "screens" / (dat.stem + ".png"))
    print(f"  Screen: {dat.stem}")

# GIF
if PV_EXE.exists():
    paths = inqview.default_density_movie(series, VISDIR / "density", pv_executable=PV_EXE)
    print(f"  GIF: {paths['gif']}")
else:
    print("  Skipping GIF: pvbatch not found")

print("Done.")
