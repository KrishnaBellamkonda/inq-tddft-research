"""analysis.py — coronene run_04_d15 (D=15.0A, sigma=0.53A, 200 eV)."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys

REPO_ROOT = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))
import inqview

RUN_DIR  = Path(__file__).parent
RESULTS  = RUN_DIR / "results"
VISDIR   = RESULTS / "visualisation"
VISDIR.mkdir(parents=True, exist_ok=True)
PV_EXE   = REPO_ROOT / "ParaView-6.1.0-MPI-Linux-Python3.12-x86_64" / "bin" / "pvbatch"

# N_elec expected = 54 coronene occupied + 1 WP = 109 (extra_states(3), 2 unoccupied)
NELEC_EXPECTED = 109

print("Loading density series...")
series = inqview.SimulationData(RUN_DIR).field_series("results/density_rt")
print(f"  {len(series.files)} frames loaded")

nelec_vals = []
for meta_path in series.files:
    f = inqview.load_real_field(meta_path=meta_path)
    dx, dy, dz = f.meta.spacing_bohr
    nelec_vals.append(float(f.array.sum() * dx * dy * dz))
nelec_arr = np.array(nelec_vals)
print(f"  N_elec mean={nelec_arr.mean():.4f} min={nelec_arr.min():.4f}")
assert np.all(np.abs(nelec_arr - NELEC_EXPECTED) < 0.1), f"N_elec deviated from {NELEC_EXPECTED}"
print(f"  PASS: N_elec = {NELEC_EXPECTED}")

# Observables
obs_path = RESULTS / "observables.csv"
if obs_path.exists():
    inqview.plot_observables_summary(obs_path, VISDIR / "observables_summary.png")
    print("  Observables summary saved")

# Screen plots (4 screens)
screens_dir = RESULTS / "screens"
(VISDIR / "screens").mkdir(exist_ok=True)
screen_patterns = {}
for dat in sorted(screens_dir.glob("*.dat")):
    pattern = inqview.load_leed_pattern(dat)
    screen_patterns[dat.stem] = pattern
    inqview.plot_leed_pattern(pattern, VISDIR / "screens" / (dat.stem + ".png"))
    print(f"  Screen: {dat.stem}  (z={pattern.z_bohr:.2f} bohr)")

# 4-panel comparison
if len(screen_patterns) == 4:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (name, pat) in zip(axes.flat, screen_patterns.items()):
        im = ax.imshow(pat.data, origin="lower", extent=pat.extent_bohr,
                       cmap="viridis", aspect="equal")
        plt.colorbar(im, ax=ax, label="rho*dt")
        ax.set_title(f"{name}  z={pat.z_bohr:.1f} bohr")
        ax.set_xlabel("x (bohr)"); ax.set_ylabel("y (bohr)")
    fig.suptitle(f"coronene run_04_d15: all screens")
    fig.tight_layout()
    fig.savefig(VISDIR / "screens" / "4panel_comparison.png", dpi=150)
    plt.close(fig)
    print("  4-panel comparison saved")

# GIF
if PV_EXE.exists():
    paths = inqview.default_density_movie(series, VISDIR / "density", pv_executable=PV_EXE)
    print(f"  GIF: {paths['gif']}")
else:
    print("  Skipping GIF: pvbatch not found")

print("Done.")
