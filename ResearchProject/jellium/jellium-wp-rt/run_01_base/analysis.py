"""analysis.py — jellium run_01_base (sigma=0.53A, 200 eV, N_elec_expected=41)."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys

REPO_ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))
import inqview

RUN_DIR  = Path(__file__).parent
RESULTS  = RUN_DIR / "results"
VISDIR   = RESULTS / "visualisation"
VISDIR.mkdir(parents=True, exist_ok=True)
PV_EXE   = REPO_ROOT / "ParaView-6.1.0-MPI-Linux-Python3.12-x86_64" / "bin" / "pvbatch"

# N_elec expected = 40 jellium + 1 WP = 41
NELEC_EXPECTED = 41

print("Loading density series...")
series = inqview.FieldSeries(RESULTS / "density_rt")
print(f"  {len(series.frames)} frames loaded")

nelec_vals = []
for frame in series.frames:
    f = inqview.load_real_field(frame.data_path)
    dx, dy, dz = f.meta.spacing_bohr
    nelec_vals.append(float(f.array.sum() * dx * dy * dz))
nelec_arr = np.array(nelec_vals)
print(f"  N_elec mean={nelec_arr.mean():.4f} min={nelec_arr.min():.4f}")
assert np.all(np.abs(nelec_arr - NELEC_EXPECTED) < 0.05), f"N_elec deviated from {NELEC_EXPECTED}"
print(f"  PASS: N_elec = {NELEC_EXPECTED}")

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
