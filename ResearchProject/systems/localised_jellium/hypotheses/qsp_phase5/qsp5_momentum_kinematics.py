"""Orbital-free lobe/TOF kinematics extractor for the qsp_phase5 WP sweep.

One 3D pass over each run's ``density_total`` VTI frames producing a compact
per-run cache (``cache/<run>_kinematics.npz``) with, per frame:

- ``rho``      : z-profile of Delta n = n_total(t) - n_gs   [electrons/Bohr]
- ``rho_wp``   : z-profile of the WP-orbital density (cross-check channel)
- ``TW_lo/hi`` : von Weizsaecker shape energy of Delta n over the vacuum lobes
                 z < -Z_B and z > +Z_B  (T_W = int |grad n|^2 / (8 n) d3r) [Ha]
- ``N_lo/hi``  : lobe norms  int_lobe Delta n d3r            [electrons]
- ``Z1_lo/hi`` : lobe first moments int_lobe z Delta n d3r   [electron Bohr]

The z-profiles alone support the full time-of-flight (TOF) flux reconstruction
(1D continuity + CAP sink); the 3D lobe integrals give the localisation (shape)
energy channel.  Everything is computed from n(r,t) only - no Kohn-Sham orbital
identification is used (density_wp is carried purely as a cross-check).

Run:  PYTHONPATH=<stack> venv/python3 qsp5_momentum_kinematics.py [run ...]
      (no args = all five runs, one worker process per run)
"""
import os
import sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(
    HERE, "..", "..", "scripts", "qsp_phase5", "wp", "results"))
CACHE = os.path.join(HERE, "cache")

RUNS = ["p5_wp_v1p3", "p5_wp_v3p0", "p5_wp_v4p0", "p5_wp_v5p0", "p5_wp_v6p0"]
DT_AU = 0.04            # propagation time step (all runs)
DX = 0.5                # grid spacing [Bohr]
Z_B = 15.5              # lobe / detector-plane boundary [Bohr] (slab 12.5 + 3)
N_FLOOR = 1e-10         # density floor for the T_W integrand [e/Bohr^3]


def _tw_lobe(dn, mask_z, dx):
    """von Weizsaecker T_W = int |grad n|^2/(8n) over a z-mask, for n = dn>=floor.

    Negative Delta-n pockets (bath polarisation leakage) are excluded: the shape
    energy is only defined for the (positive) projectile density; the excluded
    weight is small in the buffered vacuum lobes and is QC'd in the notebook.
    """
    sub = dn[:, :, mask_z]
    n = np.where(sub > N_FLOOR, sub, np.nan)
    gx, gy, gz = np.gradient(sub, dx, edge_order=1)
    g2 = gx * gx + gy * gy + gz * gz
    integrand = np.where(np.isnan(n), 0.0, g2 / (8.0 * np.where(np.isnan(n), 1.0, n)))
    return float(integrand.sum() * dx ** 3)


def extract_run(run):
    from inqview import load_vti   # import inside worker

    raw = os.path.join(RESULTS, run, "raw", "vti")
    gs = load_vti(os.path.join(raw, "density_gs_system", "density_gs_system.vti"))
    n_gs = gs.data
    z = gs.z
    lo, hi = z < -Z_B, z > Z_B

    frames = sorted(os.listdir(os.path.join(raw, "density_total")))
    steps = np.array([int(f.split("_t")[1].split(".")[0]) for f in frames])
    t_au = steps * DT_AU

    nz = len(z)
    out = {k: np.zeros(len(frames)) for k in
           ("TW_lo", "TW_hi", "N_lo", "N_hi", "Z1_lo", "Z1_hi")}
    rho = np.zeros((len(frames), nz))
    rho_wp = np.zeros((len(frames), nz))

    for i, f in enumerate(frames):
        dn = load_vti(os.path.join(raw, "density_total", f)).data - n_gs
        rho[i] = dn.sum(axis=(0, 1)) * DX * DX
        w = load_vti(os.path.join(raw, "density_wp", f)).data
        rho_wp[i] = w.sum(axis=(0, 1)) * DX * DX
        for side, m in (("lo", lo), ("hi", hi)):
            out["TW_" + side][i] = _tw_lobe(dn, m, DX)
            out["N_" + side][i] = float(rho[i][m].sum() * DX)
            out["Z1_" + side][i] = float((rho[i][m] * z[m]).sum() * DX)
        if i % 40 == 0:
            print(f"[{run}] frame {i}/{len(frames)}", flush=True)

    os.makedirs(CACHE, exist_ok=True)
    np.savez_compressed(
        os.path.join(CACHE, f"{run}_kinematics.npz"),
        z=z, t_au=t_au, steps=steps, rho=rho, rho_wp=rho_wp,
        dx=DX, z_b=Z_B, **out)
    print(f"[{run}] DONE -> cache/{run}_kinematics.npz", flush=True)
    return run


if __name__ == "__main__":
    todo = sys.argv[1:] or RUNS
    todo = [r for r in todo
            if not os.path.exists(os.path.join(CACHE, f"{r}_kinematics.npz"))]
    print("extracting:", todo or "(nothing - all cached)")
    if todo:
        with ProcessPoolExecutor(max_workers=len(todo)) as ex:
            list(ex.map(extract_run, todo))
