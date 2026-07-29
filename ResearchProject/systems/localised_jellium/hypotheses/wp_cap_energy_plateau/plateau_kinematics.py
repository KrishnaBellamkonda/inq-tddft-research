"""Orbital-free lobe/TOF kinematics extractor for the wp_cap_energy_plateau pair.

Same design as hypotheses/qsp_phase5/qsp5_momentum_kinematics.py, adapted to the
sigma=1, E=100 eV cap/nocap geometry (box 25x25x140, slab |z|<=12.5, launch
z=-20.5, CAP 60<|z|<70 in the cap run only).

One pass over each run's ``density_total`` frames -> cache npz per run with,
per frame: z-profiles of Delta n = n_total - n_gs and of the WP-orbital density
(cross-check channel), plus vacuum-lobe von Weizsaecker T_W and lobe moments.

Run:  PYTHONPATH=<stack> venv/python3 plateau_kinematics.py [cap|nocap ...]
"""
import os
import sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SYS_DIR = os.path.normpath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(SYS_DIR, "scripts", "wp_cap_energy_plateau", "wp", "results")
GS_VTI = os.path.join(SYS_DIR, "shared_gs", "slab_n102_L25x25x140_w0p5_h0p5",
                      "density_gs_system", "density_gs_system.vti")
CACHE = os.path.join(HERE, "cache")

RUNS = ["cap", "nocap"]
DT_AU = 0.02            # propagation time step
DX = 0.5                # grid spacing [Bohr]
Z_B = 16.5              # lobe boundary [Bohr] (slab 12.5 + 4 buffer)
N_FLOOR = 1e-10


def _tw_lobe(dn, mask_z, dx):
    """von Weizsaecker T_W = int |grad n|^2/(8n) over a z-mask (positive Dn only)."""
    sub = dn[:, :, mask_z]
    n = np.where(sub > N_FLOOR, sub, np.nan)
    gx, gy, gz = np.gradient(sub, dx, edge_order=1)
    g2 = gx*gx + gy*gy + gz*gz
    integrand = np.where(np.isnan(n), 0.0, g2/(8.0*np.where(np.isnan(n), 1.0, n)))
    return float(integrand.sum()*dx**3)


def extract_run(run):
    from inqview import load_vti

    raw = os.path.join(RESULTS, run, "raw", "vti")
    n_gs = load_vti(GS_VTI).data
    frames = sorted(os.listdir(os.path.join(raw, "density_total")))
    steps = np.array([int(f.split("_t")[1].split(".")[0]) for f in frames])
    t_au = steps*DT_AU

    first = load_vti(os.path.join(raw, "density_total", frames[0]))
    z = first.z
    lo, hi = z < -Z_B, z > Z_B

    out = {k: np.zeros(len(frames)) for k in
           ("TW_lo", "TW_hi", "N_lo", "N_hi", "Z1_lo", "Z1_hi")}
    rho = np.zeros((len(frames), len(z)))
    rho_wp = np.zeros((len(frames), len(z)))

    for i, f in enumerate(frames):
        dn = load_vti(os.path.join(raw, "density_total", f)).data - n_gs
        rho[i] = dn.sum(axis=(0, 1))*DX*DX
        w = load_vti(os.path.join(raw, "density_wp", f)).data
        rho_wp[i] = w.sum(axis=(0, 1))*DX*DX
        for side, m in (("lo", lo), ("hi", hi)):
            out["TW_"+side][i] = _tw_lobe(dn, m, DX)
            out["N_"+side][i] = float(rho[i][m].sum()*DX)
            out["Z1_"+side][i] = float((rho[i][m]*z[m]).sum()*DX)
        if i % 40 == 0:
            print(f"[{run}] frame {i}/{len(frames)}", flush=True)

    os.makedirs(CACHE, exist_ok=True)
    np.savez_compressed(os.path.join(CACHE, f"{run}_kinematics.npz"),
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
