#!/usr/bin/env python3
"""Free-Gaussian dispersion validation: WP in vacuum, sigma x energy sweep.

Purpose (user request, 2026-07-29): confirm the engine reproduces free-particle
Gaussian broadening after the CSD3 CUDA-12.1 rebuild, across several starting
widths and several energies.

Physics being tested — for a wavepacket amplitude
    psi(z,0) ~ exp(-(z-z0)^2 / (2 sigma0^2)) * exp(i k0 z)
the DENSITY |psi|^2 has standard deviation

    sigma_d(t) = sqrt( sigma0^2/2 + t^2 / (2 sigma0^2) )      [a.u., m_e = 1]
    R(t) = sigma_d(t)/sigma_d(0) = sqrt(1 + (t/sigma0^2)^2),  tau = sigma0^2

Derivation check (independent of the run.cpp header, which states the same):
sigma_d(0) = sigma0/sqrt(2); Delta_p = hbar/(2 sigma_d(0)) = 1/(sqrt2 sigma0);
free spreading Delta_x(t)^2 = Delta_x(0)^2 + (Delta_p t/m)^2 reproduces it exactly.

THREE independent checks come out of this sweep:
  1. sigma_d(t) tracks the analytic curve for each sigma0.
  2. Broadening is INDEPENDENT of k0 (Galilean invariance) — the sharpest check,
     free once the energy axis is swept.
  3. Centroid z_mean(t) = z0 + k0 t, and E_total conserved (kinetic only).

Measured from results/<run>/raw/observables/wp_real_space_stats.csv, column
sigma_z2 = <z^2> - <z>^2, so sigma_d = sqrt(sigma_z2).

Per .claude/rules/sigma-wp-convention.md every run is LABELLED by the wavepacket
sigma (sigma_WP = WP_SIGMA). The density width sigma_WP/sqrt(2) is derived and
never used as a label.

Run sizing, per (sigma0, k0):
  T        = sqrt(3) * sigma0^2       -> R = 2 (density width exactly doubles)
  N_STEPS  = ceil(T/dt)
  clearance= 5 * sigma_d(T)           (same 5-sigma rule the run.cpp geometry uses)
  LZ       = travel + 2*clearance     (travel = k0*T), rounded up to 5 Bohr
  LPERP    = 2*clearance              (transverse spreads by the same law)
  launch_z = -LZ/2 + clearance
The box must hold the packet for the whole run: sigma_z2 is meaningless once the
density wraps the periodic boundary, so this is a correctness constraint, not a
cosmetic one.

Usage (from this directory):
    ./dispatch_dispersion.py --dry-run     # print the matrix, run nothing
    ./dispatch_dispersion.py               # build once, then run the sweep
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

HA_TO_EV = 27.211386245988

SIGMAS_WP = (1.0, 2.0, 3.0)      # Bohr, wavepacket (amplitude) sigma
ENERGIES_EV = (1.0, 10.0, 100.0)

DT = 0.01
H = 0.4                           # grid spacing; k_max = pi/h = 7.85
R_TARGET = 2.0                    # density width doubles over the run
CLEARANCE_SIGMAS = 5.0
WF_EVERY = 5                      # sigma_z2 sampling cadence


def sigma_d(sigma0: float, t: float) -> float:
    """Density standard deviation of a free Gaussian at time t (a.u.)."""
    return math.sqrt(sigma0 * sigma0 / 2.0 + t * t / (2.0 * sigma0 * sigma0))


def plan_run(sigma0: float, energy_ev: float) -> dict:
    # T such that sigma_d(T)/sigma_d(0) = R_TARGET
    t_total = sigma0 * sigma0 * math.sqrt(R_TARGET**2 - 1.0)
    k0 = math.sqrt(2.0 * energy_ev / HA_TO_EV)
    travel = k0 * t_total
    clearance = CLEARANCE_SIGMAS * sigma_d(sigma0, t_total)

    lz = math.ceil((travel + 2 * clearance) / 5.0) * 5
    lperp = math.ceil((2 * clearance) / 5.0) * 5
    launch_z = -lz / 2.0 + clearance

    # Aliasing guard: the grid must resolve k0 plus the packet's momentum spread.
    dk = 1.0 / (math.sqrt(2.0) * sigma0)
    k_needed = k0 + 4.0 * dk
    k_max = math.pi / H

    return {
        "name": f"disp_sig{sigma0:g}_E{energy_ev:g}",
        "sigma_wp": sigma0,
        "energy_ev": energy_ev,
        "k0": k0,
        "t_total": t_total,
        "n_steps": math.ceil(t_total / DT),
        "lz": float(lz),
        "lperp": float(lperp),
        "launch_z": launch_z,
        "sigma_d0": sigma_d(sigma0, 0.0),
        "sigma_dT": sigma_d(sigma0, t_total),
        "k_margin": k_max / k_needed,
        "mpts": (lz / H) * (lperp / H) ** 2 / 1e6,
    }


def env_for(p: dict) -> dict:
    env = dict(os.environ)
    env.update(
        WP_OUT=f"dispersion/{p['name']}",
        WP_ETA="0",                      # NO CAP -> pure periodic vacuum, stock inq
        WP_SIGMA=f"{p['sigma_wp']:.10g}",
        WP_K0=f"{p['k0']:.10g}",
        WP_LZ=f"{p['lz']:.10g}",
        WP_LPERP=f"{p['lperp']:.10g}",
        WP_LAUNCH_Z=f"{p['launch_z']:.10g}",
        WP_H=f"{H:.10g}",
        WP_DT=f"{DT:.10g}",
        WP_NSTEPS=str(p["n_steps"]),
        WP_WF_EVERY=str(WF_EVERY),
        WP_MOM_EVERY=str(WF_EVERY),
    )
    return env


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print the matrix, run nothing")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    os.chdir(here)

    plans = [plan_run(s, e) for s in SIGMAS_WP for e in ENERGIES_EV]

    hdr = (f"{'run':>22} {'sig_WP':>7} {'E_eV':>6} {'k0':>6} {'T_au':>7} "
           f"{'steps':>6} {'LZ':>5} {'Lperp':>6} {'launch':>7} {'sig_d0':>7} "
           f"{'sig_dT':>7} {'Mpts':>6} {'kguard':>7}")
    print(hdr)
    print("-" * len(hdr))
    for p in plans:
        print(f"{p['name']:>22} {p['sigma_wp']:7.1f} {p['energy_ev']:6.0f} {p['k0']:6.3f} "
              f"{p['t_total']:7.2f} {p['n_steps']:6d} {p['lz']:5.0f} {p['lperp']:6.0f} "
              f"{p['launch_z']:7.1f} {p['sigma_d0']:7.3f} {p['sigma_dT']:7.3f} "
              f"{p['mpts']:6.2f} {p['k_margin']:7.2f}")

    bad = [p for p in plans if p["k_margin"] <= 1.0]
    if bad:
        print("\nABORT: grid too coarse to resolve k0 + 4*dk for: "
              + ", ".join(p["name"] for p in bad), file=sys.stderr)
        return 1

    if args.dry_run:
        print("\n--dry-run: nothing executed.")
        return 0

    # Build once with the first config, then reuse ./run for every point.
    first = env_for(plans[0])
    print(f"\n== building (inq-run, first config: {plans[0]['name']}) ==", flush=True)
    subprocess.run(["inq-run"], env=first, check=True)

    failures = []
    for p in plans[1:]:
        print(f"\n== {p['name']} "
              f"(sigma_WP={p['sigma_wp']}, E={p['energy_ev']} eV, {p['n_steps']} steps) ==",
              flush=True)
        r = subprocess.run(["./run"], env=env_for(p))
        if r.returncode != 0:
            failures.append(p["name"])
            print(f"   FAILED (exit {r.returncode})", file=sys.stderr)

    print("\n== sweep complete ==")
    print(f"   ran {len(plans) - len(failures)}/{len(plans)} runs")
    if failures:
        print("   FAILED: " + ", ".join(failures), file=sys.stderr)
        return 1
    print("   results in results/dispersion/<run>/raw/observables/wp_real_space_stats.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
