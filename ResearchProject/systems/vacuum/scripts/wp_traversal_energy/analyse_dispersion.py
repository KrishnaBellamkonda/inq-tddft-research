#!/usr/bin/env python3
"""Check measured free-Gaussian broadening against the analytic law.

Consumes the sweep produced by dispatch_dispersion.py and reports, per run, the
worst relative deviation of the measured density width from

    sigma_d(t) = sqrt( sigma0^2/2 + t^2/(2 sigma0^2) )      [a.u., m_e = 1]

Three checks, in decreasing order of how easy they are to fake:

  LONGITUDINAL  sqrt(sigma_z2) vs the law           — the packet is moving in z
  TRANSVERSE    sqrt(sigma_x2) vs the SAME law      — k0 does not enter x at all,
                so this isolates dispersion from anything velocity-dependent
  k0-INDEPENDENCE  at fixed sigma0, the three energies must broaden identically
                (Galilean invariance); reported as the spread across energies

Pure standard library on purpose — the repo venv is unusable (Python 3.6, while
inq-stack needs >=3.10) and the system pythons have no numpy.

Usage (from this directory):
    ./analyse_dispersion.py
    ./analyse_dispersion.py --tol 0.02       # pass threshold on relative error
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

HA_TO_EV = 27.211386245988
RESULTS = Path("results/dispersion")


def sigma_d(sigma0: float, t: float) -> float:
    return math.sqrt(sigma0 * sigma0 / 2.0 + t * t / (2.0 * sigma0 * sigma0))


def read_stats(path: Path) -> tuple[list[float], dict[str, list[float]]]:
    """Return (times, {column: values}) from a WPRealSpaceStats CSV.

    The first line is a '# wp_state_index=... write_every=...' comment; the
    second is the real header.
    """
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    rdr = csv.DictReader(lines)
    cols: dict[str, list[float]] = {}
    times: list[float] = []
    for row in rdr:
        times.append(float(row["time_au"]))
        for k, v in row.items():
            if k in ("step", "time_au") or v is None:
                continue
            try:
                cols.setdefault(k, []).append(float(v))
            except ValueError:
                pass
    return times, cols


def parse_run(name: str) -> tuple[float, float]:
    """disp_sig<S>_E<E> -> (sigma0, energy_eV)."""
    body = name[len("disp_sig"):]
    s_txt, e_txt = body.split("_E")
    return float(s_txt), float(e_txt)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=0.02,
                    help="max allowed relative deviation (default 2%%)")
    args = ap.parse_args()

    if not RESULTS.is_dir():
        print(f"no sweep found at {RESULTS.resolve()}", file=sys.stderr)
        return 1

    runs = sorted(RESULTS.iterdir(), key=lambda p: p.name)
    rows = []
    for run in runs:
        csv_path = run / "raw" / "observables" / "wp_real_space_stats.csv"
        if not csv_path.is_file():
            print(f"  {run.name}: MISSING {csv_path}", file=sys.stderr)
            continue

        sigma0, energy_ev = parse_run(run.name)
        times, cols = read_stats(csv_path)
        if not times:
            print(f"  {run.name}: EMPTY", file=sys.stderr)
            continue

        worst_z = worst_x = 0.0
        for i, t in enumerate(times):
            ref = sigma_d(sigma0, t)
            worst_z = max(worst_z, abs(math.sqrt(cols["sigma_z2"][i]) - ref) / ref)
            worst_x = max(worst_x, abs(math.sqrt(cols["sigma_x2"][i]) - ref) / ref)

        rows.append({
            "name": run.name, "sigma0": sigma0, "energy_ev": energy_ev,
            "t_end": times[-1], "n": len(times),
            "sd_meas_end": math.sqrt(cols["sigma_z2"][-1]),
            "sd_ref_end": sigma_d(sigma0, times[-1]),
            "worst_z": worst_z, "worst_x": worst_x,
        })

    if not rows:
        print("no usable runs", file=sys.stderr)
        return 1

    hdr = (f"{'run':>22} {'sig_WP':>7} {'E_eV':>6} {'t_end':>7} {'pts':>5} "
           f"{'sd_meas':>8} {'sd_ref':>8} {'err_z':>8} {'err_x':>8}  verdict")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        ok = max(r["worst_z"], r["worst_x"]) <= args.tol
        print(f"{r['name']:>22} {r['sigma0']:7.1f} {r['energy_ev']:6.0f} "
              f"{r['t_end']:7.2f} {r['n']:5d} {r['sd_meas_end']:8.4f} "
              f"{r['sd_ref_end']:8.4f} {r['worst_z']:7.3%} {r['worst_x']:7.3%}  "
              f"{'PASS' if ok else 'FAIL'}")

    # k0-independence: at fixed sigma0 the broadening must not depend on energy.
    print("\nk0-independence (Galilean invariance) at fixed sigma_WP:")
    print(f"  {'sig_WP':>7} {'sd_end spread':>14} {'relative':>10}  verdict")
    worst_gal = 0.0
    for sigma0 in sorted({r["sigma0"] for r in rows}):
        grp = [r for r in rows if r["sigma0"] == sigma0]
        if len(grp) < 2:
            continue
        vals = [r["sd_meas_end"] for r in grp]
        spread = max(vals) - min(vals)
        rel = spread / (sum(vals) / len(vals))
        worst_gal = max(worst_gal, rel)
        print(f"  {sigma0:7.1f} {spread:14.6f} {rel:9.3%}  "
              f"{'PASS' if rel <= args.tol else 'FAIL'}")

    worst = max(max(r["worst_z"], r["worst_x"]) for r in rows)
    print(f"\nworst deviation from analytic : {worst:.3%}")
    print(f"worst k0-dependence           : {worst_gal:.3%}")
    print(f"tolerance                     : {args.tol:.3%}")
    ok = worst <= args.tol and worst_gal <= args.tol
    print("\nOVERALL: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
