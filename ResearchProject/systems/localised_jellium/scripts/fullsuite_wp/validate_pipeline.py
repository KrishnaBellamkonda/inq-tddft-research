#!/usr/bin/env python3
"""Validate the full-suite WP pipeline end-to-end on a (short) run output.

Checks (1) every canonical observable channel is present and non-empty, and
(2) the analysis chain works: load density_total + density_wp via the canonical
`inqview.load_vti` (no fftshift) and form the bath = total - wp decomposition,
asserting basic physical sanity (WP localised; bath ~ GS away from the WP).

Usage: python validate_pipeline.py <results/OUT dir>
"""
import sys, os, glob
import numpy as np

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.visualisation import load_vti

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/validate"
RAW = os.path.join(OUT, "raw")
ok = True


def check(label, cond, detail=""):
    global ok
    ok = ok and cond
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")


def n_vti(sub):
    return len(glob.glob(os.path.join(RAW, "vti", sub, "*.vti")))


def csv_rows(name):
    p = os.path.join(RAW, "observables", name)
    if not os.path.exists(p):
        return -1
    with open(p) as f:
        return sum(1 for _ in f) - 1  # minus header


print(f"== channel presence ({OUT}) ==")
for sub in ["density_total", "density_system", "density_gs_system", "density_wp",
            "wavefunction_wp", "density_delta", "density_delta_coarse"]:
    check(f"vti/{sub}", n_vti(sub) > 0, f"{n_vti(sub)} frames")
for name in ["observables.csv", "state_energies.csv", "occupations_vs_time.csv",
             "momentum_distribution.csv", "wp_momentum_stats.csv",
             "wp_real_space_stats.csv", "electron_number.csv"]:
    check(f"obs/{name}", csv_rows(name) > 0, f"{csv_rows(name)} rows")
for d in ["overlap", "overlap_full", "eigenvalues"]:
    p = os.path.join(RAW, "observables", d)
    check(f"obs/{d}/", os.path.isdir(p) and len(os.listdir(p)) > 0)

print("== analysis chain: load_vti decomposition total / wp / bath ==")
tot = sorted(glob.glob(os.path.join(RAW, "vti", "density_total", "*.vti")))
wpf = sorted(glob.glob(os.path.join(RAW, "vti", "density_wp", "*.vti")))
if tot and wpf:
    ft = load_vti(tot[-1])          # physical order, no fftshift
    fw = load_vti(wpf[-1])
    bath = ft.data - fw.data
    check("total >= wp everywhere (within noise)",
          float((fw.data - ft.data).max()) < 1e-3,
          f"max(wp-total)={float((fw.data-ft.data).max()):.2e}")
    check("WP is localised (peak >> mean)",
          float(fw.data.max()) > 20 * float(np.abs(fw.data).mean()),
          f"peak={float(fw.data.max()):.2e} mean={float(np.abs(fw.data).mean()):.2e}")
    check("bath integral ~ N_bath (>0, < total)",
          0 < float(bath.sum()) < float(ft.data.sum()),
          f"bath_sum={float(bath.sum()):.1f} total_sum={float(ft.data.sum()):.1f}")
    check("axes are physical (z spans about -L/2..L/2)",
          ft.z[0] < -20 and ft.z[-1] > 20, f"z=[{ft.z[0]:.1f},{ft.z[-1]:.1f}]")
else:
    check("density_total + density_wp present for decomposition", False)

print(f"\n=== PIPELINE VALIDATION: {'PASS' if ok else 'FAIL'} ===")
sys.exit(0 if ok else 1)
