#!/usr/bin/env python3
"""validate_wp.py <GS_CKPT> — programmatic WP-start-position check (replaces the
manual eyeball for the autonomous chain). Reads the GS z-density profile and
asserts the WP launch site sits in vacuum, clear of the slab.

Checks (exit 0 = pass, 1 = fail):
  * launch z = -20.5 is 8 Bohr from the -12.5 slab face (config invariant);
  * GS electron density at the launch z is < 2% of the slab-centre plateau
    (i.e. the WP starts in genuine vacuum, not overlapping the slab tail);
  * the slab half-density edge is near |z| ~ 12.5 (confirms the slab is where
    the config says, so "8 Bohr from the face" is physically true).
"""
import sys
from pathlib import Path
import numpy as np

LAUNCH_Z = -20.5
SLAB_FACE = 12.5

def main() -> int:
    ck = Path(sys.argv[1])
    prof = ck / "z_density_profile.csv"
    if not prof.exists():
        print(f"FAIL: no {prof}"); return 1
    import csv
    z, ne = [], []
    with open(prof) as f:
        r = csv.DictReader(f)
        for row in r:
            z.append(float(row["z_bohr"])); ne.append(float(row["n_e_planar_avg"]))
    z = np.array(z); ne = np.array(ne)

    n_center = ne[np.argmin(np.abs(z))]
    n_launch = ne[np.argmin(np.abs(z - LAUNCH_Z))]
    frac = n_launch / n_center if n_center > 0 else 1.0

    # slab half-density edge on the -z side
    half = 0.5 * n_center
    neg = z < 0
    edge = None
    zs, ns = z[neg], ne[neg]
    order = np.argsort(zs)
    zs, ns = zs[order], ns[order]
    for i in range(1, len(zs)):
        if ns[i-1] < half <= ns[i]:
            edge = zs[i-1] + (zs[i]-zs[i-1])*(half-ns[i-1])/(ns[i]-ns[i-1]); break

    print(f"n_center={n_center:.4e}  n(launch z={LAUNCH_Z})={n_launch:.4e}  "
          f"frac={frac:.3%}  -z half-density edge≈{edge}")
    ok = True
    if frac > 0.02:
        print(f"FAIL: WP launch density {frac:.2%} of centre (>2%) — overlaps slab"); ok = False
    if edge is None or abs(abs(edge) - SLAB_FACE) > 2.0:
        print(f"FAIL: slab edge {edge} not near ∓{SLAB_FACE} (±2 Bohr)"); ok = False
    gap = abs(LAUNCH_Z) - SLAB_FACE
    print(f"launch-to-face gap = {gap:.2f} Bohr (expect 8.0)")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
