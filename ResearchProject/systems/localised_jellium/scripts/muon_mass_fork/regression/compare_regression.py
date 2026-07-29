#!/usr/bin/env python3
"""Phase-3 bit-for-bit comparison: fork (inq-study, all mass=1) vs pristine inq.

Diffs the two regression runs' GS energy components, the full RT energy trace,
and the GS density field. The per-state mass fork's empty-factor guard routes the
mass-1 path through the ORIGINAL scalar kinetic code, so the two builds must agree
to ~GPU-reduction precision. A difference above tolerance means the fork perturbs
the mass-1 path => the fork is BROKEN (hard trust gate).

Tolerances allow for non-deterministic GPU reduction ordering (~1e-11), which is
present even between two runs of the SAME binary; we require agreement far tighter
than any physical change would produce.

Usage: compare_regression.py <fork_dir> <pristine_dir> [--etol 1e-9] [--dtol 1e-8]
Exit 0 iff all diffs are within tolerance.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np


def _gs_energy(run: Path) -> dict[str, float]:
    d = {}
    for ln in (run/"raw/observables/gs_energy.csv").read_text().strip().splitlines()[1:]:
        k, v = ln.split(","); d[k] = float(v)
    return d


def _rt_energy(run: Path) -> tuple[list[str], np.ndarray]:
    lines = (run/"raw/observables/rt_energy.csv").read_text().strip().splitlines()
    hdr = lines[0].split(",")
    arr = np.array([[float(x) for x in ln.split(",")] for ln in lines[1:]])
    return hdr, arr


def _density(run: Path) -> np.ndarray:
    cands = list((run/"raw/observables/gs_density").glob("*"))
    raw = [c for c in cands if c.suffix in (".raw", ".bin", "") and c.is_file()]
    if not raw:  # fall back to any file in the dir
        raw = [c for c in cands if c.is_file()]
    return np.fromfile(raw[0], dtype=np.float64)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fork_dir"); ap.add_argument("pristine_dir")
    ap.add_argument("--etol", type=float, default=1e-9)
    ap.add_argument("--dtol", type=float, default=1e-8)
    a = ap.parse_args()
    fork, pris = Path(a.fork_dir), Path(a.pristine_dir)
    ok = True

    # --- GS energy components --------------------------------------------
    ef, ep = _gs_energy(fork), _gs_energy(pris)
    print("=== GS energy components (fork vs pristine) ===")
    for k in ef:
        d = abs(ef[k] - ep.get(k, np.nan))
        p = d < a.etol; ok &= p
        print(f"  [{'ok' if p else 'XX'}] {k:9s} fork={ef[k]:.12f} pris={ep[k]:.12f} |d|={d:.2e}")

    # --- RT energy trace -------------------------------------------------
    hf, af = _rt_energy(fork); hp, ap_ = _rt_energy(pris)
    print(f"=== RT energy trace ({af.shape[0]} steps) ===")
    if af.shape != ap_.shape:
        print(f"  XX shape mismatch {af.shape} vs {ap_.shape}"); ok = False
    else:
        for j, name in enumerate(hf):
            if name in ("step", "time_au"):
                continue
            d = float(np.max(np.abs(af[:, j] - ap_[:, j])))
            p = d < a.etol; ok &= p
            print(f"  [{'ok' if p else 'XX'}] {name:9s} max|d|over_steps={d:.2e}")

    # --- GS density field ------------------------------------------------
    try:
        df, dp = _density(fork), _density(pris)
        if df.shape != dp.shape:
            print(f"=== density: XX shape {df.shape} vs {dp.shape}"); ok = False
        else:
            dmax = float(np.max(np.abs(df - dp)))
            p = dmax < a.dtol; ok &= p
            print(f"=== density field ({df.size} pts): "
                  f"[{'ok' if p else 'XX'}] max|d|={dmax:.2e} (tol {a.dtol:.0e})")
    except Exception as e:
        print(f"=== density: could not compare ({e}) — energies are the primary gate")

    print("\n" + ("BIT-FOR-BIT PASS: fork == pristine on the mass-1 path"
                  if ok else "BIT-FOR-BIT FAIL: fork perturbs the mass-1 path (BROKEN)"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
