#!/usr/bin/env python3
"""summarise_epp_across_pairs.py — cross-pair table of the gauge-clean dE_PP.

    venv/bin/python summarise_epp_across_pairs.py

For every bulk twin pair that carries interactions.csv on BOTH halves, report

    dE_PP = E_PP(WP) - E_PP(classical)

evaluated in the SAME cell. This is the only self-Hartree quantity that is
gauge-clean: absolute E_PP carries the charged-cell G=0 gauge, so it cannot be
compared across geometries (that error invalidated a slab-vs-bulk comparison on
2026-08-01). Within one cell the gauge cancels in the difference.

WHY THIS TABLE. E_PP is the projectile's self-Hartree -- the uncancelled
self-interaction a wavepacket has in LDA and a classical point projectile does
not. It is the leading candidate for the unexplained ~2.2x residual in the
classical/WP stopping ratio. With sigma in {2,3} crossed with r_s in
{3.987, 5.702} this becomes a 2x2 design that separates a WIDTH dependence from
a DENSITY dependence -- which no single pair can do.

WINDOWING. Every quantity is evaluated on [FIT_T0, t_max] where

    t_max = min(FIT_T1 from the config, projectile-cloud clipping onset)

The clipping bound is NOT in the config: the classical projectile's Gaussian
cloud is clipped by the +z box face near the end of a run, norm_proj falls below
1 and E_PP decays. Measured on the sigma=3 r_s=5.702 run (2026-08-01): E_PP is
bit-exactly constant while norm_proj == 1, then decays from t = 21.04 a.u.
See docs/plans/bulk-jellium-ks-stopping.md section 10.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SYSTEM = Path(__file__).resolve().parents[2]
SCRIPTS = SYSTEM / "scripts"
HA_TO_EV = 27.211

# label, FIT_T0, FIT_T1 (from the config headers; sigma=1 excluded by user)
PAIRS = {
    "bulk_ks_stopping":            dict(label="r_s 5.702, sigma 2", t0=4.0, t1=18.43),
    "bulk_ks_stopping_sigma3":     dict(label="r_s 5.702, sigma 3", t0=4.0, t1=19.48),
    "bulk_ks_stopping_rs4":        dict(label="r_s 3.987, sigma 2", t0=4.0, t1=18.43),
    "bulk_ks_stopping_rs4_sigma3": dict(label="r_s 3.987, sigma 3", t0=4.0, t1=19.48),
}


def _load(variant: str, half: str):
    d = SCRIPTS / variant / half / "results/raw/observables"
    files = sorted(d.glob("interactions*.csv"))
    if not files:
        return None
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return df.sort_values("step").drop_duplicates("step", keep="last").reset_index(drop=True)


def _clip_onset(cl: pd.DataFrame):
    """Time at which the classical projectile cloud starts leaving the grid."""
    if "norm_proj" not in cl:
        return np.inf
    clean = cl["norm_proj"].to_numpy() >= 1.0 - 1e-9
    tail = len(clean)
    while tail > 0 and not clean[tail - 1]:
        tail -= 1
    return np.inf if tail >= len(clean) else float(cl["time_au"].to_numpy()[tail])


def main() -> int:
    rows = []
    print("=" * 96)
    print("Gauge-clean self-Hartree collapse  dE_PP = E_PP(WP) - E_PP(classical), same cell")
    print("=" * 96)
    for var, m in PAIRS.items():
        wp, cl = _load(var, "wp"), _load(var, "classical")
        if wp is None or cl is None:
            missing = [h for h, d in (("wp", wp), ("classical", cl)) if d is None]
            print(f"  {m['label']:<22} SKIP - interactions.csv missing on: {', '.join(missing)}")
            continue

        onset = _clip_onset(cl)
        t_max = min(m["t1"], onset)
        bound = "config" if t_max == m["t1"] else "CLIPPING"

        j = wp.merge(cl, on="step", suffixes=("_wp", "_cl"))
        if j.empty:
            print(f"  {m['label']:<22} SKIP - no shared steps between halves")
            continue
        t = j["time_au_wp"].to_numpy()
        w = (t >= m["t0"]) & (t <= t_max)
        if w.sum() < 5:
            print(f"  {m['label']:<22} SKIP - only {w.sum()} rows in window")
            continue

        epp_wp, epp_cl = j["e_pp_wp"].to_numpy(), j["e_pp_cl"].to_numpy()
        # t=0 reference: the packet and its matched classical cloud must start equal
        # (sigma_pot = sigma_WP/sqrt2). A mismatch here means the UPF and the WP
        # disagree and NOTHING downstream is comparable.
        i0 = int(np.argmin(np.abs(t)))
        d0 = (epp_wp[i0] - epp_cl[i0]) * HA_TO_EV

        # dE_PP at the END of the clean window = how much self-Hartree the packet
        # has shed relative to the rigid classical cloud.
        d_end = (epp_wp[w][-1] - epp_cl[w][-1]) * HA_TO_EV
        rows.append(dict(label=m["label"], n=int(w.sum()), t_max=t_max, bound=bound,
                         epp0=epp_wp[i0] * HA_TO_EV, d0=d0, d_end=d_end,
                         drop=(epp_wp[w][-1] - epp_wp[i0]) * HA_TO_EV,
                         cl_drift=(epp_cl[w][-1] - epp_cl[i0]) * HA_TO_EV))

    if not rows:
        print("\nNo pair has interactions.csv on both halves yet.")
        return 1

    print(f"\n{'pair':<22} {'n':>4} {'t_max':>7} {'bound':>9} "
          f"{'E_PP(0)':>9} {'dE_PP(0)':>9} {'dE_PP_end':>10} {'WP drop':>9} {'cl drift':>9}")
    print(f"{'':22} {'':>4} {'a.u.':>7} {'':>9} {'eV':>9} {'eV':>9} {'eV':>10} {'eV':>9} {'eV':>9}")
    print("-" * 96)
    for r in rows:
        print(f"{r['label']:<22} {r['n']:>4} {r['t_max']:>7.2f} {r['bound']:>9} "
              f"{r['epp0']:>9.4f} {r['d0']:>9.2e} {r['d_end']:>10.4f} "
              f"{r['drop']:>9.4f} {r['cl_drift']:>9.2e}")

    print("\nGATES")
    bad = [r for r in rows if abs(r["d0"]) > 1e-3]
    print(f"  dE_PP(0) == 0 (sigma matching, sigma_pot = sigma_WP/sqrt2): "
          f"{'PASS' if not bad else 'FAIL on ' + ', '.join(r['label'] for r in bad)}")
    bad2 = [r for r in rows if abs(r["cl_drift"]) > 1e-4]
    print(f"  classical E_PP constant in window (rigid cloud): "
          f"{'PASS' if not bad2 else 'FAIL on ' + ', '.join(r['label'] for r in bad2)}")

    if len(rows) == 4:
        by = {r["label"]: r["d_end"] for r in rows}
        print("\n2x2 DESIGN — does self-Hartree collapse track WIDTH or DENSITY?")
        try:
            w_lo = by["r_s 5.702, sigma 3"] - by["r_s 5.702, sigma 2"]
            w_hi = by["r_s 3.987, sigma 3"] - by["r_s 3.987, sigma 2"]
            d_s2 = by["r_s 3.987, sigma 2"] - by["r_s 5.702, sigma 2"]
            d_s3 = by["r_s 3.987, sigma 3"] - by["r_s 5.702, sigma 3"]
            print(f"  width effect  (sigma 3 - sigma 2):  {w_lo:+.4f} eV at r_s 5.702, "
                  f"{w_hi:+.4f} eV at r_s 3.987")
            print(f"  density effect (r_s 3.987 - 5.702): {d_s2:+.4f} eV at sigma 2, "
                  f"{d_s3:+.4f} eV at sigma 3")
            print(f"  interaction (non-additivity):       "
                  f"{(w_hi - w_lo):+.4f} eV")
        except KeyError:
            pass

    print("\nNOTE: dE_PP is gauge-clean ONLY within a cell. Do not compare these")
    print("      against a slab run's absolute E_PP (different G=0 gauge).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
