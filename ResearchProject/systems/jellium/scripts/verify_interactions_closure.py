#!/usr/bin/env python3
"""verify_interactions_closure.py — closure gates for interactions.csv.

Mandated by .claude/rules/decomposed-interaction-energies.md: the pairwise
P/S/B terms are only trustworthy if they sum back to the INQ scalars. This
checks that on a finished run and exits non-zero if they do not.

    python3 verify_interactions_closure.py <run_dir> <wp|classical>

WHAT IS AND IS NOT A REAL GATE
------------------------------
For the WP half, E_SS + E_PS + E_PP == e_hartree_check is an ALGEBRAIC
IDENTITY of compute_coulomb_wp (E_SS is *defined* as e_hartree_check - cross
+ E_PP), so it is exact by construction and proves nothing on its own. It is
still checked, because a violation means the CSV is corrupt rather than that
the physics is wrong.

The gate that carries information is the cross-check against INQ's own
scalars in observables.csv:

    WP        : e_hartree_check  == energy_hartree   (INQ)
    classical : e_ss             == energy_hartree   (INQ, written inline
                                                      as e_hartree_inq)

BULK NOTE: the background is uniform, so phi_plus == 0 identically and
E_SB = E_PB = E_BB = 0. Those columns are asserted to BE zero here; that is a
schema check, not physics. In a slab run they would be non-trivial and this
script's background gates would need extending.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Ha. INQ's Hartree energy is an O(1-100 Ha) quantity built from the same grid
# and Poisson solver, so agreement should be near machine precision. 1e-9 leaves
# room for float64 CSV round-tripping at 12 significant digits.
TOL_CLOSURE = 1e-9
TOL_ZERO = 1e-12


def _load_segments(obs_dir: Path, stem: str) -> pd.DataFrame:
    """Concatenate stem.csv + stem.fromNNN.csv in step order (resume segments)."""
    files = sorted(obs_dir.glob(f"{stem}*.csv"))
    if not files:
        raise FileNotFoundError(f"no {stem}*.csv in {obs_dir}")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return df.sort_values("step").drop_duplicates("step", keep="last").reset_index(drop=True)


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    run_dir, half = Path(sys.argv[1]), sys.argv[2]
    if half not in ("wp", "classical"):
        print(f"FATAL: half must be wp|classical (got {half!r})")
        return 2

    obs_dir = run_dir / "results" / "raw" / "observables"
    try:
        ix = _load_segments(obs_dir, "interactions")
    except FileNotFoundError as e:
        print(f"FATAL: {e}")
        return 3

    print(f"--- interactions closure: {run_dir.name}/{half} ---")
    print(f"  rows: {len(ix)}   steps {ix['step'].min()}..{ix['step'].max()}")

    fails = 0

    def gate(name, resid, tol):
        nonlocal fails
        m = float(np.nanmax(np.abs(resid)))
        ok = m <= tol
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: max|resid| = {m:.3e} (tol {tol:.0e})")
        if not ok:
            fails += 1

    # --- background terms are identically zero in BULK ---------------------
    for col in ("e_sb", "e_pb", "e_bb"):
        if col in ix:
            gate(f"{col} == 0 (bulk, phi+ = 0)", ix[col].to_numpy(), TOL_ZERO)

    if half == "wp":
        # Internal identity — corruption check, not a physics gate.
        gate("E_SS + E_PS + E_PP == e_hartree_check (identity)",
             ix["e_ss"] + ix["e_ps"] + ix["e_pp"] - ix["e_hartree_check"],
             TOL_CLOSURE)
        # The real gate: our Poisson sum vs INQ's own Hartree energy.
        try:
            obs = _load_segments(obs_dir, "observables")
            m = ix.merge(obs[["step", "energy_hartree"]], on="step", how="inner")
            if m.empty:
                print("  [WARN] no overlapping steps with observables.csv "
                      "— INQ cross-check skipped")
            else:
                print(f"  (INQ cross-check on {len(m)} shared steps)")
                gate("e_hartree_check == INQ energy_hartree",
                     m["e_hartree_check"] - m["energy_hartree"], TOL_CLOSURE)
        except FileNotFoundError as e:
            print(f"  [WARN] {e} — INQ cross-check skipped")
        norm = ix["norm_wp"].to_numpy()
        print(f"  [info] norm_wp: {norm.min():.6f} .. {norm.max():.6f} "
              f"(no CAP in this study, so it must stay ~1)")
        if norm.max() < 0.5:
            print("  [FAIL] norm_wp collapsed — E_PP is not measuring what it claims")
            fails += 1
    else:
        # Classical projectile is NOT in n, so E_SS is the whole Hartree energy.
        gate("e_ss == INQ energy_hartree", ix["e_ss"] - ix["e_hartree_inq"],
             TOL_CLOSURE)
        # E_PP of a rigid Gaussian cloud is a CONSTANT of the motion -- but only
        # while the cloud is FULLY INSIDE the box. Once the projectile nears the
        # +z face its Gaussian tail is clipped off the grid, charge is lost
        # (norm_proj < 1) and E_PP falls. Verified on the sigma=3 r_s=5.702 run
        # (2026-08-01): E_PP is constant to 2.6e-11 Ha for norm_proj == 1, then
        # decays over the final 32 of 301 rows as norm_proj -> 0.9942. That is
        # clipping, NOT egg-box error (correlation with sub-grid phase is 0.03).
        #
        # So gate constancy on the CLEAN rows and REPORT the clipping onset --
        # it is the hard upper bound on any fit window using this run.
        epp = ix["e_pp"].to_numpy()
        npj = ix["norm_proj"].to_numpy()
        # Threshold measured, not guessed (sigma=3 r_s=5.702, 2026-08-01):
        #   norm_proj == 1.0 exactly -> E_PP spread is BIT-EXACTLY ZERO (237 rows)
        #   norm_proj >= 1-1e-9      -> 1.8e-11 Ha (260 rows, to t = 20.96)
        #   norm_proj >= 1-1e-6      -> 4.5e-8 Ha  (admits the clipping shoulder)
        # 1e-9 keeps the full fit window (ends 19.48) while excluding the shoulder.
        clean = npj >= 1.0 - 1e-9
        if clean.sum() < 10:
            print(f"  [FAIL] only {clean.sum()} rows with an unclipped projectile "
                  f"cloud -- the run is too short or launched too near a face")
            fails += 1
        else:
            spread = float(np.nanmax(epp[clean]) - np.nanmin(epp[clean]))
            ok = spread <= TOL_CLOSURE
            print(f"  [{'PASS' if ok else 'FAIL'}] E_PP constant on unclipped rows "
                  f"({clean.sum()}/{len(epp)}): spread = {spread:.3e} Ha, "
                  f"value = {np.nanmean(epp[clean]):.9f} Ha "
                  f"({np.nanmean(epp[clean]) * 27.211:.4f} eV)")
            if not ok:
                fails += 1
        # Report the START OF THE TRAILING CLIPPED RUN, not the first row below
        # threshold: the first few rows sit a few 1e-9 under 1.0 purely from
        # discretising the Gaussian on the grid (0.999999996 at launch), which is
        # not clipping. Only the contiguous tail is the projectile leaving the box.
        tail = len(clean)
        while tail > 0 and not clean[tail - 1]:
            tail -= 1
        if tail < len(clean):
            tv = ix["time_au"].to_numpy()
            zv = ix["proj_z"].to_numpy()
            print(f"  [info] projectile cloud starts clipping the box face at "
                  f"t = {tv[tail]:.2f} a.u. (z = {zv[tail]:+.2f}, norm_proj -> "
                  f"{np.nanmin(npj):.6f}) -- ANY FIT MUST END BEFORE THIS")

    print(f"--- {'ALL CLOSURE GATES PASSED' if fails == 0 else f'{fails} GATE(S) FAILED'} ---")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
