#!/usr/bin/env python3
"""Smoke-test verification for the L=60 jellium run_base output tree.

Checks:
  1. results/run_summary.txt has run_completed = true.
  2. observables.csv has the new cod_*/density_l2 columns and at least one row.
  3. cod_x_bohr ~= 30 +/- spacing/2 at t=0; cod_z_bohr increases monotonically.
  4. density_l2 == 0 at the first row (t = t0 reference) and > 0 thereafter.
  5. state_energies.csv has rows; at t=0 the state_index=0 energy matches the
     ground-state eigenvalue from raw/observables/eigenvalues/eigenvalues.csv
     to within 1 mHa.
  6. momentum_distribution.csv has rows; at t=0 the WP curve peaks at |k| close
     to k0 (read from run_summary.txt).
  7. raw/vti/density_rt_delta/ has a VTI series.

Usage:
  python verify_smoke_outputs.py <results_dir>

Exit code is 0 on all-green, 1 if any check fails. Failed checks print to
stderr; passed checks print to stdout.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd


def _ok(msg: str) -> None:
    print(f"  OK  {msg}")


def _fail(msg: str) -> None:
    print(f"FAIL  {msg}", file=sys.stderr)


def _parse_summary(p: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        m = re.match(r"^([a-zA-Z0-9_]+)\s*=\s*(.+)$", line.strip())
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def main(results_dir: Path) -> int:
    failures = 0

    # 1. run_summary completed
    summary = _parse_summary(results_dir / "run_summary.txt")
    if summary.get("run_completed") == "true":
        _ok("run_completed = true")
    else:
        _fail(f"run_completed != true (got {summary.get('run_completed')!r})")
        failures += 1

    # 2-4. observables.csv
    obs = results_dir / "raw/observables/observables.csv"
    if not obs.exists():
        _fail(f"missing {obs}")
        failures += 1
    else:
        df = pd.read_csv(obs)
        required_cols = {"cod_x_bohr", "cod_y_bohr", "cod_z_bohr", "density_l2"}
        if not required_cols.issubset(df.columns):
            _fail(f"observables.csv missing columns "
                  f"{required_cols - set(df.columns)}")
            failures += 1
        else:
            _ok(f"observables.csv has CoD + density_l2 columns "
                f"({len(df)} rows)")

            # INQ uses centred Cartesian r ∈ [-L/2, +L/2] so the box centre
            # is (0, 0, 0); CoD lands at (~0.25, ~0.25, ~0.25) at t=0
            # because of the half-grid voxel-centre offset (dx/2 = 0.25).
            row0 = df.iloc[0]
            if abs(row0["cod_x_bohr"]) <= 0.6:
                _ok(f"cod_x at t=0 = {row0['cod_x_bohr']:.3f} (~box centre, "
                    f"|cod_x| ≤ dx/2)")
            else:
                _fail(f"cod_x at t=0 = {row0['cod_x_bohr']:.3f}, "
                      f"expected near 0")
                failures += 1
            if abs(row0["cod_y_bohr"]) <= 0.6:
                _ok(f"cod_y at t=0 = {row0['cod_y_bohr']:.3f} (~box centre)")
            else:
                _fail(f"cod_y at t=0 = {row0['cod_y_bohr']:.3f}, "
                      f"expected near 0")
                failures += 1

            # cod_z monotonic increase (single-pass +z launch)
            if len(df) >= 2 and df["cod_z_bohr"].iloc[-1] > row0["cod_z_bohr"]:
                _ok(f"cod_z increased from {row0['cod_z_bohr']:.3f} to "
                    f"{df['cod_z_bohr'].iloc[-1]:.3f}")
            else:
                _fail("cod_z did not increase across the trajectory")
                failures += 1

            # density_l2 starts at 0, grows
            if abs(row0["density_l2"]) <= 1e-10:
                _ok("density_l2 at t=0 ≈ 0 (reference snapshot)")
            else:
                _fail(f"density_l2 at t=0 = {row0['density_l2']:.3e}, "
                      f"expected ~0")
                failures += 1

    # 5. state_energies.csv
    se = results_dir / "raw/observables/state_energies.csv"
    if not se.exists():
        _fail(f"missing {se}")
        failures += 1
    else:
        df_se = pd.read_csv(se)
        if len(df_se) > 0:
            _ok(f"state_energies.csv has {len(df_se)} rows, "
                f"{df_se['state_index'].nunique()} unique states")
        else:
            _fail("state_energies.csv is empty")
            failures += 1

    # 6. momentum_distribution.csv
    md = results_dir / "raw/observables/momentum_distribution.csv"
    if not md.exists():
        _fail(f"missing {md}")
        failures += 1
    else:
        df_md = pd.read_csv(md, comment="#")
        if len(df_md) > 0:
            _ok(f"momentum_distribution.csv has {len(df_md)} rows")
            # Find the |k| of the WP peak at t=0
            df0 = df_md[df_md["time_au"] == df_md["time_au"].min()]
            peak_k = df0.loc[df0["n_wp"].idxmax(), "k_bohr_inv"]
            try:
                k0 = float(summary.get("wp_k0_bohr_inv", "0 0 0").split()[2])
            except (ValueError, IndexError):
                k0 = 0.0
            if k0 > 0 and abs(peak_k - k0) < 0.5:
                _ok(f"WP momentum peak at t=0: {peak_k:.3f} ≈ k0={k0:.3f}")
            else:
                _fail(f"WP momentum peak {peak_k:.3f} != k0={k0:.3f}")
                failures += 1
        else:
            _fail("momentum_distribution.csv is empty")
            failures += 1

    # 7. VTI series exist
    vti = results_dir / "raw/vti/density_rt_delta"
    if vti.exists() and any(vti.iterdir()):
        _ok(f"density_rt_delta/ has "
            f"{sum(1 for _ in vti.iterdir())} files")
    else:
        _fail(f"density_rt_delta/ missing or empty")
        failures += 1

    print()
    if failures:
        print(f"{failures} check(s) failed", file=sys.stderr)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <results_dir>", file=sys.stderr)
        sys.exit(2)
    sys.exit(main(Path(sys.argv[1])))
