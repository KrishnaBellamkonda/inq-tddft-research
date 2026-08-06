#!/usr/bin/env python3
"""Build the single-run notebooks for the CAP-free wrap-around slab KS study.

Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
Skill: .claude/skills/run-notebook/ (run_notebook_builder.py does the work; this
file is the thin per-run driver that knows THIS campaign's parameters).

16 runs: {wp, classical} x {N=100 (r_s 4.18), N=40 (r_s 5.68)} x {2.0, 2.5, 3.0, 3.5}.

WHAT IS PASSED AND WHY
----------------------
--rs           per density, so the analytical Lindhard stopping panel is drawn
               against the RIGHT electron gas (4.1815 vs 5.6751).
--proj-sigma   1.41421 = sigma_WP/sqrt2 for sigma_WP = 2. This is the CHARGE std
               and it only ever appears here and inside the binaries; every label
               stays sigma_WP = 2 (.claude/rules/sigma-wp-convention.md).
--lindhard     'both' — at sigma_pot = 1.41 the projectile is NOT point-like, so
               the finite-sigma Gaussian curve is the meaningful comparison and
               the point-charge curve is the (over-estimating) bound.
--e-gs-ha      the bare-slab GS energy at dx = 0.40, per density, enabling the
               energy-method stopping section.
--no cap-inner THERE IS NO CAP in this study, so no absorbing-band guides are
               drawn. Passing one would draw dashed lines at a fictitious z.
--twin-wp      classical runs only: adds the WP-minus-classical energy-diff bar
               GIF, which is the whole point of having matched twins.

COMPLETENESS GATE
-----------------
A run is built only if run_summary.txt says `run_completed = true`. A
still-propagating run yields plausible-looking but WRONG numbers — that exact
trap bit this project once already (deposit_stopping read 86 of 3623 steps and
returned a confident S; see docs/handovers/wavepacket-highdensity-sv-twin.md).
Incomplete runs are listed and skipped, never silently half-built.

Usage:
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 build_run_notebooks.py [--only wp|classical] [--force]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SCRIPTS = REPO / "ResearchProject/systems/localised_jellium/scripts/slab_ks_wrap"
BUILDER = REPO / ".claude/skills/run-notebook/run_notebook_builder.py"
VENV_PY = REPO / "venv/bin/python3"

# Per-density constants, measured (not assumed) from the dx = 0.40 ground states:
#   scripts/wp_highdensity_sv/gs/results/dx0p4/run_summary.txt      (N = 100)
#   scripts/slab_ks_wrap/gs/results/n40_dx0p4/run_summary.txt       (N =  40)
DENSITY = {
    100: dict(rs=4.1814717081217, e_gs_ha=207.18323030158),
    40:  dict(rs=5.6751302339093, e_gs_ha=31.529527863103),
}
VELOCITIES = (2.0, 2.5, 3.0, 3.5)

SIGMA_POT = 2.0 / (2.0 ** 0.5)   # 1.41421 — sigma_WP = 2
LAUNCH_Z = -24.0
L_SLAB = 25.0
SLAB_HALF = 12.5
GIF_SECONDS = 17.0               # readable loop pace


def run_name(n_elec: int, v: float) -> str:
    return f"n{n_elec}_v{str(v).replace('.', 'p')}"


def is_complete(results_dir: Path) -> bool:
    summary = results_dir / "run_summary.txt"
    if not summary.is_file():
        return False
    return "run_completed = true" in summary.read_text()


def build_one(half: str, n_elec: int, v: float, force: bool) -> tuple[str, str]:
    name = run_name(n_elec, v)
    results = SCRIPTS / half / "results" / name
    out = HERE / f"{half}_{name}.ipynb"

    if not results.is_dir():
        return name, "MISSING (no results dir)"
    if not is_complete(results) and not force:
        return name, "SKIP (run not complete)"

    cfg = DENSITY[n_elec]
    cmd = [
        str(VENV_PY), str(BUILDER), str(results), str(out),
        "--run-cpp", str(SCRIPTS / half / "run.cpp"),
        "--rs", f"{cfg['rs']}",
        "--e-gs-ha", f"{cfg['e_gs_ha']}",
        "--proj-sigma", f"{SIGMA_POT}",
        "--lindhard", "both",
        "--launch-z", f"{LAUNCH_Z}",
        "--v0", f"{v}",
        "--l-slab", f"{L_SLAB}",
        "--slab-half", f"{SLAB_HALF}",
        "--gif-seconds", f"{GIF_SECONDS}",
    ]
    # Matched twin: the classical notebook gets the WP-minus-classical energy GIF.
    if half == "classical":
        twin = SCRIPTS / "wp" / "results" / name
        if twin.is_dir() and is_complete(twin):
            cmd += ["--twin-wp", str(twin)]

    print(f"\n=== {half}/{name} ===\n  {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, cwd=str(HERE))
    return name, ("OK" if proc.returncode == 0 else f"FAILED (rc={proc.returncode})")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=("wp", "classical"), default=None)
    ap.add_argument("--force", action="store_true",
                    help="build even if the run is still propagating (numbers will be partial)")
    a = ap.parse_args()

    halves = [a.only] if a.only else ["wp", "classical"]
    results: list[tuple[str, str, str]] = []
    for half in halves:
        for n_elec in (100, 40):
            for v in VELOCITIES:
                name, status = build_one(half, n_elec, v, a.force)
                results.append((half, name, status))

    print("\n" + "=" * 62)
    print("RUN-NOTEBOOK BUILD SUMMARY")
    print("=" * 62)
    for half, name, status in results:
        print(f"  {half:<10} {name:<12} {status}")
    n_ok = sum(1 for _, _, s in results if s == "OK")
    n_skip = sum(1 for _, _, s in results if s.startswith("SKIP") or s.startswith("MISSING"))
    n_bad = len(results) - n_ok - n_skip
    print(f"\n  built {n_ok} | skipped {n_skip} | failed {n_bad}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
