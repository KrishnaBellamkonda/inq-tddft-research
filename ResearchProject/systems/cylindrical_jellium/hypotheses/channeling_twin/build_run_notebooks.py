#!/usr/bin/env python3
"""Build the two single-run notebooks of the annular-tube channeling twin.

Plan:  docs/plans/cylindrical-channeling-ks-stopping.md
Skill: .claude/skills/run-notebook/ (run_notebook_builder.py does the work; this
       file is the thin driver that knows THIS study's parameters).

Two runs: classical (Gaussian-charge perturbation) and wp (KS orbital), matched
in every physical parameter. The cross-run comparison lives in the SEPARATE
comparison notebook (build_comparison_notebook.py) — these two are the per-run
deep dives, each opening with the mandatory density-GIF battery
(.claude/rules/notebook-density-gif.md).

WHAT IS PASSED AND WHY
----------------------
--rs 3.000     the r_s of the tube WALL, so the analytical Lindhard panel is
               drawn against the right electron gas. NOTE this is a genuine
               approximation for a channeling geometry: Lindhard is a BULK
               response and the projectile here is in a vacuum bore, so the curve
               is an upper reference, not a prediction. Said again in the
               notebook narrative.
--proj-sigma   2.82843 = sigma_WP/sqrt2 for sigma_WP = 4. This is the CHARGE std
               and appears ONLY here and inside the binaries; every label stays
               sigma_WP = 4 (.claude/rules/sigma-wp-convention.md).
--lindhard both  at sigma_pot = 2.83 the projectile is very far from point-like,
               so the finite-sigma Gaussian curve is the meaningful comparison
               and the point-charge curve is the (badly over-estimating) bound.
--slab-half 0  THERE IS NO SLAB. The tube is uniform along z, so the medium fills
               every z the projectile visits and there are no faces to mark.
               Passing 12.5 (the default) would draw dashed lines at a fictitious
               boundary. 0 is the same switch bulk runs use.
--l-slab 60    = L_z: the medium thickness traversed in one pass, for the
               deltaE/L energy-method estimate.
--e-gs-ha      the bare-tube ground-state energy, READ FROM THE GS RUN SUMMARY
               rather than hard-coded, so it can never drift out of sync with the
               checkpoint the runs actually loaded.
--twin-wp      classical only: adds the WP-minus-classical energy-diff bar GIF,
               which is the whole point of having matched twins.
no --cap-inner THERE IS NO CAP in this study.

COMPLETENESS GATE
-----------------
A run is built only if run_summary.txt says `run_completed = true`. A
still-propagating run yields plausible-looking but WRONG numbers — that exact
trap bit this project once already (deposit_stopping read 86 of 3623 steps and
returned a confident S; docs/handovers/wavepacket-highdensity-sv-twin.md).
Incomplete runs are listed and skipped, never silently half-built.

Usage:
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 build_run_notebooks.py \
        [--only wp|classical] [--force]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SCRIPTS = REPO / "ResearchProject/systems/cylindrical_jellium/scripts/channeling_twin"
GS_SUMMARY = SCRIPTS / "gs/results/gs/run_summary.txt"
BUILDER = REPO / ".claude/skills/run-notebook/run_notebook_builder.py"
VENV_PY = REPO / "venv/bin/python3"

sys.path.insert(0, str(HERE))
import channeling_stopping as CS  # noqa: E402

RUN_NAME = {"wp": "wp", "classical": "classical"}
GIF_SECONDS = 17.0               # readable loop pace


def gs_energy_ha() -> float | None:
    """E_GS of the bare tube, from the ground state's own summary.

    Returned as None if the ground state has not been run, in which case the
    energy-method stopping section of the notebook simply does not appear —
    better than passing a plausible-looking invented number.
    """
    if not GS_SUMMARY.is_file():
        return None
    for line in GS_SUMMARY.read_text().splitlines():
        if line.strip().startswith("ground_state_energy_ha"):
            try:
                return float(line.split("=", 1)[1])
            except (IndexError, ValueError):
                return None
    return None


def is_complete(results_dir: Path) -> bool:
    summary = results_dir / "run_summary.txt"
    return summary.is_file() and "run_completed = true" in summary.read_text()


def build_one(half: str, force: bool, e_gs: float | None) -> tuple[str, str]:
    name = RUN_NAME[half]
    results = SCRIPTS / half / "results" / name
    out = HERE / f"{half}_{name}.ipynb"

    if not results.is_dir():
        return name, "MISSING (no results dir)"
    if not is_complete(results) and not force:
        return name, "SKIP (run not complete)"

    cmd = [
        str(VENV_PY), str(BUILDER), str(results), str(out),
        "--run-cpp", str(SCRIPTS / half / "run.cpp"),
        "--rs", f"{CS.RS:.6f}",
        "--proj-sigma", f"{CS.SIGMA_POT:.6f}",
        "--lindhard", "both",
        "--launch-z", f"{CS.LAUNCH_Z}",
        "--v0", f"{CS.V0}",
        "--l-slab", f"{CS.LZ}",
        "--slab-half", "0",
        "--gif-seconds", f"{GIF_SECONDS}",
    ]
    if e_gs is not None:
        cmd += ["--e-gs-ha", f"{e_gs}"]

    # Matched twin: the classical notebook gets the WP-minus-classical energy GIF.
    if half == "classical":
        twin = SCRIPTS / "wp" / "results" / RUN_NAME["wp"]
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

    e_gs = gs_energy_ha()
    if e_gs is None:
        print(f"NOTE: no ground-state energy at {GS_SUMMARY} — the energy-method "
              f"stopping section will be omitted rather than fabricated.")
    else:
        print(f"bare-tube E_GS = {e_gs:.6f} Ha (from the ground state's own summary)")

    halves = [a.only] if a.only else ["wp", "classical"]
    results: list[tuple[str, str, str]] = []
    for half in halves:
        name, status = build_one(half, a.force, e_gs)
        results.append((half, name, status))

    print("\n" + "=" * 62)
    print("RUN-NOTEBOOK BUILD SUMMARY")
    print("=" * 62)
    for half, name, status in results:
        print(f"  {half:<10} {name:<12} {status}")
    n_ok = sum(1 for _, _, s in results if s == "OK")
    n_skip = sum(1 for _, _, s in results if s.startswith(("SKIP", "MISSING")))
    n_bad = len(results) - n_ok - n_skip
    print(f"\n  built {n_ok} | skipped {n_skip} | failed {n_bad}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
