#!/usr/bin/env python3
"""Single-run notebook for the SIC-corrected channeling WP run.

Plan: docs/plans/wp-self-interaction-correction.md (reviewed 2026-08-02).

IDENTICAL to the channeling_twin per-run WP notebook (same skill builder, same
arguments — see ../channeling_twin/build_run_notebooks.py for the rationale of
every flag), pointed at scripts/channeling_sic/wp/results/wp_sic. The
cross-run three-way story (classical / wp / wp+SIC) lives in the refined
notebook; this is the per-run deep dive with the mandatory density-GIF battery
(.claude/rules/notebook-density-gif.md).

Usage:
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 build_run_notebooks.py [--force]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SIC_SCRIPTS = REPO / "ResearchProject/systems/cylindrical_jellium/scripts/channeling_sic"
TWIN_SCRIPTS = REPO / "ResearchProject/systems/cylindrical_jellium/scripts/channeling_twin"
GS_SUMMARY = TWIN_SCRIPTS / "gs/results/gs/run_summary.txt"
BUILDER = REPO / ".claude/skills/run-notebook/run_notebook_builder.py"
VENV_PY = REPO / "venv/bin/python3"

sys.path.insert(0, str(HERE))
import channeling_stopping as CS  # noqa: E402  (symlinked twin data layer)

RUN = "wp_sic"
GIF_SECONDS = 17.0


def gs_energy_ha() -> float | None:
    if not GS_SUMMARY.is_file():
        return None
    for line in GS_SUMMARY.read_text().splitlines():
        if line.strip().startswith("ground_state_energy_ha"):
            try:
                return float(line.split("=", 1)[1])
            except (IndexError, ValueError):
                return None
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    results = SIC_SCRIPTS / "wp/results" / RUN
    out = HERE / f"wp_{RUN}.ipynb"
    summary = results / "run_summary.txt"

    if not results.is_dir():
        print(f"MISSING: {results}")
        return 1
    if not a.force and not (summary.is_file()
                            and "run_completed = true" in summary.read_text()):
        print(f"SKIP: {results} not complete (use --force to override)")
        return 1

    e_gs = gs_energy_ha()
    cmd = [
        str(VENV_PY), str(BUILDER), str(results), str(out),
        "--run-cpp", str(SIC_SCRIPTS / "wp/run.cpp"),
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

    print(" ".join(cmd), flush=True)
    rc = subprocess.run(cmd, cwd=str(HERE)).returncode
    print(f"wp/{RUN}: {'OK' if rc == 0 else f'FAILED rc={rc}'}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
