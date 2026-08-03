#!/usr/bin/env python3
"""Refined-analysis notebook for the SIC-corrected channeling run.

Plan: docs/plans/wp-self-interaction-correction.md (reviewed 2026-08-02).

IDENTICAL BY CONSTRUCTION to channeling_twin/refined_analysis.ipynb: this is a
thin wrapper that re-invokes ../channeling_twin/build_refined_notebook.py
(same 29-cell battery, same fits, same figures) with the WP result tree
redirected to the SIC run through the data layer's own env knobs:

    CHAN_WP_RESULTS = scripts/channeling_sic/wp/results     (the corrected run)
    CHAN_CL_RESULTS = scripts/channeling_twin/classical/results  (unchanged ref)

refined.py / channeling_stopping.py are SYMLINKS to the channeling_twin
originals, so any fix there propagates here; the notebook cells import them
from this folder (Path.cwd()) and every figure/CSV lands here.

The user's fit windows from the uncorrected study (9-25; 21-30; 5-20) are
applied by editing the section-6 parameter cell exactly as in the original.

Usage:
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 \
        build_refined_notebook.py [--wp wp_sic] [passthrough args]
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWIN = HERE.parent / "channeling_twin"
REPO = HERE.parents[4]
SCRIPTS = REPO / "ResearchProject/systems/cylindrical_jellium/scripts"

env = dict(os.environ)
env["CHAN_WP_RESULTS"] = str(SCRIPTS / "channeling_sic/wp/results")
env["CHAN_CL_RESULTS"] = str(SCRIPTS / "channeling_twin/classical/results")

args = sys.argv[1:]
if not any(a.startswith("--wp") for a in args):
    args = ["--wp", "wp_sic"] + args

cmd = [sys.executable, str(TWIN / "build_refined_notebook.py"),
       "--out-dir", str(HERE)] + args
print("exec:", " ".join(cmd))
print("  CHAN_WP_RESULTS =", env["CHAN_WP_RESULTS"])
print("  CHAN_CL_RESULTS =", env["CHAN_CL_RESULTS"])
sys.exit(subprocess.run(cmd, env=env, cwd=str(HERE)).returncode)
