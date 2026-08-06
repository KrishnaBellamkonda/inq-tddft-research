#!/usr/bin/env python3
"""Runner for .claude/evals/clusters/cluster-r-figure-standard.md.

Pure-text structural checks (no matplotlib import — CI-safe in the ecosystem
job). The UNITS *values* + axis_label behaviour are tested by the pytest
test_theme.py in the matplotlib CI job; here we check the Cluster-R de-dup:
  1. the theme module DEFINES the canonical unit map (single source);
  2. report-figures OWNS the global standard (references the theme, holds the
     annotation rules incl. axis_label + leader-line);
  3. tufte SHED the project config (no widths table, no full units enumeration);
  4. the inqview LIBRARY (visualisation + pipeline) has no rogue styling outside
     style.py (everything goes through the theme).
    python3 .claude/evals/programmatic/run_cluster_r_eval.py
"""
# NOTE: no `from __future__ import annotations` — these evals are invoked with
# bare `python3` (see docstring above), which on CSD3/RHEL8 is 3.6.8 where that
# import is a hard SyntaxError. All annotations here are plain names; keep them
# 3.6-compatible (no list[str] / X | Y).

import os
import re

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
STYLE = os.path.join(REPO, "inq-stack/python/inqview/visualisation/style.py")
REPORTFIG = os.path.join(REPO, ".claude/skills/report-figures/SKILL.md")
TUFTE = os.path.join(REPO, ".claude/skills/tufte/SKILL.md")
LIB_DIRS = [
    os.path.join(REPO, "inq-stack/python/inqview/visualisation"),
    os.path.join(REPO, "inq-stack/python/inqview/pipeline"),
]


def main() -> int:
    fails = []
    style = open(STYLE).read()
    # 1. canonical unit map present in the theme
    for key in ("energy", "length", "time", "momentum", "stopping_power"):
        if not re.search(rf'"{key}"\s*:', style):
            fails.append(f"theme UNITS missing key {key!r}")
    if '"fs"' not in style:
        fails.append("theme UNITS does not pin time = fs")
    if "def axis_label" not in style:
        fails.append("theme missing axis_label() helper")

    # 2. report-figures owns the global standard
    rf = open(REPORTFIG).read()
    if "inqview.visualisation.style" not in rf:
        fails.append("report-figures does not reference the canonical theme")
    if "leader line" not in rf.lower():
        fails.append("report-figures missing the annotation rules (leader-line)")
    if "axis_label" not in rf:
        fails.append("report-figures does not route units through axis_label")

    # 3. tufte shed the project config
    tufte = open(TUFTE).read()
    if re.search(r"\|\s*`?single`?\s*\|", tufte):
        fails.append("tufte still contains the column-widths table")
    if re.search(r"Stopping power:\s*eV/Bohr\.\s*Energy:", tufte):
        fails.append("tufte still contains the full units enumeration")

    # 4. inqview library is theme-clean (no rogue rcParams/style.use)
    rogue = []
    for d in LIB_DIRS:
        for root, _, files in os.walk(d):
            if "__pycache__" in root:
                continue
            for fn in files:
                if not fn.endswith(".py") or fn == "style.py":
                    continue
                p = os.path.join(root, fn)
                txt = open(p).read()
                if re.search(r"plt\.rcParams\[|plt\.style\.use\(|mpl\.rcParams\.update", txt):
                    rogue.append(os.path.relpath(p, REPO))
    if rogue:
        fails.append("rogue styling outside the theme (use style.apply_theme()): "
                     + ", ".join(rogue))

    if fails:
        print(f"FAIL: {len(fails)} Cluster-R checks failed")
        for f in fails:
            print("  -", f)
        return 1
    print("PASS: Cluster-R figure standard (theme owns units; report-figures owns "
          "the standard; tufte shed project config; library theme-clean)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
