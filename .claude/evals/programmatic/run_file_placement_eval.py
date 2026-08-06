#!/usr/bin/env python3
"""Runner for .claude/evals/programmatic/file-placement-hook.md.
Feeds each spec path to classify() and asserts warn/allow. Pure stdlib.
    python3 .claude/evals/programmatic/run_file_placement_eval.py
"""
# NOTE: no `from __future__ import annotations` — these evals are invoked with
# bare `python3` (see docstring above), which on CSD3/RHEL8 is 3.6.8 where that
# import is a hard SyntaxError. All annotations here are plain names; keep them
# 3.6-compatible (no list[str] / X | Y).

import os
import sys

HOOKS = os.path.join(os.path.dirname(__file__), "..", "..", "hooks")
sys.path.insert(0, os.path.abspath(HOOKS))
from file_placement_check import classify, _to_rel, WARN, ALLOW, REPO  # noqa: E402

# (relative-or-absolute path, expected level)
WARN_CASES = [
    "inq/foo.hpp",
    "inq/src/hamiltonian/x.hpp",
    "notes.md",
    "scratch.py",
    "inq-stack/random_note.md",
]
ALLOW_CASES = [
    "inq-stack/include/inqkit/detail/x.hpp",
    "inq-stack/python/inqview/analysis/y.py",
    "docs/plans/z.md",
    "docs/notes/scratch.md",
    "docs/handovers/t.md",
    "ResearchProject/systems/jellium/run_x/run.cpp",
    ".claude/evals/programmatic/new.md",
]
# absolute-path + outside-repo plumbing
PLUMBING = [
    (f"{REPO}/inq/foo.hpp", WARN),
    (f"{REPO}/docs/plans/a.md", ALLOW),
    ("/tmp/anything.txt", None),     # outside repo -> _to_rel None -> not judged
]


def lvl(path):
    rel = _to_rel(path)
    if rel is None:
        return None
    return classify(rel)[0]


def main() -> int:
    fails = []
    for p in WARN_CASES:
        if classify(p)[0] != WARN:
            fails.append(f"expected WARN, got ALLOW: {p}")
    for p in ALLOW_CASES:
        if classify(p)[0] != ALLOW:
            fails.append(f"expected ALLOW, got WARN: {p} -> {classify(p)[1]}")
    for p, exp in PLUMBING:
        got = lvl(p)
        if got != exp:
            fails.append(f"plumbing {p}: got {got} != {exp}")

    total = len(WARN_CASES) + len(ALLOW_CASES) + len(PLUMBING)
    if fails:
        print(f"FAIL: {len(fails)}/{total} checks failed")
        for f in fails:
            print("  -", f)
        return 1
    print(f"PASS: {total}/{total} file-placement checks "
          f"({len(WARN_CASES)} warn, {len(ALLOW_CASES)} allow, {len(PLUMBING)} plumbing)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
