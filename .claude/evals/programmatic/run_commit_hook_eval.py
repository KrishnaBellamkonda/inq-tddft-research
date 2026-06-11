#!/usr/bin/env python3
"""Runner for .claude/evals/programmatic/commit-hook.md.

Feeds each spec case to the verdict function and asserts the expected verdict.
Pure stdlib; run with any python3:
    python3 .claude/evals/programmatic/run_commit_hook_eval.py
Exit 0 = all pass.
"""
from __future__ import annotations

import os
import sys

# import the hook verdict function
HOOKS = os.path.join(os.path.dirname(__file__), "..", "..", "hooks")
sys.path.insert(0, os.path.abspath(HOOKS))
from commit_message_check import check, _extract_commit_message  # noqa: E402

# (message, expected_ok, label)
REJECT = [
    ("fix(inqkit): add Claude trailer", "forbidden claude"),
    ("feature(io): wire Anthropic SDK", "forbidden anthropic"),
    ("feature(io): add AI agent parser", "standalone ai"),
    ("feature(io): add writer\n\nCo-Authored-By: Claude <x@y.z>", "forbidden trailer"),
    ("feature(io): add writer\n\n\N{ROBOT FACE} Generated with Claude Code", "forbidden attribution"),
    ("Fixed the writer bug", "no action(scope)"),
    ("feat(io): add writer", "feat not in 9-word list"),
    ("feature(io): " + "x" * 70, "subject > 72"),
    ("feature: add writer", "missing scope"),
]
ACCEPT = [
    ("feature(inqkit): add vec3 subscript operator", "valid"),
    ("chore(repo): add CI for portable tiers", "valid"),
    ("fix(jellium): correct boundary stop in raid run", "raid substring ok"),
    ("cut(inqview): remove deprecated shim", "valid cut"),
    ("sim(jellium): add E50 run defs", "valid sim"),
    ("feature(coronene+inqview): add leed overlay", "multi-scope ok"),
    ("feature(repo): add .claude/hooks/commit_message_check.py", "dot-claude path ok"),
    ("docs(repo): describe the .claude/evals layout", "dot-claude path ok 2"),
    ("docs(repo): add docs/claude-ecosystem-guide.md", "claude in slash+hyphen path"),
    ("docs(repo): polish docs/claude/skills notes", "claude mid-path (slashes)"),
]

# message-extraction sanity (hook plumbing)
EXTRACT = [
    ('git commit -m "feature(io): add writer"', "feature(io): add writer"),
    ("git commit -F - <<'EOF'\nfix(io): bug\n\nbody line\nEOF", "fix(io): bug\n\nbody line"),
    ("ls -la", None),
]


def main() -> int:
    fails = []
    for msg, label in REJECT:
        r = check(msg)
        if r.ok:
            fails.append(f"REJECT case wrongly PASSED [{label}]: {msg!r}")
    for msg, label in ACCEPT:
        r = check(msg)
        if not r.ok:
            fails.append(f"ACCEPT case wrongly REJECTED [{label}]: {msg!r} -> {r.reason}")
    for cmd, expected in EXTRACT:
        got = _extract_commit_message(cmd)
        if got != expected:
            fails.append(f"EXTRACT mismatch for {cmd!r}: got {got!r} != {expected!r}")

    total = len(REJECT) + len(ACCEPT) + len(EXTRACT)
    if fails:
        print(f"FAIL: {len(fails)}/{total} checks failed")
        for f in fails:
            print("  -", f)
        return 1
    print(f"PASS: {total}/{total} commit-hook checks "
          f"({len(REJECT)} reject, {len(ACCEPT)} accept, {len(EXTRACT)} extract)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
