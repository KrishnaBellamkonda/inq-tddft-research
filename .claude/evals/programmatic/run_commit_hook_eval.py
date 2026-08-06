#!/usr/bin/env python3
"""Runner for .claude/evals/programmatic/commit-hook.md.

Feeds each spec case to the verdict function and asserts the expected verdict.
Pure stdlib; run with any python3:
    python3 .claude/evals/programmatic/run_commit_hook_eval.py
Exit 0 = all pass.
"""
# NOTE: no `from __future__ import annotations` — these evals are invoked with
# bare `python3` (see docstring above), which on CSD3/RHEL8 is 3.6.8 where that
# import is a hard SyntaxError. All annotations here are plain names; keep them
# 3.6-compatible (no list[str] / X | Y).

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
    # BOUNDARY cases — added 2026-07-30 after a real escape. Every case above has
    # the forbidden word followed by another character, so all 22 checks passed
    # while a message ENDING in a forbidden word was silently ALLOWED: the
    # path-context exemption compared the (empty) neighbouring character with
    # `"" in "-/_"`, which is True in Python. Keep at least one word-at-the-very-end
    # and one word-at-index-0 case here forever.
    ("chore(repo): made by Claude", "forbidden claude at END of message"),
    ("feature(io): add writer\n\nthanks anthropic", "forbidden anthropic at END of body"),
    ("feature(io): add writer\n\nbuilt with ai", "standalone ai at END of body"),
    # NOTE this last one is a DEFENSIVE sentinel, not a discriminating test: a word at
    # index 0 sits in the action-word slot, so the subject-format rule rejects it
    # whichever way the before-guard behaves (verified: the pre-fix code also
    # rejected it). The `before and` guard is therefore unreachable-by-verdict today;
    # it is kept for symmetry so a future refactor cannot reintroduce the class.
    ("claude(repo): bad action word", "forbidden claude at index 0 (defensive sentinel)"),
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
    # Companions to the REJECT boundary cases: the fix must not over-block. Genuine
    # path context that happens to sit at the END of the message stays exempt,
    # because there the neighbour is a real "." or "/", not an absent character.
    ("chore(repo): tidy .claude", "dot-claude at END of message still exempt"),
    ("docs(repo): describe docs/claude/", "claude before trailing slash still exempt"),
]

# message-extraction sanity (hook plumbing)
EXTRACT = [
    ('git commit -m "feature(io): add writer"', "feature(io): add writer"),
    ("git commit -F - <<'EOF'\nfix(io): bug\n\nbody line\nEOF", "fix(io): bug\n\nbody line"),
    ("ls -la", None),
    ("git   commit -m \"fix(io): bug\"", "fix(io): bug"),
    ("git commit --amend", None),
    # FALSE-POSITIVE guard (added 2026-07-30 after it blocked a real handover write).
    # A heredoc that merely DISCUSSES git commit is not a commit message. Extraction now
    # requires -F/--file on the invocation's own first line; without this the entire
    # document body was checked and blocked on its first forbidden word. This repo
    # documents the commit rule at length, so the case is common, not hypothetical.
    ("cat >> h.md <<'EOF'\nevery `git commit` was unchecked; mentions Claude\nEOF", None),
    ("cat >> h.md <<'EOF'\ngit commit -F - shown as a docs example\nEOF", None),
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
