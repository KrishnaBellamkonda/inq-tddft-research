#!/usr/bin/env python3
"""Deterministic commit-message validator (commit-messages rule).

Two roles:
  * library  — `check(message) -> Result(ok, reason)` for the eval runner.
  * PreToolUse hook — reads the hook JSON on stdin; if the Bash command is a
    `git commit` with an extractable message, validates it and BLOCKS (exit 2,
    reason on stderr) on failure. If no message can be extracted (editor
    commit, amend without -m), it stays silent (cannot judge → allow).

Enforces the DETERMINISTIC half of the rule only (forbidden words + subject
format/length). Action-word *classification* (which of the 9) is reasoning and
stays in the rule, not here — this hook only checks the action is one of the 9.
"""
# NOTE: deliberately NO `from __future__ import annotations`. Hooks are launched
# with whatever bare `python3` is on PATH, and on this machine (RHEL8/CSD3) that is
# /usr/bin/python3 == 3.6.8, where that future import is a hard SyntaxError
# ("future feature annotations is not defined") -> exit 1 -> Claude Code reports
# "PreToolUse:Bash hook failed with non-blocking status code" and the guard is
# SILENTLY DISABLED. Every annotation here is a plain name (str/bool/Result), so
# the import buys nothing. Keep this module stdlib-only and 3.6-compatible; do not
# use PEP 585 (list[str]) or PEP 604 (X | Y) annotations. Verified 2026-07-30.
import json
import re
import sys

# The closed action-word list (commit-messages rule §5), precedence order.
ACTIONS = ["rename", "cut", "sim", "docs", "fix",
           "feature", "refactor", "add", "chore"]

# Forbidden as STANDALONE PROSE words (so "raid", "main", "detail" are fine).
# The rule targets ATTRIBUTION ("Claude Code", "Co-Authored-By: Claude", "AI
# agent"), not the word as part of a path or identifier token — we reference
# `.claude/`, `docs/claude/`, `docs/claude-ecosystem-guide.md` constantly. A
# match is EXEMPT when it sits in path/identifier context: preceded by "/" or
# "." or followed by "-", "/", or "_".
_FORBIDDEN = re.compile(r"\b(claude|anthropic|ai)\b", re.IGNORECASE)
_PATH_BEFORE = "/."
_PATH_AFTER = "-/_"


def _first_forbidden(message: str):
    """Return the first forbidden PROSE-word match, or None (path/identifier
    occurrences are exempt)."""
    for m in _FORBIDDEN.finditer(message):
        before = message[m.start() - 1] if m.start() > 0 else ""
        after = message[m.end()] if m.end() < len(message) else ""
        # The `before and` / `after and` guards are LOAD-BEARING. When the match sits
        # at the very start/end of the message there is no neighbouring character and
        # these are "", and `"" in "-/_"` is True in Python (the empty string is a
        # substring of everything) — so without the guards a message ENDING in a
        # forbidden word ("chore(repo): made by Claude") was silently ALLOWED. Absence
        # of a neighbour is prose context, not path context. Fixed 2026-07-30.
        if (before and before in _PATH_BEFORE) or (after and after in _PATH_AFTER):
            continue            # part of a path/identifier — allowed
        return m
    return None
_SUBJECT = re.compile(r"^(?:" + "|".join(ACTIONS) + r")\([^()]+\): .+")
_MAX_SUBJECT = 72
# Matches the invocation itself, tolerating extra spacing (`git   commit`).
_GIT_COMMIT = re.compile(r"\bgit\s+commit\b")


class Result:
    __slots__ = ("ok", "reason")

    def __init__(self, ok: bool, reason: str):
        self.ok = ok
        self.reason = reason

    def __repr__(self):
        return f"Result(ok={self.ok}, reason={self.reason!r})"


def check(message: str) -> Result:
    """Validate a full commit message (subject + optional body)."""
    msg = message.strip("\n")
    subject = msg.split("\n", 1)[0] if msg else ""

    # 1. forbidden PROSE words anywhere in the message (path/identifier exempt)
    m = _first_forbidden(msg)
    if m:
        return Result(False, f"forbidden word {m.group(0)!r} (commit-messages rule §1)")

    # 2. subject format: action(scope): description
    if not _SUBJECT.match(subject):
        return Result(False,
                      "subject must be 'action(scope): description' with action in "
                      + ", ".join(ACTIONS))

    # 3. subject length
    if len(subject) > _MAX_SUBJECT:
        return Result(False, f"subject {len(subject)} > {_MAX_SUBJECT} chars")

    return Result(True, "ok")


# ── PreToolUse hook entrypoint ──────────────────────────────────────────────
def _extract_commit_message(command: str):
    """Best-effort message extraction from a `git commit` shell command.
    Returns the message string, or None if it cannot be determined."""
    m = _GIT_COMMIT.search(command)
    if not m:
        return None
    # Only consider text BELONGING to the git-commit invocation, i.e. after it.
    tail = command[m.end():]
    # heredoc: git commit -F - <<'EOF' ... EOF   (also -F-)
    #
    # The -F/--file requirement is LOAD-BEARING. Previously any heredoc anywhere in a
    # command that merely CONTAINED the substring "git commit" was read as the commit
    # message, so writing documentation about commits —
    #   cat >> handover.md <<'EOF' ... every `git commit` was unchecked ... EOF
    # — had its whole body treated as a commit message and was BLOCKED on the first
    # occurrence of a forbidden word. That false positive fires constantly in this repo,
    # which documents the commit rule at length. A real heredoc commit always passes the
    # message via -F, so requiring the flag on the invocation's own first line
    # distinguishes the two. Found and fixed 2026-07-30 (it blocked this very handover).
    first_line = tail.split("\n", 1)[0]
    if re.search(r"(?:^|\s)(?:-F|--file)[=\s]*-?", first_line):
        here = re.search(r"<<-?\s*['\"]?(\w+)['\"]?\n(.*?)\n\1", tail, re.DOTALL)
        if here:
            return here.group(2)
    # -m "msg" / -m 'msg'  (concatenate multiple -m into subject\n\nbody).
    # One capturing group → findall yields the quoted string; strip the quotes.
    parts = re.findall(r"-m\s+(\"(?:\\.|[^\"])*\"|'[^']*')", tail)
    if parts:
        return "\n\n".join(p[1:-1] for p in parts)
    return None


def _hook_main():
    try:
        payload = json.load(sys.stdin)
    except Exception:
        sys.exit(0)  # not a hook invocation / no JSON → do nothing
    if payload.get("tool_name") != "Bash":
        sys.exit(0)
    command = (payload.get("tool_input") or {}).get("command", "")
    message = _extract_commit_message(command)
    if message is None:
        sys.exit(0)  # cannot judge → allow
    res = check(message)
    if not res.ok:
        sys.stderr.write(f"commit-message hook BLOCKED: {res.reason}\n")
        sys.exit(2)  # exit 2 → PreToolUse blocks the tool, reason shown to model
    sys.exit(0)


if __name__ == "__main__":
    _hook_main()
