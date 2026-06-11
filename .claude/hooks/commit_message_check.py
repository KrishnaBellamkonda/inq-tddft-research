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
from __future__ import annotations

import json
import re
import sys

# The closed action-word list (commit-messages rule §5), precedence order.
ACTIONS = ["rename", "cut", "sim", "docs", "fix",
           "feature", "refactor", "add", "chore"]

# Forbidden as STANDALONE words (so "raid", "main", "detail" are fine). The
# negative lookbehind for "." exempts dot-prefixed path tokens like ".claude/"
# (a directory we reference constantly) — the rule targets attribution
# ("Claude Code", "Co-Authored-By: Claude"), not the literal directory name.
_FORBIDDEN = re.compile(r"(?<!\.)\b(claude|anthropic|ai)\b", re.IGNORECASE)
_SUBJECT = re.compile(r"^(?:" + "|".join(ACTIONS) + r")\([^()]+\): .+")
_MAX_SUBJECT = 72


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

    # 1. forbidden words anywhere in the message (subject + body + trailers)
    m = _FORBIDDEN.search(msg)
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
    if "git commit" not in command:
        return None
    # heredoc: git commit -F - <<'EOF' ... EOF   (also -m - / commit -F-)
    here = re.search(r"<<-?\s*['\"]?(\w+)['\"]?\n(.*?)\n\1", command, re.DOTALL)
    if here:
        return here.group(2)
    # -m "msg" / -m 'msg'  (concatenate multiple -m into subject\n\nbody).
    # One capturing group → findall yields the quoted string; strip the quotes.
    parts = re.findall(r"-m\s+(\"(?:\\.|[^\"])*\"|'[^']*')", command)
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
