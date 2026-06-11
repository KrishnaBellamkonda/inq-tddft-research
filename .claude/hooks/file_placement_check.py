#!/usr/bin/env python3
"""Backstop for the file-placement rule (deterministic half).

Roles:
  * library — `classify(rel_path) -> (level, reason)` where level is "allow" or
    "warn"; tested by the eval runner.
  * PreToolUse hook — on Write/Edit, classify the target; on "warn" print a
    NON-BLOCKING nudge to stderr and exit 0 (the rule's dir table stays the
    primary guide; this only catches obvious slips like writing into upstream
    inq/ or scattering a file at the repo root).

The rule's allowlist is the source of truth; this mirrors its hard constraints.
"""
from __future__ import annotations

import json
import os
import sys

REPO = "/local/data/public/skcb2/tddft"

ALLOW = "allow"
WARN = "warn"

# Designated top-level trees (file-placement rule). `.claude/` and `shared/`
# (configs/templates) are included; inq-stack is restricted to its real subdirs.
_ALLOW_PREFIXES = (
    "docs/", "ResearchProject/", "Tutorial/", "QuantumKickExtension/",
    "shared/", ".claude/",
)
_INQSTACK_OK = ("inq-stack/include/inqkit/", "inq-stack/python/inqview/",
                "inq-stack/tests/")
_ROOT_ALLOW = {"CLAUDE.md", "CONTEXT.md", "CONTEXT-MAP.md", "README.md",
               ".gitignore"}


def classify(rel: str):
    """rel: POSIX path relative to the repo root. Returns (level, reason)."""
    rel = rel.lstrip("./") if rel.startswith("./") else rel
    parts = rel.split("/")

    if parts[0] == "inq":
        return WARN, "inside upstream INQ (inq/) — do not edit; new code goes in inq-stack/"
    if parts[0] == "inq-stack":
        if rel.startswith(_INQSTACK_OK):
            return ALLOW, ""
        return WARN, ("not a designated inq-stack subdir "
                      "(include/inqkit/, python/inqview/, tests/)")
    if rel.startswith(_ALLOW_PREFIXES):
        return ALLOW, ""
    if len(parts) == 1:                      # a file at the repo root
        if rel in _ROOT_ALLOW:
            return ALLOW, ""
        return WARN, "scattered repo-root file — use a designated dir (docs/notes/, …)"
    return WARN, f"'{parts[0]}/' is not a designated directory"


def _to_rel(file_path: str):
    """Absolute or relative path -> POSIX path from repo root, or None if the
    target is outside the repo (not our concern)."""
    if not file_path:
        return None
    ap = os.path.abspath(file_path if os.path.isabs(file_path)
                         else os.path.join(REPO, file_path))
    if ap == REPO or not ap.startswith(REPO + os.sep):
        return None
    return os.path.relpath(ap, REPO).replace(os.sep, "/")


def _hook_main():
    try:
        payload = json.load(sys.stdin)
    except Exception:
        sys.exit(0)
    if payload.get("tool_name") not in ("Write", "Edit"):
        sys.exit(0)
    fp = (payload.get("tool_input") or {}).get("file_path", "")
    rel = _to_rel(fp)
    if rel is None:
        sys.exit(0)
    level, reason = classify(rel)
    if level == WARN:
        sys.stderr.write(f"file-placement nudge ({rel}): {reason}\n")
    sys.exit(0)   # always non-blocking


if __name__ == "__main__":
    _hook_main()
