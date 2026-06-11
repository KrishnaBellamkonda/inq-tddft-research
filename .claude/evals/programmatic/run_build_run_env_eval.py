#!/usr/bin/env python3
"""Runner for .claude/evals/programmatic/build-run-env.md.
Checks settings.json env wiring + build-run skill consistency. Pure stdlib.
    python3 .claude/evals/programmatic/run_build_run_env_eval.py
"""
from __future__ import annotations

import json
import os

REPO = "/local/data/public/skcb2/tddft"
SETTINGS = os.path.join(REPO, ".claude/settings.json")
SKILL = os.path.join(REPO, ".claude/skills/build-run/SKILL.md")

WANT = {
    "INQ_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share",
    "PSEUDOPOD_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
}


def main() -> int:
    fails = []
    s = json.load(open(SETTINGS))
    env = s.get("env")
    if not isinstance(env, dict):
        fails.append("settings.json has no 'env' block")
        env = {}
    for k, v in WANT.items():
        if env.get(k) != v:
            fails.append(f"env[{k!r}] = {env.get(k)!r} != {v!r}")
    # PATH must NOT be pinned in settings env (would shadow system tools)
    if "PATH" in env:
        fails.append("settings.json env should NOT set PATH (shadows system tools)")

    skill = open(SKILL).read()
    if "settings.json" not in skill:
        fails.append("build-run skill no longer references settings.json env")
    # stale-claim guard: it must NOT claim PATH is hard-coded in settings env
    if "hard-coded in `.claude/settings.json`" in skill and "PATH" in skill.split("hard-coded")[0][-60:]:
        fails.append("build-run skill still claims PATH is hard-coded in settings env")
    if "~/.bashrc" not in skill:
        fails.append("build-run skill no longer documents the bashrc PATH mechanism")

    if fails:
        print(f"FAIL: {len(fails)} checks failed")
        for f in fails:
            print("  -", f)
        return 1
    print("PASS: build-run env wiring + skill consistency (3/3)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
