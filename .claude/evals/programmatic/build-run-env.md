# Eval: build-run env claim (LOCKED 2026-06-11)

Component (subtask 3): fix the stale `build-run` skill claim that env vars are
hard-coded in `.claude/settings.json` — add the `env` block so the claim is true.

Evaluator: **programmatic** — assert the wiring exists.

## Cases

1. `.claude/settings.json` contains an `env` block.
2. The `env` block defines `INQ_SHARE_PATH` and `PSEUDOPOD_SHARE_PATH` with the
   canonical values:
   - `INQ_SHARE_PATH = /local/data/public/skcb2/tddft/inq/install/share`
   - `PSEUDOPOD_SHARE_PATH = /local/data/public/skcb2/tddft/inq/install/share/pseudopod`
3. `PATH` (or the settings env) includes `/local/data/public/skcb2/tddft/shared/bin`.
4. Consistency: the `build-run` SKILL.md text describing the env block matches
   the actual `settings.json` keys (no stale claim).

## Pass criterion

All 4 true. This is a wiring/consistency check, not a behavioural test — it
guarantees the skill's documented contract matches reality.
