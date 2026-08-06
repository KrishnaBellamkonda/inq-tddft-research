# Eval: build-run env claim (LOCKED 2026-06-11)

Component (subtask 3): fix the stale `build-run` skill claim that env vars are
hard-coded in `.claude/settings.json` — add the `env` block so the claim is true.

Evaluator: **programmatic** — assert the wiring exists.

## Cases

1. `.claude/settings.json` contains an `env` block.
2. The `env` block defines the two additive, safe vars, pointing at the INQ install
   inside **this** repo — i.e. `<repo-root>/inq/install/share` and
   `<repo-root>/inq/install/share/pseudopod`.

   The runner DERIVES these from its own location rather than hardcoding them
   (updated 2026-07-30). They were previously pinned to the previous device's
   literal `/local/data/public/skcb2/tddft/...`, which after the CSD3 migration
   asserted a path that does not exist, so the eval failed for the wrong reason.
   `settings.json` still stores an absolute path (it is machine-specific config);
   deriving in the runner means a repo move now fails this eval LOUDLY instead of
   leaving a stale value in place unnoticed.
3. Consistency: the `build-run` SKILL.md text matches reality — it states these
   two vars are in the `settings.json` env, and that PATH entries (`shared/bin`,
   `pyenv/shims`) come from `~/.bashrc` (the Bash shell sources the user profile),
   NOT from settings.json. No stale "PATH is hard-coded in settings env" claim.

## Design note (why PATH is NOT set in settings env)

Overriding `PATH` wholesale in `settings.json` env risks shadowing system tools
(git, python3, nvcc) for every Bash call. `shared/bin` is already on PATH via
`~/.bashrc`, which the Bash tool's shell sources. So only the two additive
`*_SHARE_PATH` vars are pinned in settings; PATH stays profile-driven.

## Pass criterion

All 3 true. A wiring/consistency check: the skill's documented contract matches
the actual settings.json + the bashrc PATH mechanism.
