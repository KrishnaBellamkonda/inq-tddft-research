# Eval: file-placement hook (LOCKED 2026-06-11)

Component (subtask 3): a PreToolUse backstop hook on `Write`/`Edit` that WARNS
when a path falls outside the `file-placement` rule's allowlist or inside
upstream `inq/`. Evaluator: **programmatic** — feed each target path, assert
verdict. Backstop only (the rule's dir table stays always-on).

## Warn cases (hook must FLAG)

| Target path | Reason |
|---|---|
| `inq/foo.hpp` | inside upstream INQ (gitignored, do-not-edit) |
| `inq/src/hamiltonian/x.hpp` | inside upstream INQ |
| `notes.md` (repo root) | scattered note — belongs in docs/notes/ |
| `scratch.py` (repo root) | scattered scratch |
| `inq-stack/random_note.md` | not a designated inq-stack subdir |

## Allow cases (hook must PASS silently)

| Target path |
|---|
| `inq-stack/include/inqkit/detail/x.hpp` |
| `inq-stack/python/inqview/analysis/y.py` |
| `docs/plans/z.md` · `docs/notes/scratch.md` · `docs/handovers/t.md` |
| `ResearchProject/systems/jellium/run_x/run.cpp` |
| `.claude/evals/programmatic/new.md` |

## Pass criterion

All warn cases flagged, all allow cases silent. A WARN is non-blocking (the
model may proceed after acknowledging) — false warns are tolerable, but a
missed `inq/` write is a hard fail (protects upstream source).
