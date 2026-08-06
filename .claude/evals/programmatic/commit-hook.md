# Eval: commit-message hook (LOCKED 2026-06-11)

Component (subtask 3): a PreToolUse hook on `git commit` enforcing the
deterministic half of `commit-messages` rule. Evaluator: **programmatic** —
feed each message, assert verdict. Runner lives in `.claude/evals/programmatic/`.

## Reject cases (hook must BLOCK)

| Message (subject / body) | Reason |
|---|---|
| `fix(inqkit): add Claude trailer` | forbidden word "claude" |
| `feature(io): wire Anthropic SDK` | forbidden word "anthropic" |
| `feature(io): add AI agent parser` | standalone "ai" |
| body contains `Co-Authored-By: Claude <…>` | forbidden trailer |
| body contains `🤖 Generated with Claude Code` | forbidden attribution |
| `Fixed the writer bug` | no `action(scope):` prefix |
| `feat(io): add writer` | `feat` ∉ the closed 9-word list (must be `feature`) |
| `feature(io): <73+ char subject …>` | subject > 72 chars |
| `feature: add writer` | missing `(scope)` |

## Accept cases (hook must PASS)

| Message | Note |
|---|---|
| `feature(inqkit): add vec3 subscript operator` | valid action+scope |
| `chore(repo): add CI for portable tiers` | valid |
| `fix(jellium): correct boundary stop in raid run` | "raid"/"main" substrings are NOT the forbidden standalone words |
| `cut(inqview): remove deprecated shim` | valid action `cut` |
| `sim(jellium): add E50 run defs` | valid action `sim` |
| `feature(repo): add .claude/hooks/commit_message_check.py` | `.claude/` path token is exempt (rule targets attribution, not the dir name) |

**Path/identifier exemption:** the forbidden word is allowed when it sits in a
path or identifier token — preceded by `/` or `.`, or followed by `-`/`/`/`_`
(`.claude/`, `docs/claude/skills`, `docs/claude-ecosystem-guide.md`). It is still
blocked as a prose word ("Claude Code", "Co-Authored-By: Claude", "AI agent").
Regression-locked after the hook blocked two of its own commits (the `.claude`
path, then the `docs/claude-…` guide path).

## Closed action-word list (the 9, in precedence order)

`rename, cut, sim, docs, fix, feature, refactor, add, chore`

## Pass criterion

100% of reject cases blocked AND 100% of accept cases passed. Any
false-accept of a forbidden word is a hard fail (the rule's primary purpose).

## Note on `--no-verify`

The hook is advisory at the tool layer (it cannot intercept a raw shell
`git commit` the user runs via `!`); the eval tests the hook's verdict
function, not git internals.

## Boundary cases (added 2026-07-30 after a real escape)

The forbidden-word check exempts matches in path/identifier context by inspecting
the neighbouring character. Because `"" in "-/_"` is True in Python, a match at the
very start or END of the message compared against an ABSENT neighbour and was
wrongly exempted — so `chore(repo): made by Claude` was silently ALLOWED. All 22
pre-existing cases happened to place the forbidden word before another character,
so the suite passed throughout.

The suite therefore now pins, and must keep pinning:

- REJECT: forbidden word as the LAST token of the subject and of the body
  (`made by Claude`, `thanks anthropic`, `built with ai`).
- ACCEPT: genuine path context that ALSO ends the message stays exempt
  (`tidy .claude`, `describe docs/claude/`) — the fix must not over-block.
- A word at index 0 is a defensive sentinel only: it lands in the action-word slot
  and the subject-format rule rejects it either way.

Count is now 28 (13 reject, 12 accept, 3 extract).
