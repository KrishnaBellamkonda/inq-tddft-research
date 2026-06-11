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

## Closed action-word list (the 9, in precedence order)

`rename, cut, sim, docs, fix, feature, refactor, add, chore`

## Pass criterion

100% of reject cases blocked AND 100% of accept cases passed. Any
false-accept of a forbidden word is a hard fail (the rule's primary purpose).

## Note on `--no-verify`

The hook is advisory at the tool layer (it cannot intercept a raw shell
`git commit` the user runs via `!`); the eval tests the hook's verdict
function, not git internals.
