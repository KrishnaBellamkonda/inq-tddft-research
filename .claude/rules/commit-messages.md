# Rule: Commit Messages

Apply to: every git commit in any repo under `/local/data/public/skcb2/tddft/`.

## Rules

1. **Forbidden words.** Commit messages (subject and body) must not contain the words `claude`, `anthropic`, or `ai` (case-insensitive, including substrings only when they form the standalone word — e.g. avoid "Claude Code", "AI agent", "Anthropic SDK"; "main" or "raid" are fine).

   **Why:** the user does not want generated-by-assistant attribution or AI/vendor branding leaking into the public git history.

   **How to apply:**
   - Strip any "Co-Authored-By: Claude …" or "Generated with Claude Code" trailers — do not append them, even though they are the default in some workflows.
   - When summarising automated edits to `.claude/`, refer to them as "internal folder and environmental configurations" (or similar neutral phrasing) rather than naming the assistant.
   - If a substantive technical term legitimately requires one of these words (e.g. citing an "AI" paper), surface it to the user before committing and let them write or approve the message.

2. **Two-commit hygiene for mixed work.** When a single set of staged changes mixes (a) production research changes and (b) infrastructure/internal-config changes, split them into separate commits with clear, scoped subjects rather than a single omnibus commit.

3. **Pre-commit check.** Before running `git commit`, scan the drafted message for the forbidden words above. If any are present, rewrite before committing — never amend after the fact in a way that rewrites already-pushed history.

4. **Subject format.** Every commit subject is `action(scope): description`.

   - `action` is one of the nine action words in rule 5, classified by the
     **first matching rule, top to bottom**.
   - `scope` is the top-level component/system the change belongs to (rule 6).
   - `description` is a lowercase, imperative phrase ("add", "fix" — not
     "added", "fixes"), no trailing period, subject line ≤ 72 characters.

5. **Action words (classify by FIRST match, top → bottom).** The list is
   closed — every change maps to exactly one. If two could apply, the one
   higher in this list wins.

   | Action | Use for |
   |---|---|
   | `rename` | pure file move/rename, content unchanged |
   | `cut` | pure removal of files/content |
   | `sim` | simulation run provenance (`run.cpp`, `analyse.py`, run configs) |
   | `docs` | docs, journals, handovers, reports, plans, source/literature notes |
   | `fix` | bug fix in existing code |
   | `feature` | new code capability / behaviour |
   | `refactor` | restructure with no behaviour change |
   | `add` | new non-code, non-doc asset (catch-all for net-new content) |
   | `chore` | tooling, gitignore, build, config, repo hygiene |

   **Why precedence:** `rename`, `cut`, `sim`, `docs`, and `add` would
   otherwise overlap with `feature`/`chore`. The ordering makes the
   classification deterministic — a pure move of a run script is `rename`,
   not `sim`; a new run definition is `sim`, not `add`.

6. **Scope token.** `scope` is a component/system name, not a path:
   `inqview`, `inqkit`, `jellium`, `coronene`, `qke`, or `repo` for
   repo-wide changes (e.g. gitignore). Multi-scope is allowed with `+`
   when a change genuinely spans two components, e.g.
   `feature(coronene+inqview): …`. Prefer splitting over a third scope.

7. **Body policy.** Include a body — a `-`-prefixed bullet list of what
   changed and why, wrapped at ~72 columns, separated from the subject by a
   blank line — whenever a commit **spans more than one file** or carries
   **physics/run provenance** (run IDs, energies, parameter values). A
   subject line alone is acceptable only for a single-file, self-evident
   change. Record concrete provenance (e.g. `S(v=10.5)=0.021 eV/Bohr`,
   `run_wp_n162_L50_E100`) in the body, never invented values.

8. **Example.**

   ```
   sim(jellium): add E50–E300 classical+wp run defs

   - run.cpp + analyse.py for 14 energy-sweep runs
   - shared configs: highdens L30, sigma1 variants
   - logs/results gitignored; provenance only
   ```
