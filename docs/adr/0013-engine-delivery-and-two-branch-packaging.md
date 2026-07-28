# ADR 0013 — Engine delivery and two-branch packaging for device migration

- **Status:** accepted
- **Date:** 2026-07-28
- **Context scope:** repository packaging for moving the project to a second
  device and continuing runs there — the `quantum-stopping-power` and
  `report2/submission-package` branches, the `inq-study` submodule, and the
  `setup.sh` + `inq-local.patch` engine-bootstrap artifacts.
- **Relates to:** `inq-immutable` rule (upstream `inq/` is read-only except
  sanctioned deltas); `.gitignore` run-output policy.

## Context

The work needs to continue on a different device: build the engine, regenerate
ground states, launch new runs, and analyse them. The main repo tracks only
source; three heavy dependencies live outside git and had to be reconciled with
"download what's needed, lose nothing, ship no results":

- **`inq/` (2.5 GB)** — pristine upstream INQ at commit
  `44f73d9527ab677f38ed2138c2e83a28a5ab6c79`, plus **two uncommitted local
  edits**: the CUB fix (`external_libs/gpurun/include/gpu/reduce.hpp`, required
  for CUDA 12.5+) and a read-only `ham()` accessor
  (`src/real_time/viewables.hpp`) that `inqkit` KS-energy observables depend on.
  Public upstream is re-clonable; only the 2-file delta is precious.
- **`inq-study/` (437 MB source)** — the project-modified engine that adds
  complex-absorbing-potential (CAP) support. Stock upstream **cannot compile a
  CAP run**. Its entire delta was **uncommitted and hosted nowhere** — one disk
  away from permanent loss.
- **Ground states (112 GB)** and all run outputs — gitignored; far too large and
  explicitly out of scope ("no results of any kind").

Two consumers were needed: a full working branch for continued development
(with `docs/` and `.claude/` for context), and a clean showcase of the codebase
with neither.

## Decision

1. **`inq-study` → its own public repo, consumed as a git submodule.** Its
   uncommitted delta is first committed inside its own repo (build trees
   gitignored), then pushed to a dedicated public GitHub repo and referenced as
   a submodule (HTTPS URL, so anonymous clones of the showcase resolve). This
   backs up the irreplaceable engine *and* keeps its 437 MB out of the research
   repo's history.

2. **`inq/` → pinned clone + committed patch, not vendored.** `setup.sh` clones
   upstream at the pinned commit and applies `inq-local.patch` (both sanctioned
   edits). The patch is committed, so the 2-file delta is backed up without
   carrying 2.5 GB. This also records that `inq/` has **two** sanctioned deltas,
   not one — the `inq-immutable` rule (which names only the CUB fix) is stale and
   should be updated.

3. **Ground states are regenerated, never shipped.** The `save_gs/*/run.cpp`
   builders are part of the tracked scripts; the new device rebuilds the 112 GB
   locally.

4. **Two branches off the current work head:**
   - **`quantum-stopping-power`** — normal branch, everything (docs, `.claude`,
     scripts), no results; the active development branch.
   - **`report2/submission-package`** — an **orphan** branch (clean root, no
     history) showcasing the codebase: `inq-stack/` (inqkit + inqview + tests),
     `ResearchProject/` scripts, the `inq-study` submodule, and root
     `README.md` + `setup.sh` + `inq-local.patch`. Excludes `.claude/`, `docs/`,
     `CONTEXT.md`, `CLAUDE.md` — even from history.

## Consequences

- The irreplaceable CAP engine is committed and (once pushed) backed up; the
  central loss risk is retired.
- A new device needs a submodule init and a build step (`setup.sh`) — one extra
  action versus a vendored monorepo, in exchange for a lean, unbloated repo.
- The orphan showcase shares no history with `quantum-stopping-power`, so
  refreshing it means regenerating the snapshot rather than merging — accepted,
  because it is a presentation artifact, not a development branch.
- The submodule is only as accessible as the `inq-study` repo: it must stay
  public for the showcase to resolve for anyone (including a fresh, unauthed
  clone).

## Alternatives considered

- **Vendor `inq-study` (and/or `inq`) into the branches.** Self-contained but
  bloats every branch by hundreds of MB–GBs and mixes engine source into the
  research history permanently. Rejected in favour of submodule + patch.
- **Filtered-tip showcase** (delete `.claude`/`docs` in a commit off HEAD).
  Leaves both folders recoverable in history and pays for that history on clone.
  Rejected in favour of an orphan clean root.
- **Document-only engine delivery** (copy `inq-study` by hand). Lightest, but
  leaves the modified engine unbacked-up — the exact failure this ADR prevents.
