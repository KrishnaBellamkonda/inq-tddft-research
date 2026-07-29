# Handover — device-migration branch packaging

Rolling handover for splitting the repo into two branches for continuing work on
a second device, plus safe delivery of the INQ engine. Decisions are recorded in
`docs/adr/0013-engine-delivery-and-two-branch-packaging.md`.

## Goal

Move to a new device and continue runs + analysis. Produce two branches on
`KrishnaBellamkonda/inq-tddft-research`, no run results in either:

- **`quantum-stopping-power`** — full working branch (everything incl. `docs/`,
  `.claude/`), for active development on the new device.
- **`report2/submission-package`** — **orphan** (clean root, no history) showcase
  of the codebase: `inq-stack/` (inqkit+inqview+tests), `ResearchProject/`
  scripts, `inq-study` submodule, root `README.md`+`setup.sh`+`inq-local.patch`.
  EXCLUDES `.claude/`, `docs/`, `CONTEXT.md`, `CLAUDE.md` (even from history).

Both branch from the current `overnight-gaussian-classical` HEAD; the plain
`overnight-gaussian-classical` name is NOT pushed.

## Engine delivery (settled)

- **`inq-study`** (project-modified engine, CAP support): own **public** GitHub
  repo → git **submodule** (HTTPS URL). Its delta was uncommitted; now committed.
- **`inq/`** (upstream, pinned `44f73d9527ab677f38ed2138c2e83a28a5ab6c79`):
  `setup.sh` clones + applies `inq-local.patch`. Two sanctioned deltas: CUB fix
  (`external_libs/gpurun/include/gpu/reduce.hpp`) + read-only `ham()` accessor
  (`src/real_time/viewables.hpp`, needed by inqkit KS-energy observables).
- **Ground states (112 GB)** and outputs: NOT shipped; regenerate via
  `ResearchProject/systems/**/save_gs/*/run.cpp` on the new device.

## DONE (local, no auth) — 2026-07-28

- Committed `inq-study` CAP delta INSIDE its own repo (anti-loss). Commit
  `8c59be9` in `/local/data/public/skcb2/tddft/inq-study` (build-*/ gitignored).
  This was the critical loss risk — the modified engine had been uncommitted and
  unhosted.
- `/local/data/public/skcb2/tddft/inq-local.patch` — generated, verified it
  reverse-applies to current `inq/` (⇒ forward-applies to the pinned commit).
- `/local/data/public/skcb2/tddft/setup.sh` — engine bootstrap (clone pinned inq,
  apply patch, build; init inq-study submodule; pip install inqview).
- `docs/adr/0013-engine-delivery-and-two-branch-packaging.md`.
- Committed on `overnight-gaussian-classical` (= base for both branches):
  - `c1e6f96` docs(repo): ADR 0013
  - `e9dba4f` chore(repo): setup.sh + inq-local.patch
  - `1ca25e6` cut(coronene): remove stray tracked profile.dat
- Created local branch **`quantum-stopping-power`** (now the active branch),
  committed the handover (`95105c5`).
- Folded the mid-session background work into `quantum-stopping-power` (tree now
  fully clean): `67f6f65` chore(repo) (gitignore `*_DONE.*` + settings),
  `eb21fe4` sim(jellium) (qsp5 momentum / wp-cap plateau / vacuum pipeline),
  `e0b9976` docs(jellium) (qsp5 handover, report-2 catalogue, cap notes). All
  ~1 MB source/docs; run outputs stayed gitignored.

## DONE — 2026-07-29 (local build of submodule + orphan)

- User pushed `quantum-stopping-power` to origin (at `cf5c677`).
- Wired `inq-study` submodule into `quantum-stopping-power` (`1e8a674`):
  un-ignored `inq-study/` + `.gitmodules`, HTTPS URL, gitlink at `8c59be9`,
  branch master. Set `inq-study` remote to the (not-yet-created) GitHub URL.
- Built **`report2/submission-package`** as an ORPHAN via plumbing (no
  working-tree switch — a background run was live): `efdd52f`, 17.1 MB, 1336
  files. Includes inq-stack, ResearchProject scripts, inq-study submodule,
  root README.md + setup.sh + inq-local.patch. EXCLUDES `.claude/`, `docs/`,
  `CONTEXT.md`, `CLAUDE.md`, `ResearchProject/literature/`, all `*.ipynb`.
- Untracked `ResearchProject/literature/` on `quantum-stopping-power`
  (`29df188`) + gitignored it (files kept on disk). It held a copyrighted
  textbook epub.

### CAVEAT — copyrighted epub in already-pushed qsp history
`ResearchProject/literature/tddft/tddft-concept-and-applications.epub` (27.8 MB)
is in the history pushed to origin (`cf5c677` and ancestors). Removed from the
tree now, but it REMAINS in history. If the `inq-tddft-research` repo is public,
scrub it with `git filter-repo --path <epub> --invert-paths` + force-push.
`report2/submission-package` (orphan) is clean — no epub in its history.

## NOT DONE — pushes (user drives) + inq-study repo creation

No auth in the tool shell (no ssh-agent, no `gh`); the user drives all pushes
with their passphrase-unlocked key. Remaining, in order:

1. User creates an EMPTY public GitHub repo `KrishnaBellamkonda/inq-study`
   (no README). Remote already set locally.
2. User pushes `inq-study`: `git -C inq-study push -u origin master` (commit
   `8c59be9`). MUST happen before the submodule resolves for anyone.
3. User re-pushes `quantum-stopping-power` (now ahead of origin by `1e8a674`
   submodule + `29df188` literature-untrack): `git push origin
   quantum-stopping-power`.
4. User pushes the orphan showcase: `git push -u origin
   report2/submission-package`.
5. (Optional) delete the superseded `overnight-gaussian-classical`; scrub the
   epub from qsp history if the repo is public (see CAVEAT).

## OPEN decisions / caveats

- Surprise mid-session work: RESOLVED — folded into `quantum-stopping-power`
  (`67f6f65`/`eb21fe4`/`e0b9976`).
- **`inq-immutable` rule is stale**: it documents only the CUB fix; `inq/`
  actually has TWO sanctioned deltas (add the `ham()` accessor). Update the rule
  + its memory so a future session doesn't "restore upstream purity" and break
  KS-energy runs.
- GS-loading `run.cpp` must use RELATIVE paths (not `/local/data/public/...`) to
  resolve on the new device — verify before relying on replication.
- report2 root already has a tracked `README.md` (project readme, 5.9 KB) — for
  the orphan showcase, decide whether to reuse it or write a showcase-specific
  README.
