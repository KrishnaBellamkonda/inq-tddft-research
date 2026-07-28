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

## NOT DONE — gated on GitHub auth (SSH publickey currently failing)

The push to `git@github.com:...` fails `Permission denied (publickey)` — no
ssh-agent reachable from the tool shell; the user's key `~/.ssh/id_ed25519` is
passphrase-protected. To proceed, either the user runs the commands, or provides
`echo "$SSH_AUTH_SOCK"` from their agent so the tool shell can reach the loaded
key. Then, in order:

1. Create public GitHub repo `KrishnaBellamkonda/inq-study`; add remote in
   `inq-study/` and push commit `8c59be9`.
2. On `quantum-stopping-power`: `git submodule add
   https://github.com/KrishnaBellamkonda/inq-study.git inq-study`; commit.
3. Build **`report2/submission-package`** as an ORPHAN in ONE clean commit:
   curated file set (see Goal), with the submodule wired. (Deferred to here
   deliberately — the submodule is intrinsic to it; one clean root commit is the
   reason orphan was chosen.)
4. Push `quantum-stopping-power` and `report2/submission-package` to origin.
5. Delete/ignore the now-superseded `overnight-gaussian-classical`.

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
