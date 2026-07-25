# Plan — `inq-study`: annotated, buildable fork of the INQ engine for learning TDDFT

Status: **active** (started 2026-06-12)
Owner: user (chiddukanna). Assistant role: scaffolding + Socratic grilling (see §5).

## 1. Goal

Build a personal, buildable replica of the upstream INQ engine (`inq/`) that the
user annotates by hand, file by file, to (a) understand what each source file
does and how the engine is structured, and (b) connect each component back to
the physics of (TD)DFT. Learning method is **the user writing the comments
themselves** — the assistant must not do that learning for them.

## 2. Resolved decisions (grilling session 2026-06-12)

| Decision | Resolution | Rationale |
|---|---|---|
| Target codebase | Upstream `inq/` engine, all 180 `src/` files | Physics lives here (KS H, XC, real-time propagation, perturbations), not in inqkit/inqview which are downstream post-processing |
| Replica form | **Buildable, self-contained** fork | User wants to be able to compile/experiment, not just read |
| `external_libs/` | Physically duplicated (~427 MB) | Self-contained; survives even if `inq/` is deleted |
| Excluded from copy | `build/` (728 MB), `install/` (144 MB), `.git` (1.2 GB) | Regenerable build output / foreign upstream history |
| Name | `inq-study/` | "my study of inq"; clear to a future reader |
| Location | repo root, sibling of `inq/`, `Tutorial/`, `QuantumKickExtension/` | Matches existing big-external-tree pattern |
| Git | own fresh `git init`; gitignored by the main repo | Same pattern as `inq/`, `Tutorial/`; user's comments are the first commits in *its* history |
| First deliverable | `/understand` knowledge graph → **dashboard** → **chat** | User explicitly wants the interactive map first, before committing to a comment workflow |
| Build order | Fork first, THEN `/understand inq-study` | Graph nodes are path-bound; binding them to the commenting copy avoids a rebuild |

## 3. Deferred decisions (user will decide after seeing the dashboard)

- **The comment-writing workflow / assistant role.** Candidates discussed:
  (a) Socratic partner on pristine files; (b) seed question-scaffolds (questions,
  never answers) atop each file; (c) one fully-commented exemplar file then grill.
  User said: "Make an understand dashboard first. Then, if required I will use
  something else later." → revisit once the map is visible.
- **Graph scope:** all 180 files vs a physics-only subset for a less crowded
  first dashboard. Currently: all 180 (full engine).
- **Glossary/context strategy:** the existing root `CONTEXT.md` is scoped to the
  *inq-stack unit-testing rejuvenation* task. The engine-internals vocabulary is
  a **different context**. Likely outcome: a root `CONTEXT-MAP.md` + a new
  `inq-study/CONTEXT.md` engine glossary, created lazily as engine terms resolve.

## 4. Execution steps

1. [done] `rsync -a --exclude={build,install,.git} inq/ inq-study/` — copy (432 MB, 180 src files).
2. [done] Added `inq-study/` to the main repo `.gitignore`.
3. [done] `git init` inside `inq-study/`; pristine baseline commit `aea585e`
   (its own `.gitignore` excludes `build/`, `install/`, `.understand-anything/`).
4. [done] Staged `inq-study/.understand-anything/.understandignore` scoping
   analysis to `src/` only (excludes external_libs/build/share/examples/…).
5. [done] `/understand /…/inq-study` → knowledge graph. Scope: 196 files
   (180 src C++ + top-level/CI; 13,261 third-party filtered by `.understandignore`).
   13 semantic batches, file-analyzer subagents. Result: **414 nodes, 525 edges,
   9 layers, 13-step tour**. 0 validation issues. Graph at
   `inq-study/.understand-anything/knowledge-graph.json`.
6. [done] `/understand-dashboard` → Vite server live at
   `http://127.0.0.1:5173/?token=11bcc0b6c6e6a5c702141f66ebc13e70`
   (background task `bx0zy9en8`; token regenerates each restart).
7. `/understand-chat` available for Q&A against the graph (user requested).
8. [PENDING USER] Revisit §3 deferred decisions after exploring the dashboard.

### Knowledge-graph facts (for future sessions)
- INQ uses absolute `#include <inq/...>` paths resolved at build-configure time,
  so structural `imports` edges are sparse (only 8). Cross-file relationships are
  represented as inferred `depends_on`/`related`/`calls`/`contains`/`inherits`.
- INQ's heavy template/macro C++ defeats tree-sitter; file-analyzers read sources
  directly to produce physics-grounded node summaries.
- 9 layers: Public Interface & Entry (36) · Systems & Electronic State (13) ·
  Basis/Grids/Math (20) · Hamiltonian & Physics (15) · SCF/Eigensolvers/LinAlg (25) ·
  Real-Time TDDFT & Perturbations (18) · Observables & Field Ops (26) ·
  Parallel & Core Infra (27) · Build/CI/Docs (16).
- Re-run / update graph: `/understand /…/inq-study` (incremental on git diff).
  Restart dashboard: background the dashboard Vite command (see handover).

## 5. Assistant role (until §3 is decided)

Scaffolding + verification only:
- Stand up the fork, the knowledge graph, the dashboard, the chat.
- Provide a physics-ordered reading order (data flow:
  `interface → systems → hamiltonian → ground_state/real_time → observables`).
- As the user comments, **grill and verify** their understanding against the
  actual code and ground the TDDFT physics in trustworthy sources
  (`.claude/rules/scientific-grounding.md`). Do **not** write the learning
  comments.

## 6. Notes / invariants

- The pristine `inq/` stays untouched and remains the build target for `inq-run`.
- `inq-study/` is a learning artifact, NOT a production run target.
- Source attribution: `inq-study/` is a verbatim copy of upstream INQ
  (`inq/.git` HEAD `44f73d95`, "Merge branch 'inq_paw_dev' into 'master'").
  Keep `COPYING`/`AUTHORS` intact to preserve the GPL licence + attribution.
