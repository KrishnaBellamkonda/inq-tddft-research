# Handover — inq-study: annotated INQ engine + knowledge graph

Task plan: `/local/data/public/skcb2/tddft/docs/plans/inq-study-engine-learning.md`
Glossary context: this is a NEW context (engine internals), distinct from the
root `CONTEXT.md` (scoped to inq-stack unit-testing). Multi-context split not yet
created — see "Not done".

## 2026-06-12 — fork built, knowledge graph + dashboard live

### Done (verified)
- **Fork created:** `/local/data/public/skcb2/tddft/inq-study/` — self-contained,
  buildable copy of upstream `inq/` (rsync, excluded `build/`/`install/`/`.git`).
  432 MB, 180 `src/` C++ files + duplicated `external_libs/`. Own git, pristine
  baseline commit `aea585e`; gitignored by the main repo (`.gitignore` line added).
- **understand-anything run:** scoped to `src/` via
  `inq-study/.understand-anything/.understandignore` (excludes external_libs,
  build, share, examples, tests, python, scripts, cmake, ci, benchmarks).
  196 files analyzed, 13,261 third-party filtered. 13 semantic batches →
  file-analyzer subagents → merge → assemble-review → architecture → tour.
  **Result graph:** `inq-study/.understand-anything/knowledge-graph.json` —
  414 nodes, 525 edges, 9 layers, 13-step tour. Inline validation: **0 issues**,
  21 benign orphan warnings. `meta.json` + `fingerprints.json` written (incremental
  updates enabled).
- **Dashboard:** Vite dev server running (background task `bx0zy9en8`) at
  `http://127.0.0.1:5173/?token=11bcc0b6c6e6a5c702141f66ebc13e70`.
  Token regenerates on restart; graph dir = the project dir via `GRAPH_DIR`.

### How to restart the dashboard (token changes each time)
```bash
cd /home/raid/skcb2/skcb2/tddft/.claude/plugins/cache/understand-anything/understand-anything/2.7.6/packages/dashboard
GRAPH_DIR=/local/data/public/skcb2/tddft/inq-study npx vite --host 127.0.0.1
# then read the "🔑 Dashboard URL: …?token=…" line from stdout
```
(Deps already installed; `pnpm` is MISSING on this host but unnecessary — dashboard
node_modules + core/dist already present. `npx`/node v24 available.)

### Re-analyze / update the graph
`/understand /local/data/public/skcb2/tddft/inq-study` — incremental update keyed
on git diff vs `meta.json` commit hash. Use `--full` to force a rebuild.
`/understand-chat` answers questions against the graph (user requested this).

### Key facts discovered (survive compaction)
- INQ uses absolute `#include <inq/...>` paths resolved at build-configure time →
  structural `imports` edges are sparse (8 total); relationships are mostly inferred
  `depends_on`/`related`/`calls`/`contains`/`inherits`.
- INQ's template/macro-heavy C++ defeats tree-sitter; analyzers read sources
  directly. Node summaries are intentionally physics-grounded (KS Hamiltonian, XC,
  KB projectors, PAW, ETRS/Crank-Nicolson, kick/laser perturbations, Ewald, MPI).

### Not done / pending USER decision (deferred by user on 2026-06-12)
The user said: "Make an understand dashboard first. Then, if required I will use
something else later. /understand-chat too." So these are intentionally OPEN:
1. **Comment-writing workflow / assistant role.** User learns by writing comments
   themselves in `inq-study`. Candidate roles discussed (not chosen): (a) Socratic
   partner on pristine files; (b) seed question-scaffolds; (c) one exemplar file
   then grill. Assistant must NOT write the learning comments — verification +
   physics-grounding only. Revisit now that the dashboard exists.
2. **Glossary strategy.** Likely: root `CONTEXT-MAP.md` + new `inq-study/CONTEXT.md`
   engine glossary, populated lazily as engine terms resolve during commenting.
   Not created yet.

### Provenance / attribution
`inq-study` is a verbatim copy of upstream INQ (Andrade & Correa, LLNL) at
`inq/.git` HEAD `44f73d95`. `COPYING`/`AUTHORS` preserved (GPL). It is a learning
artifact, NOT a production run target — the pristine `inq/` remains the `inq-run`
build target.
