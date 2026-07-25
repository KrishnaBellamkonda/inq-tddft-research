# Handover: inq-study CAP enablement + engine validation (DEFERRED)

Status: **change made, UNVALIDATED — deferred by user decision 2026-06-14.**
Plan: `docs/plans/absorbing-cap-engine.md`. Parent handover:
`docs/handovers/absorbing-boundary.md`. Validation task: **#7** (task list).
Rule: `inq/` immutable (`.claude/rules/inq-immutable.md`) — engine edits in
`inq-study/` only.

## TL;DR for whoever resumes this

A **1-line `inq-study` change makes INQ's in-built `perturbations::absorbing`
(a sin² CAP) functional**, and CAP simulations run + give correct physics. But the
change is **not yet validated against INQ's own test suite**, so every CAP result
is **provisional** until Task #7 passes. The user chose to run the in-built-CAP
investigation first and defer the validation. Nothing here is blocked; it just
needs the regression run to become trustworthy.

## What was changed (inq-study ONLY — `inq/` proven pristine)

`inq-study/src/hamiltonian/self_consistency.hpp` (2 edits):
1. line ~176: `auto vscalar = vion_;` →
   `basis::field<basis::real_space, typename HamiltonianType::potential_type> vscalar(vion_);`
   (complex in real-time `ks_hamiltonian<complex>`, double in GS → GS bit-identical).
2. line ~191: `energy.external(... )` wrapped in `inq::real(...)` (CAP imaginary
   part excluded from real energetics; `energy.hpp` already `real()`s band
   energy/eigenvalues, so no further edit needed).

`shared/config.sh`: made `INQ_SOURCE`/share paths env-overridable (defaults
unchanged) so `INQ_SOURCE=…/inq-study inq-run` builds the fork. Build recipe:
```
INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study \
INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share \
PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod \
inq-run --reconfig          # first build full (~15-20 min); then incremental
```

## WHY this is a repair, not a hack (git archaeology — vindicates the team)

`perturbations::absorbing` is a real INQ feature the team built (Yifan Yao + Xavier
Andrade, Jan 2023): commits `430cd28f` create CAP, `3de075b4` wire into H,
`c390c091` make `update_hamiltonian` complex, `df2466a3` complex potential in
propagation, merged `60266e26`. **A later refactor `bd4a46fe` (2024-04-09, "Use a
field_set potential in update_hamiltonian") reintroduced a real `vscalar`
intermediate** that the perturbation writes into before the (already-complex)
`vks` copy → regressed absorbing. Because **no test drives `absorbing` through
propagation**, the regression passed CI silently. Our `inq/` HEAD (`44f73d95`,
2026-03-26) inherits it. The 1-line fix re-aligns `vscalar` with the team's own
complex `vks` design. Empirical proof it's correct: the depth sweep reproduces the
De Giovannini Fig.4 reflection U-shape (ε min ~1e-5 at η≈1 Ha, reflection turn-up
beyond).

## DEFERRED: Task #7 — engine regression (the gate)

Until this passes, CAP ε are provisional. Recipe (see task #7 for full detail):
1. **Cheap kill-switch first:** one ground-state on pristine `inq/` vs `inq-study`,
   same system → total energy MUST be identical (complexify is inert when
   pert=none). If it diverges, stop.
2. Build `inq-study` (full), run `ctest --output-on-failure --timeout 2000`
   serial + `INQ_EXEC_ENV="mpirun.openmpi -np 4" ctest …`.
3. **Differential:** also run `ctest` on pristine `inq/` as baseline; acceptance =
   zero tests passing on `inq/` but failing on `inq-study`. Watch ground_state
   (h2o, silicon), real_time, and the `perturbations::absorbing` unit test.
4. If green: lift "provisional" everywhere. If red: isolate; repair in `inq-study`
   only, never `inq/`.

## SEPARATE, pre-existing `inq/` deviation (NOT this work — user decision pending)

`inq/src/real_time/viewables.hpp` carries a 2026-05-01 in-place `ham()` accessor
(for jellium `<φ|H|φ>` observables) — a 2nd local mod beyond the allowed CUB
workaround (`reduce.hpp`). A copy is at `ResearchProject/literature/misc/viewables.hpp`.
Recommend migrating it to `inq-study` (keep `inq/` pristine) OR documenting it as a
sanctioned exception in `inq-immutable.md`. Not mine; predates this session.

## Resume checklist
- [ ] Run Task #7 (cheap GS check → full differential ctest).
- [ ] On green, remove "provisional" flags in `absorbing-boundary.md` + the CAP ipynb.
- [ ] Decide `viewables.hpp` (migrate vs sanction).
- [ ] Consider upstreaming the `vscalar` fix to INQ (it's a real upstream regression).
