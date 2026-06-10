# ADR 0006 — Minimum observable set as a declared manifest + post-run validator

Date: 2026-06-10
Status: accepted (design; implemented in phase 3)

## Context

Runs across the project (coronene LEED, jellium WP, jellium classical, free-WP)
each produce many **primary observables** written directly by `run.cpp`. There
is no enforced standard for *which* must be present, so runs silently drift:
`docs/observables/catalogue.md` §4 documents real gaps (WP momentum stats absent
from the σ=5 energy-sweep runs, density_fourier missing when the density phase is
skipped, eigenvalue/occupation inconsistencies). Cross-run comparison and the
bit-identical regression strategy both depend on every run measuring the same
core quantities the same way. We want a *deterministic* way to guarantee this.

## Decision

Define a **layered minimum observable set** (universal core ∪ per-run-type
required, plus optional) and enforce it with a two-part hook:

1. **Pre-run manifest (C++).** `run.cpp` consults a single
   `MinimumObservableSet` table in inqkit and writes
   `results/observables_manifest.json` at startup, declaring its run-type and
   the required∪optional observables (path, schema, optional physical invariant)
   it commits to produce. The run cannot under-declare its type's required set.
2. **Post-run validator (Python).** `inqview.validation.validate_run` reads the
   manifest and checks each observable at four tiers — existence, schema, finite,
   and an optional manifest-declared physical invariant. It is also runnable on
   existing runs to audit the catalogue §4 gaps retroactively.

Spec: `docs/observables/minimum-set-spec.md`.

## Alternatives considered

- **Pure C++ self-validation** (run asserts at exit). Self-contained, but
  physical-sanity checks and re-auditing old runs are awkward in C++.
- **Pure Python validator keyed off a run-type table** (no manifest). Simplest,
  audits all existing runs, but the run isn't self-describing — a run that
  silently changed what it measures wouldn't record the change.

We chose **C++ declares + Python validates** so the *contract* is owned by the
run that makes it (self-describing, travels with `results/`) while the *rich
validation* lives where the numerics/plotting stack already is (and can re-audit
history).

## Consequences

- **Hard to reverse:** once runs emit a manifest and downstream tooling/tests
  assume it, the schema becomes a stable contract (hence `schema_version`).
- A new run-type or observable is a table edit in one place (the
  `MinimumObservableSet`), not scattered across run scripts.
- `scripts/verify_smoke_outputs.py` is subsumed by the tier-4 validator.
- Promotes `energy_hartree`/`energy_xc` to the universal core (fixes the
  `ObservableSelection` default-`false` inconsistency).
- Per rule #6, both halves ship with tests (table mapping + manifest round-trip
  in inqkit; `validate_run` pass/fail cases in inqview).
- **Scope:** this ADR is design only; the migration/cleanup phase-2 work and the
  bit-identical complex-run validation do **not** depend on it landing first.
