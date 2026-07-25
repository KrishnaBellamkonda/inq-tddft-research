# Handover — Minimal observable set (define + expand)

Task: explore the .claude ecosystem to find what minimal observables every run
is required to produce, then rebuild the set by user accept/reject across all
possible primary + derived observables.

## 2026-06-15 — Set approved (spec + plan written; NO code edited)

### Done
- **Located the single source of truth:**
  `inq-stack/include/inqkit/observables/minimum_observable_set.hpp` (ADR 0006).
  Run writes `results/observables_manifest.json`; `inqview.validation.validate_run`
  is the intended 4-tier checker. Covers PRIMARY observables ONLY.
- **Mapped enforcement reality:** only hooks are `commit_message_check.py` +
  `file_placement_check.py`. C++ manifest write exists in a few run.cpp /
  `run_template.hpp`. `validate_run` wired into **zero** analyse.py. Skill Phase 7
  has an explicit TODO to link the set via deterministic hooks. Derived layer
  enforced nowhere.
- **Drift gaps found:** `density_wp` marked optional in header (skill says
  compulsory equal-cadence); `current_xyz`/`dipole_xyz` written but unrequired;
  spec's `delta_density_l2`/`step`/`time_au`/`density_delta`/overlap unimplemented.
- **Ran the full accept/reject** (8 rounds, universal core → WP → classical →
  coronene → free-WP; primary + derived each). Result recorded as the
  authoritative set.
- **Persisted:**
  - `docs/observables/minimum-set-spec.md` → new section
    "2026-06-15 — Expanded approved set" (the authoritative target).
  - `docs/plans/minimal-observable-set-expansion.md` → implementation plan
    (A: C++ primary header; B: NEW derived contract + wire validate_run into
    analyse.py; C: propagate to skill/catalogue/ADR).

### Key decisions (user)
- Scope: **primary + derived**. Structure: **universal core + per-run-type**.
- Promote `current`/`dipole`/`delta_density_l2`/`step`/`time_au` to core.
- **`density_wp` VTI → REQUIRED** for WP (closes the main drift gap).
- Complex WP wavefunction VTI + `gamma_transitions` → required for WP.
- **Rejected:** `norm_per_state`, all-orbital wavefunction dump, proxy overlap
  (WP), `knudsen_ke` (WP S(v) now from energy-decomposition + momentum-band),
  `density_fourier`, `plasmon_fft`, `planewave_decomposition`.
- Custom: `energy_bookkeeping_vs_time` (time-series, not just t_IFW bar).
- `secondary_electron_yield` accepted but **not yet implemented** → build it.

### Not done (next session)
- No code edited (user chose "spec + plan only").
- Implement plan section A (edit `minimum_observable_set.hpp` + table test),
  then B (derived contract + validate_run wiring), then C (skill/catalogue/ADR).
- Each section is a review gate; verify a real WP run PASSES the expanded set
  before claiming "enforced".

### Authoritative files
- Set: `/local/data/public/skcb2/tddft/docs/observables/minimum-set-spec.md`
- Plan: `/local/data/public/skcb2/tddft/docs/plans/minimal-observable-set-expansion.md`
- Source of truth (to edit): `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/observables/minimum_observable_set.hpp`
