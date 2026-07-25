# Plan — Expand & enforce the minimal observable set (primary + derived)

Status: **drafted 2026-06-15** (session "minimal-set-observables"). Spec-only so
far; no code edited. Authoritative set:
`docs/observables/minimum-set-spec.md` §"2026-06-15 — Expanded approved set".

## Motivation

The canonical `inq-stack/include/inqkit/observables/minimum_observable_set.hpp`
(ADR 0006) enforces only a **subset** of the observables the project treats as
minimal, and the **derived (post-processed) layer is enforced nowhere**. Result:
runs silently drift (e.g. `density_wp` saved as optional / coarse cadence;
`current`/`dipole` written but not required; no run is obliged to produce a
stopping-power, loss-function, or energy-balance output). The user walked the
full observable universe and approved an expanded two-layer set per run-type.

Current enforcement reality (verified 2026-06-15):
- Only hooks present: `commit_message_check.py`, `file_placement_check.py`.
- C++ manifest write exists in some `run.cpp` / `run_template.hpp`.
- `inqview.validation.validate_run` is BUILT but wired into **zero** `analyse.py`.
- Skill Phase 7 carries an explicit TODO to link the set via deterministic hooks.

## Scope of changes (not yet done — review gate before each)

### A. PRIMARY layer — encode in the C++ single source of truth
Edit `minimum_observable_set.hpp` (`universal_core()` + per-type blocks):
1. Promote `current_xyz`, `dipole_xyz` to `universal_core()`.
2. Add `delta_density_l2` (csv, `zero_at_t0` invariant); add `step`/`time_au`
   monotonicity invariants.
3. jellium-WP: flip `density_wp_rt` to **required**; add `density_delta_raw`,
   `density_delta_coarse`, `step_delta`, `wp_config`, `wp_injection_report`,
   full O_ij overlap, complex WP wavefunction VTI, `gamma_transitions`.
4. jellium-classical: add `density_delta_raw/coarse`, `step_delta`,
   `occupations_vs_time`, `momentum_distribution`, full O_ij overlap,
   proxy overlap, `gamma_transitions`.
5. coronene: add the three RT densities, WP initial density+WF, windowed +
   instantaneous LEED screens, GS orbital densities, WP-only overlap.
6. free-WP: add density_wp VTI, `wp_momentum_stats`, `wp_real_space_stats`,
   `momentum_distribution`.
- Update the inqkit table unit test (each run-type → expected required names)
  and the Cluster-O drift eval that pins skill Tier tables ↔ `.hpp`.

### B. DERIVED layer — new contract (the genuinely new mechanism)
The header/manifest only models raw observables. Decide & build a derived
contract so `analyse.py` is obliged to emit the approved derived set:
- **Option B1 (recommended):** extend the manifest schema (`schema_version` 2)
  with a `derived` array; `validate_run` checks each derived output's existence
  (file present under `analysis/`) + a light sanity tier. `analyse.py` calls
  `validate_run` at the end and the run is gated on it.
- **Option B2:** a separate `derived_manifest.json` produced by `analyse.py`
  itself, validated independently.
- Either way: **wire `validate_run` into the canonical `analyse.py` templates**
  (the missing enforcement half). One call, run-gating, byte-identical to raw.

### C. Propagate to docs/skill
- `tddft-simulations` skill Phase 3 (Tier 1/2 tables) + Phase 7 (resolve the
  TODO) → reference the expanded set; keep them an operational view, not a 2nd
  source of truth.
- `docs/observables/catalogue.md` §4 gap analysis → note which gaps the new
  required set now forbids.
- ADR 0006 → add a short amendment note (set expanded + derived layer added
  2026-06-15) or a new ADR if the derived contract is a material design change.

## Carry-over items / caveats
- `secondary_electron_yield` (WP derived, **required**) is **not yet
  implemented** — building it is part of this plan, not assumed to exist.
- `density_system` VTI kept for schema parity but is semantically unreliable
  (canonical bath = `total − wp` in post); document, do not trust the raw field.
- `planewave_decomposition` and the all-orbital wavefunction dump were
  **rejected** — do not re-add without revisiting the spec.
- `knudsen_ke` rejected for WP; WP S(v) is taken from `energy_decomposition_vs_z`
  + `momentum_band_free_vs_jellium`. Confirm this holds as the campaign's WP
  stopping-power method before deleting any existing knudsen tooling.

## Validation (rule #6 — both halves ship tests)
- A: extend inqkit `MinimumObservableSet` table test + manifest round-trip.
- B: `validate_run` derived-tier pass/fail cases on synthetic run dirs; a real
  completed jellium-WP run must PASS the expanded set (re-audit, like the
  phase-3 wp_momentum_stats fix).
- Do NOT claim "enforced" until a real run is validated green against the set.
