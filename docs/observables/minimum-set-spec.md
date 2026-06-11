# Minimum Observable Set — specification (phase-3 design)

Status: **design locked 2026-06-10** (grilling session). Implemented in
**phase 3** (this is the contract phase 3 builds to). Companion: ADR 0006,
`docs/observables/catalogue.md` (the full inventory this formalises),
`docs/observables_reference.md` §13.3 (new-observable rationale).

## Purpose

Make the set of **primary (direct) observables** a run produces *deterministic
and enforced*, so no run silently omits a core measurement (the §4 gap analysis
in the catalogue — WP momentum stats missing from σ=5 runs, etc. — is exactly
what this prevents). Two halves of one hook:

- **pre-run (store):** `run.cpp` writes an **observable manifest** declaring its
  run-type and the required∪optional observables it commits to produce.
- **post-run (validate):** an inqview validator reads the manifest + `results/`
  and checks each declared observable at four tiers.

## Structure (layered)

```
required(run) = UNIVERSAL_CORE  ∪  REQUIRED[run_type]
optional(run) = OPTIONAL[run_type]  (declared, validated-if-present, never fatal)
```

A **set member** is any named primary observable regardless of format (CSV
column, VTI series, `.dat` screen). Each member carries: `name`, `path` (or
CSV file+column), `format`, `cadence`, optional `schema`, optional
`invariant`.

### UNIVERSAL_CORE (every run)

| Observable | Path · column | Format | Invariant (tier 4) |
|---|---|---|---|
| step | `observables.csv:step` | CSV | monotone +WRITE_EVERY |
| time_au | `observables.csv:time_au` | CSV | monotone, t[0]=0 |
| energy_total | `observables.csv:energy_total` | CSV | drift `< 1 mHa` over run |
| energy_kinetic | `observables.csv:energy_kinetic` | CSV | finite |
| **energy_hartree** | `observables.csv:energy_hartree` | CSV | finite *(currently default-OFF — promote to core)* |
| **energy_xc** | `observables.csv:energy_xc` | CSV | finite *(currently default-OFF — promote to core)* |
| density_l2 (system) | `observables.csv:density_l2` | CSV | `density_l2(0)=0` |
| **delta_density_l2** (system) | `observables.csv:delta_density_l2` | CSV | NEW; `=0` at t0 |
| GS eigenvalues | `eigenvalues/eigenvalues.csv` | CSV | n_states rows, finite |
| GS occupations | `eigenvalues/occupations.csv` | CSV | `Σf = N_electrons` |
| GS system density | `vti/density_gs_system/…vti` | VTI | `∫n dV = N` |
| run_summary | `run_summary.txt` | text | `run_completed=true` |

Promoting `energy_hartree`/`energy_xc` to core fixes the inconsistent
`ObservableSelection` defaults (currently `false`).

### REQUIRED per run-type (∪ core)

**jellium-WP** — RT densities ×3 categories, WP stats, momentum dist, state-resolved:

| Observable | Path | Notes |
|---|---|---|
| RT density {system, wp, total} | `vti/density_{system,wp,total}/density_t*.vti` | the three **density categories** (§ below) |
| density_delta {system,wp,total} (raw + coarse) | `vti/density_delta*/…` | Δn vs t0 |
| **δn step-to-step** {system,wp,total} | `vti/density_step_delta*/…` | NEW (inqkit T16): `n(t_{i+1})−n(t_i)` |
| wp_momentum_stats | `wp_momentum_stats.csv` | ⟨p⟩,⟨p²⟩,σ_p,E_kin,norm — schema fixed |
| wp_real_space_stats | `wp_real_space_stats.csv` | ⟨r⟩,σ_r,norm |
| momentum_distribution | `momentum_distribution.csv` | n(\|k\|,t); peak at k₀ at t0 |
| state_energies | `state_energies.csv` | E_i(t) |
| occupations_vs_time | `occupations_vs_time.csv` | frozen (audit) |
| WP-only overlap | `overlap/index.csv`+snapshots | |
| norm_per_state | `norm_per_state.csv` | NEW (built): every \|ψ_i\|²≈1 |

**jellium-classical** — core ∪ {RT density system+total, electron_track, state_energies, occupations, overlap}.
- electron_track: `electron_track.csv` (step,t,x,y,z,vx,vy,vz,fx,fy,fz).

**coronene** — core ∪ {RT density {system,wp,total}, LEED screens (total + time-windowed + instantaneous + screen_config + window_ranges), WP-only overlap, WP config+injection report, GS orbital densities, initial WP density+wavefunction}.

**free-WP** — core ∪ {RT density wp, wp_momentum_stats, wp_real_space_stats, momentum_distribution} — the analytic anchor; invariants are the free-particle laws (σ(t), ⟨p⟩=k₀, centroid=k₀t).

### OPTIONAL (declared, validated-if-present)
Full O_ij overlap matrix; proxy overlap + shells; coarse density; v2 complex WP
wavefunction VTIs (`density_wp`, `wavefunction_wp`); all-orbital dump at t_f;
probability **current-density vector field on PlaneScreens** (NEW, future-todos);
gamma_transitions.

## New observables folded in (from your TODO notes)

1. **delta_density_l2** alongside density_l2 — core, all three categories.
2. **Step-to-step density difference** `n(t_{i+1})−n(t_i)` (inqkit T16) — per category.
3. **Three density categories everywhere** (system/wp/total) for every density
   metric. Requires the WP wavefunction to separate bath = `n_total − |ψ_WP|²`
   (catalogue note L100–103); for runs without saved ψ_WP, the wp/system split
   is declared OPTIONAL not REQUIRED.
4. **Vec3 current/dipole** — stored as a unit (inqkit observables_writer TODO);
   CSV columns unchanged for back-compat.
5. **Probability current density** vector field accumulated on PlaneScreens
   (future-todos.md) — OPTIONAL.
6. **Projected occupation** n_i^GS(t) (§13.3) and **energy-balance** ledger —
   derived, belong to the evaluation set, noted here for cross-reference.

## Observable manifest (pre-run, written by C++)

`results/observables_manifest.json` (written at startup, before propagation):

```json
{
  "run_type": "jellium-wp",
  "schema_version": 1,
  "write_every": 2,
  "n_steps": 190,
  "observables": [
    {"name":"energy_total","required":true,"file":"raw/observables/observables.csv",
     "column":"energy_total","format":"csv","cadence":"step",
     "invariant":{"kind":"drift_max","value_mHa":1.0}},
    {"name":"density_system_rt","required":true,
     "path":"raw/vti/density_system/density_t*.vti","format":"vti","cadence":"write_every"},
    {"name":"wp_momentum_stats","required":true,
     "file":"raw/observables/wp_momentum_stats.csv","format":"csv","cadence":"write_every",
     "schema":["step","time_au","px_mean","py_mean","pz_mean","px2_mean","py2_mean",
               "pz2_mean","sigma_px2","sigma_py2","sigma_pz2","e_kin_ha","norm_check"],
     "invariant":{"kind":"norm_band","col":"norm_check","lo":0.97,"hi":1.03}}
  ]
}
```

The run-type → required-set mapping is a single source of truth in inqkit
(a `MinimumObservableSet` table) that `run.cpp` consults to build the manifest;
the run cannot under-declare its type's required set.

## Post-run validator (inqview)

`inqview.validation.validate_run(run_dir) -> ValidationReport`. Reads the
manifest, walks `results/`, runs the four tiers per observable:

1. **existence** — file/column present.
2. **schema** — columns / array shape / cadence (row count ≈ n_steps/write_every).
3. **finite** — non-empty, no NaN/Inf.
4. **invariant** — the manifest-declared physical check (drift, norm band,
   zero-at-t0, peak-at-k₀, ∫n=N, |cod|≈0). Subsumes
   `scripts/verify_smoke_outputs.py`.

Output: a `ValidationReport` (per-observable per-tier pass/fail) + a one-line
`PASS/FAIL`; exit non-zero on any required-observable tier-1–3 failure or a
declared tier-4 invariant breach. Re-runnable on any existing run (audits the
§4 gaps retroactively — a run with no manifest is validated against the
inferred run-type set with a warning).

## Implementation status (BUILT 2026-06-11, branch observable-set/inq-stack)

- **C++ mechanism** — `inqkit/observables/minimum_observable_set.hpp` (pure,
  std-only): `RunType`, `minimum_set(type)` (the single source of truth, core ∪
  per-type), `manifest_json()` / `write_manifest()`. Pure test (45 assertions).
- **Python validator** — `inqview.validation.validate_run(run_dir)` (deps-clean):
  4 tiers + invariant registry (drift_max, norm_band, zero_at_t0, value_band,
  monotone). Test (6 cases). CLI `python -m inqview.validation <run_dir>`.
- **Verified end-to-end:** the C++ manifest is parsed + validated by the Python
  validator; a real completed jellium-WP run **PASSES** the set (the audit caught
  + fixed a wrong invariant on `wp_momentum_stats.norm_check`).

**Run adoption (rollout — one line per run template, additive ⇒ bit-identical):**
```cpp
#include <inqkit/observables/minimum_observable_set.hpp>
using namespace inqkit::observables;
// at startup, after results/ is created:
write_manifest("results/observables_manifest.json",
               RunType::jellium_wp, Cfg::WRITE_EVERY, Cfg::N_STEPS);
```
Post-run: `python -m inqview.validation <run_dir>` (or call `validate_run` in
analyse.py) gates the run. Writing a JSON file does not touch observables, so a
run that adopts this stays byte-identical.

## Tests (phase-3, rule #6)

- inqkit: `MinimumObservableSet` table unit test (each run-type → expected
  required names); manifest round-trip (write → parse → fields match).
- inqview: `validate_run` on synthetic run dirs — a complete one PASSES; one
  with a dropped required column FAILS tier-1; a NaN-injected column FAILS
  tier-3; a norm=2.0 column FAILS tier-4. Plus a real-run audit smoke test.
