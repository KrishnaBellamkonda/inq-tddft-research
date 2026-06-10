# Plan: inqkit / inqview unit-testing rejuvenation

Source prompt: `docs/prompts/codebase_rejuvination/task_unit_testing.md`
Grilling session: 2026-06-08. Glossary in `CONTEXT.md`; decisions in
`docs/adr/0001-*`, `docs/adr/0002-*`. This plan is the locked design; the
per-component test contents are locked individually during execution.

## Locked decisions (from grilling)

1. **Scope & order.** C++ `inqkit` **first**, then `inqview` core +
   `inqview.postprocess`. One-off scripts (`inqview/report1/**`,
   `inqview/scripts/**`) are **excluded** from testing and will be relocated
   during the final restructure.
2. **Test tiers (2-tier).** `unit` (single component; may link INQ on a tiny
   CPU grid) and `integration` (multi-component / round-trip / GPU). "pure"
   vs "engine-coupled" is a *property* (CI runnability), not a tier.
3. **C++ harness** (ADR 0001). External `inq-stack/tests/cpp/` mirroring the
   include tree; headers untouched; Catch2; one ctest test per file; labels
   `pure` (C++17 only) and `engine` (links INQ).
4. **CI topology** (ADR 0002). Local-first: git pre-push hook + make/ctest/
   pytest run the full suite on the GPU box every change. Hosted GitHub
   Actions = deferred pure-tier backstop. GPU never in cloud CI. Tag tests
   from day one.
5. **Verification agents (two, independent).** formula-validation agent
   (re-derives math from cited source, pre-lock) and test-validation agent
   (audits the written test, pre-suite). Separate fresh-context spawns; B
   never sees A's derivation. See `CONTEXT.md`.
6. **Mapping (subtask 1).** Plugin-first: user installs/uses
   `understand-anything`; the canonical artifact is `docs/inqkit_map.md`,
   assembled from its output.
7. **Ideas file (subtask 1).** `docs/notes/inqkit-rejuvenation-ideas.md`
   (running) → promoted to this plan in subtask 4.
8. **Error log (subtask 3).** `docs/validation/inqkit-errors.md`, structured,
   decision-pending. Source is not edited until a decision is recorded.
9. **Gate & done.** Per-component lock in subtask 3 (propose → accept → lock
   formula → write → audit → run → user "lock" → next). Per-subtask LOCKED
   summary + approval. Definition of done below.

## C++ testing surface (~24 real headers; bottom-up)

| Order | Tier | Headers |
|---|---|---|
| 1 | pure | `detail/grid_layout`, `jellium/shells`, `config/tsubonoya_2014_coronene` |
| 2 | engine | `fields/real_field_3d`, `complex_field_3d`, `orbital`, `density` |
| 3 | engine | `observables/center_of_density`, `density_delta`, `momentum_distribution`, `wp_momentum_stats`, `wp_real_space_stats`, `orbital_overlap`, `eigenvalue_dump`, `state_energy_writer`, `occupations_writer` |
| 4 | engine | `io/real_field_3d_writer`, `complex_field_3d_writer`, `vti_image_data_writer`, `observables_writer` |
| 5 | engine | `wavepacket/wavepacket`, `injection_report` |
| 6 | engine | `real_time/step_context`, `real_time_session` |
| 7 | engine | `screens/plane_screen`, `leed_pattern_accumulator` (in-progress: test only completed behavior) |

**Deferred — 11 empty/TODO placeholders** (never implemented on any branch;
handled in restructuring, not testing): `core/pipeline`, `core/session_context`,
`core/task`, `jellium/analytics`, `config/simulation_config`,
`ground_state/ground_state_tasks`, `detail/validation`, `detail/filesystem`,
`detail/text_io`, `io/manifest_writer`, `io/text_summary_writer`.

## Python phase (after C++) — recommended, confirm at phase start

- Mirror the split with a pytest marker: `@pytest.mark.engine` for tests
  needing INQ-produced data / CUDA / VTK; default = pure (numpy-only).
  Hosted CI later runs `pytest -m "not engine"`.
- Add a `[tool.pytest.ini_options]` block to `inq-stack/pyproject.toml`
  registering the `engine` marker and `testpaths = ["tests/python"]`.
- Same formula-bearing dual-agent treatment for: `postprocess/lindhard`,
  `stopping`, `fourier`/`density_fourier`, `kl_divergence`, `knudsen_ke`,
  `gamma_transitions`, `spectral_weight`.
- Relocate the stray `inqview/postprocess/test_lindhard.py` into `tests/python`.

## Definition of done

- [ ] All ~24 C++ headers have **locked** tests; `ctest` green (pure + engine).
- [ ] inqview core + postprocess at the same bar; `pytest` green.
- [ ] Every failure logged in `docs/validation/inqkit-errors.md` with a
      recorded user decision.
- [ ] Restructuring plan (subtask 4) locked, implemented, re-reviewed via the
      understand plugin (subtask 5).
- [ ] Local CI hook runs the full suite green end-to-end.
- [ ] CLAUDE.md drift resolved (no advertised-but-empty modules).

## Subtask checklist

- [ ] **1 Map** — install/use understand plugin → `docs/inqkit_map.md`;
      accumulate ideas in `docs/notes/inqkit-rejuvenation-ideas.md`.
- [ ] **2 Component→test mapping** — per header: behaviour that matters,
      expected output, failure definition, tier; user locks each plan.
- [ ] **3 Write tests** — per-component loop with both agents; document
      errors, do not edit source to pass.
- [ ] **4 Brainstorm fixes** — interview over errors + ideas → locked
      restructuring plan (promote into this file).
- [ ] **5 Implement** — per-component change + review; modular, extendable.

## Standing constraints

- Never edit a non-empty file unless the locked plan requires it.
- Never edit `runs/` under `ResearchProject/`, `Tutorial/`,
  `QuantumKickExtension/`; re-runs are new runs.
- Never modify source to make a test pass without a reviewed, accepted change.
- Never batch-generate the suite; one component at a time, reviewed.
