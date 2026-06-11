# Eval: Cluster O — minimum observable set canonicalisation (LOCKED 2026-06-11)

Canonical = `inq-stack/include/inqkit/observables/minimum_observable_set.hpp`
(executable; writes `results/observables_manifest.json` at run start). Every
other enumeration of "required observables per run-type" must agree with the
manifest it emits.

## Programmatic drift eval (.claude/evals/programmatic/)

For each `RunType` (coronene / jellium_wp / jellium_classical / free_wp):
1. Emit the manifest from the `.hpp` (build a tiny harness or call
   `manifest_json(type, …)`); parse the required observable names.
2. Assert the **spec doc** `docs/observables/minimum-set-spec.md` lists the same
   required set for that run-type (parse its table).
3. Assert the **catalogue** `scan_runs.py` `FILE_OBS`/`DIR_OBS` flag set is a
   superset of the required names (every required observable is detectable).
4. Structural: assert `tddft-simulations` SKILL Phase 3 does **not** re-enumerate
   the required Tier-1/Tier-2 lists as literal hardcoded names — it must
   reference the min-obs-set (grep for a sentinel marker the de-dup leaves).

PASS = all four agree for all four run-types. Any drift = FAIL (names the
diverging artifact + run-type).

## Functional eval (already exercised last phase — re-confirm)

A run declares its run-type → `write_manifest` emits the manifest → `validate_run`
passes its 4 tiers on a real `results/` dir → catalogue flags match the manifest.
(One jellium-wp run validated PASS last phase; re-run after Phase-3 de-dup.)

## What de-dup must NOT break

- The `.hpp` remains the only place a required observable is *defined*. The spec
  doc and skill prose are *views* of it, checked by the drift eval.
- Optional/extra observables a specific run adds (Phase 3 Tier 3) are allowed and
  not subject to the drift eval — only the required set is canonicalised.
