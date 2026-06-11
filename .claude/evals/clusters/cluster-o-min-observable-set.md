# Eval: Cluster O — minimum observable set canonicalisation (LOCKED 2026-06-11)

Canonical = `inq-stack/include/inqkit/observables/minimum_observable_set.hpp`
(executable; writes `results/observables_manifest.json` at run start). Every
other enumeration of "required observables per run-type" must agree with the
manifest it emits.

## Programmatic drift eval (.claude/evals/programmatic/run_cluster_o_drift_eval.py)

Implemented (regex over the `.hpp` + spec + skill — no compile needed):
1. Extract the canonical observable names from the `.hpp` (`csv/vti/text("…")`).
2. Assert every canonical name is **covered in the spec doc**
   `docs/observables/minimum-set-spec.md` — verbatim, else via the documented
   key→prose **bridge** in the runner (the spec currently names some observables
   in prose, e.g. `gs_eigenvalues` ↔ "GS eigenvalues"). A canonical observable
   with neither a verbatim hit nor a bridge entry is real drift (code added an
   observable the spec never documented).
3. Structural: `tddft-simulations` Phase 3 carries the canonical-reference
   **sentinel** (`min-obs-set: canonical = …`) — it defers to the `.hpp`, not a
   second source of truth.

Status 2026-06-11: **PASS** (19 canonical observables all covered; sentinel
present). The drift check first ran red — it correctly caught that the spec
names 6 observables in prose rather than by canonical key; bridged + documented.

### Remaining Cluster-O cleanup (precise)

- **Align the spec doc keys to the `.hpp` canonical names** (`gs_eigenvalues`,
  `density_total_rt`, `density_wp_rt`, `leed_screen_config`, `gs_system_density`)
  so the bridge can be deleted and the check becomes pure key-equality.
- Optional: a manifest-emitting harness for per-run-type required-set equality
  (currently the check is name-coverage, not per-run-type set equality) and a
  catalogue (`scan_runs.py`) coverage check.

## Functional eval (already exercised last phase — re-confirm)

A run declares its run-type → `write_manifest` emits the manifest → `validate_run`
passes its 4 tiers on a real `results/` dir → catalogue flags match the manifest.
(One jellium-wp run validated PASS last phase; re-run after Phase-3 de-dup.)

## What de-dup must NOT break

- The `.hpp` remains the only place a required observable is *defined*. The spec
  doc and skill prose are *views* of it, checked by the drift eval.
- Optional/extra observables a specific run adds (Phase 3 Tier 3) are allowed and
  not subject to the drift eval — only the required set is canonicalised.
