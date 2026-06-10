# inqview TODO catalogue (review extraction)

Aggregation of the 29 `TODO`/`FIXME` review comments in the inqview Python
package (`inq-stack/python/inqview/`), the backbone for the inqview test +
restructuring planning (mirrors `inqkit-todo-catalogue.md`). Grouped by theme.
Line refs are `inqview/<path>:<line>`.

**Type** key: `Q` understand/answer · `REFAC` restructure · `REDUN` possible
redundancy/removal · `METHOD` physics/method question · `VIZ` visualisation/output
· `IMPORT` import-convention · `FEATURE` new capability.

## Cross-cutting themes
- **Φ-imports** — repeated "is this import convention right / wise?" across
  modules → a package-wide import-style decision. → pipeline:22, kl_divergence:40,
  bath_energy:32, state_energies:27, lindhard import.
- **Φ-minimum-set** — "minimum set of observables/outcomes per run" — which
  phases are essential. → pipeline:8, bath_energy:39, energy_balance.
- **Φ-viz-rule** — a visualisation standard referenced as "TODO 1a/1e/1f"
  (linear+log, no clipping, per-component energy). → _common.py:36 (1a),
  plots.py:40 (1e), overlap.py:74 (1f), observables.py:111 (1e).
- **Φ-redundancy** — duplicate/unclear modules: two `screens.py`
  (screens.py:14 vs postprocess/screens.py), `vti.py` maybe redundant. → vti:4,
  screens:14.
- **Φ-cod-reuse** — wake.py should reuse the inqkit centre-of-density (ties to
  inqkit **E04** half-cell + the bath = total−wp definition). → wake:81, wake:77.

## By module
| Module:line | Type | Essence |
|---|---|---|
| lindhard.py:42 | METHOD/REFAC | analytic calc — better organised elsewhere? |
| fourier.py:8,12 | METHOD/Q | FT windowing/detrending need careful check; what is fourier.py used for now? |
| vti.py:4 | REDUN | is this file redundant / removable? |
| screens.py:14 | REDUN/Q | what does the OTHER postprocess/screens.py do? (dup) |
| plots.py:40 | VIZ | total-energy-only plot (viz rule 1e) |
| energy_balance.py:31 | METHOD | output: jellium bath vs wp orbital energy split |
| gs_projected_occupations.py:55 | VIZ | overlay heatmaps unused → need to decide/remove |
| _common.py:36 | VIZ | matches viz rule 1a |
| state_energies.py:27 | IMPORT | import statement to be corrected? |
| orbitals_per_kpoint.py:30 | FEATURE | can a band-structure plot be made here? |
| overlap.py:74 | VIZ | [0,1] normalisation hides structure (viz rule 1f) |
| pipeline.py:8,11,22 | Q/REFAC | minimum observable set; are phases independent / parallelisable; import wise? |
| kl_divergence.py:40,46,49,53 | IMPORT/METHOD/FEATURE | explain import; KL as time series (cost?); KL as WP-localisation metric; contour viz (Runfeng) |
| wake.py:74,77,81 | REFAC/METHOD | postprocess submodule sense; total−wp=system equivalence; **reuse inqkit COD** |
| density_fourier.py:28,33 | METHOD | loss function should be 3D not 1D-restricted; axial vs 1D loss modes |
| bath_energy.py:32,35,36,39 | IMPORT/VIZ/FEATURE | import convention; ensure (incomplete); jupyter notebooks for analysis; part of minimum set |
| observables.py:111 | VIZ | all_energies_vs_time vs per-component (viz rule 1e) |

## Notes for the inqview test/restructure planning
- Many TODOs are **understanding/restructuring**, not test targets — the inqview
  testable surface is the numeric post-processing (FFT/`fourier.py`, loss
  function/`density_fourier.py`, KL divergence, overlap normalisation, energy
  balance, bath = total−wp). Tier split likely: **pure-numpy** (cheap, CI) vs
  **VTK/data-dependent** (needs a run's output files).
- The **viz rule (1a/1e/1f)** is a user standard to locate + document before
  testing plot code.
- **Φ-cod-reuse + density semantics** connect inqview back to inqkit
  (E02/E04 + bath definition) — coordinate the two.
- Full grilling-based test planning (tiers, fixtures, per-module test ideas) is
  the next session's task — this catalogue is the input.
