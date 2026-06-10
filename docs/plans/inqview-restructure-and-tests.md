# Plan: inqview restructure + test suite

Status: **design locked (architecture) — physics-method TODOs still open**
Date opened: 2026-06-10
Driving review input: `docs/code-revitalisation/inqview-todo-catalogue.md`
(29 aggregated TODOs). Decisions recorded in ADRs 0003/0004/0005 and
`CONTEXT.md`.

## Goal

Restructure `inq-stack/python/inqview/` into a simple, community-releasable
library with a clean public API, and stand up a portable test suite — using
**characterization tests to guard the move** (behaviour-preserving), with
bug fixes done as separate deliberate steps afterwards (same discipline as
the inqkit round: E01–E04).

## Locked architecture (ADR 0003)

Four dependency-layered sub-packages; `import inqview.analysis`/`inqview.io`
must pull in **no matplotlib and no VTK** (enforced by a test):

| Package | Role | May import |
|---|---|---|
| `inqview.io` | loaders + field/format dataclasses | numpy only |
| `inqview.analysis` | numeric kernels → frozen dataclasses | numpy/scipy only |
| `inqview.visualisation` | all rendering (mpl, VTK/paraview, GIF) + theme | plotting libs |
| `inqview.pipeline` | thin phase orchestration (compute→plot→write) | the above |

`report1/**`, `scripts/**` = applications (relocate, not tested).

### Data contract (locked)
- Analysis kernel: `compute(...) -> FrozenDataclass` holding numpy arrays +
  metadata (units, axis names, provenance). Mirrors inqkit `Moments`/Result
  structs.
- Renderer: `render_*(result, ...) -> Figure/artefact`. Consumes the
  dataclass; never recomputes.

### Canonical visualisation theme (ADR 0004)
- Promote `report1/_shared_style.py` (+ `report-figures` skill +
  `docs/reports/report1/figures/global_style.md`) into
  `inqview.visualisation.style`. Supersedes generic `config.py`/`defaults.py`.
- **Semantic cmap roles** (no literal cmaps in phases): `sequential→inferno`,
  `diverging→RdBu_r` (zero-centred), `phase→twilight_shifted`.
- **Fixed-dimension factory**: `figure_one_col()`=3.5×3.0 in (fixed axes
  rect), `figure_two_col()`=7.0 in wide. Individual plots only; panels are a
  LaTeX concern. Guard the tight-bbox/constrained-layout pitfall (memory
  `reference_fixed_dimension_plot_pitfalls`).

## Module → package mapping (first pass — refine during migration)

- **io**: `data.py`, `fields.py`, `screens.py`→`io/leed.py` (resolve name
  collision: loader, not phase), `run_summary.py` parse half.
- **analysis**: `fourier.py`, `overlap.py` (compute), `postprocess/`
  {`wake` bath math, `density_fourier` loss fn, `kl_divergence`, `stopping`,
  `lindhard`, `_ifft`, `momentum`, `density` compute, `energy_balance`,
  `bath_energy` compute, `gamma_transitions`, `knudsen_ke`, `spectral_weight`,
  `wp_trajectory`, `state_energies`, `eigenvalues_gs`, `gs_projected_occupations`,
  `occupations`, `compare`}.
- **visualisation**: `plots.py`, `paraview.py` (verify still used), `style.py`
  (NEW), `postprocess/paraview_3d.py`, the colorbar/`shared_clim` helpers from
  `wake.py`, `_common.write_animation`, the plotting half of every phase.
- **pipeline**: `postprocess/pipeline.py` dispatcher, `_common` (ensure_dir/
  need_rebuild/title), `layout.py`, each phase's thinned `run()`.
- **internal/relocated (out of public API)**: `email.py`→`_notify`/scripts;
  `report1/**`, `scripts/**`→applications.
- **legacy verify-then-cut**: `vti.py` Python writer (C++ owns VTI writing);
  keep `convert_real_series_to_vti` only if `paraview.py` pipeline is live.

## Test strategy (ADR 0005) — whole suite pure-tier, portable

Location: `inq-stack/python/tests/` (pytest). No GPU, no INQ, fixtures < 5 MB.

| Layer | Test kind | Expected value |
|---|---|---|
| io parsing | read committed small file → compare | exact (deterministic parse) |
| analysis kernels | reduced/analytic system | analytic, `np.allclose` tolerance |
| integration | free-space-WP committed output | analytic law (σ(t), ⟨p⟩, centroid) |
| visualisation | numeric theme-config only | `figure_one_col()==(3.5,3.0)`, role→cmap, rcParams |
| deps-clean | import-graph assertion | no mpl/vtk in `analysis`/`io` (ADR 0003) |

Plus cheap **physical invariants** layered on golden (bath ≈ N_e, KL ≥ 0,
norm ≈ 1) so fixtures find bugs rather than freeze them.

### First tests to lock (order)
1. `test_deps_clean` — `import inqview.analysis`; assert no `matplotlib`/`vtk`
   in `sys.modules`. (Guards ADR 0003 from day one.)
2. `test_theme` — numeric theme constants (ADR 0004).
3. `test_fourier` — sinusoid → known peak bin.
4. io parsing tests (LEED `.dat`, meta json, VTI read) on hand-built fixtures.
5. Free-space-WP integration fixture + analytic golden.

### Fixtures
- Build a small free-space-WP output once (committed) + analytic golden.
- Hand-built tiny `.dat`/json/VTI for io parsing.
- A trimmed (grid-strided ~32³, few-frame) real WP-run snapshot for the
  pipeline orchestration smoke test — fixture-build mechanism TBD (was being
  decided; user leaning to free-space-WP analytic combo over a heavy real
  snapshot). **OPEN.**

## Method TODOs — resolved + still open

### RESOLVED (2026-06-10, recorded in `docs/validation/inqview-findings.md`)
- **density_fourier / loss function** — IV-M01 (both `axial` + `3d_binned`
  modes) + IV-M04 (renamed `PlasmonSpectrum`, peak-locator semantics).
  Formula independently validated: `|n_q(ω)|²/q²` is a plasmon **peak-locator,
  NOT** −Im[1/ε] (`docs/validation/loss-function-formula-validation.md`).
  Test = undamped-plasmon phasor → exact peak + 1/q² scaling. Bugs IV-E01
  (real-only FFT, fix), IV-E02 (relabel). Memory updated.
- **wake / bath COD** — IV-M02: inqview recomputes COD in Python (node
  convention), does NOT reuse inqkit CSV; cross-check test documents E04 as a
  dx/2 offset. (`total−wp = system` bath identity already in
  `reference_canonical_bath_density`.)
- **Φ-minimum-set** — IV-M03: one global `PHASES_MINIMUM = (summary,
  observables, density)` + per-call extras (user chose flat over per-system
  presets).

### RESOLVED (continued) — ALL method TODOs now closed
- **fourier** — IV-E03: window coherent-gain not corrected (`/n` → `/win.sum()`);
  confirmed by signal-validation agent (`fft-normalization-validation.md`); fix
  + sinusoid test during migration. Core `transform()` → analysis.
- **energy split** — IV-M07: new `energy_components` (kinetic/Hartree/xc/ext
  flow + initial-vs-final bars + ΔE(t) lines + GIF) is primary; band-sum ledger
  kept caveated; orbital-level WP×component deferred (needs saved ψ).
- **Φ-phase-independence** — IV-M08: strictly sequential; parallelism = future.
- **features** — IV-M09: band-structure deferred to multi-k (QKE/QBall) runs;
  `gs_projected_occupations` KEPT (excitation diagnostic) + t=0 identity test.

### FUTURE TODOs (parked, not this round)
- True ε-vs-k-path band structure for multi-k QKE/QBall runs (IV-M09).
- Orbital-level WP×component energy cross-decomposition (IV-M07; needs saved ψ).
- Pipeline process-pool parallelism over pure analysis phases (IV-M08).
- KL frame-to-frame drift-rate + `(k,t)` momentum-carpet contour (IV-M05).

## Cross-library integration tests + validation strategy
- **Free-WP integration tests — separate per library (IV-M11).**
  - inqkit: `test_free_wp_engine.cpp` (engine tier) — non-interacting Gaussian
    WP, assert analytic free-particle laws (σ_r(t), ⟨p⟩=k₀, ballistic centroid,
    norm, energy). First integration test on the inqkit suite.
  - inqview: independently-built free-space-WP fixture; post-process → same
    analytic laws.
- **Validation: code-review at the END (IV-M10)** — user triggers
  `/code-review ultra` once on the whole suite (assistant cannot launch it).
  Expected values are nonetheless derived analytically up front (anti-circularity).
- **FFT drift-removal (IV-M12)** — under agent review
  (`fft-drift-removal-validation.md`); fold the chosen `subtract=` default into
  the fourier kernel + fix the #2 un-subtracted-spectra correctness bug.

## TODO-list cross-reference (todo.txt / todo_later.md, examined 2026-06-10)
- Covered by our decisions: #7 projected_occupation→IV-M09, #7 energy_balance→
  IV-M07, #4 new observables→tested inqkit headers, #5 orthonormalisation→inqkit
  E03 fixed, #1/#2 FFT subtraction→IV-M12.
- Already done/stale: todo_later "MPI slice Allreduce"→inqkit E01 fixed;
  "screen-z FFT bug" + "complex orbital fftshift"→fixed Phase 3; #8 done.
- Out of scope (physics/sim): #3 QKE kicks (note: single-k → narrows IV-M09
  band-structure deferral), #6 WP-revival reading, #9 Gaussian potential;
  todo_later CAP/geometry-relax/E_cut/VESTA (coronene physics).

## Migration sequence (after physics TODOs resolved)
1. Stand up `tests/` + lock the first tests against the CURRENT (pre-move)
   code → green baseline (characterization), analytic expected values.
2. Add the inqkit free-WP engine integration test (IV-M11).
3. Create the four packages; move modules per mapping; split each phase into
   `analysis.compute` + `visualisation.render` + thin `pipeline` run.
4. Keep the suite green through every move. Cut `vti.py`/`email`/verify
   `paraview.py` once their consumers are confirmed.
5. Apply the deferred bug fixes (IV-E01 FFT real-only, IV-E03 coherent-gain,
   IV-M12 subtraction, relabels, energy_components) as separate red→green commits.
6. END: user triggers `/code-review ultra` on the full suite (IV-M10).
