# CONTEXT — Glossary for the inq-stack unit-testing rejuvenation

> A glossary only. No implementation details, no plans, no decisions —
> those live in `docs/plans/`, `docs/handovers/`, and `docs/adr/`.
> Terms are resolved during the grilling session for `task_unit_testing.md`.

## Components

- **inqkit** — the project-local C++ header library on top of INQ, at
  `inq-stack/include/inqkit/`. Header-only, C++17, engine-coupled (most
  headers take INQ `basis`/`field` objects).
- **inqview** — the Python post-processing and visualisation package at
  `inq-stack/python/inqview/`. Being restructured (2026-06-10, ADR 0003)
  into four deps-layered sub-packages:
  - **inqview.io** — loaders + field/format structures; numpy only.
  - **inqview.analysis** — numeric post-processing kernels; imports only
    numpy/scipy and returns plain dataclasses/arrays. **Never** imports
    matplotlib or VTK (the deps-clean invariant: a headless cluster node
    can compute observables without plotting deps).
  - **inqview.visualisation** — all rendering: matplotlib, VTK/paraview,
    GIF. The *only* layer allowed to import plotting/VTK libraries.
  - **inqview.pipeline** — thin phase orchestration (compute → plot →
    write artefact); calls analysis + visualisation, holds no math.
  (Supersedes the earlier **core** / **post-processing** two-way split.)
- **leed loader vs screens phase** — NOT duplicates despite the shared
  name. `inqview/screens.py` (→ `inqview.io`) is the `LeedPattern`
  dataclass + `load_leed_pattern` loader; `inqview/postprocess/screens.py`
  (→ `inqview.pipeline`) is the LEED render phase that imports that loader.
  The naming collision is a usability smell to fix in the restructure,
  not a redundancy to delete. Same primitive-vs-phase pairing for
  `overlap.py`.
- **one-off scripts** — single-use presentation/figure code:
  `inqview/report1/**` and `inqview/scripts/**`. Out of scope for testing;
  to be relocated during the final restructuring step, not tested.

## Visualisation standard (2026-06-10, ADR 0004)

- **canonical theme** — the single library-wide plotting standard, promoted
  from `report1/_shared_style.py` + the `report-figures` skill +
  `docs/reports/report1/figures/global_style.md` into
  `inqview.visualisation.style`. Supersedes the generic `config.py`/
  `defaults.py` values (cividis/coolwarm/(6.4,4.2) → the designed standard).
- **semantic cmap role** — phases request a *role*, never a literal cmap:
  `sequential → inferno`, `diverging → RdBu_r` (zero-centred),
  `phase → twilight_shifted`. Keeps every phase visually consistent and
  makes the cmap choice testable/centralised.
- **fixed-dimension figure factory** — `figure_one_col()` (3.5×3.0 in, with
  a fixed axes rectangle so every one-column panel shares an identical data
  box) and `figure_two_col()` (7.0 in wide). The library emits **individual
  plots only**; panel composition is a downstream LaTeX concern. Known
  fragility (see memory `reference_fixed_dimension_plot_pitfalls`):
  tight-bbox + constrained-layout silently break the fixed width — the
  geometry test must guard this.

## Observables (2026-06-10, ADR 0006)

- **primary (direct) observable** — a quantity written *directly by the
  simulation* (`run.cpp` / inqkit writers): energy components, current,
  dipole, the density fields, WP stats, momentum distribution, occupations,
  eigenvalues, LEED screens, the projectile track. Contrast **derived
  observable** — anything computed *afterwards* by inqview from primary ones
  (spectra, COD, stopping power, LEED IFFT, KL, loss function).
- **minimum observable set** — the standardised set of primary observables a
  run is required to produce. **Layered**: a small **universal core** every
  run must emit, plus a **per-run-type required set** (coronene adds LEED
  screens; WP runs add momentum/real-space stats + momentum distribution;
  classical adds the projectile track), plus **optional** extras. A member is
  any named primary observable regardless of format (CSV column / VTI series /
  `.dat` screen). Distinct from the **minimum *evaluation* set** (the required
  *derived* artefacts of `analyse.py`, catalogue §2.1).
- **run-type** — the classification a run declares (coronene / jellium-WP /
  jellium-classical / free-WP) that selects which required set applies.
- **observable manifest** — a machine-readable file the run writes at startup
  declaring its run-type and the required∪optional observables it commits to
  produce (with expected path + schema + optional physical invariant). The
  pre-run half of the enforcement hook; the contract a run is held to.
- **density category** — every density-derived observable exists in three
  flavours: **system** (target electrons, WP excluded), **wp** (the injected
  packet alone), **total** (system + wp). In current jellium raw output
  `density_total` and `density_system` are byte-identical (both WP-included);
  the bath/system-only field is reconstructed as `n_total − |ψ_WP|²`.
- **validation tier** — the post-run validator checks each declared observable
  at four tiers: (1) **existence**, (2) **schema** (columns/shape/cadence),
  (3) **finite** (non-empty, no NaN/Inf), (4) **physical invariant** (optional,
  manifest-declared — e.g. norm∈[0.97,1.03], density_l2(0)=0, |cod_x|≈0).
  Tiers 1–3 mandatory; tier 4 per-observable opt-in.

## Test tiers (this project — 2-tier scheme)

- **unit test** — exercises a **single** component (one inqkit
  function/struct, or one inqview function). May link the INQ engine on a
  minimal hand-built CPU grid. Asserted against a **known/analytic expected
  value** defined and accepted *before* the test is written — never
  retrofitted to current code output. No multi-component chaining.
- **integration test** — exercises **multiple** components assembled
  end-to-end: chained inqkit modules, a C++→Python round-trip (write then
  read back), or anything requiring a real ground-state / GPU run. A
  separate, explicit stage begun only once the relevant units are locked.

### inqview suite portability (2026-06-10, ADR 0005)

- The **whole inqview suite is pure-tier** — pure-Python/numpy on committed
  fixtures, **no GPU, no INQ engine, no multi-GB data**. The engine/GPU is
  used once *by us* to generate fixtures, never re-run in the suite.
- **portable** — a test must pass on any user's machine. Bit-exact golden is
  banned (BLAS/FFT reorder floats across backends); golden is compared with
  physical tolerances (`np.allclose`). Exact comparison is allowed **only**
  for deterministic I/O parsing.
- **free-space-WP fixture** — a free-space wave-packet propagation output
  (committed, small) is the integration anchor: its physics is analytic
  (σ(t) spreading, ⟨p⟩ conserved, centroid = k₀·t/m), so golden is paired
  with a known law, not raw capture. Renderers are **not** tested; only a
  numeric theme-config test guards the designed proportions.

### Property of a unit test (not a tier)

- **pure** — needs no INQ engine; compiles with a C++17 compiler alone
  (or pure Python). Runs on any CI runner.
- **engine-coupled** — links the INQ engine (or imports a CUDA/VTK-backed
  Python path); needs the INQ build (CPU or GPU). Determines where in CI it
  can run.

## Verification agents (independent, fresh-context)

- **formula-bearing component** — a function whose output is defined by a
  mathematical formula with a citable source (e.g. `center_of_density` =
  ∫r·n/∫n, `momentum_distribution` = |ψ̃(k)|², `jellium/shells` magic
  numbers, `lindhard`, `stopping`). Triggers the formula-validation agent.
- **formula-validation agent** — a fresh-context subagent given ONLY the
  formula-as-implemented and its cited source. Re-derives or sanity-checks
  the math independently of the code and of the main session. The formula is
  **locked** only when this agent and the user independently agree.
- **test-validation agent** — a separate fresh-context subagent given the
  written test and the already-locked expected value. Audits the test for
  circularity (asserting code output vs the verified value), correct
  tolerances/units, and isolation. Runs before the test enters the suite.
  Never sees the formula-validation agent's derivation — independence is the
  point.

## Test lifecycle states

- **under-review** — proposed candidate; expected result and (for
  formula-bearing functions) the verified formula are still being agreed.
  Not in the suite.
- **locked** — user-accepted case, classification fixed, expected result
  accepted. Only locked tests enter the suite.
