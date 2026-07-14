# CONTEXT — Glossary for the inq-stack unit-testing rejuvenation

> A glossary only. No implementation details, no plans, no decisions —
> those live in `docs/plans/`, `docs/handovers/`, and `docs/adr/`.
> Terms are resolved during grilling sessions (originally `task_unit_testing.md`;
> the absorbing-boundary vocabulary was added for
> `docs/campaigns/absorbing_boundary/benchmarks_and_parameter_search.md`).

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
- **scientific-figures skill** (2026-06-25) — the figure RULE SET layered on the
  canonical theme; single source of truth referenced by `report-figures` and
  `notebook-making`. Holds the hard rules: density-based system-design plots,
  minimal legends, colorbar-outside-same-height, titles present-for-presentation
  /cropped-for-report, captions in the slide spec, `.drawio`+matplotlib workflows,
  header-only table colour.
- **"linear response"** (2026-06-25, user) — the project-canonical DISPLAYED label
  for the Lindhard / RPA reference curve on every figure and in prose. The code
  module name is unchanged (`inqview.analysis.lindhard_elf` / `pipeline.lindhard`);
  only the rendered label is "linear response", never "Lindhard".
- **system-design plot** (2026-06-25, user) — a figure proving simulation geometry
  is correctly placed, built from the run's REAL density (total-density xz slice
  preferred, or Δn) with dashed slab-extent and CAP-extent lines. NOT a hand-drawn
  cartoon. Run parameters (N, dx, r_s, v…) go on the slide, never in the figure.

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

## Absorbing boundaries (2026-06-13)

> Vocabulary for the De Giovannini–Larsen–Rubio (2014, arXiv:1409.1689)
> reflection-error study, resolved while grilling
> `docs/campaigns/absorbing_boundary/benchmarks_and_parameter_search.md`.

- **absorber / absorbing boundary (AB)** — a device at the box edge that damps
  the outgoing wavepacket to emulate an open boundary. Three families in the
  source paper: **mask function absorber (MFA)**, **complex absorbing potential
  (CAP)**, **smooth exterior complex scaling (SES)**. This work treats the MFA
  (implementation + parameter study) and the CAP (feasibility analysis only).
- **mask function absorber (MFA)** — an absorber applied as a real multiplicative
  mask `M(x) ∈ [0,1]` once per step: `ψ(t+dt) = M(x)·U(t+dt,t)·ψ(t)` (paper
  Eq. 12). The canonical mask is the **sin² mask** (paper Eq. 13): `M(x)=1` for
  `x<0`, `1 − sin²(xπ/2L)` for `0≤x≤L`.
- **complex absorbing potential (CAP)** — an absorber realised as a non-Hermitian
  imaginary potential added to H over `[0,L]`, making `U` non-unitary. Subject of
  the Task-1 feasibility analysis only (not implemented this run).
- **box geometry** — propagation domain `[−X, L]` with `X = 6σ`; the **inner
  region** is `[−X, 0]` (absorber-free), the **absorber region** is `[0, L]`. The
  WP starts at `x₀ = −3σ` (halfway between the left wall and the absorber).
- **Gaussian WP convention** — `σ = 4√2/k₀`, so the momentum width is `k₀/8` and
  the spreading velocity `k₀/8 ≪ k₀` (packet stays coherent). The inqkit
  `WavePacket.sigma` argument *is* this σ (density std = σ/√2).
- **propagation time τ** — `τ = 2(3σ+L)/k₀`, one full round trip; the run stops
  before the reflected wave bounces off the left wall at `−X`.
- **reflection error ε(E,L)** — the survival/reflection measure: the WP norm
  remaining in the **inner region** at `t=τ`, `ε = ∫_{[−X,0]} |ψ_WP(τ)|² dV`
  (paper Eq. 7). Faithful because the ideal free packet has exited the box by τ,
  so the inner region holds only the reflected wave. Needs **no** per-point free
  reference run. `E = k₀²/2` Hartree `= 13.6 k₀²` eV (the paper's `E = 5k₀²/4`
  prefactor is a horizontal-axis relabel only).
- **anchor run** — a no-absorber (hard-wall) propagation scattered across a few
  energies; it should yield `ε ≈ 1` (full reflection), pinning the
  high-reflection asymptote and validating the pipeline end-to-end.
- **two-sided absorber** — an absorber placed on BOTH z-boundaries, total width L
  **split L/2 at each end**. A forward-moving packet meets only the far half, so a
  two-sided L absorbs the beam like a single-sided L/2; the near half catches
  back-scattered / k<0 flux. (CAP form = sum of two `absorbing` slabs; mask form =
  a symmetric two-sided `MaskAbsorber`.) Contrast the earlier **single-sided**
  sweeps. Plan: `docs/plans/twosided-cap-vs-mask.md`.
- **comfortable region** — the (width L*, CAP depth η*) operating point that holds
  the reflection error ε low across the target energy band; the deliverable of the
  two-sided study. The **user** reads it off the ε(E,L) maps — never auto-declared.
- **survival fraction (vs reflectivity)** — for a quasi-monochromatic packet
  (`σ ∝ 1/k₀`, ~12 % spread) ε is true beam **reflectivity**. For a sharply
  localised production packet (`σ = 0.5 Bohr`) the packet disperses to tens–hundreds
  of Bohr before reaching the edge, so ε is instead the **net surviving un-absorbed
  fraction** of a broad dispersing packet — an absorber-applied-to-workload metric,
  not monochromatic reflectivity. The σ=0.5 baselines are a deferred successor task.
- **quasi-1D emulation** — the 1D paper problem run in INQ's 3D engine with a
  **minimal transverse box**. Justified by free-particle factorization
  (`ψ = φ_x·φ_y·φ_z`; the mask depends only on x; transverse norm = 1), so the
  3D inner-region integral equals the 1D ε exactly and transverse extent does
  not affect ε.
- **mask mechanism (in-callback mutation)** — the MFA is applied entirely in the
  inq-stack wrapper by multiplying `electrons.kpin()` by `M(x)` inside the
  per-step `real_time::propagate` callback (the callback captures the outer
  non-const `electrons`; INQ's `viewables` observer is const and cannot be used).
  inq/ and inq-study stay byte-identical for Task 2.

- **showcase run** — a small hand-picked subset (~6) of the ε(E,L) sweep that, in
  addition to the minimum set, writes density VTI frames (all variants) and the
  final WP wavefunction (complex field) for deep post-analysis. Contrast the bulk
  runs, which write only the reduced-cadence CSVs + `epsilon.txt`.
- **excitation metric (vacuum)** — in the vacuum free-WP study there is no Fermi
  sea to excite, so "KS occupation/excitation" is realised as the WP's own
  momentum-space redistribution: `|ψ̃_WP(k)|²` at t=0 vs t=τ (which k-components
  the absorber removed/reflected) plus the surviving WP occupation. The vacuum
  analogue of the jellium GS-decomposition metric.

## System folder structure (ADR 0007)

- **run-set** — a group of runs produced together for one purpose (e.g. the 72
  MFA runs that make the ε(E,L) reflectivity curves). The unit a `hypotheses/`
  subfolder analyses.
- **`shared_gs/`** — the converged ground state(s) reused across all of a
  system's runs (so SCF is not re-run per job). New unifying name; legacy systems
  keep `checkpoints/` (jellium) / `save_gs/` (coronene).
- **`shared/`** — shared config headers / `Common_`-derived cfg structs.
- **`scripts/`** — how runs are *produced*: the build-once binary, the GPU
  dispatcher, `gpu_probe`, the per-run `analyse.py` template.
- **`<sweep_name>/<run_name>/`** — runs are **grouped by sweep** (ADR 0007
  amendment 2026-06-15), one folder per sweep holding its run subdirectories;
  outputs only. Supersedes the original flat top-level `run_*`. The `<sweep_name>`
  matches the analysis folder `hypotheses/<sweep_name>/` exactly. (jellium/coronene
  stay grandfathered-flat; vacuum is migrated as the reference instance.)
- **`hypotheses/<sweep_name>/`** — what a run-set *means*: combined CSVs,
  aggregation/plotting scripts, the study `.ipynb`, `README.md` + figures, and a
  `tests/` subfolder for task-specific implementation/mechanism checks (distinct
  from library-generic `inqkit` unit tests, which live in `inq-stack/tests/`).
  Bare sweep names (no `NN_` prefix) so runs and analysis share one `<sweep_name>`.
- **Study notebook** — the executed `<sweep>_study.ipynb` that narrates a run-set
  (or a *consequential* individual run that has a `docs/plans/` entry): context →
  formulas (every term defined) → fully-reconstructable setup → linked source
  files → results → takeaway. Generated (never hand-edited) by
  `build_<sweep>_report.py` and **auto-rebuilt by the run machinery** (dispatcher
  tail, once per batch; planned single run via its `analyse.py`), not by a
  Claude-Code hook. Authoring contract: the `notebook-making` skill.
- **run-notebook** — the deep **single-run** analysis notebook: the full
  standardised plot battery for ONE run, for analysing it individually in depth.
  Sibling of the **Study notebook** (which narrates what a run-SET *means*); the
  run-notebook narrates what *one run shows*. Reuses the house-narrative context
  (title+question, fully-reconstructable config, linked source files) and
  specialises the results into the **standard battery** — density matrix
  (carpets + one lead GIF), energetics (system components + WP-orbital total
  energy), projectile/transport (COD, stopping, current), collective response
  (dipole→plasmon, E-field), momentum (1D n(k) before/after; 2D k_z–k_⊥ scattering
  is a *future observable*), KS excitation (`gs_projected_occupations`,
  `gamma_transitions`, eigenvalues), loss function (with a low-resolution note),
  and integrity (WP norm, N(t)). **Auto-gated**: a section appears only when its
  observable is present (so it adapts to WP / classical / baseline runs). An
  **assembler over `inqview.pipeline`** (it runs the phases, then embeds their
  figures), NOT a reimplementation. Authoring contract: the **`run-notebook`
  skill** (skill-local, shippable builder); the generated `.ipynb` lives in the
  run's `hypotheses/` folder per ADR 0007.
- **CAP thin-absorber tuning (`cap_thin_L5`, 2026-06-15)** — a vacuum sweep using
  the **in-built inq CAP** (`perturbations::absorbing`, not the MFA) at fixed thin
  width `L=5` Bohr, sweeping depth `η ∈ {−0.01,−0.05,−0.30}` Ha across 11 energies
  (1–100 eV) to find the ε(E) curve whose minimum sits near 10 eV. Plan:
  `docs/plans/cap-thin-absorber-tuning.md`. ε provisional until Task #7.

## CAP in jellium baselines (2026-06-17)

> Vocabulary for `docs/campaigns/cap_in_jellium/baseline_runs.md` (the first use of
> the inq-study CAP in an *interacting* jellium bath). All absorption numbers stay
> PROVISIONAL until the inq-study engine regression (Task #7) passes.

- **bath drainage** — the continuous absorption of *equilibrium* electron-gas
  density by the CAP. Because the CAP is absent from the (Hermitian) ground state
  and switches on abruptly at t=0, the filled Fermi sea inside the absorbing slabs
  is no longer stationary and is steadily removed; current flows from the free
  region toward the slabs. Characterising this (NOT making it negligible) is the
  point of Baseline 1.
- **plane / flux screen (jellium)** — a constant-z (xy) monitor plane emitting the
  **planar density** ρ(x,y;z,t) AND the **integrated z-flux** ∮J_z dA. Built on
  `inqkit::screens::PlaneScreen` (density) + a new flux reducer over the
  `observables::current_density` field. **NOT** a LEED screen — there is no
  diffraction physics here; the LEED accumulator (coronene far-field FFT) is the
  wrong primitive and the term "LEED screen" is retired for jellium.
- **continuity / CAP-sink check** — with a CAP the continuity equation gains a sink
  ∂ₜn + ∇·J = −2η·sin²·n. Comparing dN_free/dt (from density) against the boundary
  flux ∮J·dA (from the CAP-edge flux screens) isolates the CAP-removed term as the
  residual — the quantitative backbone tying density, current, and screens together.
- **Baselines 0–3** — 0: the plain Hermitian GS as the t=0 reference (CAP-free,
  reused checkpoint); 1: CAP on, no projectile (bath-drainage reference); 2: CAP +
  classical σ=0.5 electron at 100 eV; 3: CAP + σ=0.5 Gaussian WP at 100 eV. Runs
  1–3 share one 140-a.u. window so 1 is the exact subtraction reference for 2–3.
- **transit window** — t ∈ [0, t⋆] with t⋆ = (z_edge − z₀)/v₀ ≈ 10.3 a.u. at 100 eV:
  launch to the projectile (at its initial speed) reaching the far free-region edge
  z=+15. The only window where the projectile is cleanly pre-absorber; all clean
  stopping read-off happens here.
- **stopping observable** — the projectile's *mechanical* (drift) kinetic-energy
  loss ΔKE over the transit, S = ΔKE/x. Classical run: −dKE/dz from the Ehrenfest
  ion velocity (clean, unambiguous). Quantum WP: the **coherent-peak** momentum of
  n_wp(k,t), NOT the second moment ½⟨k²⟩ (the second moment is inflated by the
  scattering tail + coarse k-binning → non-monotonic, unusable). **NOT** the
  total-energy difference E_baseline−E_run (during the Hermitian transit, stopping
  energy is *redistributed* into bath excitations that stay in the total, so it
  doesn't appear there; the residual is dominated by differential CAP drainage).
- **projectile edge-zeroing (REQUIRED for classical runs)** — the classical
  Ehrenfest projectile has a *long-range* erf-Coulomb tail `C·erf(r/(σ_pot√2))/r`
  and is never absorbed by the CAP, so after traversing it **wraps the periodic box
  and re-enters the slab** (seen in `p2_classical`), contaminating `E_total(t_f)`.
  Fix: **zero the ion's radial potential once it is about to leave the box** (`z_ion
  ≥ z_edge`), so the late-time state is a clean relaxing bath with no projectile.
  Neutrality-safe because the projectile UPF has `z_valence=0` (no valence charge).
  run.cpp-only mutation; never edit `inq/`. Procedure: `tddft-simulations` skill §2d″.
  **Not yet implemented** — existing classical runs lack it. (Distinct from the
  `p3_classical` mid-slab *trapping* anomaly, which edge-zeroing does NOT address.)
- **pseudopotential (ghost) projectile** — the classical projectile represented as an
  INQ *ion* carrying a UPF whose local potential is the erf-Coulomb `erf(r/(σ_WP/√2·√2))/r`,
  `z_valence=0`. Because that potential is *long-range* (≈1/r) and INQ, seeing
  `z_valence=0`, places the whole tail as a truncated short-range local potential, the
  projectile↔background diagnostic aliases and diverges with the UPF radial cutoff r_cut
  (linearly, unbounded, past the slab). The ledger residual it yields is r_cut- and
  grid-contaminated (see `reference_ghost_upf_tail_aliasing`). Contrast the **perturbation
  projectile**.
- **perturbation (Gaussian-charge) projectile** — the classical projectile represented NOT
  as an ion but as a stationary Gaussian CHARGE added to the KS potential via its Poisson
  potential (`gaussian_projectile_perturbation`, composed with the background by
  `perturbations::sum`). No UPF, no r_cut, no aliasing: the projectile↔electron term and the
  projectile↔background diagnostic `U_proj_bg = −∫n_proj·φ₊` are computed in one Poisson
  convention. This is the *accurate* representation — the residual `d(E_H+E_ext) − U_proj_bg`
  comes out clean and equals the wavepacket Hartree self-energy (≈20.8 eV in-cell p2, r=12,
  σ=0.5). The `pseudopotential (ghost) projectile`'s r_cut-dependent value was the artifact.
- **U_proj_bg** — the classical projectile↔positive-background Coulomb energy, absent from
  INQ's 8-term total (the projectile is not an INQ ion coupled to the background by Ewald).
  Tracked as a diagnostic: the accurate value is the r_cut-invariant "ideal" `−∫n_proj·φ₊`.
- **ledger residual** — `d(E_H+E_ext) − U_proj_bg` (WP minus classical, projectile-background
  added back). With the perturbation projectile it equals the WP Hartree self-energy; combined
  with dXC (self-exchange-correlation) it yields the LDA **self-interaction error** — the only
  genuinely unaccounted energy after localisation (dKin = 3/(4σ_WP²)).
- **quantum component of stopping** — S_WP − S_classical at *matched* σ and v: the
  genuinely quantum effects the classical particle can't have — projectile
  zero-point momentum spread Δp=1/2σ, Pauli exchange with the (identical) bath
  electrons, diffraction. Spreading biases S_WP *low* (charge dilutes, form factor
  e^{−q²σ²/2} cuts coupling), so the full-transit S_WP−S_classical is a *lower*
  bound, NOT the upper bound originally sought. The dispersion-free estimate is the
  **t→0 initial-rate difference** (both projectiles share the rigid σ shape before
  spreading acts; zero-point + exchange are already present at t=0).
- **projectile "width" (user term — ALWAYS the wavepacket σ).** When the user says the
  **width** of a projectile (e.g. "width 0.5 Bohr"), it means the **wavepacket σ_WP**:
  ψ ∝ exp(−r²/2σ_WP²). The classical twin's Gaussian charge is **auto-adapted** so its
  *effective* width is identical: charge std = σ_WP/√2 (so "width 0.5" ⇒ classical
  Gaussian charge std = 0.354 **automatically**). Operationally: ALWAYS build the
  classical UPF via `inqview.io.gaussian_psp.generate_gaussian_psp(sigma_wp=<width>)`
  (it sets charge std = width/√2), so the classical and WP projectiles present the
  **identical** cloud exp(−r²/σ_WP²) to the bath. **Never** hand-pick or round the
  classical charge std. See the √2 validation + exact-match mandate below.
- **σ-convention UNIFICATION (2026-06-21 — wavepacket is the single source of truth).**
  `σ` now means the **wavepacket σ** everywhere: ψ∝exp(−r²/2σ²) ⇒ charge/density std
  = σ/√2. `inqview.io.gaussian_psp.generate_gaussian_psp(sigma_wp=…)` was CHANGED to
  take this unified σ and build its erf charge at std σ/√2, so a classical projectile
  and a WP given the SAME σ present the **identical** cloud exp(−r²/σ²) to the bath.
  (Code: `sigma_wp`/`sigma_charge` fields on `GaussianPspResult`; tests updated +
  passing.) New unified UPFs are named `electron_gaussian_wpsigmaXpY.upf` (e.g.
  `…wpsigma3p0.upf` = σ_wp 3, charge std 2.121 — the matched companion for the WP σ=3
  run; `…wpsigma0p5.upf` = σ_wp 0.5, charge std 0.35355 — the exact-matched companion
  for the WP σ=0.5 slab runs).
  - **√2 VALIDATED + EXACT-MATCH MANDATE (2026-06-23, independent agent).** Confirmed
    from `wavepacket.hpp:234,254` + `gaussian_psp.py:126`: WP **density std = σ_WP/√2**
    exactly (0.5→0.35355). **The localised-jellium slab pair is a ~1 % mismatch, not
    exactly matched:** the WP run used σ_WP=0.5 (density std 0.354) but the classical
    run loaded the **legacy** `electron_gaussian_sigma0p35.upf` whose charge std is
    **0.350** (rounded). Mandate for ALL future quantum-vs-classical runs: generate the
    classical UPF with `generate_gaussian_psp(sigma_wp=…)` so charge std = σ_WP/√2
    **exactly** — never reuse a rounded legacy `sigmaXpY` file for a matched pair. The
    slab plots/notebooks were relabelled to the actual widths (0.354 WP vs 0.350
    classical) and compare to **point-charge Lindhard only** (the σ was chosen from the
    `06_sigma_convergence` sweep to sit in the linear-response regime; measured
    0.706 ≈ point-charge 0.719 eV/Bohr).
  - **LEGACY registry (old convention: filename/run σ = CHARGE STD = unified σ_wp/√2).**
    These predate the unification and are LEFT AS-IS (option (a)); read their σ as
    charge std and multiply by √2 for the unified σ_wp:

    | legacy file / run | charge std (old σ) | unified σ_wp |
    |---|---|---|
    | `electron_gaussian_sigma0p15.upf` + `run_…sv_sigma0p15` | 0.15 | 0.212 |
    | `…sigma0p25.upf` + `…sv_sigma0p25` | 0.25 | 0.354 |
    | `…sigma0p35.upf` + `…sv_sigma0p35` | 0.35 | 0.495 |
    | `…sigma0p4.upf` | 0.40 | 0.566 |
    | `…sigma0p5.upf` + `run_sv_sigma0p5`, **cap_baselines B2** | 0.50 | 0.707 |
    | `…sigma3p0.upf` + `…sv_sigma3p0` | 3.0 | 4.243 |
    | graphene `…sigma1p47*.upf` | 1.47 | 2.079 |

    The `06_sigma_convergence` S(v) study keeps its charge-std labels (no WP there);
    treat its σ as charge std. **WP runs** (`run_wp_*_sigma{0p5,1,3,8}`) already use the
    WP σ ⇒ their labels ARE unified σ_wp, unchanged.
  - **⚠ cap_baselines B2-vs-B3 is NOT width-matched.** B2 classical used charge std 0.5
    (= unified σ_wp 0.707); B3 WP used σ_wp 0.5 (density std 0.354). The notebook labels
    both "σ=0.5" but they are √2-different clouds — the B2/B3 wake/stopping comparison
    conflates a width difference with the classical-vs-quantum difference. Flagged in
    the notebook + handover; do not read it as a matched pair.
- **quantum S(E) sweep (`qsp_phase5`)** — the localised-jellium WP stopping power
  S vs drift energy E=½k₀² (`hypotheses/qsp_phase5`). Its classical + Lindhard
  overlays are **bulk** references (σ_WP=0.5 ⇔ the bulk σ_q=0.354 `sigma0p35` set);
  slab-WP-vs-bulk is a labelled **geometry estimate**, not a matched pair (ADR 0010).
  WP points carry a convergence flag — high-v are true values, slow/54 eV are
  **upper bounds**. σ_WP=0.5 also imposes a velocity floor v>σ_p=1.0 (k₀>σ_p).
- **WP self-interaction error (SIE)** — the WP is one electron in the KS determinant,
  so it feels its own Hartree (E_H=1/(√(2π)·s)≈21.7 eV for density-std 0.354), only
  partly cancelled by LDA exchange → **residual ≈7 eV**, with no classical
  counterpart. This contaminates S_WP−S_classical (it is artifact, not quantum
  physics) and drives extra self-spreading. Must be bounded by a **vacuum-WP control**
  (same WP, no bath) before any "quantum component" is reported.
- **matched-pair quantum-vs-classical experiment** — the salvageable design: classical
  (regenerated at σ_pot=σ_WP/√2) and WP run at IDENTICAL charge clouds, in the
  non-relativistic non-spreading window (density std s₀≈1.5–2, σ_label≈2.1–2.8,
  E≈0.8–2.5 keV, spread<10%, v<0.1c). ΔS=S_WP−S_classical there is the quantum
  component at *finite* width; a σ-scan {1.5,2.0,2.5} extrapolated to s₀→0 reaches
  toward the point/Lindhard limit (which is itself unmeasurable as a non-spreading WP).
- **spreading limit (why σ<0.5 fails)** — free Gaussian spreads on τ=2σ², so
  "no appreciable spread over the transit" needs 2σ² ≫ t⋆ ⇒ σ ≳ 2.3 Bohr at 100 eV;
  σ<0.5 gives τ=0.5 a.u. ≪ t⋆ (≈21× spread). σ and no-spread are uncertainty-
  conjugate (small σ is *worse*); and σ<0.5 is under-resolved at dx=0.40
  (≲1.25 pts/σ; floor ≈1.5–2 dx). "Bind σ to E to suppress spreading" is therefore
  impossible at the target — abandoned in favour of σ≈3 (rigid, resolved) for the
  full-transit measurement + the t→0 method for the small-σ quantum component.

## Jellium electron-gas analytics (2026-06-17)

Vocabulary for the analytical electron-gas reference notebook
(`ResearchProject/systems/jellium/hypotheses/00_jellium_reference/`).

- **the density knob (`RS`)** — the single notebook control parameter: the
  Wigner-Seitz radius `r_s` of the homogeneous electron gas. The electron count
  `N=162` is held fixed (always a closed shell); the box side `L` is *derived*
  (`L=(N·4π r_s³/3)^{1/3}`, giving `L=50` Bohr at `r_s=5.69`). Every analytical
  quantity recomputes from this one knob. NOT to be confused with a literal
  density `n` input — `n`, `L`, `kF`, `ω_p` are all outputs of `RS`.
- **jellium shell / magic number / closed shell** — free electrons in a cubic
  periodic box are plane waves indexed by integer triples `(nₓ,nᵧ,n_z)`; states
  with equal `|G|²=nₓ²+nᵧ²+n_z²` are exactly degenerate and form a *shell*.
  Cumulative fills at shell closures are the jellium *magic numbers*
  `2,14,38,54,66,114,162,…`; `N=162` closes the `|G|²=6` shell (a *closed
  shell*). Canonical enumeration: `inqkit::jellium::shells` (`shells.hpp`).
- **loss function L(q,ω)** — the *analytical* RPA energy-loss function
  `Im[−1/ε_RPA(q,ω)]`, the quantitative spectral weight feeding the f-sum rule
  and Lindhard stopping. Canonical implementation: `inqview.analysis.lindhard_elf`
  (full complex Lindhard argument; f-sum verified). **Distinct from** the
  **plasmon (peak-)locator** `|n_q(ω)|²/q²` of `inqview.analysis.plasmon_spectrum`,
  which only locates plasmon *peak positions* from rt-TDDFT and is NOT the
  quantitative ELF (the naming was deliberately separated — see that module's
  docstring and `docs/validation/loss-function-formula-validation.md`).
- **discrete box mode `q_m`** — a wavevector the finite `L`-box actually
  supports, `q_m = 2π m / L` (`m=1,2,3,…`); the axial set `F[0,0,m]` that
  `plasmon_spectrum.extract_axial_nq` measures. The analytical continuous
  `L(q,ω)` / Bohm-Gross dispersion is overlaid with these `q_m` markers to link
  bulk RPA to what the simulated box can resolve.
- **Rayleigh resolution time `T_min`** — the minimum rt-TDDFT propagation
  *duration* for an FFT to separate the loss-function frequencies:
  `T_min = 1/min(Δf)` (Rayleigh frequency-resolution criterion; FFT bin width is
  `1/T`). The frequency set is the per-mode spectral features — plasmon lines
  `ω_pl(q_m)` and e-h band edges `ω_±(q_m)` (NOT the continuum interior, nor every
  discrete e-h microline). Energies in Ha are angular `ω`; ordinary `f=ω/2π`, so
  `T_min = 2π/min(Δω)`. Sets *resolution* (total time) — distinct from the timestep
  `dt`, which sets the *Nyquist ceiling* (max resolvable `ω`).

## Graphene CAP (2026-06-18)

Scattering of an electron projectile off graphene with a complex absorbing
potential, replicating Yao & Schleife. Glossary only; locked decisions live in
`docs/plans/graphene-cap.md`.

- **Feasibility replica** — a deliberately reduced version of a published setup
  (here: 4×4 / 32-atom graphene, 60 Bohr z-cell, 3-trajectory ensemble) chosen to
  fit a one-GPU-overnight budget. Its numbers are **methodology-faithful but NOT
  the paper's converged values**; every deviation is tabulated and labelled as
  such in the deliverables.
- **Channeling (A) vs centroid (O) trajectory** — the two perpendicular (+z)
  impact geometries of Yao & Schleife: **channeling** aims the projectile through
  a graphene **hexagon hollow** site, **centroid** through a carbon **atom**.
  Distinct from a *grazing* trajectory (projectile skimming parallel to the
  sheet), which is incompatible with periodic graphene and is not used.
- **Planar-integrated charge-density difference Δn(z,t)** — the projectile-induced
  density change integrated over the two transverse directions, `∫∫ Δn dx dy`,
  plotted vs z and time (the paper's Fig. 1 trace); the most directly
  literature-comparable observable for this system.
- **Whole-system field (KS)** — the gauge-invariant single fields representing all
  electrons: total density `n(r,t)=Σᵢfᵢ|ψᵢ|²` (real) and current density
  `j(r,t)` (vector). There is **no meaningful "sum of orbitals" wavefunction** —
  the many-electron state is a Slater determinant and `Σᵢψᵢ` is gauge-dependent.
  The only single complex wavefunction is the **WP (projectile) orbital**.

## Density decomposition — WP run analysis (2026-06-22)

Canonical vocabulary for decomposing a wavepacket run's density. Applies to every
WP run, not just localised jellium. Glossary only.

- **Total** — all electron density, bath + projectile. Channel `density_total`,
  `n_total(t)`.
- **Wavepacket (WP)** — the injected projectile orbital density `|ψ_wp|²`. Channel
  `density_wp`.
- **Bath** — the responding electron gas *without* the projectile, defined
  run-independently as `n_bath = n_total − n_wp`. The on-disk `density_system`
  label is **not** trusted (WP-included in some runs, bath-only in others).
- **Jellium / background** `n₊` — the *positive*, static neutralising charge. **Not**
  an electron density and never moves; do not conflate with the bath. ("Jellium-only
  system" in loose speech means the **bath** electrons.)
- **Density views (per system):** *absolute* `n(t)`; *Δ-vs-first* (induced)
  `n(t) − n(0)` using the first RT frame per channel — plus, for the bath only, a
  true-induced `n_bath(t) − n_bath^{GS}` panel referenced to the ground state
  (`density_gs_system`); *Δ-vs-previous* (flux, `∝ −∇·j`) `n(t) − n(t−Δt_frame)`.
- **Coordinate order** — inqkit VTIs are in **physical order** (`Origin=−L/2`);
  load via `inqview.io.load_vti` and **never** `np.fft.fftshift` them. Only LEED
  screen `.dat` files are FFT-natural and need a shift.

## Localised jellium (2026-06-21)

A jellium target confined to a finite region of the INQ cell (sphere/slab/box),
into which a projectile is fired from surrounding vacuum. Glossary only; decisions
live in `docs/plans/` + `docs/adr/`; theory in
`docs/notes/localised-jellium-theory.md`.

- **Delocalised (whole-cell) jellium** — the *existing* jellium in the repo. It is
  not an explicit background at all: `extra_electrons(N)` with no ions, and the
  periodic Poisson solver **drops the G=0 component**, which is identical to a
  uniform positive background spread over the *entire* cell. The status quo, and
  the reason today's jellium cannot be localised by parameter tweaks.
- **Localised background** `n₊(r)` — an explicit confined positive charge density,
  e.g. `n₀·Θ(R_cl − |r − r₀|)` (sphere) or a slab/box region, with
  `n₀ = 3/(4π r_s³)` and `∫n₊ = N` (charge-neutral). Replaces the whole-cell
  uniform background with a finite one.
- **Background perturbation** — the chosen *mechanism* (Option A): an `inqkit`
  class implementing INQ's perturbation duck-type whose `.potential(t, v)` adds a
  **static** `v_bg(r) = −poisson(n₊)` to the KS potential. Because both
  `ground_state::calculate` and `real_time::propagate` take the perturbation, the
  background confines the electrons during the SCF *and* persists through the
  projectile flight, with **no edit to `inq/`**. (Rejected alternative: smeared
  positive pseudo-ions, which cannot make a flat-top interior.)
- **Confinement radius `R_cl`** — the sphere radius / slab half-width of the
  background. For a sphere, neutrality fixes `R_cl = r_s · N^{1/3}`.
- **`E_self`** — the classical electrostatic self-energy of the background
  (`(3/5)N²/R_cl` for a uniform sphere); must be tracked explicitly for the
  cluster-energy → HEG-limit benchmark to be meaningful.
- **Spherical-cluster magic numbers (2, 8, 18, 20, 34, 40, 58, …)** — the
  closed-shell electron counts of a *localised spherical* jellium (angular-momentum
  shells 1s,1p,1d,2s,…). **Distinct from** the *periodic-box* shell table in
  `inqkit/jellium/shells.hpp` (2, 14, 38, …, 162), which belongs to delocalised
  jellium. Validating the localised cluster against the box table would be wrong.
- **ΔE_total anomaly (energy-oscillation phenomenon)** — in many localised-jellium
  RT runs, `ΔE_total(t) = E_total(t) − E_ref` does **not** decay monotonically once
  the CAP begins absorbing; it *oscillates*, and in several runs rises **above 0**,
  which is unphysical: the closed system has no energy influx and a CAP can only
  *remove* energy. Seen across WP, effective-mass, heavier-electron, and some
  truncated-classical runs; contrast the `p3_wp` run, where ΔE_total decays to a
  stable plateau as expected. Cause **under investigation** — candidate mechanisms:
  the CAP acting as a non-Hermitian energy *source* in the reported ledger; the
  static background `v_bg` contribution being absent from the reported energy
  functional; a wrong subtracted `E_ref`; propagator/grid numerics; or a
  density-dependent (time-dependent) KS Hamiltonian. Note:
  `docs/notes/localised-jellium-energy-oscillation-investigation.md`; diagnosis
  campaign: `docs/campaigns/localised_jellium/energy-oscillation-diagnosis.md`.
  UPDATE 2026-07-13: diagnosis campaign CONFIRMED (conf 0.90) the CAP-gated
  bookkeeping mechanism (dominant channel: norm-divided kinetic filtering);
  the cap_fix campaign then found setup configs with no observed turn-up —
  but see *period-lengthening reading* below.
- **Period-lengthening reading (of the cap_fix results)** — the user's 2026-07-14
  re-interpretation: the "improvement ladder" of cap_fix may not remove the
  oscillation but *lengthen its period* — the time of the E_total minimum drifts
  monotonically later along the ladder (t_min 21.6 → 27.8 → ~33 → 36.4 → >48),
  every sufficiently long two-sided run turns up, and configs with "no rise" have
  **less than one period of data**. "Monotone to t=T" therefore only means "no turn
  observed yet", never "no oscillation".
- **PBC-vs-open-z channels (Arm A / Arm B)** — "periodic boundary conditions cause
  the oscillation" splits into two distinct mechanisms. **Arm A, density
  recirculation:** the FFT-grid wavefunction always wraps in z (regardless of the
  cell's periodicity setting), so unabsorbed density re-enters and re-interacts;
  only absorbers stop it; testable only indirectly (period vs L_z scaling).
  **Arm B, electrostatic periodicity:** the Hartree/Poisson interaction with
  periodic images along z plus the charged-cell G=0 convention; directly
  switchable in INQ via cell `periodicity(2)` (slab-truncated Poisson) vs
  `periodicity(3)`. A p2 RT run must load a GS converged at p2 (a p3 GS is not an
  eigenstate of the p2 Hamiltonian → spurious t=0 kick). Absolute energies are
  convention-dependent across p2/p3 — compare only ΔE_total(t) shapes. Campaign:
  `docs/campaigns/localised_jellium/pbc-open-z-oscillation.md`.

## Campaigns (2026-06-22)

> Vocabulary for the `campaigns` skill (`.claude/skills/campaigns/`) and the
> `docs/campaigns/` tree (renamed from `docs/prompts/`).

- **campaign** — one prompt `.md` file under `docs/campaigns/<area>/`. It defines
  a single hypothesis-testing run-set that a **fresh agent executes autonomously,
  end-to-end, with no user in the loop**. The *file* is the unit (not the folder);
  "each prompt is a campaign". Authored via the five gated stages of the
  `campaigns` skill and shaped by `docs/campaigns/template.md`.
- **`campaign_autorun` is NOT a "campaign" (naming clash, 2026-07-06).** The
  localised-jellium GS ladder `scripts/campaign_autorun/` is a **run-set / sweep**
  (hypotheses H0–H5, ~90 GS + frozen single-point runs), analysed under
  `hypotheses/campaign_autorun_study/`. It predates the `campaign` glossary term
  above and keeps its folder name; when someone says "the campaign_autorun
  campaign" they mean this run-set, not a `docs/campaigns/` prompt file.
- **area** — the folder grouping campaigns by theme (e.g. `cap_in_jellium`,
  `ml-patterns`, `td-hf`). Purely organisational; carries no status of its own.
- **campaign status** — the `status` frontmatter field, agent-set end-to-end:
  `draft → ready → running → blocked → paused → done` (`blocked` = waiting on a
  dependency; `paused` = deliberately stopped). Together with the `tasks:` list
  (each `{name, done}`) it yields the **`x/N`** progress (`x` = tasks done).
  Frontmatter is the single source of truth; the executing agent flips `done`
  flags and bumps `status` as it runs.
- **campaign INDEX** — `docs/campaigns/INDEX.md`, a status-grouped table
  regenerated (never hand-edited) by `.claude/skills/campaigns/build_index.py`
  from every campaign's frontmatter. The portfolio view of all campaigns at once.
- **autonomy-readiness** — the Stage-5 gate: a checklist whose every box must be
  answerable from the prompt text alone before `status` may become `ready`. A
  compact echo ships in each prompt as the `<preflight>` block for the executing
  agent to re-verify before burning GPU.

## Run inventory: catalogue vs database (2026-06-30)

> Vocabulary resolved while grilling the `ml-patterns/pattern-finding-in-wp-classical-runs`
> campaign. Two distinct artefacts with similar names — do NOT conflate.

- **run catalogue** — the **thin, flag-based** inventory `docs/runs_catalogue.csv`,
  maintained by the `tddft-run-catalogue` skill (`scan_runs.py`). One row per run,
  observables recorded as 0/1 **presence flags**, a small metadata subset, jellium
  + coronene only. Purpose: answer "which runs have observable X?". Deliberately
  coarse; complements (does not feed) campaigns.
- **run database** — the **rich, reproducibility-grade** inventory built for the
  ml-patterns campaign: a **new, from-scratch** artefact (wide canonical CSV +
  nested JSON mirror), one row per run across **all six systems**, holding *every*
  parameter needed to reproduce a run (parsed from `run_summary.txt`, with
  config-header / `run.cpp` fallback) **plus file PATHS** (not flags) to each
  observable, **plus** derived physics features and ML-shape columns. Missing
  values are the literal token `NULL`. Distinct from — and not coupled to — the
  thin **run catalogue**.
- **twin_run_id / matched pair** — the run-database column linking a classical
  projectile run to its wavepacket twin at **matched (σ_WP, v)**. The
  classical-vs-quantum induced-density question is only defined *between* matched
  twins; `pair_width_matched` (bool) guards the documented √2 mismatches (B2/B3,
  legacy `sigma0p35`) so non-matched "pairs" are not silently compared.
- **relevant_to_induced_density** — the run-database boolean that flags whether a
  run carries an electron-gas + induced-density signal (true for jellium/localised/
  cylindrical/coronene projectile runs; false for vacuum free-WP absorber runs and
  WP-only graphene), so the off-topic vacuum sweep can be filtered out downstream.

## Bulk-jellium PDE-discovery redo (2026-07-03)

> Vocabulary resolved while grilling the bulk-only PDE-discovery extension of the
> `ml-patterns/pattern-finding-in-wp-classical-runs` campaign. This redo narrows
> scope to **pure bulk jellium only** (no slab, no other systems) and re-centres
> on discovering a governing differential equation for the induced density.

- **two-track design** — the redo runs two parallel tracks judged together at the
  end. **Track A** = the falsifiable gates (form-factor + wake) re-run clean on
  bulk-only. **Track B** (headline) = open governing-equation discovery. A
  synthesis judge reports where the two tracks agree/diverge. Neither track's
  outcome is retried into the other's.
- **induced-density PDE / field PDE** — the Track-B headline object: a governing
  partial differential equation for the bath response `n_bath(r,t) − n_bath^GS`,
  discovered by **weak-form SINDy / PDE-FIND** from a **broad, agnostic** operator
  library (minimal physics priors). Physical meaning of surviving terms
  (restoring `−ω_p²n`, dispersion `∇²n`, advection `v·∇n`, projectile source) is
  assigned **post-hoc** by interpretation, not baked into the library. A
  **latent ODE** (SINDy on POD modes) is the reduced-order support, not the
  headline.
- **separate-then-compare** — `PDE_classical` and `PDE_WP` are each discovered on
  their **own** runs (classical = coulombic point sweep; WP = σ-fixed velocity
  sweep), *not* on the difference field. Comparison happens only afterwards
  (shared vs unique active terms, coefficient ratios). Contrast with the prior
  campaign's difference-first `Δn = n_WP − n_classical` framing.
- **three validation walls** — a discovered PDE term counts as "physics" only if
  it survives all three: (1) a **pinned** calibration/held-out cell split (library
  + sparsity tuned on calibration, coefficients reported from held-out; extends
  ADR 0011); (2) **temporal forward-prediction** (fit on a held-out trajectory's
  early window, integrate forward, score the later window); (3) **bootstrap
  coefficient stability** (term persists across resampled calibration subsets).
- **aim-reached gate** — the autonomous orchestrator's stop condition: both
  deliverables exist and validate (Track-A held-out verdicts rendered; Track-B
  PDEs pass the three walls for classical AND WP + comparison done), OR the
  validation metric plateaus (K rounds no gain), OR the **12 h wall-clock cap** is
  hit — then it stops and reports best-so-far honestly. A refute/inconclusive/
  partial is a valid reported outcome, never retried into a confirm.

## Nazarov–Gross mass sweep (2026-07-11)

- **NG theorem (working name)** — the Nazarov & Gross 2025 (arXiv:2510.26222)
  claim under test: projectiles of the same charge but different mass, moving at
  the same velocity, experience *different* friction in an electron liquid — a
  purely quantum effect tied to the projectile's spatial extent. Derived in the
  **slow-projectile (friction) limit**; in the fast limit same-charge projectiles
  converge to mass-independent (Lindhard/classical) stopping.
- **null branch** — a fixed-velocity mass sweep run in the *fast* regime
  (v ≫ v_F), where NG predicts **no** appreciable S difference across masses.
  Serves as the control/noise-floor half of the theorem test: an observed
  splitting there is an artefact channel (packet spreading, SIE, grid), not NG
  friction.
- **slow branch / discriminating regime** — the sub-Fermi-velocity companion
  sweep (v ≲ 0.5·v_F) where NG's mass-dependent friction S = Q(m)·v is the
  predicted, measurable signal. Only this branch can positively validate the
  theorem; the null branch alone can merely fail to falsify it.
- **fixed-velocity mass sweep** — the sweep design in which the projectile
  VELOCITY is the invariant across rungs (k₀ = m·v and E = m·v²/2 both scale
  with m). Chosen because stopping theory parameterises S by velocity; contrast
  fixed-k₀ or fixed-energy sweeps, where v differs per rung and S(v) curvature
  contaminates the mass comparison.
- **aliasing mass ceiling** — in a fixed-velocity sweep the packet's momentum
  envelope (centre k₀ = m·v, fixed width σ_p = 1/(√2σ_WP)) slides toward the
  grid Nyquist edge as m grows; the ceiling is the largest m whose aliased tail
  passes the cutoff guard. It is a *fast-branch* constraint — at slow v the
  ceiling sits far above any mass of interest.
- **spreading systematic** — the mass-dependent free-packet dispersal rate
  (∝ 1/m at fixed σ_WP) that can imprint an apparent mass dependence on S even
  where NG predicts none; any flatness/violation verdict must be checked against
  the measured packet width channel before being attributed to NG physics.
