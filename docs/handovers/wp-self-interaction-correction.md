# Handover: WP self-interaction correction (SIC) — review, implement, validate, run

**Rolling file. Latest milestone at top.**
**Repo:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research` (branch `quantum-stopping-power`)
**Machine:** CSD3, login node `login-q-2`, GPU partition `ampere` (A100, sm_80)
**Plan:** `docs/plans/wp-self-interaction-correction.md` (REVIEWED + AMENDED 2026-08-02, §0)
**Parent study:** `docs/handovers/cylindrical-channeling-ks-stopping.md` (2026-08-02 5th:
the ~20 % WP stopping deficit is dominated by LDA self-interaction error)

---

## 2026-08-05 — REPORT-2 SUBSECTION WRITTEN: "Self-Hartree quantification"

The self-interaction work of 2026-08-02 is now a written subsection of report 2
(`sec:results-self-hartree`), sitting immediately after the bulk-jellium
chapter. Six figures + three paragraphs + a generated numbers file. **Complete
and verified; nothing pending except three bibliography keys the user must add.**

### Where everything lives (absolute paths)

| What | Path |
|---|---|
| section folder | `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/reports/report2/drafts/draft1/figures/self_hartree/` |
| shared loader | `.../figures/self_hartree/_selfhartree_data.py` |
| 6 generators | `.../figures/self_hartree/make_sh_{epp_scaling,xc_cancellation,residual_sie,excess_spreading,excess_vs_sigma,epp_in_medium}.py` |
| panel copies | `.../figures/self_hartree/panel/` (`PANEL=1 python make_sh_*.py`), all 2.948 × 2.53 in |
| every quoted number | `.../figures/self_hartree/self_hartree_numbers.txt` ← `make_sh_numbers.py` |
| tests (19, passing) | `.../figures/self_hartree/tests/test_selfhartree_data.py` |
| prose + 2 figure envs | `.../draft1/template/report2.tex`, subsection `sec:results-self-hartree` |
| compilable panel test doc | `.../draft1/template/self_hartree_panel_final.tex` |
| production log | `.../draft1/plots_draft1_log.md` (top entry) |

### The argument the subsection makes, in three paragraphs

1. **Magnitude.** E_PP = A/σ_WP with A = 1/√(2π) Ha·Bohr = 10.86 eV·Bohr;
   measured 10.84 (−0.18 %). At σ_WP = 2 that is 4.82 eV, ~5 % of the 100 eV
   projectile, and it has no classical counterpart. E_PP(0) is identical across
   twin halves to <1 meV, which validates σ_pot = σ_WP/√2 from the data.
2. **Cancellation and persistence.** LDA xc cancels 79.7 → 97.7 % over
   σ_WP = 1 → 8; the residual is 0.2–0.8 eV and falls as ~σ^−2 to σ^−3. It is
   not static: in vacuum the residual self-field makes the packet spread 5.3 %
   (LDA) / 22.9 % (Hartree-only) faster than free at σ_WP = 2, τ = 1.875, and
   SIC-PZ returns it to free to 1e-8.
3. **The bound.** Three independent measurements: (i) in the bulk run the
   residual drifts −0.45 eV over the fit window against a 30.0 eV
   classical-minus-WP drift-energy gap = **1.5 %** (σ_WP = 3: −0.32 vs 25.8 eV,
   1.2 %); (ii) the direct SIC run of 2026-08-02 recovers ~9 % of the ~20-point
   cylindrical deficit; (iii) across the 2.92× density lever the stopping ratio
   moves 13 % while ΔE_PP moves 0.2 %.

### THE NEW RESULT (structural, not fitted) — carry this forward

**For a Gaussian one-electron density, U and the Dirac exchange energy both
scale as 1/a, so |E_x|/U is a PURE NUMBER: 0.67840955 at every width.** LDA
exchange therefore removes a fixed 67.84 % of the self-Hartree and *all* width
dependence of the cancellation comes from correlation. This is why the residual
dies far faster (~σ^−2 to σ^−3) than E_PP itself (1/σ), and it is what makes the
bound quotable rather than hand-waved. Confirmed two ways (closed form and
radial quadrature of ∫n^{4/3}, agreeing to 4e-16) and asserted in the suite.

**This also refines, rather than contradicts, the intuition "use a big enough
packet".** E_PP does fall as 1/σ_WP, but at fixed DIMENSIONLESS time the excess
spreading GROWS with σ (the coupling ~ σ, since KE ~ σ^−2 while U ~ σ^−1). It is
only at fixed PHYSICAL duration that widening helps — 7.9 % → 2.3 % from
σ_WP = 4 to 8 — because a wider packet evolves more slowly. Both readings are in
figure `fig:self-hartree-bound`(a); quoting only one of them is how this gets
stated wrongly.

### Validation actually performed (not assumed)

- **The analytic kernel is not a model of the functional — it IS the
  functional.** It reproduces INQ's own `exc_self_ha` from all six vacuum
  `sweep_s*_sic_pzrun/raw/observables/sic.csv` to **<1e-9 relative**, asserted at
  load in `vacuum_exc_measured()`, and independently reproduces the cylindrical
  study's measured σ_WP = 4 terms (2.7139 / −1.8412 / −0.6320 eV).
- Vacuum protocol gates re-asserted at load: free reference vs closed-form
  spreading 4.8e-12, var(p) drift 8.0e-13, E_PP(0)·σ constant to 1.2e-12 Ha.
- Twin-pair gates re-asserted at load: E_PP(0) equal across halves, a(0) = σ/√2,
  classical E_PP rigid.
- **19 tests, all passing**; catalogue row added to
  `docs/validation/test-catalogue.md`. One test caught a mistyped literal in its
  own expected value (0.678416 vs the true 0.67840955) — the derived-value rule
  working as intended.
- Margin check: all 12 PNGs clear the canvas edge by ≥15 px.
- `template/self_hartree_panel_final.tex` compiles; both figures render with
  complete captions on one page each (pages checked as images, not just exit
  codes). `report2.tex` compiles once `placeins.sty` is available — it is NOT
  installed in this environment, which is pre-existing and unrelated.

### Traps hit, all now documented in the plots log

1. **Six half-slot panels do not fit one figure** — 3 rows leave ~1.3 in for the
   caption and it ran off the page mid-sentence. Split into two figures by
   argument (`fig:self-hartree` a–d, `fig:self-hartree-bound` a–b).
2. **`style.figure_two_col()` clips its x-label at 3.0 in** — it uses
   `add_subplot` default margins, not the theme's fractional rect. Caught by the
   ink-bbox margin check, invisible by eye. Use
   `fig.add_axes(style._ONE_COL_AXES_RECT)`.
3. **|x| on a log axis hides a sign change.** The residual changes sign between
   the isolated and periodic bases; plotted as |·| each crossing became a spike
   that reads as a numerical glitch. Signed + linear.
4. **Two box conventions give two DIFFERENT laws.** Scaled sweep: E_PP(0)·σ
   constant. Fixed production box: A/σ − C, C = 0.6016 eV. Merging them is
   wrong by the Madelung offset; SH(a) draws the fixed-box line only over the
   σ = 1.6–3.8 span its two anchor points support.

### Outstanding — needs the user

1. **Three bibliography keys** are cited and do not yet exist in
   `references.bib`: `perdew1981sic` (Perdew & Zunger, PRB 23, 5048 (1981)),
   `leslie1985` (J. Phys. C 18, 973), `makov1995` (PRB 51, 4014). They are listed
   in a comment block at the head of the subsection in `report2.tex`.
2. The subsection forward-references `sec:results-cylindrical-jellium` for the
   direct SIC test; that section is still `% TODO.` in the skeleton.
3. Nothing committed to git (not asked).

---

## 2026-08-02 (5th) — Report-ready figure set produced (user request)

`docs/reports/report_figures_self_hartree_removal/` — two subfolders, all PNGs
at 600 dpi, canonical theme (ADR 0004), time axes in fs, individual one-column
figures (no titles; captions belong to LaTeX), scripts alongside the figures:

- `xc_correction_analysis/` (vacuum five-run study; script
  `make_xc_correction_figures.py`, reuses the hypotheses data layer):
  `00_setup` (initial-density x–z slice, linear|log, σ_WP/cell annotated),
  `01a` width vs analytic, `01b` excess spreading, `01c` var(p) conservation,
  `02a` E_PP, `02b` var(p)/2m, `02c` released-vs-absorbed, `03a` log-distance
  from the closed form, `03b` subtracted self-terms, `03c` E_corr vs E_tot drift.
- `cylindrical_jellium/` (uncorrected twin WP vs SIC replica; script
  `make_cylindrical_comparison_figures.py`): `01a/01b/01c` E_PP / E_PS / ΔE_SS,
  `02a–02d` ΔE components (total/kinetic/hartree/xc), `03a` σ_perp with free
  dispersion + bore radius 10 Bohr, `03b` σ_z.
- Colour contract across both folders: LDA-uncorrected red, full SIC green.
- Two render defects found by visual inspection and fixed: bore-radius label /
  legend collision in `03a_sigma_perp`; 3-decimal ticks clipping the y-label in
  the ΔE component plots (MaxNLocator(5)).

---

## 2026-08-02 (4th) — Vacuum (Tier V) notebook extended to all FIVE runs, on user request

The vacuum study notebook
`ResearchProject/systems/vacuum/hypotheses/wp_selfinteraction/selfinteraction.ipynb`
previously covered only the three-theory difference measurement
(noninteracting/hartree/lda); the two SIC runs existed on disk but were verified
only through in-run gates. Extended:

- `selfinteraction.py`: `THEORIES` → 5 (+`INTERACTING` tuple, labels/colors);
  `load()` reads `sic.csv` when populated (header-only ⇒ `u_self is None`);
  `load_all` skips absent theories but requires the reference (`_smoke` suffix
  still works); `summary_table` rows for all interacting theories, xc-difference
  row preserved.
- `build_selfinteraction_notebook.py`: five-curve spreading/var(p)/energy plots;
  NEW section 3 "The correction, verified against the closed form (Tier V)" —
  `numerics_gate` applied to `sic_pzrun` (the reference's own gate reused as the
  Tier V criterion), fig `si_figs/03_sic_verification.png`: (a) log-distance
  from the closed form, (b) subtracted self-terms u_self/exc_self, (c)
  E_corrected vs E_total drift. Sections renumbered (energy→4, answer→5, GIF→6,
  now 5 GIF batteries).
- Tests: +3 in `tests/test_selfinteraction.py` (load_all tolerance, header-only
  sic.csv, planted-perfect-correction passes the gate + ~0 % excess row) —
  **13/13 PASS**. Catalogue row appended (`docs/validation/test-catalogue.md`,
  SIC section).
- Also FIXED the (pre-existing) density-GIF cell: it called
  `make_density_gif_battery(total_subpath=...)`, a keyword the function does not
  have, so every earlier build silently fell back to "battery unavailable" and
  the notebook shipped without its mandatory animations. Now calls the real
  signature (`run_label`/`slab_face`/`cap_inner`, `cap_lines=()` for the vacuum
  box — no slab/CAP guide lines) and unpacks the `(gifs, vmax)` return.
- FINAL notebook: 17 cells, **0 errors, 15 embedded GIFs** (5 runs × 3 kinds,
  total-density only — vacuum has no bath), 67 MB, figures in `si_figs/`
  (01_spreading, 02_energy, 03_sic_verification + 5 battery dirs). Tier V gate
  reprinted in-notebook: sic_pzrun |σ/analytic−1| = 1.6e-8, var(p) drift
  8.3e-5, PASS with the reference's own gate; sic_h −44.5 % width (self-bound),
  closure 0.0 Ha exact for all four interacting runs.

---

## 2026-08-02 (3rd) — PRODUCTION COMPLETE. RESULT: the deficit is ~91 % GENUINE, only ~9 % SIE — the parent study's "substantially an artefact" reading is OVERTURNED

Chain finished green end-to-end: 32615190 (tierb) 7 m 24 s, 32615191 (prod)
52 m 36 s, exit 0. Production SIC integrity was essentially perfect:
`max_overlap_pre` ≤ **4.1e-9** all run (the projection never had real work to
do), `cum_norm_removed` ~1e-14, `E_corrected` drift **−1.1e-5 eV over 1500
steps** (the D2 drift channel exists but is numerically negligible here).
`E_total` fell 0.178 eV — the correction doing its intended work, not an error.
u_self(t=0) = 1.9361 eV matches the independent `interactions.csv` E_PP exactly
(two code paths, one number).

### The three-way table (fits via the SAME refined.fit_in_window code path,
### classical always over the same window; S in eV/Bohr)

| window | est. | S_wp_SIC | S_cl | ratio SIC | ratio UNCORRECTED |
|---|---|---|---|---|---|
| **9–25** | T1 | **0.08962** | 0.10956 | **0.818** (r² 0.9995) | 0.801 |
| 21–30 | T2 | 0.07767 | 0.13217 | 0.588 | 0.634 |
| 5–20 | T2 | 0.02645 | 0.09791 | 0.270 | 0.132 |

Whole-run impulse ratio 0.784 (was 0.764). ⟨r_perp⟩/free at t=30: **1.399**
(was 1.467). var(p)/2m growth **+1.988 eV** (was +2.139). var(p_z) growth
+29.0 % (was +44.5 %). f_bore end 0.496 (was 0.457). E_PP 1.936 → 0.339 eV
(was → 0.292) — still decays, now purely as a SIZE diagnostic since the term
is absent from the WP's Hamiltonian.

### Reading (the §5 conditional, answered)

- **The headline T1 ratio moved 0.801 → 0.818: removing the self-interaction
  recovers only ~1.7 points of the ~19.9-point deficit (~9 %).** The remaining
  ~18 points are genuine quantum kinematics — a real WP-vs-point-charge
  difference, qualitatively consistent with Nazarov & Gross (2025).
- **The excess transverse expansion is NOT SIE-driven.** With the self-field
  removed (verified exact to 1e-8 in vacuum), the packet still expands 1.40×
  faster than free and var(p_perp)/2m still gains ~93 % of its former growth —
  while its former "reservoir" (E_PP) no longer even appears in the WP's
  Hamiltonian. The parent study's energy-balance argument (E_PP release 1.64 ≈
  var gain 1.95 eV) was a coincidence of scale, not causation. The driver is
  the wall/bath interaction itself.
- The parent handover's claim "a KS-orbital stopping power measured this way
  carries an SIE contamination at the tens-of-per-cent level" is corrected to:
  **~2 points out of 20 in the T1 window** (uncertainty: the in-medium PZ xc
  subtraction over-removes attraction where bath density dominates — an O(bath
  overlap) ambiguity, inference — but the bore window has f_bore ≥ 0.95).
- T2 remains window-dominated (0.270 vs 0.588 across windows) — unchanged
  conclusion: T2 is not a usable estimator here.

### Deliverables — DONE and verified

- `hypotheses/channeling_sic/wp_wp_sic.ipynb` — 53 cells, **0 errors** (full
  run-notebook battery, same builder + args as the twin's wp notebook).
- `hypotheses/channeling_sic/refined_analysis.ipynb` — 29 cells, 2.58 MB,
  **0 errors** — cell-for-cell identical structure to the twin's (same builder,
  windows baked in); 13 figures in `refined_figs/`;
  `refined_stopping_summary.csv` written (path defs agree to 0.05 % —
  Ehrenfest consistency held).
- Ledger closure on the SIC run: |E_SS+E_PS+E_PP − E_hartree| ≤ 5.0e-11 Ha,
  |E_SB+E_PB − E_external| ≤ 5.0e-11 Ha; WP norm drift 2.7e-9.

### Still pending (deliberately)

- Run catalogue: `scan_runs.py` expects the old `<run>/results/` layout and
  holds NO channeling rows at all (parent study left the same item open) —
  extend the scanner before upserting, don't hack it per-run.
- Journal entry: needs the user's own observation text (rule).
- σ_WP sweep hedge (plan §6): now LOW value — the direct SIC test already
  separated the artefact from the physics.
- Nothing committed to git (user has not asked).

---

## 2026-08-02 (2nd) — ALL GATES GREEN through Tier V; correction verified against the closed form at 1e-8; Tier B running

Chain: 32615188 (chan-tests) COMPLETED 48 s — **all 8 engine tests passed**,
incl. `test_wp_sic_engine` (kick semantics, exact Q re-orthogonalisation, D1
run-consistency). 32615189 (wp-si Tier V) COMPLETED 34 m 47 s — the three
uncorrected runs were already complete (job 32615079 from the earlier session;
skip logic worked) so it ran only sic_h + sic_pzrun.

### Tier V verdict: SIC-PZrun passes its closed-form gates near machine precision

| gate (sic_pzrun) | measured | tolerance |
|---|---|---|
| var(p_z) vs 1/(2σ²) | dev **5.6e-8** rel | 1e-3 |
| σ_dens(t_end) vs √(σ²/2+t²/2σ²) | dev **1.6e-8** rel | 5e-3 |
| ⟨p_z⟩ − k0 (zero-force) | −2.5e-18 | 1e-4 |
| E_corrected drift | 2.9e-5 eV (WARN band; split-operator residual) | PASS<1e-5 / WARN<1e-3 |

SIC-H: zero-force PASS (−3.6e-17), E_corr drift 1.9e-4 eV (WARN band).

### The five-run σ_final(t=30) ordering — every sign as predicted (plan §4)

free reference 6.0104: hartree **8.63** (+44 %, Hartree self-repulsion) >
lda **6.48** (+7.9 %, xc partially cancels Hartree in vacuum) >
**sic_pzrun 6.0104 (free to 1.6e-8)** > sic_h **3.34** (−45 %: only the
attractive xc self-term remains → the packet self-binds). NOTE: in vacuum the
xc self-term largely CANCELS the Hartree one (lda ≪ hartree), so SIC-H alone
would over-correct badly — D1/D3 of the review were load-bearing.

**Decision (plan §4 rule): production variant = SIC-PZrun** (already the
chain's default). 32615190 (tierb) RUNNING; 32615191 (prod) pending afterok.
Background monitor b8s40vqxb watches until prod resolves.

---

## 2026-08-02 — Plan reviewed against literature, TWO real defects found and fixed; SIC implemented wrapper-only; chain being submitted

### User directive (this session)

Review the plan against papers that have done this, list the downsides clearly,
amend, implement, validate thoroughly, produce run notebooks identical to the
channeling_twin refined-analysis + wp-run notebooks in a NEW hypotheses folder,
submit to CSD3 autonomously, analyse when results arrive.

### Literature digested (source notes written, all in docs/sources/)

| Note | What it grounds |
|---|---|
| `messud-2008-tdsic.md` (PRL 101, 096404 — full text read) | variational TDSIC: Lagrange multipliers + symmetry condition ⟨ψβ\|Uβ−Uα\|ψα⟩=0; our Q-projected kick IS the one-sided multiplier scheme → orthonormality EXACT, energy conservation NOT guaranteed |
| `perdew-zunger-1981-sic.md` | SIC functional; one-electron exactness; the polarised-vs-unpolarised xc consistency caveat (D1) |
| `mundt-2007-zero-force-tdkli.md` (PRA 75, 050501(R)) | simplified TD-SIC schemes can violate zero-force/energy conservation SECULARLY → measure, never assume |
| `nazarov-2025-quantum-projectile-stopping.md` (arXiv:2510.26222, Nazarov & Gross) | a genuine quantum-kinematic stopping difference EXISTS (mass-dependent friction) → post-SIC residual deficit is a RESULT, not failure |

### The review verdict (plan §0, D1–D8) — the downsides, explicitly

1. **D1 (would have failed its own Tier V):** plan specified the canonical PZ
   *polarised* xc self-term `v_xc[n_wp,0]`; INQ runs are spin-restricted, so the
   Hamiltonian contains the *unpolarised* `v_xc[n_wp]` — subtracting the
   polarised form leaves a ~26 % exchange remnant and the free-dispersion gates
   fail for a CORRECT implementation. Fixed: subtraction evaluated
   run-consistently (unpolarised `XC_LDA_X`+`XC_LDA_C_PZ` through INQ's own
   `xc_term::evaluate_functional`). Variants renamed **SIC-H** / **SIC-PZrun**.
2. **D2 (wrong primary gate):** `E_corrected` exact conservation is impossible
   for the projected one-sided scheme in the jellium
   (dE/dt = 2 Im Σ_j ⟨wp|h|j⟩⟨j|v|wp⟩; grounded in Messud 2008). Now: hard gate
   in vacuum only; measured + soft-gated (<0.1 eV, non-secular) in jellium.
3. **D3:** SIC-H and SIC-PZrun *must* differ in vacuum (xc self-binding) — the
   old "if indistinguishable choose SIC-H" rule was vacuous; Tier V now
   *separates* Hartree and xc shares of the SIE.
4. **D4 (stated limitation):** bath still feels `v_xc[n_S+n_wp]` vs classical's
   `v_xc[n_S]` — not removed, arguably physics; caveat for the notebook.
5. **D5:** tagged-orbital identity relies on INQ propagating columns
   independently (it does); full double-set TDSIC / ADSIC / TD-KLI all rejected
   because they modify bath dynamics or mix the projectile away.
6. **D6:** Tier V never exercises Q (no bath) → engine test + Tier B carry that.
7. **D7:** the disk blocker is stale — 127 GB free measured; nothing deleted.
8. **D8:** expectation reframed as a SPLIT of the deficit (SIE artefact vs
   genuine quantum kinematics), not a collapse to the classical value.

### Implemented (all compiles clean on the login node, CUDA 12.1 toolchain)

| Path | What |
|---|---|
| `inq-stack/include/inqkit/wavepacket/self_interaction_correction.hpp` | NEW. `SelfInteractionCorrection{Mode::none/hartree/pz_run}`: per-step kick `ψ_wp ← N·Q·exp(+i dt v_SIC)ψ_wp`; v_SIC = poisson(n_wp) [+ unpol LDA vxc via INQ's own evaluate_functional]; modified-GS projection over occupied states returns `max_overlap_pre`/`norm_removed`; `measure()`, `corrected_energy()`. Wrapper-only, single-rank/gamma. Note: internal helpers are public — nvcc forbids device lambdas in private methods. |
| `inq-stack/tests/include/inqkit/wavepacket/test_wp_sic_engine.cpp` (+ CMake reg) | NEW engine test, 3 cases: kick semantics (density/⟨p⟩ invariant, var(p) responds); exact Q re-orthogonalisation with leak reported + bath columns untouched; **D1 run-consistency** (`u_self==energy_hartree`, `exc_self==energy_xc` at 1e-9 for a 1-e system, polarised PZ asserted OUT of tolerance). Compiled OK; RUNS in the chan-tests job (login node has no GPU). |
| `ResearchProject/systems/vacuum/scripts/wp_selfinteraction/run.cpp` | EXTENDED (was the 3-theory difference measurement from earlier today, smokes all passed, `noninteracting` prod complete in 854 s). Adds `WP_SIC=none/h/pzrun`, Strang half/full/half kicks in the callback, `sic.csv`, END-OF-RUN CLOSED-FORM GATES (exit 4 on failure). Rebuilt OK. |
| `shared/bin/run-wp-si.slurm` | now runs 5 configs (3 theories + sic_h + sic_pzrun), skips completed outputs, prints gate summaries |
| `ResearchProject/systems/cylindrical_jellium/scripts/channeling_sic/wp/run.cpp` | NEW production run: clone of channeling_twin/wp + SIC (CH_SIC=pzrun default, CH_SIC_TIER=b/prod), sic.csv (segment-suffixed), resume-aware Strang boundary (`sic_boundary=` in rt_state.txt; closed→opening half-kick on load; mixed-mode resume refused), Tier-B hard gates / prod soft gates. E_total explicitly NOT a gate under SIC. **Building in background on login node** (fresh CMake tree; ~74 % through libxc when last checked). |
| `shared/bin/run-chan-sic.slurm` | stages tierb (200 steps, builds, hard gates, ckpt cleanup) / prod (1500 steps, execs ./run) |
| `shared/bin/submit-channeling-sic.sh` | the chain: chan-tests → wp-si(Tier V) → chan-sic tierb → chan-sic prod, each afterok |
| `shared/bin/run-chan-tests.slurm` | + `test_wp_sic_engine` in the gate list |
| `ResearchProject/systems/cylindrical_jellium/hypotheses/channeling_sic/` | NEW analysis folder: `refined.py`/`channeling_stopping.py` SYMLINKED to channeling_twin (identical by construction); `build_refined_notebook.py` wrapper (env-redirects CHAN_WP_RESULTS to the SIC run, CHAN_CL_RESULTS to the twin classical, `--wp wp_sic`, `--out-dir` here); `build_run_notebooks.py` (per-run wp_sic deep dive via the run-notebook skill builder, same args as twin) |
| `.../channeling_twin/build_refined_notebook.py` | + `--out-dir` flag (backwards-compatible; default unchanged) |
| `docs/validation/test-catalogue.md` | + SIC section (engine + run-level gate rows) |

### Validation state

- Compile: header + engine test + vacuum binary VERIFIED on login node.
  channeling_sic binary: background build in progress at last check.
- Engine test execution, Tier V gates, Tier B gates: PENDING (SLURM chain).
- The chain is self-gating: any failure blocks production (afterok).

### What is NOT done

- Chain not yet submitted at the time of this handover write (submission is the
  immediate next step once the channeling_sic build finishes — the build is not
  strictly a prerequisite since tierb rebuilds via inq-run, but a login-node
  compile error is cheaper to fix before queueing).
- Notebooks: built only AFTER prod completes (completeness gate in builders).
- Fit windows for the refined notebook: apply the user's windows (9–25; 21–30;
  5–20) in the section-6 parameter cell after the build, then re-execute.
- σ_WP sweep hedge (plan §6): deferred deliberately.
- Journal entry: needs the user's own observations (rule), not written.
- Run catalogue upsert: after prod completes.

### Resume instructions (any session)

1. `squeue -u skcb2` — look for chan-tests / wp-si / chan-sic jobs; logs are
   `<name>-<jobid>.out` in the repo root.
2. If the chain has not been submitted: check `/tmp/chansic_build.log` or just
   run `./shared/bin/submit-channeling-sic.sh` (tierb rebuilds regardless).
3. If wp-si FAILED: read the `[FAIL]` lines in `wp-si-*.out` — a sic_pzrun gate
   failure means the implementation, not the physics (plan §4 decision rule).
4. After chan-sic prod completes: from `hypotheses/channeling_sic/`,
   `PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3
   build_run_notebooks.py && ... build_refined_notebook.py`, then edit the
   section-6 window cell (`T_WIN_CL`/`T_WIN_WP`) and re-execute.
5. Compare S ratios against the uncorrected study's table (handover
   `cylindrical-channeling-ks-stopping.md`, 2026-08-02 4th: T1 9–25 ratio
   0.801) — the headline is how much of the 20 % deficit the correction
   removes.
