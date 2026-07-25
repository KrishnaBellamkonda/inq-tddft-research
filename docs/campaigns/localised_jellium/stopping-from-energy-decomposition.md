---
id: stopping-from-energy-decomposition
area: localised_jellium
title: Stopping power from the decomposed energy ledger (classical + wavepacket)
status: draft
hypothesis: "The decomposed energy ledger (pairwise Coulomb E_PP/E_PS/E_SS/E_SB/E_PB/E_BB + kinetic split + E_xc) admits a stopping-power definition that (i) for a classical Gaussian-charge projectile reproduces the headline deposit S=E_absorbed/L_slab (Definition 2) to within a stated tolerance, and (ii) generalises to the wavepacket without invoking an ill-defined projectile KE — yielding a consistent, gauge-clean S across the select plateaued+decomposed localised-jellium runs and, aggregated over same-system/different-energy runs, a coherent S(E0)/S(v0) curve."
handover: docs/handovers/stopping-from-energy-decomposition.md
tasks:
  - { name: "Phase 1 — notebook §A/B/C0/C/D: enumerate every recorded energy; verify the WP/classical closure relations numerically; MEASURE E_sp(0) vs E_sp(t_final) on a reference run and decide the E_sp treatment (§C0); propose the ranked S-from-decomposition formulae C1/C2/C3 with caveats (Fable 5, §C); brainstorm energy-structure-simplifying setup modifications (§D, Fable 5). MANUAL GATE = user approves formulae + E_sp treatment + setup ideas", done: true }
  - { name: "Phase 2 — lock the approved formula into a tested extraction kernel (E_sp term exposed switchably); formula-validation agent + code-test known-case (classical decomposition-S == Definition-2 deposit within tol; ledger closure ~1e-9; synthetic fixture). Done: tests green + user+agent agree", done: false }
  - { name: "Phase 3 — USER-ASSISTED run pick; apply C1 (classical + WP-with-caveat), C2 (full plotted chain: ΔE(t), ΔKE_proj(t), fitted region, best-fit line, value), C3 (ZPE/localisation subtracted) per run + C1-vs-C2 classical agreement; per-run sanity (energy conservation, gauge test, plateau). Done: S computed + sane for every selected run", done: false }
  - { name: "Phase 4 — A1 DISCREPANCY INVESTIGATION (Fable 5): why is deposit-based S(v) >> Lindhard bulk? find the plateaued-WP runs; rough E_sp-subtraction recompute of E_absorbed vs Lindhard; test (a) localised!=bulk (compare classical), (b) decomposition-sum vs E_total baseline, (c) CAP total-energy anomaly. Done: each hypothesis tested with an experiment + a stated conclusion", done: false }
  - { name: "Phase 5 — aggregate by system, sweep projectile energy → S(E0) headline + S(v0) companion (classical & WP overlaid), Lindhard eyeball overlay (non-gating), write S_decomposition_summary.csv", done: false }
  - { name: "Phase 6 — synthesis: fold the A1 conclusion + the §D setup ideas into a short next-batch recommendation for clean S extraction. Done: one notebook holds all phases, executed end-to-end; handover + INDEX + frontmatter done flags updated", done: false }
blocked_reason: "Phase-4 critical analysis (§E) done EARLY on existing data (2026-07-22, user-directed, no re-runs): the WP deposit-based S ~8x Lindhard is a WP-method artifact (deposit exceeds drift KE at low v = WP internal energy leaking through an energy-lossy CAP); (a) refuted, (b)+(c) one mechanism; classical S=0.25 eV/Bohr the only trustworthy number. Phases 2-3,5-6 BLOCKED on a rebuilt observable (bath-only decomposed / vacuum-referenced deposit) carried through a dynamic, fully-absorbing, ledger-equipped run (does not exist yet). User deferred the re-run; awaiting go."
---

# Stopping power from the decomposed energy ledger (classical + wavepacket)

<identity>
You are a scientific computing researcher running first-principles rt-TDDFT in INQ.
You adhere to this repository's rules, skills, and workflows (twin-run-analysis,
stopping-power-extraction, notebook-making, code-test, formula-validation,
simulation-validation, literature-review, handover-update, scientific-grounding).
σ-convention is UNIFIED: σ means the wavepacket σ_WP; the classical Gaussian
charge std = σ_pot = σ_WP/√2 (`.claude/rules/sigma-wp-convention.md`, CONTEXT.md).

**Model routing (user decision, this campaign):** use **Fable 5** for the
conceptual brainstorm / formula-suggestion work (Phase 1 §C, any later definition
refinement); use **Sonnet** for data-fetch and simple per-run analysis
(column enumeration, CSV loading, per-run number-crunching). Reserve the strongest
reasoning for the formula-derivation and formula-validation steps.
</identity>

<description>
This campaign **derives a stopping-power formula from decomposed energy** and
**applies it** to the localised-jellium runs that both (a) carry the full pairwise
energy ledger and (b) reach a clean energy plateau. It is the **analysis dual** of
the `classical-highdensity-sv` campaign, which is the *data-generation* half: that
campaign records the full ledger with clean-exit plateaus and explicitly reserves
**"Definition 1 — energy-component decomposition (formula still being derived by
the user; DATA-COLLECT ONLY)"**. THIS campaign realises Definition 1.

**No GPU runs are launched here.** The data is the plateaued+decomposed runs from
`classical-highdensity-sv` (once its Phase 3 sweep completes) plus any other select
runs the user identifies at the start of Phase 3. The only new code is a **derived
observable** — the Definition-1 extraction kernel — which is `code-test` +
`formula-validation` gated before it produces any headline number.

**Why this exists.** Stopping power S = energy transferred from the projectile to
the target electrons per unit projectile path length. Two dual views: the
projectile's energy LOSS per unit length, and the target's energy ABSORPTION per
unit length. For the **classical** projectile (an external potential) these are
exact duals by construction — `−ΔKE_proj = ΔE_electronic + ΔU_proj_bg`. For the
**wavepacket** the projectile is one more Kohn–Sham electron, so "the projectile's
energy" is NOT a partition-free observable; the project convention forbids
`−dKE_proj/ds` for the WP and reads stopping from the total electronic deposit
(`feedback_quantum_stopping_not_from_projectile_ke`,
twin-run-analysis rule 6). The decomposed ledger opens a middle path: build S from
named, attributable channels (KE_slab, E_SS, E_SB, E_xc, and — where defensible —
the projectile-partition terms) so that classical and WP are computed by the *same*
functional of the *same* columns and are therefore comparable apples-to-apples.

**Decision it informs.** Which decomposition-based S definition to adopt as the
project's Definition 1, and what the classical-vs-WP stopping comparison actually
is on the high-density benchmark slab.

**Success criterion.** (1) The classical decomposition-S reproduces the headline
Definition-2 deposit `E_absorbed/L_slab` to within a stated tolerance on a reference
run (the built-in self-test). (2) The approved WP formula is well-defined and
gauge-clean (its gauge-invariant terms Δ(E_SS,E_SB,E_BB)≈0 across the twin; it does
not rest on an individual gauge-dependent Hartree/external term). (3) Aggregated
over a same-system energy sweep it yields a coherent S(E0)/S(v0) curve.
**Failure.** No admissible channel combination reproduces the classical deposit, OR
the WP formula is irreducibly gauge-dependent (individual pair terms do not cancel
in the twin difference) — reported honestly as a refutation, never retried into a
confirm.
</description>

<observables_set>
**Consumes existing observables only — emits NO new simulation observable.** Per
selected run, already on disk (ADR-0006 set + full decomposition, see
`classical-highdensity-sv` `<observables_set>`):
- **Lumped stores:** `energy_total, energy_kinetic, energy_hartree, energy_xc,
  energy_external` (+ `energy_nonlocal, energy_ion, energy_ion_kinetic` where
  present) — `raw/observables/observables.csv`.
- **Pairwise Coulomb ledger:** `e_ss, e_pp, e_ps, e_sb, e_pb, e_bb`
  (+ `e_hartree_check, e_external_check, norm_*`) — `raw/observables/interactions.csv`
  (`inqkit/jellium/interaction_energies.hpp`).
- **Projectile track (classical):** `proj_z, proj_vz, energy_proj_ke,
  energy_proj_bg_ideal` — `raw/observables/projectile.csv`.
- **WP diagnostics (WP):** centroid/σ_z where recorded (else reconstructed in
  post from density frames per the canonical-bath-density / VTI-coordinate rules).

**NEW derived observable (pre-gated):** the **Definition-1 stopping kernel** — the
approved decomposition-S functional. It is DERIVED (post-hoc from primary columns),
so it is routed through `code-test` + `formula-validation` + a catalogue row in
Phase 2 BEFORE it produces any headline S in Phases 3–4. No new simulation kernel.
</observables_set>

<resolved_decisions>
All locked via grill-with-docs (2026-07-21). Engine/closure facts carry file refs.

**Framing.**
- New standalone campaign, bound to `classical-highdensity-sv` as its primary data
  dependency (that campaign's reserved "Definition 1"). Deliverable = a formula +
  ONE executed Jupyter notebook holding all four phases. No GPU runs here.

**Candidate formula shortlist (seeds Phase 1 §C; user approves the final set).**
Ranked from the brainstorm (Fable 5, full A–E menu kept as a notebook appendix):
- **Headline — matched-estimator total-deposit** (both runs, identical functional):
  `S = ΔE_target / Δs`, with `ΔE_target(t) = E_electronic(t) − E_electronic(0)`
  and the **channel split** `ΔE_target = ΔKE_slab + ΔE_SS + ΔE_SB + ΔE_xc`
  (each Δ relative to t=0). For the classical run this equals the Definition-2
  deposit `E_absorbed`; the notebook proves that equality as the self-test.
- **Classical conservation anchor (B1):** `S(v0) = −d(½·m·proj_vz²)/ds` over the
  early v≥0.85·v0 window (light-projectile rule) — sanity channel, never the WP
  headline.
- **WP projectile-partition (B3′, exploratory, gauge-checked):**
  `E_proj^WP = KE_proj + E_PP + E_PS + E_PB`, corrected by a WP-in-vacuum baseline
  to remove free-dispersion drift of E_PP (∝1/σ) and KE_internal;
  `S = −d[E_proj^WP − E_proj^vac]/ds`. Reported only if it passes the gauge test.
- **Target-absorption (A1)** and **irreversibility qualifier (D1/D2, matched-face
  / hysteresis-loop)** as accompanying diagnostics.
- **Path length ds:** projectile arc-length ∫|v|dt (classical); slab thickness
  L_slab for the aggregate deposit; WP centroid arc-length ONLY while the density
  is unimodal (flag bimodality). Every reported S carries its ds choice.

**Aggregation.** Group by *system* (same σ_WP / r_s / slab geometry); sweep
*projectile energy* → **S(E0)** headline + **S(v0)** companion; classical and WP
overlaid on one axis (same estimator). Lindhard / linear-response overlay is an
eyeball comparison only (NEVER a gate — `feedback_fourier_loss_function_gate`).

**Gauge discipline.** Individual `E_hartree`/`E_external` differences are Poisson
G=0-convention dependent in the net-charged WP cell; only gauge-clean combinations
are physical (`reference_charged_cell_hartree_convention`). The gauge test
(Δ(E_SS,E_SB,E_BB)≈0 between twins) must pass before any pair-term-based WP number
is quoted.

**Closure relations (code-verified ~1e-9/1e-10, asserted in
`inqkit/jellium/interaction_energies.hpp:17-18` and
`reference_twin_pairwise_decomposition`):**
- classical: `energy_hartree = E_SS` ; `energy_external = E_SB + E_PS`
- WP: `energy_hartree = E_SS + E_PS + E_PP` ; `energy_external = E_SB + E_PB`
- `energy_total = energy_kinetic + energy_hartree + energy_xc + energy_external
  (+ energy_nonlocal + energy_ion)`

**File placement (ADR-0007).**
- Notebook + builder + tests + combined CSVs + kernel:
  `ResearchProject/systems/localised_jellium/hypotheses/stopping_from_decomposition/`.
- If the kernel is promoted to a library function: `inqview.analysis` (numpy-only,
  deps-clean) with its unit test in `inq-stack/python/tests/`.
- One executed `stopping_from_decomposition.ipynb` (all four phases) is the spine.

**Data sources (Phase 3, user-assisted at that phase's start).**
- Primary: `classical-highdensity-sv` Phase-3 sweep (6 velocities, classical, full
  ledger, clean plateau) + its later matched WP run(s).
- Plus any other select decomposition+plateau runs the user names. Run
  identification is DEFERRED to the start of Phase 3 (user in the loop there only).
</resolved_decisions>

<guard_rails>
- **Formula human-gate (Phase 1, hard).** Phases 2–4 do NOT begin until the user
  has approved the decomposition-S formula(e) in the notebook. This is the
  campaign's one mandatory human checkpoint (plus the Phase-3 run pick).
- **Self-test gate (Phase 2, numeric).** The classical decomposition-S MUST equal
  the Definition-2 deposit `E_absorbed/L_slab` on a reference run to within a stated
  tolerance (target ≲ few %); ledger closure MUST hold to ~1e-9. If not, the
  formula/implementation is wrong — fix before Phase 3, do not proceed.
- **Gauge gate (per WP run).** The gauge test Δ(E_SS,E_SB,E_BB)≈0 must pass before
  any pair-term-based WP number is reported; if a real gauge is present, restrict
  to gauge-clean combinations and say so.
- **Plateau requirement.** Only runs whose `E_electronic` has plateaued (deposit
  complete; |dE/dt| over the final 15% < 5% of the deposit) are admissible for the
  aggregate; a non-plateaued run's S is a lower bound and is flagged, not averaged
  in silently.
- **NEVER** report the projectile-KE-drag S as the WP headline or transfer it to
  the WP (`feedback_quantum_stopping_not_from_projectile_ke`). B4 (WP centroid-drift
  KE) is a labelled diagnostic overlay only, never a results-table "S_WP".
- **NEVER** `np.fft.fftshift` a VTI (`vti-coordinate-mapping` rule); load via
  `inqview.load_vti`. Density GIF at the TOP of any per-run notebook produced
  (`notebook-density-gif` rule).
- **Lindhard/linear-response is comparison-only, NEVER a gate**
  (`feedback_fourier_loss_function_gate`).
- **E_sp treatment (Phase 1 §C0, then Phase 4).** C1 drops `E_sp`; that is exact
  only if `E_sp(t_final)≈E_sp(0)`. MEASURE both ends before asserting it; the
  correction convention (drop / offset-remove initial E_sp / keep) is user-approved
  at the gate, never silently chosen. The Phase-4 E_sp-subtraction recompute is a
  ROUGH probe — its verdict (does it move S toward Lindhard?) is reported honestly,
  confirm OR refute, never retried into a confirm.
- **PROVISIONAL caveats:** the approved formula is not locked until Phase 1's human
  gate; the select run-set is user-identified at Phase 3; the classical benchmark
  data depends on `classical-highdensity-sv` completing its Phase 3 (open
  dependency — if that sweep is not yet on disk, Phase 3 proceeds on whatever
  select plateaued+decomposed runs the user names, and the S(E)/S(v) aggregate is
  built from those).
</guard_rails>

<tasks>
1. **Phase 1 — understand + propose formulae (MANUAL GATE).** In the one notebook:
   - **§A (Sonnet):** enumerate EVERY energy recorded in a selected reference run —
     lumped stores, pairwise ledger, projectile/track columns — grouped by CSV file,
     noting which exist in the WP run vs the classical run.
   - **§B (Sonnet):** show how they compose — load the reference run and VERIFY the
     closure relations numerically (classical: E_hartree=E_SS, E_external=E_SB+E_PS;
     WP: E_hartree=E_SS+E_PS+E_PP, E_external=E_SB+E_PB; and E_total sum), printing
     the residuals (~1e-9 expected).
   - **§C0 — the E_sp question (Sonnet measure + Fable 5 interpret).** C1 excludes
     the slab–projectile interaction `E_sp ≡ e_ps`. That exclusion is EXACT only when
     `E_sp(t_final) ≈ E_sp(0)`. MEASURE both ends on a reference run (and plot E_sp(t))
     to show whether it holds; the user's premise is E_sp(0) is substantial while
     E_sp(t_final) ≈ 0 (projectile absorbed). State the treatment (drop / offset-remove
     the initial E_sp / keep) that the user approves at the gate. Feeds C1 and §D and
     Phase 4.
   - **§C (Fable 5):** present the ranked candidate S-from-decomposition formulae
     (C1 target-absorption `T_slab+E_ss+E_sb+E_xc`; C2 classical projectile-KE;
     C3 WP projectile-loss `T_WP+E_pp+E_ps+E_pb` with ZPE/localisation subtracted),
     each with a plain-text formula, which representation it applies to, its physical
     meaning, gauge/conservation status, failure modes, the numeric extraction recipe,
     and the explicit WP caveat (KS orbitals cleanly separate slab vs WP); full A–E
     menu as an appendix.
   - **§D — setup-simplification brainstorm (Fable 5).** Propose setup modifications
     that simplify the energy structure for cleaner S (e.g. a long run that fully
     absorbs the projectile so `E_total(t_final)=E_slab(t_final)` and `E_sp(final)≈0`;
     other analogous simplifications). Feeds Phase 6.
   - **Done:** the user reviews §A/§B/§C0/§C/§D and APPROVES the formula(e) + the E_sp
     treatment + the setup ideas. (composes: twin-run-analysis, literature-review; no GPU.)
2. **Phase 2 — implement + validate the kernel (code-test + formula-validation).**
   Lock the approved formula into a reusable extraction function (skill-local under
   the hypotheses folder, or `inqview.analysis` if promoted), **exposing the E_sp term
   switchably** so C1's include/exclude convention is a parameter. Write known-case
   tests: (a) classical decomposition-S == Definition-2 deposit `E_absorbed/L_slab` on
   a reference run within tolerance; (b) ledger closure ~1e-9; (c) a synthetic fixture
   with an analytic answer. Dispatch a `formula-validation` agent (given ONLY the
   formula + its source) and reconcile with the user. Add a catalogue row.
   **Done:** tests green, formula-validation + user agree → formula locked.
3. **Phase 3 — apply to the select runs (USER-ASSISTED, then autonomous).** At the
   START of this phase, work WITH the user to identify the select decomposition+
   plateau runs (classical-highdensity-sv sweep + WP + any others named). For each:
   apply **C1** (classical + WP-with-caveat), **C2** with the FULL plotted chain
   (ΔE(t) plot → ΔKE_proj(t) plot → stated fitted region → best-fit line shown →
   value), and **C3** (ZPE/localisation subtracted). Report the **C1-vs-C2 classical
   agreement** explicitly (they should coincide). Run per-run sanity — energy
   conservation (E_conserved flat), the gauge test, plateau confirmation — recording
   S with its ds, mean velocity, and flags. **Done:** S computed and sane for every
   selected run; anomalies surfaced to the user.
4. **Phase 4 — A1 discrepancy investigation (Fable 5 + experiments).** Diagnose why
   the existing deposit-based S(v) (Definition 2, `E_absorbed/L_slab`) comes out
   **much larger than the Lindhard bulk-jellium prediction**. Steps: (i) FIND the
   well-behaved plateaued-WP localised-jellium runs where this S(v) was extracted
   (A1-T1); (ii) the ROUGH E_sp probe — recompute `E_absorbed` with the initial
   slab–projectile interaction removed (`D_corr = D_raw + E_sp(0)`, since
   `E_sp(t_final)≈0`) and compare the resulting S to Lindhard; (iii) test the three
   hypotheses with an experiment each — **(a)** localised jellium ≠ bulk at the same
   density (so Lindhard is the wrong reference → compare against a classical
   projectile instead), **(b)** the baseline is wrong (use the decomposition-sum,
   not `E_total`, whose difference may carry terms we should not track), **(c)** the
   CAP distorts the total-energy curve. Use Fable 5 to reason through (a)/(b)/(c).
   **Done:** each hypothesis has an experiment and a stated conclusion; the corrected
   vs Lindhard comparison is reported honestly (confirm or refute).
5. **Phase 5 — aggregate + plot (autonomous).** Group the per-run S by system
   (σ_WP/r_s/slab); build the headline **S(E0)** and companion **S(v0)** figures with
   classical and WP overlaid (canonical theme, `report-figures`/`scientific-figures`);
   overlay the Lindhard/linear-response curve as a non-gating eyeball; write
   `S_decomposition_summary.csv`.
6. **Phase 6 — synthesis + setup recommendations (autonomous).** Fold the Phase-4
   conclusion and the §D setup ideas into a short "how to run the next batch for a
   clean S" recommendation at the end of the notebook. **Done:** the single notebook
   executes end-to-end (0 errors), figures + summary CSV in the hypotheses folder;
   handover + frontmatter `done` flags + INDEX updated. (composes: notebook-making.)
</tasks>

<rules>
- ALWAYS compute classical and WP S by the SAME functional of the SAME columns
  (apples-to-apples); if a definition cannot be applied identically to both, label
  it single-representation and do not overlay it on the shared S axis without a caveat.
- ALWAYS carry the ds (path-length) choice with every reported S.
- ALWAYS pass the gauge test before quoting any pair-term-based WP number.
- NEVER launch a GPU run from this campaign — it is analysis-only; if new data is
  needed, that is `classical-highdensity-sv`'s job or a separate campaign.
- NEVER let the Lindhard reference become the reported S; comparison only.
- NEVER report projectile-KE drag as the WP headline.
- Use Fable 5 for formula brainstorming/derivation; Sonnet for data-fetch + simple
  per-run analysis (identity block).
- Everything lands in ONE notebook (`stopping_from_decomposition.ipynb`).
</rules>

<preflight>
- [ ] Intent self-contained: falsifiable hypothesis + success (classical
      decomposition-S reproduces Definition-2 deposit within tol; WP formula
      gauge-clean; coherent S(E0)/S(v0)) / failure (no admissible combination, or
      irreducible gauge dependence) criteria; each task has a done-criterion.
- [ ] Setup reproducible, zero guessing: consumes existing plateaued+decomposed
      runs (classical-highdensity-sv Phase-3 + user-named select runs); closure
      relations + file refs stated; formula shortlist + ds choices + aggregation
      axes locked; ONE notebook; file placement per ADR-0007. NO GPU runs.
- [ ] New code pre-gated: the Definition-1 extraction kernel is DERIVED → routed
      through code-test + formula-validation + catalogue row in Phase 2 BEFORE it
      produces any headline S in Phases 3–4.
- [ ] Validation & guard rails: Phase-1 formula HUMAN GATE; Phase-2 numeric
      self-test (classical decomposition-S == E_absorbed/L within tol; closure
      ~1e-9); per-WP-run gauge gate; plateau requirement; Lindhard comparison-only;
      PROVISIONAL caveats (formula not locked pre-gate; run-set user-identified at
      Phase 3; classical data depends on classical-highdensity-sv) named.
- [ ] Autonomous mechanics: Phases 2, 4, 5 & 6 fully autonomous; Phase 1 (formula +
      E_sp treatment + setup ideas) and Phase 3 (run pick) are the two human
      touchpoints; model routing Fable/Sonnet per identity; single notebook
      auto-built/executed; handover pointer present; agent updates handover + flips
      frontmatter done/status. (No dispatcher/GPU probe needed — analysis-only.)
- [ ] Grounding: stopping-power definitions = Correa 2018
      (docs/sources/correa-2018-electronic-stopping-power.md) + twin-run-analysis;
      closure/gauge = interaction_energies.hpp:17-18 +
      reference_twin_pairwise_decomposition + reference_charged_cell_hartree_convention;
      WP-KE prohibition = feedback_quantum_stopping_not_from_projectile_ke.
</preflight>
