---
id: tdhf-wp-orbital-approx
area: td-hf
title: "Is the KS-WP orbital a good approximation to the HF orbital?"
status: blocked
hypothesis: "The KS wavepacket orbital approximates the corresponding HF orbital; the fidelity F(t)=|<psi_HF|psi_KS>|^2 and physical-observable gaps quantify the limit of physical interpretability of the KS-WP orbital and the size of its self-interaction error (SIE) — tested in free propagation, a high-density jellium slab, and coronene."
handover: docs/handovers/td-hf-wp-orbital-approx.md
tasks:
  - { name: "orbital_fidelity kernel (inqview.analysis) + code-test + formula-validation + catalogue row", done: false }
  - { name: "TD-HF RT feasibility + cost pilot (size the small high-density slab)", done: false }
  - { name: "LDA exchange-dominance check at slab density (r_s~=4)", done: false }
  - { name: "Phase A - free propagation: 3 arms + analytic, metrics, notebook", done: false }
  - { name: "Phase B - small high-density slab: 3 arms, metrics, notebook (gated on slab validation)", done: false }
  - { name: "Phase C - coronene: 3 arms, metrics, notebook (caveated)", done: false }
  - { name: "Synthesis notebook - pin thresholds, SIE/interpretability bounds, verdict", done: false }
blocked_reason: "Phase B/C gated on locjel slab validation + TD-HF pilot (no in-repo precedent); thresholds deferred. Phase A + orbital_fidelity kernel + pilot are runnable now."
---

# Is the TDDFT KS Wavepacket orbital a good approximation to the HF Orbital?

<identity>
You are a scientific computing researcher working on first-principles simulations.
You understand the first-principles domain, write scientific-standard code, and
adhere to the rules, principles, and workflows established in this repository.
`inq/` is immutable; engine changes live only in `inq-study/`; new library code
lives in `inqkit`/`inqview`.
</identity>

<!-- ============================================================== -->
<!-- STAGE 1 — FRAME (locked <date>). Stages 2–5 below are PENDING.  -->
<!-- ============================================================== -->

<frame>

## The question

Kohn–Sham orbitals are mathematical constructs with **no individual physical
meaning**; Hartree–Fock orbitals, being the single-particle orbitals of a Slater
determinant, **do** carry physical meaning. When we inject a Gaussian wavepacket
as an *extra KS orbital* and evolve it under mean-field TDDFT, **how well does
that KS-WP orbital approximate the corresponding HF orbital?**

The gap between them bounds (a) **how far the KS-WP orbital may be physically
interpreted**, and (b) **the size of the self-interaction error (SIE)** in the
mean-field KS treatment of the injected electron — directly informing a later
SIE-removal refinement of the projectile simulations.

## Hypothesis (falsifiable)

The KS-WP orbital is a good approximation to the HF orbital — near-exact in free
propagation, degrading under interaction — and the measured gap is the SIE /
interpretability bound.

## Decision this informs

Whether KS-WP observables (centroid track, momentum-loss / stopping, density,
wake) can be quoted as physically meaningful as-is, or whether an SIE correction /
exact-exchange treatment is required first — and roughly how large that correction
is.

## Comparison metrics (LOCKED; thresholds deferred)

Orbitals are not observables (arbitrary global phase), so all metrics are
**phase-invariant**. Compute and report each as a function of time:

1. **Orbital fidelity** — `F(t) = |⟨ψ_HF(t) | ψ_KS(t)⟩|²` (headline; 1.0 = identical).
2. **Density agreement** — L1 and L2 norm of `|ψ_KS(t)|² − |ψ_HF(t)|²`.
3. **Physics gap** — divergence in the **centroid track** and the **mean
   momentum-loss** (the stopping-relevant quantity), reported in eV against the
   ~7 eV SIE scale.

> **OPEN (must be pinned before the campaign runs):** the pass/fail thresholds
> (e.g. free-propagation sanity band, "good-approximation" band under
> interaction). Deferred by user 2026-06-22; compute first, set bands later.

## Scope — three phases

- **Phase A — Free propagation** (WP alone, vacuum). The *cleanest SIE probe*: a
  single electron has HF = exact one-electron orbital, while a single KS electron
  in LDA still feels its **own** Hartree+xc, so any `F(t)<1` here is **pure SIE**
  (no bath, no correlation). Also the build/validation testbed for the TD-HF tool.
- **Phase B — High-density jellium slab** (~Na, r_s≈4). Primary interaction case,
  chosen **because exchange dominates correlation** there, so HF (exchange-only,
  no correlation) is physically appropriate.
- **Phase C — Coronene** (central molecular target). Included per user; **carry
  the caveat** that a molecular target is *not* exchange-dominated, so HF's missing
  correlation is most questionable here — coronene results are interpreted with
  that flag (separates "is the KS orbital good?" from "is HF even valid here?").

**Excluded:** pure jellium N=234 (L=50 cubic) — low density ⇒ correlation too
large for an exchange-only HF reference to be meaningful.

## Hard dependencies / blockers (this campaign must NOT run until resolved)

- **Jellium-slab validation** — Phase B reuses the slab from `locjel-campaign`;
  its implementation must be validated first (`locjel-campaign` Phase 2). Until
  then this campaign is `blocked` on it.
- **TD-HF tooling** — INQ is a DFT engine; a TD-HF reference orbital must be
  obtained (external package or implementation). Choosing/integrating it is a
  Stage-2/3 task and a feasibility gate.
- **Exchange-vs-correlation check** — verify LDA exchange dominates correlation at
  the slab density (~Na, r_s≈4) before trusting Phase B.

</frame>

<!-- ====================== PENDING STAGES ======================= -->

<observables_set>
Per **arm** (full-LDA, exchange-only-LDA, HF), in EACH phase:
- **Complex WP orbital `ψ_WP(r,t)`** dumped via inqkit `ComplexField3DWriter` at a
  **matched 150-frame cadence**, on the **identical grid** across all three arms
  (arms differ only in theory ⇒ grids coincide ⇒ overlaps need no interpolation).
  This is the data the fidelity metric consumes.
- Standard **ADR-0006 jellium-WP set**: energy components (per step), density VTI,
  current, WP momentum distribution `n(k)`, COD track.

Derived metrics (post-processing, all phase-invariant), per arm-pair:
- **`F(t) = |⟨ψ_a(t)|ψ_b(t)⟩|²`** — orbital fidelity (NEW `orbital_fidelity` kernel).
- **density L1/L2** of `|ψ_a|² − |ψ_b|²` (same kernel).
- **centroid-track** and **mean-momentum-loss** divergence (existing
  `center_of_density` / `momentum_distribution` kernels), reported in eV vs the
  ~7 eV SIE scale.
Pairwise reads: HF−(x-only-LDA) = exchange-treatment/SIE; full-LDA−(x-only-LDA) =
correlation; full-LDA−HF = total gap carried by production KS-WP.
</observables_set>

<resolved_decisions>
<!-- STAGE 3 verified facts (line-refs from inq/; numbers from literature). -->

<engine_feasibility kind="verified">
RT-TDHF is supported IN-ENGINE (inq-study); no external code needed.
- Theory selector `theory.hartree_fock(coeff=1.0)` (inq/src/options/theory.hpp:74),
  exchange-only LDA via `functional(XC_LDA_X)` (theory.hpp:161), full LDA `lda()`.
- Exact exchange in RT uses **ACE acceleration** (exchange_operator.hpp:49,88,105).
- **Spin-polarised exact exchange SUPPORTED**: `spin_config::POLARIZED`
  (states/spin_config.hpp:20), requested via `electrons.spin(POLARIZED)`; exchange
  is same-spin via the kpin structure (exchange_operator.hpp:133-160).
- **WP injection is theory-agnostic**: `inqkit::WavePacket` +
  `extra_electrons(N)` (wavepacket.hpp:119-172) — no LDA-only / ETRS-only / real-
  potential assumption; works under HF + CN.
- **NO precedent anywhere in this repo** for hartree_fock()/pbe0()/b3lyp() in any
  run/test (only `.lda()` is ever used; sole RT-CN example is inq/tests/h2o.cpp:103).
  ⇒ the TD-HF feasibility+cost pilot (task 2) is a HARD GATE, not a formality.
</engine_feasibility>

<propagator kind="locked">
**Crank-Nicolson for ALL THREE arms.** Forced: ETRS asserts no exact exchange
(real_time/etrs.hpp:26 `assert(not sc.has_exact_exchange())`); only CN handles it
(crank_nicolson.hpp:125,173). LDA arms also use CN for apples-to-apples. This
INVERTS the usual jellium WP guidance (ETRS, to protect CAP absorption,
[[reference_inq_propagator_mask_absorber]]) — benign here: no CAP, so CN's
norm-conserving unitary step is correct. dt = 0.02 a.u. (matches production jellium
runs); pilot-confirm CN-TDHF stability. See memory
[[reference_inq_rt_tdhf_requires_crank_nicolson]].
</propagator>

<slab kind="grounded">
Phase B reuses the localised-jellium slab DENSITY (r_s=3.996, Na-like) but SMALL.
Reference (shared/configs/slab_n234_L50.hpp): L=50 cubic, slab 25 Bohr thick
(half-width 12.5, full 50x50 face), N=234, n0=3.744e-3 a0^-3, dx=0.50 Bohr.
Small slab = shrink the volume at FIXED density so N_occ ~ few dozen; exact size
set by the TD-HF cost pilot (task 2). Gated on locjel slab-implementation
validation (`locjel-campaign` Phase 2).
</slab>

<coronene kind="grounded">
Phase C reuses the converged GS checkpoint
`ResearchProject/systems/coronene/checkpoints/tsubonoya_2014_paper_replica_gs`
(C24H12, ALDA, cutoff 54 Ha, box 34.77x34.77x59.90 Bohr, 8 extra states). Load via
`electrons.load(...)` with matching cutoff/extra_states. OVERRIDE its native WP
(sigma=1.0/E=200) with this campaign's locked sigma=0.5 / E=100 eV.
</coronene>

<exchange_correlation kind="grounded">
At r_s=4 (HEG, per electron, Ha): Dirac exchange eps_x = -0.4582/r_s = -0.115;
PW92 correlation eps_c ~= -0.045 (pin exact value in a docs/sources note — PW92
evaluated directly). |eps_x|/|eps_c| ~= 2.5: exchange dominates but correlation is
~28% of exchange, NOT negligible. This is WHY the three-arm decomposition is
needed: full-LDA - exchange-only-LDA MEASURES the correlation contribution rather
than assuming it away. (CAUTION on units: the textbook -0.9163/r_s is RYDBERG, not
Hartree.)
</exchange_correlation>

<scientific_grounding kind="verified">
- It is ESTABLISHED that KS and HF *occupied ground-state* orbitals are very close
  (Della Sala & Goerling). The OPEN question — this campaign's contribution — is
  whether that near-equivalence survives for a **real-time-propagated injected
  scattering orbital** where SIE acts dynamically. Reframes the campaign from
  "confirm the obvious" to "test the static result in the dynamic non-equilibrium
  regime".
- HF is the **SIE-free reference, NOT the 'more accurate' theory** (it omits
  correlation; TDDFT can outperform TDHF). The campaign measures a *gap*, not a
  quality ranking.
- SIE grounding: Perdew-Zunger SIC; "Orbital Anatomy of Self-Interaction in KS-DFT"
  (arXiv:2407.09680). ACE cost: Lin Lin, JCTC 2016.
- `docs/sources/` notes to author (via `literature-review`) for: KS-vs-HF orbital
  near-equivalence; PZ-SIC + SIE-in-RT-TDDFT; PW92 eps_c(r_s=4); ACE.
</scientific_grounding>

<file_placement kind="locked">
Runs live in their HOME systems (ADR-0007 sweep layout), NOT a new system:
- Phase A (free WP) -> `systems/vacuum/tdhf_free/<arm>/`, scripts in
  `systems/vacuum/scripts/tdhf_free/`; per-phase notebook in
  `systems/vacuum/hypotheses/tdhf_free/`.
- Phase B (slab)   -> `systems/localised_jellium/tdhf_slab/<arm>/` (+ scripts/,
  hypotheses/tdhf_slab/). Reuses the locjel slab GS.
- Phase C (coronene) -> `systems/coronene/tdhf_coronene/<arm>/` (+ scripts/,
  hypotheses/tdhf_coronene/). Reuses the coronene GS checkpoint.
- `<arm>` in {full_lda, xonly_lda, hf}.
New kernel `orbital_fidelity` -> `inqview/analysis/` (deps-clean) + test in
inq-stack/tests/python/.
CROSS-SYSTEM synthesis notebook (spans all three systems) -> `docs/reports/
td-hf-orbital/` (manuscript-level cross-system home per ADR-0007 / file-placement).
</file_placement>
</resolved_decisions>

<guard_rails>

**HARD GATES (ordering — no expensive HF run starts before these pass):**
1. `orbital_fidelity` kernel written + tests green + `formula-validation` agent +
   catalogue row — BEFORE any HF run.
2. TD-HF feasibility + cost pilot PASS (below) + slab sizing fixed — BEFORE Phase B.
3. locjel slab-implementation VALIDATED (`locjel-campaign` Phase 2) — BEFORE Phase B.
4. exchange-dominance ratio recorded (`docs/sources/` note) — BEFORE trusting Phase B.
5. coronene GS confirmed loadable — BEFORE Phase C.
6. Order: Phase A → B → C.

**PILOT GATE (task 2 — first-ever HF+CN+WP+spin-polarised run in this repo; Phase A
free, ~50–100 CN steps).** PASS iff ALL:
- completes, no NaN / non-finite values;
- total energy REAL and stable; exact-exchange energy term finite and nonzero
  (`has_exact_exchange()` true ⇒ ACE active);
- CN conserves norm: |‖ψ_WP‖ − 1| < 1e-3 across the pilot;
- s/step measured ⇒ size the Phase-B small slab so a full HF run fits the GPU-time
  budget (target ≤ ~24 GPU-h/arm; tune).
FALLBACK if HF is infeasible at a physically-meaningful slab: shrink further; if
still infeasible, restrict to Phase A and report HF-in-slab as COST-BLOCKED.

**ORBITAL-IDENTITY RULE (keeps F(t) well-defined under interaction).** The WP is the
designated extra-state slot; compare that slot across arms (identical ordering).
Identity guard each frame = max single-orbital overlap of the live WP with the t=0
injected Gaussian: ≥ 0.7 ⇒ single-orbital F(t) meaningful; < 0.7 (WP hybridises /
scatters into several orbitals) ⇒ single-orbital fidelity UNRELIABLE there — fall
back to the WP DENSITY comparison (|ψ|² via the slot) + subspace fidelity, and FLAG
the window in the notebook.

**ABORT / SAFETY:** NaN / non-finite or complex total energy beyond tol → abort+log;
norm drift |‖ψ‖−1| > 1e-2 → abort/flag (CN should conserve); GS-HF SCF
non-convergence → abort; GPU occupied by another user → warn (`cudaMemGetInfo`
probe; NVML broken); GPU is the default.

**BOUNDARY / CADENCE:** 4σ/1σ launch-stop rule (free + slab); complex WP orbital +
density at matched **150-frame** cadence; identical grid across arms (assert box/dx
equal before comparing).

**PROVISIONAL:** no in-repo precedent for `hartree_fock()`+`crank_nicolson()`+
injected-WP+spin-polarised ⇒ ALL HF numbers PROVISIONAL until the pilot validates
the path; verify `inq-study ≡ inq` (diff) on non-exchange files for attribution.
Thresholds DEFERRED ⇒ verdict (task 7) is descriptive until pinned.
</guard_rails>

<run_matrix>
**3 phases × 3 arms = 9 runs**, identical setup within a phase, only the theory
(`inq-study` `theory.lda()` / `lda` exchange-only / `theory.hartree_fock()`)
changing between arms. σ=0.5 Bohr, E=100 eV, single energy, everywhere.
Spin-polarised + collinear, WP in the majority (↑) channel, identical across arms.

| Phase | System | Arms | Notes |
|---|---|---|---|
| A — free | 1-electron WP in vacuum, box per 4σ/1σ rule | full-LDA, x-only-LDA, HF | + **analytic free-Gaussian** cross-check (HF of 1 e = exact). Pure-SIE decomposition; cheapest; TD-HF testbed. |
| B — slab | **small** high-density jellium slab, r_s≈4, N_occ ≈ few dozen (sized by the cost pilot) | full-LDA, x-only-LDA, HF | **Gated on `locjel` slab-implementation validation**. Exchange-dominated ⇒ HF appropriate. |
| C — coronene | reuse coronene GS + WP params | full-LDA, x-only-LDA, HF | Run **last**; HF on C₂₄H₁₂ heavy; molecular ⇒ HF-correlation caveat carried. |

Cost: RT-TDHF ~O(N_occ²) FFT-pairs/step. The **HF arm's pilot s/step sets the
feasible system size/time** per phase. GPU default; `cudaMemGetInfo` probe.
</run_matrix>

<tasks>
1. **`orbital_fidelity` kernel** — `inqview.analysis` (deps-clean numpy):
   `F(t)=|Σ ψ_a* ψ_b dV|²` + density L1/L2. Formula-bearing ⇒ `code-test`
   (identical→1, orthogonal→0, global-phase invariance) + `formula-validation`
   agent + catalogue row. **Done-criterion:** kernel + tests green + catalogue row,
   BEFORE any HF run.
2. **TD-HF feasibility + cost pilot** — short inq-study RT run with
   `theory.hartree_fock()` (Phase A first): energy real, no NaN, exchange operator
   active; measure s/step ⇒ size the Phase-B small slab. **Done:** pilot passes +
   slab size fixed.
3. **LDA exchange-dominance check** — confirm |E_x| ≫ |E_c| at r_s≈4 (numeric +
   `docs/sources/` note). **Done:** ratio recorded and gate passed.
4. **Phase A** — 3 arms + analytic; compute F(t)/density/centroid/momentum; per-phase
   notebook. **Done:** notebook executed (0 errors), metrics reported.
5. **Phase B** — *gated on `locjel` slab validation*; 3 arms; metrics; notebook.
6. **Phase C** — *gated on coronene GS*; 3 arms; metrics; notebook; HF caveat stated.
7. **Synthesis notebook** — cross-phase SIE / interpretability bounds; **pin the
   deferred pass/fail thresholds**; verdict on the hypothesis. **Done:** final
   notebook executed, thresholds stated, verdict recorded.

Each completed task ⇒ flip its frontmatter `done` flag + update the handover.
</tasks>

<rules>
ALWAYS:
- All three arms identical except `theory` (and exchange coeff); assert
  box/grid/dt/WP identical before any cross-arm comparison.
- Crank-Nicolson for ALL arms; spin-polarised collinear, WP in the majority channel.
- Flip the frontmatter `done` flag + update the handover at each task completion.
NEVER:
- Never use ETRS with exact exchange (engine asserts; see memory
  `reference_inq_rt_tdhf_requires_crank_nicolson`).
- Never start an expensive HF run before the kernel + pilot gates pass.
- Never edit `inq/`; any engine work goes only in `inq-study/`.
- Never quote single-orbital F(t) inside a window where the identity guard failed
  without the flag.
</rules>

<preflight>
Re-verify from THIS prompt before burning GPU:
- [ ] **Intent** — hypothesis + metrics {F(t), density L1/L2, centroid, momentum-loss}
  present; thresholds DEFERRED (verdict is descriptive until pinned in task 7).
- [ ] **Setup** — σ=0.5 / E=100 eV; CN all arms; dt=0.02; spin-polarised collinear,
  WP majority channel; grids identical across arms; slab = locjel density, small &
  pilot-sized; coronene GS path set.
- [ ] **New code pre-gated** — `orbital_fidelity` kernel + tests + formula-validation
  + catalogue row BEFORE any HF run.
- [ ] **Validation/gates** — pilot PASS (no NaN, real energy, ACE active,
  |‖ψ‖−1|<1e-3) + slab sizing BEFORE Phase B; locjel slab validated BEFORE Phase B;
  exchange-dominance recorded; orbital-identity guard active; abort conditions set.
- [ ] **Autonomous mechanics** — GPU default (`cudaMemGetInfo` probe; NVML broken);
  per-phase Gmail; per-phase + cross-system synthesis notebooks auto-built
  (`analyse.py`/dispatcher tail); agent updates handover + frontmatter done/status.
- [ ] **Grounding** — HF = SIE-free reference (NOT "more accurate"); ε_x=−0.4582/r_s
  verified; KS≈HF *static* known, the *dynamic injected-orbital* regime is the open
  question; `docs/sources/` notes owed.
- [ ] **HARD** — no repo precedent for HF+CN+WP+spin-pol ⇒ the pilot is the go/no-go.
</preflight>
