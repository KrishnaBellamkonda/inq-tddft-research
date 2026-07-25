---
# ROUGH DRAFT — authored interactively 2026-06-27 to capture intent only.
# NOT autonomy-ready. To be built out carefully later via the `campaigns` skill.
# This campaign is the GATE that unblocks `classical-projectile-fix` (Campaign 1)
# via threads B (SIE decomposition) + D (long-range cutoff).
id: localised-jellium-gs-study
area: localised_jellium
title: "Localised jellium GS — parameter study + analytical mental models"
status: paused
hypothesis: "The localised jellium slab's ground state and its response to a projectile can be understood from simple building blocks: (i) a charged-plates/capacitor mental model of the GS density; (ii) a clean decomposition E_total(t=0)−E_GS−T_WP = E_SIE + E_WP–jellium-repulsion, with E_SIE the part that survives box elongation; (iii) a defensible long-range cutoff for the classical projectile; and (iv) parameter sweeps (Lx, Ly, Lz, w) that map how surface effects change — together giving the intuition needed to fix the classical projectile (Campaign 1)."
handover: docs/handovers/localised-jellium-gs-study.md
tasks:
  - { name: "Worksheet B+C+D inputs (material + .docx + arXiv PDFs + md/XML plan) → external agent builds worksheet", done: false }
  - { name: "Pre-gate code: Phi/sigma_s extraction + ghost-background integral (code-test + formula-validation + catalogue)", done: false }
  - { name: "H0 base WP-vs-classical E_total(0) gap (matched sigma, r=4 & 40, L_z=120 w=0); notebook", done: true }
  - { name: "H1 edge model — clean w (Gibbs vs Friedel); notebook", done: true }
  - { name: "H2 GS convergence + open-z viability (Lz sweep, periodicity-2 GS smoke); notebook", done: true }
  - { name: "H3 surface energetics — sigma_s, e_bulk, charged-plates (thread A); notebook", done: true }
  - { name: "H4 WP energetics — PBC-vs-open-z verdict + E_SIE route 1 + k0 control; notebook", done: true }
  - { name: "H5 classical subtraction — route 2 + thread D cutoff; pins E_SIE/reference for Campaign 1; notebook", done: true }
  - { name: "Cumulative campaign notebook (highlight plots, overall story) + synthesis (charged-plates, zero-point E)", done: true }
blocked_reason: "H0-H5 all executed + emailed (audit clean, path bug fixed); accepted as-is 2026-06-27. PROVISIONAL follow-ups documented (not execution errors): (1) open-z periodicity-2 net-charge G=0 reference -> per-2 E_SIE biased (PBC E_SIE=4.3 eV is sound); (2) H2 work-function Phi; (3) H3 sigma_s E_self correction; (4) H5 ghost-background integral. Plus worksheet (external agent), pre-gate code, cumulative notebook + synthesis."
---

# Localised jellium GS — parameter study + analytical mental models

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository.
</identity>

<rough_draft_banner>
This is a ROUGH DRAFT capturing the user's intent + the simple questions resolved
2026-06-27. The purpose of the run-sets is to BUILD UNDERSTANDING (analytical
mental models), per the campaign title. Do NOT execute as-is. "(deferred to the
/campaigns pass)" marks intent, not locked decisions.
</rough_draft_banner>

<gate>
**This campaign unblocks Campaign 1** (`docs/campaigns/localised_jellium/classical_projectile_fix.md`).
Its critical-path deliverables for that are **thread B** (the energy reference
`E_jellium(t=0)=E_total(0)−⟨T_WP⟩−SIE`, properly decomposed) and **thread D** (the
long-range cutoff for the classical projectile). The other threads deepen
intuition but do not gate Campaign 1.
</gate>

<existing_material>
- **A 472-line theory worksheet already exists:** `docs/notes/localised-jellium-theory.md`
  — HEG reference; background electrostatics (G=0/neutrality/self-energy);
  KS energy decomposition + sign audit; spherical-cluster + **Lang–Kohn slab**
  benchmarks (work function, surface energy 86.4 erg/cm²); interior-density
  convergence; brief projectile section; reference list. **EXTEND it — do not
  restart.** Missing: charged-plates mental model (A), SIE decomposition (B),
  Coulomb-vs-distance (C), long-range cutoff (D), run-mining (H).
- Existing source notes: `correa-2018-…`, `quijada-2007-…`,
  `nazarov-gross-2025-…`, `stopping-power-formulae`, `stopping-power-jellium-anchors`.
</existing_material>

<description>
**Purpose (user's words, lightly edited).** Start from small building blocks and
gain a strong intuition for what is happening in this localised jellium slab
system — making analytical **mental models** (e.g. charged plates / capacitors).
Then understand the dynamics of a classical projectile going through, again
building from small controllable pieces. Build simple systems, understand them
well, then add complexity step by step.

**Character (resolved 2026-06-27): a production parameter study** whose runs serve
the understanding. Two run-sets form the backbone; the analytical threads are the
interpretation layer wrapped around them.

**The threads (all from the mind-dump):**
- **A — charged-plates mental model.** Take the GS charge distribution; by its mean,
  assign +/− "plates": a net-positive smoothing region just outside the slab; a
  tall negative charge peak at the boundary; a net-positive interior; a net-negative
  boundary plate; a net-positive region beyond. Interpret the projectile's
  interaction as interaction with this alternating plate stack (the picture evolves
  in time but is a useful anchor).
- **B — SIE decomposition (critical path).** The current SIE estimate
  `E_total(jellium+WP) − E_GS(jellium) − T_WP` is **not quite right** — it omits the
  classical repulsion between projectile and jellium. Correct frame:
  `E_total(t=0) − E_GS − T_WP = E_SIE + E_WP–jellium-repulsion`, where **E_SIE is the
  part that does NOT vanish as the box is elongated**.
- **C — Coulomb-vs-distance.** Use classical electrostatics to estimate the
  WP–jellium Coulomb repulsion; then measure `E_total(WP) − E_total(classical)` for
  several WP–jellium distances r, in an **elongated-z box (same x,y)**, to map the
  interaction vs r and isolate the box-elongation-invariant SIE (this IS run-set 2).
- **D — long-range cutoff (critical path, owned here).** Decide a defensible cutoff
  for the classical projectile's radial Coulomb potential to avoid loop-around / PBC
  self-interaction; grounded in literature (below). Hand the prescription to
  Campaign 1.
- **E — localisation cost / zero-point significance.** Does the WP localisation cost
  (zero-point error) materially change the classical-vs-WP comparison? Assumed no —
  verify deeply.
- **F — config-parameter / surface effects.** Vary Lx, Ly, Lz (one at a time) and w
  (the erfc boundary-smoothing width); for each, see how surface effects change
  (this IS run-set 1).
- **G — literature + worksheet.** Collect localised/surface-jellium literature
  (user + agent read independently); extend the worksheet.
- **H — mine existing classical runs.** Scout existing localised-jellium classical
  VTI/observables for behaviour clues.
</description>

<literature>
Found 2026-06-27 (→ `docs/sources/` notes in task 1):
- **Cutoff / finite-size (thread D):** "Trajectory sampling and finite-size effects
  in first-principles stopping power calculations", *npj Comput. Mater.* (2023),
  arXiv:2307.03213 — periodic-image re-crossing into excited density; plasmon-cutoff
  finite-size error ~8%. **Coulomb-cutoff (truncated kernel)** + **Martyna–Tuckerman**
  nonperiodic Poisson solver — standard image-artifact removal (ORNL "Coulomb finite
  size effects"). "Examining RT-TDDFT for stopping power", arXiv:1805.01377.
- **Jellium-surface foundations (threads A/F/G):** Lang–Kohn (work function, surface
  energy, Friedel oscillations) — already cited in the worksheet Part 5; pin the
  finite/localised-jellium specifics in the careful pass.
</literature>

<observables_set>
(deferred to the /campaigns pass) Reuse the full suite + current cadence (ADR-0006).
Especially for understanding: GS **density profiles** (axial + radial), the
**energy decomposition** (Hartree / XC / kinetic / external / background self-energy
— worksheet Part 3.4), and density VTI for the charged-plates picture. For run-set 2:
E_total at t=0 for matched WP vs classical at each distance.

Note (2026-07-07): in the ΔE-component plots each energy term is plotted **relative
to the far reference**, e.g. the U_H curve is `U_hartree(r) − U_hartree(r=40)` (not the
absolute Hartree) — this isolates the r-dependent WP–slab part; the large r-independent
offset (slab self-energy, plus the periodicity-2 Poisson G=0 reference `0.5·rc²` that
flips the *absolute* Hartree/external signs vs periodicity 3) drops out.
</observables_set>

<run_sets>
- **Run-set 2 — WP-distance E_total(r) sweep (threads B+C+D).** σ_WP=0.5, jellium
  centered at z=0, baseline slab geometry (slab_n82: 50×50 face, 25 Bohr thick,
  N=82, r_s≈5.67) in a **bigger box L_z=120** (z∈[−60,+60], spacing 0.5). Measure
  `excess(r,BC) = E_total(t=0) − E_GS(BC) − ⟨T_WP⟩` (= E_SIE + E_cross(r)) vs WP–slab
  distance r, comparing **periodicity 3 (PBC, production default) vs periodicity 2
  (open-z)** — INQ supports open-z natively (poisson.hpp:184-213, 2d kernel opens z,
  Rozzi cutoff = thread D's cutoff in-engine).

  **Rung 1 (COARSE, locked 2026-06-27, session 3):** fresh slab GS in L_z=120 under
  periodicity 3 AND 2; r ∈ {4, 16, 28, 40} Bohr from the near face (r_min=8σ_WP);
  stationary WP k₀=0 (⟨T_WP⟩=81.6 eV zero-point), plus a k₀=2.71/100 eV control at
  r=4 under both BCs. 2 GS + 8 + 2 = 12 INQ runs (WP are 1–5-step). **Done:** verdict
  PBC-vs-open-z for production. Smoke-test the periodicity-2 GS FIRST (no project GS
  has used periodicity 2). E_total via WP inject + ≥1 step; excess is post-processing
  (no new kernel).

  **Route 2 — classical subtraction (LOCKED 2026-06-27 session 3):** matched ghost
  `electron_gaussian_wpsigma0p5.upf` (charge std σ_WP/√2=0.354), k₀=0, t=0 energetics,
  **full mirror of the coarse WP grid: r∈{4,16,28,40} × periodicity {3,2} = 8 classical
  runs.** **Ghost–background term `∫v_ghost·n₊` is MANDATORY** (re-added in analysis;
  launch-far alone fails — bare ghost–electron Coulomb decays only as ~N/r ≈ 56 eV at
  r=40; the chargeless ghost omits the background attraction). Then
  `E_SIE = E_WP − [E_cl + ∫v_ghost·n₊] − ⟨T_WP⟩`. **Three deliverables:** (i) E_SIE(r)
  cross-check vs the route-1 plateau (≈4.5 eV); (ii) E_cross^WP(r) == corrected
  E_cross^cl(r) (WP–slab interaction is purely classical); (iii) classical
  periodicity-3-vs-2 image error → feeds Campaign 1's cutoff (thread D). Ghost–bg
  integral pre-gated (formula-validation: numeric vs closed-form erf).

  **Rung 2+ (deferred):** finer r-sweep once the coarse BC verdict is in.

  **Coarse run-set 2 total:** 2 GS + 8 WP + 2 WP(k₀ control) + 8 classical = 20 INQ
  runs (WP/classical are 1–5-step; GS is the cost).

- **Run-set 1 — surface-physics sweep (thread F, locked 2026-06-27 session 3).**
  GS-only (no projectile), baseline slab_n82 (Lx=Ly=50, a=12.5, Lz=90, n₀=1.31e-3,
  spacing 0.5). **Metrics per config:** profile (interior n₀ vs HEG, Friedel
  amplitude+wavelength π/k_F≈9.3, spill-out length, surface dipole, E_GS) **+ work
  function Φ = V_vac−E_Fermi + surface energy σ_s** (liquid-drop fit over the
  thickness series; 86.4 erg/cm² target, worksheet Part 5).

  **Sweeps (one-at-a-time; Lx/Ly DROPPED as in-plane convergence checks):**
  - **edge-width w (FIRST):** {0, 0.5, 1.0, 1.5, 2.0} Bohr at baseline a,Lz →
    Gibbs-vs-Friedel discriminator + dipole/Φ sensitivity → **lock the smallest w
    that kills Gibbs while keeping π/k_F Friedel**, use it as baseline w for the rest.
  - **thickness a:** {7.5, 12.5, 17.5, 22.5} (2a=15/25/35/45, N≈50/82/114/148 at
    fixed n₀, N rounded even) → σ_s + e_bulk (liquid-drop) + interior-n₀ breakdown.
  - **Lz (vacuum):** {50, 70, 90, 120} (vacuum/side 12.5/22.5/32.5/47.5) → Φ plateau.

  **Pre-gate (required before the sweep):** Φ and σ_s extraction helpers via
  code-test + formula-validation + catalogue row (NEW code). ~11–13 GS runs total
  (baseline/Lz=90 + Lz=120 reuse existing/run-set-2 GS where applicable).
</run_sets>

<hypothesis_ladder>
**Execution as an ordered ladder of falsifiable hypothesis-sets (locked 2026-06-27
session 3, GS-foundation-first).** Simplest/most-informative first; each rung's result
feeds the next. The clean w (H1) is upstream of every later GS (incl. run-set 2's box);
the GS surface picture (H3) feeds the interpretation of E_cross in H4/H5. This refines
the earlier "critical-path-first" spine — the Campaign-1 gate (E_SIE + cutoff) now lands
at H4–H5 on top of a validated GS.

| # | Falsifiable hypothesis | Runs | Decides / feeds |
|---|---|---|---|
| **H0** | the base WP−classical `E_total(0)` gap equals the WP localisation energy (≈81.6 eV) | WP & matched classical (k₀=0) at r∈{4,40}, **L_z=120 box, w=0**, periodicity 3 (1 GS + 4 tiny runs) | **predicted FALSE — raw gap is ghost-bg-artifact-dominated (~−650 eV near), NOT localisation** → motivates the ghost-bg correction (H5); seeds run-set 2 |
| **H1** | finite erfc `w ≳ grid` kills Gibbs while keeping Friedel (λ=π/k_F≈9.3); below it ringing tracks the grid | w∈{0,.5,1,1.5,2}, GS-only, periodicity 3 | **clean w** for all later GS |
| **H2** | neutral-slab interior n₀/E are box- & BC-independent (periodicity 3≈2; Φ plateaus by Lz~90) ⇒ open-z usable | Lz∈{50,70,90,120} + periodicity-2 GS smoke, clean w | Lz, Φ, **open-z viability** |
| **H3** | E(N) liquid-drop-linear ⇒ σ_s (~86 erg/cm²), e_bulk→HEG; thin slabs lose the bulk | thickness a∈{7.5,12.5,17.5,22.5} | σ_s, e_bulk, **charged-plates (A)** |
| **H4** | image error excess(r,3)−excess(r,2) is sig/negligible ⇒ production BC; excess(r)→E_SIE≈4.5 eV; k₀-indep | WP E_total(r), k₀=0, {4,16,28,40}×{3,2} + k₀ control | **BC verdict** + **E_SIE route 1** |
| **H5** | corrected route-2 E_SIE matches route-1 ∀r; E_cross^WP=E_cross^cl; classical image mirrors WP | classical mirror {4,16,28,40}×{3,2} + ghost-bg | **pins E_SIE/reference for Campaign 1**; thread-D cutoff |

**Pre-gate positions:** Φ/σ_s code before H2/H3; ghost-bg integral before H5;
periodicity-2 GS smoke is the H2 gate for all open-z runs.
</hypothesis_ladder>

<notebook_contract>
**Per-hypothesis notebooks + one cumulative campaign notebook (locked; user is the
primary reader — "looking closely to learn as much as possible").** Per
`notebook-making` house narrative (context → formulas, every term defined →
full reconstructable setup → linked source files → results → takeaway).

- One `.ipynb` per rung: `H0_base_difference`, `H1_edge_model`, `H2_gs_convergence`,
  `H3_surface_energetics`, `H4_wp_energetics`, `H5_classical_subtraction` — each STATES
  its hypothesis up front, shows the outcome plots that prove/disprove it, and ends with
  an explicit verdict.
- `campaign_cumulative.ipynb` — stitches the **highlight plot** from each rung into the
  overall story (GS edge → bulk/surface → open-z → E_SIE → reference for Campaign 1).
- Placement per ADR-0007: `hypotheses/<NN_purpose>/` with the run-set's combined CSVs,
  builders, and figures; canonical theme (report-figures). Auto-build via the
  dispatcher / `analyse.py` tail where applicable.
</notebook_contract>

<autonomous_execution>
**Running autonomously, no user in the loop (user decision 2026-06-27).** Execute the
ladder strictly in order **H0 → H1 → H2 → H3 → H4 → H5 → cumulative**. Each phase:
build machinery (just-in-time) → run on GPU (probe first; GPU0 may be busy — use a
free card, warn if a run is another user's) → analyse → make the phase highlight plot
→ **email the user** → flip the frontmatter task `done: true` + update the handover →
launch the next phase.

**Executor = a PYTHON orchestrator, NOT bash** (user decision 2026-06-27; bash
autonomous scripts are brittle). Headless, idempotent-resume (skips completed runs),
per-phase try/except + failure emails, one-shot sim retry:
`scripts/campaign_autorun/orchestrate.py` (parametrised run.cpp in `gs/`,`wp/`,
`classical/`; per-phase analysis+email in `analyse_phase.py`). Launch:
`GPU=1 nohup venv/bin/python3 orchestrate.py &`.

**Per-phase email is MANDATORY — invoke the `email-notifications` skill.** Every phase
email MUST carry the highlight plot and the four parts: (1) the hypothesis being
tested, (2) what was done, (3) what the plot shows, (4) the conclusion from the phase
results. No phase is "done" until its email is sent.

**Pre-gates still gate their phases:** Φ/σ_s extraction before H2/H3; ghost-bg integral
before H5 (and used in H0's analysis). periodicity-2 GS smoke is the H2 gate.
Abort a phase on NaN / complex energy / GPU unavailable, email the failure, and stop
the chain for that branch.
</autonomous_execution>

<resolved_decisions>
Locked 2026-06-27:
- **Character:** production parameter study serving an understanding goal.
- **Gate role:** delivers threads B + D that unblock Campaign 1.
- **Worksheet:** extend `docs/notes/localised-jellium-theory.md`, don't restart.
- **Owns** the long-range cutoff (thread D).
- **Placement:** `docs/campaigns/localised_jellium/`.

Locked 2026-06-27 (campaigns grill, session 2):
- **Spine = "critical path first, one campaign":** sequence (1) cutoff (D) →
  (2) run-set 2 elongated-z SIE (B+C) → hand B+D to Campaign 1; THEN (3) run-set 1
  (F) → (4) charged-plates (A) + zero-point (E) → (6) synthesis.
- **E_SIE extraction (working design, reflected in the worksheet; formal frontmatter
  lock still pending):** BOTH routes cross-checked — route 1 far-launch/vacuum
  plateau + route 2 Hartree-matched classical subtraction; their agreement is the
  falsifiable test. All t=0 static energetics (single-point energy evals).
- **B+C+D worksheet inputs produced** for the user's external worksheet agent:
  `docs/notes/localised-jellium-sie-worksheet/` — `reference_pack.md` (material),
  `worksheet_plan.md` (md+XML spec), `resources/*.docx` (authored, 248 eqns) +
  two fetched arXiv PDFs (2307.03213, 1805.01377). See handover.

Deferred to the /campaigns pass: baseline geometry; the Lx/Ly/Lz/w grids; the
elongated-box dims + distance grid r; the cutoff prescription itself; what
"surface-effect metric" to tabulate; whether any new analysis kernel is needed
(pre-gate if so).
</resolved_decisions>

<guard_rails>
(deferred to the /campaigns pass)
- **Build simple → one parameter at a time** (vary a single config knob per run).
- **SIE extraction needs box-elongation convergence** — report E_SIE only once the
  residual is stable vs further elongation.
- Charge neutrality preserved at every config; N(t)≈const.
- Any new analysis kernel pre-gated (code-test + formula-validation + catalogue).
- Abort on NaN / complex energy / GPU occupied; 300-frame VTI cadence; results
  PROVISIONAL until cross-checked against the worksheet self-test.
</guard_rails>

<tasks>
(rough — done-criteria sharpened via /campaigns)
1. **Literature + worksheet extension** — pin the cutoff/finite-size + surface-jellium
   refs into `docs/sources/`; add charged-plates (A), SIE-decomposition (B), and
   cutoff (D) sections to the worksheet. *Done when:* sources + extended worksheet.
2. **Mine existing classical runs (H)** — scout VTI/observables for behaviour clues.
   *Done when:* a short findings note exists.
3. **Run-set 1 (F)** — config-parameter sweep → surface-effect table/notebook.
   *Done when:* parameter→surface-effect mapping recorded.
4. **Run-set 2 (B+C)** — WP-distance elongated-box sweep → E_SIE + WP–jellium
   repulsion decomposed. *Done when:* the Campaign-1 energy reference is pinned.
5. **Long-range cutoff (D)** — decide + validate the projectile cutoff prescription.
   *Done when:* prescription documented + handed to Campaign 1.
6. **Synthesis** — charged-plates mental-model writeup (A) + zero-point significance
   verdict (E) + worksheet self-test. *Done when:* the understanding doc exists.
</tasks>

<rules>
- ALWAYS keep `inq/` immutable; analysis/config work in `inqkit`/`inqview`.
- Ground every physical/mental-model claim per literature-review; label inferences
  explicitly ("Inference: …").
- Report energies at 2–3 s.f. per the rounding rule; carry Ha/eV units explicitly.
- New analysis kernel → code-test + formula-validation + catalogue row first.
</rules>

<preflight>
(rough draft — NOT yet autonomy-ready; reminder of what /campaigns must satisfy)
- [ ] Baseline geometry + Lx/Ly/Lz/w grids + elongated-box dims + distance grid set.
- [ ] SIE box-elongation-convergence criterion defined (numeric).
- [ ] Cutoff prescription chosen + validated; surface-effect metric defined.
- [ ] Worksheet extension + docs/sources notes planned; any new kernel pre-gated.
- [ ] Notebook/output contract + handover pointer; cross-link to Campaign 1.
- [ ] Grounding: cutoff/finite-size + surface-jellium refs cited; inferences labelled.
</preflight>
