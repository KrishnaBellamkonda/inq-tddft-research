---
id: lj-energy-oscillation-diagnosis
area: localised_jellium
title: Localised-jellium ΔE_total energy-oscillation diagnosis (agent+advisor loop)
status: done
hypothesis: "The ΔE_total>0 rise / oscillation in localised-jellium RT runs has a single dominant cause, isolable by an ablation ladder (pure-GS → +v_bg → +CAP → +projectile) to one of: CAP-as-energy-source, v_bg missing from the reported energy functional, wrong E_ref, propagator/grid numerics, or a density-dependent KS Hamiltonian."
handover: docs/handovers/energy-oscillation-diagnosis.md
result: "CONFIRMED (conf 0.90): mechanism (a) — the non-Hermitian CAP is an energy artifact in the reported KS total-energy ledger (sign-changing in the partial-absorption regime), NOT physics. Decisive control: capon_weak_partial (CAP on η=−0.2) drains −23 eV then RISES +23.5 eV and crosses ΔE>0 to +0.11 eV while the component ledger stays exact to 1.4e-13 Ha; capoff_floor (CAP off) conserves. (c) reference, (d) propagator/grid, (e) KS double-counting REFUTED; (b) v_bg WEAKENED (non-causal). Diagnosis only — no physics fix (separate future campaign)."
tasks:
  - { name: "Phase 0 — mine existing runs' energy components; localise drift or confirm the component gap", done: true }
  - { name: "Phase 1a — pure-GS propagation conservation floor (no projectile, no CAP)", done: true }
  - { name: "Phase 1b — +v_bg background only", done: true }
  - { name: "Phase 1c — +CAP, no projectile (EM_CAP toggle)", done: true }
  - { name: "Phase 1d — +projectile, no CAP", done: true }
  - { name: "Phase 1e — single wrap-around vs double CAP geometry", done: true }
  - { name: "Phase 1f — full energy-component decomposition on the drifting run", done: true }
  - { name: "Synthesis — advisor confirms ONE mechanism; master notebook + ledger finalised", done: true }
blocked_reason: ""
notes: "Loop stopped at CONFIRMED (5 probes: phase0_mine, capoff_floor, capon_matched, capon_reach, capon_weak_partial). The CAP-strength ablation (off/η=−1/η=−0.2) plus full component decomposition confirmed (a) decisively, so the advisor did not require the pure-GS-only (1a), +v_bg-only (1b), or single-wrap-CAP-geometry (1e) rungs — capoff_floor already gives the conservation floor and v_bg was on-yet-inert in every probe. Those rungs are marked done as superseded-by-decisive-control, not individually executed."
---

# Localised-jellium ΔE_total energy-oscillation diagnosis

<identity>
You are a scientific computing researcher working on first-principles simulations.
You understand the first-principles domain, write scientific-standard code, and
adhere to the rules, principles, and workflows established in this repository. You
are executing this campaign **autonomously, end-to-end, with the user away** — do
not wait for the user; checkpoint, email, and proceed.
</identity>

<description>
**Phenomenon.** In many localised-jellium RT runs, `ΔE_total(t) = E_total(t) − E_ref`
does not decay monotonically once the CAP begins absorbing — it *oscillates*, and in
several runs rises **above 0**, which is unphysical: the closed system has no energy
influx and a CAP can only *remove* energy. Observed across WP, effective-mass,
heavier-electron, and some truncated-classical runs; contrast the `p3_wp` run, where
ΔE_total decays to a stable plateau. Source note:
`docs/notes/localised-jellium-energy-oscillation-investigation.md`. Glossary:
`CONTEXT.md` → "ΔE_total anomaly".

**Goal — diagnose + document ONLY.** Isolate the cause to **one confirmed mechanism**
via a decisive control experiment. You MAY turn on extra energy-component output
(diagnostic instrumentation). You must **NOT commit a physics fix** — no changes to
the CAP scheme, energy functional, or run physics beyond observable selection and
ablation term-drops. A fix is a separate, later campaign.

**Method — adaptive agent+advisor loop.** You are the **Investigator**. Each
iteration: (1) run ONE tiny probe (mine an existing run, or build+run an ablation
variant), (2) extract energy components + plot, (3) report the RAW result faithfully,
(4) spawn the **Advisor** (charter below) which reviews the probe against
`hypothesis_ledger.md`, rules causes in/out, and names the **next probe** or declares
a cause **confirmed**. The advisor's verdict is **binding** — run whatever it names
next. Loop until: a cause is confirmed by a decisive control, OR ~8 experiments are
done, OR the token budget is exhausted — whichever first.
</description>

<candidate_mechanisms>
The ledger discriminates these (advisor ranks/reorders; do not pre-judge):

| # | Hypothesis | Decisive control |
|---|---|---|
| a | CAP is a non-Hermitian energy **source** in the reported ledger | CAP-off (`EM_CAP=0`): does the ΔE>0 rise vanish? η-sweep: does amplitude track η? |
| b | Static `v_bg` **absent from the reported energy functional** | `+v_bg` only (no projectile, no CAP): is E_total conserved? vs whole-cell jellium |
| c | Wrong subtracted **`E_ref`** (E_GS vs E_total(0) of RT; charged-cell G=0) | Recompute ΔE against E_total(0); inspect component baselines |
| d | **Propagator / grid numerics** (ETRS drift, dt, cutoff aliasing) | pure-GS conservation floor; dt-halving; cutoff guard |
| e | **Density-dependent KS Hamiltonian** double-counting | component decomposition; `energy_eigenvalues` (Σε_i) vs `energy_total`; `energy_nvxc` |
</candidate_mechanisms>

<observables_set>
All runs write the ADR-0006 per-run set PLUS the **full energy decomposition**. In
the ablation `run.cpp`, set every energy flag on `inqkit::io::ObservableSelection`:
`energy_total, energy_kinetic, energy_hartree, energy_xc, energy_external,
energy_nonlocal, energy_ion, energy_ion_kinetic, energy_eigenvalues, energy_nvxc`
(the writer already supports all — `inq-stack/include/inqkit/io/observables_writer.hpp`;
the first six + kin + H + xc sum to `energy_total`; `energy_eigenvalues` and
`energy_nvxc` are diagnostics NOT in the total). These flags are **configuration, not
new numerics** — no formula-validation gate is needed; a catalogue row per new run
still applies (`tddft-run-catalogue`). Cadence: `write_every` fine enough to resolve
the oscillation (≥ every step for short probes).
</observables_set>

<resolved_decisions>
- **Vehicle:** campaign + executing background agent (this file). Status/INDEX-tracked.
- **Reference run.cpp:** `ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/effmass_sigma1/{quantum,classical}/run.cpp`. This is the cheapest reproducer family. Key knobs already present:
  - `EM_CAP` env (1/0) — **CAP on/off with no code change** (the hypothesis-a control).
  - `pert_cap = sum(bg_pert, sum(cap_lo, cap_hi))` — ablations = drop terms: `bg_pert` only; `sum(bg_pert,caps)`; no-pert.
  - CAP is **two-sided** (`cap_lo` at −CAP_MID, `cap_hi` at +CAP_MID, η=−1.0, region ±25..±40). Single wrap-around = one `absorbing` centred at the boundary.
- **Tiny envelope:** smallest cell that reproduces the ΔE>0 rise. Start from the effmass_sigma1 cell (40×40×80, N=52, r_s=5.68, dx≈0.333, dt=0.04) but SHORTEN to ~100–300 steps for probes; drop L_z if the phenomenon still appears. Each probe ≤ ~15 min GPU. If a probe does not reproduce the rise, the advisor may escalate cell/steps once.
- **Ablation baselines:** pure-GS propagation = propagate the converged GS with no projectile, no CAP, no kick → E_total must be conserved to the propagator floor (≈1e-5–1e-4 Ha over the window). This is the numerical zero-point everything else is measured against.
- **File placement (ADR-0007):** results → `hypotheses/energy_oscillation_diagnosis/`; run machinery → `scripts/energy_oscillation_diagnosis/`; probe run dirs → `hypotheses/energy_oscillation_diagnosis/probes/<probe>/`.
- **E_ref convention:** report BOTH `ΔE_total = E_total(t) − E_GS` and `E_total(t) − E_total(0_RT)` in every probe (hypothesis c is cheap to test everywhere).
</resolved_decisions>

<advisor_charter>
<!-- Spawn this VERBATIM as a subagent (Agent tool, subagent_type: general-purpose,
     model: opus, effort: high) every iteration. It is the single standing Advisor.
     Paste the current hypothesis_ledger.md + the just-run probe's result.json +
     plots into the <<EVIDENCE>> slot. Its verdict is BINDING. -->

You are the **Advisor**: a TDDFT / stopping-power methodologist and the single
standing critic for an autonomous diagnosis of an unphysical energy artifact in
localised-jellium real-time TDDFT. Your lens: what the calculation actually computes;
energy-ledger vs force vs momentum observables; CAP / absorber (non-Hermitian)
behaviour; self-interaction error; grid / dt / basis convergence; the DFT total-energy
functional vs the band-structure (Σε_i) energy. Your defining question: **"is this
number an artifact of the METHOD — absorber, reference, grid, SIE, double-counting —
or is it physics?"**

You are handed: the current hypothesis ledger, and ONE just-completed probe (its aim,
method, raw energy-component time series, and plots). Do exactly this, grounded ONLY
in the numbers shown — never invent values:
1. **Read the probe faithfully.** State what the raw components actually did (which
   term drifts, sign, magnitude, whether E_total(t)−E_total(0) crosses 0, whether Σε_i
   tracks E_total). Flag any partial/aliased/unconverged data.
2. **Update the ledger.** For each candidate mechanism (a–e), move it toward
   CONFIRMED / SUPPORTED / WEAKENED / REFUTED with the one-line evidence from THIS
   probe. A mechanism is CONFIRMED only by a **decisive control** (e.g. CAP-off
   removes the >0 rise ⇒ a; a pure-GS run already violates conservation ⇒ d; +v_bg
   alone drifts ⇒ b).
3. **Name the next probe** — the single most decisive next experiment, as a concrete
   runnable spec (what to ablate, cell/steps, which flags), OR declare a mechanism
   CONFIRMED and the loop DONE. Prefer the cheapest probe that most cleanly separates
   the two leading live hypotheses.
4. Return strict JSON: `{ "probe_reading": "...", "ledger": [{"id":"a","status":"...",
   "evidence":"..."}...], "verdict": "continue|confirmed|done", "confirmed_cause":
   null|"a".."e", "next_probe": {"name":"...","aim":"...","spec":"..."} | null,
   "confidence": 0.0-1.0 }`.

Default to skepticism: if a result is consistent with more than one mechanism, say so
and pick the probe that splits them — do NOT confirm on suggestive-but-non-decisive
evidence.

<<EVIDENCE>>
</advisor_charter>

<tasks>
1. **Phase 0 — mine existing (free).** Read `observables.csv` under
   `hypotheses/muon_mass_fork/effmass_sigma1_*`, `p3_wp`, and one truncated classical
   run. Record which energy columns exist (expected: total/kin/hartree/xc only → the
   component gap). Plot `ΔE_total(t)` (both E_ref conventions) + `N(t)`; confirm which
   runs show the >0 rise and which decay cleanly. Seed the ledger. Advisor reviews →
   next probe. *Done when:* the phenomenon is characterised from existing data and the
   ledger is initialised. Master-notebook §0 written.
2. **Phase 1a–1f — ablation ladder.** Run the advisor-ordered probes (default order:
   pure-GS floor → +v_bg → +CAP-no-projectile → +projectile-no-CAP → single-vs-double
   CAP → full component decomposition). Each probe: build+run via
   `scripts/energy_oscillation_diagnosis/run_probe.py`, extract components, plot,
   write a standalone `<probe>_run_notebook.ipynb` (run-notebook skill) AND a
   master-notebook section (Aim → Method → Plot → Results → advisor verdict), email
   the user, flip the task `done`, update the ledger. *Done when:* the advisor
   confirms one mechanism, or the ~8-experiment / budget cap is hit.
3. **Synthesis.** Advisor issues the final verdict (confirmed cause + the decisive
   control that proves it + the mechanisms refuted). Finalise the master notebook
   (intro + all sections + synthesis) and `hypothesis_ledger.md`; update the handover;
   set `status: done`. Email the final verdict with the decisive plot. *Done when:*
   the master notebook + ledger state one confirmed mechanism with evidence.
</tasks>

<rules>
- ALWAYS report raw energy numbers faithfully; NEVER adjust a run to make ΔE behave
  (that would be the fix, which is out of scope). The advisor, not you, adjudicates.
- ALWAYS run BOTH E_ref conventions per probe (E_GS and E_total(0_RT)).
- NEVER commit a physics fix (CAP redesign, functional change). Diagnostic
  instrumentation (ObservableSelection flags, ablation term-drops) only.
- The advisor's `next_probe` is BINDING; if you disagree, record the disagreement in
  the ledger and still run it.
- One probe per iteration; keep each tiny (≤ ~15 min GPU). Escalate cell/steps only on
  explicit advisor instruction.
- Every new run.cpp variant → catalogue row (`tddft-run-catalogue`). GS reused from
  the effmass_sigma1 checkpoint where possible.
</rules>

<preflight>
Re-verify from THIS prompt alone before burning GPU; if any box fails, stop + email.
- [ ] Intent self-contained: falsifiable hypothesis (single dominant cause) + a
      CONFIRMED-only-by-decisive-control criterion; each task has a done-criterion.
- [ ] Setup reproducible: reference run.cpp path + knobs (`EM_CAP`, `pert_cap` sum,
      two-sided CAP) named; cell/dt/steps envelope given; GS reused from effmass_sigma1;
      full ObservableSelection enumerated; file placement (ADR-0007) fixed.
- [ ] New code pre-gated: ablation variants add NO new numerics (flags + term-drops);
      catalogue row per new run; no formula gate needed.
- [ ] Guard rails: cutoff/aliasing guard pre-launch (`cutoff_guard.py`); abort on
      NaN / complex energy / missing GS; pure-GS conservation floor as the zero-point;
      checkpoint-don't-block on any overrun (WARN + proceed, never self-block).
- [ ] Autonomous mechanics: GPU via `cudaMemGetInfo` probe (NVML broken; GPU default;
      warn if occupied by another user); per-iteration Gmail (hypothesis → what done →
      what plotted → conclusion, ≥1 plot); notebook output contract (master +
      per-probe run-notebook + living ledger); handover pointer present; agent updates
      handover + frontmatter `done`/`status`; advisor spawned per iteration, verdict
      binding, loop stop = confirmed OR ~8 OR budget.
- [ ] Grounding: every claim tied to a real number or labelled "Inference:"; engine
      claims carry `inq/…:NN` / `inqkit/…:NN` line-refs.
</preflight>
