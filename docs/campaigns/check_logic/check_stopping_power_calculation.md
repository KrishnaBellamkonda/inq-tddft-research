---
id: check-stopping-power
area: check_logic
title: "Stopping-power + Fourier-analysis training"
status: done
hypothesis: "A hand-built, fully-understood stopping-power + Fourier-analysis workflow reproduces the outsourced pipeline on known-answer cases and can be locked into deterministic, self-contained skills."
handover: docs/plans/stopping-power-fourier-training.md
tasks:
  - { name: "Stopping-power notebook: sections 2-6 + critique + rebuild; hand-S reproduces stopping_vs_v on v3p0; stress-tested v0p8/v0p6; method locked in README", done: false }
  - { name: "Fourier-audit notebook: synthetic -> QKE v0p0626 energy -> E15 n_q; every fourier.py method validated or flagged", done: true }
  - { name: "Fill verdicts on 3 validation dossiers (loss-function-formula / fft-drift-removal / fft-normalization) + fix density_fourier BUG A/B via code-test + catalogue row", done: true }
  - { name: "Encode stopping-power-extraction skill under .claude/skills/ + README", done: true }
  - { name: "Encode fourier-analysis skill under .claude/skills/ (the named gating deliverable)", done: true }
blocked_reason: "DONE 2026-06-25 by user (gate released): strict gate-release condition MET — both skills exist (stopping-power-extraction + fourier-analysis) AND all three dossier verdicts filled (loss-function/fft-drift-removal/fft-normalization). density_fourier BUG-A/B fixed via code-test. Residual non-gating follow-up (task 1): the stopping_power_extraction.ipynb energy-regression rewrite + S(v)-vs-Lindhard notebook + README — folded into the encoded skill; does NOT block the gate."
---

# Proper training: stopping power for classical projectiles + Fourier-transform analysis

<campaign_type>
**INTERACTIVE TRAINING CAMPAIGN — gate node. NOT an autonomous run-set.**

This campaign is deliberately human-in-the-loop: the user instructs, the agent
performs, the user critically evaluates and gives feedback, the agent improves —
repeated until the user is confident the method is correct. It therefore does
**not** follow, and is **exempt from**, the campaigns autonomy-readiness checklist
(there is no `<preflight>` / no fresh-agent autonomous execution). It launches no
GPU runs; it consumes the outputs of runs that already exist.

Its role in the campaign system is twofold:
1. **Progress tracker** — the 5 `tasks:` above are the checkpoints.
2. **Gate node** — `status: done` is the signal that RELEASES the gates on two
   other campaigns:
   - `quantum-kick-extension` (status `blocked` until the **fourier-analysis
     skill** here exists), and
   - `cap-jellium-loss-function` (tasks 5–7 already gated on this campaign).

   The specific named deliverable those campaigns wait on is the
   **`fourier-analysis` skill** (task 5). Gate-release condition, locked with the
   user 2026-06-22: **both skills exist AND the three validation-dossier verdicts
   are filled** (the strict "skill exists + verdicts filled" criterion).
</campaign_type>

<identity>
You are a research assistant in a lab who is currently being trained by the user
to properly do these tasks. Specifically: correctly calculating stopping-power
for classical projectiles, and using Fourier transforms properly on the energy /
density data with sanity checks that confirm it is being done correctly. You
write scientific-standard code and adhere to the repository's rules and workflows.
</identity>

<description>
For training, we take concrete examples and the user hand-holds the agent through
calculating the stopping power, then through Fourier-transform analysis. For both,
a tight feedback loop: the user instructs, the agent performs, the user evaluates
and gives feedback, the agent improves — until the user is confident it is right.
All work lives in `.ipynb` files so the user sees code and visualisations
together.

**Why this campaign exists.** The out-of-box defaults for Fourier transforms and
peak analysis in Fourier space, and the extraction of stopping-power metrics, were
largely outsourced. The user wants hands-on control and to understand each step
(what Hann windowing is, what detrending is, how S(v) is extracted), then to lock
each into a replicable, deterministic workflow encoded as a skill.

**The decision it informs.** Until this is complete the project may NOT build any
new loss-function `L(q,ω)` or Fourier/spectral analysis (the standing
loss-function gate). Completing it unblocks the QKE peak-frequency analysis and
the loss-function feasibility campaign.

**Success / failure (per task, see `<tasks>`).** Success = each method is
understood, reproduces the existing kernel on a known-answer case, is critically
stress-tested, and is locked into a self-contained skill; the three validation
dossiers carry the user's verdicts. Failure / stop = a method cannot be made to
reproduce its known-answer case, or a dossier verdict rejects the formula — in
which case surface it rather than papering over.
</description>

<context>
The user believes the Fourier-transform defaults and Fourier-space peak analysis
were largely outsourced, and wants to regain control and understand each step
(Hann windowing, detrending, etc.) to build a replicable workflow. Likewise the
user is not 100% sure how the stopping-power metrics are extracted and wants to be
hands-on. Work one problem at a time, build up by understanding every step, then
encode everything in deterministic workflows and skills.
</context>

<artefacts_on_disk>
Grounded inventory (2026-06-22) so a session resuming this knows what already exists:

- **Plan (live doc):** `docs/plans/stopping-power-fourier-training.md` — grilling
  complete; method shape, ladder, example all agreed.
- **Stopping notebook:** `docs/validation/stopping-power-extraction/stopping_power_extraction.ipynb`
  — **Section 1 only** (load & look). Sections 2–6 + critique + rebuild remain.
- **Fourier notebook:** `docs/validation/stopping-power-extraction/fourier_analysis.ipynb`
  — **not created**.
- **README:** `docs/validation/stopping-power-extraction/README.md` — **not created**.
- **Production code (audit targets, do NOT treat as black boxes):**
  - `inq-stack/python/inqview/analysis/stopping_extract.py` — `Track`, `load_track`,
    `stopping_vs_v(track, transient_bohr=3.0, window=11)`. Production-ready.
  - `inq-stack/python/inqview/analysis/fourier.py` — `FourierTransform`,
    `WindowSpec`, `FourierResult`; carries the embedded TODOs (lines ~8–13)
    doubting windowing/detrend/convenience methods, "especially for the
    QuantumKickExtension run".
  - `inq-stack/python/inqview/pipeline/density_fourier.py` — n_q(ω) loss-function
    path. **BUG A (~line 313):** takes `sig.real` before the time-FFT (folds ±ω,
    halves amplitude, loses direction). **BUG B (~line 344):** plots `|n_q(ω)|`,
    not the intended `|n_q(ω)|²/q²`.
- **Validation dossiers awaiting the user's verdict lines:**
  - `docs/validation/loss-function-formula-validation.md` (peak-locator vs true
    −Im[1/ε]; BUG A/B; proposes 3 synthetic numpy tests)
  - `docs/validation/fft-drift-removal-validation.md` (baseline subtraction:
    initial vs mean vs detrend)
  - `docs/validation/fft-normalization-validation.md` (window coherent-gain,
    interior-bin doubling)
</artefacts_on_disk>

<resolved_decisions>
Locked with the user during the 2026-06-22 grill:

- **campaign_type** = interactive training / gate node; exempt from autonomy
  checklist (see `<campaign_type>`).
- **gate_release** = strict: both skills exist AND all three dossier verdicts
  filled. `status: done` then releases `quantum-kick-extension` and
  `cap-jellium-loss-function`.
- **skill_count** = **two separate skills** — `fourier-analysis` (the named
  gating deliverable) and `stopping-power-extraction`. Each self-contained under
  `.claude/skills/<skill>/` (skills ship all their own artefacts; nothing in
  `docs/`).
- **method_shape** = understand-in-context from first principles → the existing
  kernel is the *destination*, not a black box → critically stress-test →
  rebuild a clean deterministic version. Independent re-derivation keeps judgement
  honest (mirrors the formula/test-validation independence ethos).
- **stopping example** = `run_sv_sigma0p5/results/v3p0/` (v₀=3.0, ~constant v,
  clean single slope). Stress-test held in reserve: `v0p8` (Barkas into noise),
  `v0p6` (sub-v_F, noisiest). r_s=5.69 jellium, m=m_e, free Ehrenfest, σ=0.5 erf
  electron; track columns `step,time_au,x,y,z,vx,vy,vz` (motion along z).
- **fourier example** = three-stage, known-answer first:
  1. **Synthetic** signals with analytic-truth peaks — validate window
     coherent-gain, baseline subtraction, and expose BUG A (real-before-FFT
     folding). Use the 3 synthetic tests proposed in
     `loss-function-formula-validation.md` (undamped plasmon; real-part folding;
     damped oscillator line-shape).
  2. **`fourier.py` on the QKE Li v0p0626 multi-k energy series** — known answer:
     plateau-detrend FFT → **6.480 eV vs paper 6.5 eV** (the run the TODO names).
  3. **`density_fourier.py` on the E15 jellium long run** n_q — known answer:
     ω_p ≈ 3.473 eV, Δω ≈ 0.09 eV.
- **output contract** = both notebooks executed against the `inqview-venv` kernel,
  outputs embedded; figures as `.png` with the canonical theme; the README ties
  the two notebooks and records the LOCKED conventions (which baseline mode for
  which observable, window/normalization convention, transient cutoff, the loss-
  function peak-locator caveat).
- **file_placement** = notebooks + README under
  `docs/validation/stopping-power-extraction/`; skills under `.claude/skills/`;
  any code fix to `density_fourier.py` is real library code → `code-test` +
  catalogue row in `docs/validation/test-catalogue.md`.
</resolved_decisions>

<tasks>
The agent flips the matching frontmatter `done` flag and updates the live plan as
each completes. Tasks 1–3 are interactive (user in the loop); 4–5 encode the
locked result.

1. **Stopping-power notebook.** Continue `stopping_power_extraction.ipynb` from
   Section 1 through the first-principles ladder:
   (2) define S = −dE_proj/dx = −dKE/ds, one global linear fit → single number
   (Ha/Bohr and eV/Å); (3) motivate + apply the transient cut (`transient_bohr`),
   re-fit; (4) cross-check S′ = +dE_electrons/ds on the same window, use E_total
   drift as integrator health; (5) why bin by v(t) (v≈const here → one slope;
   decelerating runs need a *local* slope — bridge to the kernel design);
   (6) call `inqview.analysis.stopping_extract.stopping_vs_v` and show it
   reproduces the hand number. Then **critique** (window size, transient length,
   finite-diff vs fit, KE-vs-E_elec disagreement, low-v breakdown on v0p8) →
   **rebuild** clean. *Done = notebook executed end-to-end; hand-derived S matches
   the kernel on v3p0; method locked in the README.* Uses `notebook-making`,
   `code-test`.

2. **Fourier-audit notebook.** Create `fourier_analysis.ipynb` following the
   three-stage example above. Audit every `fourier.py` method against first
   principles: window coherent-gain normalization, interior-bin doubling, the
   three baseline modes (initial / mean / detrend) and which is correct for which
   observable, transient cutoff, zero-padding (interpolation vs information),
   Gaussian frequency smoothing, and the convenience methods
   (`transform_energy/current/dipole`). *Done = each method validated or flagged
   with reasoning; synthetic stage reproduces analytic truth; QKE stage recovers
   6.48 eV; E15 stage recovers ω_p; notebook executed.* Uses `notebook-making`,
   `literature-review` (grounding window/baseline conventions), `code-test`.

3. **Verdicts + density_fourier bug fixes.** Fill the user's verdict lines on the
   three dossiers. Apply the fixes the audit confirms in `density_fourier.py`:
   **BUG A** (use the complex signal, not `.real`, before the time-FFT) and
   **BUG B** (compute/plot `|n_q(ω)|²/q²` per the documented intent, with the
   peak-locator caveat in bold). *Done = three verdicts filled; fixes applied;
   `code-test` known-case passes (the synthetic loss-function tests); a row added
   to `docs/validation/test-catalogue.md`.* Uses `code-test`, `formula-validation`.

4. **Encode the stopping-power-extraction skill.** Ship a self-contained
   `.claude/skills/stopping-power-extraction/` capturing the locked workflow
   (load track → S = −dKE/ds → transient cut → local slope vs v → kernel call →
   stress-test checklist). *Done = skill exists, self-contained, references the
   notebook as the worked example.* Uses `notebook-making` conventions.

5. **Encode the fourier-analysis skill (gating deliverable).** Ship a
   self-contained `.claude/skills/fourier-analysis/` capturing the locked Fourier
   methodology: window + coherent-gain, baseline-per-observable rule, transient
   cutoff, zero-pad/smoothing guidance, peak-attribution checklist, and the
   loss-function `|n_q|²/q²` peak-locator caveat. *Done = skill exists,
   self-contained → this campaign's `status` may flip to `done`, RELEASING the
   gates on `quantum-kick-extension` and `cap-jellium-loss-function`.* Uses
   `notebook-making` conventions.
</tasks>

<rules>
- **ALWAYS** keep the user in the loop on tasks 1–3 — instruct → perform →
  evaluate → improve. Do not "run ahead" and declare a method correct without the
  user's confirmation.
- **NEVER** build any new loss-function / Fourier spectral analysis for OTHER work
  while this campaign is incomplete — the standing loss-function gate
  (`feedback_fourier_loss_function_gate`) remains in force until task 5 ships.
- **ALWAYS** validate against a known-answer case before trusting a method
  (synthetic analytic truth, then the QKE 6.48 eV / E15 ω_p anchors).
- **NEVER** treat the existing kernels as black boxes — reach them by independent
  derivation (the whole point of this training).
- Skills are self-contained — bundle ALL artefacts skill-locally; nothing in
  `docs/`.
- `density_fourier.py` is real library code — any fix goes through `code-test` +
  a catalogue row, never "compiles ⇒ works".
</rules>

<gate_relationship>
Downstream consumers of this campaign's completion (do not edit them from here;
they read this campaign's `status`):
- `docs/campaigns/quantum_kick_extension/quantum_kick_extension.md` — `blocked`
  until the **fourier-analysis skill** (task 5) exists.
- `docs/campaigns/cap_in_jellium/loss_function_hypothesis_checking.md` — tasks 5–7
  already gated on `check-stopping-power`.
When this campaign flips to `status: done`, a session may unblock those two.
</gate_relationship>
