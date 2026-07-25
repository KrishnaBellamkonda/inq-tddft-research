---
id: ml-pattern-finding-wp-classical
area: ml-patterns
title: "Pattern-finding in WP/classical runs — quantum signatures in induced density"
status: ready
hypothesis: "In bulk jellium, the bath induced-density wake of a finite-σ quantum wavepacket projectile differs from that of a point classical projectile at matched velocity in two interpretable, SIE-controlled ways the scalar S(v) cannot capture — a form-factor softening of the q-space induced-density ratio n_WP(q)/n_classical(q) that follows F_WP/F_ONCV (≈exp(−q²σ_pot²/2)) across the σ-sweep, and a collective wake whose DMD frequency obeys λ(v)=2πv/ω_p across the velocity-sweep — both recovered after the GS→rigid-motion→Lindhard→vacuum-WP subtraction ladder, with verdicts read on a held-out cell split."
handover: docs/handovers/ml-pattern-finding-wp-classical.md
tasks:
  - { name: "T0 — run database built + validated (data inventory)", done: true }
  - { name: "T1 — pre-gate new kernels (induced-density normaliser/subtraction ladder, POD, DMD): code-test + formula-validation + catalogue row", done: false }
  - { name: "T2 — Rung 1 (bulk jellium) form-factor cut: POD on Δn, σ-sweep @E=100 vs exp(−q²σ_pot²/2); done = held-out agreement within ±20% reported (confirm/refute/inconclusive)", done: false }
  - { name: "T3 — Rung 1 wake gate: DMD wake cut @σ=5, λ(v)=2πv/ω_p; done = held-out agreement within ±20% reported", done: false }
  - { name: "T4 — Rung 1b (localised slab): transfer frozen pipeline, wake @σ=0.5; geometry comparison vs bulk", done: false }
  - { name: "T5 — Rung 2 (dynamics): DMD/Koopman + SINDy on bulk-jellium induced density", done: false }
  - { name: "T6 — Exploratory: exchange/xc-hole + diffraction on vacuum-WP-subtracted field (caveated, not headline)", done: false }
  - { name: "T7 — Synthesis notebook + handover + frontmatter updates", done: false }
  - { name: "T8 — [REDO] bulk-only scope + cell pinning: freeze form-factor cut, classical (coulombic) + WP velocity sweeps, pin Track-A + Track-B calib/held-out splits", done: true }
  - { name: "T9 — [REDO] new weak-form SINDy/PDE-FIND kernel, PRE-GATED (formula-validation CONFIRM + 6/6 code-tests recovering known PDEs + catalogue rows)", done: true }
  - { name: "T10 — [REDO] Track A gates re-run clean on bulk-only: form-factor + wake held-out verdicts", done: false }
  - { name: "T11 — [REDO] Track B: discover PDE_classical (coulombic point sweep); three walls. VERDICT: REFUTE (0/3). Panel: real 'no resolvable wake PDE' but cause confounded (nonlinear screening vs transverse-mean dilution — needs on-axis line cut)", done: true }
  - { name: "T12 — [REDO] Track B: discover PDE_WP (sigma=5). Pipeline VERDICT: CONFIRM (u_tt=c^2 u_xx, c≈v). *** PANEL OVERTURN (2026-07-04): ARTIFACT — velocity-sweep used WP-INCLUDED density (dN=+1e blob translating at v); bath-only _wf gives u_tt=0. c≈v is blob kinematics, NOT collective physics. See handover panel milestone. ***", done: true }
  - { name: "T13 — [REDO] compare: 'ratio 5' is a cross-velocity artifact (WP E100 / classical E25). Panel-invalidated alongside T12.", done: true }
  - { name: "T14 — [REDO] synthesis + notebook + handover + frontmatter flips", done: true }
blocked_reason: ""
---

# Pattern-finding in WP/classical runs — quantum signatures in induced density

> Authored interactively (grill-with-docs + campaigns Mode A), 2026-06-30/07-01.
> Every substantive decision is locked in `<resolved_decisions>`. T0 (DB) is built +
> validated (2 rounds, PASS), the `<preflight>` checklist is green, and the user
> reviewed + approved the prompt — `status: ready` set 2026-07-01. A fresh agent may
> execute T1–T7 autonomously; re-verify `<preflight>` before any heavy work.

<identity>
You are a computational condensed-matter physicist with ML/AI expertise, working
on first-principles rt-TDDFT simulations. You understand jellium stopping-power
physics, write scientific-standard code, and adhere to this repository's always-on
rules (σ_WP convention, never-fftshift-a-VTI, inq/ immutable, validation gates,
number rounding, scientific grounding).
</identity>

<description>
**Problem.** The project's core target is the electronic stopping power of a light
projectile in jellium. Prior work reduced the classical-vs-quantum comparison to a
single scalar — the "quantum component of stopping" S_WP − S_classical — which is
contaminated by the WP self-interaction error (SIE, ~7 eV) and dispersion, and
throws away all spatial structure.

**This campaign** takes the novel lens: study the **induced electron-density field**
n(r,t) (the VTI series already on disk) instead of the scalar S, and use
**interpretable ML** to find and *explain* the spatial/dynamical differences
between a **point classical** projectile and a **finite-σ quantum wavepacket** at
matched velocity. It is an **analysis-only** campaign — no new INQ runs.

**Decision it informs.** Whether the induced-density field carries a measurable,
defensible quantum signature beyond linear response and beyond the SIE artefact —
i.e. whether "quantum stopping" has a spatial fingerprint the scalar S(v) misses.
This feeds the thesis (a real physics result) or bounds the quantum component
spatially (a publishable null).

**Gap (grounded).** ML scalar-stopping surrogates are mature (Ward et al., npj
Comput. Mater. 2024 — density → scalar force); induced-density analysis to explain
stopping is done by hand (Kononov et al. 2025). Nobody applies interpretable ML to
the induced-density *field* to characterise the classical↔quantum difference. See
`docs/campaigns/ml-patterns/research/ml_induced_density_research.md` + `docs/sources/`.

**Success / failure (falsifiable).** CONFIRM if, on **held-out** cells, the
form-factor q-ratio n_WP(q)/n_classical(q) matches F_WP/F_ONCV (≈exp(−q²σ_pot²/2)
where F_ONCV≈1) and the DMD wake frequency matches λ(v)=2πv/ω_p, both within ±20%.
REFUTE if a method-valid analysis misses both.
INCONCLUSIVE if method validity is unreachable in ≤4 tries. A refute/inconclusive
is a legitimate, reported outcome — never retried into a confirm.
</description>

<observables_set>
**Consumed (existing primary observables, from the run database):** the density
VTI *series* per run — `density_total`, `density_wp`, `density_system` — plus the
GS density (`density_gs_system`) for the induced reference. Cadence per run =
`frame_dt_au = dt_au × write_every` (in the DB). The campaign's working field is
the **bath response** n_bath(t) − n_bath^GS, n_bath = n_total − n_wp (WP) /
n_total (classical) — the canonical run-independent bath density.

**New derived kernels (PRE-GATED in T1 before any headline analysis):**
- induced-density normaliser / subtraction ladder (reuses `lindhard_elf`,
  `fourier`, `pipeline/wake`, `bath_energy`);
- **POD** (SVD truncation) and **DMD/Koopman** (windowed) — formula-bearing →
  each gets a `formula-validation` agent + `code-test` + a catalogue row.
Persistent-homology descriptors and SINDy are pre-gated likewise before T5/T6.
All new code is **campaign-local** under `docs/campaigns/ml-patterns/` (see
ADR — not promoted to `inqview` unless later proven useful).
</observables_set>

<resolved_decisions>
All locked in the 2026-06-30 grill; each ended in an explicit user lock.

- **Aim** — physics DISCOVERY via an interpretable ML lens (the physics result is
  the product; ML is the representation/discovery tool).
- **Primary target** — the *defensible core*: (iii) collective **wake**
  phase/wavelength + (iv) **form-factor softening**. Discovery axis = the nonlinear
  residual vs Lindhard linear response. Exchange/xc-hole (i) + diffraction (ii) are
  **exploratory only** (T6), on the vacuum-WP-subtracted field with explicit
  artifact caveats (research: SIE-confounded).
- **Substrate** — **point-classical (ONCV, `classical_potential_form=coulombic`)
  vs finite-σ WP at matched velocity** (point→smeared wake *is* exp(−q²σ_pot²/2)).
  **No new runs.** DB confirms ~362 directed `point_vs_wp` matches in jellium.
- **System scope (complexity ladder)** — Rung 1 = **bulk jellium** (consistent
  grid, abundant); Rung 1b = **localised slab** (thesis geometry, σ_WP=0.5,
  `sigma_matched_gauss`); Rung 2 = dynamics on bulk. Form-factor *requires* bulk
  (only bulk has the σ-sweep).
- **Rung-1 grid (two orthogonal cuts):** form-factor cut = fixed **E=100 eV**,
  σ_WP ∈ {0.5,1,3,5,8}; wake cut = fixed **σ_WP=5**, E ∈ {20…1500} eV. Both have
  matched point-classical partners with hundreds–thousands of density frames.
- **Observable** — bath response n_bath(t) − n_bath^GS (above).
- **Methods** — Rung 1 form-factor (T2): the **headline metric is the q-space
  induced-density ratio** R(q) = n_WP(q)/n_classical(q) (azimuthally + temporally
  reduced to a robust scalar per q per σ), compared to the analytic form-factor
  prediction. **POD** on Δn (= n_WP − n_classical, matched v) is the *real-space
  structural* support (mode geometry of the difference), NOT the form-factor metric
  itself; persistent-homology descriptors support it. Rung 2: **DMD/Koopman** (modes
  tagged with frequency+decay), **windowed over the early near-constant-velocity
  stretch** (light-projectile deceleration rule), with `dt < π/ω_p`; SINDy in latent
  space.
- **Falsification** — **parameter-free** overlays (no fitting of σ_eff or ω):
  - **Form-factor (T2):** R(q) = n_WP(q)/n_classical(q) vs the KNOWN
    **F_WP(q)/F_ONCV(q)**, where F_WP(q)=exp(−q²σ_pot²/2) (σ_pot from the WP UPF) and
    **F_ONCV(q) is the ACTUAL classical-projectile form factor** (the ONCV electron
    is NOT a true δ — F_ONCV is computed from its UPF charge density in T1, never
    assumed = 1; the prediction reduces to exp(−q²σ_pot²/2) ONLY where F_ONCV≈1 over
    the resolved q-range, which T1 must establish).
  - **Wake (T3):** DMD dominant ω vs the KNOWN ω_p=√(4πn); λ(v) vs 2πv/ω_p.
  - **±20% agreement band**; every verdict read on **held-out** cells (pinned split
    below).
- **Pinned calibration/held-out split (pre-registered; ADR 0011)** —
  - *Form-factor cut* (E=100 eV): **CALIBRATION** σ_WP ∈ {1, 5}; **HELD-OUT**
    σ_WP ∈ {0.5, 3, 8}.
  - *Wake cut* (σ_WP=5): sort the matched energies ascending by velocity;
    **CALIBRATION** = even-index energies (0-indexed 0,2,4,…), **HELD-OUT** =
    odd-index energies. (Deterministic; no hardcoded energy list.)
  - The ≤4-try loop tunes the **shared pipeline config** to maximise ±20% agreement
    on the **CALIBRATION** cells of both cuts; freezes it; reports the T2 verdict on
    form-factor HELD-OUT and the T3 verdict on wake HELD-OUT. All ≤4 attempts logged.
- **Autonomy** — **fully autonomous Python orchestrator** (NOT bash), idempotent
  resume, per-phase Gmail (4-part + ≥1 plot). After T1, the Rung-1 analysis runs a
  **bounded agentic validation+retry loop (≤4 tries)** that maximises ±20%
  agreement **on the pinned CALIBRATION split (above)**, **freezes** the winning
  pipeline config, and reports the verdict from the **pinned HELD-OUT cells**; **all
  ≤4 attempts are logged**. After ≤4 it proceeds with the best config to T4–T7.
  (See ADR `docs/adr/0011-held-out-split-anti-phacking.md`.)
- **Code placement** — everything **campaign-local** under
  `docs/campaigns/ml-patterns/` (kernels + glue + tests); `inqview` untouched
  (reversible; promote later if proven). Validation gates still apply in place.
- **Notebooks** — one auto-built study notebook **per rung** + a **final synthesis**
  notebook, under `docs/campaigns/ml-patterns/notebooks/`; each records the exact
  DB cells, the calibration/held-out split, the frozen config, all ≤4 attempts, and
  the verdict. Per-phase email carries that rung's key figures.
- **Input** — the reproducibility-grade **run database** (`docs/run_database.csv` +
  `.json`, 581 runs × 137 cols), built + independently validated (T0). The campaign
  selects cells from it; it is the single source of run truth.
</resolved_decisions>

<guard_rails>
- **Subtraction ladder is mandatory before ANY ML:** GS → rigid-projectile-motion
  → Lindhard linear-response → vacuum-WP/SIE. Raw n(r,t) ML is forbidden (it
  "discovers" the GS, the trivial translation, and the SIE).
- **Never `np.fft.fftshift` a VTI** (physical order; load via `inqview.load_vti`).
  Only LEED `.dat` are FFT-natural (not used here).
- **Anti-p-hacking (ADR 0011):** the retry loop tunes config on the **calibration**
  split only; the verdict is **always** read from **held-out** cells; the split is
  **PINNED** (form-factor: calib σ_WP∈{1,5}, held-out σ_WP∈{0.5,3,8}; wake:
  even/odd-velocity-index energies) — not chosen at runtime; **all ≤4 attempts are
  logged** in the notebook. Never tune on, or report agreement from, held-out cells.
- **Honesty:** CONFIRM / REFUTE / INCONCLUSIVE are all valid outcomes. A refute is
  reported, not retried away. INCONCLUSIVE = method validity unreachable in ≤4.
- **Abort/skip conditions:** a cell whose density series is missing/short, NaN/Inf
  in a field, or grid-shape mismatch within a cut → skip + log (never silently
  drop; report what was dropped). DMD requires `dt < π/ω_p` and a coherent dominant
  mode, else flag method-invalid.
- **Grounding:** every scientific/numerical choice cited or labelled "Inference:".
  The research's `[abstract-only]` method citations (DMD-on-electron-phonon,
  POD/DMD & persistent-homology reviews, contrastive-physics, β-VAE) must be
  verified in T1 before they enter a notebook/manuscript. The SIE≈7 eV figure is
  from the project brief, not externally verified — label as such.
- **Number rounding:** 2 s.f. default in all reported numbers / captions / tables.
</guard_rails>

<tasks>
- **T0 — run database (data inventory).** Built (581×137) + independently validated
  (graphene σ fix, classical_potential_form by V(r), reworked twins + match_type,
  filled σ/velocity). *Done-criterion:* round-2 validation PASS. Skills: (built via
  the campaign's `build_run_database.py` + `validate_run_database.py`).
- **T1 — pre-gate new kernels.** Build the induced-density normaliser, POD, DMD
  (campaign-local). *Done:* each formula-bearing kernel passes `formula-validation`
  + `code-test` (known-case) + has a catalogue row, BEFORE T2. Verify the
  `[abstract-only]` method citations (`literature-review`). **Also: compute the
  ACTUAL ONCV projectile form factor F_ONCV(q) from its UPF charge density and
  establish the q-range where F_ONCV≈1** — so the T2 prediction (F_WP/F_ONCV) is
  honest, never an assumed point charge.
- **T2 — Rung 1 form-factor cut.** Split is PINNED (calib σ_WP∈{1,5}, held-out
  σ_WP∈{0.5,3,8} @E=100). Headline metric = the **q-ratio R(q)=n_WP(q)/n_classical(q)**
  vs the KNOWN **F_WP(q)/F_ONCV(q)** (F_WP=exp(−q²σ_pot²/2)); POD on Δn gives the
  real-space structural support. ≤4-try loop maximises ±20% agreement on the
  calibration σ; freeze; *Done:* **held-out** verdict (confirm/refute/inconclusive)
  on σ_WP∈{0.5,3,8} + all attempts logged in `rung1_bulk_…ipynb`.
- **T3 — Rung 1 wake gate.** DMD (windowed) on the σ_WP=5 velocity-sweep; split is
  PINNED (calib = even-velocity-index energies, held-out = odd). ±20% to
  λ(v)=2πv/ω_p and ω vs ω_p on the **held-out** energies. *Done:* held-out verdict +
  attempts logged.
- **T4 — Rung 1b localised slab.** Transfer the FROZEN pipeline to the σ_WP=0.5
  slab (wake only); compare geometry vs bulk. *Done:* slab verdict + bulk-vs-slab
  panel in `rung1b_slab.ipynb`.
- **T5 — Rung 2 dynamics.** DMD/Koopman + SINDy on bulk-jellium induced density;
  mode spectrum vs Bohm-Gross. *Done:* `rung2_dynamics.ipynb` with mode table.
- **T6 — Exploratory exchange/diffraction.** On the vacuum-WP-subtracted field
  only; *Done:* `exploratory_…ipynb` with explicit SIE-artifact caveats; NOT a
  headline claim.
- **T7 — Synthesis.** Cross-rung `synthesis.ipynb`; update handover + flip
  frontmatter `done`/`status`.
Each Tn: orchestrator runs it, emails the 4-part result (+plot), flips `done`.
</tasks>

<rules>
- ALWAYS pre-register the calibration/held-out split before analysing a cut; read
  every verdict from held-out cells; log all ≤4 attempts.
- NEVER optimise for, or report, agreement on the held-out cells (p-hacking).
- NEVER claim a confirmed signature without the full subtraction ladder applied.
- NEVER `np.fft.fftshift` a VTI; ALWAYS load via `inqview.load_vti`.
- A refute or inconclusive is a valid result — report it; do not retry it into a
  confirm.
- All new kernels pre-gated (formula-validation + code-test + catalogue) before use.
- This campaign launches NO INQ runs and needs NO GPU scheduling.
</rules>

<preflight>
Analysis-only campaign — adapted from the standard checklist. Re-verify from this
prompt alone before running:
- [ ] **Intent self-contained:** falsifiable hypothesis + ±20% held-out
  confirm/refute/inconclusive criteria; every Tn has an unambiguous done-criterion.
- [ ] **Inputs reproducible, zero guessing:** the run database (581×137,
  validated) is the source of run truth; the (σ,E) cells per cut are named; the
  bath-response observable + cadence are defined; file placement is campaign-local.
- [ ] **New code pre-gated:** POD, DMD, normaliser → formula-validation + code-test
  + catalogue row BEFORE T2; `[abstract-only]` citations verified.
- [ ] **Validation & guard rails:** mandatory subtraction ladder; pre-registered
  split + held-out verdict + logged attempts (anti-p-hacking ADR); never-fftshift;
  skip+log on missing/NaN/shape-mismatch; DMD `dt<π/ω_p` + coherent mode; ≤4-try
  cap then proceed with best.
- [ ] **Autonomous mechanics:** Python orchestrator (not bash), idempotent resume,
  per-phase Gmail (4-part + plot), per-rung auto-built notebooks + synthesis,
  handover pointer present; agent updates handover + frontmatter done/status. NO
  GPU/dispatcher/pilot-gate (analysis-only).
- [ ] **Grounding:** every scientific/numerical choice cited or "Inference:";
  parameter-free predictions (σ_pot, ω_p) sourced; SIE figure labelled unverified.
</preflight>

---

# 2026-07-03 — Bulk-only PDE-discovery redo (two-track, autonomous, 12 h)

> Authored via grill-with-docs, 2026-07-03. Extends this campaign in place (user
> lock). Reuses the validated T0 run database + T1 kernels; adds Track-B
> governing-equation discovery as the headline. Every decision below ended in an
> explicit user lock during the grill. **Not yet executed — awaiting user go.**

## Why redo
The prior autonomous run (T1–T7, 2026-07-01) landed **T2 form-factor CONFIRM**
(2/3 held-out, weak calibration), **T3 wake INCONCLUSIVE** (1/3), **T4 slab
INCONCLUSIVE**, and a **thin, unvalidated T5 SINDy** (a 2-mode latent ODE with
fitted coefficients, no forward-prediction test, no physics interpretation). The
part the user most wants — *governing differential equations suggestive of
physics* — was the weakest. This redo narrows scope and makes equation discovery
rigorous.

## Goal string (user)
> "identify the patterns in the induced density in bulk jellium caused due to
> classical and quantum projectiles. Find similarities and differences. Try to
> find differential equations that might [be] suggestive (of some physics)."

## Locked decisions (grill 2026-07-03)
- **Scope** — **pure bulk jellium only**; DROP the localised slab (old T4) and all
  other systems. Bulk has the σ-sweep + a clean velocity sweep for both projectiles.
- **Backbone** — **two-track**, judged together at the end:
  - **Track A** (falsifiable spine) — re-run the form-factor + wake gates *clean on
    bulk-only*; render held-out verdicts. Reuses T2/T3 logic + frozen kernels.
  - **Track B** (headline, discovery) — discover a **governing field PDE** for the
    induced bath density.
- **DE form** — **field PDE headline** (`∂ₜ²n = L[·]` via **weak-form
  SINDy / PDE-FIND** on `n_bath(r,t) − n_bath^GS`) **+ latent ODE support** (SINDy
  on POD modes as a reduced-order-model cross-check).
- **Library** — **broad, agnostic, minimal priors**; sparsity decides. Physical
  names (`−ω_p²n`, `∇²n`, `v·∇n`, projectile source) assigned **post-hoc** by a
  physics-interpreter, never seeded. (ADR 0012.)
- **Separate-then-compare** — discover `PDE_classical` and `PDE_WP` on their **own**
  runs (NOT the difference field); compare terms afterwards. Matching only needed
  at the comparison step, so discovery is data-robust.
- **Three validation walls** (ADR 0012, extends 0011) — a term is "physics" only if
  it survives: (1) pinned calibration/held-out **cell split** (tune on calib,
  report from held-out); (2) **temporal forward-prediction** (fit early window,
  integrate forward, score later window); (3) **bootstrap coefficient stability**.
- **Engine** — a **Python orchestrator** extending
  `docs/campaigns/ml-patterns/orchestrate.py`: idempotent/resumable, per-phase
  Gmail (4-part + ≥1 plot), spawning **discovery + adversarial-skeptic +
  physics-interpreter + synthesis-judge** subagents. **Hard 12 h wall-clock cap.**
- **Termination (aim-reached gate)** — stop when BOTH deliverables validate
  (Track-A held-out verdicts rendered; Track-B PDEs pass all three walls for
  classical AND WP + comparison done), OR the validation metric **plateaus**
  (K rounds no gain), OR **12 h** elapses — then report best-so-far honestly. A
  refute / inconclusive / partial is a valid, reported outcome.
- **New code pre-gated** — the weak-SINDy/PDE-FIND kernel is formula-bearing →
  `formula-validation` + `code-test` (recover a KNOWN PDE from synthetic data:
  heat/advection) + catalogue row BEFORE any headline use (T9). POD/DMD/normaliser/
  form-factor kernels are already validated (T1); reused unchanged.
- **Code placement** — campaign-local under `docs/campaigns/ml-patterns/`;
  `inqview` untouched. Notebooks per track + a synthesis notebook.

## Pinned bulk-jellium cells (from the validated run DB, verified 2026-07-03)
All at r_s≈5.69, L50, grid **125³** (dx≈0.4), ω_p≈3.5 eV (plasma period ≈49 au;
frame_dt 0.02–0.06 au is Nyquist-safe).
- **Form-factor cut** (E=100 eV, σ-sweep): σ_WP ∈ {0.5, 1, 3, 5, 8}, bath-only
  `density_system`, 230–330 frames. *Track-A split (ADR 0011): calib σ∈{1,5},
  held-out σ∈{0.5,3,8}.*
- **Classical velocity sweep** (coulombic **point**): E ∈ {20,25,50,100,300,600} eV,
  190–457 bath frames. (Classical *gaussian* runs are unusable — 6 coarse frames.)
- **WP velocity sweep** (σ_WP=5 primary, σ_WP=1 secondary): E ∈ {20,25,50,100,300,
  600} eV, matching the classical sweep. σ=3 @E=25 has ~10 001 frames → prime cell
  for the temporal forward-prediction wall.
- **Wake split** (Track A, σ=5): calib = even-velocity-index E, held-out = odd
  (ADR 0011). Skip+log E=700/1500 (1–4 frames) and E=15 (frame_dt=4.0).
- **Track-B split** (pinned): calibration E ∈ {20,50,300}, held-out E ∈ {25,100,600}
  (deterministic; identical for classical and WP so the comparison is on the same
  held-out velocities).

## Redo tasks (T8–T14)
T8 scope+cell pin → T9 pre-gate PDE-FIND kernel → T10 Track-A gates (bulk) →
T11 discover PDE_classical → T12 discover PDE_WP → T13 compare + latent-ODE +
interpret → T14 synthesis judge + notebooks + handover/frontmatter. The
orchestrator wraps T10–T13 in the aim-reached gate; T11/T12 each log all refine
attempts and read the verdict from held-out only.

## Guard rails (in addition to the originals above)
- Broad agnostic library is defensible ONLY behind the three walls — never report a
  term from calibration or from a pointwise `∂ₜn` fit that fails forward-integration.
- "No stable interpretable PDE survives" is a valid, reported outcome — do NOT tune
  the library until something appears.
- Subtraction ladder (GS → rigid-motion → Lindhard → vacuum-WP) still mandatory
  before ANY discovery; never `np.fft.fftshift` a VTI; load via `inqview.load_vti`.
- Number rounding 2 s.f.; every physical claim cited or "Inference:".
