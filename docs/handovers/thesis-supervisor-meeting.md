# Handover — Emilio thesis-supervisor meeting deck (26 Jun 2026)

Rolling handover for building the supervisor presentation. Workspace:
`/local/data/public/skcb2/tddft/docs/reports/26-06-2026-meeting-emilio/`.

## 2026-06-25 — draft-1 build underway (Sections 1–2 done)

### Goal
A fresh assertion-evidence talk: "A localised jellium slab makes the stopping
power of a quantum wavepacket and a charge-matched classical projectile directly
measurable, revealing the difference between them." (thesis-claim FIRM, user
2026-06-25 — sign/magnitude/cause of the difference left open.)

### Source documents (all in the workspace dir)
- `00_simulation_inventory.md` — 6-thread inventory.
- `initial_idea.md` — user's idea dump (the spec).
- `02_storyline_skeleton.md` — XML storyline, 22 message-beats, thesis-claim firm.
- `03_nuts_and_bolts_plan.md` — 25-slide build spec (idea/plots/equations per slide).
- `plots_log.md` — per-figure production log (appended by each build agent).
- `figures/` — output PNGs/GIFs (600 DPI); `figures/equations/` — equation PNGs.
- `build/build_sectionN.py` — per-section build scripts.

### Build loop (the agreed pattern)
For each section: one general-purpose **build agent** (precise spec, venv python,
canonical `inqview.visualisation.style` theme, σ_WP labelling, .png only, no image
preview) → one **haiku validator agent** (checks code/data/method/conventions, NOT
visuals — user owns visual review). Continue until the deck is built.

### DONE + validated
- **Section 1** (Slides 2–4): `fig_s1_potential.png`, `fig_s1_system_design.png`,
  `fig_s1_stopping_sv.png`, `fig_s1_stopping_se.png`. Method A (fixed-20%-time,
  free-intercept slope, se error bars) on classical r_s=5.69 runs; Lindhard overlay.
  Validator PASS 8/8.
- **Section 2** (Slides 5–9): `fig_s2_workflow.png`, `equations/eq_mask.png`,
  `equations/eq_cap.png`, `fig_s2_cap_overlay.png`, `fig_s2_cap_absorb.gif`,
  `fig_s2_mask_reflect.gif`, `fig_s2_reflectivity_L.png` (+`_linear`),
  `fig_s2_eta_sweep.png`, `fig_s2_eta_L_grid.png`. Validator PASS 10/10.
- **Section 3a** (Slides 10–13, S3.0/S3.1/S3.1b/S3.2): `fig_s3_injection_workflow.png`,
  `fig_s3_gs_baseline.png`, `fig_s3_sie_table.png`, `fig_s3_energy_bookkeeping.png`,
  `equations/eq_sie.png`, `eq_stopping.png`, `eq_energy_decomp.png`. Validator PASS 10/10.
  - **EXACT SIE def found** (`qsp_phase1/sie/run.cpp:9-10`): `E_SIE = E_total(0) − E_GS −
    KE_WP` = 0.1618 Ha = 4.40 eV. KE_WP zero-point ≈ 80.8 eV (theory 3/(4σ²)=81.6),
    dominates and exceeds SIE — annotated. No longer "TODO/invent".
  - **GS baseline (N=82 box):** N=82 (∫n), r_s=5.667, n₀=1.312e-3, k_F=0.339, L_z=25.
    GS VTI `qsp_phase1/gs/results/density_gs_system/density_gs_system.vti`.
  - **Twin run = qsp_phase2** (`p2_wp`/`p2_classical`, N=82, 50×50×70, σ_WP=0.5, E=100eV)
    — this is the genuine 50×50×70 box (the user's "most recent test run"), NOT
    `fullsuite_*` (N=234, 50³). qsp_phase2 also = the stopping-power skill's worked
    examples. **Use qsp_phase2 for Section 3b** unless user says otherwise.

### Verified facts / corrections (must survive)
- **kF = 0.337 Bohr⁻¹** for r_s=5.69 (n=1.296e-3), NOT 0.473. (An Explore agent gave
  0.473; build+validator both re-derived 0.337 = 1.9192/r_s. Load-bearing for
  Lindhard.) Used correctly in Section 1.
- **ω_p = 3.473 eV** at n₀=1.296e-3 (lowest absorbable energy; jellium-application
  bound — absorber benchmark itself was vacuum/free-WP. Cross-system caveat printed.)
- **Reflectivity ε** = surviving inner-region WP fraction ∫|ψ|²/N₀ (from
  `twosided_cap_vs_mask/build_twosided_report.py`). Chosen L=20 Bohr.
- **CAP ε values are PROVISIONAL** until the inq-study engine regression (Task #7).
- **Fourier gate LIFTED** (user 2026-06-25): build Slides 17/18 with the
  `fourier-analysis` skill; keep the Δω≈9 eV resolution flag.
- **slab≈bulk citation validated**: Quijada et al., PRA 75, 042902 (2007); source note
  `docs/sources/quijada-2007-cluster-bulk-stopping.md`. Velocity-gated condition; the
  preceding mechanism sentence is user-paraphrase, not verbatim (flagged).
- σ-WP convention is now an always-on rule: `.claude/rules/sigma-wp-convention.md`.

### NOT done yet (remaining build)
- **Section 3** (Slides 10–23, S3.0–S3.12 + S3.1b): injection workflow, GS baselines,
  SIE table, energy-bookkeeping system design (Quijada cite), WP-vs-classical 2×2
  density gifs, energy curves, 2D momentum, **loss function + FFT spectra (Fourier
  skill)**, spreading problem, broadening plot, muon outlook (FUTURE — workflow only,
  no run), wide-electron + quantum-S(v)-sweep placeholders.
  - Data threads: `localised_jellium` (`qsp_phase1`, `fullsuite_wp`,
    `fullsuite_classical`), `jellium_wp_stopping`. SIE=4.40 eV, E_GS=−45.759 Ha.
  - Suggested sub-batches: (3a) S3.0–S3.5,S3.8,S3.9 data figs; (3b) S3.6–S3.7 Fourier;
    (3c) S3.10–S3.12 outlook/placeholders.
- **Conclusion** (Slide 24) + **Appendix** (A1 cap_in_jellium, A2 quantum_classical_nocap
  → real dir `jellium/hypotheses/qvc_nocap_sigma3`).
- **Equation PNG consolidation**: global list in `03_nuts_and_bolts_plan.md` (13 eqns).
  Still TODO: reflectivity exact form, SIE exact form (do NOT invent — confirm).
- **Deck assembly** (pptx) once figures approved — reuse prior `build_*_deck.py` /
  pptxgenjs from earlier `*-meeting-emilio` folders.

### Open decisions (user)
- Visual review of PNGs (user is the only previewer).
- KE_WP in the energy decomposition (Slide 13) must be the MEASURED WP KE (~80 eV
  zero-point dominates, exceeds SIE) — confirm before transcribing to a slide.

## 2026-06-25 (later) — figure-standard overhaul + Section 1 REMADE

### New skill + conventions (user feedback)
- **`scientific-figures` skill** created (`.claude/skills/scientific-figures/SKILL.md`
  + bundled `workflow_render_template.py`); cross-referenced from `report-figures`
  and `notebook-making`; terms added to `CONTEXT.md`. It is the figure RULE SET:
  - No analytical text / arrowheads on canvas (→ caption/slide).
  - Legends name series identity ONLY.
  - **System-design plots = REAL density** (total-density xz via `load_vti`, dashed
    slab+CAP extents) — never cartoons. Params → slide, not figure.
  - **Colorbar outside the axes, same height as panel** (`make_axes_locatable`).
  - Tables: header row coloured only.
  - **Workflows = `.drawio` source + matched matplotlib PNG** (no CLI; pattern =
    `docs/diagrams/build_contribution_page.py` + `render_contribution_png.py`, idiom:
    Carlito/Calibri, muted palette, containers/lanes, orthogonal edges, 1920×1080).
  - Presentation → titles PRESENT (report → cropped). **Captions live in the slide
    spec** (`03_nuts_and_bolts_plan.md` `- **caption:**` lines), never on canvas.
  - "linear response" is the label for the Lindhard/RPA curve (never "Lindhard").

### σ CONVENTION √2 FINDING (critical — affects all σ labels)
- The classical **σ-sweep run-dir names are charge-std σ_pot**, NOT σ_WP. Verified:
  `electron_gaussian_sigma0p5.upf` V(0)=1.596 Ha ⇒ σ_charge=0.5; config says
  "sigma=0.5 Bohr Gaussian **charge**". So **true σ_WP = √2 × dir number**.
- User decision: **label by TRUE σ_WP**. Slide-4 curves are now
  σ_WP = 0.35 (run sigma0p25), 0.50 (sigma0p35), 0.71 (sigma0p5), 4.24 (sigma3p0).
- **Working width = σ_WP = 0.5** (the WAVEPACKET runs, e.g. qsp_phase2, are genuinely
  σ_WP=0.5 — zero-point 3/(4·0.5²)=3 Ha=81.6 eV confirms). Its **matched classical**
  is σ_charge=0.354 = the `sigma0p35` run = the **σ_WP=0.50 curve**. Do NOT call
  σ_WP=0.71 the working width.
- WP runs are named by σ_WP; CLASSICAL runs are named by σ_charge. Always reconcile
  via σ_WP = √2 σ_charge before labelling.

### Section 1 — REMADE + validated under new standard (PASS 11/11)
- `fig_s1_potential.png`: annotations + arrowheads stripped; legend = identities only.
- `fig_s1_system_design.png`: REAL total-density xz slice from
  `run_classical_n162_L50_E100_v2/results/raw/vti/density_total/` via `load_vti`
  (no fftshift); dashed box edges + projectile trajectory; NO fictitious CAP (bulk
  S(v) runs have no engine CAP, only a one-traversal path cap); colorbar outside,
  same height; no param legend.
- `fig_s1_stopping_sv.png` / `_se.png`: multi-σ convergence, **log-x**, "linear
  response" reference, true σ_WP = 0.35/0.50/0.71/4.24 (Gaussian ladders only,
  ONCV `_v2` runs discarded), Method-A error bars; σ=0.71 series psp-gated to
  `sigma0p5.upf` (5 pts; `sig0p4_v1p0` excluded).
- Captions added to Slides 2–4 in the plan.

### PENDING remakes (Sections 2 & 3a were built pre-standard — must be redone)
- **Section 2:** workflow `fig_s2_workflow` → `.drawio`+matplotlib (was cartoon);
  `fig_s2_eta_sweep` is the WRONG plot → use the correct reflectivity curve from
  `twosided_cap_vs_mask_study.ipynb`; `fig_s2_eta_L_grid` colorbar → outside +
  same height; minimal legends; "linear response" n/a here.
- **Section 3a:** `fig_s3_injection_workflow` → `.drawio`+matplotlib; `fig_s3_gs_baseline`
  → SPLIT into two separate plots (no number pile-up, drop legend); SIE table →
  header-row-only colour; `fig_s3_energy_bookkeeping` system-design → REAL density
  (Δn or total) of qsp_phase2 with dashed slab+CAP extents.
- Then the still-unbuilt: Section 3b (twin comparison, qsp_phase2), 3c (Fourier via
  `fourier-analysis` skill), 3d (spreading/broadening/muon/placeholders),
  conclusion, appendix, equation consolidation, pptx assembly.

## 2026-06-25 (later still) — DRAFT 1 COMPLETE

**Deliverable:** `docs/reports/26-06-2026-meeting-emilio/emilio_meeting_draft1.pptx`
(28 slides, 16:9; builder `build/build_deck.py`; python-pptx). Captions in speaker
notes (AE: clean canvas). All figures built; sections 3b/3c/3d added this round.

### Built this round
- **3b (stills, not gifs — gif render timed out):** `fig_s3_wake_compare.png` (2×2
  total+Δn, WP|classical, frame t=0.198 fs), `fig_s3_norm_vs_time.png`,
  `fig_s3_energy_curves.png`, `fig_s3_momentum_2d.png`. Builder
  `build/build_section3b_stills.py`.
- **3c (Fourier skill):** `fig_s3_loss_function.png` from the RESOLVED long run
  `run_plasmon_n162_L50_E15` (T=2000 au, Δω=0.086 eV; low-q peaks on ω_p≈3.47 eV,
  locator only); `fig_s3_response_spectra.png` from p2_wp (T=40 au, Δω=4.27 eV —
  UNRESOLVED, bold-flagged). p2_wp has no n_q(t), hence E15 for the loss function.
- **3d + 13 equation PNGs:** `fig_s3_spreading.png` (measured in-medium σ(t) +
  analytic free-spread overlay), `fig_s3_broadening.png` (analytic), 
  `fig_s3_muon_spread.png` (analytic outlook). Equations in `figures/equations/`.

### PHYSICS FLAGS surfaced in the headline twin (qsp_phase2) — for user, NOT style
1. **Norm loss NOT comparable:** ΔN_WP=0.865 vs ΔN_classical=0.037 (~24×). This
   CONTRADICTS the idea-dump's "WP norm-loss ≈ classical norm-loss coincidence"
   (Slide 14 callout). Flagged to-verify in caption/log.
2. **ΔE_WP = −117.3 eV (NEGATIVE)** vs ΔE_classical = +70.6 eV. WP twin shows net
   energy DECREASE over the run — threatens the S=ΔE/L_z bookkeeping for the WP.
   Surface before any slide claim rests on it.

### Draft-1 pending (known, for the feedback pass)
- **Figure-style problems remain** (user: fix later). Esp.: Section 2 & 3a NOT yet
  remade to the `scientific-figures` standard; workflow diagrams (`fig_s2_workflow`,
  `fig_s3_injection_workflow`) are still matplotlib CARTOONS → convert to
  `.drawio`+matplotlib idiom; `fig_s2_eta_sweep` is the wrong plot; `eta_L_grid`
  colorbar; `gs_baseline` split; SIE table header-only colour.
- **Appendix A1/A2 figures pending** (placeholders in deck): export cap_in_jellium
  B0–B3 and quantum_classical_nocap (`jellium/hypotheses/qvc_nocap_sigma3`).
- **Slides 22/23** intentional planned-result placeholders (wide-electron; quantum
  S(v) sweep — runs not done).
- The two physics flags above need user resolution before the bookkeeping slides
  are trusted.

---

## 2026-06-26 — Figure-correction pass (grill-with-docs, FB-001…FB-018)

Full review pass on `emilio_meeting_draft1.pptx`. Running good/bad record in
`docs/reports/26-06-2026-meeting-emilio/feedback_log.md` (FB-001…FB-018, all
APPLIED). Deck rebuilt → **32 slides** (added 3 section dividers + 1 absorptive-
loss slide).

**Done + regenerated (all build scripts re-run clean):**
- Systemic: tight-bbox save (`style.save_presentation`, FB-001); `×10⁻ⁿ` notation
  via `use_mathtext` + `style.sci_notation` (FB-003).
- §1 system_design: dashed start-plane line, no trajectory, snapshot time `t=4.30 a.u.` (FB-002).
- §2: cap_overlay per-band sin² SHAPE fixed vs `absorbing.hpp` (peak at band
  centre, was wall) + de-cluttered (FB-005); eta_sweep → two-sided L=20 (FB-007);
  eta_L_grid square + matched colorbar (FB-006); reflectivity de-texted (FB-008);
  gif y-labels `⟂`→`n_WP`, gifs regenerated (FB-004).
- §3a: gs_baseline split into two panels, metrics → slide (FB-014).
- §3b: energy_curves → `ΔE(t)` + annotated refs (FB-013/015); momentum →
  `|k⊥| vs |k_z|` map from FFT of `wavefunction_wp` (FB-017).
- §3c: loss function log + linear + **absorptive `|Im n_q|/q²`** + plasmon/e-h
  overlays (FB-016); −Im[1/ε] normalisation FLAGGED uncertain (Fourier gate).
- §3d: broadening legend outside (FB-011); muon recast as energy×σ broadening,
  comparable to electron (FB-018).
- Workflows: Graphviz `.dot`→PNG, role-coloured nodes, arrows kept
  (`build/build_workflows.py`, FB-010).
- Deck: section dividers + energy-bookkeeping **table** (FB-009/012); captions in
  black speaker notes.

**Skill refined:** `.claude/skills/scientific-figures/SKILL.md` — added the save
rule, arrowhead scope, Graphviz §7, and §10 (×10ⁿ, no exotic glyphs, verify-shape
-vs-source, signal-vs-baseline, provenance) + checklist.

**Still open (not deck-blocking):**
- Two physics flags persist and need user resolution: `ΔE_WP = −117.3 eV`
  (negative) vs `ΔE_cl = +70.6 eV`; WP norm drop 0.865 vs classical 0.037 (~24×,
  contradicts the "comparable norm-loss" claim on slide 14).
- −Im[1/ε] absolute normalisation (FB-016) — needs the external-perturbation
  factor; only lineshape/peak claimed.
- Doc bug: `build_twosided_report.py:64` LaTeX has the same `/(2·L_half)` CAP
  shape error as FB-005 — fix the study doc later.
- FB-010 node-styling note ("…exact neutrality and same object…") never parsed;
  implemented as role-coloured nodes pending user restatement.

---

## 2026-06-26 — DRAFT 2 built (additive plots + gifs + section names)

`emilio_meeting_draft2.pptx` — **38 slides** (draft 1 = 32; +6 new §3 slides),
10.3 MB, zero missing-image warnings. Draft 1 untouched. Plan:
`docs/plans/emilio-draft2.md`. Builders: `build/build_section3_new.py` (plots),
`build/build_section3b_gifs.py` (gifs), `build/build_deck2.py` (deck).

**Task 1 — 4 new §3 plots:** `fig_s3_efield.png` (|E| xz, mid frame step 400),
`fig_s3_kl_divergence.png` (KL→~20), `fig_s3_ks_energies.png` (ΔE_i bar, max
193 eV), `fig_s3_momentum_kt.png` (|k|-vs-time carpet).

**Task 2 — gifs:** BUILT 4 strided (every-4th, 51-frame) xz comparison gifs
`fig_s3_{total,delta}_{wp,cl}.gif` (2×2 slide 141; still kept) + animated
`fig_s3_norm_anim.gif` (dN_wp=0.865 vs dN_cl=0.037). REUSED `fig_s3_ks_delta_wp/cl.gif`
+ `fig_s3_momentum_anim.gif`. CAP/mask 1D gifs kept (not swapped).

**Task 3 — dividers retitled** to exactly: "Gaussian Potentials", "Absorbing
boundary conditions", "Localised Jellium slab".

**New slides:** 141 (2×2 animated wakes), 142 (norm anim), 161 (momentum
|k|-vs-time carpet+gif), 181 (KS ΔE bar+gif), 182 (E-field), 183 (KL).

**User will now edit draft 2 manually.** Carried-over flags unchanged (ΔE_WP
negative; WP norm-loss ~24×; −Im[1/ε] normalisation).
