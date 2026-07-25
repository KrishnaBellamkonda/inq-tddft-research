---
# ROUGH DRAFT — created 2026-07-06 during a grill-with-docs session as a sibling
# idea to the muon per-orbital-mass fork. NOT autonomy-ready, NOT to be run now.
# Its first task is literature research + a worksheet; then simple INQ GS demos.
# Harden later via /campaigns (gated authoring + INDEX registration).
id: mass-tuned-bands
area: mass_tuned_bands
title: "Effective-mass tuning of jellium band structure — emulating materials via m*"
status: draft
hypothesis: "Because the free-particle dispersion is E(k) = ħ²k²/2m, rescaling the
  per-particle mass in the Kohn–Sham kinetic operator (-ħ²/2m ∇²) rescales the
  band curvature and density of states of a jellium / simple periodic system. A
  tunable per-state effective mass in inq-study can therefore emulate different
  materials' band dispersions within one engine — light mass → wide, highly
  dispersive bands; heavy mass → narrow, flat bands (heavy-fermion-like). This is
  a MODEL that imposes m* by hand, not a first-principles derivation of m* from a
  periodic potential; jellium is the clean testbed because its bands are exactly
  parabolic so m sets the curvature directly."
handover: docs/handovers/mass-tuned-bands.md
tasks:
  - { name: "Literature research — collect + download all relevant papers (effective-mass approximation, k·p theory, band-structure engineering, heavy-fermion/flat-band systems, any DFT/TDDFT work using mass rescaling); build a reading worksheet + docs/sources notes", done: false }
  - { name: "Simple INQ GS demonstrations — GS of jellium / a simple periodic system at several global masses; extract eigenvalue spectrum E(k)/DOS; show the parabola scales as 1/m and matches ħ²k²/2m analytically", done: false }
blocked_reason: ""
---

# Effective-mass tuning of jellium band structure

<identity>
You are a scientific computing researcher working on first-principles simulations.
σ always means the wavepacket width σ_WP.
</identity>

<rough_draft_banner>
ROUGH DRAFT. Two tasks only, as specified by the user (2026-07-06): (1) research
+ worksheet, (2) simple INQ GS demonstrations. Do NOT execute as-is; harden via
/campaigns first (done-criteria, k-point/box/functional locks, validation gates).
</rough_draft_banner>

<motivation>
The per-state inverse-mass fork designed for the muon campaign
(`docs/campaigns/muon_projectile/inq_study_engine_notes.md` §2) is GENERAL: it can
set an arbitrary mass on any orbital (or all orbitals). Setting a *global* mass on
a jellium bath turns the engine into a band-curvature knob:

- Free-particle / jellium dispersion: **E(k) = ħ²k²/2m**. Halving m doubles the
  band width and curvature; multiplying m by 200 flattens the band ~200×.
- This is exactly the **effective-mass approximation** of solid-state physics
  (Ashcroft & Mermin; Kittel) run in reverse: instead of *reading* m* off a
  computed band curvature, we *impose* m to emulate a target curvature/DOS.
- **Materials modelling angle:** a single tunable-mass jellium can stand in for
  systems with very different band masses (light-mass semiconductors vs
  heavy-fermion / flat-band correlated systems) at the level of the free-particle
  dispersion and DOS.
</motivation>

<engine_dependency>
Uses the SAME `inverse_mass_` fork as the muon campaign (§2 of the muon notes) —
NO new engine machinery. Two differences from the muon-projectile use:
1. **Mass is GLOBAL** (all orbitals share m), not per-orbital.
2. **Band structure is a GROUND-STATE property**, so the mass MUST be applied in
   the GS path (`ground_state/calculator.hpp` + `ground_state/initial_guess.hpp`),
   not only in `real_time/propagate.hpp`. The generalisable design threads
   `electrons.inverse_mass()` through ALL ks_hamiltonian construction sites,
   defaulting to all-ones (backward compatible). See muon notes §2 + §5.
</engine_dependency>

<task_1_literature>
**Research + worksheet.** Collect and DOWNLOAD the relevant literature; turn it
into a reading worksheet.
- Topics: effective-mass approximation & its validity; k·p perturbation theory;
  band-structure engineering / strain-tuned masses; heavy-fermion & flat-band
  physics; any DFT/TDDFT work that rescales particle mass or uses fictitious
  masses (e.g. Car–Parrinello fictitious electron mass — DISTINGUISH: that is a
  dynamics trick, not a physical band mass); mass-dependent HEG / jellium scaling
  (r_s in effective atomic units — cross-ref muon notes §3 point 3).
- Deliverables: `docs/sources/<author-year>.md` notes per key paper
  (literature-review skill) + a single worksheet (a notebook or md) summarising
  the technique, its assumptions, and where imposing-m is/ isn't legitimate.
- Ground every claim (scientific-grounding rule); label inferences.
</task_1_literature>

<task_2_gs_demo>
**Simple INQ GS demonstrations that the mechanism works.**
- Build GS of jellium (and/or a simple periodic test system) at several global
  masses m ∈ {e.g. 0.5, 1, 2, 10, 207}·m_e using the inverse-mass fork.
- Extract the eigenvalue spectrum and, with **k-point sampling** (NOTE: this needs
  a multi-k-point run — the existing gamma-only jellium WP machinery does NOT
  resolve E(k); a band path / k-grid is required), plot E(k) and the DOS.
- **Validation:** show E(k) scales as 1/m and matches the analytic ħ²k²/2m for the
  free/jellium case (formula-validation + code-test). Confirm the m=1 case
  reproduces the unforked engine bit-for-bit.
- Done when: a figure of E(k) vs m + a table of fitted curvatures vs 1/m, with the
  analytic overlay, in an executed notebook.
</task_2_gs_demo>

<open_questions>
- Which "simple periodic system" (bare jellium is parabolic-trivial; a weak
  periodic potential would show real band folding/gaps whose curvature THEN
  depends on both the potential and m). Decide the demonstrator system.
- k-point strategy (band path vs uniform grid) and functional (LDA for jellium).
- Does imposing global m and using electron-parameterised LDA-XC give a physically
  meaningful DOS, or only a kinetic-parabola demo? (cross-ref muon notes §3.)
- Relationship to the muon campaign's GS-mass requirement — shared engine work.
</open_questions>

<rules>
- NEVER edit `inq/`; the mass fork lives in `inq-study` only.
- Reuse the muon `inverse_mass_` fork; do not build a parallel mechanism.
- Ground the effective-mass framing in literature; label the "impose-m" modelling
  choice explicitly (it is not first-principles m*).
- Report numbers at 2 s.f. (3 s.f. for near-equalities).
</rules>
