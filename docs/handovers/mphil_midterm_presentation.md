# Handover: MPhil mid-term presentation (assertive-evidence)

## Current status
Infrastructure done; in the **brainstorm phase** (user-led loop). Workspace, AE
skill, requirements, rubric all written. Brainstorm task outputs in
`brainstorms/tasks/`.

**2026-05-31 milestone:**
- New skill `tddft-run-catalogue` (scanner + CSV at **docs/runs_catalogue.csv**,
  80 runs), hooked into `tddft-simulations`. AE skill gained a phase-loop
  protocol, a `brainstorms/tasks/` rule, and a doc-conventions rule (define
  abbreviations; explain+trust-label new papers).
- `brainstorms/tasks/literature-key-messages.md`: **Nazarov–Gross & Yao–Schleife
  read in full and summarised** (key ideas + evidence); exact-factorisation +
  IMFP + abbreviations added. `jellium-density-motivation.md`: r_s=5.69≈Cs,
  r_s=3.41 actually ≈Li (not Na — flagged).
- **Coronene run `run_broadening_35x35x80` set up & launch-ready** (config
  `broadening_35x35x80.hpp`, run.cpp, analyse.py; GS gs_35x35x80_cut40). Launch
  z=+30, N_STEPS=848, WRITE_EVERY=1; centroid-gate PASS (ends −35, inside ±40).
  **BLOCKED: GPU driver/library mismatch (NVRM 535.288.01 vs lib 535.309.01) —
  needs nvidia module reload / reboot (root). Not launched.**
- Pending: T3 loss-function stopping (method explained, awaiting approval);
  T5 coronene pre-collision extraction (run blocked); exact-factorisation
  deep-PDF read optional.

**2026-05-31 (cont.) — WP broadening validated, 1D + 3D:**
- Run `run_broadening_35x35x80` completed; `density_rt_wp` VTI = WP orbital
  density alone (t=0 z-centroid = +30 confirms, no molecular density mixed in).
- **CANONICAL free-spread law (pin this — factor-2 error recurred twice):**
  `density` width, not wavefunction width. dens0 = sigma0/sqrt2;
  spreading time **tau = 2*dens0^2 = sigma0^2**;
  `sigma_dens(t) = dens0 * sqrt(1 + (t/tau)^2)`. With sigma0=1.0016 →
  dens0=0.708, tau=1.003, giving 5.57 Bohr at collision (t=7.82). WRONG version
  used tau=2*sigma0^2 (=2.006) → predicts ~2.85, half-rate. Cohen-Tannoudji free
  Gaussian, m=hbar=1.
- **Pre-collision validation PASS:** z-only mean|fwhm-analytic| 0.136 Bohr;
  **3D: sigma_x=sigma_y=sigma_z to ~0.13 Bohr each (perfect isotropy), r_rms
  residual 0.227 Bohr.** Transverse centroid frozen (cx,cy drift 0.031 Bohr).
  Slight +0.25-0.32 Bohr excess only at the last point (t=8.40, collision onset)
  — consistent with WP starting to feel coronene and/or transverse box (35 Bohr,
  half=17.5; 3-sigma~18 approaches boundary near collision).
- Assets: `docs/presentations/assets/precollision_broadening{,_3d}.{py,png,csv}`.
  Both scripts now use the correct law; PADded x-axis (T_COLL+0.75 a.u.).
- **Open recommendation:** a dedicated vacuum free-propagation run (large cubic
  box, no molecule, no orthogonalisation) would give a caveat-free isotropic-
  spreading panel AND double as the methods "running-race" GIF source
  (brainstorm L26). Coronene pre-collision is already sufficient to make the
  point; vacuum run is optional polish — awaiting user decision.

**2026-06-01 — GS-occupations (e-h decomposition) plot: saved-orbital blocker:**
- Goal: correct the mechanisms-slide excitation plot by filling the **grey
  "untracked" region** of `fig_gs_decomposition.png` (orbitals > 100) so the
  charge discrepancy → 0, using analytic jellium plane-wave GS orbitals
  101–150 overlapped with the final KS orbitals. Run = **E25**
  `run_wp_n162_L50_E25_sigma1_v2`; plasmon stays on E15 (the two panels are
  necessarily different runs — short snapshot vs 2000 a.u. for ω resolution).
- **Physics confirmed:** the +0.497 e "untracked" charge is **100 % bath
  electron–hole excitation** (the WP's initial occupation is f≈0 in
  `occupations.csv`, so the GS-projected occupation tracks only the 162 bath
  electrons). User's premise was correct; the grey region genuinely is e-h.
- **BLOCKER (code-level proof):** the corrected plot **cannot be produced from
  saved data — a re-run is required.** `inqkit::observables::OrbitalOverlapMatrix`
  (`inq-stack/include/inqkit/observables/orbital_overlap.hpp:60-61, 90-122`)
  computes `|⟨ψ_i^GS|ψ_j(t)⟩|²` in-loop from GPU memory and writes ONLY the
  squared scalars onto the truncated 100-orbital GS basis; it never `.save()`s
  any orbital. The only real-space orbital persisted for E25 is the **WP**
  (`results/raw/vti/wavefunction_wp/`). No bath orbitals, no INQ state-save, no
  restart checkpoint exist (searched exhaustively; GS save dir holds only
  `density_gs_system.vti` + observables). From squared, basis-truncated overlaps
  one cannot recover the complex evolved orbitals nor re-project them onto new
  analytic orbitals 101–150 — so the grey region is unrecoverable in post.
- The momentum-distribution `n_total − n_wp` tail recovers the e-h charge only
  to ~10 % (radial |k| binning smears the Fermi-shell boundary; raw per-G
  coefficients not saved), so it cannot close the discrepancy to *exactly* zero.
- **Fix (awaiting user go-ahead):** short E25 re-run with the overlap observer's
  `n_ref` raised 100 → ~200 (the GS orbitals ARE the analytic plane waves), so
  the bath e-h excitation above the old cutoff is captured per-orbital and the
  discrepancy closes to ~0. Edit `run_template.hpp:269`, ensure GS carries
  ≥200 states, rebuild, re-propagate (GPU ~30–60 min), rerun inqview `overlap`.
- **Lesson recorded in skill** `tddft-simulations` Phase 3d: any GS-occupations
  graph MUST persist the GS + final-state KS orbital wavefunctions AND the
  overlaps — not just the in-loop overlap scalars — or run with `n_ref` wide
  enough to cover the expected excitation.

**2026-06-02 — storyline3 deck edits + read-aloud script:**
- `drafts/build_storyline3_deck.py` edited + rebuilt (`drafts/storyline3.pptx`,
  21 slides). Changes:
  - **S12 mechanisms:** swapped the E15 binned excitation plot → **E25
    `fig_gs_decomposition.png`** (GS-projected occupation = e-h channel) beside
    the E15 loss function (plasmon). Equations now GS-projected occupation (left)
    + loss function & plasmon frequency (right). NOTE in speaker notes that the
    grey untracked region still needs the staged E25 re-run to fill (see
    2026-06-01 entry).
  - **S10 wave nature:** replaced the combined GIF with 3 panels — animated
    `assets/xz_density.gif` + static report screens
    `fig_leed_transmission_fft.png` (transmission, structure-factor overlay) +
    `fig_leed_backscatter_centre.png` (backscatter).
  - **References (A3):** rewrote all 20 to **full APA**, alphabetical, two
    columns, **8.5 pt**. 15 sourced from `draft5/references.bib`; 5 not in the
    bib added from literature (Ashcroft & Mermin, Giuliani & Vignale, Boudaïffa
    2000, Alizadeh 2015, Pimblott 2007). **Alizadeh-2015 and Pimblott-2007 page
    ranges flagged — user to confirm.**
  - **Appendix A4:** replaced the orphaned E15 charge-recovery slide with a new
    **wake** slide — assertion "deeper hole behind / larger accumulation ahead",
    `wake_2d.gif` + `wake_1d.gif` (copied into `assets/`), plus 5 qualitative
    observations (sweeps not robust; feature persists toward classical/low-σ;
    expected lower quantum SP; higher-σ hole moves faster; does density agree
    with stopping power?).
- **Read-aloud script:** `drafts/storyline3_script.md` — word-for-word, per slide,
  paced to ~10 min, with stage directions; appendix marked Q&A-only.
- **A5 appendix (deck now 22 slides):** added `fig_master_stopping_vaxis_lowdens.png`
  (copied to `assets/`) with a **placeholder TBD headline** — user will write the
  assertion later.
- **Co-supervisor: Runfeng Zhou — CONFIRMED** (user, 2026-06-02). Spelling kept as
  a **practice-check flag** in the title + A1 speaker notes and the script (user
  asked to verify it while practicing).
- **Still open:** confirm the two flagged references (Alizadeh-2015, Pimblott-2007
  page ranges); write the A5 assertion; the exact-zero E25 decomposition still
  needs the staged re-run.

## What changed
- Built `docs/presentations/` workspace (full clean reorg + lazy phase folders).
- Renamed/rewrote the method doc: `sample-structure.md` →
  `method/assertive-evidence-presentation.md` (skill-shaped, with frontmatter).
- Updated `requirements.md` (venue-specific) and created
  `evaluation/rubric-checklist.md` (from the marking sheet).
- Recorded specs cross-checked against the **official AE PowerPoint template**
  in `templates/` and Alley 2013 / PSU checklist.

## Files touched (absolute)
- /local/data/public/skcb2/tddft/docs/presentations/method/assertive-evidence-presentation.md (new skill)
- /local/data/public/skcb2/tddft/docs/presentations/requirements.md (rewritten)
- /local/data/public/skcb2/tddft/docs/presentations/evaluation/rubric-checklist.md (new)
- Moved: 3 source PDFs → `reference/`; AE source PDFs+txt → `method/assertive-evidence/`.
- Memory: feedback_presentation_skill_building.md (+ MEMORY.md index line).

## Workspace layout
```
docs/presentations/
  requirements.md                        venue constraints + maths/lit mandate
  method/assertive-evidence-presentation.md   the SKILL (general, no rubric)
  method/assertive-evidence/             source PDFs + assertion-evidence.txt
  brainstorms/brainstorm1.md             user's (problem statement, in progress)
  evaluation/rubric-checklist.md         scorable marking criteria
  drafts/  assets/                       (decks, regenerated figures)
  reference/                             handbook, marking sheet, Report1 PDFs
  templates/                             two official AE .pptx templates
```
Other phase folders (`storylines/`, `outlines/`, …) created lazily.

## Key decisions locked
- Skill is **venue-agnostic**; rubric/limits live only in requirements + evaluation.
- 7 workflow phases (added optional **Preempting questions**, phase 6).
- 7 slide types (added optional **Background/importance** type 2; **References**
  + **anticipated-question** slides are *backup/appendix after the conclusion* —
  talk ends on conclusion which stays up during Q&A).
- Typography: title 36 pt bold, headline 28 pt, body 18–24, refs 14; bold
  sans-serif (Calibri/Arial); white bg, black text.
- ≤20 wpm is a **warning, not a blocker**.
- **Citation rule (personal):** Physics-APA footer ref, numbered by first
  occurrence, References slide after conclusion.
- **Maths + numerics + background-literature slides are mandatory** (rubric).

## Trusted sources used
- Alley, *The Craft of Scientific Presentations* 2nd ed. (2013), Ch.4.
- Penn State AE checklist (writing.engr.psu.edu/AE_checklist.pdf); the official
  PSU AE PowerPoint template (theme + slides + notes parsed directly).
- Report1 PDF (content), Assessors_Marking_Sheet.pdf, student-handbook.pdf.

## Known issues / open points
- "Physics-APA" citation style under-specified (numeric vs author-date) — to
  refine with user.
- Slide **toolchain not chosen** (Beamer vs PowerPoint/Keynote) — deferred to the
  slides phase; PDF is the final artefact.
- Folder structure flagged as possibly over-engineered — **review & prune at end
  of workflow** (user's explicit request — remind them).

## Assumptions still in play
- Pure white background accepted (cream parked).
- Content brainstorming is user-led; agent supports method/infrastructure.

## Exact next steps
1. Continue content via phases: brainstorm → storyline → outline.
2. When concrete artefacts appear (strong assertions, working figures), **prompt
   the user to document them** in the skill's "Sample assertions" section.
3. At slides phase: choose toolchain; apply spec table; wire citations.
4. Score each draft with `evaluation/rubric-checklist.md`.
5. End of workflow: review/prune folder structure and finalise the skill.

---

## Milestone 2026-06-02 — Draft-2 figure fixes A1–A7 (complete)

Disambiguated via `/grill-with-docs`; spec in
`docs/plans/presentation_figure_fixes_A1_A8.md`. All seven figure tasks done.

**Output (all in `docs/presentations/assets/draft2_fixes/`):** one `aN_*.py`
script + PNG per task; `INDEX.md` maps task → figure → known-case result.
Report1 pipeline left untouched (new scripts are self-contained adaptations).

**Done & known-case-verified:**
- A1 `a1_free_wp_combined.png`: ΔE replaces abs-energy; σ(t)+ΔE stacked shared-x;
  |ΔE|<1e-4 meV in IFW. σ_r(0)=SIG0/√2 ✓.
- A2 `a2_coronene_setup.png`: k₀ arrow at WP entry (x=0, z=+12), label outside.
- A3 `a3_gs_decomposition.png`: charge balance moved BELOW (horizontal bars,
  +496.8 untracked unclipped); Σ=−5.6e−17 ✓.
- A4 `a4_coronene_target{,_linear}.png`: report style, fixed dims; ∫n_2D=108.000 e ✓.
- A5 `a5_jellium_schematic.png` + `a5_metals_rs_table.png` (Cs bold, footnote
  r_s=5.69); n_sim→r_s=5.69 consistency ✓.
- A6 `a6_momentum_{1d,2d}.png`: signed-k_z marginal + full negative-k_z 2D; norm
  conserved ✓. Backward weight t0=2.1% (Gaussian tail) → t_end=3.1%, net +0.95%.
- A7 `a7_stopping_{KSenergy,expectation}.png`: two defs. **KS-energy S2>0 all E
  (0.64→0.26); expectation S1<0 all E (−0.072→−0.034)** — WP kinetic energy rises
  via σ_p² spreading. Loss-fn(red)+classical+Lindhard box-q on both.

**Awaiting USER VERDICT** (verification mode — agent does not judge): which WP
stopping definition is "closer to truth" (A7); whether to keep Lindhard on A7
(kept by default). 

**Slide-level, NOT figure (user action in PPTX):** A2 floating formula box +
width/energy boxes are PowerPoint shapes — delete on the slide. A8 resolved:
Lindhard box-q is the analytical anchor (no SRIM/Bethe/Echenique).

**Unverified:** figures not visually previewed by agent (per user rule — user previews).

---

## Milestone 2026-06-02 (late) — draft3_freewp batch

New folder `docs/presentations/assets/draft3_freewp/` (own INDEX.md). Plan:
`docs/plans/presentation_freewp_stopping_equations.md`. Grounding:
`docs/sources/stopping-power-formulae.md`.

**Done & known-case verified (run venv python from repo root):**
- 5 equation PNGs (`render_freewp_equations.py`, true usetex): eq_spreading, eq_bethe,
  eq_bloch (full Bethe+Bloch digamma), eq_lindhard, eq_ks_stopping. LHS=S(v),
  explicit constants. Smoke test PASSED.
- `freewp_sigma_t.png`: σ_r(0)=3.5355=σ₀/√2 ✓ (IFW t≤3.5; barely spreads, σ₀=5 large).
- `freewp_energy.png`: absolute E(t), E(0)=100.8 eV, max|ΔE|=1e-4 meV ✓.
- `freewp_xz_density.gif`: xz slice, **z horizontal / x vertical**, density_wp norm=1 ✓,
  fixed log clim, 18 frames.
- `freewp_total3d.gif`: pvbatch volume render (log), density_wp, 18 frames (EGL warns benign).
- `stopping_KSenergy_sigmasweep.png`: adapted from make_fig_master_stopping_vaxis.py —
  REMOVED red loss-fn + grey analytical Lindhard; KEPT classical Ehrenfest; WP σ=0.5/1/3/5/8
  under KS-energy (S2) def; legend ascending + Classical. S2: σ0.5=1.223, σ1 0.637→0.256,
  σ3=0.086, σ5 0.130→0.012, σ8=0.004; classical 0.720→0.053.
- RECREATED coronene `../xz_density.gif` with z-horizontal/x-vertical orientation
  (`make_coronene_anim_gifs.py` build_xz edited; in place, decks still reference it).

**Decisions:** energy=absolute E(t); xz/3D=GIFs; eq LHS=S(v) explicit constants
(switchable to a.u. on request); Bloch=full correction; orientation z-horizontal.

**Field note:** free-WP density_total/_system integrate to 3.0 (WP + faint ~2-e diffuse
background artifact); used density_wp (norm=1) for clean gifs/3D.

**Unverified:** figures not visually previewed by agent (user previews per rule).
**Default flagged for user:** equation RHS uses explicit constants, not atomic units.

### 2026-06-02 (latest) — revisions + concentrated σ=1 + coronene xz draft3

- `freewp_sigma_t.png` / `freewp_energy.png`: removed green overlay + "interference-free"
  annotation. `freewp_energy.png` now plots total + directed (⟨p⟩²/2m) + spread
  (σ_p²/2m); closure |total−(dir+spr)|=5e-10 eV ✓ (t=0: 100.82 = 100.00 + 0.82 eV).
- `freewp_sigma1_xz_density.gif` + `freewp_sigma1_total3d.gif`: concentrated σ=1 packet
  EXPANSION (run_free_wp_L50_E25_sigma1_v2, density_wp norm=1 ✓, t≤5 a.u., σ_d 0.71→3.8,
  26 frames each). E=25 (only σ=1 free run with clean density_wp; E100 σ=1 has only
  density_rt_total).
- `coronene_xz_density.gif` (draft3, styled w/ colorbar) + in-place `../xz_density.gif`:
  coronene xz in new schema (z horizontal / x vertical), 61 / 21 frames.

### 2026-06-03 — INQ workflow diagram: contribution page + Bohr eq

- `eq_bohr.png` added to `draft3_freewp/` (render_freewp_equations.py); Bohr 1913
  classical stopping; recorded in docs/sources/stopping-power-formulae.md.
- `Misc/INQ-flow-chart.drawio`: spliced a 5th page `contribution` (id
  contrib_inqkit_inqview) via `Misc/build_contribution_page.py`. 16:9 (1920x1080),
  slide-level. GREY spine = INQ base (Input → GS SCF → real_time::propagate);
  GREEN = my work in two lanes — inqkit (in-run C++: WP injection/orthonormalisation/
  validation, per-step callback, field/VTI/CSV/manifest I/O) and inqview (post-run
  Python: direct observables densities/wavefunctions + derived loss fn/momentum/
  stopping). Legend grey vs green. Existing 4 pages untouched. Backup at
  `Misc/INQ-flow-chart.drawio.bak`. XML validated (5 pages, well-formed); not
  PNG-rendered (no drawio CLI) — export the `contribution` page from the drawio app.
- Encoding gotcha (fixed): drawio stores HTML labels XML-escaped (`&lt;font&gt;`),
  and HTML entities must be literal chars (—, spaces) not `&nbsp;`/`&mdash;` to avoid
  over-escaping. Generator handles this via xml_esc().

### 2026-06-03 (cont.) — magic table, induced-Δn gifs, workflow PNG report-grade

- `magic_numbers_table.png` (magic_numbers_table.py): report tab:magic-numbers
  (|G|²=0..6, cumulative N, N=162 highlighted). cum N=[2,14,38,54,66,114,162] ✓.
- `wp_induced_density_diff.gif` + `classical_induced_density_diff.gif`
  (induced_density_diff_gifs.py): report-standard induced-Δn (E25 σ1). WP =
  density_delta_jell − density_delta_free (removes WP self-density); classical =
  density_delta_classical. Shared symmetric clim 1.01e-3 (99.5 pct) across BOTH;
  RdBu_r, TwoSlopeNorm, gaussian σ=1, z-horizontal, apply_style serif, 28 frames,
  t=0..11 a.u. t≈0 frames exactly 0 ✓. Method = report draft5 make_fig_density_diff_2d.py.
  The "WP − classical" (_2d_diff) panel intentionally NOT produced (user request).
- `INQ-flow-chart-contribution.png` redone to report standard (serif Computer-Modern
  family, balanced full-canvas 1920×1080, callback→{wavepacket,I/O} tree, refined
  arrows/palette). Iterated with preview (user granted preview permission for this PNG).
  Renderer: Misc/render_contribution_png.py.

### 2026-06-03 (cont.) — draft3_wp rename, gif extension + stutter fix, workflow text

- `draft3_freewp` -> `draft3_wp` (all results; paths sed-updated in scripts +
  Misc/render_contribution_png.py).
- σ=1 free gifs extended t≤5 -> t≤11 a.u. (freewp_sigma1_xz/_3d; 56 motion frames).
- Loop stutter (abrupt last->first jump) fixed on ALL evolution gifs via
  draft3_wp/_gifpad.py write_gif: duplicated start/end frames (~0.4s/0.9s holds);
  GIF writer merges them into longer-duration end frames (PowerPoint-safe). Verified:
  sigma1 xz first=500ms last=1000ms mid=120ms, total 7.98s.
- INQ-flow-chart-contribution.png: enlarged fonts (TITLE_SZ 24, SUB 17, MOD 13; grey
  titles tsize 21 + "Real-time TDDFT" shortened; container 19.5; main 28; legend own
  row). Box titles ~17pt on-slide (was ~13). Iterated with preview.

### 2026-06-03 (cont.) — sigma1 gif trim + workflow font/annotations

- freewp_sigma1_xz_density.gif & freewp_sigma1_total3d.gif: T_START=0.5 a.u. — first 3
  frames (t=0/0.2/0.4) dropped; the sub-2-grid-cell packet rendered as a blocky
  "initial disturbance" (data was clean, norm=1; purely a rendering artifact). Gifs
  now start t=0.60; colour scale taken from first kept frame.
- INQ-flow-chart-contribution.png: font -> Calibri (Carlito, metric-compatible clone;
  Calibri not installed). Removed the inqkit::/inqview module-name annotation lines
  from every block (user: no value). Renderer Misc/render_contribution_png.py.
