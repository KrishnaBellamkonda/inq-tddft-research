# Handover: Localised jellium GS — parameter study + analytical mental models

Campaign: `/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium/ground_state_parameter_study.md`
Gates Campaign 1: `/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium/classical_projectile_fix.md`

---

## Milestone: 2026-06-28 — dense re-run + notebooks + EXPERT PANEL + summary (all 5 user stages done)

**Stages 1-5 complete.** (1) 7 phase/campaign notebooks built; (2) denser grids re-run
(no failures); (3) notebooks refreshed; (4) scientific-panel run (4 opus experts +
judge, ~448k tok); (5) `campaign_summary.ipynb` built (verdict + plots + conclusions).
Notebooks: `hypotheses/campaign_autorun_study/` (H0-H5 + campaign + summary; builders
`build_notebooks.py`, `build_summary.py` — re-run to refresh).

**PANEL VERDICT (provisional, user owns it) — the science conclusion:**
- **Use PBC (periodicity-3). Energy reference for Campaign 1: E_SIE ≈ 4.4 eV** (2 s.f.;
  three-way agreement: PZ model 4.4 / far-launch 4.55 / PBC-at-r=40 4.26), quoted as a
  fixed-geometry, fixed-BC difference. r→∞ may drop ~0.3-0.8 eV (not yet converged).
- **DISCARD the open-z (periodicity-2) SIE for the charged WP.** Negative SIE is excluded
  by positivity of E_H+E_xc; it is the `0.5·rc²` (rc=L_z) G=0 monopole term applied to
  the net-charged WP cell (poisson.hpp:49), absent from the neutral GS. The open-z
  "more-negative-with-r" trend is real E_cross(r) (slope ≈ −0.046 eV/Bohr, present in
  BOTH BCs) sitting on the spurious offset — NOT a self-energy.
- Absolute GS energies (−161/−109/+60 Ha) are a box/BC gauge, not physics (interior n₀,
  N correct). Open-z is clean for the NEUTRAL GS; wrong only for the charged difference.
- **Unresolved:** 3× monopole magnitude gap (~2 eV envelope vs ~6 eV observed; likely INQ
  FFT normalisation, unverified); 4.3 eV convergence (r=40 only 7.5 Bohr from edge).

**DECISIVE NEXT TEST (panel, cheap, not yet run):** isolated σ_WP=0.5 electron in an
EMPTY box, {periodicity 3,2}×{L_z=90,120} (4 configs, seconds each) + in-plane face scan
(L_xy=50/70/100). PBC→≈4.4 eV & L_z-flat confirms the SIE; open-z→∝L_z² proves the
monopole. Locks the number before it goes to Campaign 1.

**Open questions for the user (panel §7):** pure SIE vs full t=0 excess E_SIE+E_cross(r);
launch r; restore ghost ∫v_ghost·n₊ before WP-vs-classical comparison; resolve the 3×
gap analytically or accept the L_z² test; run the 4-config vacuum control to lock the number.

---

## Milestone: 2026-06-28 — denser re-run launched + notebooks built; panel+summary pending

User asked: denser data for every phase, re-run, refresh notebooks, expert panel, summary notebook.

**Denser grids (user-approved) — orchestrator + analyse made grid-agnostic (dynamic r):**
- H0 (NEW phase in orchestrator): r={4,12,20,28,36,40}, wp+cl, periodicity 3, L_z=120.
- H1 w={0,0.25,0.5,0.75,1,1.5,2,3}; H2 Lz={50,60,70,80,90,110,130,150}+open-z@{90,120};
  H3 a={5,7.5,10,12.5,15,17.5,20,22.5,27.5} (N scaled); H4/H5 r={2,6,10,14,18,22,26,30,34,40}×{per3,per2}.
- **RE-RUN LAUNCHED** headless (orchestrate.py, PID 1570086, GPU1, `orchestrate_dense.log`).
  Idempotent: skips the ~25 existing runs, computes ~14 new GS + fast projectile points (~3-5h).
  H0 done (4.1 min, emailed). User GPU0 run untouched.

**Notebooks BUILT (step 1):** `hypotheses/campaign_autorun_study/` — `build_notebooks.py`
generates 7 .ipynb (H0-H5 + `campaign_autorun_study.ipynb`) with hypothesis→setup→
results→takeaway + **embedded PNG plots** (base64; no kernel needed — venv jupyter CLI
is pyenv-shadowed, so notebooks embed plots rather than execute). Re-run the builder to REFRESH.

**Continuation:** waiter (bg id b4vesq4v9) fires when PID 1570086 exits → then:
(1) refresh notebooks (`build_notebooks.py`), (2) **expert panel** via `scientific-panel`
skill — 3 experts: computational condensed-matter, TDDFT/jellium, codebase+this/related
campaigns, (3) **summary notebook** (ideas/conclusions/problems + key plots).

Usable result so far: PBC E_SIE ≈ 4.3 eV (open-z still needs the net-charge G=0 fix).

---

## Milestone: 2026-06-27 (session 3) — ALL PHASES H0–H5 executed + emailed; H4/H5 path bug fixed

### Autonomous run complete (H0–H5 + cumulative emailed). Frontmatter H0–H5 = done.
- **Pipeline bug found+fixed:** projectile `run.cpp` writes `results/<LJ_OUT>/raw/...`
  (nested), but the orchestrator `_done()` + `analyse_phase.py` expected
  `rundir/results/` and `rundir/raw/`. Sims SUCCEEDED (rc=0) but were mis-marked FAIL
  → re-run → analyse crash (`FileNotFoundError observables.csv`). **Fix:** both
  `_done()` and `analyse._find_obs()` now **glob** (`**/run_summary.txt`,
  `**/observables.csv`) → layout-tolerant. No sim rerun needed (data was intact);
  idempotent resume then SKIPPED all 16 projectile sims and re-ran only the analyses.

### Headline results (PROVISIONAL — flags below)
- **H4 critical path: PBC (periodicity-3) E_SIE plateau = 4.3 eV ✓** (matches known ~4.5).
- **H4 open-z (periodicity-2) E_SIE = −2.1 eV — UNPHYSICAL.** New finding: a net-charged
  WP cell under periodicity-2 carries a G=0 compensation term (`0.5·rc²`, 2D kernel)
  absent from the NEUTRAL GS, so naive `E_wp − E_GS` is biased for open-z. **The open-z
  energy reference needs the net-charge G=0 correction before the PBC-vs-open-z verdict
  is trustworthy.** Trust the PBC E_SIE (4.3 eV) for now.
- H0: raw gap artifact-dominated (done earlier). H1: edge model. H2: GS conv + open-z
  GS usable. H3: thickness/σ_s (RAW, E_self caveat). H5: classical image error −1.2 eV;
  route-2 ghost-bg flagged.

### Open follow-ups (flagged in emails, for review / next session)
1. **Open-z net-charge G=0 reference** (H4/H5 periodicity-2) — fix before the BC verdict
   and before handing the energy reference to Campaign 1.
2. Φ extractor (H2), σ_s E_self correction (H3), ghost-background integral (H5) — the
   pre-gates, still best-effort/flagged.

### Pipeline now robust (Python orchestrator)
`orchestrate.py` (glob `_done`) + `analyse_phase.py` (glob `_find_obs`) are
layout-tolerant + idempotent-resumable. Re-run any subset: `orchestrate.py H4 H5 CUMULATIVE`.

---

## Milestone: 2026-06-27 (session 3) — switched to PYTHON orchestrator (bash master retired)

**User decision:** autonomous executor must be Python, not bash (bash autonomous
scripts are brittle). Bash `master.sh` stopped (mid binary-build, before H1; user's
GPU0 qsp_phase5 run protected/untouched). Replaced by:
**`scripts/campaign_autorun/orchestrate.py`** — RUNNING headless (GPU1, nohup),
currently in H1 (gs w=0). Features: structured logging (`orchestrate.log`),
**idempotent resume** (skips runs with `run_completed = true` — safe to re-launch),
per-phase `try/except` + full-traceback **failure emails** (chain continues),
one-shot sim retry, reuses `analyse_phase.py` (+ `email-notifications`).
Run with: `GPU=1 nohup venv/bin/python3 orchestrate.py &` (phase filter: pass phase
names as argv, e.g. `orchestrate.py H4 H5`).

Pipeline files (all built/validated): parametrised `gs/run.cpp`, `wp/run.cpp`,
`classical/run.cpp` (env-driven geometry + periodicity 2/3); `analyse_phase.py`
(H1–H5 + per-phase email); binaries `gs/run`,`wp/run`,`classical/run` compiled.
periodicity-2 (open-z) VALIDATED via smoke (N=42, interior n₀ ok).

**Preference codified** in `campaigns` (autonomy checklist) + `tddft-simulations`
(Phase 6) skills: Python orchestrator over bash for autonomous multi-phase runs.

Order now running: H1 (edge-w) → H2 (Lz + open-z) → H3 (thickness) → H4 (WP
PBC-vs-open-z, E_SIE) → H5 (classical mirror) → cumulative; email per phase.

---

## Milestone: 2026-06-27 (session 3) — H0 COMPLETE + emailed; machinery + email channel validated

**H0 done (frontmatter `done: true`).** 4 runs (WP/classical × r=4,40) in L_z=120 box,
k₀=0; analysed + plotted + **emailed** (`hypotheses/h0_base_difference/analyse_h0.py`,
`H0_base_difference.png`). Email channel confirmed working (creds at
`~/.config/inqview/gmail_credentials.json`).

**H0 RESULT (GS = −108.534 Ha):**
| run | excess vs GS | reading |
|---|---|---|
| WP r=4 | +87.6 eV | localisation+SIE |
| WP r=40 | +85.9 eV | localisation+SIE (distance-STABLE) |
| cl r=4 | +187.5 eV | unscreened ghost (near) |
| cl r=40 | +11.2 eV | unscreened ghost (far) |
- WP zero-point KE measured = `E_kin(WP)−E_kin(bath) = 3.00 Ha = 81.7 eV` ✓ (=3/4σ²).
- **E_SIE (route 1, k₀=0, r=40) ≈ 4.2 eV** ✓ (matches known ~4.5 eV).
- Raw WP−cl gap flips sign: r=4 **−100 eV**, r=40 **+75 eV** → artifact-dominated,
  NOT the localisation energy. Hypothesis answered: NO. Motivates ghost-bg (H5).

**Validated for the headless run:** GS pattern, WP injection+energy (core of H4),
classical ghost, observables step-0 readout, plotting, email. **Still UNTESTED:
periodicity-2 (open-z)** — used by H2/H4/H5; must smoke-test before those run headless.

**False alarm corrected:** a mid-run −108.12 value I flagged as a WP bug was the
classical r40 run, not the WP. WP machinery is correct.

---

## Update: 2026-06-27 (session 3) — L_z=120 GS validated; absolute-energy box-dependence CONFIRMED (affects H2/H3)

H0 GS (slab_n82_L50x50x120, periodicity 3, w=0) built + converged. Validation:
- **interior n₀(|z|<6) = 1.317e-3** vs target 1.312e-3 (0.4%) ✅
- **∫n d³r = 82.00** ✅ (correct count); peak n(z)=1.48e-3 at z≈−7.2.
- **absolute E_GS = −108.53 Ha vs −160.99 Ha at L_z=90** — a 52 Ha shift. This is the
  **box-dependent background self-energy** (charged slab under G=0 uniform compensation;
  worksheet Part 3.4 "absolute total energy is meaningless, only fixed-geometry
  differences survive"). **NOT a bug** — interior density proves the GS is sound.

**CONSEQUENCE (important, must handle in H2/H3):** absolute E_GS is NOT box- or
geometry-independent. So:
- **H0/H4/H5** (fixed L_z=120): `E_total(0) − E_GS` differences are CLEAN (E_self cancels). ✅
- **H2 (Lz sweep):** do NOT compare absolute E across L_z — use Φ (a potential difference)
  for the work-function plateau, not total energy.
- **H3 (thickness/σ_s):** the liquid-drop fit E(N) across thicknesses changes the background
  slab ⇒ E_self differs per thickness; the σ_s extraction must account for the
  box/background-dependent E_self (worksheet E_self term) or use a consistent reference,
  else σ_s is corrupted. Pre-gate the σ_s extraction with this in mind.

Dispatcher continuing: WP{r4,r40} then classical{r4,r40}. H0 email on completion.

---

## Milestone: 2026-06-27 (session 3, cont.) — AUTONOMOUS execution armed + email-notifications skill created

### Decisions (user)
- **Run the whole ladder H0→H5 autonomously**, no user in the loop, in the locked
  order. Each phase: build → run (GPU) → analyse → highlight plot → **email** →
  flip task done + handover → launch next. Campaign `status: running`. Documented in
  the campaign `<autonomous_execution>` section.
- **Per-phase email is mandatory** with a plot and the four parts: hypothesis being
  tested, what was done, what the plot shows, conclusion. → new **`email-notifications`
  skill** (`.claude/skills/email-notifications/SKILL.md`) owns this contract + the send
  mechanism (`inqview.email.send_run_email`, Gmail SMTP — MCP only drafts).
  `tddft-simulations` now points to it.

### State of the autonomous run
- **H0 is RUNNING** (background dispatcher `run_h0.sh` on GPU1; GPU0 = user's
  qsp_phase5 run, untouched). Sequence: GS(L_z=120) → wp{r4,r40} → cl{r4,r40}.
- On H0 completion (background notification): validate GS energy (≈ −160.99 Ha),
  read the four E_total(0), compute raw gap + ghost-bg-corrected gap, make the H0
  highlight plot, **email per email-notifications**, flip H0 task done, then build +
  launch H1 (edge-w sweep). Continue the chain to H5.

### Email-sending mechanism (reference)
`inqview.email.send_run_email(subject, body, attachments=[png], to="chiddukanna@gmail.com",
html_body=...)` run via `/local/data/public/skcb2/tddft/venv/bin/python3`. Needs ≥1 plot.

---

## Update: 2026-06-27 (session 3, cont.) — H0 added as the motivating rung-0

**H0 (base WP-vs-classical E_total(0) gap).** User addition: the simplest projectile
question — is `E_WP(0)−E_classical(0)` equal to the WP localisation energy (81.6 eV)?
- WP (σ_WP=0.5, k₀=0) vs matched classical (`electron_gaussian_wpsigma0p5`, k₀=0), same
  position, t=0 energetics; compare E_WP(0), E_classical(0), E_GS.
- **Distances r=4 AND r=40** (near+far) → r=40 (z=−52.5) exceeds the L_z=90 box, so
  **H0 runs in the L_z=120 box (w=0, periodicity 3)**: 1 new GS + 4 tiny runs. H0 thus
  seeds run-set 2 (first slice of the H4/H5 matrix at k₀=0).
- **Predicted finding (from Phase-3 numbers):** raw gap is NOT the localisation energy —
  it's ghost-bg-artifact-dominated (E_cl−E_GS ~ +700–800 eV near from the unscreened
  `∫v_ghost n_GS`; E_WP−E_GS ~ +86 eV), so E_WP−E_cl ~ −650 eV near. Only after the
  ghost-bg correction does the gap reduce to localisation+SIE ≈ 86 eV. → motivates H5.
- Ladder now **H0 → H1 → H2 → H3 → H4 → H5 + cumulative** (frontmatter tasks 0/9).
- Note: H0's L_z=120 GS is at w=0; H4/H5 use the clean-w L_z=120 GS (from H1) — at most
  one extra GS if H1 finds w=0 unacceptable for energetics.

**Next:** build H0 = L_z=120 config header + GS run.cpp (w=0, periodicity 3) + WP &
classical projectile runs at r=4,40 + `H0_base_difference.ipynb`.

---

## Milestone: 2026-06-27 (session 3, cont.) — campaign reframed as an ordered hypothesis-ladder + notebook contract

### Decision (user): organise as a logical complexity ladder of falsifiable hypothesis-sets, GS-foundation first; each set its own notebook + a cumulative campaign notebook
Locked order (`<hypothesis_ladder>` + `<notebook_contract>` in the campaign file;
frontmatter `tasks` reframed to match):

1. **H1 edge model** — w∈{0,.5,1,1.5,2}, GS-only, periodicity 3 → clean w (Gibbs vs
   Friedel). *Upstream of every later GS.* → `H1_edge_model.ipynb`
2. **H2 GS convergence + open-z** — Lz∈{50,70,90,120} + periodicity-2 GS smoke, clean w
   → Lz, Φ, open-z viability. → `H2_gs_convergence.ipynb`
3. **H3 surface energetics** — thickness a∈{7.5,12.5,17.5,22.5} → σ_s, e_bulk,
   charged-plates (thread A). → `H3_surface_energetics.ipynb`
4. **H4 WP energetics** — WP E_total(r) {4,16,28,40}×{3,2} + k₀ control → PBC-vs-open-z
   verdict + E_SIE route 1. → `H4_wp_energetics.ipynb`
5. **H5 classical subtraction** — classical mirror {4,16,28,40}×{3,2} + ghost-bg →
   route-2 E_SIE cross-check + thread-D cutoff; **pins E_SIE/reference for Campaign 1**.
   → `H5_classical_subtraction.ipynb`
6. **Cumulative** `campaign_cumulative.ipynb` — highlight plots → overall story.

Notebook house style per `notebook-making` (hypothesis stated up front → outcome plots
→ explicit verdict); ADR-0007 placement `hypotheses/<NN_purpose>/`; canonical theme.

### Why GS-first (refines the earlier critical-path-first spine)
Clean w (H1) feeds run-set 2's GS; the surface dipole/charged-plates picture (H3) feeds
the E_cross interpretation in H4/H5. Gate (E_SIE) still delivered, at H4–H5, on a
validated GS.

### Exact next steps (updated)
1. **Stage-5 autonomy checklist** (deferred by user until the ladder ordering was set —
   now it is) OR start building H1 (cheapest, upstream): config header + GS run.cpp +
   w-sweep dispatcher + `H1_edge_model` notebook.
2. periodicity-2 GS smoke test is the H2 gate for all open-z runs.
3. Pre-gate Φ/σ_s (before H2/H3) and ghost-bg integral (before H5).

---

## Milestone: 2026-06-27 (session 3, cont.) — route-2 classical subtraction LOCKED; reference pack corrected

### Locked route-2 spec (run-set 2 critical path — pins E_SIE for Campaign 1)
- Matched ghost **`electron_gaussian_wpsigma0p5.upf`** (charge std σ_WP/√2=0.354;
  prefer over legacy `sigma0p35`), k₀=0, t=0 energetics.
- **Full mirror: r∈{4,16,28,40} × periodicity {3,2} = 8 classical runs.**
- **MANDATORY ghost–background correction** (re-added in analysis):
  `E_SIE = E_WP − [E_cl + ∫v_ghost·n₊] − ⟨T_WP⟩`.
- **Three deliverables:** (i) E_SIE(r) cross-check vs route-1 plateau ≈4.5 eV;
  (ii) E_cross^WP == corrected E_cross^cl (WP–slab interaction purely classical);
  (iii) classical periodicity 3-vs-2 image error → Campaign 1 cutoff (thread D).
- **Pre-gate:** ghost–bg integral helper (formula-validation: numeric vs closed-form erf).
- **Coarse run-set 2 total = 20 INQ runs** (2 GS + 8 WP + 2 WP control + 8 classical).

### PHYSICS CORRECTION (important — reference pack + plan + docx updated this session)
Earlier I wrote "two equivalent fixes (re-add ghost–bg OR launch far)". **WRONG —
they are NOT equivalent.** Launch-far alone fails for the classical: the bare
ghost–electron Coulomb decays only as ~N/r (≈56 eV even at r=40), because the
chargeless ghost feels bare electrons but its `∫v_ghost·n₊` background-attraction is
dropped (z_valence=0). The WP, a real electron, sees the *neutral* slab (gets both
terms) → ~0.5 eV. So **ghost–bg re-add is mandatory at every r.** Fixed in
`reference_pack.md` C.6, `worksheet_plan.md` C6, and regenerated the docx (254 eqns).

### Run-set 2 now fully designed (coarse). Still open: rung-2 finer sweep (after BC
verdict); synthesis (threads A/E); Stage-5 autonomy checklist; build machinery +
periodicity-2 GS smoke test.

---

## Milestone: 2026-06-27 (session 3, cont.) — run-set 1 (surface-physics sweep) fully designed

### Locked run-set 1 spec (thread F, GS-only, baseline slab_n82)
- **Metrics:** profile (interior n₀, Friedel amp+λ=π/k_F≈9.3, spill-out, surface
  dipole, E_GS) + **work function Φ** (V_vac−E_Fermi) + **surface energy σ_s**
  (liquid-drop over the thickness series; 86.4 erg/cm² target).
- **Sweeps (one-at-a-time; Lx/Ly dropped = in-plane convergence, not surface physics):**
  - **w FIRST:** {0, 0.5, 1.0, 1.5, 2.0} Bohr → lock smallest w killing Gibbs while
    keeping π/k_F Friedel; use as baseline w for the rest (avoids chicken-and-egg).
  - **thickness a:** {7.5, 12.5, 17.5, 22.5} (N≈50/82/114/148 at fixed n₀, even) →
    σ_s + e_bulk + interior-n₀ breakdown.
  - **Lz:** {50, 70, 90, 120} → Φ plateau (reuse Lz=90 + run-set-2 Lz=120 GS).
- **Required pre-gate:** Φ and σ_s extraction helpers (NEW code) → code-test +
  formula-validation + catalogue BEFORE the sweep. Profile metrics need no new code.
- ~11–13 GS runs (GS-only, cheaper than run-set 2's WP runs).

### Still deferred (next grills)
- Run-set 2 rung 2+: finer r-sweep (after BC verdict) + **route-2 classical
  subtraction** (Hartree-matched / ghost–background term) ⇒ pin E_SIE for Campaign 1.
- Synthesis (thread A charged-plates writeup, thread E zero-point verdict).
- Stage 5 autonomy checklist + frontmatter `tasks` refinement.
- Build machinery (config headers, run.cpp, dispatcher, analyse.py); periodicity-2
  GS smoke test gates run-set 2.

---

## Milestone: 2026-06-27 (session 3) — run-set 2 COARSE experiment fully specified (E_total vs r, PBC vs open-z)

### Current status
The grill resumed and fully specified **rung 1 of the complexity ladder**: a coarse
`E_total(t=0) vs r` energetics experiment that doubles as a **method-selection test**
(full PBC vs open-z) for production. Not yet built or run. Engine capability for
open-z **verified in INQ source** (big de-risk). Remaining: build run machinery, then
the finer sweep + route-2 classical subtraction.

### Engine verification (INQ supports open-z natively — source line-refs)
- `inq/src/systems/cell.hpp:129-142` — `.periodic()`=periodicity 3, `.periodicity(2)`=slab,
  `.finite()`=0; periodicity 1 throws "not implemented".
- `inq/src/solvers/poisson.hpp:184-213` — `solve`/`in_place` dispatch on
  `cell().periodicity()`: 3→3d FFT, 2→`poisson_solve_2d`, else→0d.
- `poisson.hpp:41-53` `poisson_kernel_2d` opens the **z-axis** (`gpar=hypot(g0,g1)`,
  `gz=|g2|`, `rc=rlength()[2]=L_z`, box doubled only in z via `enlarge({1,1,2})`) →
  aligns with `SLAB_AXIS=2`. **This is thread D's cutoff, already in-engine** (Rozzi
  PRB 73, 205119 2006 for the 0d kernel, `poisson.hpp:57-67`).
- **periodicity 0 is WRONG for this geometry:** its `rc=min_rlength/2=50` Bohr < the
  ~67 Bohr slab+WP z-extent → would clip. Use **periodicity 2** (open z only).

### Locked coarse-run spec (rung 1)
- σ_WP = **0.5**; jellium **centered at z=0**; slab faces ±12.5; spacing 0.5.
- **Box L_z = 120** (z∈[−60,+60], 240 z-pts). Fresh slab GS converged in this box
  under **periodicity 3 AND periodicity 2** (2 GS).
- **r-grid (from near face): {4, 16, 28, 40} Bohr** (r_min = 8σ_WP = 4 Bohr;
  WP centroid z = −16.5, −28.5, −40.5, −52.5).
- **BCs compared: periodicity 3 (PBC, production default) vs periodicity 2 (open-z).**
- **WP momentum: k₀=0 (stationary)** for the sweep → ⟨T_WP⟩ = zero-point only
  = 81.6 eV (3 Ha). **k₀-independence control at r=4:** also run k₀=2.71 (100 eV)
  under both BCs (+2 runs) → confirm excess identical.
- Quantity: `excess(r,BC) = E_total(t=0) − E_GS(BC) − ⟨T_WP⟩` (= E_SIE + E_cross(r));
  `E_total` read after ~1–5 steps once stable (ETRS, dt≈0.01). E_GS = bare slab.
- **Run count:** 2 GS + 8 sweep WP + 2 control WP = 12 INQ invocations (WP runs are
  tiny 1–5-step; GS is the only substantial cost).
- **Done-criterion:** a verdict — keep PBC or switch to open-z for production —
  from the `excess(r,3) − excess(r,2)` curve; THEN a finer r-sweep.
- **No new analysis kernel** (excess is post-processing arithmetic) → no pre-gating.

### Validation plan (complexity-ladder smoke tests BEFORE the sweep)
1. **periodicity-2 GS smoke test** — no project GS has used periodicity 2 before;
   confirm the slab GS SCF converges under periodicity 2 and the interior n₀/energy
   match the periodicity-3 / L_z=90 GS (neutral-slab interior is box- and BC-independent).
2. **WP-under-periodicity-2 smoke test** — confirm WP injection + the net −1 charge
   give a sensible finite E_total under the 2d kernel (G=0 term handled by the kernel,
   `poisson.hpp:49`).
3. Abort on NaN / complex energy / GPU occupied.

### Files touched this session
- `/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium/ground_state_parameter_study.md` — run_sets updated with the coarse-run spec.
- (this handover)

### Assumptions still in play
- periodicity-2 GS converges for this slab (verify via smoke test 1 — not yet run).
- A fresh L_z=120 GS reproduces the L_z=90 interior n₀/energy (verify).

### Exact next steps
1. Decide: **build the coarse-run machinery now** (config header for L_z=120 +
   run.cpp supporting periodicity 3/2 + k₀ + r, dispatcher, analyse.py) — OR keep
   grilling (finer sweep; route-2 classical subtraction; run-set 1).
2. Smoke-test the periodicity-2 GS (validation step 1) before the 12-run matrix.
3. After the coarse verdict: lock the finer r-sweep and the route-2 classical
   subtraction (Hartree-matched / ghost-background term), then run-set 1.

## Milestone: 2026-06-27 — spine locked; B+C+D worksheet inputs produced for external worksheet agent

### Current status
Campaign still `status: draft` (not autonomy-ready). This session did two things:
(1) **locked the campaign spine** via the `/campaigns` grill, and (2) on user
redirect, **produced the inputs for the SIE/Coulomb/cutoff worksheet** that the
user's *separate* worksheet-building agent will turn into a polished worksheet.
The run-matrix numbers (grids, distances, cutoff value) are still deferred — the
session pivoted from locking run decisions to generating worksheet material before
those were set.

### Decisions locked this session
- **Campaign spine = "critical path first, one campaign"** (user-selected). Order:
  (1) cutoff prescription (thread D) → (2) run-set 2 elongated-z SIE (threads B+C)
  → hand B+D to Campaign 1; THEN (3) run-set 1 Lx/Ly/Lz/w sweep (F) → (4)
  charged-plates (A) + zero-point (E) → (6) synthesis. Respects the earlier
  "two run-sets = backbone" lock.
- **E_SIE extraction strategy (working assumption, reflected in the worksheet, NOT
  yet formally locked as a run decision):** BOTH routes, cross-checked —
  route 1 = far-launch/vacuum plateau; route 2 = Hartree-matched classical
  subtraction; their *agreement* is the campaign's falsifiable test. User
  redirected to the worksheet before selecting the AskUserQuestion option, but the
  reference pack + plan commit to this both-routes design.
- **Worksheet scope = full thread B + C + D** (SIE decomposition + Coulomb-vs-distance
  + long-range cutoff), one worksheet. **Resource format = Word (.docx)** for the
  authored pack; fetched sources stay PDF. **Resource content = authored reference
  pack + fetched open-access arXiv PDFs.**

### What changed (artefacts produced)
- New folder `docs/notes/localised-jellium-sie-worksheet/` holding the three
  worksheet-build inputs the user's downstream agent consumes:
  - `reference_pack.md` — authored "material": full B+C+D derivations, boxed
    results, in-repo numeric anchors, citations. (LaTeX math.)
  - `worksheet_plan.md` — the **markdown + XML-tag plan/spec** the worksheet agent
    executes (per-Part learning objectives, DERIVE-THIS problems mapped to
    reference-pack sections, VC callouts, self-test, output contract).
  - `resources/localised_jellium_sie_reference_pack.docx` — pandoc render of the
    reference pack, **248 native Word (OMML) equations** verified present.
  - `resources/arxiv_2307.03213_finite_size_stopping.pdf` — fetched (14 pp).
  - `resources/arxiv_1805.01377_rttddft_stopping.pdf` — fetched (8 pp).

### Files touched
- `/local/data/public/skcb2/tddft/docs/notes/localised-jellium-sie-worksheet/reference_pack.md`
- `/local/data/public/skcb2/tddft/docs/notes/localised-jellium-sie-worksheet/worksheet_plan.md`
- `/local/data/public/skcb2/tddft/docs/notes/localised-jellium-sie-worksheet/resources/localised_jellium_sie_reference_pack.docx`
- `/local/data/public/skcb2/tddft/docs/notes/localised-jellium-sie-worksheet/resources/arxiv_2307.03213_finite_size_stopping.pdf`
- `/local/data/public/skcb2/tddft/docs/notes/localised-jellium-sie-worksheet/resources/arxiv_1805.01377_rttddft_stopping.pdf`
- `/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium/ground_state_parameter_study.md` — resolved_decisions updated with the spine lock + worksheet-inputs pointer (this session).

### Commands run
```bash
# fetch open-access sources
curl -sL -o resources/arxiv_2307.03213_finite_size_stopping.pdf https://arxiv.org/pdf/2307.03213
curl -sL -o resources/arxiv_1805.01377_rttddft_stopping.pdf https://arxiv.org/pdf/1805.01377
# render docx with native equations
pandoc reference_pack.md -o resources/localised_jellium_sie_reference_pack.docx \
  --from markdown+tex_math_dollars --to docx --toc --toc-depth=2
# verify OMML equations: 248 blocks present
unzip -p resources/...docx word/document.xml | grep -o '<m:oMath>' | wc -l
```

### Key physics captured (repo-grounded; any worksheet must reproduce)
- Decomposition: `E_tot(0) − E_GS − ⟨T_WP⟩ = E_SIE + E_cross(r)`, with
  `E_SIE = E_H[n_WP] + E_xc[n_WP]` (Perdew–Zunger one-electron SIE, r-independent)
  and `E_cross(r) → 0` for a neutral slab.
- `⟨T_WP⟩ = ½k₀² + 3/(4σ_WP²)`; σ_WP=0.5,100eV → 181.6 eV (run 180.8). Using
  "+100 eV" overcounts SIE by ~82 eV.
- Net SIE ≈ **4.5 eV** (far-launch, σ_WP=0.5) vs bare self-Hartree
  `q²/(2s√π)` ≈ **22 eV** (XC cancels most). SIE ∝ 1/σ_WP → re-measure per σ.
- **Crux (thread C):** classical ghost t=0 excess = **+798 eV** (bare unscreened
  `∫v_ghost n_GS`; chargeless ghost omits the compensating `∫v_ghost n₊`), NOT
  comparable to the WP's +185.9 eV → naïve WP−classical subtraction fails; fix by
  re-adding ghost–background term OR launching both far.
- Cutoff window: `λ_TF ≈ 1.5 Bohr ≪ R_c < L_min/2 = 25 Bohr` (baseline); options:
  geometric (elongate+stop), truncated kernel, Martyna–Tuckerman. Finite-size
  plasmon-cutoff error ~8% (arXiv:2307.03213).

### Tests and validation
- Proposed: (worksheet self-test embedded in `worksheet_plan.md`); route1-vs-route2
  E_SIE cross-check (the campaign's falsifiable test).
- Approved/Run: none (no simulations this session; worksheet inputs only).
- Outcomes: docx equation render verified (248 OMML blocks). arXiv PDFs verified
  (valid PDF, correct page counts).
- Remaining gaps: docs/sources notes for the NEW external refs (Perdew–Zunger 1981,
  Martyna–Tuckerman 1999, arXiv:2307.03213) not yet written; the actual worksheet
  (built by the user's agent) does not yet exist.

### Trusted sources used
- Perdew & Zunger 1981 (PRB 23, 5048) — one-electron SIE.
- Lang & Kohn 1970/71/73 — jellium slab, work function, image plane.
- Martyna & Tuckerman 1999 (JCP 110, 2810) — nonperiodic Poisson.
- arXiv:2307.03213 (npj 2023), arXiv:1805.01377 — finite-size stopping.
- In-repo: `docs/sources/{nazarov-gross-2025,correa-2018,quijada-2007,stopping-power-*}.md`;
  `docs/notes/localised-jellium-theory.md` Parts 2/3/5/8.

### Known issues / blockers
- Campaign not autonomy-ready: run grids, elongated-box dims, distance grid r, the
  chosen cutoff value, and the surface-effect metric are all still deferred.

### Assumptions still in play
- The both-routes-cross-checked E_SIE design (reflected in the worksheet) is the
  intended run-set-2 strategy; needs a formal lock when the /campaigns grill resumes.
- The user's external worksheet agent consumes (a) the md+XML plan and (b) the
  .docx + PDF resources; format chosen accordingly (.docx for authored pack).

### Exact next steps
1. **User:** feed `docs/notes/localised-jellium-sie-worksheet/worksheet_plan.md` +
   `resources/` to the worksheet agent; review the produced worksheet.
2. Resume the `/campaigns` grill on `ground_state_parameter_study.md`: lock run-set 2
   geometry (50×50 in-plane to match Campaign 1), the distance grid r + elongation
   L_z axis, and the cutoff value (thread D).
3. Formally lock the E_SIE extraction strategy (both routes) in the campaign frontmatter.
4. Write docs/sources notes for Perdew–Zunger 1981, Martyna–Tuckerman 1999,
   arXiv:2307.03213 (campaign task 1 remainder).
5. Then proceed to run-set 1 (Lx/Ly/Lz/w sweep) design and the autonomy checklist.
