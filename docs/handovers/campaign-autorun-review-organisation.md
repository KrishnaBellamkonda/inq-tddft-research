# Handover: organise `campaign_autorun` run-set for review

**Task:** reorganise every notebook of the localised-jellium `campaign_autorun`
run-set (H0–H5) so the user can read each result and independently confirm it, with
NO assistant interpretation. Scope + design fixed via `/grill-with-docs` (2026-07-06).
Plan: `/local/data/public/skcb2/tddft/docs/plans/campaign-autorun-review-organisation.md`.

## 2026-07-06 — DONE (built + executed, 0 errors)

All artefacts live in
`/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/campaign_autorun_study/`.

### Deliverables (all executed to 0 errors, verified)
- **`00_index.ipynb`** — single canonical entry point: neutral ladder table
  [question | what was done | links], per-hypothesis highlight plot, one merged
  *⚠ Provisional (you own the verdict)* box. Supersedes both old aggregators.
- **`H0…H5_*.ipynb`** (6 study notebooks) — rebuilt by `build_notebooks.py` in the
  neutral order: **Question you were aiming to answer** (user's own `hyp` wording) →
  **What was done** (`setup`+`method`) → **Results** (plot + recomputed numbers) →
  **⚠ Provisional** box (the old `take`, quarantined). Plus an "Independently confirm"
  links block.
- **`runs/H0…H5_runs.ipynb`** (6 run-evidence notebooks) — new `build_run_evidence.py`.
  One per hypothesis; a table with EVERY run in the sweep (config + step-0/converged
  energy + excess-over-GS + interior n0) and the absolute `run_summary.txt` path per
  row. Verified: H0 shows classical excess 187→11 eV (decays), WP 87.6→85.9 eV (flat);
  H4 shows excess−ZP ≈ 5–6 eV per run.
- **`runs/rep_*.ipynb`** (6 representative single-run deep-dives) via the extended
  run-notebook assembler: `rep_H0_wp_r28_p3`, `rep_H0_cl_r28_p3` (H0 pair),
  `rep_H4_wp_r28_p2`, `rep_H4_wp_r28_p3` (H4 pair), `rep_H2_gs_lz120`,
  `rep_H3_gs_a15_N98`. GS reps show converged density slice + n(z); single-point reps
  show step-0 energy + pipeline energetics. Figures beside each nb in `*_figs/`.

### Builder change (additive, shared skill-local)
`/home/raid/skcb2/skcb2/tddft/.claude/skills/run-notebook/run_notebook_builder.py`:
- New helpers `load_static_density` / `gs_density_fig` / `gs_profile_fig` /
  `single_point_energy_md`; a `static={}` branch (runs only when `ser is None`) and
  two body sections (static density; "Result — single-point energy decomposition").
- Name heuristic: `name = rd.parent.name if rd.name=="results"` (GS results dir).
- **Additive only** — trajectory runs unaffected; also benefits jellium GS baselines.
- Verified on both a GS run (`gs_lz120`) and a single-point run (`wp_r28_p3`,
  E_tot(0)−E_GS = +86.19 eV).

### Docs / memory updated
- `CONTEXT.md` — glossary bullet: `campaign_autorun` is a run-SET, not a "campaign".
- Memory `feedback_verification_user_owns_verdict` — extended with the results-review
  notebook rule (Question/What/Results + quarantined provisional box).

## 2026-07-07 — energy-component streaming + H0 periodicity-2 re-run (DONE)

Follow-up requested by the user: stream **all** INQ energy components each step,
re-run the H0 insertion-energy experiment at **periodicity 2**, update the notebook.
Plan: `/local/data/public/skcb2/tddft/docs/plans/h0-energy-components-periodicity2.md`.

- **inqkit (additive, header-only; no `inq/`/`inq-study/` edit):**
  `inq-stack/include/inqkit/real_time/step_context.hpp`,
  `.../real_time/real_time_session.hpp`,
  `.../io/observables_writer.hpp` now stream `energy_{external,nonlocal,ion,ion_kinetic,exact_exchange,nvxc,eigenvalues}`
  from INQ's public `energy` accessors. New `ObservableSelection` flags default
  **false** → existing runs' CSV schema unchanged.
- **Run drivers:** `scripts/campaign_autorun/{classical,wp}/run.cpp` enable all new
  `sel.energy_*`. Both binaries rebuilt via `inq-run` (INQ_SOURCE=inq-study, GPU 1);
  compiled clean, smoke runs passed.
- **Re-run:** `scripts/campaign_autorun/rerun_h0_p2.py` (new dispatcher) → 12
  single-point runs in `runs/h0_p2/{wp,cl}_r{4,12,20,28,36,40}_p2`, periodicity 2,
  off the open-z GS `runs/h2/gs_p2_lz120/checkpoint`. All completed.
- **Validated:** per-row `total == Σ(kinetic+external+non_local+hartree+xc+exact_exchange+ion+ion_kinetic)`
  to ~1e-13 Ha (machine precision) for all 12 runs. Recorded in
  `docs/validation/test-catalogue.md` (new "inqkit full energy-component streaming").
- **Data (neutral):** WP p2 excess 81.2→79.6 eV vs r (flat, near 81.6 eV zero-point);
  classical p2 excess 185→12 eV (decays). `E_ext` now MEASURED — classical
  r-dependence sits in `E_ext` (153.6→147.2 Ha); WP shows `E_ext`↑ vs `U_H`↓.
- **Notebooks updated + executed (0 errors):** `H0_base_difference.ipynb` gained a
  "periodicity-2 measured decomposition" section (builder `build_notebooks.py`,
  `p2full_md`/`p2full_code`); new `runs/H0_p2_runs.ipynb` evidence table (12 rows,
  `E_ext_Ha`/`sum_minus_total` cols) via `build_run_evidence.py` (`H0_p2` key);
  `00_index.ipynb` re-executed. Kept neutral — interpretation quarantined.
- **Not committed** (user has not requested). Two-commit split when they do:
  inqkit change + run drivers/notebooks (production research) vs any `.claude`.

## 2026-07-07 (later) — interpretation aids + extended-r sweep (DONE)

User began interpreting the p2 decomposition; asked for (a) an individual-run bar
chart summing to total, (b) GS electron-vs-background charge-distribution plots,
(c) larger-r runs to test whether the classical excess reaches 0.

- **New notebook `hypotheses/campaign_autorun_study/H0_p2_interpretation.ipynb`**
  (builder `build_h0_p2_interpretation.py`), executed 0 errors, 4 figures:
  1. **Waterfall** energy decomposition (logical order T, E_ext, E_nl, U_H, E_xc,
     E_ion, E_ion,kin, E_xx) for wp_r4_p2 + cl_r4_p2; running sum == E_total (1e-13).
     Only T/E_ext/U_H/E_xc non-zero (LDA jellium; E_nl=E_ion=E_xx=0).
  2. **GS charge distribution** n_-, n_+, n_- - n_+ in xz/yz/xy + planar n(z) profile
     (both integrate to 82 e; electron spill-out + Friedel = surface dipole).
  3. **Extended-r** excess vs r.
- **Physics verified for the answer:** classical Δ has T/U_H/E_xc frozen, dE_ext ≡
  dTotal (E_ext overlaps total). WP: dU_H=+199.5 eV, dE_ext=-197.9 eV cancel to
  +1.65 eV. ABSOLUTE component signs are reference-dominated (U_H=-83, E_ext=+150
  Ha) — background enters as v_bg=-poisson(n₊) with Hartree on neutralised n_e-n₊
  (`background_perturbation.hpp`); read physics from Δ, not absolute sign.
- **Extended-r sweep (`rerun_h0_p2_far.py`, Lz=200 box, own p2 GS E_GS=60.25 Ha):**
  r={4..76}, wp+cl, all sum-checked 1e-13. **Classical excess: 185→21(r36)→4.1(r44)
  →0.7 eV(r52 minimum)→rises 1.0/1.4/1.8 (r60/68/76).** So it bottoms near r≈52,
  NOT r=40 (where it was still 12 eV) — confirmed the user's skepticism. WP flat
  75.5–77.4 eV. Data presented neutrally; the r≈52 minimum + slight rise is shown,
  not interpreted (user owns verdict; slight rise past r=52 is a candidate follow-up
  — possible far-edge/box effect, unexamined).

## 2026-07-07 (later still) — electrostatic sheet/slab model notebook (DONE)

User supplied their own theoretical modelling (infinite sheet / slab electrostatics,
SI) and asked for a notebook of the *expected* behaviour to compare with the KS runs.
- **`hypotheses/campaign_autorun_study/theoretical_slab_model.ipynb`** (builder
  `build_theoretical_model.py`), executed 0 errors, 4 figures. Implements verbatim:
  single sheet `φ=-σ/(2ε₀)|z_q|`, slab `φ=-ρ₀L/(2ε₀)z_q` (+ exact parabolic interior),
  collapsed sheet. SI constants; slab = same params as runs (n0=1.312e-3, L=25 Bohr).
- **Key model output:** uniform field ⇒ U LINEAR in distance, slope **5.61 eV/Bohr**,
  unbounded (does NOT decay) — the "expected" curve. Slab vs collapsed sheet differ
  only by a constant (φ offset 35 V ⇒ −35 eV energy), same slope outside.
- **Extension (labelled):** two coincident opposite sheets (neutral slab) cancel
  exactly (bridge to the neutral-slab behaviour); assumptions table lists model-vs-KS
  differences (infinite/rigid/single-sheet/point vs periodic/screened/neutral/finite).
  Presented neutrally — no verdict on the model-vs-sim discrepancy (user owns it).

## 2026-07-07 (sign flip + empirical model + UPF-cutoff finding)

- **Sign flip (p2 vs p3) explained + code-grounded:** `poisson.hpp` — 3D kernel G=0→0
  (zero cell-average), 2D kernel G=0→`0.5·rc²` (large, ∝box²). Net-charged n_e picks up
  this offset Φ; `E_hartree=½∫n_e·poisson(n_e)` gets +½ΦN, `E_external=∫n_e·(−poisson(n₊))`
  gets −ΦN → opposite sign, external shift = −2× Hartree shift (data: ΔHartree≈−188,
  Δexternal≈+351≈−2×). Absolute signs are a reference artifact; Δ and neutral total are
  invariant. One-line note added to `docs/campaigns/localised_jellium/ground_state_parameter_study.md`.
- **KEY FINDING — classical decay is (largely) a UPF radial-cutoff artifact.** The ghost
  UPF `electron_gaussian_wpsigma0p5.upf` has z_valence=0 and radial grid **r_max=50 Bohr**
  (1/r tail in Rydberg = +2/r, ends at 50). A point charge above an infinite sheet with
  Coulomb truncated at r_cut has kernel `(r_cut−|dz|)`, decaying to 0 as the slab passes
  beyond r_cut. An unfitted model (r_cut=50) reproduces the classical excess: 209 eV at
  r=4 (sim 185), →0 by r≈40–50 (sim bottoms 0.7 eV at r≈52). So the classical decay is
  dominated by the projectile potential's finite range, NOT physical screening. REVISES
  the earlier "image/screening" reading. To confirm: regenerate UPF with larger r_max, re-run.
- **`theoretical_slab_model.ipynb`** (builder `build_theoretical_model.py`) now has 7 figs:
  single-sheet, slab-vs-collapsed, expectation, neutral superposition, **empirical-density
  plate model** (real n(z) sheet stack → NET residual ±2 eV surface dipole), **cutoff test**
  (above), **image-potential** `−q²/(4z)` (physical expectation once cutoff removed).

## 2026-07-07 (cutoff sweep — DECISIVE) + E_ion answer

- **4-UPF cutoff sweep (CPU) — the classical decay IS the projectile cutoff.** Generated
  4 ghost UPFs truncated at r_cut={10,20,30,40} Bohr (`cutoff_test/make_cutoff_upfs.py`);
  parametrised `classical/run.cpp` with `LJ_PROJ_UPF`; built a CPU binary (`inq-run --cpu`,
  ~6.6 min/run); ran 24 classical single-points (`cutoff_test/run_cutoff_sweep.py`, 8
  concurrent × 6 threads, N_STEPS=1, periodicity 2, Lz=120 GS). RESULT — each ΔE_total(r)
  curve dies at its OWN cutoff:
  ```
  r_cut  dE at r=2,8,16,24,32,40 (eV)
   10    6.7  0.6  0.0  0.0  0.0  0.1   dead ~r=10
   20   35.4 15.6  1.7  0.0  0.1  0.3   dead ~r=20
   30   85.4 53.4 21.4  4.0  0.5  0.7   dead ~r=30
   40  140.7 107.7 64.0 28.6 7.8 1.3    dead ~r=40
  ```
  Conclusive: the classical projectile–slab interaction range is set entirely by the UPF
  radial cutoff (an appreciable, controlling effect), NOT physics. Added as the 9th figure
  in `theoretical_slab_model.ipynb` (auto-built on sweep completion).
- **E_ion answer:** `E_ion` and `E_ion_kinetic` ARE part of `total()` and ARE streamed —
  identically 0 for the ghost (z_valence=0; jellium background is external, not Ewald). The
  ghost–background interaction the user expects is the separately-omitted `∫v_ghost·n₊`.
- **Key disentangling (deliberation):** WP has NO radial cutoff (full `poisson` solve), so
  WP flatness = neutral-slab cancellation + screening (physics); classical decay = cutoff
  artifact. E_xc(WP)−E_xc(classical) isolates the WP–slab Pauli/correlation term. Next
  proposed: (1) full-component WP-vs-classical overlay + E_xc diff (no runs); (2) save n_WP,
  poisson(n_WP) vs v_ghost + induced density; (3) corrected classical (+ghost-bg) vs model.

## 2026-07-09 — Periodicity-3 full-component mirror + Part III of the model notebook (DONE, 0 errors)

Scoped via `/grill-with-docs`. User wanted, for ONE run (r=28): (1) the WP-vs-classical total
magnitude difference, (2) each run's individual component decomposition, (3) the component-wise
difference, (4) an analytic check that the "missing energy" (E_ion-ion Ewald + projectile-ion
Coulombic) matches. Ground-up, no assumptions.

### Grill outcomes (physics corrected before building)
- **E_ion (Ewald) is NOT missing/unstored.** `classical/run.cpp:69` streams `energy_ion` +
  `energy_ion_kinetic`; both are **identically 0** (ghost z_valence=0; jellium background is an
  *external* potential, not Ewald). The classical run's only omitted piece is `∫v_ghost·n₊`
  (`ghost_background_term_omitted=true`, re-added in H5). Retargeted the "missing energy" to the
  **WP self-Hartree** (Gaussian Coulomb self-energy).
- **Individual E_hartree/E_external are NOT readable at periodicity 2.** `poisson.hpp:49` sets the
  p2 G=0 term to `0.5·rc²` (rc=Lz=120) → each of Hartree/external carries a large opposite-sign
  offset on a net-charged cell; only their SUM is physical. Confirmed with data (cl E_ext=+154 Ha,
  E_H=−81 Ha for a 60 Ha system). Fix = periodicity 3 (`poisson.hpp:31`, G=0→0).

### New runs — `runs/h0_p3/{wp,cl}_r{4,12,20,28,36,40}_p3` (12 single-point, DONE)
- Dispatcher **`scripts/campaign_autorun/rerun_h0_p3.py`** — mirror of `rerun_h0_p2.py` at
  **periodicity 3**, p3 GS `shared_gs/slab_n82_L50x50x120`, ALL energy components streamed.
  Backends: **wp→GPU** (CUDA binary), **cl→CPU** (binaries were inconsistent; physics is
  backend-identical, verified per-row by sum(parts)==total to ~1e-13). Idempotent.
- **Key p3 result (physical signs, offset-free):** E_hartree>0 (e–e repulsion), E_external<0
  (electrons in +background well) — OPPOSITE sign to p2 (the sign-flip finding, now clean).
  Classical kinetic/hartree/xc are FROZEN GS values (ghost = pure external ⇒ step-0 density is the
  GS density); ALL classical r-dependence is in E_external.

### Notebook Part III (`build_theoretical_model.py` → 18 code cells, 0 errors, 17 figs)
4 panels (r=28 single run + r-sweep): **P1** total-magnitude ΔE=+32 eV; **P2** per-run waterfall
(each sums to total 1e-13); **P3** component difference vs r; **P4** self-energy ledger + caveat.

### Findings (recorded; user owns verdicts)
- **Zero-point KE:** d(kinetic) = +3.004 Ha (81.7 eV) at every r = exact `3/(4σ²)` = 3.000 Ha (0.1%).
- **Self-XC:** d(xc) = −0.605 Ha (−16.5 eV), r-independent (WP LDA self-XC).
- **Self-Hartree:** closed form `1/(2σ_ρ√π)` = 0.798 Ha (21.7 eV) ↔ numeric FFT-Poisson on saved
  `n_WP` = 0.774 Ha (21.1 eV), agree ~3%.
- **KEY CAVEAT (charged-cell convention):** raw `E_hartree(WP)−E_hartree(CL)` = −29 eV (p3) vs
  −274 eV (p2) vs physical +22 eV — matches NEITHER. Inserting the WP makes the cell net −1
  charged, so the Poisson G=0 convention injects `~N_slab·mean(V[n_WP])` into the Hartree/external
  split. The self-Hartree is recovered from `n_WP` directly, NOT from a component subtraction.
  See [[reference_charged_cell_hartree_convention]].

### Open / not done
- Optional (declined by default, offered): extend p3 sweep past the ghost UPF cutoff (~r=55) to
  show `dtotal(r)` plateau to the +3.2 Ha quantum floor.
- Not committed to git (user has not requested).

## 2026-07-08 — Screening / WP-potential test + classical-vs-WP deconstruction (DONE, 0 errors)

User asked (notes l.109–112) to run the two "Potential Learning" tests autonomously on
CPU and extend `theoretical_slab_model.ipynb`. Plan:
`/local/data/public/skcb2/tddft/docs/plans/screening-wp-potential-test.md`.

### Code (run machinery)
- **`scripts/campaign_autorun/wp/run.cpp`** — added env-gated `LJ_SAVE_DENSITY=1` block that,
  at t=0 (after WP injection, before propagate), writes VTIs `density_wp` (=|ψ_WP|²),
  `density_total`, `density_bath` via `RealField3DWriter`. Rebuilt as a **CPU** binary
  (`inq-run --cpu`, `ENABLE_CUDA=OFF`, INQ_SOURCE=inq-study); compiles + runs (exit 0).
- **`scripts/campaign_autorun/screening_wp_test.py`** — driver; ran WP insertion at r=12 (clean)
  and r=4 (near surface), p2, 1 step, CPU. Output `runs/screening_wp/wp_r{4,12}_p2/` — 2/2 done,
  6 VTIs written, `∫n_WP = 1.0000` both.

### Notebook (`build_theoretical_model.py` → 13 figs, 0 errors)
Added **Part II** (4 sections): LEDGER (full deconstruction, sum(parts)=total to 4.6e-13),
XCDIFF (E_xc diff + total-diff split), WPPOT (★ the screening/WP-potential test), BATH
(screening baseline). Neutral; σ=σ_WP labels; figures inline (matches Part I).

### Findings (recorded; user owns verdicts)
- **Deconstruction exact.** WP−CL = r-**independent** quantum self-energy (dKin=+81.6 eV =
  zero-point 3/(4σ²); dXC=−16.5 eV = WP self-XC) + r-**dependent** electrostatic d(Hartree+external)
  −169 eV (r=4) → +2 eV (r=40). WP excess FLAT ~80 eV; classical decays 185→12 eV (Part-I cutoff term).
- **E_xc(WP)−E_xc(CL) = −16.5 eV, r-independent** → WP self-XC, not a slab-proximity screening signature.
- **WP-potential test (Learning #2):** n_WP tracks the ideal Gaussian charge (σ_ρ=σ_WP/√2=0.354) at
  the core; `poisson(n_WP)` (Python FFT-Poisson, validated vs analytic erf to RMS 3e-4 Ha) overlays the
  analytic ghost potential `erf(r/(√2σ_ρ))/r`. Same source ⇒ any energy gap is quantum, not a different
  potential. Tail deviations = orthogonalisation vs occupied slab states.
- **Screening baseline = 0 EXACTLY.** `n_slab(t=0) − n_GS` bit-identical (max|Δ|=0). Genuine screening
  is dynamical (~T_plasmon ≈ 4893 steps) — flagged as GPU follow-up, not faked.
- **GOTCHA:** `inqkit::fields::density::total()` returned **82 e (slab only, WP excluded)** in this
  config — the exact ambiguity flagged in the `density.hpp` TODO. WP captured via `density::orbital`
  (=1 e). Screening baseline therefore uses `density_total − gs` (NOT the saved `density_bath`, which
  subtracts a WP that `total` never contained). See [[reference_inqkit_density_total_excludes_wp]].

### Open / follow-up (GPU)
- Static screening: re-converge SCF with a fixed external Gaussian ghost (GS run) → induced density.
- Dynamical screening: ~T_plasmon of propagation (~5000 steps) — GPU only.
- Not committed to git (user has not requested).

## Constraints honoured
- No assistant interpretation; questions taken from the user's own `hyp` wording;
  existing verdicts kept but quarantined in marked provisional boxes.
- No new INQ runs / no GPU — reads existing results only. venv python throughout.
- Figures travel beside notebooks; VTIs read via `inqview.load_vti` (no fftshift).

## Not done / open (for the user, not the assistant)
- Deriving learnings + further experiments — the user's job (by design).
- **Stale aggregators left on disk, OFF the reading path** (not deleted, not
  regenerated): `campaign_autorun_study.ipynb` and `campaign_summary.ipynb`. Superseded
  by `00_index.ipynb`. Decide whether to delete.
- Nothing committed to git yet (user has not requested a commit).
- Cosmetic: `rep_H3_gs_a15_N98` density fig shades the slab at ±12.5 (default) though
  a=15; guide-only, not a data error.
