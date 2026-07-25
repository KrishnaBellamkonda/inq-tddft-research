# Handover: Absorbing-boundary benchmarking & parameter search

Plan: `/local/data/public/skcb2/tddft/docs/plans/absorbing-boundary.md` (authoritative, all decisions LOCKED).
Glossary: "Absorbing boundaries" section of `/local/data/public/skcb2/tddft/CONTEXT.md`.
Prompt: `/local/data/public/skcb2/tddft/docs/prompts/absorbing_boundary/benchmarks_and_parameter_search.md`.

---

## Milestone: 2026-06-17 — two-sided notebook figures revised (linear companions + 1% line)

Edited the builder `ResearchProject/systems/vacuum/hypotheses/twosided_cap_vs_mask/build_twosided_report.py`
and re-ran it (venv python, `PYTHONPATH=inq-stack/python`) — notebook rebuilt,
**0 error cells**, 105 runs.

- The three **reflectivity (ε) curves** (§1 ε(E)-by-η, §2&3 ε(E)-by-L CAP|mask,
  §4 ε-vs-L) now each show a **log panel beside a linear panel (y 0–1)**;
  §2&3 became a 2×2 (top log / bottom linear). ε now plotted as a **fraction**
  (was %). Every panel carries a **dashed 1 % line (ε=0.01)** labelled "1%".
- §6 t_absorb (time, not reflectivity) left unchanged.
- **§5 ε(E,L) heatmap** now a 2×2: top row log₁₀ ε, bottom row **linear ε**
  (CAP η=−0.5 │ mask), one shared colorbar per row.
- **New §8 — density carpets (z–t)**: 2D maps `x=z (Bohr), y=time (a.u.),
  colour=⟂-summed WP density` for CAP η=−0.5, L=30 Bohr at **2 / 10 / 100 eV**
  (lowest / 10 eV / closest-to-100). Absorber edges marked (dashed inner edge,
  dotted box wall). Reuses §7's VTK `zprof`; time axis exact via step·dt
  (dt=0.01). Output `fig_density_carpet_L30.png`.

Second 2026-06-17 pass (3 more requests):
- **§5 heatmap colormap REVERSED** (`viridis`→`viridis_r`): bright/yellow now =
  low reflectivity (comfortable), per user "we want lower reflectivity".
- **New §1b — ε vs η at L=20** (one curve per energy): the orthogonal cut of §1
  (η on x-axis). DATA NOTE: the η sweep exists **only at L=20** (η∈{−0.3,−0.5,
  −0.7,−1.0}); energies with all 4 η = 2/10/32/100/300 eV. log │ linear(0–1) +
  1% line. Output `fig_eps_vs_eta_L20.png`. (§1 already had the ε(E)-by-η view.)
- **New §9 — xz-plane density GIFs**: animated mid-y slice (x horizontal, z
  vertical, colour=density) for CAP η=−0.5, L=30 Bohr at 2/10/100 eV. Frames
  strided to ≤120. Outputs `fig_xz_density_E{2,10,100}.gif`. Complements the §8
  static carpets. Notebook timeout bumped 1800→3600 s for the extra GIFs.
- **Anchor energy = 10 eV confirmed** (`ANCHOR_E=10`, `ETA_STAR=−0.50`,
  `ANCHOR_L=20`). Abstract prose ("anchor energy"/"anchor width") replaced with
  the concrete `10 eV` / `L=20 Bohr` / `η=−0.5` in titles + markdown.
- Codified the principle in the **report-figures skill** (new production rule 6):
  prefer the concrete value over an abstract label; reflectivity curves get a
  linear 0–1 companion + explicit threshold line.

---

## Milestone: 2026-06-16 (eve) — two-sided CAP vs mask: implemented, validated, sweep LAUNCHED

Plan: `docs/plans/twosided-cap-vs-mask.md` (DESIGN LOCKED via grill). Successor
prompt (σ=0.5 baselines): `docs/prompts/absorbing_boundary/sigma0p5_baselines_with_locked_params.md`.

**Implemented (inq/ AND inq-study UNTOUCHED — all new code is inq-stack + scripts):**
- Two-sided **mask**: `inq-stack/include/inqkit/absorbers/mask_absorber.hpp`
  (`TwoSidedMaskAbsorber`, `inner_region_norm_twosided`) + INQ-free shape header
  `mask_shape.hpp` (`sin2_mask_value_twosided`).
- Two-sided **CAP**: NO new engine code — `run.cpp` composes
  `absorbing(η,+mid,w) + absorbing(η,−mid,w)` via existing `perturbations::sum`.
- Run binary: `scripts/twosided_cap_vs_mask/run.cpp` (env mode=cap|mask), built
  (101 MB) against inq-study. Dispatcher: `scripts/twosided_cap_vs_mask/dispatch.py`.
- Notebook builder: `hypotheses/twosided_cap_vs_mask/build_twosided_report.py`
  (notebook-making skill; 7 results + 2 GIFs).

**Validation (all PASS):**
- Pure mask-shape test `inq-stack/tests/include/inqkit/absorbers/test_mask_shape.cpp`
  (symmetry, inner=1, walls=0, mid=½) — ctest PASS.
- Engine mechanism check `hypotheses/twosided_cap_vs_mask/tests/mechanism_check.py`:
  **cap absorbed 99.2% (ε=0.78%), mask absorbed 76.5% (ε=23.3%)** at anchor
  (E≈10, L=20) — both PASS (mask is the softer absorber by design; CAP stronger at
  equal width — the study's central comparison, already visible).
- Catalogued in `docs/validation/test-catalogue.md`.

**Geometry (LOCKED, energy-scaled, quasi-monochromatic):** σ=4√2/k0; two-sided,
L split L/2 per end; z_in=6σ, Lcell=12σ+L; WP launched z0=−z_in+4σ (≥4σ from near
absorber, negligible CAP density), moving +z; ε=∫_{|z|<z_in}|ψ|²/N0; τ=2(z_in+L)/k0.

**Sweep RUNNING (launched 2026-06-16 ~17:05, dispatcher PID 1565208, nohup, 2 GPUs):**
`scripts/twosided_cap_vs_mask/dispatch.out`. 3 phases, EMAILS user per phase, then
auto-builds the notebook + final email.
- Grids: E={2,4,10,16,32,64,100,300,1000} (1 eV DROPPED — σ-scaled box too costly);
  L_total={10,16,20,26,30}; η={−0.3,−0.5,−0.7,−1.0}; anchor E=10, L=20, η*=−0.5.
- Phase 1 CAP η-sweep @ L=20 (E_ETA×η=20). Phase 2 CAP L-sweep @ η* (E×L=45).
  Phase 3 mask L-sweep (E×L=45). ~110 runs; resumable (skips existing epsilon.txt).
- Run dirs `twosided_cap_vs_mask/run_{cap,mask}_E{E}_L{L}[_eta{|η|}]/`. Smoke dirs
  `scripts/.../smoke_{cap,mask}/` feed the mechanism check (keep).

**⚠ FLAGGED to user (pre-existing, NOT mine, NOT touched):** `inq/src/real_time/
viewables.hpp` carries a project-local `ham()` edit ("jellium KS-energy observables,
2026-05-01") — a violation of the inq-immutable rule predating this task. Does NOT
affect this sweep (builds vs inq-study; vacuum runs don't use it). Governance item
for the user; per the rule I did not edit/revert it.

**Open:** sweep completing overnight (phase emails are the live updates); the user
reads the **comfortable region** (L*, η*) off the notebook §4–5 (not auto-declared).
All CAP ε PROVISIONAL until Task #7. After lock, run the σ=0.5 baseline successor.

---

## Milestone: 2026-06-16 — notebook-making skill + cap_Lopt_E10 notebook (gap filled)

**Done.**
- Built the missing **`cap_Lopt_E10` study notebook** (the widen-L sin² baseline
  had results but no notebook): `ResearchProject/systems/vacuum/hypotheses/cap_Lopt_E10/cap_Lopt_E10_study.ipynb`
  (+ `build_Lopt_report.py`, `cap_Lopt_combined.csv`, `fig_Lopt_eps_vs_L.png`,
  `fig_Lopt_eps_vs_eta.png`). Executed, 0 errors. **Result (E=10 eV, sin², ETRS,
  PROVISIONAL):** depth optimum **η=−0.5** (beats −0.3 and −1.0 at every width);
  ε = 2.74% (L8) → 1.13% (L10) → 0.49% (L12) → 0.15% (L15). η=−0.3 L6 = 14.3%
  (worst). Reaching ε≈1% costs ≥ L≈10 Bohr.
- The four current notebooks (decision-support reminder):
  `hypotheses/cap_real/cap_real_study.ipynb`,
  `hypotheses/cap_thin_L5/cap_thin_L5_study.ipynb` (current),
  `hypotheses/cap_monomial/cap_monomial_study.ipynb` (current),
  `hypotheses/cap_Lopt_E10/cap_Lopt_E10_study.ipynb` (new). All live.
- **New skill `notebook-making`** (`.claude/skills/notebook-making/SKILL.md`):
  house narrative (context → formulas-with-every-term-defined → fully
  reconstructable setup → linked source files → results → takeaway), grounding/
  PROVISIONAL rules, and the auto-build convention.
- **Auto-build design (user-locked):** per-**sweep** notebook · trigger in the
  **run machinery** (dispatcher/`analyse.py` tail), NOT a Claude-Code hook ·
  rebuilt **once at end of batch**. PLUS a *consequential individual run* (one with
  a `docs/plans/` entry) auto-produces its own single-run notebook via its
  `analyse.py`. Worked wiring added to `scripts/cap_Lopt_E10/dispatch.py`
  (`autobuild_notebook()` tail).
- CONTEXT.md glossary: added **"Study notebook"** term.

**Open / not done.** Existing dispatchers (`cap_thin_L5`, `cap_monomial`,
`mfa_sweep`, jellium) not retrofitted with the tail call — convention applies
going forward; retrofit only if re-run. Task #7 still pending → all CAP ε
PROVISIONAL.

---

## Milestone: 2026-06-15 (eve) — Monomial CAP (inq-study) + widen-L baseline

Plan: `/local/data/public/skcb2/tddft/docs/plans/cap-monomial-inq-study.md`.
Driver: thin L=5 sin² CAP floors at ε≈0.20 (10 eV); user wants ~0. Confirmed
tree-wide that INQ's built-in CAP is **sin²-only, no monomial/shape knob**
(`inq/src/perturbations/absorbing.hpp`, the ONLY CAP class) → a monomial must be
a NEW perturbation in inq-study.

**Part A — built-in widen-L baseline [RUNNING, bg bqm7sfsg9]:** `scripts/cap_Lopt_E10/`
2D sweep L∈{6,8,10,12,15}×η∈{−0.3,−0.5,−1.0} at E=10 eV, built-in sin² CAP (reuses
cap_sweep binary). Finds smallest L reaching target ε. ~13/15 done at last check.

**Part B — `absorbing_monomial` in inq-study [DONE: header+test+binary src]:**
- `inq-study/src/perturbations/absorbing_monomial.hpp` (NEW; inq/ untouched) —
  ramp `V=i·eta·s^n`, s∈[0,1] across slab (0 inner edge → max at wall, unlike the
  sin² hump). Ctor `(amplitude, mid_pos, width, order)`. Cites De Giovannini 2014
  §IV / Riss-Meyer 1996. In-header Catch2 test (construct + has_potential).
- Known-case test `hypotheses/cap_monomial/tests/monomial_shape_check/run.cpp`:
  checks (1) it absorbs, (2) ε(n=1)<ε(n=4) — the order-monotonicity signature
  unique to s^n (sin² can't have it). **BUILDING vs inq-study [bg b8gapatl8].**
- `scripts/cap_monomial/run.cpp` = cap_sweep/run.cpp + 3 edits (include header,
  read CAP_ORDER, use absorbing_monomial); full obs set + manifest. Not yet built.
- `scripts/cap_monomial/dispatch.py` — benchmark n∈{1,2,3,4}×η∈{−0.1..−0.5} at
  L=5,E=10 (16 runs) vs sin² baseline. Ready, gated on the test passing.

**RESULTS [DONE 2026-06-16]:**
- Known-case test PASSED (absorbs; ε(n=1)=0.223 < ε(n=4)=0.548). Catalogued.
  inq/ verified pristine (diff shows only absorbing_monomial.hpp new).
- Built cap_monomial binary; ran 16-run benchmark (n∈{1,2,3,4}×η∈{−0.1..−0.5},
  L=5, E=10 eV). **Monomial n=1 (linear ramp) beats sin² 2.5×: ε=8.3% (η=−0.5)
  vs sin² L=5 = 20.9%.** Lower n + deeper η monotonically better (ramp peaks at
  the wall where sin² hump is zero). BUT 8.3% ≈ sin² at only ~L6 — buys ~1 Bohr
  effective width, not ε→0.
- Built `hypotheses/cap_monomial/cap_monomial_study.ipynb` (0 errors) + figs;
  emailed comparison to user ([cap-thin-L5] thread). ε PROVISIONAL until Task #7.
- **Part A baseline (cap_Lopt_E10, 15 runs):** sin² needs L≈8→2.7%, L10→1.1%,
  L12→0.49% at 10 eV; η=−0.5 sweet spot. (No analysis notebook built yet.)

**OPEN LEVERS toward ε→0 at L=5 (both untested):**
1. **Deeper η at n=1** — benchmark stopped at η=−0.5, trend still improving;
   reuse cap_monomial binary with η∈{−0.7,−1.0,−1.5,−2.0} at n=1, L=5, E=10
   (~quick). Watch for the over-absorption U-turn (steep ramp → wall reflection).
2. **Transmission-free CAP** (Manolopoulos, JCP 117 9552, 2002) — new inq-study
   perturbation `absorbing_tf`, same scaffolding as absorbing_monomial; the form
   designed for ε→0 at SHORT L. The principled route if thinness is hard.

---

## Milestone: 2026-06-15 — Thin-CAP (L=5) reflectivity tuning + vacuum reorg

Plan: `/local/data/public/skcb2/tddft/docs/plans/cap-thin-absorber-tuning.md`.
Grilled 2026-06-15 (`/grill-with-docs`). Goal: find the in-built CAP params giving
low ε across 1–100 eV with the **minimum near 10 eV**, under a **thin L=5 Bohr**
absorber and **shallow** η ∈ {−0.01,−0.05,−0.30} Ha.

**Folder reorg (DONE) — amends ADR 0007 (flat `run_*` → grouped-by-sweep).**
Runs now group as `systems/vacuum/<sweep_name>/<run_name>/`; analysis in
`hypotheses/<sweep_name>/` (bare names, NN_ dropped). Migrated on disk:
- `run_cap_*` (17) → `vacuum/cap_real/run_cap_*`
- `runs/run_mfa_*` (72) → `vacuum/mfa_sweep/run_mfa_*`; `epsilon_grid.csv` → `hypotheses/mfa_sweep/`
- top-level `mfa_sweep/` machinery → `vacuum/scripts/mfa_sweep/`
- `hypotheses/01_cap_real/` → `hypotheses/cap_real/`; fixed `build_cap_report.py`
  glob (`cap_real/run_cap_*`) + re-ran → 17 runs, 0 errors (notebook intact).
Ecosystem reconciled: ADR 0007 Amendment (2026-06-15), `CONTEXT.md` folder section,
`.claude/rules/file-placement.md`. (`vacuum/tests/` left in place — note for later.)

**Study design (33 CAP runs = 3 η × 11 E).** L=5 fixed; E = {1.01,1.87,3.48,6.46,
7,10,15,22.24,41.28,76.62,100} eV (MFA-comparable ladder + densified at 10 eV);
ETRS; full free-WP minimum observable set per run; all runs write density_wp VTI
(grids tiny: 8×8×Nz, ~250 MB total). **No rebuild** — reuses the env-parameterised
`scripts/cap_sweep/run` binary (CAP_L=5).

**Machinery (DONE):**
- `scripts/cap_thin_L5/dispatch.py` — new menu, points BINARY at `scripts/cap_sweep/run`,
  run dirs `cap_thin_L5/run_cap_k{k0:.2f}_L5_eta{|η|:.2f}`, GPUs 0,1, optional
  substring filter for smoke.
- `hypotheses/cap_thin_L5/build_thin_report.py` — combined CSV + overlaid ε(E)
  curves (3 η, log-log, 10 eV guide, min-markers) + absorbed-fraction + dynamics +
  density GIF + executed `cap_thin_L5_study.ipynb`. READY (not yet run).

**Smoke (DONE):** 10 eV/η=−0.05 → 52 s (5784 steps, ~9 ms/step). ε=0.729,
absorbed=0.255 → thin+shallow UNDER-absorbs (expected high floor); η=−0.30 will do
better. Manifest validates: tiers 1–3 ok, tier-4 energy-drift + norm-band fail by
design.

**STATUS (in flight):** 33-run sweep launched (detached PID; `scripts/cap_thin_L5/sweep.log`).
Waiter bg-task `b93qk4u89` notifies at completion. Firm ETA ~60–70 min compute on
2 GPUs (1.01 eV run dominates ~21 min each). On completion: run `build_thin_report.py`,
then **email ε(E) curves + density GIF to chiddukanna@gmail.com** via
`inqview.email.send_run_email` (subject family `[cap-thin-L5]`). ε **PROVISIONAL**
until Task #7.

---

## Milestone: 2026-06-14 (autonomous) — In-built CAP investigation, FULL observable set

User redirect: investigate the **in-built** `perturbations::absorbing` thoroughly,
run with the **full free-WP minimum observable set** (ADR 0006), document in ipynb;
DEFER inq-study validation (Task #7) → `docs/handovers/inq-study-cap-deferred.md`.
Running AUTONOMOUSLY (ETA ~1.5 h).

**Skepticism resolved + vindicated (git archaeology):** the absorbing CAP is a real
team feature (Yao+Andrade 2023); a 2024 refactor `bd4a46fe` silently regressed it
(real `vscalar`, no test drives it through propagation). My 1-line inq-study fix is
a REGRESSION REPAIR, not "fixing broken team code". `inq/` confirmed pristine
(only pre-existing reduce.hpp + viewables.hpp mods, neither mine).

**Augmented `scripts/cap_sweep/run.cpp`** now emits the full free-WP set + manifest:
observables.csv (energy_total/kinetic/**hartree/xc**, current, dipole, **density_l2**
via DensityDelta on system density), wp_momentum_stats, wp_real_space_stats,
momentum_distribution, **gs_eigenvalues + occupations** (final-state), **gs_system
density VTI**, density_wp VTI (61 frames), run_summary — all at manifest paths
`results/raw/observables/...`. Smoke validated: tiers 1-3 PASS; tier-4
energy-drift + wp norm-band FAIL **by design** (CAP is non-Hermitian → energy+norm
decrease = the absorption signature; documented, not a defect).

**2D sweep launched** (`dispatch.py`, bq6pxrrq0): depth η{0.01..4.0}@L20 ∪ width
L{5..50}@η0.5 ∪ energy k0{0.86,1.28,2.0,2.71}@η0.5 ≈ 17 flat `run_cap_*`, 2 GPUs.
Report: `hypotheses/01_cap_real/build_cap_report.py` (rewritten for 2D + full obs +
validation summary + density gif). Validator path: `python -m inqview.validation
<run_dir>` reads `<run_dir>/results/observables_manifest.json`. STILL PROVISIONAL
pending Task #7.

**AUTONOMOUS RUN COMPLETE (2026-06-14 ~02:30).** 17/17 runs OK. All validate as
`free-wp` with **tier 1-3 failures = 0** (observable set complete) and exactly the
2 expected tier-4 violations (energy drift + wp norm-band = absorption signature).
Results: depth U-shape (min ε=1.2e-5 @η=1, L20), width monotonic (L5→0.14 …
L50→3e-7), energy band (best ~mid-E). Deliverable executed (8 cells, 0 errors):
`hypotheses/01_cap_real/cap_real_study.ipynb` + figs fig_cap_{depth,width_energy,
dynamics}.png + fig_cap_density.gif. Each run ~19 MB full observables. Task #7
(engine regression) remains the ONLY outstanding gate → results provisional.

---

## Milestone: 2026-06-14 — Folder standard (ADR 0007) + CAP route decided

### Folder structure standardised (ADR 0007)
Grilling session locked the canonical `systems/<name>/` layout, ALIGNED to the
established jellium/coronene convention (NOT the cleaner-sounding alternative):
- `shared_gs/` (new unifying GS name), `shared/`, `scripts/`, FLAT top-level
  `run_<type>_<params>/`, plural `hypotheses/<NN_purpose>/` (numbered, with a
  `tests/` subfolder for task-specific checks).
- Two-tier tests: library-generic `inqkit` → `inq-stack/tests/`; task-specific →
  `hypotheses/<NN>/tests/`.
- Existing systems **grandfathered** (not migrated); vacuum is the reference
  instance (rearrange still PENDING — see below).
- Encoded in: `docs/adr/0007-system-folder-structure.md`, `.claude/rules/file-placement.md`,
  `CONTEXT.md`, skills (`tddft-simulations`, `build-run`, `literature-review`,
  `report-writing`), memories (`reference_system_folder_structure`,
  `feedback_ipynb_run_reports_direction`). Agents/hooks needed no change.

### KEY FINDING — INQ's built-in CAP is non-functional (engine route required)
`perturbations::absorbing` (inq/inq-study `src/perturbations/absorbing.hpp`) IS a
region-restricted sin² CAP `V=+i·η·sin²` in a FRACTIONAL z-slab. BUT a probe
(`vacuum/hypotheses/01_cap_real/tests/cap_probe/run.cpp`) FAILED TO COMPILE:
```
absorbing.hpp:45  error: no operator "+=" : double += inq::complex
  self_consistency hands perturbations a field<real_space, DOUBLE> scalar potential
```
→ INQ's local scalar potential is **real-typed**, so an imaginary CAP has nowhere
to live. `absorbing` is dead code (its unit test never calls `potential()`).
**Decision (user): implement a TRUE integrated CAP in `inq-study`** by
**complexifying the scalar potential** so the in-built `perturbations::absorbing`
works (user chose: use the built-in perturbation, complexify, inq-study only).
The callback route (multiply `exp(η·sin²·dt)`) was REJECTED as a re-skinned mask.
Plan: `docs/plans/absorbing-cap-engine.md`.

#### Phase 1 IN PROGRESS (engine change made; build running 2026-06-14)
Discovered the change is SURGICAL, not a multi-file refactor: `scalar_potential_`
is already `field_set<…, PotentialType>` (`ks_hamiltonian.hpp:47`) and `vks` is
`field_set<…, HamiltonianType::potential_type>` — so the KS potential is ALREADY
complex in real-time (`ks_hamiltonian<complex>`). Only `vscalar` was the real
bottleneck. Edits (inq-study ONLY):
- `src/hamiltonian/self_consistency.hpp:176` — `auto vscalar = vion_;` →
  `field<real_space, typename HamiltonianType::potential_type> vscalar(vion_);`
  (complex in RT, double in GS → GS bit-identical). Uses field's converting ctor.
- `src/hamiltonian/self_consistency.hpp:191` — wrap `energy.external(...)` in
  `inq::real(...)` (CAP imaginary part excluded from energetics).
- `shared/config.sh` — made `INQ_SOURCE`/share paths env-overridable (defaults
  unchanged) so `INQ_SOURCE=…/inq-study inq-run` builds the fork. Build with:
  `INQ_SOURCE=…/inq-study INQ_SHARE_PATH=…/inq/install/share PSEUDOPOD_SHARE_PATH=…/inq/install/share/pseudopod inq-run --reconfig`.
**Phase 1 GREEN (2026-06-14 01:02).** cap_probe compiled + ran against inq-study
(was `double += complex` against inq/) and ABSORBED: η=−0.5 Ha, E=13.6 eV, L=20 →
`total_wp_norm(τ)=4e-6, absorbed=0.9999, inner_eps=5e-5`. ETRS stable (no NaN,
"ended normally"); energy stayed real (~67.3 Ha) — `energy.hpp` already `real()`s
band energy/eigenvalues (`:45,:91`), so NO further engine edit needed. First-ever
run of `perturbations::absorbing` in this codebase.
Build recipe (works): `INQ_SOURCE=…/inq-study INQ_SHARE_PATH=…/inq/install/share
PSEUDOPOD_SHARE_PATH=…/inq/install/share/pseudopod inq-run --reconfig`.

**Phase 2/3 staged:** production binary `vacuum/scripts/cap_sweep/run.cpp` (engine
CAP + full min-observable suite + showcase density VTI) + `dispatch.py` (η-depth
menu, 2-GPU, flat `run_cap_*`). η=−0.5 already near-total absorption → menu widened
to span regimes (weak incomplete-absorption → strong CAP-reflection turn-up).
NEXT: build cap_sweep against inq-study, run depth sweep, build report ipynb
(density gif + energetics + ε-vs-depth + data-path links) in `hypotheses/01_cap_real/`.

**Phase 2 DONE (2026-06-14) — depth sweep, textbook CAP U-shape.** 11 runs
`vacuum/run_cap_*` (flat, ADR 0007), full min-observable suite each, showcase
(k1.28,η0.25) has density VTI+raw frames. ε(|η|) at E=22 eV, L=20:
0.01→0.241, 0.02→0.197, 0.05→0.108, 0.10→0.041, 0.25→2.4e-3, 0.50→3.0e-5,
**1.00→1.2e-5 (min)**, 2.00→6.1e-5, 4.00→1.7e-4 — i.e. weak under-absorbs →
sweet spot η≈1 Ha → strong REFLECTS (turn-up), matching De Giovannini Fig.4.
Both GPUs were free (25 GB). Report builder: `hypotheses/01_cap_real/build_cap_report.py`
(vtk reads density_wp/*.vti for the gif). ⚠️ **ALL CAP ε PROVISIONAL until Task #7**
(inq-study ctest regression) confirms the complexify didn't break the engine —
user: "determines if my simulation results are meaningful at all".

**Phase 3 DONE (2026-06-14 01:20).** `hypotheses/01_cap_real/cap_real_study.ipynb`
(executed, 7 code cells, 0 errors) + figures: `fig_cap_density.gif` (61-frame WP
z-profile meeting the CAP), `fig_cap_eps_vs_depth.png` (the U-shape),
`fig_cap_survival_vs_time.png`, `fig_cap_energetics.png`, `fig_cap_eps_vs_energy.png`.
Notebook carries the governing eq + method header + data-path links + observations.
Engine-CAP sub-project (Phases 1–3) COMPLETE; only Task #7 validation outstanding.
OPEN (user decision, not auto): `inq/src/real_time/viewables.hpp` carries a
pre-existing 2026-05-01 `ham()` accessor — a 2nd undocumented in-place inq/ edit
beyond the allowed CUB workaround; recommend migrating it to inq-study or
documenting as a sanctioned exception. Also still deferred: vacuum rearrange,
MFA mask weight/depth benchmark, notebook headers/Q&A on the 3 feasibility ipynbs.
Coordinate convention resolved: `rvector`=contravariant/fractional,
`rvector_cartesian`=Bohr; fractional CAP slab `[0.5−L/Lcell, 0.5]` == MFA absorber
`[(6σ−L)/2,(6σ+L)/2]`, so engine-CAP ε is directly comparable to MFA ε.

---

## Milestone: 2026-06-13 (complete) — Overnight run DONE; both tasks delivered

### Current status: COMPLETE
72/72 ε(E,L) runs OK (0 aborts) in ~3.4 h on 2 GPUs. Both notebooks built+executed.
Summary email sent to chiddukanna@gmail.com (msg-id 178133386721…). Nothing pending.

### Results
- **ε(E,L) surface**: ε ∈ [5.9e-10, 0.806] over L∈{5,10,20,30,40,50} × 12 energies
  E=0.54–489.8 eV. Clean monotonic roll-off (high-reflection plateau ~0.5–0.8 at
  low E → ~1e-9 at high E), matching the paper's Fig. 3 family.
- **Cross-validated** vs independent NumPy 1D split-operator `cap_toy` (no shared
  code): e.g. (E=13.6 eV, L=20) INQ ε=0.086 vs toy ~0.08.
- **6 showcase runs** carry density VTI (total/system/wp) + final WP wavefunction.

### Final config (what produced the data)
ETRS propagator (NOT CN — CN renormalises the WP, undoes absorption); dt=0.01
(paper value — ε is dt-sensitive: 0.086/0.047/0.030 at dt=0.01/0.02/0.03, a
first-order mask-method property, documented in the notebook); adaptive
dx=clamp(0.75/k0, 0.18, 0.30); NPERP=8 (ε transverse-insensitive, validated);
σ pinned to 4√2/k0; ε = inner-region survival / N0 (cancels transverse truncation).
Dispatcher: cheapest-first, incremental CSV, cudaMemGetInfo GPU probe.

### Deliverables (all on disk)
- `docs/reports/absorbing-boundary/mfa_reflectivity_study.ipynb` (Task 2, executed).
- `docs/reports/absorbing-boundary/feasibility_cap.ipynb` (Task 1, executed).
- `docs/reports/absorbing-boundary/{cap_toy.py,build_*.py,send_summary.py,*.png}`.
- `ResearchProject/systems/vacuum/runs/epsilon_grid.csv` + 72 run dirs.
- `ResearchProject/systems/vacuum/{mfa_sweep/run.cpp, dispatch_sweep.py, gpu_probe.cu, tests/gate1_mask_absorber/}`.
- `inq-stack/include/inqkit/absorbers/mask_absorber.hpp`.

### Added 2026-06-13 (later): SES feasibility + inq-immutable rule
- **SES (smooth exterior complex scaling) feasibility study** —
  `docs/reports/absorbing-boundary/feasibility_ses.ipynb` (+ `ses_toy.py`,
  `build_ses_feasibility.py`). Verdict: FEASIBLE but the HARDEST of the three
  families in INQ — SES modifies the KINETIC operator (complex metric), not just
  adds a potential, and INQ's FFT kinetic assumes a uniform real grid. Practical
  route = SES-CAP (paper Eq. 27, `V0+V1∂x+V2∂x²`) via FFT-applied derivatives +
  Crank–Nicolson. Difficulty order: MFA (done) < CAP < SES. NumPy FD+CN toy
  reproduces the paper's Fig. 9 (largest absorption region; ε≲0.01 at L=30).
- **New rule** `.claude/rules/inq-immutable.md`: `inq/` is strictly read-only;
  engine changes go in `inq-study/`, wrapper code in `inq-stack/`.

### Optional follow-ups (NOT blocking; for a future session)
- If paper-exact ε calibration is wanted, run a dt→0 convergence study (expensive)
  or match dx=0.1 (needs a non-renormalising fine-grid propagator).
- Task-1 CAP is "GO" as a separate engineering task (CN + complex PotentialType in inq-study).
- `inq/` and `inq-study/` remain byte-identical (verify before any commit).

---

## Milestone: 2026-06-13 (later) — Implementation built; gates 1 & 2 passing; key propagator findings

### Current status
Mask absorber implemented, gate-1 PASSED (all 9 checks), gate-2 (pilot) passing
after fixing two propagator issues found by running. The validated production
binary is `ResearchProject/systems/vacuum/mfa_sweep/run` (env-driven). The
72-run sweep is ready to dispatch. Task-1 feasibility notebook built+executed.

### Critical findings (must survive)
- **Mask mechanism = in-callback mutation, verified** (gate-1 fidelity bit-identical, feedthrough norm-drop).
- **Propagator MUST be ETRS, NOT Crank–Nicolson.** INQ's CN **renormalises the WP
  orbital to unit norm every step**, which silently UNDOES the mask's absorption
  (observed: norm jumped 0.12→1.0 at step 1, ε came out 8.2). ETRS does not
  renormalise → the mask's norm reduction (absorption) is preserved.
- **ETRS Taylor needs dt·E_max ≲ 2** → grid spacing is clamped: `dx=clamp(0.75/k0, 0.18, 0.30)`.
  dx=0.15 (E_max≈219, dt·E_max≈2.2) ABORTS; floor 0.18 (E_max≈152, 1.5) is safe.
  Adaptive dx (≈0.75/k0) also cuts low-k0 cost since wide packets are over-resolved.
- **Transverse box can be tiny (NPERP=12, Lperp=12·dx).** The injected WP is
  transverse-truncated (norm_after≈0.35, N0≈0.12) but ε=inner/N0 **cancels the
  transverse factor exactly**: NPERP=24 gave ε identical to NPERP=12 (0.0856481).
  σ MUST stay pinned to 4√2/k0 (ε-formula condition).
- **Anchors DROPPED.** INQ finite-cell hard-wall reflection is messy (ε≈0.43, not
  1); the paper's ε plateaus below 1 anyway. Low-E masked runs give the high-ε
  end; the independent `cap_toy` 1D split-operator reference validates the pipeline.
- **Pilot ε cross-validated:** k0=1 (E=13.6 eV), L=20 → ε=0.086, matching the
  cap_toy mask value (~0.08) with NO shared code. Norm trajectory shows correct
  physics: 0.12 → packet travels +z → absorbed (norm→0.011, 91%) → 8.6% reflects.

### Files (new/changed this milestone)
- `inq-stack/include/inqkit/absorbers/mask_absorber.hpp` — MaskAbsorber + inner_region_norm (ε reducer).
- `ResearchProject/systems/vacuum/mfa_sweep/run.cpp` — env-driven production binary (ETRS, adaptive dx).
- `ResearchProject/systems/vacuum/tests/gate1_mask_absorber/run.cpp` — gate-1 (4 checks, PASS).
- `ResearchProject/systems/vacuum/dispatch_sweep.py` — 2-GPU dispatcher (72 jobs, anchors dropped).
- `ResearchProject/systems/vacuum/gpu_probe.cu` — cudaMemGetInfo probe (NVML broken).
- `docs/reports/absorbing-boundary/{cap_toy.py, build_feasibility.py, build_reflectivity.py, feasibility_cap.ipynb}`.
- `docs/validation/test-catalogue.md` — gate-1 + agent rows.

### Exact next steps
1. Confirm final gate-2 (dx-convergence k0=1 dx0.20 vs 0.30; high-k0 k0=6 stable).
2. `cd ResearchProject/systems/vacuum && venv/bin/python3 dispatch_sweep.py` (72 runs, 2 GPUs).
3. `cd docs/reports/absorbing-boundary && venv/bin/python3 build_reflectivity.py` (Task-2 notebook).
4. Summary email; final handover update.

---

## Milestone: 2026-06-13 — Grill complete, plan crystallised, mask mechanism verified on GPU

### Current status
Grill-with-docs session finished; the overnight plan is fully crystallised and
user-approved. The single make-or-break feasibility item — applying the mask
each step from the inq-stack wrapper without editing the engine — is **verified
on GPU** (not just reasoned). No implementation of the mask absorber, ε reducer,
configs, or notebooks has started yet. Ready to execute the gated pipeline in
`docs/plans/absorbing-boundary.md` §10.

### Key verified result (the foundation)
In-callback mutation IS Eq. 12 and works. Standalone test
`ResearchProject/systems/absorbing-boundary/mask_mechanism_check/run.cpp`
built+ran on GPU (exit 0):
- **FIDELITY PASS**: M≡1 callback ⇒ `|ΔN|=0.00e+00`, `|Δz|=0.00e+00` vs no-callback baseline (bit-identical).
- **FEEDTHROUGH PASS**: sin² absorber callback ⇒ WP norm 0.999999 → 0.018999 (drop 0.981); surviving remnant at the absorber edge.
Mechanism: the `real_time::propagate` callback captures the outer **non-const**
`electrons` (propagate holds it by reference; callback fires after each ETRS
step) → multiply `electrons.kpin()` by M(z) in the callback. INQ's `viewables`
observer is const and is NOT used for mutation. `inq/` + `inq-study` stay byte-identical.

### What changed (this session)
- `CONTEXT.md`: added "Absorbing boundaries" glossary section (ε, MFA/CAP, box
  geometry, τ, anchor run, quasi-1D emulation, mask mechanism, showcase run,
  vacuum excitation metric) + generalised the header note.
- Created `docs/plans/absorbing-boundary.md` (the full overnight plan).
- Created the mask-mechanism verifier (run.cpp) and ran it (built artefacts present).

### Files touched
- `/local/data/public/skcb2/tddft/CONTEXT.md` — glossary section appended.
- `/local/data/public/skcb2/tddft/docs/plans/absorbing-boundary.md` — plan (new).
- `/local/data/public/skcb2/tddft/docs/handovers/absorbing-boundary.md` — this file.
- `/local/data/public/skcb2/tddft/ResearchProject/systems/absorbing-boundary/mask_mechanism_check/run.cpp`
  — mechanism verifier (TO BE RELOCATED to `ResearchProject/systems/vacuum/tests/mask_mechanism_check/`).

### Commands run
```bash
# build + run the mask mechanism verifier on GPU (from its dir)
cd ResearchProject/systems/absorbing-boundary/mask_mechanism_check && inq-run
# verdict: FIDELITY PASS, FEEDTHROUGH PASS (exit 0)
```

### Tests and validation
- Proposed (gate-1): mechanism (done), mask-shape unit, ε-reducer engine unit.
- Approved: gated pipeline (Task 1 autonomous; Task 2 sweep blocked on gates).
- Run: mechanism check — **PASS** (fidelity + feedthrough).
- Outcomes: mechanism locked; mask-shape + ε-reducer tests NOT yet written.
- Remaining gaps: ε formula must clear the formula-validation agent; mask-shape
  + ε-reducer tests + test-validation agent; gate-2 pilot (incl. transverse-
  insensitivity check) before the sweep.

### Trusted sources used
- De Giovannini, Larsen & Rubio 2014 (arXiv:1409.1689) — method, Eq. 12 (mask),
  Eq. 13 (sin² mask), Eq. 7/8 (ε), geometry (X=6σ, x₀=−3σ, τ=2(3σ+L)/k₀),
  σ=4√2/k₀. PDF validated page-by-page against the prompt's extraction — accurate.

### Attribution notes
- `mask_mechanism_check/run.cpp` mask + ε loops adapted from the injection/norm
  loops in `inq-stack/include/inqkit/wavepacket/wavepacket.hpp:213–257`.
- Free-WP recipe (non_interacting, ghost via extra_electrons(2.0)) from
  `inq-stack/tests/include/inqkit/wavepacket/test_free_wp_engine.cpp`.

### Known issues / blockers
- **NVML/`nvidia-smi` broken** (driver mismatch) — compute works; dispatcher must
  use `cudaMemGetInfo` probe, NOT nvidia-smi. Other users (lm2153, fb638) had
  long-lived python procs at grill time → single-GPU fallback + WARN if occupied.

### Assumptions still in play
- Free-particle factorization (D5): 3D inner-region ε = 1D paper ε for any
  transverse box. EMPIRICAL — gate-2 pilot must confirm ε-vs-Lperp insensitivity.
  If it fails, STOP and reconsider.
- Periodic cell is equivalent to the paper's hard walls because M=0 at the right
  cell edge and the run stops before the reflected wave reaches the left wall.
- Ghost 2e occupied state never perturbs the WP (non-interacting). Density
  "total"/"system" variants are dominated by the flat ghost; "wp" is the physics.

### Exact next steps
1. Relocate `mask_mechanism_check/` → `ResearchProject/systems/vacuum/tests/` (delete the `absorbing-boundary/` scratch dir).
2. Write `inq-stack/include/inqkit/absorbers/mask_absorber.hpp` + the
   `inner_region_norm` ε reducer. Send ε = ∫|ψ|²·𝟙[z<z_abs0] to the
   **formula-validation agent**; lock only on agreement.
3. Write gate-1 unit tests (mask-shape pure; ε-reducer engine known-case) →
   **test-validation agent** → run → record in `docs/validation/test-catalogue.md`.
   Gate-1 = mechanism PASS ∧ shape PASS ∧ ε PASS.
4. Build `systems/vacuum/shared/configs/mfa_mask_sin2.hpp` + run generator
   (geometry formulas: plan §2). Run gate-2 pilot (k₀≈1.5, L=20): norm∈[0.95,1.05],
   |ΔE|<1 mHa, ε∈(0,1), ε insensitive to Lperp.
5. Gates pass ⇒ dispatch 76 runs (L∈{5,10,20,30,40,50} × 12 k₀ giving E≈0.5–490 eV
   + 4 anchors) on 2 GPUs; catalogue each (`tddft-run-catalogue`).
6. `docs/reports/absorbing-boundary/mfa_reflectivity_study.ipynb` (Task 2 curves).
7. `docs/reports/absorbing-boundary/feasibility_cap.ipynb` (Task 1, Level 2:
   source-grounded analysis + NumPy 1D split-operator sin² CAP toy demo).
8. Summary email to chiddukanna@gmail.com (ε(E,L) PNG, counts, wall-time, failures);
   failure emails on gate-1/gate-2 halts.

### Do-not-repeat / guardrails
- NEVER edit `inq/` or `inq-study/` for Task 2 (mechanism needs no engine edit).
- Do NOT poll nvidia-smi (broken). Do NOT lower k₀_min below 0.2 (cost ~1/k₀³).
- venv python for all post-processing; figures `.png` via canonical theme.
