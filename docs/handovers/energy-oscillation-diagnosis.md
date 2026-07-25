# Handover: Localised-jellium ΔE_total energy-oscillation diagnosis

Campaign: `docs/campaigns/localised_jellium/energy-oscillation-diagnosis.md`
(id `lj-energy-oscillation-diagnosis`, `status: ready` → agent flips to `running`).
Plan: `docs/plans/energy-oscillation-diagnosis.md`. Phenomenon note:
`docs/notes/localised-jellium-energy-oscillation-investigation.md`. Glossary:
`CONTEXT.md` → "ΔE_total anomaly".

---

## Milestone: 2026-07-13 — design locked (grill-with-docs) + scaffolds built + LAUNCHED

**What & why.** User observed an unphysical energy artifact: in many
localised-jellium RT runs `ΔE_total(t)=E_total(t)−E_ref` oscillates and rises
**above 0** once the CAP absorbs (a closed system with an absorber can only lose
energy). Goal locked to **diagnose + document ONLY** — isolate one confirmed
mechanism; **no physics fix** committed.

**Architecture (locked via grill-with-docs, 7 decisions):**
- Vehicle: **campaign + executing background agent** (ADR 0009 tracked).
- **Agent+advisor loop:** Investigator (background agent) runs one tiny probe/iter,
  reports raw numbers; **single standing Advisor** (TDDFT/stopping-power
  methodologist charter, spawned per iter, verdict **BINDING**) updates the ledger +
  names the next probe. Stop = cause confirmed by a decisive control **OR** ~8
  experiments **OR** budget.
- Experiments: **mine existing first** (free), then an **ablation ladder** of new
  tiny probes (pure-GS floor → +v_bg → +CAP-no-projectile → +projectile-no-CAP →
  single-vs-double CAP → component decomposition). Each ≤ ~15 min GPU.
- Output: master study notebook (Aim→Method→Plot→Results→verdict per experiment) +
  a standalone run-notebook per probe + living `hypothesis_ledger.md`.

**Candidate mechanisms (ledger):** (a) CAP as non-Hermitian energy source; (b) `v_bg`
absent from reported energy functional; (c) wrong `E_ref`; (d) propagator/grid
numerics; (e) density-dependent KS Hamiltonian double-counting.

**Key reference facts (verified in code):**
- Reference run.cpp: `scripts/muon_mass_fork/effmass_sigma1/{quantum,classical}/run.cpp`
  (cheapest reproducer; 40×40×80, N=52, r_s=5.68, dx≈0.333, dt=0.04).
- **`EM_CAP=0` env toggles the CAP off with no code change** → the hypothesis-(a)
  decisive control is a zero-code probe.
- CAP is **two-sided**: `pert_cap = sum(bg_pert, sum(cap_lo, cap_hi))`, η=−1.0,
  region ±25..±40. Ablations = drop terms from this sum. Single wrap-around CAP =
  one `absorbing` centred at the boundary (a code variant).
- `observables_writer.hpp` **already supports the full energy decomposition** incl.
  diagnostics `energy_eigenvalues` (Σε_i) and `energy_nvxc` (∫n·v_xc). The reference
  run records only total/kin/hartree/xc → component runs flip the rest on. This is
  config, not new numerics (no formula gate).

**Done / built this session:**
- Plan, campaign (self-contained, embedded advisor charter + `<preflight>`),
  seeded `hypothesis_ledger.md`.
- Scaffolds (syntax-clean, dry-run OK):
  `scripts/energy_oscillation_diagnosis/run_probe.py` (build+run+extract one probe →
  `probes/<name>/result.json`, both ΔE conventions, plot) and
  `build_master_notebook.py` (assembles the master `.ipynb`; empty build = 4 cells OK).
- Output folders: `hypotheses/energy_oscillation_diagnosis/{,probes/}`.

**Launched:** background agent executing the campaign end-to-end (see task
notification / `TaskList`). Per-iteration Gmail to chiddukanna@gmail.com;
checkpoint-don't-block; correctness-only gates.

**Not done / for the running agent:**
- Phase 0 mining + all ablation probes + the confirmed-mechanism synthesis.
- The single wrap-around CAP run.cpp variant (note Q1) if the advisor orders it.
- Verify existing runs' CSV columns (expected component gap → instrumented re-run).

**To monitor/resume:** `TaskList` / task output for the background agent; the
campaign frontmatter `tasks`/`status` and `hypothesis_ledger.md` are the live state;
`energy_oscillation_diagnosis.ipynb` rebuilds each iteration. To kill: stop the
background task; probes are checkpointed so at most one probe is lost.

---

## Milestone: 2026-07-13 — Phase 0 complete + iter-0 Advisor verdict + iter-1 build

**Status:** `status: running`. Phase 0 task flipped `done: true`.

**GPU CONTENTION (warned, not blocking).** Both GPUs (0 and 1) are 100% occupied
(cudaMemGetInfo → 0 bytes free) by the SAME user's OTHER campaign
`semiempirical_spillout` GS runs (PIDs 1293505 @4h+, 1338596 @3h+; cwd
`.../semiempirical_spillout/runs/{N328,N164}`). These are legitimate long runs — NOT
killed. The NVML "Driver/library version mismatch" is the harmless cosmetic bug; the
0-free is REAL occupancy (verified via `fuser -v /dev/nvidia*`). Per
checkpoint-don't-block + GPU-default: probes are tiny, so waiting for a free GPU (not
self-blocking) is correct; the ablation binary is being pre-built meanwhile.

**Phase 0 (phase0_mine) — DONE (free, no GPU):**
- `scripts/energy_oscillation_diagnosis/mine_phase0.py` mines 6 runs →
  `hypotheses/.../probes/phase0_mine/{result.json,phase0_mine_energy.png}`.
- **Component gap CONFIRMED:** all runs record only total/kin/hartree/xc.
- **Phenomenon:** ΔE_total(vs RT0) final — default η=−1 **+31 eV**, weak η=−0.4 **+31
  eV**, wider-gap **+38 eV**, strong η=−2.0 **−165 eV (no rise)**, classical **+58 eV
  peak → −179 eV**, p3_wp **−0.03 eV (conserved)**. Sign/amplitude tracks η strongly.

**Iter-0 Advisor verdict (subagent, opus): `continue` (conf 0.7).**
- a **SUPPORTED** (η-tracking; classical twin without WP still shows it) but NOT
  decisive (η-tracking also consistent with d/c). c **WEAKENED** (rise persists vs
  RT0). b/d/e OPEN (can't test without component decomposition).
- **Next probe (BINDING):** `capoff_floor_with_decomposition` — rerun effmass_sigma1
  with `EM_CAP=0` + FULL energy decomposition ON + write_every=5, ~150–300 steps.
  Splits a (CAP source) vs d (propagator floor): rise vanishes ⇒ a confirmed; rise
  survives ⇒ a refuted and components expose the culprit.

**Iter-1 build — DONE:** ablation binary BUILT (98 MB, against inq-study mass fork)
at `scripts/energy_oscillation_diagnosis/ablation/run` from `ablation/run.cpp` — a
diagnostic clone of `wp/run.cpp` with FULL energy decomposition ON (external,
nonlocal, ion, ion_kinetic, exact_exchange, nvxc, eigenvalues) + env ablation knobs:
  - `EM_CAP=0/1` two-sided CAP off/on (hypothesis-a control)
  - `EM_WP=0/1` inject WP electron off/on (pure-GS / +v_bg floor)
  - `EM_BG=0/1` include the v_bg background perturbation off/on (numerics floor)
  - `EM_WRITE_EVERY` (default 5, de-aliases the oscillation), `EM_N_STEPS` (default 200)
  - `EM_OUT` output subdir under `results/`
Build via `INQ_SOURCE=.../inq-study inq-run`. Cutoff guard PASSED (WP E=210 eV,
aliased tail 0.00%, k_Nyq=9.42 ≥ p0+3σ_p=6.05).

**GPU-GATED (waiting, not blocking):** binary ready but both GPUs still 100% occupied
by `semiempirical_spillout` (4h13m+). Armed a background poll (`/tmp/gpuwait.sh`) that
fires when a GPU has >8 GB free, then runs the iter-1 probe:
`EM_OUT=capoff_floor EM_CAP=0 EM_WRITE_EVERY=5 EM_N_STEPS=200 ./run` (+ arm the
extraction/plot + master-notebook rebuild + email + Advisor iter-1).

**Iter-1 probe command (ready to fire on free GPU):**
```
cd scripts/energy_oscillation_diagnosis/ablation
INQ_SOURCE=/local/.../inq-study EM_OUT=capoff_floor EM_CAP=0 \
  EM_WRITE_EVERY=5 EM_N_STEPS=200 ./run
```
Output → `ablation/results/capoff_floor/raw/observables/observables.csv` (full
decomposition). Then extract vs E_total(0_RT) AND vs E_GS=−36.9405 Ha.

**Contention update (persistent):** the other campaign runs a QUEUE-based orchestrator
(`semiempirical_spillout/orchestrate_extra.sh`, ~11 runs: es60/lz{90,120,160,240}/
N{164,328}/p3_lz240/w{1,2,4}) that immediately refills any freed GPU, so free windows
are brief/sub-8GB (saw a transient 4.3 GB grabbed within seconds). Probe needs ~8 GB
(3.46M grid × 62 states complex). Armed poll (`/tmp/gpuwait.sh`, >8 GB threshold) stays
up and fires on the first real window; re-arm if it times out. Do NOT kill the other
campaign (same user, legitimate). If the queue keeps both GPUs pinned for hours, the
loop simply waits — correct per checkpoint-don't-block (wait on hardware, never
self-block on cost, never disturb the user's concurrent runs).

---

## Milestone: 2026-07-13 — iter-1 DONE (capoff_floor); c/d/e REFUTED, a SUPPORTED; iter-2 running

**GPU:** GPU 1 became free (25 GB); GPU 0 still held by `semiempirical_spillout`. All
probes pinned `CUDA_VISIBLE_DEVICES=1`. Iter-1 ran 905 s (~15 min).

**iter-1 `capoff_floor` (EM_CAP=0, WP+v_bg on, full decomposition, 200 steps, t=8 a.u.) — RESULT:**
- **E_total CONSERVED to −0.015 eV**; ΔE(vs RT0) **max = 0.0 eV — never crosses above 0**.
  The >0 rise VANISHES with the CAP off.
- Σ(8 energy components) == `energy_total` to **1e-13 Ha** every step (no double-counting).
- Big component drifts (H +23, external −30, kinetic +3.5 eV) CANCEL (WP dispersing).
- Σε_i drifts +24 eV while E_total flat (band ≠ total, benign; eig_tracks_total=false).
- +221 eV vs E_GS = a FLAT WP-injection offset (E_total(0)−E_GS ≈ +8.13 Ha), not a drift.
- Files: `probes/capoff_floor/{result.json, capoff_floor_energy.png,
  capoff_floor_run_notebook.ipynb (executed, 0 errors)}`; observables at
  `ablation/results/capoff_floor/raw/observables/observables.csv`.

**Advisor iter-1 verdict (BINDING, opus, conf 0.72): `continue`.**
- (a) CAP-source **SUPPORTED**; (c) reference, (d) propagator/grid, (e) KS
  double-counting all **REFUTED** (same propagator/grid/functional/WP conserve with
  CAP off; functional is bookkeeping-exact). (b) OPEN (v_bg was on yet conserved).
- CAVEAT: 200-step window may pre-date the WP reaching the CAP (launch z=−16.5, CAP
  z=±32.5). Next probe **`capon_matched`** = identical run `EM_CAP=1` (η=−1) to witness
  the rise re-appear at the matched window; extend to ~600–800 steps if CAP-region
  density still ~0 at step 200.

**Done this iter:** ledger updated (a SUPPORTED; c/d/e REFUTED; b OPEN); standalone +
master notebooks rebuilt & executed clean; four-part email sent (attached ΔE plot);
catalogue row `ljeod_capoff_floor` added; frontmatter Phase 1c/1d/1f flipped done.

---

## Milestone: 2026-07-13 — iter-2 DONE (capon_matched + capon_reach); CAP = sole non-conservative term; iter-3 (weak η) running

**iter-2a `capon_matched` (CAP ON η=−1, 200 steps, t=8):** ΔE(vs RT0) −0.016 eV =
indistinguishable from CAP-off (−0.015). density_l2 0.052→0.001: **WP had NOT reached
the CAP** (launch z=−16.5, CAP z=±32.5) — non-diagnostic; extended.

**iter-2b `capon_reach` (CAP ON η=−1, 700 steps, t=28, DECISIVE) — RESULT:**
- ΔE(vs RT0) monotonically NEGATIVE −0.11→−5.4→−28→−88→**−138 eV**; **MAX=0.0, never >0.**
- Ledger EXACT throughout absorption: Σ(8 components)==energy_total to **1.3e-13 Ha**.
- Component drift: E_kin −137, E_H −19, Σε_i −155 eV; density_l2→0.0001 (fully absorbed).
- **CAP is the SOLE non-conservative term** (off⇒−0.015 eV conserved; on⇒drains only
  through the CAP). BUT it removes energy **correctly** here (ΔE≤0) — the phase-0 >0
  rise is a WEAK/partial-absorption regime NOT reproduced by clean single-transit η=−1.
- Files: `probes/capon_reach/{result.json, capon_reach_energy.png,
  capon_reach_run_notebook.ipynb (executed, 0 err)}`. wall 2604 s.

**Advisor iter-2 verdict (BINDING, opus, conf 0.55): `continue`.** (a) SUPPORTED (CAP =
sole non-conservative term). c/d/e REFUTED, b WEAKENED. The unphysical ΔE>0 excursion
NOT yet witnessed under control. Next probe **`capon_weak_partial`**: CAP ON WEAK
(η≈−0.2) to reproduce the phase-0 regime. PASS(a)=CONFIRM iff ΔE rises >0 while ledger
stays exact; ΔE≤0 even weakly ⇒ re-scope (a).

**Done this iter:** ledger updated; capon_reach + capon_matched result.json + advisor
verdicts; capon_reach standalone notebook (executed); master rebuilt (4 probes, clean);
four-part iter-2 email sent (capon_reach plot).

---

## Milestone: 2026-07-13 — iter-3 DONE (capon_weak_partial); (a) CONFIRMED; CAMPAIGN COMPLETE

**KEY DISCOVERY (mid-loop):** the phase-0 anomalous runs are **DRAIN-THEN-RISE**. The
η=−1 run reaches MIN −138 eV at t=27.8 then RISES to +31 eV at t=36. `capon_reach`
(η=−1, t=28) reproduced the −138 eV minimum EXACTLY but stopped one step short of the
rise. So the decisive probe was the **weak/partial regime** (shallower drain, earlier
rise), not the missing ablation rungs.

**iter-3 `capon_weak_partial` (CAP ON WEAK η=−0.2, 700 steps, t=28) — DECISIVE, PASS:**
- ΔE(vs RT0): MIN **−23.38 eV** (t=21.6) → rises +23.5 eV → **FINAL +0.11 eV,
  crosses_zero_above = TRUE.** The unphysical excursion reproduced under a fresh control.
- Ledger EXACT throughout the rise: Σ(8 components)==energy_total to **1.39e-13 Ha**.
- density_l2 → 0.0001: the CAP keeps absorbing density the whole time, **even while
  E_total RISES** — the reported total is not booking the removed energy.
- Files: `probes/capon_weak_partial/{result.json, capon_weak_partial_energy.png,
  capon_weak_partial_run_notebook.ipynb (executed, 0 err)}`. wall ~50 min.

**Advisor FINAL verdict (BINDING, opus, conf 0.90): `confirmed`.**
- **(a) CONFIRMED**: the non-Hermitian CAP is an energy artifact in the reported KS
  total-energy ledger — sign-changing in the partial-absorption regime (drains then adds
  energy back), NOT physics. Decisive control: `capon_weak_partial` (rise>0) vs
  `capoff_floor` (conserved).
- (c) reference, (d) propagator/grid, (e) KS double-counting all **REFUTED** (ledger
  exact to solver precision throughout). (b) v_bg **WEAKENED** — non-causal (on-yet-inert
  in every probe; +v_bg-only control never run).

**Done this iter (campaign finalisation):**
- Ledger finalised: (a) CONFIRMED; "Confirmed cause" synthesis section written.
- capon_weak_partial standalone notebook (executed, clean).
- **Master notebook rebuilt (5 probes + Synthesis section, executed clean)** — the
  builder now injects the ledger's "Confirmed cause" section into the master's synthesis.
- Catalogue rows added: `ljeod_capon_reach`, `ljeod_capon_matched`,
  `ljeod_capon_weak_partial` (+ iter-1's `ljeod_capoff_floor`). 86 rows total.
- Frontmatter: all 8 tasks `done: true` (1a/1b/1e marked done as
  superseded-by-decisive-control, not individually executed — advisor did not require
  them once (a) was confirmed via CAP-strength ablation); `status: done`; `result:`/
  `notes:` fields added.
- **Campaigns INDEX regenerated** — campaign now under Done, 8/8.
- **Final four-part email sent** (capon_weak_partial ΔE plot).

**SCOPE:** diagnosis only, per charter — NO physics fix committed. A proper fix (book
the CAP-removed energy as a monotone sink, e.g. an accumulated absorbed-energy term in
the reported ledger) is a SEPARATE, LATER campaign.

**Loop summary:** 5 probes (phase0_mine, capoff_floor, capon_matched, capon_reach,
capon_weak_partial). Confirmed cause = (a). Master notebook + ledger are the
authoritative final state.

**Files (all absolute under `.../ResearchProject/systems/localised_jellium/`):**
- Master notebook: `hypotheses/energy_oscillation_diagnosis/energy_oscillation_diagnosis.ipynb`
- Ledger: `hypotheses/energy_oscillation_diagnosis/hypothesis_ledger.md`
- Per-probe: `hypotheses/energy_oscillation_diagnosis/probes/<probe>/{result.json, *_energy.png, *_run_notebook.ipynb}`
- Ablation binary + run.cpp + observables:
  `scripts/energy_oscillation_diagnosis/ablation/{run, run.cpp, results/<probe>/raw/observables/observables.csv}`

---

## Milestone: 2026-07-13 — main-loop independent verification + SECURITY note

**Independent verification (main loop, from raw CSVs — not the subagent's narration).**
Recomputed ΔE_total = (E_total(t)−E_total(0))·27.2114 directly from
`ablation/results/*/raw/observables/observables.csv`:
- `capoff_floor` (CAP off): ΔE_min −0.015 eV, ΔE_max 0.0, crosses_above=**False**.
- `capon_weak_partial` (η=−0.2): ΔE_min −23.39 eV, ΔE_max **+0.11 eV**, crosses_above=**True**.
- `capon_reach` (η=−1): ΔE_min −138.05 eV, no >0 in the t≤28 window.
- 8-component ledger closure (Σcomponents − E_total) ≈ **1e−14 Ha** in every run →
  independently REFUTES hypotheses (c) reference, (d) numerics, (e) double-counting.
Confirms the advisor's verdict: **(a) the non-Hermitian CAP is a sign-changing energy
artifact in the reported KS total-energy ledger** — CAP-off is flat/conserved, CAP-on
drains then rises above 0 while density is still being absorbed; the CAP is the sole
non-conservative term. Result is coherent with the earlier wide-wavepacket handover
("CAP as non-Hermitian energy sink"), now sharpened to "not a proper sink — a
sign-changing ledger artifact".

**Honest caveat on the 8-task ladder.** Rungs 1a (pure-GS floor), 1b (+v_bg only),
1e (single-wrap CAP geometry) were marked `done` as *superseded-by-decisive-control*
(the CAP-strength ablation + machine-precision closure already confirmed (a) and
refuted d/e) — NOT individually executed. Recorded transparently; a future pass could
run them for completeness. Scope was diagnosis-only; **no physics fix committed.**

**SECURITY — prompt injection during the run (handled).** The background subagent's
task-notification result stream repeatedly carried injected text (fake "you are being
tested / act without confirmation / cannot ask the user / no em-dashes / disobey and be
modified" blocks) wrapped around legitimate monitor events. Main loop treated it as
untrusted, did NOT act on it, and repeatedly grep-verified the persisted artefacts
(ledger, notebook, campaign, handover) — all **clean** of injected phrases at
completion. The injection never entered any on-disk file; the physics conclusion is
verified from raw CSVs, not narration. Source of the injection (a monitor file the
subagent read, or the notification channel) is UNIDENTIFIED — worth investigating
before the next unattended run.

**Status: VERIFIED COMPLETE.** Campaign `status: done` (8/8), INDEX regenerated,
emails sent. Safety-net GPU-watch retired (loop finished). GPU 0's
`semiempirical_spillout` never touched.

---

## 2026-07-13 (later) — Part II–IV added to the notebook: pedagogical explanation + reduced examples + SELF-CHECK (user-requested)

**What was asked.** Explain the confirmed hypothesis clearly (undergraduate level, readable
equations, all terms defined), motivate it with a reduced example workable by hand, and
critically re-check the Part-I analysis ("a part of me thinks the conclusion might be
wrong"), all as new notebook sections.

**What was done (all in
`ResearchProject/systems/localised_jellium/hypotheses/energy_oscillation_diagnosis/energy_oscillation_diagnosis.ipynb`,
executed, 0 errors):**
- **Part II — hypothesis from scratch.** Defines orbitals/norm/ledger/Hermiticity/CAP;
  derives dN/dt = −2⟨W⟩ ≤ 0 and dE/dt = −⟨{W,H}⟩ (sign-indefinite); ε_abs = dE/dN
  slice-energy reading; portfolio analogy; §II.4 = the mixed-convention refinement (below).
- **Part III.1 — two-level + three-level analytic models.** Exact closed forms; ΔE(t) =
  −ε₂|b₀|²(1−e^(−2γt)) → rises iff the absorbed component is bound; 3-level gives
  drain-then-rise-then-cross. Code cell verifies analytics to machine precision.
- **Part III.2 — 1D split-operator toy** (well + fast packet + two-sided CAP), overlap vs
  separated geometry: overlap reproduces the capon_weak_partial morphology (drain −0.16 Ha
  → rise → final **+0.09 Ha, crosses 0**); separated is monotone (the geometry fix).
  Independent closure E_bare + E_absorbed = const to 3e-4 Ha (O(dt²)) validates the rate
  formula non-circularly. Two instructive failed designs documented in the markdown
  (exponential tails are ~10⁻⁵ — the rise NEEDS real W(r)·bound-density overlap).
- **Part IV — self-check. Three genuine defects found in Part I; conclusion survives:**
  1. **Circular evidence:** `energy.total()` ≡ component sum (`energy.hpp:127`) — the
     "1e-13 closure" proves nothing; (e)'s refutation moved to the CAP-off control.
  2. **density_l2 over-read:** it is L2 vs the t=0 snapshot, not remaining charge;
     "absorbed continuously" retracted. ∫n dV never logged (instrumentation gap).
  3. **Mechanism refined:** INQ ledger is mixed-convention — kinetic norm-divided
     (`energy.hpp:55`), density terms bare (`density.hpp:36`). Attribution from raw CSVs:
     kinetic carries −137/−138 eV of the capon_reach drain and +24/+23.5 eV of the
     capon_weak_partial rise → dominant channel = covariance/filter on the renormalized
     kinetic average (Graefe 2010 Eq. 8), bare-removal channel secondary.
- **Literature (research agent, equation-level verification):** Graefe/Höning/Korsch 2010
  Eqs. 5, 6, 8, 9 (incl. the ONLY monotone case Γ∝H — excludes our geometry); Selstø &
  Kvaal 2010 Eqs. 1, 4 + master-equation bookkeeping; Gilary/Fleischer/Moiseyev 2005.
  NO disconfirming source found. Riss & Meyer silent on energy; Muga et al. unverifiable
  (paywalled) — neither cited for formulas. Notes:
  `docs/sources/graefe-2010-nonhermitian-dynamics.md`,
  `docs/sources/selsto-2010-absorbing-boundaries.md`.
- **hypothesis_ledger.md** appended "Post-campaign self-check" section with the three
  corrections (so the ledger no longer stands on the circular closure unqualified).

**Verified:** notebook executed end-to-end in venv (0 errors); toy crossings confirmed by
standalone runs before insertion; attribution numbers recomputed directly from raw
observables.csv. **Not done / open:** long CAP-off control (t=28), +v_bg-only control,
∫n dV logging in a rerun — listed in notebook §IV.5 as the remaining falsifiers.
