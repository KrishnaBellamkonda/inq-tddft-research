---
id: cap-jellium-loss-function
area: cap_in_jellium
title: "Loss-function feasibility for wavepacket stopping power"
status: ready
hypothesis: "The loss-function route — extract an interference-free L(q,ω) from a CAP'd jellium TDDFT run (no wrap, B1-subtracted) and integrate to S(v) — is a FEASIBLE way to obtain the stopping power of a WAVEPACKET projectile: S_LF(v) reproduces the directly-measured WP stopping power S_WP = -d<T_WP>/dz (peak-driven structure + v-trend; absolute magnitude provisional), with the classical projectile and analytic Lindhard as references."
handover: docs/handovers/cap-in-jellium-loss-function.md
tasks:
  - { name: "GS reference reuse + load check (gs_L50_cubic_N162_dx0p40 == B0)", done: true }
  - { name: "run.cpp velocity-gauge kick mode + build-once (modes classical/wp/kick vs inq-study CAP)", done: true }
  - { name: "pilot @ E15 (~100 steps, all 3 modes): s/step, real energy, CAP absorbs, kick rings at omega_p", done: false }
  - { name: "9 production runs (classical/wp/kick x E15/E20/E30), T~2000 a.u.", done: false }
  - { name: "[UNGATED 2026-06-25] loss-function extraction via fourier-analysis skill (|n_q|²/q² peak-LOCATOR, complex FFT, Δω=2π/T flagged) + reliability gate R1/R2/R3 + classical-vs-WP L (H1, Aim A)", done: false }
  - { name: "[UNGATED 2026-06-25] direct WP stopping S_WP=-d<T_WP>/dz at E15/E20/E30 vs S_LF(v) from medium L; classical+Lindhard refs (H2 feasibility)", done: false }
  - { name: "[UNGATED 2026-06-25] study notebook (per-velocity comparison rows + feasibility verdict) + handover + frontmatter flips", done: false }
  - { name: "[decision gate, post-verdict] surface to user: replicate this campaign for the localised jellium system? (do not auto-launch)", done: false }
blocked_reason: ""
---

# Loss-function hypothesis check (CAP'd projectile, low velocity)

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository.
</identity>

<description>
This campaign is the **successor to the abandoned 2026-06-01 task**
`docs/handovers/projectile_loss_function_interference.md`. That task tried to
obtain a loss function L(q,ω) from a moving projectile and was rejected with a
hard conclusion: in a periodic box a projectile **must wrap**, and the wrap
produces interference that is "very evidently visible and visually mimics
plasmons" — so **no periodic-box projectile run is interference-free**. The escape
routes proposed then were (1) remove the projectile before it wraps, (2) enlarge
the box, or (3) accept a kick. It also noted a CAP "cannot be used because the
uniform jellium bath fills the boundary too" — **but that was before the CAP
baselines existed.**

The CAP baselines (B0–B3, `docs/campaigns/cap_in_jellium/baseline_runs.md`, DONE)
changed the game: a two-sided sin² CAP **absorbs the projectile at the far
boundary so it never wraps**, the B1 bath-drainage is now characterised and
subtractable, and — crucially — once the projectile is absorbed we can keep
propagating an **isolated, interference-free ringing bath** to resolve ω_p(q).

**Core question (what this campaign is FOR):** is the loss-function route a
**feasible way to obtain the stopping power of a WAVEPACKET projectile**? The
loss function is a *medium* property, so integrating it gives the stopping of an
idealised *point charge* of velocity v (Lindhard). The WP projectile is neither a
point (extended Gaussian σ=0.5) nor classical (momentum spread, quantum
spreading) — so whether the loss-function S(v) actually reproduces the WP's real
stopping power is an open, falsifiable question. We answer it against the WP's
**own directly-measured** stopping power as ground truth.

We test two nested hypotheses:

- **H1 (enabling).** With the projectile absorbed by the CAP (no wrap) and the B1
  bath-drainage subtracted, a **reliable, interference-free** loss function
  L(q,ω) can be extracted from a low-velocity jellium projectile-transit run of a
  reasonable duration (R1–R3 gate). *Precondition for H2.*
- **H2 (core feasibility).** The stopping power
  S_LF(v) = (2/πv²)∫(dq/q)∫₀^{qv} ω L(q,ω) dω derived from that L(q,ω)
  **reproduces the directly-measured WP stopping power** S_WP = −d⟨T_WP⟩/dz
  (⟨T_WP⟩ = ⟨ψ_WP|−½∇²|ψ_WP⟩, `docs/sources/stopping-power-formulae.md:50`), with
  the **classical-projectile** S(v) and **analytic Lindhard** S(v) as references.
  **Two-tier verdict (the deliverable, not a defect):** because the |n_q|²/q²
  method is a validated *peak-locator* (peak positions reliable; absolute
  magnitude NOT — quadratic-vs-linear, missing 4π/q²), feasibility is judged
  separately for (a) **structure + v-trend** — does S_LF track S_WP in shape and
  velocity dependence (achievable), and (b) **absolute magnitude** — is absolute
  S(v) recoverable without FDT calibration (expected NOT). Reporting *which tier
  is feasible* is the campaign's headline result.

**Stated comparison aims (equal-standing probes).** At each velocity we run a
**classical** projectile, a **wavepacket (WP)** projectile, and a **quantum
kick**, all on equal footing. Two explicit aims beyond H1/H2:
- **(Aim A) Classical-vs-WP loss-function comparison at each energy** — how the
  classical and quantum projectiles excite the medium differently is itself a
  result, not a means to an end.
- **(Aim B) The quantum kick is the per-velocity baseline** the classical and WP
  loss functions are compared against (the q=0 plasmon / Santervás-Arranz
  reference).

The output is an executed `.ipynb` per the `notebook-making` house style, with
the loss functions plotted **classical | WP | kick side-by-side, one row per
velocity**, and the setup + key quantities + their values clearly discussed.
Email the results at the end of each phase.

**HARD GATE (see <guard_rails>):** the Fourier / loss-function analysis
(temporal FFT of n_q(t): windowing, detrending, peak extraction) **must not be
performed** until the user-led training campaign
`docs/campaigns/check_logic/check_stopping_power_calculation.md`
(id `check-stopping-power`, tasks 2–3: "Fourier-analysis check" + "encode
deterministic workflow") is complete. The 9 simulations may run before then; the
analysis waits.
</description>

<observables_set>
Reuse the **maximal CAP-baseline observable suite** (already built in inqkit /
inqview during `cap-jellium-baselines`): energy components per step,
total current/dipole/density_l2, density VTI, current-density-field VTI, 9
plane/flux screens, region N(t), per-orbital energy, total-system n(k), E-field
post-kernel. On top of that, the loss-function-critical observables:

- **n_q(t) — the primary loss-function input.** Discrete induced-density Fourier
  amplitudes δn_q(t) = FFT[n(r,t) − n(r,0)] at the box modes
  q_m = m·(2π/L) ẑ, m = 1..6 (q = 0.126..0.754 Bohr⁻¹), columns
  `time_au, m, q_au, re_n_q, im_n_q, abs_n_q` (same schema as the existing
  `run_*/results/analysis/observables/n_q_vs_time.csv`). **Sampled at
  Δt_save ≤ 4 a.u.** (≥ 500 samples over T≈2000 a.u. — matches the proven E15
  cadence; Nyquist ω_max = π/Δt_save ≈ 21 eV ≫ ω_p, FFT bin Δω = 2π/T ≈ 0.086 eV).
  This dedicated lightweight dump is independent of the 300-frame VTI cadence.
- **ΔE(t) — kick-run excess energy.** ΔE(t) = (E_total(t) − E_GS)·Ha→eV, per
  step, for the kick runs — the Santervás-Arranz diagnostic input
  (`docs/sources/santervas-arranz-prr-2025.md`).
- WP runs additionally need the **bath isolation** n_bath = n_total − n_wp (the
  canonical bath density) before n_q is taken; classical runs are already
  pure-bath (the Gaussian projectile is an external potential, not electron
  density).
- **Direct WP stopping power — the H2 ground truth.** S_WP = −d⟨T_WP⟩/dz with
  ⟨T_WP⟩ = ⟨ψ_WP|−½∇²|ψ_WP⟩ (Yao–Schleife drift+spread ⟨T⟩=(⟨p⟩²+σ_p²)/2,
  `stopping-power-formulae.md:50`), measured at each velocity from the WP orbital
  kinetic energy vs the WP centroid z(t) (`center_of_density` track, already
  produced for B3 as `cap_b3_wp_centroid.csv`). Cross-check against −dE_proj/dz.
  Derived in post; **no new in-run code** — but its extraction is GATED (it is a
  stopping-power extraction, trained in `check-stopping-power` task 1).

**No new C++ kernel is introduced** (CAP, n_q extraction, and the velocity-gauge
kick all already exist — see <resolved_decisions>). The only new code is the
**temporal-FFT loss-function analysis**, and that is the GATED item routed
through the `check-stopping-power` training (it IS the code-test +
formula-validation for this method) — never an ad-hoc reimplementation here.
</observables_set>

<resolved_decisions>

<run_matrix kind="locked">
3 velocities × 3 probes = **9 new runs**, all in the CAP box, all reusing the
validated GS `gs_L50_cubic_N162_dx0p40`.

Velocities (low-velocity window only; v = √(2E/27.2114) a.u.; v_F = 0.337):
| Energy | v (a.u.) | v/v_F | Note |
|---|---|---|---|
| E15 | 1.050 | 3.1 | resonant: ω_kin(m=1) ≈ ω_p → strongest plasmon coupling |
| E20 | 1.212 | 3.6 | just above resonance |
| E30 | 1.485 | 4.4 | upper end of the low-v window |

Probes at each velocity:
- **classical** — Gaussian projectile (σ=0.5 Bohr erf charge), CAP ON, .ehrenfest()
  ionic propagator. Pure-bath n_q. Also yields the direct classical S(v) via the
  stopping force.
- **wp** — Gaussian wavepacket (σ=0.5 Bohr), CAP ON. Bath isolated as
  n_total − n_wp before n_q.
- **kick** — velocity-gauge uniform kick, A ∝ v·ẑ (see <quantum_kick>),
  **CAP OFF**. The q=0 plasmon / medium reference. *(Inference: CAP-off chosen
  because a kick has no projectile to absorb and a CAP would only drain the
  reference bath; matches how the historical clean reference run_plasmon_E15 was
  done. Flag at review.)*

Run dirs (ADR-0007 grouped-by-sweep):
`cap_loss_function/run_{classical,wp,kick}_E{15,20,30}/`.
</run_matrix>

<quantum_kick kind="locked-mechanism">
INQ's `perturbations::kick` (`inq/src/perturbations/kick.hpp:60-73`) applies
`exp(i·k·r)` to every orbital — a **uniform momentum boost**, in **velocity
gauge** for periodic cells (`kick.hpp:44-47`: efield→0, vpot = −kick_field). In
jellium (electron sea displaced relative to the FIXED positive background) this
drives the **q=0 plasmon** (the uniform sloshing mode) at ω_p = 3.47 eV — the
standard RT-TDDFT optical kick. The kick amplitude ∝ v, so **one velocity-gauge
kick with A ∝ v·ẑ per velocity**. **No custom code** is required (the deferred
single-q kick — plasmon-detection plan "Run A", task #16 — is NOT needed here).
Diagnostic = Santervás-Arranz protocol (`docs/sources/santervas-arranz-prr-2025.md`):
ΔE(t) → detrend second-half plateau → Hann → 8× zero-pad → |rfft|² → normalise
0–20 eV; "per velocity" tests the α(v) countercurrent / linear-response check.
**Inference to verify in the pilot:** that a uniform velocity-gauge kick excites
a clean q=0 plasmon ring-down at ω_p in INTERACTING jellium (confirm via the
pilot ΔE(t) FFT, do not assume).
</quantum_kick>

<sv_route kind="locked">
**One L → all v (medium-property route).** L(q,ω) is a property of the medium, so
S(v) for any v follows from the single integral
S(v) = (2/πv²)∫(dq/q)∫₀^{qv} ω L(q,ω) dω (Lindhard 1954,
`docs/sources/stopping-power-formulae.md`; Im[−1/ε] = L). Reconstruct S_LF(v) from
the **single best CAP'd-projectile L(q,ω)**; the ≥2 velocities exist to prove
probe-independence (R3), not to build a per-v curve. The **H2 feasibility target
is the direct WP stopping power** S_WP = −d⟨T_WP⟩/dz measured at E15/E20/E30 (3
points): does the single medium-L S_LF(v) curve pass through them? The
**classical S(v)** is **REUSED** from the existing classical stopping-force sweep
(`run_sv_sigma0p5/` + classical energy-sweep runs), and analytic **Lindhard** is
the theory line — both are references, not the target.
CAVEAT (`docs/validation/loss-function-formula-validation.md`, memory
`reference_loss_function_method`): the |n_q(ω)|²/q² method is **quadratic** in
n_q (a plasmon PEAK-LOCATOR), whereas the true loss function −Im[1/ε] is
**linear** in Im χ (missing the 4π/q² Coulomb factor). Peak POSITIONS ω_p(q) are
trustworthy; absolute lineshape / area / cross-q intensity are NOT — so the S(v)
amplitude is shape-normalised unless FDT-calibrated, and H2 is judged on
**peak-driven structure + relative trend**, with absolute magnitude labelled
provisional.
</sv_route>

<reliability_gate kind="locked">
H1 "reliable" is operationalised by three checks (R1+R2 are the gate; R3 is the
strongest evidence; peak SHARPNESS is secondary, not a gate):
- **R1 — physically correct peaks.** ω_p(q) from the CAP'd-projectile L(q,ω)
  lands on the RPA/Bohm-Gross dispersion ω(q)² = ω_p² + (3/5)v_F²q² + q⁴/4 within
  a stated tolerance, with q→0 → ω_p0 = 3.47 eV. (The only claim the method is
  validated for.) Reference numbers per mode m=1..6 are in
  `docs/plans/jellium_plasmon_detection.md` §0.
- **R2 — wrap-line gone.** The CAP run's L(q,ω) shows **no** spurious
  wrap-frequency peak (ω_kin(m) = m·v·q₁) that the **matched no-CAP run does**.
  Negative control = the existing no-CAP `n_q_vs_time.csv` at matched v
  (`run_plasmon_n162_L50_E15`, `run_classical/wp_n162_L50_E20`) — analysis-only,
  no new run.
- **R3 — probe-independence (benchmark-free).** Peak positions ω_p(q) agree
  across the E15/E20/E30 CAP runs where modes overlap (a medium property must be
  v-independent). Agreement → reliable; disagreement → still probe-contaminated.
- Secondary (report, do not gate): Δω = 2π/T per run and peak FWHM; flag any mode
  with FWHM < Δω as resolution-limited rather than failing it.
</reliability_gate>

<geometry kind="locked-inherited">
Inherited verbatim from `cap-jellium-baselines` (NOT re-litigated): cubic L=50
Bohr, N=162, r_s=5.69; ω_p=3.473 eV, T_p=49.2 a.u., k_F=v_F=0.3374 a.u.;
n=1.296e-3 a₀⁻³. CAP = 20 Bohr total, two-sided, 10 Bohr/side, slabs
[−25,−15]∪[+15,+25], free region [−15,+15], η=−0.5 Ha, sin² (inq-study
`perturbations::absorbing` composed two-sided via `perturbations::sum`).
Projectile launch z0 ≈ −13 (4σ inside the −z CAP edge), moving +z, exits through
the +z CAP. Boundary rule references the **CAP edges (±15)**, not the box faces.
</geometry>

<propagator kind="locked-inherited">
Electronic propagator = **ETRS** (INQ default). NEVER crank_nicolson with the CAP
(CN renormalises each step → undoes absorption; see memory
`reference_inq_propagator_mask_absorber`). `.ehrenfest()` is the ionic propagator
for the classical projectile.
</propagator>

<duration_and_io kind="locked">
**Long ring-down, T ≈ 2000 a.u.** (≈ 40 plasmon periods; Δω = 2π/2000·Ha2eV ≈
0.086 eV), dt = 0.020 a.u. ⇒ ~100,000 steps per run. Rationale: after the CAP
absorbs the projectile (~3× t_cross: E15 36 a.u., E20 31 a.u., E30 26 a.u.) the
bath rings interference-free, and R1 (peaks-on-Lindhard) needs the frequency
resolution — which the CAP route now provides at no physics cost (only GPU
time). This supersedes the "3–5× traversal" length used for the baselines.
- n_q(t) dumped every ≤ 4 a.u. (≥500 samples). ΔE(t), scalar observables every
  step. Density + current-density VTI @ **300 frames** (WRITE_EVERY ≈ 333).
- GPU scheduling via `cudaMemGetInfo` probe (NVML/nvidia-smi broken — memory
  `reference_gpu_driver_mismatch`); GPU is the default; warn if a GPU is occupied
  by another user. Dispatcher runs the 9 jobs across available GPUs and emails
  per run (Gmail, `tddft-simulations` skill).
- **COST NOTE:** 9 × ~100k steps is heavy (the E15 plasmon run was ~12 h for 100k
  steps). The pilot gates the real s/step before the full set is launched; if
  s/step is prohibitive, surface it and propose a shorter T (with the Δω penalty
  stated) rather than silently truncating.
</duration_and_io>

<pilot kind="locked">
**PILOT FIRST @ E15, ~100 steps, all 3 modes** before the long runs. Gate:
(i) real s/step on the current GPU (feeds the cost decision);
(ii) energy stays real, total energy + norm decrease smoothly for the CAP
projectile modes (absorption signature), no NaN/complex energy;
(iii) GS loads + slabs absorb;
(iv) **kick mode: ΔE(t) over the pilot window shows the onset of a q=0 ring at
ω_p** (confirms the <quantum_kick> Inference) — if it does not ring, STOP and
surface (the velocity-gauge kick may not couple as expected in interacting
jellium).
This is still early CAP-in-interacting-jellium territory and inq-study Task #7
(engine regression) is OPEN → ALL absorption/loss numbers PROVISIONAL until #7.
</pilot>

<file_placement kind="locked">
ADR-0007 grouped-by-sweep (jellium otherwise flat):
- `scripts/cap_loss_function/` — ONE env-driven run.cpp (modes classical/wp/kick;
  kick mode adds the velocity-gauge kick + CAP-off branch), built ONCE against
  inq-study (CAP engine); `dispatch.py` (multi-GPU, emails per run); `analyse.py`
  (per-run pipeline → REPORT.md, per `feedback_per_run_analyse_py`).
- `cap_loss_function/run_{classical,wp,kick}_E{15,20,30}/` — outputs (logs
  gitignored; provenance only).
- `hypotheses/cap_loss_function/` — study `.ipynb` (auto-built by dispatcher
  tail), combined CSVs, `tests/` (task-specific checks). Figures as PNG, canonical
  theme.
</file_placement>

</resolved_decisions>

<guard_rails>
- **HARD GATE — analysis tasks (Fourier AND stopping-power extraction).** Tasks
  5–7 **MUST NOT run** until `check-stopping-power` (id `check-stopping-power`) is
  DONE: task 1 = stopping-power extraction (covers the direct WP S_WP =
  −d⟨T_WP⟩/dz), tasks 2–3 = Fourier-analysis workflow (covers the temporal FFT of
  n_q(t)/ΔE(t): windowing, detrending, peak extraction, L(q,ω), S_LF(v)). Both the
  stopping-power extraction and the loss-function analysis are the methods the user
  is explicitly reclaiming control of. The 9 simulations (tasks 1–4) MAY run
  before then. If you reach task 5 and the training is incomplete, **STOP and
  surface it** — do not fall back to ad-hoc / out-of-box defaults.
- **Abort conditions:** NaN, complex total energy, GS fails to load, or a GPU
  occupied by another user (warn; do not evict).
- **Pilot-first numeric gate** (see <pilot>): no long run launches until the pilot
  passes, including the kick-rings-at-ω_p check.
- **Boundary + cadence:** 4σ/1σ rule referenced to the CAP edges (±15); density
  + current VTI @ 300 frames; n_q(t) at ≤ 4 a.u.
- **PROVISIONAL:** all absorption/loss numbers are provisional pending inq-study
  Task #7 (engine regression). The |n_q|²/q² loss function is a peak-LOCATOR;
  absolute lineshape/area/cross-q intensity and the S(v) magnitude are NOT trusted
  (peak positions + relative trend only).
- **Open dependency named:** the classical S(v) comparison reuses the existing
  classical sweep; if that sweep is unavailable/incompatible, surface it rather
  than launching a new classical-with-CAP sweep on your own initiative.
</guard_rails>

<tasks>
1. **GS reference reuse + load check.** Reuse `gs_L50_cubic_N162_dx0p40`; confirm
   it loads and its density / energy components match the B0 reference. *(no new
   run; simulation-validation Tier-0.)* Done when the load is confirmed and logged.
2. **run.cpp kick mode + build-once.** Extend the CAP run.cpp (modes
   classical/wp/kick) with the velocity-gauge kick (A ∝ v·ẑ) + CAP-off branch for
   the kick mode; build ONCE against inq-study. *(config addition, not a new
   kernel; build-run skill.)* Done when all three modes dispatch and the binary
   builds clean.
3. **Pilot @ E15 (~100 steps, 3 modes).** Run the <pilot> gate. Done when (i)–(iv)
   pass (incl. kick rings at ω_p) and s/step is recorded for the cost decision.
4. **9 production runs.** classical/wp/kick × E15/E20/E30, T≈2000 a.u., dispatched
   across GPUs, Gmail per run. Done when all 9 report `run_completed = true` with
   n_q(t)/ΔE(t) + VTI written. *(tddft-simulations skill.)*
5. **[GATED on check-stopping-power] Loss-function extraction + reliability gate
   (H1, Aim A).** Using the trained deterministic Fourier workflow only: extract
   L(q,ω) per run (bath-isolate the WP runs first), evaluate **R1, R2, R3** with
   explicit numbers per velocity, render an H1 verdict per velocity, and compare
   the **classical-vs-WP** loss functions at each energy (Aim A). Done when R1–R3
   are reported with tolerances + pass/fail per velocity and the classical/WP L
   comparison is rendered.
6. **[GATED] WP-stopping feasibility (H2 — the core test).** (a) Extract the
   **direct WP stopping power** S_WP = −d⟨T_WP⟩/dz at E15/E20/E30 (trained
   stopping-power workflow; cross-check −dE_proj/dz). (b) Reconstruct S_LF(v) from
   the single best medium L(q,ω). (c) **Feasibility verdict, two-tier:** does
   S_LF(v) pass through the 3 direct WP points in *structure + v-trend* (tier a),
   and is *absolute magnitude* recoverable without FDT calibration (tier b)?
   Overlay classical sweep + analytic Lindhard as references. Done when the
   S_LF(v) curve, the 3 S_WP points, and the explicit per-tier feasibility verdict
   exist.
7. **[GATED] Study notebook + handover + frontmatter.** Executed `.ipynb` in
   `hypotheses/cap_loss_function/` per `notebook-making`: context → formulas (every
   term defined) → full reconstructable setup → linked source files → results →
   takeaway; loss functions **classical | WP | kick in a row per velocity** (Aim
   A + Aim B), plus the S_LF(v)-vs-S_WP feasibility figure + verdict. Update
   `docs/handovers/cap-in-jellium-loss-function.md` and flip the frontmatter
   `done`/`status`. Done when the notebook is executed and the handover +
   frontmatter are current.
8. **[DECISION GATE — post-verdict, do not auto-launch] Replicate for localised
   jellium?** After tasks 6–7 (the H2 feasibility verdict + notebook) are complete,
   surface to the user the decision: *should this entire loss-function-feasibility
   campaign be replicated for the **localised jellium** system* (cf. the active
   `locjel-campaign`, `docs/campaigns/localised_jellium/`)? Present the bulk-jellium
   verdict (which feasibility tier held — structure/trend vs absolute magnitude) as
   the basis, note that the localised geometry changes the boundary/CAP and q-mode
   structure (so feasibility does not transfer for free), and **WAIT for the user's
   go/no-go** — do NOT author or launch a localised-jellium campaign on your own
   initiative. *Done = the decision has been put to the user with the verdict
   summary, and their answer is recorded in the handover.*
</tasks>

<follow_up_gate>
**Replicate for localised jellium? (user decision, post-verdict — task 8).** This
campaign tests loss-function feasibility in BULK periodic jellium. Whether the same
route works for a **localised jellium** geometry (finite slab/cluster — cf. the
active `locjel-campaign`, `docs/campaigns/localised_jellium/`) is a separate,
important question the user wants to decide deliberately rather than by default. On
completing task 7, the agent MUST put this decision to the user (task 8) with the
bulk verdict as evidence and WAIT — it must NOT auto-author or auto-launch a
localised-jellium replication. If the user says yes, a new campaign is authored
then (via the `campaigns` skill), gated like this one on the `fourier-analysis`
skill. This gate is a follow-up/exit decision; it does not block tasks 1–7.
</follow_up_gate>

<rules>
- ALWAYS: GPU by default; ETRS (never CN) with the CAP; reuse the existing GS;
  reuse the existing classical S(v) sweep; keep all loss numbers PROVISIONAL
  (Task #7) and the |n_q|²/q² quantity labelled a peak-locator.
- NEVER: perform the temporal Fourier / peak analysis OR the stopping-power
  extraction (S_LF or direct S_WP) before `check-stopping-power` training is done;
  modify `inq/` (use inq-study for the CAP engine); reimplement the Fourier or
  stopping workflow ad-hoc; invent peak positions, tolerances, or S(v) values;
  claim absolute-magnitude feasibility without FDT calibration.
- The classical-vs-WP comparison (Aim A) and the kick baseline (Aim B) are
  first-class results — present them, do not bury them under H1/H2.
</rules>

<preflight>
- [ ] Intent self-contained: CORE = is the loss function feasible for WP stopping
      power? H1 (reliable L, R1–R3 gate) → H2 (S_LF(v) reproduces direct WP
      S_WP=−d⟨T_WP⟩/dz; two-tier verdict: structure/trend vs absolute magnitude;
      classical + Lindhard as references); every task has a done-criterion.
- [ ] Setup reproducible: L=50/N=162/r_s=5.69; GS = gs_L50_cubic_N162_dx0p40 (task 1
      load-check); ETRS, dt=0.020, T≈2000 a.u. (~100k steps), E∈{15,20,30}; n_q ≤4
      a.u., VTI @300; file placement scripts|cap_loss_function|hypotheses (ADR-0007).
- [ ] New code pre-gated: no new C++ kernel (CAP/n_q/velocity-gauge kick all exist,
      kick.hpp:60-73); the ONLY new code is the temporal-FFT analysis, gated on the
      check-stopping-power training (which IS its code-test + formula-validation).
- [ ] Validation & guard rails: pilot-first gate incl. kick-rings-at-ω_p; abort on
      NaN/complex energy/GPU-occupied; boundary ±15 + 300-frame cadence; PROVISIONAL
      (Task #7); peak-locator caveat; classical-sweep-reuse dependency named.
- [ ] Autonomous mechanics: GPU via cudaMemGetInfo (NVML broken, warn if occupied);
      dispatcher multi-GPU + per-run Gmail; notebook = per-velocity comparison rows,
      auto-built via dispatcher/analyse.py tail; handover pointer present; agent
      flips frontmatter done/status.
- [ ] Grounding: kick.hpp:44-73, santervas-arranz-prr-2025, stopping-power-formulae
      (Lindhard), jellium_plasmon_detection §0 (peak numbers), loss-function-formula
      -validation (quadratic-vs-linear caveat); Inferences (kick CAP-off; q=0 ring)
      labelled and pilot-verified.
- [ ] HARD GATE re-checked: check-stopping-power tasks 2–3 DONE before any task 5–7.
</preflight>
