---
# ROUGH DRAFT — authored interactively on 2026-06-27 to capture intent only.
# NOT autonomy-ready. Consolidates the existing sketch in `draft_campaigns.md`
# ("Campaign 2 — large rigid σ") + task P2.3 of the quantum-stopping-power
# campaign + the user's 2026-06-27 mind-dump. To be built out carefully later via
# the `campaigns` skill. Numbers below are starting points, not locked decisions.
id: wide-wavepacket-lowspread
area: localised_jellium
title: "Wide low-spread wavepacket — isolating purely quantum stopping (localised slab)"
status: running
hypothesis: "With a wide, near-rigid wavepacket (σ chosen via the Gaussian spreading law so it does not disperse appreciably over transit) and a matched-σ classical projectile through the localised jellium slab, any WP−classical difference in stopping is attributable to purely quantum effects (Pauli + interference) rather than to dispersion or an interaction-range mismatch."
handover: docs/handovers/wide-wavepacket-lowspread.md
tasks:
  # --- Design (this /campaigns session) ---
  - { name: "D  Design locked: sigma_WP=3.5, E=300, box 50x50x101, dx=0.40, CAP=phase-5 (eta-0.7/10-side), 4-sigma launch, extensive observable suite, Phase-0/1, 6-energy grid", done: true }
  # --- Setup: build-once validated blocks ---
  - { name: "S1 Matched classical UPF electron_gaussian_wpsigma3p5.upf (sigma_pot=2.475); verified by V(r) + known-case generator test", done: true }
  - { name: "S2 Slab GS (LZ=101, dx=0.40, N=82, periodicity 3): build+run+validate (interior n0~1.31e-3, SCF converged, energy stable)", done: true }
  # --- Phase 0: human-in-the-loop design gate ---
  - { name: "P0a Wide-sigma CAP-completeness smoke (WP-only): residual inner norm <2%, N_total drift <2%", done: false }
  - { name: "P0b First matched WP+classical run (E=300); quantum S extracts cleanly + energy/centroid agree -> user sign-off", done: false }
  # --- Phase 1: autonomous production ---
  - { name: "P1a Sweep machinery: env-driven WP+classical run.cpp (full observable suite) + Python orchestrator (resumable, 2-GPU, per-run email, retry)", done: true }
  - { name: "P1b Autonomous E sweep: 6 energies {200,280,360,440,520,600} x (WP+classical) = 12 runs + 1 vacuum-WP SIE control", done: false }
  - { name: "P1c S(E) overlay (WP vs classical vs Lindhard vs sigma=0.5 reference) + WP-classical quantum component (SIE bounded) + executed study notebook", done: false }
  # --- Gate ---
  - { name: "R  Autonomy-readiness checklist pass -> flip status draft->ready", done: false }
blocked_reason: ""
---

# Wide low-spread wavepacket — isolating purely quantum stopping (localised slab)

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository. σ always means the wavepacket width σ_WP.
</identity>

<rough_draft_banner>
This is a ROUGH DRAFT. Its job is to capture the user's intent and the simple
questions already resolved (2026-06-27), so a later `/campaigns` session can turn
it into an autonomy-ready prompt. Do NOT execute it as-is. Sections marked
"(open — finalised in the /campaigns pass)" are intent, not locked decisions.
</rough_draft_banner>

<prior_art>
This campaign already existed in sketch form — this file consolidates it; it does
not invent it:
- `docs/campaigns/jellium_wp_stopping/draft_campaigns.md` → **"Campaign 2 — large
  rigid σ"** (operating point, caveats).
- Task **P2.3** in `docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md`
  ("large-σ low-spread WP + matched classical … isolates quantum-vs-classical
  stopping from dispersion") — **retained there and cross-linked here** (user
  chose a standalone file, not retiring P2.3).
- Supporting analysis: `docs/campaigns/jellium_wp_stopping/brainstorming-jellium-campaigns.ipynb`;
  Campaign-1 restrictions: `notes_campaign1_sigma05_restrictions.md`.
- Literature anchor: `docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`
  (quantum vs classical projectile gives different stopping; classical
  point-charge limit recovers Lindhard; the projectile **width** matters).
- **Related but separate dataset:** the executed σ_WP=3 *no-CAP* twins in
  `docs/campaigns/quantum_classical_nocap/` are on **bulk** cubic jellium
  (50³, N=162), NOT this localised slab — do not conflate.
</prior_art>

<description>
**Why this campaign exists (user's words, lightly edited).** The aim is a
simulation system capable of identifying the artefacts that arise due to **purely
quantum effects** that are missing in classical projectiles. To do this we select
a wavepacket and a classical projectile with the **same Gaussian width**, and we
ensure that width is chosen so the **wavepacket does not change its width
appreciably** over the transit. That makes it convenient to isolate the purely
quantum effects.

**The problem with the existing σ=0.5 run.** In the quantum-stopping-power run the
WP width is σ_WP=0.5 Bohr, so by the time it reaches the jellium slab it has
**expanded appreciably**. Its comparison is against a classical projectile whose
radial-potential σ does **not** change. So (a) the Coulomb interaction between the
two differs, and (b) the WP's effective interaction radius grows relative to the
fixed classical particle — both can contaminate the comparison.

**The fix (spreading law).** In free propagation the final σ depends only on the
initial σ and the elapsed time:
`σ(t) = σ₀ · √(1 + (ħ t / (2 m σ₀²))²)`.
Choose σ₀ and the total time so the packet stays near-rigid. This generally points
to a **large σ**, which means the cell sizes, the WP energy (hence total
simulation time), and other parameters must be adjusted accordingly. Then the
interaction range of the quantum and classical projectiles is the same (or on a
similar scale), and any remaining WP−classical difference is a purely quantum
effect.

**Decision this informs.** Whether a clean, dispersion-free quantum stopping power
S(E) can be measured on the localised slab and compared, like-for-like, to the
matched classical projectile — the central quantum-vs-classical result of the
jellium stopping programme.
</description>

<workflow>
**Two-phase (user decision, 2026-06-27).** The campaign mixes a human-in-the-loop
design phase with an autonomous production phase — capturing the user's "one
experiment at a time" preference (skill-change request #1/#2):

- **Phase 0 — iterative design (human-in-the-loop).** Start from a single system
  setup, run it, the user reviews the results and gives feedback, then the design
  is improved considering that feedback and retried. This iteration continues
  until the system is believed to work — i.e. **the quantum stopping power can be
  calculated clearly from the data**. The autonomous phase does NOT start until
  the user signs off the design.
- **Phase 1 — autonomous sweep.** Once the design works, run a sweep over the WP
  energy (with varying σ sizes) and produce the **S(E)** plot for the quantum case
  plus the matched classical, end-to-end with no user in the loop.
</workflow>

<observables_set>
**LOCKED (2026-06-30): the EXTENSIVE suite used by the localised-jellium phase-5
runs** (verbatim from `scripts/qsp_phase5/wp/run.cpp`), at the phase-5 cadence
(`WRITE_EVERY=4`; state energies/occupations every 5·WRITE_EVERY; WF every
`WF_EVERY`). This is the full raw+processed set the user is happy with.

**Scalars — `raw/observables/observables.csv`:** `step`, `time_au`,
`energy_total`, `energy_kinetic`, `energy_hartree`, `energy_xc`,
`current_x/y/z`, `dipole_x/y/z`, `density_l2`.

**Per-state / electronic structure:** `state_energies.csv` (per-orbital energies
vs t, `StateEnergyWriter`); `occupations_vs_time.csv` (`OccupationsWriter`);
`eigenvalues/` (copied from the GS checkpoint).

**WP-resolved (the quantum-stopping observables):**
- `momentum_distribution.csv` — **n(k,t)**, 64 bins → stopping via the
  **coherent-peak centroid** (sub-bin fit), NOT ½⟨k²⟩ (scattering-inflated).
- `wp_momentum_stats.csv` — mean k, spread (`WPMomentumStats`).
- `wp_real_space_stats.csv` — WP centroid z(t) **and width σ(t)** (`WPRealSpaceStats`)
  → **directly measures the spreading**, so the run validates the §1 analytical
  prediction (2.6 % at the slab).
- `overlap/`, `overlap_full/` — orbital-overlap (autocorrelation) matrices.
- `electron_number.csv` — **N_total(t)** (bath / CAP-absorption guard).

**Density fields (VTI, physical order — load via `inqview.load_vti`, never
fftshift):** `density_total`, `density_system`, `density_gs_system`,
`density_wp`, `density_delta` + `density_delta_coarse` (L2, 3.0-Bohr coarse bin);
`wavefunction_wp` (complex field, every `WF_EVERY`).

**Provenance:** `run_summary.txt` (full config snapshot).

**Classical run differs:** drop the WP-specific block (momentum/overlap/
wavefunction_wp); **add `electron_track.csv`** (z, vz, F) → S = −dKE/dz. Keep the
scalars, per-state, density VTI, and the N_total guard.
</observables_set>

<resolved_decisions>
Locked in this rough-draft session (2026-06-27):
- **Structure:** standalone dedicated file (this one); consolidates the
  `draft_campaigns.md` "Campaign 2" sketch; P2.3 retained + cross-linked.
- **Workflow:** two-phase — iterative human-in-the-loop design → autonomous
  E(×σ) sweep (see <workflow>).
- **System:** the **localised jellium slab** (matches the σ=0.5 production run),
  two-sided sin² CAP; NOT the bulk σ=3 no-CAP dataset.
- **σ-convention:** σ = σ_WP everywhere; matched classical UPF at σ_pot=σ_WP/√2.

**LOCKED in the 2026-06-30 /campaigns deliberation** (planning notebook:
`wide_wavepacket_planning.ipynb` + `build_wide_wp_planning_notebook.py`, same
folder; figs in `wide_wp_planning_figs/`):
- **Operating point: σ_WP = 3.5 Bohr, E = 300 eV** (v = 4.70 a.u.). Spreading law
  in the planned **50×50×100 box** gives spread **2.6 % at the slab** (matched-σ
  essentially perfect) and **17.6 % at the far CAP**. σ₀=0.5 reference spreads
  ~1020 % at the slab — the production-run contamination this campaign fixes.
- **Box: 50×50×101 Bohr** (z∈[−50.5,+50.5]), **spacing dx = 0.40** (refined from
  phase-5's 0.50 to keep the cutoff well within 4σ_p up to E=600 — see Phase-1
  grid; clear of the dx=0.30 WP-init deadlock). **LZ=101 chosen so the
  equidistant launch gives EXACTLY 4σ clearance** (user decision 2026-06-30: CAP↔WP
  gap = 4σ). Launch z₀ = −26.5; clearance to BOTH the slab face and the CAP inner
  face = 14.0 Bohr = 4σ at σ=3.5 (σ=3.5 sits exactly on the 4σ boundary rule).
- **Slab unchanged:** |z|<12.5, N=82, r_s≈5.67 (matched to production).
- **CAP = same as the previous campaign (phase-5), user decision 2026-06-30:**
  two-sided sin² absorbing potential, **η = −0.7 Ha, 10 Bohr/side**. In the LZ=101
  box, region **[±40.5, ±50.5]** (inner faces ±40.5). Real INQ (narrow packet) →
  ε≈0.2 % at these params; wide-σ adequacy verified by the Phase-0 gate (below).
- **Sim-time estimate:** wall ≈ 0.060 h/a.u. (phase-5's 0.054 ×1.11 for LZ 90→100);
  geometry-minimal τ≈30 a.u. → ~1.8 h/run at E=300.
- **CAP-completeness finding (REAL INQ vs toy, planning notebook §4 —
  `cap_reflectivity_real_vs_toy.png`).** Real INQ two-sided-CAP data
  (`systems/vacuum/hypotheses/twosided_cap_vs_mask/twosided_combined.csv`) has
  E=300 eV points. At phase-5's per-side width **Lhalf=10** (= dataset `L20`) and
  η=0.7 the CAP reflects only **~0.2 %** for the *narrow* benchmark packet
  (σ=4√2/k₀≈1.2); raising |η|→1.0 → ~0.03 %. The 1D `cap_toy` forward model is
  **~5× pessimistic** vs these real runs. **Caveat:** the real runs use the narrow
  k₀-tied packet, NOT our σ≈4 — the toy (upper bound) predicts higher wide-σ
  reflection, so **wide-σ CAP adequacy is the one open risk**, to be closed by a
  Phase-0 real wide-σ run.

Open — finalised in the next /campaigns prompts:
- **Phase-0 CAP-completeness gate:** confirm the phase-5 CAP (η=−0.7, 10/side)
  absorbs the *wide* σ≈4 packet (residual norm → 0). Fallbacks if it fails: raise
  |η|→1.0 (cheaper) or widen to Lhalf=15 (needs LZ≈110).
- **Spread tolerance** as a formal Phase-1 gate (the slab-spread is ≤~3 % across the
  intended E range at σ₀=3.5 — effectively already satisfied).
- **Matched classical UPF** at σ_pot = 3.5/√2 ≈ 2.47 (generate) + new slab GS for
  the LZ=100/110 box.
- ~~The E×σ grid~~ — **LOCKED below (<phase1>): one width σ=3.5, E∈{200,300,400,
  500,600} eV, dx=0.40.**
</resolved_decisions>

<phase0>
**PHASE 0 — iterative human-in-the-loop design (concrete, 2026-06-30).** The gate
before any autonomous sweep. Run one experiment, user reviews, refine, retry —
until quantum S extracts cleanly. Locked geometry/CAP/observables are above.

**P0.1 — Operating point — ✅ DONE.** σ_WP = 3.5 Bohr, E = 300 eV (v=4.70 a.u.),
box 50×50×101, CAP η=−0.7 / 10 Bohr-side / [±40.5,±50.5], launch z₀=−26.5 (4σ).
Analytical spread: 2.6 % at the slab, 17.6 % at the far CAP (planning notebook §1).

**P0.2 — Setup + GS + CAP-completeness smoke.**
- (a) **Matched classical UPF** at σ_pot = σ_WP/√2 = **2.475 Bohr** (electron-Gaussian
  UPF). Verify by the actual V(r) data, NOT the header/filename
  ([[reference_inq_ignores_is_coulomb_upf_flag]]).
- (b) **Slab GS for the LZ=101 box** (periodicity 3, slab |z|<12.5, N=82, dx=0.5).
  Validate: SCF converged, interior n₀≈1.31e-3 a₀⁻³, energy stable, interior
  matches the production slab (BC-/box-independent per GS study).
- (c) **Wide-σ CAP-completeness pilot** — short WP run (σ=3.5, E=300, locked CAP)
  propagated until the packet reaches/passes the far CAP. **GATE (numeric):**
  residual inner-region WP norm after absorption **< 2 %**, and N_total drift
  **< 2 %** (no bath drain). *If it fails:* raise |η| → 1.0 (cheap) or widen to
  Lhalf=15 (→ LZ≈111), re-test. *Done when:* both criteria pass.

**P0.3 — First matched WP + classical experiment.**
- (a) **WP run** (σ=3.5, E=300, full observable suite above) **+ matched classical
  run** (σ_pot=2.475 UPF, identical box/CAP/GS), launched concurrently on 2 GPUs.
- (b) **Extract quantum S** by the energy method `S=[E_total(t_f)−E_ref]/L_z` with
  `E_ref = E_total(0) − ⟨T_WP⟩ − E_SIE`, **cross-checked** against the
  momentum-centroid stopping (n(k,t) coherent peak). **Classical S** = −dKE/dz from
  `electron_track.csv`.
- (c) User reviews → refine → retry. *Done when:* the **user signs off** that
  quantum S extracts cleanly, i.e. all gates below hold.

**Phase-0 numeric gates (all must hold before Phase 1):**
- spread@slab ≤ ~3 % (predicted 2.6 %; read from `wp_real_space_stats.csv` σ(t));
- CAP residual WP norm < 2 %; N_total(t) drift < 2 %;
- E_total(t) plateaus before E_final is read (else result = upper bound);
- quantum S finite + stable; two channels (energy vs centroid) agree;
- abort on NaN / complex energy / GPU occupied by another user.

**Validated-block ladder (cheap → full):** UPF V(r) check → GS SCF smoke → CAP
wide-σ smoke (P0.2c) → single matched WP+classical (P0.3). Each rung gates the
next; the autonomous Phase-1 sweep starts only after P0.3 sign-off.
</phase0>

<phase1>
**PHASE 1 — autonomous E sweep at one width (LOCKED 2026-06-30).** Starts ONLY
after the P0.3 sign-off. Produces the S(E) curve for the wide WP + matched
classical, and the WP−classical "quantum component" with SIE bounded.

- **Width:** σ_WP = 3.5 (single width; the σ-axis sweep is a deferred follow-on).
- **Grid spacing dx = 0.40 Bohr** — set by the **cutoff/aliasing guard** (mandatory):
  `k_max = π/dx ≥ k0 + 4σ_p`, σ_p = 1/(√2·σ_WP) = 0.202 a.u. Phase-5's dx=0.50
  fails this at E≥500 (1.1σ_p); dx=0.40 gives ≥6σ_p margin to 600 eV.
- **Energy grid (5 runs):** E ∈ **{200, 300, 400, 500, 600} eV**
  (v = 3.83–6.64 a.u.). Cutoff margin {20, 16, 12, 9, 6}σ_p — all well within 4σ.
  spread@slab {3.9, 2.6, 2.0, 1.6, 1.3}% — all rigid. (All high-v tail, v/v_F≈11–20;
  the Bragg peak at ~2–6 eV is unreachable by a rigid packet — expect a
  monotonically falling S(E).)
- **Per E:** WP run + matched classical run (σ_pot=2.475 UPF) → **10 production
  runs**. Plus **1 vacuum-WP SIE control** (bounds the quantum component). σ_WP=0.5
  phase-5 S(E) reused as the spreading-contaminated reference; Lindhard
  (point + finite-σ) analytical.
- **Cost:** wall ≈ 0.118 h/a.u. (phase-5 0.054 × grid factor 2.19 for LZ 101 + dx
  0.40); ~3–4 h/run → ~35 h compute, ~17 h on 2 GPUs.
- **Per-run S:** energy method `S=[E_total(t_f)−E_ref]/L_z`,
  `E_ref=E_total(0)−⟨T_WP⟩−E_SIE`, cross-checked vs n(k,t) coherent-peak centroid;
  classical S = −dKE/dz. Convergence gate + N-guard as Phase-0.
- **Output:** per-run REPORT + S(E) overlay (WP vs classical vs Lindhard vs the
  σ=0.5 reference) updated/emailed after each run; final executed study notebook
  with the WP−classical quantum component (SIE bounded).
</phase1>

<caveats>
(from the existing draft — bake these in)
- **"Matched σ" does not FULLY isolate quantum-ness.** The WP carries Pauli +
  self-interaction error (SIE); the classical ghost carries neither. Quantify and
  bound the SIE with a per-σ **vacuum-WP control**; report any WP−classical
  "quantum component" only with the SIE flagged + bounded.
- **CAP completeness.** The CAP was tuned for v≈2.7; a faster rigid packet absorbs
  less per length, so pilot-check absorption completeness (residual norm → 0) and
  raise η if needed.
- **PROVISIONAL** w.r.t. the inq-study engine regression (cf. the nocap campaign's
  Task #7) until that is cleared.
</caveats>

<guard_rails>
(open — to be locked via /campaigns)
- **Spreading check:** σ(t_final)/σ₀ ≤ the chosen tolerance (the whole point).
- **CAP completeness:** residual WP norm → 0 over the run; else raise η.
- **N(t) ≈ const** outside the absorber window (bath conservation).
- **Wrap truncation:** truncate before the packet/projectile re-enters the box.
- Boundary 4σ/1σ launch-stop rule + 300-frame VTI cadence (always-on rules).
- SIE bounded via vacuum-WP control before any "quantum component" is reported.
- Abort on NaN / complex energy / GPU occupied by another user.
- Phase 1 does not start until Phase 0 is signed off by the user.
</guard_rails>

<tasks>
(rough — done-criteria sharpened via /campaigns)
**Phase 0 — iterative design (human-in-the-loop)**
1. **P0.1** Pick σ₀ and transit time from the spreading law for the chosen spread
   tolerance; land an operating point (start from σ≈4, E≥300 eV). *Done when:*
   σ₀/E/time chosen with the spreading-law calc shown.
2. **P0.2** Size the cell/energy/time for that σ; generate the matched σ UPF if
   absent; pilot-check CAP completeness. *Done when:* a smoke run exits clean,
   spread ≤ tolerance, CAP residual → 0.
3. **P0.3** Run ONE WP+classical matched experiment; user reviews; refine; retry.
   *Done when:* the user signs off that quantum S extracts cleanly.
**Phase 1 — autonomous sweep**
4. **P1.1** Run the E sweep (over varying σ) → S(E) for WP + matched classical.
   *Done when:* S(E) curves exist for the grid.
5. **P1.2** Build the S(E) overlay + WP−classical "quantum component" with SIE
   bounded; executed study notebook. *Done when:* notebook + verdict recorded.
</tasks>

<rules>
- ALWAYS keep `inq/` immutable; engine work (if any) goes in `inq-study`.
- σ always means σ_WP; surface σ_pot=σ_WP/√2 only in a methods footnote.
- Report S numbers at 2 s.f. (3 s.f. only for genuine near-equalities).
- NEVER report a WP−classical "quantum component" without the SIE bounded.
- Phase 1 (autonomous) starts ONLY after Phase 0 user sign-off.
</rules>

<preflight>
(rough draft — NOT yet autonomy-ready; reminder of what /campaigns must satisfy)
- [ ] Spread tolerance fixed; σ₀/E/time derived from the spreading law with values.
- [ ] Box/cell sized; matched σ UPF generated; CAP completeness pilot passed.
- [ ] E×σ sweep grid enumerated; observables per run + cadence locked.
- [ ] SIE-bounding vacuum-WP control defined; coherent-peak S estimator fixed.
- [ ] Phase-0 → Phase-1 sign-off gate explicit; notebook contract + handover set.
- [ ] Grounding: spreading law + Nazarov–Gross premise cited; engine claims with
      source line-refs.
</preflight>
