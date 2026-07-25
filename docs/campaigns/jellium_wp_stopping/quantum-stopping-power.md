---
id: jwps-quantum-stopping-power
area: jellium_wp_stopping
title: "Quantum (wavepacket) vs classical stopping power — localised jellium slab"
status: paused
hypothesis: "In a localised jellium slab (r_s≈5.67), the electronic stopping of a quantum electron wavepacket (σ_WP=0.5, 100 eV) — from the converged E_total energy balance — differs from the matched classical projectile (Ehrenfest ΔKE_ion); Phase 1 first fixes the GS and quantifies the wavepacket self-interaction (SIE)."
handover: docs/handovers/localised-jellium.md
tasks:
  - { name: "P1.1 GS of localised jellium slab (82 e, r_s≈5.67, |z|<12.5, box 50×50×70) — converged + validated", done: true }
  - { name: "P1.2 WP+slab short RT, WP launched far (z≈−32, CAP off) — record E_total(0) + KE_WP(⟨p²⟩/2)", done: true }
  - { name: "P1.3 SIE estimate BOTH ways; cross-check SIE_a−SIE_b ≈ zero-point KE 3/(4σ²)=81.6 eV", done: true }
  - { name: "P2.1 WP+classical convergence/CAP test (40 a.u.) — DONE: CAP good (refl 1%), WP NOT converged (slope -0.9 eV/au), classical KE conservative dip-recovery (S_face=0.507 eV/Bohr), spreading ×41", done: true }
  - { name: "P2.2 WP longer-τ run — extend sim time until E_total equilibrates (slope→0); estimate from P2.1 absorption tail; PRECONDITION for any WP stopping number", done: false }
  - { name: "P2.3 large-σ low-spread WP + matched classical — σ chosen so the packet does not appreciably disperse over transit; isolates quantum-vs-classical stopping from dispersion", done: false }
  - { name: "P2.4 classical with projectile REMOVAL at box edge (z≥+35) — zero the ion's Gaussian charge (z_valence=0 ⇒ neutrality-safe) or park; same τ as WP; makes runs comparable by stopping periodic re-entry. Probe ion-mutation path first", done: false }
  - { name: "P3.1 BIG-BOX production (σ=0.5 WP + matched classical): box 50×50×90, TWO-SIDED CAP [±35,±45] η=−0.7 10 Bohr/side, equidistant launch z=−23.75, k₀=2.711 (100 eV), τ=100 a.u. (dt 0.04 pending smoke), energy-method stopping S=[E_total(t_f)−E_GS]/L_z with E_GS=−70.22568 Ha (GS=shared_gs/slab_n82_L50x50x90). Supersedes P2.2 (fixes init-absorption + adds longer τ + corrected energy bookkeeping)", done: false }
blocked_reason: "P3.1 EXECUTING (2026-06-25): big-box pair built against inq-study; dt=0.04 smoke running, production launches on confirmation. P2.2 folded into P3.1; P2.3 (large-σ) + P2.4 (classical removal) still deferred."
---

# Quantum (wavepacket) vs classical stopping power — localised jellium slab

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository.
</identity>

<description>
**Question.** Does treating the projectile *quantum mechanically* (an electron
wavepacket) change the electronic stopping power vs a *classical* projectile of
the same charge, width and velocity, in a localised jellium slab? The classical
number comes from the Ehrenfest projectile's own kinetic-energy loss; the quantum
number must come from the **total-energy balance** (the WP is a quantum orbital
with no classical trajectory). That asymmetry is why the WP is the hard case — its
number needs the total energy to **converge** (full absorption), which is the
make-or-break test of the later phases.

**This campaign is authored in phases.** **Phase 1 (locked here)** is the
foundation: converge + validate the ground state at the chosen density, and
quantify the wavepacket **self-interaction energy (SIE)** — the spurious
one-electron self-repulsion LDA does not cancel, which would otherwise contaminate
any WP−classical comparison. **Phase 2+ (deferred)** — the real-time production
runs (CAP, box length, total time, steady-state test for both runs, the stopping
comparison) — will be specified *after* analysing Phase-1 results.

**Decision Phase 1 informs.** The SIE magnitude bounds how much of any future
WP−classical difference is artifact vs physics; the GS feasibility/cost informs
the Phase-2 run plan. Success = a converged closed-shell-ish GS + a self-consistent
SIE number (the two estimates differing by exactly the zero-point KE).
</description>

<observables_set>
Phase-1 minimal set (ADR-0006), all cheap:
- **P1.1 GS:** total energy + components, eigenvalues/occupations (closed-shell
  check), converged density VTI. No dynamics.
- **P1.2 short RT (few steps, CAP off):** `observables.csv` (energy_total at t=0),
  `wp_momentum_stats.csv` (⟨p²⟩, e_kin_ha at t=0 → KE_WP), `wp_real_space_stats`
  (confirm launch position/σ). Minimal IO; no VTI cadence needed.
No NEW observable/kernel is introduced in Phase 1 (nothing to pre-gate). The
Phase-2 stopping extraction will use the **`stopping-power-extraction` skill**
(localised-slab branch: ΔE_total/L_z with a convergence gate for the WP; Ehrenfest
ΔKE_ion slope-fit for the classical) — specified in Phase 2.
</observables_set>

<resolved_decisions>
**Density & geometry (LOCKED).**
- Density **r_s ≈ 5.67** — chosen to match the long-standing **n162-in-50³**
  jellium (r_s=5.69) so Phase-2 stopping is directly comparable to the existing
  S(v) convergence runs, and to gauge computational cost/feasibility.
- **82 background electrons** in the slab (even; nearest practical to the 81.2 the
  density implies; pick the nearest **closed-shell** count at GS validation if 82
  is open-shell — not a blocker). WP run = 83 e (background+WP); classical run = 82.
- Slab **|z| < 12.5 Bohr** (25 Bohr thick), unchanged from the baseline.
- Box **50 × 50 × 70 Bohr** — x,y=50 preserves the in-plane density (r_s 5.67);
  z extended to 70 purely as vacuum, giving room for the CAP **and** a far WP launch
  (slab volume is unchanged by the z-extension, so density is unchanged).
- slab volume = 2·12.5·50·50 = 62500 Bohr³; n = 82/62500 = 0.001312 ⇒ r_s = 5.665.
- **Region layout (z):** slab [−12.5,12.5] · free [±12.5,±25] (12.5 Bohr each) ·
  CAP [±25,±35] (10 Bohr each) — see Phase-2 CAP block.

**SIE diagnostic (LOCKED) — run BOTH references.** The WP carries
`KE = drift (½k₀²) + zero-point (3/4σ² = 81.6 eV at σ=0.5)`; "+100 eV" is the drift
only, so it omits the zero-point.
- `SIE_a = E_total(0)[WP far] − (E_GS_slab + 100 eV)`   (user reference; = SIE + zero-point)
- `SIE_b = E_total(0)[WP far] − E_GS_slab − KE_WP`       (KE_WP = ⟨p²⟩/2 measured at t=0; = SIE)
- **Cross-check:** `SIE_a − SIE_b` must ≈ zero-point KE `3/(4σ²) = 81.6 eV`. Report
  **SIE_b** as THE SIE. (At r_s=4 the old `p3_wp` gave SIE_b ≈ 4.5 eV; expect a
  similar weakly-density-dependent value.)

**GS-vs-RT (LOCKED).** The slab-alone energy is a **pure GS**. The WP+slab energy
needs **RT injection** — the WP is a *moving* Gaussian (carries k₀), not an
energy-minimizing eigenstate, so a GS solver would collapse it. Inject and read
**E_total at t=0** (run a few steps only, CAP off) — this also yields KE_WP.

**WP / classical parameters (LOCKED for Phase 1).** σ_WP = 0.5 (density std 0.354);
matched classical charge std 0.354 (`electron_gaussian_sigma0p35.upf`); nominal
100 eV drift (k₀ from ½k₀²·27.2114 = 100). WP launch z ≈ −32 (far; ≈19.5 Bohr from
the slab face −12.5; 3 Bohr from the −35 box edge; CAP off in Phase 1, so no
absorption at launch).

**File placement (ADR-0007).** New config header + run machinery under
`ResearchProject/systems/localised_jellium/scripts/<sweep>/`; GS checkpoint under
`shared_gs/slab_n82_L50x50x60_rs5p67`; Phase-1 analysis +
`brainstorming-jellium-campaigns`-linked notebook under
`hypotheses/<sweep>/`. Logs gitignored; provenance only.

**Phase 2 (DEFERRED — open).** CAP = **10 Bohr each** (regions [±25,±35], inner
edge ±25, **η = −0.7**) — with the 70-Bohr z-box this leaves a **12.5 Bohr free
region** each side, ample far-launch room (conflict resolved). Still to specify
after Phase-1 analysis: total sim time, **steady-state (energy) test for BOTH WP
and classical**, periodic-wrap avoidance, and the stopping comparison.

**Phase 3 — BIG-BOX PRODUCTION (LOCKED 2026-06-25; σ=0.5 WP + matched classical,
energy-method stopping).** Fixes the three P2.1 issues — (a) WP norm absorbed at the
70-box launch, (b) E_total NOT converged at τ=40 (13.6% WP unabsorbed), (c) the
retained-energy stopping definition.
- **Box 50×50×90** Bohr (z∈[−45,45]), spacing 0.50 (100×100×180 grid); slab |z|<12.5 unchanged.
- **Region layout (z):** slab [−12.5,12.5] · free [±12.5,±35] (22.5 Bohr each) · CAP [±35,±45].
- **CAP = TWO-SIDED** (the benchmarked "known devil"): `cap_lo+cap_hi`, **η=−0.7 Ha,
  10 Bohr/side**, region [±35,±45] (inner faces ±35). Same η + per-side width as the P2.1
  CAP ⇒ the ~1.3% reflection benchmark carries over. Engine `inq-study/absorbing.hpp` is
  PRISTINE (a seam-centred CAP variant was built + Python-validated, then REVERTED). `inq/`
  never touched.
- **Launch z = −23.75** — EQUIDISTANT: 11.25 Bohr (≈22.5σ) to BOTH the slab face (−12.5) and
  the CAP inner face (−35) ⇒ kills the P2.1 init-absorption.
- **Projectile:** σ_WP=0.5 (charge std 0.354), k₀=v₀=2.711 (100 eV drift); matched classical
  Gaussian-electron ion (`electron_gaussian_sigma0p35.upf`) at the same v.
- **dt = 0.04 a.u.** (pending the 150-step stability smoke; fall back to 0.02 if unstable),
  **τ = 100 a.u.** (2500 steps at dt=0.04), **ETRS** (CAP ⇒ non-Hermitian), LDA, built vs inq-study.
- **GS:** `shared_gs/slab_n82_L50x50x90`, **E_GS = −70.22568 Ha**. The −24.5 Ha shift vs the
  70-box GS is a charged-slab-in-PBC electrostatic constant (kinetic/xc unchanged) that
  CANCELS in E_total−E_GS; never compare absolute E_GS across box sizes.
- **Stopping (energy method, user-defined):** S = [E_total(t_f) − E_GS]/L_z, L_z = 25 Bohr.
  The jellium SYSTEM ≡ density REMAINING in the box; CAP-absorbed = transmitted/reflected WP +
  secondaries (ledgered as a DIAGNOSTIC, NOT added back). INQ `energy_total` is the clean real
  ⟨H₀⟩ (CAP does not contaminate it). t=0 sanity: `E_total(0) − ⟨T_WP⟩ − SIE ≈ E_GS` — subtract
  the WP KINETIC ⟨T_WP⟩ = ½k₀²+3/(4σ²) ≈ 6.675 Ha (NOT the eigenvalue ε_WP → double-counts SIE);
  SIE (4.40 eV) only at t=0, expect E_GS + few-eV cross-Hartree (not exact). **Convergence gate:
  WP norm < 0.02 AND E_total plateau** before any S is quoted. ⟨p_z⟩-loss (Method A) rejected
  (KS orbital ≠ physical WP).
- **New diagnostics (post-hoc):** CAP-absorbed-energy ledger (WP vs bath split), convergence
  triple (residual norm, |dE/dt|, plateau width), measured t=0 WP↔slab cross-term — all from
  the saved density VTIs + energy trajectory; P2.1's full diagnostic suite retained.
- **Files:** `scripts/qsp_phase3/{gs,wp,classical}/run.cpp`, `run_production.sh`;
  `shared/configs/slab_n82_L50x50x90.hpp`; schematic `docs/notes/qsp_bigbox_run_schematic.png`.
</resolved_decisions>

<guard_rails>
- **GS validation gate (P1.1):** SCF converged (energy + density tolerances);
  record whether 82 e is closed-shell (eigenvalue gap at the Fermi level); if badly
  open-shell, try the nearest closed-shell count and note it. Abort on NaN / complex
  energy / non-convergence.
- **Short-RT gate (P1.2):** E_total(0) and KE_WP stable over the first few steps
  (no immediate blow-up); WP norm = 1.000 at injection (not clipped by the box edge
  at z≈−32 — 3 Bohr = ~8σ_dens from the −35 edge, tail negligible).
- **SIE sanity (P1.3):** `|（SIE_a − SIE_b) − 81.6 eV|` small (≲ a few eV); if it
  is NOT, STOP — either KE_WP or the energy reference is mis-handled, surface it.
- **GPU is the default** — schedule via the `cudaMemGetInfo` probe (NVML broken);
  warn if a GPU is occupied by another user. Phase 1 is cheap (GS + few-step RT).
- **Phase-2 geometry (resolved):** the 10-Bohr CAP each side ([±25,±35]) in the
  50×50×70 box leaves a **12.5 Bohr free region** each side ([±12.5,±25]) — ample
  room to launch the WP far from the slab and outside the CAP.
- **PROVISIONAL** until Phase 2 defines the stopping runs; this file is Phase-1-only.
</guard_rails>

<tasks>
**P1.1 — GS of the localised slab.** Build the new config (82 e, r_s 5.67, |z|<12.5,
box 50×50×60), converge the GS, validate (energy components, eigenvalues/occupations,
density VTI), save the checkpoint. Done = converged GS checkpoint + a one-line
closed-shell note. Skills: `tddft-simulations`, `simulation-validation`.

**P1.2 — WP+slab short RT (WP far, CAP off).** Inject the σ_WP=0.5 / 100 eV WP at
z≈−27, propagate a few steps, record E_total(0) and KE_WP (⟨p²⟩/2 from
`wp_momentum_stats`). Done = E_total(0) + KE_WP recorded, WP norm=1.000 confirmed.
Skill: `tddft-simulations`.

**P1.3 — SIE estimate (both ways) + cross-check.** Compute SIE_a, SIE_b, and the
zero-point cross-check; report SIE_b. Done = both numbers + the cross-check in the
handover and a short notebook section (link to `brainstorming-jellium-campaigns`).
Skills: `notebook-making`, `literature-review` (SIE / zero-point grounding).

**P2 — RT production (DEFERRED).** Specify after Phase-1 analysis (CAP, box,
sim-time, steady-state for both runs, stopping via `stopping-power-extraction`).
</tasks>

<rules>
- ALWAYS keep x,y=50 (in-plane density) when sizing the box; only z varies.
- ALWAYS subtract the FULL measured KE_WP (drift + zero-point) for the physical SIE
  (SIE_b); the "+100 eV" reference (SIE_a) is kept ONLY as the zero-point cross-check.
- NEVER run Phase-2 production from this file — it is Phase-1-only until extended.
- NEVER edit `inq/` (engine immutable); new code → `inq-stack/`, engine experiments
  → `inq-study/`.
</rules>

<preflight>
- [ ] Intent self-contained: Phase-1 hypothesis (GS + SIE) with explicit
      done-criteria; SIE cross-check is the falsifiable gate.
- [ ] Setup reproducible: r_s 5.67 / 82 e / |z|<12.5 / box 50×50×70; GS = task P1.1
      (no pre-existing checkpoint at this density); σ_WP=0.5, 100 eV, launch z≈−32;
      CAP off in Phase 1; file placement per ADR-0007.
- [ ] New code pre-gated: NONE in Phase 1 (no new observable/kernel).
- [ ] Validation & guard rails: GS convergence gate; SIE cross-check ≈81.6 eV;
      abort on NaN/complex/non-convergence; Phase-2 CAP-geometry conflict recorded.
- [ ] Autonomous mechanics: GPU via cudaMemGetInfo (NVML broken; warn if occupied);
      handover = docs/handovers/localised-jellium.md; agent flips frontmatter
      done/status; Phase-1 notebook section linked to brainstorming notebook.
- [ ] Grounding: zero-point KE 3/(4σ²) and SIE cited (literature-review); density
      r_s from n=82/62500; engine GS-vs-RT claim grounded (moving WP is not a GS
      eigenstate).
</preflight>
