# Plan: Nazarov–Gross mass ladder — dense jellium slab, wide wavepacket

**Status:** DESIGNED, not yet implemented. Branch `quantum-stopping-power`, device CSD3.
**Handover:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/handovers/nazarov-gross-slab-mass-ladder.md` (to be created at first milestone)
**Source note:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/sources/nazarov-2025-quantum-projectile-stopping.md`
**Supersedes:** `docs/campaigns/nazarov_gross_comparison/nazarov-gross-comparison.md`
(2026-07-11 design: r_s ≈ 4 slab, σ_WP = 0.5, high-v null branch + slow pilots —
staged but never launched; its `scripts/nazarov_gross/` machinery is reusable).

---

## 1. Aim (one sentence)

Measure whether the electronic stopping power of a jellium slab depends on the
**mass** of a same-charge, same-velocity quantum projectile — the central claim of
Nazarov & Gross (arXiv:2510.26222, 2025) — using a **bath-deposit** definition of
S that means the same thing for a classical perturbation and for a Kohn–Sham
wavepacket of any mass.

## 2. What the theory actually claims (and what we can and cannot test)

Nazarov & Gross factorise the many-body wavefunction as Ψ(**R**, r, t) =
χ(**R**,t)Φ_**R**(r,t) (Eq. 7), reduce to mean-field TDSCF (Eqs. 17–23), and derive
the friction coefficient

    Q = ∫ (ê·∇)V₀⁽ᵉ⁾(r) ∂_ω Im χ₁⁽ᵉ⁾(r,r′,ω)|_{ω=0} (ê·∇′)V₀⁽ᵉ⁾(r′) dr dr′     (Eq. 41)

which is **identical to the classical result except that the point-charge Coulomb
potential is replaced by V₀⁽ᵉ⁾**, the potential of the projectile's ground-state
density n₀⁽ⁿ⁾. The projectile is *bound* in the well its own screening cloud digs
(Eq. 32), and its ground state solves −ħ²/2M ∇²χ₀ + V₀⁽ⁿ⁾χ₀ = E₀χ₀ (Eq. 31).

**Mass enters through exactly one door: the width of χ₀.** The authors say so
(Sec. VII): "a result of the differences in the wave packets' sizes of particles
with different masses." Fig. 3 (Z = +1): Q = 0.45 (M = ∞), 0.385 (p), 0.27 (μ⁺),
≈ 0 (e⁺) at r_s = 0.5.

### Scope honesty — three things this campaign is NOT

1. **Not their Q.** Q ≡ lim_{v→0} S(v)/v is a strict friction limit. We run at
   v = 1.4 v_F, below the Bragg peak but not in the linear regime (measured with
   `inqview.analysis.lindhard_elf`: Q is within 5 % of Q(0) only up to ~0.42 v_F
   at r_s = 4, ~0.65 v_F at r_s = 2.07; the Bragg peak sits at ~1.6 v_F at both
   densities). We test the *physical claim*, not reproduce the *coefficient*.
2. **Not their mechanism's equilibrium half.** Their width comes from a bound
   ground state; ours comes from the injected σ₀ plus dispersion during transit.
   Same chain (mass → width → coupling), different link.
3. **Not their heavy rungs.** See §4 — μ and p are reachable in principle here
   (unlike at high v) but are outside the user-chosen mass range; the classical
   twin is the M→∞ anchor.

**What it IS:** a direct, mass-resolved measurement of the deposit into the
electron liquid, plus — and this is the part the analytic theory structurally
cannot do — an independent **σ sweep at fixed mass** that separates the width
effect from the mass effect.

## 3. Definition of S (locked)

**Primary — bath deposit.** `S = dE_bath/ds` fitted over an early in-slab window
with v ≥ 0.85 v₀ (`.claude/rules/light-projectile-stopping.md`). Chosen over the
KS-orbital −dT₁/ds because:

- it is **representation-independent** — the identical observable for the
  classical perturbation and for every mass, so the ladder sits on one axis;
- it measures **dissipation**, not bookkeeping: energy stored in the packet's own
  spreading is internal and reversible and correctly never appears (in the
  channeling twin, 54.2 % of the WP's drift loss was re-absorbed by its own
  spreading and never reached the bath);
- the slab is transversely **translation-invariant**, so transverse spreading
  does not convert into a geometric loss the way a finite cylindrical bore does
  (in `channeling_twin`, f_bore fell 0.997 → 0.506 and the impulse ratio tracked
  it at r = +0.98 — a confound with the same sign as the effect being sought).

**Secondary — KS-orbital, free of extra cost.** T₁ = ⟨p⟩²/2M and T₂ = ⟨p²⟩/2M are
written by the same runs. `ΔE_bath − (−ΔT₁)` is the packet's internal-energy
uptake: a quantum diagnostic with no classical counterpart.

**Not the headline:** −dKE/ds is the energy-conservation cross-check only
(`.claude/skills/stopping-power-extraction`, user decision 2026-06-30).

## 4. The three numerical constraints, and how each was calibrated

### 4.1 Aliasing — `π/h ≥ M·v + 3/(2σ_WP)`

The packet carries k₀ = M·v with momentum std σ_k = 1/(2σ_WP); the grid holds
k_max = π/h. Calibrated against this repo's own accepted guard (`nazarov_gross`
campaign: M = 2.2 at v = 2.711, h = 0.35, σ_WP = 0.5): 5.964 + 3.000 = 8.964 vs
π/0.35 = 8.976 — reproduced exactly.

Coarsest allowed h at v = 1.075:

| σ_WP | 3/(2σ) | M=1 | M=1.2 | M=3 | M=5 |
|---|---|---|---|---|---|
| 0.5 | 3.00 | 0.77 | 0.73 | 0.49 | 0.36 |
| **4.0** | **0.375** | **2.15** | **1.84** | **0.87** | **0.55** |

At σ_WP = 4 every entry is far above the production grid h = 0.50. **With a wide
packet at below-peak velocity the projectile no longer sets the grid — the bath
does — and aliasing ceases to be a constraint.**

### 4.2 Timestep — `dt ≤ 0.08 · min(M, 1) · h²`

From `dt · k_max²/(2M) ≤ 0.395`, calibrated on two working runs, both of which
turn out to have been sitting exactly on the ceiling:

| run | h | M | dt used | ceiling |
|---|---|---|---|---|
| p3 (`fullsuite_wp`) | 0.50 | 1 | 0.02 | 0.0200 |
| `sigma1_masspair` | 0.50 | 2 | 0.04 | 0.0400 |

**The `min(M, 1)` is the part that is easy to get wrong** (and I got it wrong at
first draft). One `dt` advances **all 124 orbitals**, and the 103 bath states have
m = 1. A *heavy* projectile therefore buys nothing — its own kinetic operator is
gentler, but the bath's is not — while a *light* one tightens the ceiling for
everybody. `sigma1_masspair` could use dt = 0.04 at M = 2 only because it is a
different, lighter-state system; in this slab the bath pins dt at 0.02 for every
rung with M ≥ 1.

**Consequence — cost is FLAT for M ≥ 1 and scales as 1/M below it.** Transit time
L/v is mass-independent, so step count ∝ 1/(v·min(M,1)·h²):

| rung | dt | steps | ≈ wall @ 5 s/step |
|---|---|---|---|
| classical, M = 3, 1.2, 1.0 | 0.020 | 2560 | 3.6 h each |
| M = 0.5 | 0.010 | 5120 | 7.1 h |

Ladder ≈ 25 h, σ sweep ≈ 11 h, pilots ≈ 2.5 h, GS ≈ 2 h → **≈ 41 GPU-hours**.
The dt rule is an INFERENCE from two points; the binaries enforce it as a hard
refusal and the Phase-2 drift gate confirms it before the ladder runs.

### 4.3 Traversal — the wide packet is what makes below-peak possible

At v = 1.4 v_F, r_s = 2.5, 25 Bohr slab, M = 1 (linear-response Lindhard;
`stopping_power_sigma` is documented in this repo as over-suppressing, so these
are the optimistic end of a bracket whose other end is the point-charge column):

| σ_WP | S/S_point | S (eV/Bohr) | deposit (eV) | deposit/KE |
|---|---|---|---|---|
| 0.5 | 0.830 | 3.99 | 99.7 | 6.3 — stops |
| 1.0 | 0.548 | 2.63 | 65.8 | 4.2 — stops |
| 2.0 | 0.200 | 0.960 | 24.0 | 1.5 — stops |
| **3.0** | 0.066 | 0.317 | 7.9 | **0.50 — crosses** |
| **4.0** | 0.021 | 0.101 | 2.5 | **0.16 — crosses easily** |

Widening the packet cuts the coupling by 1–2 orders of magnitude. This is not a
trick: the weak coupling of a delocalised projectile *is* the effect being
measured, so the mechanism under test is also what makes the test feasible.

Density then supplies the kinetic-energy budget (KE ∝ v_F² ∝ 1/r_s² rises steeply
while the σ-suppressed deposit stays nearly flat), 15 Bohr slab:

| r_s | v = 1.4 v_F | deposit (eV) | KE(M=1) (eV) | ratio |
|---|---|---|---|---|
| 2.0 | 1.343 | 4.32 | 24.6 | 0.18 |
| **2.5** | **1.075** | **4.76** | **15.7** | **0.30** |
| 4.0 | 0.672 | 5.25 | 6.1 | 0.86 |
| 5.665 (existing slab) | 0.474 | 5.00 | 3.1 | 1.63 — stops |

**The existing r_s = 5.665 slab cannot do this run. The campaign requires a new,
denser ground state.**

## 5. Locked parameter table

| | value | rationale |
|---|---|---|
| Density | **r_s = 2.50** (n = 0.01526 a₀⁻³, v_F = 0.767, E_F = 8.0 eV, ω_p = 0.438, T_p = 14.3 a.u.) | Bragg peak ≈ 1.23 a.u.; dense enough that M = 1 survives, not so dense the electron count explodes |
| Slab | L_z = **15 Bohr**, face 30 × 30 | N = 206 e⁻ → 103 occupied + 20 extra + WP ≈ **124 states** |
| Box | 30 × 30 × 120 Bohr, **h = 0.50** | 60 × 60 × 240 = 864 k points; points × states ≈ p3 → ~5 s/step |
| CAP | two-sided, **η = −1.0 Ha**, region [±45, ±60] | η = −1.0 is the validated value (0.13 % reflection, no bath drain, `muon-mass-fork` CAP study 2026-07-11) |
| Launch | z = −25 | 20 Bohr = 5σ clearance from the CAP face; 17.5 Bohr run-up to the slab |
| Velocity | **v = 1.075 a.u. = 1.40 v_F**, E = 15.7 eV at M = 1 | below the Bragg peak with margin |
| σ_WP | **4.0**; classical twin at σ_pot = 2.828 | σ = 3 gives more signal but 30 % KE loss and wrap risk; σ = 4 gives 16 % loss and a 2.5 eV deposit ≈ 2500× the energy drift |
| Mass ladder | classical (M→∞), **M = 3, 1.2, 1.0, 0.5** | ordered cheapest-first; M = 0.25 excluded — it reaches σ ≈ 12 Bohr and would wrap a 30 Bohr cell |
| dt | 0.02·min(M,1) → 0.02 / 0.02 / 0.02 / 0.01 | **bath-limited** (§4.2): one dt advances all 124 orbitals and the bath is m = 1, so heavy rungs get no speed-up. The binaries refuse to start above the ceiling. |
| Theory | **LDA (Hartree + xc), no SIC** | user decision; residual bounded by the Phase-1 vacuum control (§6, step 6) |
| Disk cadences | density VTI ÷24, wavefunction ÷4, checkpoint ÷3 of N_STEPS | Measured on this grid: density VTI 6.9 MB, wavefunction VTI 13.8 MB, RT checkpoint **1.7 GB**. Interior checkpoints ROLL into one directory, so ÷3 buys 2 interior writes + the mandatory final one for 1.7 GB total, not 3×. ≈ 390 MB frames + 1.7 GB checkpoint per WP run; pilot checkpoints are pruned. Campaign total ≈ 24 GB against 220 GB free. |

Estimated cost: 55 Bohr of path ≈ 51 a.u. → 2560 steps (3.6 h) for every rung
with M ≥ 1, 5120 steps (7.1 h) for M = 0.5. **Ladder + σ sweep + pilots + GS
≈ 41 GPU-hours.**

## 6. Steps

### Phase 0 — Build (no GPU physics)

1. **Shared config header** `shared/configs/slab_n206_L30x30x120_rs2p5.hpp` with
   every constant of §5. All downstream code reads it, so no rung can drift.
2. **Ground state + validation.** Converge once, reuse everywhere. Gates: r_s
   recovered from ∫n = 2.50; E_F within a few % of 8.0 eV; flat interior density;
   no Gibbs ringing at the erfc-softened faces (EDGE_WIDTH ≥ h).
3. **Two `run.cpp` binaries.**
   - **WP half** — mass fork on `inverse_mass[0][wp_idx]`, re-applied on resume
     since `save`/`load` does not persist it.
   - **Classical half** — the M→∞ anchor as a **perturbation, NOT a UPF/ion**
     (user decision 2026-08-05): `inqkit::dynamics::moving_gaussian_projectile_perturbation(proj, SIGMA_POT)`
     tracking a live `inqkit::dynamics::Projectile` advanced by velocity-Verlet in
     the RT callback, exactly as `scripts/localised_jellium_dynamics/proj_dyn/run.cpp`
     does. σ_pot = 2.828427 is a runtime argument, so **no UPF generation is
     needed** and the σ sweep of step 10 costs nothing extra on the classical side.
     (The alternative UPF/ion route of `fullsuite_classical` would need one
     `.upf` per σ and is deliberately not used.)
   Both halves carry the two-sided CAP, final-timestep checkpoint/resume
   (`.claude/rules/final-timestep-checkpoint.md`), and `interactions.csv`
   (`.claude/rules/decomposed-interaction-energies.md`). Build once, drive by env.
4. **Pre-flight guards** per rung: §4.1 aliasing, §4.2 dt, and `check_twin.py`
   parity → `twin_manifest.json` (`projectile` field will read
   `gaussian_charge_perturbation` vs `wavepacket_orbital`).

### Phase 1 — Cheap controls (vacuum, ~1 h)

5. **Free-dispersion validation of the mass fork.** One non-interacting vacuum
   run per mass; σ(t) must match σ₀√(1+(t/2Mσ₀²)²). End-to-end proof that the
   fork changes only the one wavepacket. Reuses `vacuum/scripts/wp_selfinteraction`.
6. **LDA self-interaction control.** Same masses at LDA; excess width over the
   non-interacting reference measures the SIE contamination.
   **The question is whether that excess is mass-DEPENDENT**, since only a
   mass-dependent error contaminates a mass ratio.
   *Context:* the 24-run vacuum sweep shows the net LDA width error GROWS with σ
   (3.1 % at σ = 1 → 7.9 % at σ = 4 → 9.7 % at σ = 8) — a wide packet is not
   automatically safe. But free spreading ∝ t/(2Mσ₀) and self-repulsion
   acceleration ∝ 1/(Mσ²) both carry 1/M, so the fractional effect should largely
   cancel in a ratio (INFERENCE — this step measures it). It matters because
   d ln S/d ln σ ≈ 4.8 here, so a 7 % width error becomes ~34 % in S if it does
   not cancel.

### Phase 2 — Pilot (~6 h, ≈5 % of campaign cost)

7. **Three short runs** — classical, M = 1, M = 0.5 at σ_WP = 4, 600 steps
   (enough to clear the slab).
8. **Pilot gate.** Pass requires: energy drift < 1e-3 Ha; CAP reflection < 1 %;
   no transverse wrap; ΔE_bath resolved well above drift; M = 1 and M = 0.5
   deposits separating. Fail → revise σ_WP or density BEFORE the ladder.
   Per `.claude/rules/checkpoint-dont-block.md`, a *cost* overrun is a WARN and
   proceeds; only these correctness gates block.

### Phase 3 — Production (~1 day)

9. **Mass ladder** — classical + M ∈ {3, 1.2, 1.0, 0.5}, fixed v and σ₀, full
   length. Cheapest-first so an overrun sacrifices the least informative rung.
10. **σ_WP sweep at fixed mass** — M = 1 at σ_WP ∈ {2, 3, 4, 6}. Decouples width
    from mass; the strongest novel result in the campaign. Costs nothing extra on
    the classical side now that σ_pot is a runtime argument (step 3).
11. **Thickness check** — one rung repeated at L_z = 25 Bohr, separating the bulk
    deposit (∝ L_z) from the fixed entrance/exit surface term.

### Phase 4 — Analysis

12. **Primary S extraction** — `S = dE_bath/ds` over the early in-slab window
    (v ≥ 0.85 v₀; widen to 0.70/0.50 if sparse). KE channel reported only as the
    conservation cross-check.
13. **Secondary KS-orbital channel** — T₁, T₂ and the gap against ΔE_bath (the
    packet's internal-energy uptake).
14. **Mechanism verification** — S vs *measured* mid-transit σ_eff rather than
    mass. If mass acts only through width, every rung collapses onto the step-10
    σ curve. This is the direct test of the Nazarov–Gross mechanism.
15. **Deliverables** — study notebook under
    `ResearchProject/systems/localised_jellium/hypotheses/ng_mass_ladder/` with
    the mandatory density-matrix GIF (`.claude/rules/notebook-density-gif.md`),
    report figures, run-catalogue rows (`tddft-run-catalogue`),
    `docs/validation/test-catalogue.md` entries, and a rolling handover.

## 7. Known risks / open uncertainties

| # | risk | mitigation |
|---|---|---|
| 1 | Traversal estimates use the packet's INITIAL width; the packet spreads and S falls exponentially with σ, so the real deposit is smaller and traversal easier — by an amount not computable statically | Phase-2 pilot measures it directly |
| 2 | The light arm self-rescues (spreads → decouples → survives) but may then deposit almost nothing, leaving little to distinguish from zero | 2.5 eV deposit is ~2500× the energy drift; a near-zero light-arm deposit IS the NG result, provided the gate in step 9 confirms separation |
| 3 | `stopping_power_sigma` over-suppresses (repo-documented), so every deposit above is a lower bound; the true value lies in a wide bracket | All sizing quoted as a bracket; the pilot replaces it with a measurement |
| 4 | dt ceiling is an inference from two runs | Step 9 energy-drift gate; back off 2× on failure |
| 5 | LDA SIE may impart a mass-dependent width error with the same sign as the effect | Step 7 measures it; if mass-dependent, either re-run with SIC-PZ (`inqkit/wavepacket/self_interaction_correction.hpp`, validated to width-ratio exactly 1.0 at σ = 1…8) or quote it as a stated systematic |
| 6 | CAP contaminates `E_total` once absorption starts (`muon-mass-fork` finding) | Analysis window must close before the packet's tail reaches the CAP; the 120 Bohr box and z = −25 launch give 5σ clearance |
| 7 | Referee objection: "just a Gaussian form factor, reproducible classically with a wider cloud" | The width is not a free parameter — it is fixed by M through ħ, and the classical twin at identical σ₀ shows no mass dependence. Step 10 + step 14 make this the explicit argument |
| 8 | **The σ = 6 sweep point is containment-marginal by construction**: 4σ = 24 Bohr inside a 30 Bohr transverse cell at t = 0, so its periodic images overlap from the start. σ = 2 also reaches 4σ ≈ 25 Bohr by exit (it spreads ×3) | Both are FLAGGED automatically — `pilot_gate` and every run notebook print `4·σ_perp` against the 30 Bohr cell. Treat σ = 6 as an indicative endpoint, not a quantitative one; if the collapse test hinges on it, re-run that point in a 40 Bohr cell (N = 275, ~160 states) rather than trusting it |

## 8. Defect found while designing (not fixed here)

`ResearchProject/systems/localised_jellium/hypotheses/classical_slab_stopping/analyse_classical_baseline.py:98-101`
calls `stopping_power_point(V0, RS)` and `stopping_power_sigma(V0, RS, …)`, but
the signature is `(v, kF, …)` and `RS` (line 34) is r_s, not k_F. At r_s = 4,
v = 1.3 this makes the Lindhard reference **13.17 eV/Bohr instead of 1.95 —
6.75× too high**. Should be `L.kF_from_rs(RS)`. Flagged, outside this campaign's
scope.

## 9. Cross-references

- Paper: arXiv:2510.26222 — note `docs/sources/nazarov-2025-quantum-projectile-stopping.md`
- Superseded design: `docs/campaigns/nazarov_gross_comparison/nazarov-gross-comparison.md`
- Competing geometry (rejected, §3): `docs/handovers/cylindrical-channeling-ks-stopping.md`
- Vacuum SIE sweep (steps 6–7 machinery): `ResearchProject/systems/vacuum/scripts/wp_selfinteraction/`
- Mass fork provenance: `docs/handovers/muon-mass-fork.md`
- Reference run.cpp to clone — WP: `scripts/nazarov_gross/wp/run.cpp`,
  `scripts/sigma1_masspair/wp/run.cpp`; classical (perturbation route):
  `scripts/localised_jellium_dynamics/proj_dyn/run.cpp`
- Classical projectile header: `inq-stack/include/inqkit/dynamics/moving_gaussian_projectile_perturbation.hpp`
- Rules in force: `light-projectile-stopping`, `final-timestep-checkpoint`,
  `decomposed-interaction-energies`, `sigma-wp-convention`, `notebook-density-gif`,
  `checkpoint-dont-block`, `validation-gates`
