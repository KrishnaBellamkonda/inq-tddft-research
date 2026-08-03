# Plan — Cylindrical proximity ladder: weak → strong coupling KS stopping

Status: **DESIGN LOCKED 2026-08-02**, implementation starting.
Predecessor: `docs/plans/cylindrical-channeling-ks-stopping.md` (the R_in = 10 rung,
complete — its twin IS rung 1 of this ladder and is NOT re-run).
Handover: `docs/handovers/cylindrical-proximity-ladder.md` (to be created).

---

## 1. Aim

The channeling twin (R_in = 2.5 σ_WP) validated a KS-orbital definition of stopping
power against the classical ΔE/Δs definition **in the weak-coupling limit**: the
projectile flew down an almost-empty bore and the T₁ estimator matched the classical
twin to ~20 %. This campaign walks that same twin from grazing to fully immersed by
bringing the tube wall inward, and asks **where — and how — the classical/quantum
agreement breaks**.

Every rung is a matched pair: a wavepacket (σ_WP = 4 Bohr, occupied KS orbital) and a
classical Gaussian projectile at the matched σ_pot = σ_WP/√2 = 2.828 Bohr
(`.claude/rules/sigma-wp-convention.md`). Projectile energy is held at 50 eV
(v₀ = 1.917, v/v_F = 3.00) and the jellium density at r_s = 3.000 throughout.

## 2. What this parametrisation does and does not reach

**Reaches.** A ~10× sweep in the bath density the projectile samples, taking
fractional energy loss over the run from ~13 % to an estimated ~66 %. That is
perturbative → strongly non-linear, which is the stated goal. The projectile stays
fast throughout (v/v_F falls 3.0 → ~2.0), so it never crosses the Bragg peak and the
regime stays one-sided and interpretable.

**Does not reach.** The textbook "full stopping power" regime. The projectile's
Gaussian form factor is exp(−q²σ_pot²/2), which is 0.37 at q = 0.5, 0.018 at q = 1,
and 3×10⁻²⁶ at q = 2v₀ = 3.83. The plasmon pole sits at q_min = ω_p/v = 0.174 and the
electron–hole continuum runs to q = 2v. **This projectile couples to the collective
response and essentially nothing else**, at every rung. Moving the wall in scales how
much medium responds; it does not harden the projectile.

> The ladder is therefore **weak-collective → strong-collective**. Reaching the pair
> channel is a *σ_WP* axis, not an *R_in* axis. The vacuum σ-sweep
> (`systems/vacuum/hypotheses/wp_selfinteraction/sigma_sweep.py`) already maps the
> self-interaction cost of moving along it. Out of scope here; recorded so the
> campaign's conclusions are not over-claimed.

## 3. The ladder

Fixed: r_s = 3.000000, L = 40×40×60, spacing 0.5 Bohr, edge_width 0.5, L_z tube axis.
N_e is chosen EVEN and R_out is then solved so n₀ = N/V is exactly the r_s = 3 value
(the ∫n₊ = N exactness the G=0 cancellation needs). R_out moves by ≤ 0.014 Bohr, far
inside the 0.5 Bohr erfc edge — physically the same wall.

| rung | R_in | R_in/σ_WP | N_e | R_out | n_states | shape |
|---|---|---|---|---|---|---|
| 1 (done) | 10 | 2.50 | 160 | 14.000 | 104 | annulus |
| 2 | 8 | 2.00 | 220 | 14.000 | 143 | annulus |
| 3 | 6 | 1.50 | 266 | 13.986 | 172 | annulus |
| 4 | 4 | 1.00 | 300 | 14.000 | 195 | annulus |
| 5 | 0 | filled | 326 | 13.986 | 211 | **cylinder** |

`n_states = N/2 + extra`, extra ≈ 30 % of occupied (matching rung 1's 24/80). The
smearing is 0.00862 eV, which is cold for the dense subband structure of the filled
nanowire — **if SCF struggles at rungs 4–5, raise `extra_states` before raising T.**

### 3.1 The rungs are exponentially spaced and they merge in time

WP charge inside the wall region, exp(−R_in²/2σ_d(t)²) with σ_d(t) the free spreading:

| R_in | t = 0 | t = 13 | t = 30 |
|---|---|---|---|
| 10 | 0.19 % | 2.3 % | 25 % |
| 8 | 1.8 % | 9.0 % | 41 % |
| 6 | 10.5 % | 26 % | 61 % |
| 4 | 37 % | 55 % | 80 % |
| 0 | 100 % | 100 % | 100 % |

Labels "2.5σ, 2σ, 1.5σ, 1σ" read linear; the coupling spans a factor of 190 at t = 0
and only 3 by t = 30. **The rungs are distinct only early** — which is where S is
fitted, so this is workable, but it fixes the analysis protocol (§5).

### 3.2 "Distance to the jellium" is not "distance to the electrons"

Rung 1's own ground state: `n_bore_electrons = 16.02` of 160, `nbar_bore = 9.7 % of n₀`,
`bore_depletion_ratio = 0.129`. The bore is not empty — that spill-in density is *why*
rung 1 registers any stopping at all. The campaign's x-axis is therefore the **measured**
mean bath density sampled by the packet, not the nominal R_in (§4.1).

## 4. Implementation

### 4.1 New / changed code

- **`inqkit::jellium::background_shape::cylinder`** — DONE. A filled tube is its own
  shape, not `annulus` with R_in = 0: the erfc step is centred on its nominal edge, so
  `background_mask(0,0,w) = ½` would put n₊ = n₀/2 exactly on the tube axis, relaxing
  to n₀ only by d ≈ 2w — silent, and maximal precisely where the projectile flies.
  `cylinder_mask()` has one erfc factor and no degenerate edge.
  Tests T0.7/T0.8/T0.9 added (§6). Purely additive: nothing else switches on the enum
  and no existing run passes `inner_radius = 0`.
- **`⟨n_bath⟩_WP` observable** — TODO. `∫n_wp·n_bath d³r / ∫n_wp d³r`, written every
  step alongside `wp_radial_occupancy`. This is the campaign's coupling coordinate and
  the only variable that is well defined at *every* rung including the filled one
  (`f_bore` is meaningless at R_in = 0). Two field products per step — negligible.
- **Config header** `shared/configs/proximity_ladder_rs3.hpp` — TODO. One struct per
  rung, same derivation-not-assumption style as `channeling_tube_rs3.hpp`.

### 4.2 Reused unchanged

Both `run.cpp` (WP and classical) from `scripts/channeling_twin/`, parameterised by
rung. Interaction energies (`.claude/rules/decomposed-interaction-energies.md`),
final-timestep checkpointing, momentum stats, real-space stats, density frames — all
already wired.

### 4.3 The same-N control (one extra pair)

The ladder moves three things together: proximity, N_e (160 → 326), and the target's
mode spectrum (thin annulus with two coupled surfaces → solid nanowire). In this
geometry they are inseparable, because the electrons added *are* the close ones.

Control: **R_in = 4 with R_out = 10.583**, giving the same annulus volume and hence
N = 160 and n₀ unchanged. Compare its S against rung 4 (R_in = 4, R_out = 14, N = 300).
Agreement ⇒ the ladder is a proximity sweep and the far material is inert. Disagreement
⇒ S must be reported per-electron or against ⟨n⟩_WP, not against R_in. One pair, ~2 GPU-h,
and it decides how every other rung may be read.

## 5. Analysis protocol

**Per rung — one notebook** (`hypotheses/proximity_ladder/rung_<R_in>/`), containing,
per the user's list:

1. interaction energies, classical and WP (E_SS, E_PS, E_PP, E_SB, E_PB)
2. projectile position / trajectory
3. T₁ = ⟨p⟩²/2m, T₂ = ⟨p²⟩/2m, and the var(p) term T₂ − T₁
4. classical ½mv²
5. classical vs WP energy-loss definitions overlaid
6. WP momentum-loss distribution at several times
7. S from T₁, T₂ and the classical definition, with uncertainties
8. density-matrix GIF at the top (`.claude/rules/notebook-density-gif.md`)

Reuse `hypotheses/channeling_twin/build_report_figures.py` — the panel set already
covers 1–8; parameterise it by rung rather than rewriting.

**Fit window — the one protocol change.** Do NOT fit S over a fixed time window across
rungs; by §3.1 that compares different couplings. Define the window by a **common band
in ⟨n_bath⟩_WP**, plus the light-projectile constraint v ≥ 0.85 v₀
(`.claude/rules/light-projectile-stopping.md`). Report the window's mean v and point
count with every S, as rung 1 already does.

**Cross-rung comparison** (`hypotheses/proximity_ladder/`): S_T1, S_T2, S_classical vs
⟨n_bath⟩_WP (primary) and vs R_in/σ_WP (secondary, for readability). The headline
result is the **ratio** S_WP/S_classical as a function of coupling — where the KS
definition stops tracking the classical one, and which term (var(p), E_PP, xc coupling)
accounts for the divergence.

### 5.1 Two asymmetries that grow along the ladder — measure, do not assume

- **Pauli blocking is negligible, and it is worth logging the proof.** The WP is
  Gram–Schmidt orthogonalised against every occupied bath orbital
  (`channeling_twin/wp/run.cpp:288`). At 37 % or 100 % spatial overlap that sounds
  destructive, but the overlap is set in *momentum* space: k₀ = 1.917 sits 7.2 σ_p above
  k_F = 0.640, so ⟨φ_bath|ψ_WP⟩ ~ 10⁻¹¹. Record `report.max_overlap` per rung; it should
  stay ≲ 10⁻¹⁰. It would only start to matter if a rung decelerated below v ≈ 1.2.
- **xc coupling is the asymmetry that does grow.** The WP adds ~32 % to the on-axis
  density (peak n_wp = 2.81×10⁻³ vs n₀ = 8.84×10⁻³) and LDA xc is non-linear in n, so
  E_xc[n_bath + n_wp] ≠ E_xc[n_bath] + E_xc[n_wp]. The classical twin contributes
  nothing to xc at all. Diagnostic: E_xc[n_tot] − E_xc[n_bath] − E_xc[n_wp]. Part of
  the WP−classical gap at the strong rungs is this, not dispersion.
- **Self-interaction is bounded and shrinking.** E_PP is a property of the packet alone,
  so its absolute contribution is roughly rung-independent while the real stopping grows
  ~10×; the SIE *fraction* should fall along the ladder. Rung 1 measured it at 20.9 % of
  the excess spreading (SIC-PZ twin, agreeing with the vacuum prediction to 0.4 %).
  Re-run SIC on the **filled rung only** to bound the other end — 2 points bracket the ladder.

## 6. Validation

| ID | Tier | Check | Status |
|---|---|---|---|
| T0.7 | pure | filled cylinder carries n₀ **on the axis** at w > 0 (the n₀/2 trap) | written |
| T0.8 | pure | hollow tube still carves its bore at w > 0 (no regression) | written |
| T0.9 | engine | `cylinder` builder ∫n₊ = n₀πR²L_z, and ignores `inner_radius` | written |
| — | engine | full inqkit suite green after the enum addition | pending GPU |
| G1 | per-rung GS | ∫n₊ = N to 1e-9; r_s = 3.000; SCF converged | pending |
| G2 | per-rung smoke | cutoff guard PASS; t=0 injection gates pass; `max_overlap` ≲ 1e-10 | pending |
| G3 | per-rung run | interaction-energy closure vs INQ scalars (`decomposed-interaction-energies`) | pending |
| G4 | per-rung run | clean initial-drag S exists (finite, ≥30 early-window points) | pending |

T0.7/T0.8 arithmetic pre-verified on the host against the production functions
(scratchpad `mask_check.cpp`, all pass; the old composition returns 0.5000 on the axis
where 1.0 is required). Engine tier needs an A100 — login node has CUDA 11.4 and cannot
build the suite.

## 7. Cost

Measured rung 1: WP 3089 s, classical 2339 s at 104 states. Scaling ∝ n_states:

| rung | n_states | WP | classical |
|---|---|---|---|
| 2 | 143 | 1.18 h | 0.89 h |
| 3 | 172 | 1.42 h | 1.07 h |
| 4 | 195 | 1.61 h | 1.22 h |
| 5 | 211 | 1.74 h | 1.32 h |

**RT total for 4 new rungs = 10.5 GPU-h.** Plus 4 ground states (est. 4–6 GPU-h,
rung 5 the worst), the same-N control pair (~2 GPU-h), and one SIC run at rung 5
(~1.7 GPU-h). **Campaign ≈ 19–21 GPU-h.**

## 8. Execution order

0. inqkit test suite on GPU — validates `cylinder` before anything depends on it.
1. 4 ground states (rungs 2–5). Gate G1.
2. Smoke runs, 20 steps, both twins, all rungs. Gate G2.
3. Production twins, rungs 2–5. Gates G3/G4.
4. Same-N control pair at R_in = 4; SIC at rung 5.
5. Per-rung notebooks, then the cross-rung comparison.

Steps 1–3 are checkpointed and resumable
(`.claude/rules/final-timestep-checkpoint.md`); a projected overrun WARNs and proceeds
rather than blocking (`.claude/rules/checkpoint-dont-block.md`).
