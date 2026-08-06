# Plan: bulk-jellium KS-orbital stopping power (classical + wavepacket twin)

**Created:** 2026-07-30
**Branch:** `quantum-stopping-power`
**Machine:** CSD3, `ampere` partition (A100), account `mphil-nikiforakis-skcb2-sl2-gpu`
**Status:** design locked with user (2026-07-30); implementation not started

---

## 1. Objective

Extract the electronic stopping power S = −dT/ds of a 100 eV electron traversing a
**fully periodic bulk jellium** bath, using **KS-orbital-dependent definitions** of
both the projectile kinetic energy T and its position s, and compare a **classical
point-like projectile** against a **quantum wavepacket** projectile that is identical
in every other physical parameter.

The scientific question is what the *quantum* representation of the projectile
changes: specifically, how much of the apparent stopping is drift-momentum loss
versus momentum-width broadening (angular scattering + localisation energy).

---

## 2. Locked configuration (user decisions, 2026-07-30)

| Quantity | Value | Source |
|---|---|---|
| Cell | 46 × 46 × 80 Bohr, orthorhombic, **periodic in x, y, z** | user (wider xy, shorter Lz) + §3 analysis |
| Bath | N = 218 electrons, 109 occupied spatial states | shell-closure search, §3.3 |
| r_s | 5.70 Bohr (n = 1.3283e-3 e/Bohr³) | matches legacy L50/N162 r_s = 5.69 |
| ħω_p | 3.46 eV; v_F = 0.337; E_F = 0.0568 Ha | free-electron-gas formulae at that n |
| Shell gap at N=218 | 0.24 eV | closed-shell ⇒ clean SCF |
| Grid spacing | dx = 0.40 Bohr → 115 × 115 × 200 = 2.65 M points | user (2026-07-30, "dx=0.4, no check") |
| Cutoff | 30.84 Ha (π²/2dx²) | same as all prior jellium runs |
| XC | LDA ground state, ALDA in real time | jellium convention |
| Temperature | 0.00862 eV (≈100 K) smearing | `Base` default; shell gap ≫ smearing |
| Projectile energy | 100 eV = 3.6749 Ha → k₀ = v = 2.7111 a.u. | user |
| σ_WP | **2.0 Bohr** (ψ-width: ψ ∝ exp(−r²/2σ²)) | user |
| σ_d(0) | 1.4142 Bohr (density std = σ_WP/√2) | |
| σ_pot (classical) | **1.4142 Bohr** = σ_WP/√2 | `.claude/rules/sigma-wp-convention.md` |
| Classical mass | m = m_e = 1.0 a.u., **free Ehrenfest** | user |
| Classical charge | −1 e (Gaussian electron UPF) | parity with WP |
| Launch z₀ | −32.0 Bohr (= −L_z/2 + 4σ_WP) | `boundary_rule.hpp::launch_z` |
| dt | 0.04 a.u. | user |
| N_STEPS | **646** (centroid → +38 Bohr = stop_z; 70 Bohr traversal; t = 25.84 a.u.) | `boundary_rule.hpp::n_steps_for` |
| Fit window | t ∈ [4.0, 19.0] a.u. = steps 100–475 (40.7 Bohr of path) | §3.2 |

v/v_F ≈ 8.0 — well above the Fermi velocity, plasmon excitation kinematically
allowed. Both runs share ONE ground state.

### Naming

- GS checkpoint: `ResearchProject/systems/jellium/checkpoints/gs_L46x46x80_orth_N218_dx0p40/`
- Runs: `ResearchProject/systems/jellium/bulk_ks_stopping/{wp,classical}/`
- Machinery: `ResearchProject/systems/jellium/scripts/bulk_ks_stopping/`
- Analysis: `ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/`
- Config header: `ResearchProject/systems/jellium/shared/configs/bulk_ks_stopping_L46x46x80.hpp`

(ADR 0007 + 2026-06-15 sweep-grouping amendment; jellium is grandfathered-flat for
existing runs, but this new sweep uses the canonical layout.)

---

## 3. Design derivation — why this box

### 3.1 Free-dispersion law

A Gaussian with ψ ∝ exp(−r²/2σ²) has density std

    σ_d(t) = sqrt(σ²/2 + t²/(2σ²))

verified on this engine by the vacuum sweep (`shared/bin/run-dispersion.slurm`,
`ResearchProject/systems/vacuum/scripts/wp_traversal_energy/`). At σ = 2 Bohr,
σ_d grows 1.41 → 7.2 Bohr over 19 a.u. of flight.

**Minimum achievable arrival width:** minimising σ_d(t) over σ gives σ_opt = √t and
σ_d,min = √t. For an 80-Bohr flight (t ≈ 26 a.u.) no Gaussian can arrive narrower
than ≈5.1 Bohr. Shortening L_z is therefore the only lever with unbounded returns —
hence L_z = 80 rather than the originally proposed 110.

### 3.2 Two independent clean-window constraints

**Longitudinal** — the leading 3σ_d(t) tail must not reach the +z face:

    z₀ + v·t + 3·σ_d(t) < L_z/2   ⇒   t_IFW = 19.0 a.u.  (step 475)

**Transverse** — periodic images must not overlap. Criterion L_xy ≥ 6·σ_d(t), i.e.
≤0.54 % of the packet norm outside the transverse Wigner–Seitz cell
(1 − erf(3/√2)² = 0.0054):

| L_xy | t_transverse | binding |
|---|---|---|
| 35 | 16.0 | transverse |
| 40 | 18.4 | transverse |
| 45 | 20.8 | longitudinal (19.0) |
| **46** | **21.3** | **longitudinal (19.0)** |
| 50 | 23.2 | longitudinal (19.0) |

L_xy = 46 saturates the constraint with ~12 % transverse headroom over the
longitudinal bind. Widening beyond that buys nothing at σ_WP = 2, L_z = 80.

**Fit window** = [4.0, 19.0] a.u.: the lower edge drops the launch/injection
transient (the WP is orthogonalised against the occupied manifold at t=0, and the
bath needs ~1 plasma period, 2π/ω_p = 49 a.u.... see §3.4 caveat), the upper edge is
t_IFW. 40.7 Bohr of path — ample for a dT/ds slope.

### 3.3 Shell closure

Free-electron levels of the orthorhombic box enumerated and sorted; N = 218 sits in a
0.24 eV gap at r_s = 5.70 — the best (widest-gap, closest-to-legacy-density) closure
in the 44–50 Bohr transverse range. Verified numerically in-session, not assumed.

### 3.4 KNOWN CAVEAT — the plasma period exceeds the run

2π/ω_p = 2π/0.1272 Ha = **49.4 a.u.**, versus a 25.8 a.u. run and a 15 a.u. fit
window. The bath cannot complete a single plasma oscillation during the run. This is
an unavoidable consequence of the low jellium density (r_s = 5.70) chosen for
continuity with the existing body of work — the wake criterion
(5·2π/ω_p ≈ 247 a.u.) is geometrically impossible for a light projectile that
traverses the box in 26 a.u.

Per `.claude/rules/light-projectile-stopping.md` this is *expected and acceptable*:
we extract S as the **initial drag** over a near-constant-velocity window, not as a
steady-state wake result. It must be stated as a limitation in the notebooks: the
extracted S is a transient/initial-drag stopping power, not a converged
steady-state S(v). A denser bath (r_s ≈ 4, N ≈ 500) would shorten the plasma period
to 31 a.u. at ~2.5× the cost — recorded as a follow-up option, not in scope.

### 3.5 Rule amendment flagged

`ResearchProject/systems/jellium/shared/configs/boundary_rule.hpp::ifw_end_z` is a
**static-σ** rule (+L/2 − 3σ_launch) that ignores wavepacket spreading. For this run
it predicts t_IFW = 24.3 a.u. against the dispersion-aware 19.0 a.u. — a 22 %
over-estimate. Short legacy flights masked this. **Action:** add a
`ifw_end_t_dispersive(σ, L, v)` helper (solve z₀ + vt + 3σ_d(t) = L/2 numerically or
by a closed-form bound) and note the limitation in the header. Do NOT silently change
the existing helpers — legacy runs' recorded N_STEPS depend on them.

---

## 4. Stopping-power definitions

Projectile kinetic energy, two definitions (both already emitted per step by
`inqkit::observables::WPMomentumStats`):

| # | T | Column(s) | Meaning |
|---|---|---|---|
| **1** | ⟨p²⟩/2m | `e_kin_ha` | Full orbital kinetic energy (INQ's native measure) |
| **2** | ⟨p⟩²/2m | ½(`px_mean`²+`py_mean`²+`pz_mean`²) | Drift-only; discards localisation + scattering |

T₁ − T₂ = (3/2)·σ_p² summed over axes = **5.10 eV at t = 0** (= 3/(4σ_WP²)); its
*change* over the trajectory is the momentum-broadening contribution to apparent
stopping.

Projectile position, two definitions:

| # | s | Source |
|---|---|---|
| **3** | density centroid ⟨z⟩ of the WP KS orbital | `WPRealSpaceStats` (`z_mean`) + new circular-mean column |
| **4** | ∫₀ᵗ ⟨p_z⟩(t′)dt′ + z₀ | cumulative trapezoid of `pz_mean` |

**Exact identity:** the WP run has no ions, so the KS Hamiltonian is purely local
(kinetic + Hartree + ALDA) and Ehrenfest gives d⟨z⟩/dt = ⟨p_z⟩/m *exactly*. With no
CAP in this run, defs 3 and 4 MUST agree to numerical precision. Their comparison is
therefore a **validation check**, not an independent physics channel; deviation
localises to periodic wrap or WP-orbital norm leakage. (Contrast qsp5, where CAP
non-unitarity broke the identity at t≈5 — see `docs/handovers/qsp5-momentum-stopping.md`.)

**Extracted quantities:**

- WP: S_ij = −d(T_i)/d(s_j) for i∈{1,2}, j∈{3,4} — four numbers, of which S_13≈S_14
  and S_23≈S_24 by the identity above; the physics contrast is S_1* vs S_2*.
- Classical: S_cl = −d(½ m v_z²)/dz from the Ehrenfest ion trajectory, same fit window.

Fit: OLS slope of T vs s over t ∈ [4.0, 19.0] a.u., with (a) the residual plot, (b) a
bootstrap CI, and (c) a window-sensitivity scan (vary both edges ±3 a.u.) reported as
the systematic. Per `.claude/rules/number-rounding.md`, headline S to 2 s.f.

### 4.1 Relation to the `stopping-power-extraction` skill (added 2026-07-31)

The skill locks the headline method by geometry: for a **continuous traversal** it is
**Method A** — free-intercept slope of the deposited energy ΔE_total(x) after a fixed
20%-of-time transient cut — and −dT/dx is demoted to a conservation cross-check.

**That hierarchy inverts for the WP half, and the reason is structural, not a
preference.** The wavepacket *is* an occupied KS orbital, so its energy lives inside
`energy_total`. The system is closed, nothing absorbs, and E_total is constant to
**2.6e-4 eV** over the whole run. Method A's fit target is identically zero: the method
is **undefined here**, not merely imprecise. The four −dT_i/ds_j slopes are therefore
the measurement for this half — which is exactly why KS-orbital-dependent definitions
were requested.

The external reference comes from the **classical twin**, where the projectile sits
outside `energy_total` and Method A applies normally. Computed with the skill's own
kernels (`.claude/skills/stopping-power-extraction/stopping_power.py`):

| Channel | S (eV/Bohr) | r² |
|---|---|---|
| Method A, skill default (fixed 20% time cut) | 0.41 | 0.998 |
| Method A restricted to the WP window [4, 18.97] | **0.38** | 0.994 |
| Sanity: −dKE_ion/dx, same window | 0.38 | 0.994 |

Deposit and kinetic channels agree to **0.18%** — independent CSVs, so this is a direct
energy-conservation confirmation and the classical reference is trustworthy. Use the
**window-restricted 0.38** for like-for-like comparison against the WP numbers.

---

## 5. Work items

### 5.1 New code (inqkit) — needs tests before use

**W1. Circular-mean centroid in `WPRealSpaceStats`** (user-approved 2026-07-30).
Add columns `x_mean_circ`, `y_mean_circ`, `z_mean_circ` computed as

    ⟨z⟩_circ = (L_z/2π) · arg ⟨ψ| e^{i 2π z / L_z} |ψ⟩

(Resta/Mermin phase estimator — exact in a periodic cell; reduces to the naive
centroid for a well-localised packet). Post-processing unwraps the branch.
Existing naive columns are KEPT so the two can be cross-checked.

*Known-case test* (`code-test` skill, `inq-stack/tests/include/inqkit/observables/`):
analytic Gaussian placed at several offsets including straddling the boundary;
circular estimator must recover the offset to <1e-3 Bohr while the naive one fails.

**W2. Dispersion-aware IFW helper** in `boundary_rule.hpp` (§3.5), with
`static_assert` smoke values.

### 5.2 Run machinery

**W3. Gaussian UPF at σ_pot = 1.4142 Bohr** via `inqview.io.gaussian_psp` →
`shared/pseudopotentials/electron_gaussian_wpsigma2p0.upf`. Validate against the
existing `test_gaussian_psp.py` invariants (charge normalisation, r→0 limit,
asymptotic −1/r).

**W4. Config header** `bulk_ks_stopping_L46x46x80.hpp` — `Common_` base + `_WP` and
`_Classical` derived structs, following `electron_proj_E1000_L40x40x150.hpp`
(the orthorhombic precedent).

**W5. Ground state** `run_gs.cpp` → `checkpoints/gs_L46x46x80_orth_N218_dx0p40/`.
Shared by both runs.

**W6. `run.cpp` × 2** (WP + classical), cloned from
`run_wp_n162_L50_E100_sigma1_v2/run.cpp` with:
- orthorhombic cell (`systems::cell::orthorhombic`), not `cubic`
- **full energy decomposition ON**: `energy_{total,kinetic,hartree,xc,external,
  nonlocal,ion,ion_kinetic,exact_exchange,nvxc,eigenvalues}` + `proj_bg`
- `WPMomentumStats` / `WPRealSpaceStats` at `write_every = 1` (these ARE the
  stopping data; one extra FFT per step against ~109 in the propagator — negligible)
- WP-orbital density VTI at `WRITE_EVERY = 2` (323 frames — the user asked for the
  WP KS-orbital density at a cadence); total/system density at `WRITE_EVERY = 2`;
  complex wavefunction at `WF_WRITE_EVERY = 6` (108 frames)
- classical run additionally writes `projectile.csv` (z, v_z, KE per step)
- **final-timestep checkpoint + `LJ_RESUME` branch** per
  `.claude/rules/final-timestep-checkpoint.md`, plus interior checkpoints every 200
  steps per `.claude/rules/checkpoint-dont-block.md`

Estimated storage: ~14 GB per run, ~28 GB for the pair (966 GB free — fine).

**W7. SLURM dispatcher** `shared/bin/run-bulk-ks-stopping.slurm` (stock `inq`, not
`inq-study` — no CAP in this run), 1 GPU per run, both runs submittable in parallel.

### 5.3 Analysis

**W8. `analyse.py`** per run (energy ledger, kinematics, S extraction).

**W9. Run notebooks × 2** (`run-notebook` + `notebook-making` skills), each containing:
1. Density-matrix GIFs at the TOP, displayed inline
   (`.claude/rules/notebook-density-gif.md`): xz mid-y slices of **WP orbital
   density**, **total electron density**, and **induced Δn = n(t) − n(0)**, linear +
   log panels, fixed limits. For the classical run the "WP" panel is replaced by the
   projectile-localised induced density.
2. Step-by-step derivation of how S was computed (formulas, every term defined).
3. Each KE term (T₁, T₂, T₁−T₂) vs time, individually.
4. Each position term (s₃ naive, s₃ circular, s₄) vs time, individually, plus the
   Ehrenfest-identity residual s₃ − s₄.
5. T vs s for all four combinations, fit window shaded, slope + CI annotated.
6. Full energy decomposition vs time and the conservation check.
7. Takeaway + stated limitations (§3.4).

**W10. Cross-run comparison** in `hypotheses/bulk_ks_stopping/` — classical vs WP S,
and the WP−classical induced-density difference GIF.

---

## 6. Validation plan (`simulation-validation` / `validation-gates`)

**Tier A — cheap, must pass before production (no user approval needed):**
- A1. `WPRealSpaceStats` circular-centroid known-case test (W1) — analytic Gaussian.
- A2. `gaussian_psp` UPF invariants for σ_pot = 1.4142 (W3).
- A3. `boundary_rule` static_asserts for the new dispersive helper (W2).
- A4. GS sanity: N_electrons = 218 recovered by ∫n dV; occupations = 2.0 for the
  lowest 109 states and ≈0 above (closed shell); E_GS/N compared with the uniform-gas
  LDA value at r_s = 5.70. **OUTCOME (job 32400615):** E_kinetic = 7.396581 Ha vs
  the exact plane-wave sum 7.395907 (6.7e-4 Ha); E_total = −15.848701 Ha vs the
  analytic −15.830 (0.12 %); gap 0.2518 vs 0.2435 eV; occupations 1.9999988 / 2.3e-6.
  > **UNITS TRAP (cost a spurious gate failure).** ε_x = −0.9163/r_s is in
  > **RYDBERG**; in Hartree it is −0.4582/r_s. Same for the PZ81 correlation
  > constants (γ = −0.1423 Ry). Using the Rydberg form as Hartree doubles E_xc and
  > makes the predicted E_total ≈ 2× too negative (this plan's first draft said
  > −33.37 Ha). Gate hard only on unit-safe structural quantities — the discrete
  > kinetic sum, occupations, electron count, the gap — and treat
  > functional-dependent totals as reported, not gated.
- A5. WP injection report: norm_after ≈ 1, max_overlap with occupied manifold < 1e-3.
- A6. t=0 momentum gate: `pz_mean` = 2.7111, `sigma_pz2` = 1/(2σ²) = 0.125,
  `e_kin_ha` = ½(k₀² + 3σ_p²) = 3.8624 Ha (= 105.10 eV), with
  σ_p² = 1/(2σ²) = 0.125. Exact analytic values.

  > **CORRECTED 2026-07-30 (post-run).** This plan originally used
  > σ_p = 1/(2σ) ⇒ σ_p² = 1/(4σ²) = 0.0625, giving T₁ = 3.7686 Ha (102.55 eV)
  > and T₁−T₂ = 3/(8σ²) = 2.55 eV. **Those values are wrong** — they were copied
  > from an incorrect docstring in `wp_momentum_stats.hpp` (since fixed). For
  > ψ ∝ exp(−r²/2σ²) the correct momentum variance is σ_p² = 1/(2σ²) = 0.125.
  > The check that settles it: the real-space density std is σ/√2 = 1.41421, so
  > the wrong value implies σ_d·σ_p = 0.354 < ½ — a violation of the Heisenberg
  > bound. A Gaussian is minimum-uncertainty and must give exactly ½.
  > Confirmed against the live run to 5 decimal places.

**Tier B — short GPU jobs, ~20 min each, recommended:**
- B1. **Smoke run**: 20 steps of each run.cpp. Confirms it builds against stock `inq`
  on CSD3, fits in GPU memory, writes every expected artefact, and gives a
  measured cost/step → a real wall-time projection before committing.
- B2. **Vacuum dispersion control**: the same WP (σ=2, 100 eV) in an EMPTY 46×46×80
  box for 646 steps. Verifies σ_d(t) follows the free law on THIS grid, and gives the
  no-bath baseline that S must be measured against. Also directly tests the §3.2
  clean-window prediction.
- B3. **Energy conservation**: with no CAP the total energy must be conserved;
  |ΔE_total| over the run is the numerical-quality figure of merit (target ≲0.01 eV,
  cf. 0.0017 eV over 100 a.u. in the nocap qsp5 pair).

**User decision 2026-07-30: B1, B2 and B3 all APPROVED. C1 declined for now.**
Cadence locked at `WRITE_EVERY = 2` (323 density frames), `WF_WRITE_EVERY = 6`
(108 complex-wavefunction frames), per-step momentum/real-space stats.

**Tier C — expensive, user decides:**
- C1. dt convergence: repeat the WP run at dt = 0.02 (1292 steps) and compare S.
  *(Declined 2026-07-30; revisit if S becomes a headline report number.)*
- C2. dx convergence: GS + short RT at dx = 0.35 (user declined the GS-only check;
  offering it again here as an explicit Tier C item).
- C3. Denser bath r_s ≈ 4 (N ≈ 500) to bring the plasma period inside the run (§3.4).

**Nothing is declared correct on a green compile.** Every claim in the notebooks
cites which of the above produced it.

---

## 7. Open items / risks

| Risk | Mitigation |
|---|---|
| A100 memory: 4.6 GB orbitals × propagator workspace | B1 smoke run measures it; fall back to 2-GPU state-decomposition MPI if needed |
| Classical projectile may stop before t = 19 a.u. (light-projectile rule) | Extract S as initial drag over v ≥ 0.85·v₀; if it decelerates hard, that IS the result — report the deceleration sweep |
| WP orbital norm leaking into the bath (orthogonality violation over time) | `norm_check` column monitored per step; Ehrenfest residual s₃ − s₄ is the second detector |
| Plasma period > run length (§3.4) | Stated as a limitation, not a defect; C3 is the fix if wanted |
| `inq-study` not needed — confirm stock `inq` suffices | No CAP is constructed anywhere in these run.cpp files (contrast `wp_traversal_energy/run.cpp`, cause 8 in the CSD3 handover) |

---

## 8. Sources

- Wavepacket dispersion law: standard free-particle Gaussian result; verified
  empirically on this engine (`shared/bin/run-dispersion.slurm`).
- σ_WP / σ_pot √2 convention: `.claude/rules/sigma-wp-convention.md`.
- Light-projectile S extraction: `.claude/rules/light-projectile-stopping.md`.
- Prior momentum-based stopping methodology and its pitfalls:
  `docs/handovers/qsp5-momentum-stopping.md`.
- Free-electron-gas magic numbers: `docs/sources/free-electron-gas-magic-numbers.md`.
- Resta phase-estimator for a periodic position operator: R. Resta,
  Phys. Rev. Lett. 80, 1800 (1998) — **to be verified and written up as a source note
  in `docs/sources/` before W1 is marked complete** (currently cited from memory;
  treat as unverified until then).

---

## 9. Density-replica pair at r_s ≈ 4 (added 2026-07-31)

**Question.** Does the 6.6× S_classical/S_WP gap measured at r_s = 5.702 close
when the bath is denser? This is the discriminating test for the
initialisation density-clearing hypothesis (handover, 2026-07-31): clearing is a
screening-length effect and must be density-dependent; free-packet dispersion is
not.

**Design principle — one variable.** σ_WP = 2, E = 100 eV, launch z = −32,
dt = 0.04, N_STEPS = 646 and L_z = 80 are held **bit-identical** to the
r_s = 5.702 pair. Only the bath density moves.

| | r_s = 5.702 | r_s = 3.987 |
|---|---|---|
| Cell / N / n | 46×46×80, 218, 1.2878e-3 | 40×40×80, **482**, **3.7656e-3** |
| dx / grid / states | 0.40, 2.65 M, 130 | **0.50**, 1.024 M, 262 |
| ħω_p / plasma periods | 3.46 eV / 0.52 | 5.92 eV / **0.89** |
| fit window | [4, 18.97] | **[4, 18.43]** |

**Why not r_s = 3 (the original target).** In the 46×46×80 box it needs N = 1514
→ 778 states → ~99 GB of orbital storage for ETRS's three copies, against an
80 GB A100. Infeasible, not merely expensive. 42×42×80 lands at ~74 GB (8 %
headroom) and was also rejected.

**dx = 0.50 is verified, not assumed.** WP momentum moments on the exact
80×80×160 grid reproduce ⟨p_z⟩, T₁ and T₁−T₂ to machine precision; k_max = 3.461
against Nyquist 6.283. **Declare** in any cross-pair comparison that the two
pairs use different dx — each pair's internal ratio is unaffected.

**Window inversion.** At L_xy = 40 the transverse limit (18.43) binds before the
longitudinal one (18.97) — the reverse of the L=46 pair. FIT_T1 takes the min, so
the analysis is correct, but a strict pair-to-pair comparison should re-fit the
r_s = 5.702 run on [4, 18.43]. Pinned by `static_assert`.

**Validation.** Same Tier A4 analytic GS gates, recomputed for this system:
E_kinetic = 33.593666 Ha (exact plane-wave sum), E_hartree = 0, E_total in
[−38.30, −36.06], n_occupied = 241, gap 0.252 eV. GS exits 3 on failure and the
chain is `afterok`, so a bad ground state cannot seed production.

**Autonomous chain.** 32439807 (GS) → 32439808 (wp) + 32439810 (classical) →
32439811 (both notebooks). Projected ~2.4 h / ~3.2 h from the 1.57× cost factor;
12 h walls give ~4× headroom, which is what keeps the non-resumable classical
half safe.

---

## 10. Interaction-energy retrofit of the remaining twin pairs (added 2026-08-01)

**Goal.** Give every non-sigma=1 bulk twin pair the P/S/B pairwise decomposition
(`.claude/rules/decomposed-interaction-energies.md`) so classical and WP halves
are comparable in representation-independent terms. sigma=1 excluded by user
instruction.

**Scope.** 3 pairs = 6 runs: `bulk_ks_stopping` (sigma 2, r_s 5.702),
`bulk_ks_stopping_rs4` (sigma 2, r_s 3.987), `bulk_ks_stopping_rs4_sigma3`
(sigma 3, r_s 3.987). `bulk_ks_stopping_sigma3` already carries it.

**Both halves, not just the classical twin.** Absolute E_PP carries the
charged-cell G=0 gauge, so only `dE_PP = E_PP(WP) - E_PP(classical)` measured in
the SAME cell is meaningful. A pair with only one instrumented half yields
nothing.

**Chain.** 32526619 (preflight build, all 6) -> 32526620 (array 0-5, afterok).
Each task runs `scripts/verify_interactions_closure.py` and exits non-zero if
the pairwise terms do not sum back to INQ's scalars.

### NEW HARD CONSTRAINT — projectile cloud clipping bounds the fit window

Discovered 2026-08-01 while validating the closure verifier on the completed
sigma=3 r_s=5.702 classical run. E_PP of a rigid Gaussian cloud must be a
constant of the motion, and it is — **bit-exactly** (spread 0.000e+00 Ha) while
the cloud lies fully on the grid. It then DECAYS over the final rows as the
cloud is clipped by the +z box face and `norm_proj` falls 1.000000 -> 0.994247.

- Onset for that run: **t = 21.04 a.u.** (z = +27.15, i.e. ~6 sigma_pot from the
  face). Its fit window ends at 19.48 (common cross-sigma 9.37), so **no existing
  result is contaminated** — but the margin is only ~1.6 a.u.
- This is NOT egg-box/grid-phase error: correlation of E_PP with sub-grid phase
  is 0.03. It is charge leaving the grid.
- `verify_interactions_closure.py` now reports the onset per run as the hard
  upper bound on any fit window. **Any future S(v) fit or dE_PP integration must
  end before it.**

Consequence for analysis: when fitting, take `t_max = min(FIT_T1_AU,
clipping_onset)` rather than trusting the config's window alone, since the
config's window was derived from interference/transverse limits and does not
know about projectile-cloud clipping.
