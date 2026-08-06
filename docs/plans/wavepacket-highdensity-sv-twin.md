# Plan: wavepacket twin of the high-density classical S(v) benchmark

**Created:** 2026-07-30
**Branch:** `quantum-stopping-power`
**Machine:** CSD3, `ampere` partition (A100), account `mphil-nikiforakis-skcb2-sl2-gpu`
**Status:** design drafted; three decisions pending with the user (§5)

Twin of campaign `classical-highdensity-sv`
(`docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md`,
handover `docs/handovers/classical-highdensity-sv-benchmark.md`).

---

## 1. Objective

Re-run the six-point high-density S(v) benchmark of the localised jellium slab with
a **quantum electron wavepacket** replacing the classical Gaussian-charge projectile,
holding **every other physical parameter fixed**, so the classical curve
(S = 1.09 → 0.28 eV/Bohr over v = 2.0 → 4.5) becomes the like-for-like reference for
a quantum stopping-power curve.

Stopping is extracted with the KS-orbital definitions locked in
`docs/plans/bulk-jellium-ks-stopping.md` §4 (T₁/T₂ × s₃/s₄), plus the pairwise
Coulomb ledger as the bath-energy channel (§4 below).

---

## 2. The width mapping (the one parameter that is NOT copied verbatim)

Per `.claude/rules/sigma-wp-convention.md`, σ always means **σ_WP**, the wavepacket
ψ-width, and the classical Gaussian *charge* std is the derived internal quantity
σ_pot = σ_WP/√2.

The classical campaign ran `PROJ_SIGMA_POT = 0.5/√2 = 0.35355` Bohr
(`shared/configs/slab_n100_L35x35x85.hpp:64`), which is by construction the σ_pot of
**σ_WP = 0.5**. The matched wavepacket therefore uses

    sigma_WP = 0.5 Bohr        (psi ~ exp(-r^2 / 2 sigma_WP^2))
    sigma_d(0) = sigma_WP/sqrt(2) = 0.35355 Bohr   == the classical charge std

so the WP's **density** std at t = 0 equals the classical Gaussian charge std exactly.
Both halves are labelled **σ = 0.5** in every table, axis and caption.

`inqkit::WavePacket{}.sigma(...)` takes σ_WP directly, and the classical `run.cpp`
takes `LJ_SIGMA` = σ_WP and forms `SIGMA_POT = LJ_SIGMA/√2` internally — so **both
halves are driven with the same env value `LJ_SIGMA = 0.5`**. No conversion is applied
at the call site; the √2 lives inside each binary.

---

## 3. Parameter table — classical (as executed) vs wavepacket (proposed)

### 3.1 Shared, held fixed

| Quantity | Value | Source |
|---|---|---|
| Cell | 35 × 35 × 85 Bohr, orthorhombic | `slab_n100_L35x35x85.hpp:40-42` |
| Boundary | `periodicity(2)` — x,y periodic; z open **for electrostatics only** | campaign `<rules>`; §5.2 |
| Grid spacing | dx = 0.50 Bohr → 70 × 70 × 170 = 833 k points | `:43` |
| Cutoff | 19.74 Ha = 537 eV | `:12` |
| Slab | 25 Bohr thick (half-width 12.5), edge_width 1.0, centred z = 0 | `:47-49` |
| Bath | N = 100 electrons, n₀ = 3.2653e-3 → **r_s = 4.183** | `:52-54` |
| Extra states | 24 (⇒ 74 states total) | `:55` |
| Temperature | 0.00862 eV (~100 K) | `:56` |
| XC | LDA ground state, ALDA in real time | campaign |
| CAP | **none** | campaign `<rules>` |
| dt | 0.04 a.u. | orchestrate.py |
| Launch z₀ | −24.0 Bohr (11.5 Bohr standoff from the −12.5 slab face) | orchestrate.py `LAUNCH_Z` |
| σ (label) | **0.5** (σ_WP; classical σ_pot = 0.35355) | §2 |
| Charge | −1 e | both |
| Ground state | `slab_n100_L35x35x85_dx0p5_per2`, E_GS = 207.18322156141 Ha | handover |

### 3.2 What necessarily differs

| | Classical | Wavepacket |
|---|---|---|
| Representation | rigid Gaussian **charge**, moving external perturbation | extra **KS orbital**, part of the system |
| Mass | m = 1 (explicit `PROJ_MASS`) | m = 1 implicitly (it *is* an electron) |
| Drive | Ehrenfest velocity-Verlet, Hellmann–Feynman force | none — the TDDFT Hamiltonian propagates it |
| Velocity control | v₀ = K0/m | group velocity ⟨p_z⟩/m = k₀ ⇒ **k₀ = v₀** |
| Energy bookkeeping | E_total **changes** (external perturbation does work) | E_total **exactly conserved** (time-independent H) |
| Can leave the box | yes (z-open Poisson clips the charge) | **no** — orbital wraps on the FFT grid (§5.2) |

### 3.3 Per-run grid (the six points)

`N_steps(classical) = ceil(1.4 × 69 / (0.5 v dt))` (orchestrate.py `n_steps_for`);
WP runs are **1.5×** that, at the same dt, per user instruction 2026-07-30.

| v₀ = k₀ | N_cl | t_cl (a.u.) | **N_wp** | **t_wp (a.u.)** | KE_cl (eV) | T₂ = ⟨p⟩²/2 (eV) | T₁ = ⟨p²⟩/2 (eV) | ckpt every | frame every |
|---|---|---|---|---|---|---|---|---|---|
| 2.0 | 2415 | 96.60 | **3623** | 144.92 | 54.4 | 54.4 | 136.1 | 725 | 12 |
| 2.5 | 1932 | 77.28 | **2898** | 115.92 | 85.0 | 85.0 | 166.7 | 580 | 10 |
| 3.0 | 1610 | 64.40 | **2415** | 96.60 | 122.5 | 122.5 | 204.1 | 483 | 8 |
| 3.5 | 1380 | 55.20 | **2070** | 82.80 | 166.7 | 166.7 | 248.3 | 414 | 7 |
| 4.0 | 1208 | 48.32 | **1812** | 72.48 | 217.7 | 217.7 | 299.3 | 362 | 6 |
| 4.5 | 1074 | 42.96 | **1611** | 64.44 | 275.5 | 275.5 | 357.1 | 322 | 5 |

Total WP work: 14 429 steps.

**Classical reference results** (`sv_sweep/S_summary.csv`, for overlay):

| v₀ | v_final | v_mean | S (eV/Bohr) | E_absorbed (eV) |
|---|---|---|---|---|
| 2.0 | 1.42 | 1.82 | 1.087 | 27.2 |
| 2.5 | 2.11 | 2.41 | 0.970 | 24.3 |
| 3.0 | 2.77 | 2.95 | 0.709 | 17.7 |
| 3.5 | 3.36 | 3.47 | 0.509 | 12.7 |
| 4.0 | 3.91 | 3.98 | 0.374 | 9.3 |
| 4.5 | 4.44 | 4.49 | 0.283 | 7.1 |

S ∝ v^−1.72 (Bethe tail).

### 3.4 The localisation-energy asymmetry (physical, unavoidable)

σ_p² = 1/(2σ_WP²) = 2.0 Bohr⁻², so

    T1 - T2 = 3/(4 sigma_WP^2) = 3.0 Ha = 81.6 eV

A σ_WP = 0.5 packet carries **81.6 eV of localisation (zero-point) energy** on top of
its drift energy. At v = 2.0 the WP's total kinetic energy is 136 eV against the
classical projectile's 54 eV — 2.5×. This is not a modelling error; it is what
"same width, quantum" means, and it is precisely the T₁ vs T₂ contrast the study is
designed to measure. It must be stated in every notebook.

---

## 4. Stopping-power extraction

Four combinations, per `docs/plans/bulk-jellium-ks-stopping.md` §4:

| T | column | s | source |
|---|---|---|---|
| T₁ = ⟨p²⟩/2m | `e_kin_ha` (`WPMomentumStats`) | s₃ = density centroid | `z_mean` / `z_mean_circ` (`WPRealSpaceStats`) |
| T₂ = ⟨p⟩²/2m | ½(px_mean²+py_mean²+pz_mean²) | s₄ = ∫⟨p_z⟩dt | cumulative trapezoid of `pz_mean` |

S_ij = −dT_i/ds_j. **Use `z_mean_circ`** — the naive centroid is discontinuous across
the periodic z face (§5.2) and the circular estimator is exact in a periodic cell.

**Bath-energy channel (the WP analogue of the classical E_absorbed).** A WP run has a
time-independent Hamiltonian, so `energy_total` is *exactly conserved* and the
classical Definition-2 (`E_total(plateau) − E_GS`) has **no WP analogue**. The energy
deposited in the bath must instead be isolated from the pairwise Coulomb ledger
(`inqkit::jellium::interaction_energies`), whose documented WP closure is

    E_hartree = E_SS + E_PS + E_PP
    E_external = E_SB + E_PB

with n_P = |ψ_WP|² (the WP orbital density) and n_slab = n_total − |ψ_WP|². The bath's
internal energy change ΔE_SS + ΔE_SB (+ Δ of the bath's kinetic and xc share) is the
quantity to compare against the classical 27.2 → 7.1 eV. Constancy of `energy_total`
is a hard correctness gate.

---

## 5. Open decisions (user)

### 5.1 σ_WP = 0.5 gives a very short clean window

Free-Gaussian dispersion (verified on this engine, `run-dispersion.slurm`):

    sigma_d(t) = sqrt(sigma_WP^2/2 + t^2/(2 sigma_WP^2))

At σ_WP = 0.5 the spreading rate is **1.414 Bohr per a.u.** — 4× the σ_WP = 2.0 case
of the bulk study. Consequences:

| t (a.u.) | 0 | 1 | 2 | 4 | 10 | 20 | 96.6 | 145 |
|---|---|---|---|---|---|---|---|---|
| σ_d (Bohr) | 0.35 | 1.46 | 2.85 | 5.67 | 14.2 | 28.3 | 137 | 205 |

- **Transverse images overlap at t = 4.12 a.u. (step 103)**, where 6σ_d = 35 Bohr =
  L_xy (the ≤0.54 %-outside-cell criterion of the bulk plan §3.2). This bound is
  **velocity-independent** and binds hardest.
- The leading 3σ tail reaches the +z face at t = 7.6–10.7 a.u. (step 190–266).
- By the end of even the *classical* duration the packet is smeared over ~1.6 box
  lengths; it is not a localised projectile.

So the WP is "projectile-like" for roughly the first **100–260 steps of 1611–3623**.
Lengthening the runs 1.5× does not extend the physics window — it extends the record
of a delocalised packet. Options in §6.

### 5.2 The wavepacket cannot leave the box (verified)

`periodicity(2)` changes **only** the Poisson kernel (`solvers/poisson.hpp:189,206`
→ Rozzi slab kernel), ionic replicas (`ionic/periodic_replicas.hpp:39`), the Ewald
branch (`ionic/interaction.hpp:282-312`), the kick gauge, and spatial partitions.
It is **never** consulted by the wavefunction basis or the kinetic operator:
`basis/fourier_space.hpp` builds G-vectors unconditionally over all three axes and
`hamiltonian/ks_hamiltonian.hpp:200-204` applies −½∇² as a plain 3-D FFT multiply.

⇒ **A KS orbital travelling in +z wraps around and re-enters at −z.** It does not
leave the box and does not reflect. Confirmed empirically in-repo:
`docs/handovers/pbc-open-z-oscillation.md:20` ("wavefunction always wraps on the FFT
grid (p2 switches electrostatics only)") and
`docs/handovers/wp-localised-jellium-solving-cap.md:45` ("longer no-CAP run (t=16,
~2 wraps) conserves energy to 0.15 meV → periodic re-entry is a clean unitary op").

This breaks the classical campaign's central mechanism — projectile exits ⇒
E_electronic plateaus ⇒ E_absorbed — for the quantum half. The wrap is unitary and
energy-conserving, but the re-entering packet arrives *behind* the slab and
contaminates the bath. Options in §6.

### 5.3 The classical per-step data did not survive the machine migration

Only `sv_sweep/S_summary.csv` (6 rows) and per-run `REPORT.md` / `result.json` are in
git. The raw `dyn/results/vXpY/` outputs (observables.csv, projectile.csv,
interactions.csv, density VTIs) lived on `/local/data/public/skcb2/tddft`, which does
not exist on CSD3. The slab ground state is likewise absent and **must be recomputed
here** (gate: reproduce E_GS = 207.183 Ha).

Without re-running the classical half there can be no per-step overlays, no shared
energy-ledger comparison, and no WP−classical induced-density difference GIF — only
a six-number S(v) comparison.

---

## 6. Resolutions — LOCKED (user, 2026-07-30)

1. **Width — σ_WP = 0.5, unchanged.** The exact match to the classical
   σ_pot = 0.35355. The 4.12 a.u. transverse window and the 81.6 eV localisation
   energy are reported as stated limitations, not designed away.
2. **Boundary — TWO CAPs.** 12.5 Bohr per z face, |η| = 1 Ha: sin² absorbing bands
   over z ∈ ±[30, 42.5]. This DEPARTS from the campaign's CAP-free rule,
   necessarily — that rule relies on the projectile being an external charge that
   can leave the box, which an orbital cannot do (§5.2). Applied as η = −1.0 Ha,
   the sign that absorbs (+1.0 would be a gain medium). Consequence:
   `energy_total` is no longer conserved, so the correctness gate becomes
   norm/absorption monitoring plus the ledger closure checks.
3. **Classical — NOT re-run.** Comparison is against the published
   `S_summary.csv` values only (S = 1.087 → 0.283 eV/Bohr). No per-step overlays
   and no WP−classical difference GIFs are therefore possible; the synthesis
   notebook compares S(v) curves.
4. **Grid — dx = 0.40, velocity grid CUT to four points.** v = 2.0, 2.5, 3.0, 3.5.
   Their σ_pz² errors are +0.05 / +0.26 / +1.24 / +5.06 % and ⟨p_z⟩ errors ≤0.47 %.
   v = 4.0 (+17.9 %) and v = 4.5 (+55.1 %) are EXCLUDED rather than caveated
   (§6b item 6). The quantum curve is thus a 4-point curve against the classical
   6-point curve. Recoverable later at dx = 0.30 (every moment ≤0.11 %).
5. **Vacuum controls — ADDED** (W7): one CAP-only free-WP run per velocity, same
   grid and step count, to separate real stopping from CAP attrition (§6b item 4).

**Fit window: t ∈ [~0.5, 4] a.u. (steps ~12–100)** — set by two independent limits
that agree: transverse image overlap at 4.12 a.u., and the onset of CAP attrition.

### Final per-run table (as submitted)

| idx | v = k₀ | N_steps | t (a.u.) | density/ | wavefn/ | ckpt/ | σ_pz² err |
|---|---|---|---|---|---|---|---|
| 0 | 2.0 | 3623 | 144.9 | 12 | 36 | 724 | +0.05 % |
| 1 | 2.5 | 2898 | 115.9 | 10 | 30 | 579 | +0.26 % |
| 2 | 3.0 | 2415 | 96.6 | 8 | 24 | 483 | +1.24 % |
| 3 | 3.5 | 2070 | 82.8 | 7 | 21 | 414 | +5.06 % |

11 006 production steps at dx = 0.40 (1.64 M grid points, 74 states).

---

---

## 6b. Results of the CAP validation replica (2026-07-30, jobs 32416846 / 32417361)

`scripts/wp_highdensity_sv/cap_check/run.cpp` — free WP (σ_WP = 0.5, k₀ = 2,
launch z = −24) in an EMPTY 35×35×85 `periodicity(2)` box, production CAP
(two sin² bands, 12.5 Bohr/face, |η| = 1 Ha), `non_interacting`, 1200 steps.
Figure: `hypotheses/wp_highdensity_sv/cap_check/cap_validation.png`.

**1. CAP geometry CONFIRMED.** The binary reports
`mid=0.42647059 width=0.14705882 -> +z band [30, 42.5] Bohr, -z band [-42.5, -30]`.
The fractional-coordinate mapping is right.

**2. The orbital really does wrap without a CAP — measured.** In the CAP-off
control the circular centroid runs −24 → +41.4 (t = 32) → **−28.2 (t = 40)**:
it crossed the +z face and re-entered at −z, with norm conserved to 0.998 over
48 a.u. (a clean unitary wrap). This is the direct empirical confirmation of §5.2.

**3. The CAP absorbs with NO reflection.** With the CAP on, norm falls
1 → 0.226 over 48 a.u. and `min ⟨p_z⟩ = +0.61` — never negative, so nothing
bounces off the absorbing band. |η| = 1 Ha over 12.5 Bohr is adiabatic, as the
vacuum study found at η = −1.0/W = 15.

**4. ⚠️ THE CAP ITSELF PRODUCES A SPURIOUS DECELERATION.** This is the important
finding. In a run with *no bath, no forces, nothing but free dispersion*, the
surviving packet's mean momentum falls

| t (a.u.) | 0 | 4 | 8 | 12 | 16 | 24 | 48 |
|---|---|---|---|---|---|---|---|
| norm/norm(0) | 1.000 | 0.998 | 0.984 | 0.953 | 0.844 | 0.568 | 0.226 |
| ⟨p_z⟩ | 1.99 | 2.00 | 2.05 | 2.03 | 1.81 | 1.30 | **0.61** |

The CAP-off control holds ⟨p_z⟩ = 1.985 flat over the same interval, so the
entire drop is CAP attrition: because σ_WP = 0.5 spreads at 1.414 Bohr/a.u., the
packet's *leading* edge reaches the +z band first and is removed preferentially,
dragging ⟨p_z⟩ of what remains downward. (Note this is the OPPOSITE sign to the
qsp5 expectation that a CAP biases ⟨p_z⟩ upward by eating the slow tail.)

**Consequence:** any S = −dT/ds fitted over a window where the CAP is active
measures the CAP, not the jellium. Two independent constraints now agree on the
same window:

- transverse periodic images overlap at **t = 4.12 a.u.** (§5.1)
- CAP attrition is < 0.3 % of norm and ⟨p_z⟩ is unbiased for **t ≲ 4 a.u.**

⇒ **Fit S over t ∈ [~0.5, 4] a.u. (steps ~12–100).** Everything later is recorded
(density GIFs, absorption physics) but is not slope data.

**5. ACTION: a vacuum control per velocity.** The cheapest rigorous fix is to run
the same free-WP-plus-CAP replica at each of the six k₀ (1 state instead of 74,
so ~50× cheaper than a production run). That gives the CAP-only ⟨p_z⟩(t) and
T(t) baseline at each velocity, so the bath's true contribution is the
*difference* between the slab run and its vacuum twin. Added as W7.

**6. dx = 0.50 fails badly at high velocity — MEASURED.** The `cap_v4p5` variant
aborted on its own t=0 gates, which is the gate working as intended:

| quantity at t=0, v=4.5, dx=0.5 | analytic | measured | error |
|---|---|---|---|
| ⟨p_z⟩ | 4.5 | 3.44 | **−23.6 %** |
| σ_pz² | 2.0 | 9.05 | **+353 %** |
| T₁ − T₂ | 3.0 Ha | 6.53 Ha | +118 % |

Real-space quantities (norm, centroid, density std) are all fine — the corruption
is purely in momentum space, from the k-distribution folding at k_Nyq = π/dx.

A fold model (wrap N(k₀, σ_p²) into [−k_Nyq, k_Nyq)) reproduces this (predicts
⟨p_z⟩ = 3.20, σ_pz² = 10.27) and therefore predicts the other grids:

| dx | ⟨p_z⟩ err @v=4.0 | σ_pz² err @v=4.0 | ⟨p_z⟩ err @v=4.5 | σ_pz² err @v=4.5 | cost |
|---|---|---|---|---|---|
| 0.50 | −16.7 % | +205 % | −29.0 % | +414 % | 1.0× |
| 0.40 | −1.3 % | **+17.9 %** | −3.1 % | **+55.1 %** | 2.0× |
| 0.35 | −0.10 % | +1.4 % | −0.31 % | +5.7 % | 2.9× |
| 0.30 | −0.00 % | +0.02 % | −0.01 % | +0.11 % | 4.7× |

**This supersedes the dx = 0.40 choice of 2026-07-30**, which was made on a
*tail-fraction* estimate (0.89 % of weight beyond Nyquist at v = 4.5) that badly
understated the moment error — the aliased weight lands near k ≈ −6, so it has
enormous leverage on both ⟨p_z⟩ and ⟨p²⟩. σ_pz² sets T₁ − T₂, a headline
observable, so +55 % at v = 4.5 is not acceptable. Re-decision required.

---

## 7. Work items

- **W1.** Recompute the slab GS on CSD3 → `shared_gs/slab_n100_L35x35x85_dx0p5_per2/`;
  gate on E_GS = 207.183 Ha, ∫n dV = 100, r_s = 4.18.
- **W2.** `scripts/classical_highdensity_sv/wp/run.cpp` — WP fork: inject σ_WP = 0.5,
  k₀ = v at z = −24; static slab background perturbation; full energy decomposition;
  `WPMomentumStats` + `WPRealSpaceStats` every step; pairwise ledger every step with
  n_P = WP orbital density; total / WP-orbital / induced density VTIs; complex WP
  wavefunction at a coarser cadence; **≥4 retained numbered checkpoints** + final,
  `LJ_RESUME` branch.
- **W3.** t=0 analytic gates (abort before burning GPU): ⟨p_z⟩ = k₀,
  σ_pz² = 1/(2σ²) = 2.0, T₁ = ½(k₀²+3σ_p²), T₁−T₂ = 3.0 Ha, centroid = −24,
  σ_d = 0.354, norm ≈ 1, max_overlap < 1e-3.
- **W4.** SLURM dispatcher, six jobs submittable in parallel (1 A100 each, dependency
  on the GS job).
- **W5.** Per-run `analyse.py` + run notebook each (density-matrix GIF at top per
  `.claude/rules/notebook-density-gif.md`; step-by-step S derivation; T₁/T₂/T₁−T₂;
  s₃ naive/circular/s₄ + Ehrenfest residual; all four S_ij; full ledger +
  conservation; stated limitations).
- **W6.** Synthesis notebook in `hypotheses/classical_highdensity_sv/wp_sv_sweep/`:
  quantum S(v) vs the classical curve.

## 8. Cost

Reference: bulk WP run (2.65 M points, 129 states, 646 steps) = 5531 s, 40 GB.
This slab is 833 k points (3.2× smaller) and 74 states ⇒ roughly 2–3 s/step.
14 429 WP steps ⇒ ~8–12 GPU-hours total; longest single run (3623 steps) ~2–3 h,
inside the 24 h partition limit. Storage ~10–14 GB per run (~70 GB for six), against
895 GB free. Re-running the classical half adds ~6 GPU-hours and ~25 GB.

## 9. Sources

- Campaign + handover: `docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md`,
  `docs/handovers/classical-highdensity-sv-benchmark.md`.
- KE/centroid definitions: `docs/plans/bulk-jellium-ks-stopping.md` §4.
- σ convention: `.claude/rules/sigma-wp-convention.md`.
- Ledger closure: `inq-stack/include/inqkit/jellium/interaction_energies.hpp` header.
- periodicity(2) semantics: `inq/src/solvers/poisson.hpp:189,206`,
  `inq/src/basis/fourier_space.hpp:60-151`, `inq/src/hamiltonian/ks_hamiltonian.hpp:200-204`
  (verified 2026-07-30, independent read).
