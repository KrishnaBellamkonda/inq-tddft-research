# Stopping power in r_s = 5.69 jellium — regime classification and simulation plan

This document is a reference for the Claude Code agent driving the next round of simulations. Its purpose is threefold:

1. Codify the *regime classification* of electronic stopping in jellium in a way that the agent can use to make decisions about new simulations and to interpret existing ones.
2. Lay out the *operational definitions* of stopping power for classical-projectile and wave-packet simulations, so the agent can write or extend postprocessing pipelines consistently.
3. Propose two to three *concrete simulations* fitting a 9-hour budget on two GPUs (one for classical, one for wave packet), with clear scientific motivation.

The audience is a Claude Code agent that has direct access to the `ResearchProject/systems/jellium/` tree and the `inq-stack/python/inqview/postprocess/` pipeline. The agent should treat every claim here as a starting point and verify it against the actual repository before launching simulations.

---

## 1. Resolved questions and shared run configuration

The three blocking questions are resolved as follows. Each shapes the rest of the plan.

**Q1. inqkit per-element bulk-copy "bug" is actually a VRAM ceiling, NOT a software defect.** Earlier WP runs at `dx ≤ 0.30 / 102 states / 24 GB` failed because the smaller grid spacing inflates the orbital count and density buffers beyond what the GPU memory can hold. The conclusion: all WP runs in this campaign must stay at `dx = 0.40` Bohr, and **WP-1500 eV is not feasible** — it requires a finer dx to keep Nyquist headroom (WP at 1500 eV has `k_max = k₀ + 3 σ_k = 11.10` Bohr⁻¹, against `k_Nyquist = π/dx = 7.854` Bohr⁻¹ at dx=0.40 → 41 % over-Nyquist, heavily aliased). The matched WP-vs-classical comparison therefore cannot extend above ~600 eV; the existing classical-1500 eV run remains a useful *unpaired* high-v anchor, but no WP counterpart exists or can be made on this hardware.

**Q2. Wave-packet width σ = 5.0 Bohr fixed across the velocity sweep.** Inherited from `Base_N162_L50_E1p5` and used at both 100 eV and 600 eV WP runs. The corresponding `k₀σ` values are:

| E (eV) | k₀ (Bohr⁻¹) | σ (Bohr) | k₀σ | Packet vs de Broglie wavelength (σ / λ_dB = k₀σ / 2π) |
|---|---|---|---|---|
| 50 (proposed) | 1.92 | 5.0 | 9.6 | 1.53 |
| 100 (done) | 2.71 | 5.0 | 13.6 | 2.16 |
| 300 (proposed) | 4.69 | 5.0 | 23.5 | 3.74 |
| 600 (done) | 6.64 | 5.0 | 33.2 | 5.29 |

At every energy in the proposed sweep, k₀σ ≫ 1 and σ ≥ λ_dB. **The WP runs are therefore deep in the classical-packet limit, not the quantum-packet limit.** Each WP "projectile" is a localized object many de Broglie wavelengths wide. Two consequences:

- Projectile-wavefunction-spread effects (packet broadening, diffraction, single-collision quantum reflection) are negligible by design at this σ.
- The genuine differences between WP and classical-projectile runs must come from the host side: exchange and Pauli-blocking against the WP electron (which is identical to bath electrons, while the classical projectile uses a custom pseudopotential and is distinguishable), and the finite spatial extent of the source (a 5-Bohr Gaussian creates a different wake than a point charge).

**Q3. `observables.csv` is generated for every run** (classical and WP) at the same cadence with the same column schema (`step, time_au, energy_total, energy_kinetic, energy_hartree, energy_xc, current_x/y/z, dipole_x/y/z, density_l2`). No format conversion needed. The existing 1500 eV CSV is representative.

### 1.1 Shared run configuration (read from the user-supplied config files)

All five existing runs (and any new ones in this campaign) share the configuration below. New runs differ from this template only by `WP_EKIN_EV` and the consequent `WP_K0`, `PROJ_VEL_Z`, and `N_STEPS`.

| Parameter | Value | Notes |
|---|---|---|
| Cell | Cubic 50 × 50 × 50 Bohr, periodic | INQ centred Cartesian: `r ∈ [−25, +25]` Bohr |
| Jellium electrons | N = 162 | Closed-shell with `EXTRA_STATES = 20` (or 4 in `Base`) |
| Density | n = N / V = 1.296 × 10⁻³ Bohr⁻³ | r_s = 5.69 |
| Grid spacing | dx = 0.40 Bohr | k_Nyquist = π / 0.40 = 7.854 Bohr⁻¹ |
| Time step | dt = 0.020 a.u. | ETRS stable while `dt · v < dx`; at v=10.5, dt would have to drop to 0.005 |
| Total steps | N_STEPS = ceil((L_box/2 − z_launch) / v_initial / dt) = ceil(35 / v_initial / dt) | "Maximum traversal" rule: centroid reaches +z box face at t_end |
| WP envelope width | σ = 5.0 Bohr | Inherited from `Base_N162_L50_E1p5` |
| WP central wavenumber | k₀ = √(2 · E_eV / 27.2114) Bohr⁻¹ | `k0_from_ev` in `base.hpp` |
| Projectile launch | (0, 0, −10) Bohr | Same for WP and classical, identical across the E-sweep |
| Projectile direction | +z | All projectiles move toward +z box face |
| Classical projectile species | "H" with custom UPF + mass override | `PROJ_MASS_AMU = 1.0 / 1822.8885` → m = m_e |
| Classical projectile UPF | `shared/pseudopotentials/electron-ONCV-1.2.upf` | Pseudo-hydrogen tuned to match an electron |
| Total propagation time | t_total = 35 / v_initial a.u. | Single-pass; centroid reaches +25 Bohr exactly at t_end |

The propagation-time rule deserves emphasis: the simulation ends *exactly* when the projectile centroid reaches the +z box face at its initial velocity. This means the total wall-time per run scales as `1/v_initial`. Low-v sims are necessarily long; high-v sims are short.

### 1.2 Why WP-1500 cannot run at dx = 0.40

For completeness, the Nyquist diagnostic at each existing or proposed energy:

| E (eV) | k₀ (Bohr⁻¹) | k_max = k₀ + 3σ_k (Bohr⁻¹) | k_max / k_Nyquist | Status |
|---|---|---|---|---|
| 50 (proposed) | 1.92 | 2.52 | 0.32 | very clean |
| 100 (done) | 2.71 | 3.31 | 0.42 | very clean |
| 300 (proposed) | 4.69 | 5.29 | 0.67 | clean |
| 600 (done) | 6.64 | 7.24 | 0.92 | 8 % headroom, marginal but ok |
| 1500 (classical only) | 10.5 | 11.10 | 1.41 | 41 % over Nyquist — WP impossible |

`σ_k = 1 / (2σ) = 0.10` Bohr⁻¹ at σ = 5 Bohr. The "3σ_k" is the standard outer extent of the Gaussian momentum distribution captured before aliasing. At dx = 0.40, the 600 eV run is the highest-energy WP that stays clean.

---

## 2. Physical parameters for r_s = 5.69 jellium

All quantities given in atomic units (a.u.) where length is in Bohr (a₀) and energy in Hartree (Ha), unless otherwise stated. 1 Ha = 27.2114 eV. Electron mass m_e = 1, ℏ = 1, |e| = 1.

| Quantity | Symbol | Value (a.u.) | Value (eV, Bohr) |
|---|---|---|---|
| Wigner-Seitz radius | r_s | 5.69 | — |
| Electron density | n = 3/(4π r_s³) | 1.295 × 10⁻³ | — |
| Fermi wavevector | k_F = (3π²n)^(1/3) | 0.337 | — |
| Fermi velocity | v_F = k_F | 0.337 | — |
| Fermi energy | E_F = k_F²/2 | 0.0568 | 1.55 eV |
| Plasma frequency | ω_p = √(4πn) | 0.1276 | 3.47 eV |
| Thomas-Fermi screening length | λ_TF ≈ 1/√(4 k_F/π) | 1.71 | — |
| Bohr velocity | v_Bohr | 1.000 | — |
| Box length (current simulations) | L | 50.0 | — |
| Number of jellium electrons (current sims) | N | 162 | — |
| Smallest |q| in box | q_min = 2π/L | 0.1257 | — |
| Plasmon Landau-damping cutoff | q_c ≈ ω_p/v_F | 0.379 | — |
| Maximum energy transfer in single collision (small v limit) | T_max ≈ 2 m_e v² | — | E_kin × 2 |

The smallest box momentum 0.126 a.u. is comparable to the plasmon-cutoff momentum 0.38 a.u. There are only about 118 modes with |q| < q_c in the cell. This is the dominant systematic uncertainty in the stopping-power signal — see §6 for the v-dependent box deficit.

---

## 3. Regime classification framework

### 3.1 Why two dimensionless parameters are sufficient

For a charged projectile passing through a homogeneous electron gas, the dimensional analysis says the stopping power can depend on at most four independent quantities: the projectile charge Z₁, its mass m_p, its velocity v, and the electron density n (which determines all host scales via k_F, v_F, ω_p). Forming dimensionless ratios:

- v/v_F — sets the *host response regime*. Below 1, the projectile is slower than the Fermi sea and only low-energy electron-hole pairs near E_F can be excited; the system responds quasi-statically. Above 1, the projectile outruns the Fermi sea and excites the full spectrum of e-h pairs plus collective plasmon modes. Above ~5, the projectile is so fast that the dipole approximation holds and the loss function is exhausted at I ≈ ℏω_p — the Bethe regime.
- κ = 2|Z₁|/v (in atomic units; equivalent to κ = 2|Z₁|e²/ℏv in SI) — sets the *projectile-scattering regime*. This is the Sommerfeld / Bohr parameter, the ratio of the Coulomb collision diameter b = 2|Z₁|e²/(m_e v²) to the projectile de Broglie wavelength λ_dB = ℏ/(m_e v). When κ ≪ 1, the projectile's wavefunction barely deflects on the scale of the Coulomb interaction and the Born approximation (perturbation theory in the projectile-electron interaction) holds. When κ ≫ 1, the projectile-electron collision must be treated by classical impact-parameter mechanics — quantum diffraction is irrelevant. When κ ~ 1, neither limit applies and the full Bloch interpolation formula or numerical methods are needed.

For a fixed projectile species and a fixed jellium target, κ and v/v_F together parameterize the configuration space of stopping experiments completely. Mass m_p enters only through the kinetic-energy / velocity conversion (E_kin = m_p v²/2) and doesn't appear in any cross section formula; the projectile's relevant property is its charge and its speed.

### 3.2 Why κ specifically

The Sommerfeld parameter κ is the *ratio* of two length scales that the projectile carries with it:

- Classical collision diameter `b = 2|Z₁ Z₂|e² / m_e v²` — the impact parameter at which a free Coulomb encounter classically scatters by 90°.
- Quantum de Broglie wavelength `λ_dB = ℏ / m_e v` — the spatial extent over which a quantum wavepacket of the projectile is uncertain in position.

`κ = b / λ_dB` (up to numerical factors of order 1).

When κ ≫ 1, the classical impact parameter is well-defined relative to the projectile's quantum uncertainty — the projectile follows a classical trajectory through a sequence of distinguishable encounters. The stopping cross section computed classically (Bohr) and quantum-mechanically (Bethe) agree to leading order in 1/v, but differ in their logarithmic correction factor:

- Bohr: `S = (4π Z₁² e⁴ n / m_e v²) × ln(C m_e v³ / |Z₁ e| ω₀)` with C = 2e^(-γ_Euler) ≈ 1.1229. The argument of the logarithm is the ratio of two impact parameters: the maximum impact parameter where energy transfer happens (Bohr's adiabatic cutoff `v/ω₀`) over the classical collision diameter `b`.
- Bethe: `S = (4π Z₁² e⁴ n / m_e v²) × ln(2 m_e v² / I)`. Here the argument is the ratio of two momentum transfers: the maximum quantum momentum transfer `2 m_e v` over the mean excitation momentum `m_e √(I/E)`.

When κ ≪ 1, the projectile is a delocalized wave on the scale of its interaction with the electron, and the Born approximation gives the Bethe form. When κ ~ 1, both pictures partially apply, and the Bloch correction interpolates:

`L_Bloch = ln(2 m_e v² / I) − Re[ψ(1 + iκ/2) − ψ(1)]`

where ψ is the digamma function. For κ → 0, L_Bloch → L_Bethe. For κ → ∞, L_Bloch → L_Bohr.

This is why κ is the right axis: it is exactly the parameter that interpolates the two limits of stopping theory. For an electron projectile (|Z₁| = 1) in jellium with n = 1.295 × 10⁻³, the line κ = 1 sits at v = 2 a.u., E = 54.4 eV. Below 54 eV, every simulation must reckon with the Born approximation's failure; above 54 eV, the classical and quantum scattering theories agree.

### 3.3 Why v/v_F specifically

The Fermi velocity v_F sets the scale at which Pauli blocking matters. For v < v_F, only electrons within energy ≲ v · v_F · m_e of the Fermi surface can absorb energy from the projectile, because all lower-momentum final states are occupied. The available phase space for excitation goes as `(v/v_F)²` at small v, and the stopping power becomes linear in v: `S(v) = Q(r_s, Z₁) · v` where Q is the friction coefficient.

For v > v_F, the projectile excites a fully thermalized spectrum. The full Lindhard dielectric function `ε(q, ω) = 1 − v_q χ⁰(q, ω)` with the Lindhard susceptibility χ⁰ captures the complete linear-response physics. The stopping integrand is the loss function Im[−1/ε(q, ω)], and its weight is split between two channels:

- The plasmon pole near ω ≈ ω_p (small q, coherent collective mode) — dominates at low v, gives a sharp absorption edge.
- The electron-hole continuum at ω = ε(k+q) − ε(k) for k within the Fermi sphere (all q, incoherent particle-hole pairs) — dominates at high v, gives a smooth contribution.

For v ≫ v_F, the plasmon dispersion ω_pl(q) ≈ ω_p[1 + (3/10)(qv_F/ω_p)²] is sampled below its Landau-damping cutoff for all relevant q (because q · v ≫ ω_p), and the integral collapses to the Bethe formula with I = ℏω_p (exact equipartition between collective and single-particle channels per Lindhard 1954, recovered by the f-sum rule).

### 3.4 The two-axis classification — graphical

Use this ASCII layout as a reference. Horizontal axis: log(v/v_F). Vertical axis: log(κ).

```
   κ                                                                  
   ↑                                                                  
50 |  ┌──────────┬────────────────┬───────────────────────┐           
   |  │          │                │                       │           
   |  │ Friction │   Bragg peak   │   Bohr classical      │           
   |  │ + quantum│   nonperturb.  │   regime              │           
   |  │ scattering   "phase-shift"│   (heavy ions, etc.)  │           
   |  │ (Pauli)  │   Echenique    │   κ≫1, Born fails     │           
   |  │          │   Pitarke      │   classical ε valid   │           
   |  │          │   Bloch        │                       │           
 1 |  ├─ Pauli ──┼─ Lindhard ─────┼──── Bloch ───────────-│ — κ=1 line─
   |  │ blocking │   full RPA     │   interpolation       │           
   |  │          │   dielectric   │   κ ~ 1               │           
   |  │ + κ small│   ε(q,ω)       │                       │           
   |  │ but no   │                │                       │           
   |  │ host     │ ω_p + ω_eh     │   Bethe formula:      │           
   |  │ classical│                │   S=(4πZ²n/mv²)ln(2v²/I)
   |  │ regime   │                │                       │           
0.05|  └──────────┴────────────────┴───────────────────────┘           
       0.3           1                5                30          100 → v/v_F
                                                                    
                  ←─── host response classical ───→                  
       host quantum                                                  
```

For the electron projectile (Z₁ = 1, m_p = m_e), the configuration space reduces to a *line* in this 2D plot: κ · (v/v_F) = 2 / v_F = 5.93. The line slopes down at 45° on log-log axes. The relevant energies along this line are:

| E (eV) | v (a.u.) | v/v_F | κ = 2/v | Regime (this work) |
|---|---|---|---|---|
| 1.55 | 0.337 | 1.00 | 5.93 | Pauli boundary, friction |
| 13.6 | 1.000 | 2.97 | 2.00 | Bragg peak, deep quantum scattering |
| 27 | 1.41 | 4.19 | 1.42 | Bragg-Bethe transition |
| 50 | 1.92 | 5.69 | 1.04 | Bethe onset, κ borderline |
| 100 | 2.71 | 8.04 | 0.74 | Lower Bethe, κ marginally classical |
| 300 | 4.69 | 13.9 | 0.43 | Mid Bethe |
| 600 | 6.64 | 19.7 | 0.30 | Mid Bethe |
| 1500 | 10.5 | 31.2 | 0.19 | Deep Bethe |
| 3000 | 14.85 | 44.1 | 0.135 | Deep Bethe, severe box artifact |

### 3.5 Why this 2D space is enough

For a heavier or higher-charge projectile (e.g. proton, Z₁ = 1 but m_p ≈ 1836; or He²⁺, Z₁ = 2), the same (v, κ) coordinates work — only the mapping from kinetic energy to velocity changes (E_kin = m_p v² / 2). A proton at 1500 eV has v = √(3000 / 1836) × 0.0214 = 0.244 a.u., not 10.5 a.u. — far below v_F. A 100 keV proton has v ≈ 2 a.u., comparable to a 54.4 eV electron. So the regime an experiment occupies is set by velocity, not by kinetic energy per se. The electron-projectile case studied here is special only in that velocity is high at modest kinetic energy, putting the deep-Bethe regime in reach with small E.

A heavier projectile (Z₁ > 1) shifts the κ values upward at the same v: for a Z₁ = 5 projectile, κ = 10/v in a.u., so the κ = 1 line crosses E = 2 × 54.4 × 25 = … the boundary moves up. This is why heavy ions in solids are firmly in the Bohr classical regime (κ ≫ 1) at most experimental velocities, while light projectiles transition through the Bloch regime.

---

## 4. The regimes in detail

For each regime, I give: the equation for stopping power, the assumptions baked in, the validity boundaries in (v/v_F, κ) space, and the failure modes when those assumptions are violated.

### 4.1 Friction regime (v < v_F, any κ)

**Equation.** `S(v) = Q(r_s, Z₁) × v`, where Q is the friction coefficient. To lowest order in Z₁:

`Q ≈ (4π² Z₁² n) / (m_e v_F³ k_F²) × ln(1 + (2 k_F / k_TF)²)`

For r_s = 5.69 jellium and Z₁ = 1, the nonlinear Echenique-Nieminen-Ritchie (ENRA) DFT calculation gives Q ≈ 0.5–1.5 a.u. depending on the sign of the charge (Barkas-Andersen Z₁³ effect).

**Assumptions.** (i) v ≪ v_F so Pauli blocking restricts excitations to a thin shell around E_F. (ii) Linear response (or nonlinear DFT) is valid for the *static* screened projectile potential. (iii) The projectile is treated as a static scattering center with screened phase shifts {δ_l} evaluated at k = k_F. The friction coefficient is then `Q = (m_e v_F² / π) × ∑_l (l+1) sin²(δ_l − δ_{l+1})` — Lindhard's friction sum.

**Validity.** v/v_F < 0.5–1. Above v/v_F ≈ 1, the projectile starts to outrun the screening cloud and the Pauli-blocking approximation breaks down. Below ~0.1 v_F, the projectile is so slow that lattice (or for jellium, density-wave) effects can matter — but these are absent in a uniform jellium.

**Failure modes.** Underestimates S at v ≈ v_F (where the friction line and the Bragg peak intersect). Misses the entire plasmon contribution because plasmons need v · q ~ ω_p, which requires v at least of order v_F.

**Theory pointers.** Echenique, Nieminen, Ritchie, *Solid State Commun.* 37, 779 (1981); Puska & Nieminen, *PRB* 27, 6121 (1983); Echenique, Pitarke, Chulkov, Rubio, *Chem. Phys.* 251, 1 (2000).

### 4.2 Bragg peak (1 ≲ v/v_F ≲ 5, any κ)

**Equation.** No closed form. Stopping power is found from the full Lindhard linear-response formula:

`S(v) = (2 Z₁² / π v²) × ∫₀^{q_max} (dq/q) × ∫₀^{qv} dω × ω × Im[−1/ε(q, ω)]`

with q_max = 2 m_e v + k_F (kinematic limit) and ε = ε_RPA or ε_RPA-with-local-field-correction. The peak of S(v) sits near v ≈ 1.5 v_F at S_peak ≈ 0.3–1 eV/Bohr for r_s = 5.69 (this is a literature estimate; the present simulation campaign should refine it).

**Assumptions.** (i) Linear response in the projectile-electron interaction — i.e., the induced density is proportional to the bare projectile potential. (ii) RPA for the dielectric function (or with a static local-field correction G(q) like Corradini et al.). (iii) No memory in xc: adiabatic LDA or static xc kernel.

**Validity.** Whenever the dielectric formalism converges — formally any v, but in practice the Lindhard formula loses accuracy when κ is large (perturbation breaks down) and when v ≪ v_F (nonlinear screening dominates and ENRA is more accurate).

**Failure modes.** Underestimates S at large κ (need Bloch corrections or beyond-Born scattering). Misses Barkas-Andersen Z₁³ effect (asymmetry between Z₁ > 0 and Z₁ < 0). At low v, RPA gives the wrong friction coefficient by a factor of 2–3 because it omits nonlinear screening.

**Theory pointers.** Lindhard, *Mat. Fys. Medd. Dan.* 28, no. 8 (1954); Mahan *Many-Particle Physics* ch. 5; Giuliani & Vignale ch. 4.

### 4.3 Bethe regime (v/v_F ≳ 5, κ ≲ 1)

**Equation.** `S(v) = (4π Z₁² e⁴ n / m_e v²) × L(v)` with the stopping number

`L(v) = ln(2 m_e v² / I) − v²_F / v² × <shell correction> − v²_F / (2 v²) × <Lindhard −1/2 correction> − higher-order Z₁ terms`

For pure jellium (no atomic shells), the mean excitation energy is `I = ℏω_p` (Lindhard's equipartition theorem), and the Lindhard −1/2 correction is the leading O(v_F²/v²) term:

`L_jellium(v) = ln(2 v² / ω_p) − 1/2 + O(v_F⁴/v⁴)`

The minus-one-half is the next-to-leading-log correction from the kinematic integration limits at q → 2 m_e v. It pulls S down by about 10% in the Bethe regime.

**Assumptions.** (i) First Born approximation for projectile-electron scattering (κ ≪ 1). (ii) Dipole approximation for the response of the host — valid when the projectile's deflection in any single collision is small compared to its velocity, which holds for v ≫ v_F. (iii) Plane-wave projectile (delocalized over the system).

**Validity.** v/v_F ≳ 5 and κ ≲ 0.5. The closer to those boundaries, the more important higher-order corrections become.

**Failure modes.** Diverges at v → 0 (since log argument goes negative — the formula gives unphysical negative stopping). Underestimates S at high κ — needs Bohr or Bloch corrections. Misses Barkas Z₁³ effect.

**For r_s = 5.69 jellium, numerical Bethe-Lindhard predictions:**

| E (eV) | v (a.u.) | 2v²/ω_p | L = ln − 1/2 | S (Ha/a₀) | S (eV/Bohr) |
|---|---|---|---|---|---|
| 50 | 1.92 | 57.8 | 3.56 | 0.01566 | 0.426 |
| 100 | 2.71 | 115 | 4.24 | 0.00940 | 0.256 |
| 300 | 4.69 | 346 | 5.35 | 0.00395 | 0.107 |
| 600 | 6.64 | 691 | 6.04 | 0.00223 | 0.0606 |
| 1500 | 10.5 | 1728 | 6.95 | 0.00103 | 0.0279 |
| 3000 | 14.85 | 3456 | 7.65 | 5.65 × 10⁻⁴ | 0.0154 |

### 4.4 Bohr classical regime (κ ≫ 1)

**Equation.** `S(v) = (4π Z₁² e⁴ n / m_e v²) × ln(C m_e v³ / |Z₁| e² ω₀)` with C = 2 e^(−γ) = 1.1229 and ω₀ = ℏω_p for jellium.

**Assumptions.** (i) Classical orbit picture for the projectile-electron Coulomb encounter. (ii) Target electrons modeled as classical harmonic oscillators with resonance ω₀. (iii) Long-range cutoff at Bohr's adiabatic radius v/ω₀. (iv) Short-range cutoff at the classical collision diameter b = 2|Z₁|e²/(m_e v²).

**Validity.** κ ≫ 1, which for unit-charge projectile means v ≲ 1 a.u. = 27 eV. Most light-projectile experiments at relevant energies do NOT lie here. Heavy ions (Z₁ ~ 10–90) lie deep in this regime.

**Failure modes.** At κ ~ 1 the classical picture loses validity. Use Bloch interpolation. At low v (entering friction regime) Pauli blocking takes over.

### 4.5 Bloch interpolation (κ ~ 1)

**Equation.** Replace the stopping-number logarithm by

`L_Bloch(v) = ln(2 m_e v² / I) − Re[ψ(1 + iκ/2) − ψ(1)]`

For κ → 0, ψ → ψ(1) = −γ, so L_Bloch → L_Bethe. For κ → ∞, ψ → ln(iκ/2), so L_Bloch → L_Bohr.

**Numerical values of the Bloch correction `−Re[ψ(1+iκ/2) − ψ(1)]`** (subtracted from L_Bethe):

| κ | Correction (negative) |
|---|---|
| 0.1 | −0.0016 |
| 0.3 | −0.0144 |
| 0.5 | −0.0395 |
| 0.74 | −0.0850 |
| 1.0 | −0.150 |
| 1.5 | −0.317 |
| 2.0 | −0.526 |

So at the 100 eV point (κ = 0.74), the Bloch correction reduces L by 0.085 out of L_Bethe ≈ 4.24, a 2% reduction — small but not negligible. At 50 eV (κ = 1.04), the correction is 0.16 out of 3.56, a 4.5% reduction. This is the *expected size of the κ-driven deviation* you should look for in the 50 eV simulation.

### 4.6 Beyond perturbation theory

For *any* combination of (v/v_F, κ) where the perturbation expansion in Z₁ breaks down, one needs:

- Nonlinear DFT (static screening): Echenique-Nieminen-Ritchie, gives the proper Q at low v.
- Full quantum scattering (phase shifts on the screened potential): valid at all v, gives the proper transport cross section.
- Real-time TDDFT (this project): gives the *dynamical* nonlinear response — captures wake formation, plasmon emission, screening dynamics, all in one framework. This is the gold-standard approach for the Bragg-peak region and is what the present simulation campaign is designed to provide.

---

## 5. Existing simulations — regime classification

Five runs exist, all on the same `r_s = 5.69, L = 50³ Bohr, N = 162, dx = 0.40, dt = 0.020` cell. (The `run_classical_e1000_L40x40x150` directory exists but never executed; exclude it from analysis.)

| Path | Projectile type | E (eV) | v (a.u.) | v/v_F | κ | k₀σ | Regime |
|---|---|---|---|---|---|---|---|
| `run_classical_n162_L50_E100` | classical (modified PP + Ehrenfest) | 100 | 2.71 | 8.04 | 0.74 | — | lower Bethe, κ borderline |
| `run_classical_n162_L50_E600` | classical | 600 | 6.64 | 19.7 | 0.30 | — | mid Bethe, κ small |
| `run_classical_e1500_L50_cubic` | classical | 1500 | 10.5 | 31.2 | 0.19 | — | deep Bethe, κ small |
| `run_wp_e100_*` | wave packet | 100 | 2.71 | 8.04 | 0.74 | 13.6 | lower Bethe, classical-packet limit |
| `run_wp_e600_*` | wave packet | 600 | 6.64 | 19.7 | 0.30 | 33.2 | mid Bethe, classical-packet limit |

All five WP and classical projectiles launch from `(0, 0, −10)` Bohr and propagate to `z = +25` Bohr — the +z box face — at their initial velocity. The simulation `N_STEPS` is set per-run from `ceil(35 / v_initial / dt) = ceil(1750 / v_initial)` at dt = 0.020 a.u.

**WP-1500 eV is not feasible at this grid** (see §1.2 Nyquist). The classical-1500 eV run will remain unpaired for the comparison plot. This is a hardware-driven asymmetry, not a planning choice.

All five runs are firmly in the *Bethe regime* for host response (v/v_F ≥ 8). All five are *marginally to deeply classical* in projectile scattering (κ ≤ 0.74). None probes the Bragg peak, none probes the friction regime, none unambiguously enters the κ > 1 quantum-scattering regime.

The 1500 eV classical run has the cleanest re-analyzed stopping power signal: `S = 0.0249 ± 0.0002 eV/Bohr` (from a windowed fit in t ∈ [1.5, 3.5] a.u., which avoids both the initial bath-relaxation transient and the periodic-image wake-overlap region near t = 4.0). This is 89% of the Bethe-Lindhard prediction 0.0279 eV/Bohr — a 11% box-truncation deficit, consistent with the q_min · v / ω_p ≈ 10.5 cutoff ratio.

**The agent must redo this windowed analysis for the 100 and 600 eV classical runs.** The original linear-fit-over-full-run analysis is contaminated by:

1. *Initial transient* — first ~0.5 a.u. of bath relaxation as the projectile's screening cloud establishes itself. Excludes the first ~5 Bohr of projectile travel.
2. *Periodic-image wake-overlap* — last fraction of the run, when the projectile sits within ~50 Bohr of its own periodic image. The exact end-of-clean-window depends on v: at v = 10.5 it's around t = 3.5 (z = 36.75 Bohr); at v = 6.64 (600 eV) it would be around t ≈ 5.5 (z ≈ 37 Bohr) if the run goes that far; at v = 2.71 (100 eV) the periodic-image overlap is a much smaller effect because the wake wavelength λ_wake = 2π v / ω_p = 133 Bohr is only 2.7× the box, vs 515 Bohr at 1500 eV.

The agent should produce a *windowed stopping power* for each existing run with proper bounds and error bars.

---

## 6. Operational definitions of stopping power for the postprocess

### 6.1 Classical-projectile stopping power

Use the *bath energy gain method*. By energy conservation, in a closed simulation with the projectile driven externally at fixed velocity (Ehrenfest with custom-mass UPF), the total Kohn-Sham energy `energy_total` from `observables.csv` rises at a rate equal to the negative of the work done on the projectile per unit time:

`d E_bath / d t = − d E_projectile / d t = F_drag · v`

So `S(v) = (1/v) × d E_bath / d t`, evaluated as the slope of `delta_E_bath = energy_total(t) − energy_total(0)` against `t` in a clean window.

**Launch geometry and window construction.** From the config files (`electron_proj_E*_L50_cubic.hpp`), all classical projectiles launch at `z₀ = −10` Bohr in the centred-Cartesian cell `r ∈ [−25, +25]`. They propagate in +z at fixed initial speed `v = √(2 · E_eV / 27.2114)` Bohr/atu. The simulation ends when the centroid reaches `z = +25` Bohr, so `t_end = 35 / v`. Project displacement at any time is `Δz(t) = v · t`, ranging from 0 at launch to 35 Bohr at the end of the run.

**Recommended windowing rule:**

1. `t_start`: time at which `Δz` first reaches ~3 Bohr (clear of the initial wake-formation transient). Concretely `t_start ≈ 3 / v`.
2. `t_end_clean`: time at which `Δz` reaches 28 Bohr (the projectile is then 17 Bohr ahead of its launch and 7 Bohr from the +z box face). At this point the projectile is still 7+ Bohr from the periodic image of its launch position, so the wake has minimal self-interaction. Concretely `t_end_clean ≈ 28 / v`.
3. Linear-regression `delta_E_bath` vs `time_au` on `t ∈ [t_start, t_end_clean]`. Report slope ± standard error.
4. `S(v) = (slope_in_Ha_per_atu × 27.2114 / v)` in eV/Bohr.

For the 1500 eV run, `t_start = 0.286` a.u. and `t_end_clean = 2.667` a.u. give roughly the same window I used by hand in the empirical re-analysis (which found S = 0.0249 eV/Bohr over t ∈ [1.5, 3.5]). The user-supplied 1500 eV config also has `N_STEPS = 860` at `dt = 0.005`, giving `t_total = 4.30` a.u. — the last 1.6 a.u. (z ∈ [28, 45] Bohr) is where the periodic-image-wake overlap progressively dominates and the slope falls off.

For the 100 eV existing run, `t_total = 12.92` a.u. (from `N_STEPS = 646` at `dt = 0.020`), `t_start ≈ 1.1` a.u., `t_end_clean ≈ 10.3` a.u. — a 9.2 a.u. window, much wider than the 1500 eV case, so the linear fit should be statistically very tight. But note: at v = 2.71 and λ_wake = 2π v / ω_p = 133 Bohr (vs box = 50 Bohr), the wake is wrapped onto itself 2.7 times within the box — the "image overlap" effect is qualitatively different from high-v cases, and the running slope may show plasmon-period oscillations rather than monotonic degradation.

### 6.2 Wave-packet stopping power

The wave packet is one Kohn-Sham orbital initialized as `φ_proj(r, t=0) = (2πσ²)^(−3/4) exp(−|r − r₀|² / 4σ²) exp(i k₀ · r)` with σ = 5.0 Bohr, k₀ = (0, 0, k₀ᶻ), and `r₀ = (0, 0, −10)` Bohr. It evolves under the time-dependent Kohn-Sham Hamiltonian alongside the 162 bath orbitals. Its centroid `⟨z_proj⟩(t) = ⟨φ_proj | ẑ | φ_proj⟩` is not constrained to follow a fixed velocity — it evolves self-consistently and slows down as energy transfers to the bath.

**Note on the regime:** At every energy in the proposed sweep, `k₀σ ≥ 9.6` (see §1, Q2). The packet is therefore essentially a classical localized projectile, not a quantum-mechanical few-wavelength wave packet. The WP-vs-classical comparison is therefore *not* testing projectile-wavefunction-spread effects (those are tiny by design). It is testing two host-side differences:

- **Indistinguishability with bath electrons.** The WP electron is the same fermion species as the 162 bath electrons, so the full Slater determinant must antisymmetrize between them. The classical-projectile sim uses a separately-specified pseudo-hydrogen with custom mass `m = m_e` — antisymmetrization with bath electrons does not apply because the projectile is treated as a distinct particle species.
- **Finite spatial extent of the source.** A σ = 5.0 Bohr Gaussian creates a wake that differs from a point charge's wake in the small-q (long-wavelength) region of Im[−1/ε(q,ω)]. Specifically, the source spectrum has a Gaussian cutoff at `|q| ~ 1/σ = 0.2` Bohr⁻¹, comparable to `q_min = 0.126` Bohr⁻¹ of the box. The WP simulation therefore couples to only the few `|q| < 0.2` Bohr⁻¹ modes that the box supports, while the classical-projectile sim's point charge couples to all modes up to dx-grid cutoff.

**Primary operational definition (bath-energy method):**

`S_WP(v_initial) ≡ d E_bath(t) / d ⟨z_proj⟩(t)`

evaluated as a parametric slope in `t`. This is the *direct analog* of the classical definition — bath energy gain per unit projectile displacement — and is robust to packet dispersion or partial reflection because it only requires `⟨z_proj⟩(t)` to be monotonic (not constant in velocity).

**Secondary operational definition (kinetic-energy method):**

`S_WP,kin(v_centroid) ≡ − d ⟨T_proj⟩(t) / d ⟨z_proj⟩(t)`

where `⟨T_proj⟩(t) = ⟨φ_proj | −∇²/2 | φ_proj⟩` and `v_centroid(t) = d ⟨z_proj⟩ / dt`. In the classical-packet limit (k₀σ ≫ 1, no dispersion, no reflection), `S_WP,bath` and `S_WP,kin` should agree to leading order.

The agent must implement the following for each WP run:

1. Extract `⟨z_proj⟩(t)` from the saved wavefunction snapshots (cadence: every `WRITE_EVERY = 2` steps, so 0.040 a.u. between saves).
2. Extract `⟨T_proj⟩(t)` from the same snapshots via the kinetic-energy operator.
3. Extract `σ_z(t) = √(⟨z²⟩ − ⟨z⟩²)` and `σ_x(t)`, `σ_y(t)` as packet-spread diagnostics. Initial value is σ = 5.0 Bohr.
4. Compute `S_WP,bath` per the primary definition and `S_WP,kin` per the secondary.
5. Apply the *same windowing rule* as the classical case: clip the early ~3-Bohr-traversal transient and the late ~7-Bohr-from-boundary wake-overlap window. Convert `Δ⟨z_proj⟩(t)` to time bounds.
6. Report both `S_WP,bath` and `S_WP,kin` with error bars; flag any discrepancy > 5% as a dispersion / reflection signal.

**Required new postprocess module** (suggested name `inqview/postprocess/wavepacket_observables.py`). For each WP run, emit a `wp_observables.csv` with columns `step, time_au, z_centroid, x_centroid, y_centroid, T_expectation, sigma_z, sigma_x, sigma_y, packet_norm`. The `packet_norm = ⟨φ_proj | φ_proj⟩` is a sanity-check column — should stay at 1.0 to machine precision throughout, and any drift signals numerical issues.

### 6.3 Cross-comparison metric

For each energy E with both a classical and a WP run, define:

`Δ(E) ≡ S_WP,bath(E) / S_classical(E) − 1`

The plot of Δ(E) vs E (or vs v/v_F or κ) is the central scientific result of this comparison.

**What Δ(E) measures in this campaign, given σ = 5 fixed:**

- Δ → 0 in the high-v limit (large k₀σ, where exchange and source-extent effects are smallest compared to the dominant Bethe stopping mechanism). At 600 eV (k₀σ = 33.2), Δ should be very small — a percent or two at most.
- Δ becomes non-zero at low v due to:
  - *Exchange with bath* — the WP electron's Pauli interaction with the bath wavefunctions reduces the available phase space for excitation; the classical projectile does not have this constraint.
  - *Wake spectrum cutoff* — the source's Gaussian momentum cutoff at `q ~ 0.2` Bohr⁻¹ removes coupling to plasmon-like modes at `q > 0.2`. At low v this matters more because the Bethe integrand sits at smaller q.
  - *Reflection / packet structure* — at very low v, even a wide packet partly back-scatters from the bath response. Captured by `σ_z(t)` and `packet_norm(t)` diagnostics, not by Δ directly.
- The classical-100-eV vs WP-100-eV pair gives a baseline for Δ in the lower-Bethe regime (v/v_F = 8, κ = 0.74).
- The classical-50-eV vs WP-50-eV pair tests Δ in the Bragg-peak-onset regime (v/v_F = 5.7, κ = 1.04), where both κ-driven (Bloch) and source-cutoff-driven (q-truncation) effects should be largest.

The sign of Δ is informative. If WP < classical (Δ < 0), the most likely cause is the wake-source spectrum cutoff (the WP couples to fewer host modes than a point source). If WP > classical (Δ > 0), exchange effects with bath are enhancing energy loss — unusual but physically possible.

---

## 7. Postprocess deliverables

The agent should produce a single combined-analysis output for the campaign. Minimum content:

1. `S_table.csv` with columns: `run_path, projectile_type, E_eV, v_au, v_over_vF, kappa, S_measured (eV/Bohr), S_uncertainty (eV/Bohr), S_Bethe_pure, S_Bethe_Lindhard, S_Bloch_corrected, box_deficit_fraction, window_t_start, window_t_end, window_z_start, window_z_end, packet_sigma_initial, packet_sigma_final, v_centroid_initial, v_centroid_final`. For classical runs, the packet-σ columns are NaN; for WP runs, all classical-only columns are NaN as appropriate.
2. `S_vs_E_plot.pdf` (or PNG): a log-log plot of S(E) with the Bethe-Lindhard curve, Bloch-corrected curve, all measured classical points, all measured WP points, and the proposed new points. Use the `chart`-module conventions from this conversation (eV/Bohr for S, eV for E).
3. `delta_vs_E_plot.pdf`: a log-linear plot of `Δ(E) = S_WP/S_classical − 1` vs E, showing the classical-to-WP deviation as a function of energy.
4. `box_deficit_diagnostic.pdf`: for the 1500 eV run, plot the running stopping-power estimate `dE_bath/dz` smoothed over rolling 0.4 a.u. windows, vs `z_proj`, to visualize the periodic-image-wake-overlap effect. Repeat for 100 and 600 eV runs once their CSVs are extracted.

---

## 8. Recommended new simulations under 9-hour, 2-GPU budget

**Constraints:**
- Same density (r_s = 5.69) and box (50³ Bohr) as existing runs.
- 9 hours of wall-clock time.
- One GPU dedicated to classical-projectile simulations, the other to WP simulations.
- Same N = 162 electron count.

## 8. Recommended new simulations under 9-hour, 2-GPU budget

**Constraints (from the user):**
- Same density (r_s = 5.69) and box (50³ Bohr) as existing runs.
- Same N = 162, same dx = 0.40 Bohr, same dt = 0.020 a.u.
- Same launch geometry: z₀ = −10 Bohr, traversal to z = +25 Bohr.
- 9 hours of wall-clock time.
- One GPU dedicated to classical-projectile simulations, the other to WP simulations.
- WP-1500 eV impossible at this dx (see §1.2); deep-Bethe regime stays unpaired.

### 8.1 Cost model

The user-supplied configs make N_STEPS deterministic:

`N_STEPS = ceil(35 / v_initial / dt) = ceil(1750 / v_initial)` at dt = 0.020 a.u.

The empirical per-step cost from existing runs (omitting the unusually-fast 1500 eV which ran at dt = 0.005 and used a longer-than-single-pass `N_STEPS = 860`):

| Run | Type | N_STEPS | Wall time | s / step |
|---|---|---|---|---|
| classical 100 eV | classical | 646 | 198 min | 18.4 |
| classical 600 eV | classical | 264 | 63 min | 14.3 |
| WP 100 eV | wave packet | 646 | 110 min | 10.2 |
| WP 600 eV | wave packet | 264 | 57 min | 13.0 |

Classical per-step cost is mildly higher at lower v (likely due to slower self-consistent convergence in the slower-projectile regime), but typically 13–20 s. WP per-step cost is 10–13 s. The agent should be conservative and budget `classical: 18 s/step, WP: 12 s/step` for new runs.

For the proposed energies:

| E (eV) | v (a.u.) | N_STEPS | classical wall time | WP wall time |
|---|---|---|---|---|
| 50 | 1.92 | 912 | ~4.6 h | ~3.0 h |
| 300 | 4.69 | 374 | ~1.9 h | ~1.3 h |
| 25 (optional stretch) | 1.36 | 1287 | ~6.4 h | ~4.3 h |

### 8.2 GPU 0 (classical) — recommended pair

**Sim 1: classical 300 eV.** Estimated cost: ~1.9 h.

- Config: clone `electron_proj_E600_L50_cubic.hpp` to `electron_proj_E300_L50_cubic.hpp`. Change `WP_EKIN_EV = 300.0` and re-derive `N_STEPS = 374`. Everything else identical.
- Projectile velocity: v = 4.69 a.u.
- Regime: mid Bethe, κ = 0.43, v/v_F = 13.9.
- Bethe-Lindhard prediction: S = 0.107 eV/Bohr.
- Expected measured: S ≈ 0.094–0.100 eV/Bohr after ~10 % box deficit.
- Scientific purpose: fills the gap between the 100 eV and 600 eV classical points. With four classical Bethe-regime points (100, 300, 600, 1500), the agent can fit `S × v² = A × [ln(2v²/I_eff) − 1/2]` with A and I_eff free. Expected `I_eff ≈ ℏω_p = 3.47 eV` if Bethe holds; any systematic departure quantifies finite-box and finite-time bias as a function of v.

**Sim 2: classical 50 eV.** Estimated cost: ~4.6 h.

- Config: clone same template to `electron_proj_E50_L50_cubic.hpp`. Change `WP_EKIN_EV = 50.0`, `N_STEPS = 912`. Verify the Nyquist check still passes (k_max = 2.52 Bohr⁻¹ vs k_Nyquist = 7.85; 32 % usage, very clean).
- Projectile velocity: v = 1.92 a.u.
- Regime: top of Bragg peak, κ = 1.04, v/v_F = 5.69.
- Bethe-Lindhard prediction: S = 0.426 eV/Bohr. With Bloch correction at κ = 1.04: S = 0.407 eV/Bohr.
- Expected measured: 0.32–0.41 eV/Bohr depending on the strength of Bragg-peak vs Bethe-asymptote deviation.
- Scientific purpose: the quantum-dominated point. Even the classical-projectile method shows deviations from Bethe here because the *host response itself* is in the Bragg-peak regime where the full Lindhard dielectric matters, not the asymptotic ln-dependence. The deviation `(S_measured − S_Bethe_Lindhard) / S_Bethe_Lindhard` at this energy quantifies the size of beyond-asymptotic corrections (Bloch, Bragg-peak Lindhard).

**Total GPU 0 wall time: ~6.5 h.** ~2.5 h margin against the 9 h budget.

### 8.3 GPU 1 (wave packet) — recommended pair

**Sim 3: WP 300 eV.** Estimated cost: ~1.3 h.

- Config: same 300 eV file as Sim 1, but instantiate `Electron_Proj_E300_L50_cubic_WP_dx0p40` struct (with `WP_ENABLED = true`, `N_IONS = 0`). σ = 5 Bohr inherited.
- Scientific purpose: matched WP partner for classical 300 eV. With WP-100, WP-300, WP-600 the agent has three matched (classical, WP) Bethe-regime pairs for the Δ(E) plot.
- Expected: Δ(300) ≈ 0 within the percent-level statistical noise (k₀σ = 23.5, deep in the classical-packet limit, exchange effects suppressed by phase-space factor).

**Sim 4: WP 50 eV.** Estimated cost: ~3.0 h.

- Config: 50 eV file as Sim 2, WP struct.
- Scientific purpose: the matched WP partner for the most interesting classical point. At κ = 1.04 and v/v_F = 5.69, this is where both κ-driven (Bloch / Born) and host-side (exchange, source-extent) deviations between WP and classical are largest.
- Diagnostics to watch: σ_z(t) trajectory (initial 5 Bohr; predicted final ≤ 5.3 Bohr from dispersion physics; any major spread signals real packet structure or back-scattering); `packet_norm(t)` (should stay = 1.0 to machine precision).
- Cost note: even though k₀σ = 9.6 is in the classical-packet limit, the packet's spatial footprint (3σ ≈ 15 Bohr) is comparable to the launch-offset (z₀ = −10 Bohr means the leading edge of the packet is already at z = +5 at t = 0). The early-time bath-energy signal will reflect the packet's pre-existing overlap with the bath. Use `t_start ≈ 1.6` a.u. (Δz = 3 Bohr) for the windowed slope.

**Total GPU 1 wall time: ~4.3 h.** ~4.7 h margin against the 9 h budget.

### 8.4 Optional Sim 5 (GPU 1 surplus): WP 25 eV (stretch)

If GPU 1 finishes the matched pair early, an unmatched WP 25 eV can extend the quantum-regime sweep into the Bragg peak proper.

- v = 1.36 a.u., v/v_F = 4.0, κ = 1.47.
- N_STEPS = 1287, estimated cost ~4.3 h on GPU 1. Possible after WP 300 + WP 50 = 4.3 h.
- Combined GPU 1 total: 8.6 h. Just fits.
- Caveats: no matched classical 25 eV (would cost ~6.4 h on GPU 0, blowing the budget). The WP-25 measurement therefore stands alone as a "reach" point in the κ > 1 regime, not as a paired comparison.
- Scientific value: tests whether the WP-vs-Bethe-Lindhard deviation grows monotonically with κ across the κ = 1 boundary, useful for tuning future classical sims at the same energy.

### 8.5 What NOT to do

- **Do not attempt WP at 1500 eV.** The Nyquist diagnostic gives a 41 % aliasing margin at dx = 0.40 (see §1.2). The result would be unphysical regardless of how much wall time is invested.
- **Do not attempt classical at 25 eV in the current budget.** It would cost ~6.4 h on GPU 0, leaving no room for the 50 eV sim that the user has flagged as the priority quantum-dominated point. If a future budget allows, the 25 eV pair makes a clean κ = 1.5 anchor point.
- **Do not change dt, dx, σ, N, L, or z₀ between any of these runs.** Cross-comparison only works if every cell/projectile/grid parameter is held fixed and only `WP_EKIN_EV` (and the derived `N_STEPS`) varies. The config-inheritance pattern in `Common_E*_L50_cubic` enforces this; the agent should reuse it verbatim.

### 8.6 Budget summary

| GPU | Sim | Energy | Wall time est. |
|---|---|---|---|
| 0 (classical) | classical 300 eV | new | 1.9 h |
| 0 (classical) | classical 50 eV | new | 4.6 h |
| **GPU 0 total** | | | **6.5 h** |
| 1 (wave packet) | WP 300 eV | new | 1.3 h |
| 1 (wave packet) | WP 50 eV | new | 3.0 h |
| 1 (wave packet, optional) | WP 25 eV | new | 4.3 h |
| **GPU 1 total (with optional)** | | | **8.6 h** |
| **GPU 1 total (without optional)** | | | **4.3 h** |

The minimum two-pair plan (300 eV + 50 eV, both projectile types) costs 6.5 h on the cost-driving GPU (classical), well within 9 h. The agent has 2.5 h of safety margin on GPU 0 and 4.7 h on GPU 1.

---

## 9. Comparison plan

After all sims complete, the agent should produce the following final analysis package:

### 9.1 Master stopping-power plot

Log-log axes: E (eV) horizontal from 30 to 2000, S (eV/Bohr) vertical from 0.01 to 1. Show:

- Bethe-pure curve (solid line).
- Bethe-Lindhard curve (dashed).
- Bloch-corrected curve (dotted).
- Measured classical points: 100, 300, 600, 1500 eV with error bars.
- Measured WP points: WP-100, WP-300, WP-600, WP-50 (and WP-1500 if available) with error bars.
- The Bragg-peak indicator (a shaded band 1.5 v_F ≤ v ≤ 5 v_F).

### 9.2 Δ(E) plot

Linear vertical axis, log horizontal axis. Show `Δ(E) = S_WP(E)/S_classical(E) − 1` with error bars at all matched-pair energies. Theoretical guides:

- The pure-Bloch correction line: predicts the Δ from κ alone, ignoring projectile-wavefunction effects.
- A horizontal zero line: classical-limit prediction.
- Markers labeled with κ values to help interpret which deviations are κ-driven vs σ-driven.

### 9.3 Box-deficit diagnostic

For each classical run, plot `dE_bath/dt` smoothed over 0.4 a.u. windows, vs `z_proj`. Identify the clean window. Compute and report the *windowed* stopping power and the *full-range* stopping power separately. The difference is the systematic from periodic-image wake overlap.

### 9.4 Quantitative tests of the regime hypotheses

**Test 1 — Bethe scaling.** With four classical points in the Bethe regime (100, 300, 600, 1500 eV), fit:

`S(v) × v² = A × [ln(2 v² / I_eff) − 1/2]`

with `A` and `I_eff` free. Expected: `A ≈ 4π n × Z₁² = 0.01627` a.u., `I_eff ≈ ℏω_p = 0.1276` Ha = 3.47 eV. Any departure quantifies systematic finite-box and finite-time effects.

**Test 2 — Bloch correction visible at 50 eV.** Compare classical-50 eV measurement to Bethe-Lindhard prediction. The classical-projectile run uses Ehrenfest, which implicitly captures any Born-vs-Bohr-mixed projectile-electron physics through the full real-time Coulomb dynamics. Difference from Bethe-Lindhard signals the size of higher-order corrections at κ ≈ 1.

**Test 3 — Δ(E) vanishes in the classical limit.** Plot Δ(E) for the matched pairs. If Δ → 0 at high E (small κ, large k₀σ), the WP framework reproduces classical physics in the classical limit — a basic consistency check.

**Test 4 — Δ(E) becomes large near κ = 1.** Δ(50 eV) is expected to be significantly larger in magnitude than Δ(100, 300, 600 eV). The sign of Δ(50) is informative: a negative Δ means WP physics reduces stopping below classical (consistent with packet partially passing through coherently); a positive Δ means WP physics enhances it (consistent with packet exchange / Pauli with bath electrons that the classical sim cannot capture).

---

## 10. Assumptions and caveats

1. **The classical projectile velocity is held fixed by the Ehrenfest framework with custom-mass UPF** (`PROJ_MASS_AMU = 1.0 / 1822.8885` → m = m_e). With the user's setup, the projectile *does* respond to forces from the bath, so v decreases slightly during the run (~1 % over a full traversal, per the config comments). The bath-energy method remains valid because it identifies `dE_bath/dt = F_drag · v_actual(t)`, not `F_drag · v_initial`. The agent should optionally extract `v_actual(t)` from the projectile trajectory (if logged) to refine the S calculation; otherwise using `v_initial` introduces at most a 1 % error.

2. **The wave-packet width σ = 5.0 Bohr is fixed by `Base_N162_L50_E1p5` across the entire velocity sweep.** This is a deliberate choice that puts every WP run deep in the classical-packet limit (k₀σ ≥ 9.6). The WP-vs-classical comparison therefore probes host-side physics (exchange, source extent), not projectile-wavefunction structure. Any future campaign that wants to test projectile-wavefunction effects must vary σ as well, ideally taking σ × k₀ down toward 1.

3. **The classical-projectile-with-Ehrenfest method is a quantum mechanical method on the host side.** It solves TDDFT for the 162 bath electrons; only the projectile is treated as a classical point charge with a custom pseudopotential. Both classical and WP runs share the same host treatment. The "classical vs quantum" comparison in this campaign is between two different projectile descriptions of the *same* host TDDFT calculation.

4. **The Bethe formula assumes infinite homogeneous jellium.** The finite 50 Bohr box truncates the loss-function integral at q_min = 0.126 a.u., losing weight in the small-q (long-wavelength plasmon) channel. Empirically the 1500 eV run shows an 11 % deficit relative to Bethe-Lindhard. Scaling: deficit grows roughly with `q_min · v / ω_p`. For the proposed runs:
   - 50 eV: q_min · v / ω_p = 1.92 → expected deficit ≤ 5 %.
   - 300 eV: 4.7 → expected deficit ≈ 9–12 %.
   - The agent should use the 300 eV measurement to refine the deficit-vs-v scaling and apply it as a systematic correction to all Bethe-regime points before extracting `I_eff`.

5. **The Lindhard −1/2 correction is the leading next-to-leading-log term**, valid to O(v_F²/v²). At v/v_F = 5.7 (50 eV), the next-order correction is `(v_F/v)² ≈ 0.031`, so the Lindhard formula has residual ~3 % error from neglected higher-order terms. At 1500 eV this is < 0.001 and negligible.

6. **The Bloch correction in §4.5 is computed for *isolated* atoms** with discrete excitation spectra (Bohr's harmonic oscillator picture). For jellium, the analogous correction is computed from the same digamma function but with v_F replacing some orbital velocity scales. The numerical correction sizes in §4.5 are atomic-target estimates and may be off by a factor of 2 for jellium. The 50 eV simulation is what will pin down the true magnitude.

7. **Mass and charge convention.** S is energy lost *per unit length* of projectile travel, in eV/Bohr. Positive sign means energy lost to the bath. The Bethe formula has Z₁² = +1 for an electron. The Bohr formula has the same Z₁², but the Barkas-Andersen Z₁³ correction (present in nonlinear theory) flips sign between electrons and protons. The current simulations use electron projectiles only, so Z₁³ effects are not directly testable here.

8. **The classical projectile uses a custom UPF pseudopotential** (`electron-ONCV-1.2.upf`) which is a pseudo-hydrogen tuned to behave like an electron. It is *not* an exact electron — it has a non-zero core size and is distinguishable from bath electrons (no antisymmetrization with bath orbitals). This is the fundamental asymmetry between the classical and WP runs and is the source of any non-zero Δ(E). The agent should report the pseudopotential parameters (cutoff radius, projector channels) in the final analysis so readers can assess the implicit smoothing of the projectile-electron Coulomb interaction.

9. **WP-1500 eV is impossible at the current dx = 0.40 Bohr** because the packet's momentum spectrum extends to k_max = k₀ + 3σ_k = 11.10 Bohr⁻¹, 41 % above the grid's Nyquist limit. To get WP-1500 eV would require dx ≤ 0.275 Bohr, which the available GPU VRAM does not support at N = 162 + WP state count. The classical-1500 eV anchor therefore stays unpaired in the Δ(E) plot, and the agent should make this asymmetry explicit in the comparison.

---

## 11. Quick-reference table for the agent

| Operation | Where | How |
|---|---|---|
| Read existing run output | `ResearchProject/systems/jellium/run_*/results/analysis/observables/observables.csv` | `pandas.read_csv` |
| Compute windowed S(v) for a classical run | new postprocess `inqview/postprocess/stopping_classical.py` | (1) `delta_E_bath = (energy_total − energy_total[0]) * 27.2114` in eV. (2) Set `v_initial = sqrt(2 * E_eV / 27.2114)`. (3) `t_start = 3 / v_initial`, `t_end = 28 / v_initial`. (4) Linear regression `delta_E_bath` vs `time_au` in `[t_start, t_end]`. (5) `S = slope / v_initial` in eV/Bohr. (6) Report ±σ from regression covariance. |
| Compute WP `⟨z⟩(t)`, `⟨T⟩(t)`, `σ_z(t)`, `packet_norm(t)` | new postprocess `inqview/postprocess/wavepacket_observables.py` | Loads the projectile-orbital snapshot files from `results/snapshots/`. For each save: `z̄ = ∫ |φ|² z dV`, `T = (1/2) ∫ |∇φ|² dV` (negative-Laplacian form, avoids the second derivative), `σ_z = √(z̄² − ⟨z⟩²)`, `‖φ‖² = ∫ |φ|² dV`. Save as `wp_observables.csv` at the same cadence as the bath CSV (every 2 steps). |
| Compute WP stopping power | new postprocess (extends the above) | `S_WP,bath(t) = (delta_E_bath(t) − delta_E_bath(0)) / (⟨z⟩(t) − ⟨z⟩(0))` in eV/Bohr, plotted parametrically. Take the average slope over the clean window `Δz ∈ [3, 28]` Bohr. Cross-check against `S_WP,kin = −d⟨T⟩/d⟨z⟩` parametrically. |
| Bethe-Lindhard predicted S | (formula) | `S = 4π × 0.001295 × [ln(2 v² / 0.1276) − 0.5] / v²` in Ha/Bohr, then × 27.2114 for eV/Bohr. |
| Bloch correction at given κ | (formula) | `import scipy.special as sp; correction = -sp.digamma(1 + 1j*kappa/2).real + sp.digamma(1).real`. Subtract from the Bethe stopping number L. |
| Configure new classical run at E eV | clone `electron_proj_E600_L50_cubic.hpp` → `electron_proj_E{E}_L50_cubic.hpp` | Edit `WP_EKIN_EV = E`. `N_STEPS` auto-derives from formula. Instantiate `Electron_Proj_E{E}_L50_cubic_Classical_dx0p40` struct in the run. |
| Configure new WP run at E eV | same config file as classical | Instantiate `Electron_Proj_E{E}_L50_cubic_WP_dx0p40` struct in the run. σ inherited as 5.0 Bohr. |
| Reference dt, dx, σ, launch z₀ | from `shared/configs/base_n162_L50_E1p5.hpp` and `electron_proj_E*_L50_cubic.hpp` | dt = 0.020, dx = 0.40, σ = 5.0, z₀ = −10. |

---

## 12. Decision tree for the agent

If the agent encounters issues, follow this flow.

1. **Existing 100 eV / 600 eV classical CSVs not in the format of the 1500 eV CSV** → if columns differ, write a wrapper that maps to the canonical column set `(step, time_au, energy_total, energy_kinetic, energy_hartree, energy_xc, current_*, dipole_*, density_l2)`. If `energy_total` is missing (or labeled differently), reconstruct it as `energy_kinetic + energy_hartree + energy_xc + ion_potential_term` from whatever fields exist. Flag any reconstruction.

2. **Postprocess for windowed S(v) does not exist yet** → write `inqview/postprocess/stopping_classical.py` following §6.1 and the quick-reference table. Test on the 1500 eV CSV; expected S = 0.0249 ± 0.0002 eV/Bohr (the value used as ground truth in this plan). If your code gives a different value, debug before applying to 100 and 600 eV.

3. **Wavepacket-observables module does not exist** → write `inqview/postprocess/wavepacket_observables.py` per §6.2 step 6 and the quick-reference table. Test on the WP 600 eV run first (highest k₀σ, packet closest to classical limit) — `S_WP,bath` should agree with `S_classical(600 eV)` to within a few percent. If it doesn't, there's a normalization or sign error.

4. **Classical 50 eV sim overruns 5.5 h** → kill and reduce traversal: drop `t_end` from `35 / v_initial` to `25 / v_initial`, giving `N_STEPS = 651` and ~3.3 h cost. This loses 10 Bohr of clean window but still allows a windowed slope over `Δz ∈ [3, 18]` Bohr (15 Bohr of clean signal). Re-launch.

5. **WP 50 eV shows packet norm drift > 0.5 %** → numerical issue, not physics. Tighten the ETRS or RK4 time-step control. If `packet_norm` drops below 0.99, the simulation has lost projectile probability to the absorbing boundary (or grid noise) and the S extraction is unreliable.

6. **WP 50 eV shows σ_z growing past 8 Bohr** → real packet dispersion. Record it; the WP stopping power is still extractable via the bath-energy method (`S_WP,bath`), but the kinetic-energy method (`S_WP,kin`) will give a smaller value because some kinetic energy is going into packet spread rather than into the bath. Report both numbers and explain the difference in the writeup.

7. **The matched-pair Δ(E) is large at 600 eV (k₀σ = 33.2)** — i.e. Δ(600) > 5 % — → systematic issue, not physics. The classical-packet limit should give Δ → 0 there. Likely culprits: the classical projectile's UPF cutoff is too soft (smoothing the Coulomb singularity differently than the WP feels it), or the WP `packet_norm` is drifting. Diagnose with the `σ_z` and `packet_norm` traces.

8. **No bath-energy CSV for one of the existing classical runs** → re-run the analysis pipeline (the user has confirmed `observables.csv` is generated for every run; if missing, the run output may need to be regenerated, but this is unusual). If the pipeline never ran post-simulation, run it now: typically `python -m inqview.postprocess.observables <run_dir>`.

---

End of plan. The agent should treat this as a starting point and update the operational details based on what the codebase reveals.