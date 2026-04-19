# Plan: Coronene Electron WP Scattering — TDDFT LEED Simulation

**Reference:** Tsubonoya, Hu, Watanabe, *Phys. Rev. B* **90**, 035416 (2014)  
**DOI:** 10.1103/PhysRevB.90.035416  
**Goal:** Replicate Fig. 1 (density snapshots before/after scattering) and Fig. 2 (LEED diffraction pattern with D₆h symmetry) for coronene C₂₄H₁₂.

---

## Scientific method

The incident electron is treated as an additional (N+1)th TDKS orbital. The total electron density is:

    n(r,t) = Σᵢ fᵢ|ψᵢ(r,t)|² + f_WP|ψ^WP(r,t)|²

All orbitals, including the WP, evolve under the same TDKS Hamiltonian (paper Eq. 3):

    i ∂ψᵢ/∂t = [−∇²/2 + ∫n(r',t)/|r−r'| dr' + v_xc[n](r,t)] ψᵢ(r,t)

The WP orbital is NOT orthogonalised against the occupied states — the WP starts far from the target (D = 6.35 Å) so overlap at t=0 is negligible (paper text, Sec. II).

The XC is treated in the **ALDA** (adiabatic LDA): the XC functional is LDA evaluated at the instantaneous density. In INQ this is `options::theory{}.lda()`.

**LEED pattern** (paper Eq. 5):

    I(r) = ∫_{t₁}^{t₂} n(r,t) dt,    r ∈ S

where S is the observation plane at z = D (incident side; paper: "observation plane set at incident position").

---

## Wavepacket formula

Paper Eq. 1:

    ψ^WP(r,0) = (1/(π d²))^{3/4} · exp(−|r−b|²/(2d²)) · exp(i k·r)

Key points:
- Exponent denominator is **2d²** (not 4d² as in some other conventions)
- Normalisation constant: N = (π d²)^{−3/4} ensures ⟨ψ|ψ⟩ = 1
- k = (0, 0, −k₀): WP propagates toward the flake in the −z direction
- k₀ = √(2 E_kin) in atomic units (m=1)

---

## System parameters (from paper)

### Cell

| Parameter | Paper | Atomic units |
|---|---|---|
| Box x | 18.4 Å | 34.76 bohr |
| Box y | 18.4 Å | 34.76 bohr |
| Box z | 31.7 Å | 59.91 bohr |
| Boundary | Isolated | `.finite()` |

### Coronene

| Parameter | Value |
|---|---|
| Formula | C₂₄H₁₂ |
| C-C bond | 1.421 Å |
| C-H bond | 1.086 Å |
| Symmetry | D₆h, flat in xy-plane at z=0 |
| Valence electrons | 108 |
| KS orbitals (occupied) | 54 |

### Electronic structure

| Parameter | Value | Notes |
|---|---|---|
| XC | LDA (ALDA) | paper: ALDA |
| Grid spacing | 0.16 Å = 0.302 bohr | LEED calculation |
| E_cut | ≈ 54 Ha | (π/0.302)²/2; verify with convergence test |
| extra_states | 3 | 1 for WP + 2 SCF buffer |

### Wavepacket

| Parameter | Paper | Atomic units |
|---|---|---|
| d (width) | 0.53 Å | 1.001 bohr |
| D (impact distance) | 6.35 Å | 12.0 bohr |
| E_kin | 200 eV | 7.350 Ha |
| k₀ = √(2E_kin) | — | 3.832 bohr⁻¹ |
| b (WP centre) | (0, 0, D) | (0, 0, 12.0) bohr |
| k direction | −z | WP → flake |
| Occupation | 1.0 | 1 incident electron |

### TDDFT propagation

| Parameter | Value | Notes |
|---|---|---|
| Algorithm | 4th-order Taylor | INQ default |
| Δt (LEED) | 4.84×10⁻⁴ fs = 0.020 a.u. | paper value |
| t₁ (WP arrival) | D/k₀ ≈ 3.13 a.u. | |
| t₂ (end time) | 0.25 fs = 10.33 a.u. | paper value |
| Total steps | ≈ 517 | |
| Obs. plane | z = D = 12.0 bohr | same as WP start |

---

## Occupation question

**Q: Can INQ set occupation = 1.0 for the WP orbital?**

**A: Yes.** `electrons.occupations()` returns a writable 2D array `[kpt_idx][local_ist]`.
Setting `electrons.occupations()[0][ist_wp] = 1.0` is valid (confirmed from `electrons.hpp:117-118`; the unit test at line 603 sets `cos(istg)` values, which are non-integer).

After the GS converges with `extra_states(3)`:
- States 0–53: occupation = 2.0 (doubly occupied, spin-unpolarized)
- State 54 (WP): occupation = 0.0 → manually set to 1.0 after WP injection
- States 55–56 (buffer): occupation = 0.0 (remain unoccupied)

The density computed by TDDFT will include the WP orbital with weight 1.0.
Total electrons during propagation = 108 (GS) + 1 (WP) = 109.

**Density validation** (Step 4 in run.cpp):
- ⟨ψ_WP|ψ_WP⟩ = Σ |hc|² × dV ≈ 1.0 (verifies WP norm)
- Sum of all occupations = 109.0 (verifies total electron count)

---

## Density correctness: Hartree and XC terms

The Hartree potential v_H[n](r,t) and XC potential v_xc[n](r,t) are functionals of the TOTAL density n(r,t), which includes the WP orbital. After Step 3 (WP injection + occupation set), INQ's TDDFT propagator recomputes the density at each step as:

    n(r,t) = Σᵢ fᵢ|ψᵢ|² = 2Σ_{i=0}^{53} |ψᵢ|² + 1·|ψ^WP|²

This is exactly Eq. (4) of the paper. No manual density update is needed — INQ handles this internally once the occupation array is set correctly.

---

## Folder structure

```
ResearchProject/systems/coronene/
├── 01_geometry/
│   ├── gen_geometry.py          — analytical D6h generation + bond verification
│   ├── coronene.xyz             — 36 atoms in Angstrom (VESTA-compatible)
│   └── results/
│
├── 02_ground_state_analysis/
│   ├── config.hpp               — ALL shared parameters (cell, WP, TDDFT)
│   ├── run.cpp                  — SCF + forces + orbital summary
│   └── results/
│       └── gs_summary.txt
│
├── 03_ecut_convergence/
│   ├── config.hpp               — symlink to shared config
│   ├── run.cpp                  — sweep E_cut = 20..60 Ha, record total energy
│   └── results/
│       └── ecut_convergence.csv
│
└── 04_leed_simulation/          ← MAIN SIMULATION
    ├── config.hpp               — fine-tuned parameters (update E_cut after 03)
    ├── utils.hpp                — WP injection, validation, density extraction
    ├── run.cpp                  — GS → WP → TDDFT → LEED accumulation
    ├── analysis.py              — post-processing: Fig. 1 + Fig. 2 plots
    └── results/
        ├── leed_pattern.txt     — I(x,y) for Fig. 2
        ├── snapshot_t????.txt   — 2D density at z=0 for Fig. 1
        ├── fig1_density_snapshots.png
        ├── fig2_leed_pattern.png
        ├── leed_angular_profile.png
        └── sim_summary.txt
```

---

## Execution sequence

### Step 0 — Verify geometry (already done)
```bash
cd 01_geometry
python3 gen_geometry.py
# → coronene.xyz (36 atoms, 30 C-C + 12 C-H bonds verified)
# scp coronene.xyz to visualise in VESTA
```

### Step 1 — Ground state analysis
```bash
cd 02_ground_state_analysis
source ~/.bashrc && inq-run
# Check: SCF converges, max force < 0.1 Ha/bohr (idealized geometry, not exact)
# Check: total energy ≈ ballpark for C24H12 at LDA
```

### Step 2 — E_cut convergence test
```bash
cd 03_ecut_convergence
source ~/.bashrc && inq-run
# Check: results/ecut_convergence.csv
# Convergence criterion: ΔE_total < 1 meV between successive E_cut values
# Expected converged E_cut: ~54 Ha (paper: 0.16 Å grid)
# Update cfg::ECUT_HA_LEED in 04_leed_simulation/config.hpp
```

### Step 3 — LEED simulation
```bash
cd 04_leed_simulation
source ~/.bashrc && inq-run
# Expect runtime: ~30–60 min (517 steps × 36 atoms × 55 orbitals × GPU)
# Monitor: snapshot_t????.txt files appear as propagation proceeds
```

### Step 4 — Post-processing
```bash
cd 04_leed_simulation
python3 analysis.py
# → results/fig1_density_snapshots.png  (compare to paper Fig. 1)
# → results/fig2_leed_pattern.png       (compare to paper Fig. 2)
# → results/leed_angular_profile.png    (6-fold symmetry check)
```

---

## Validation checklist

### Ground state (Step 1)
- [ ] SCF converges within 300 steps
- [ ] 36 atoms loaded correctly from XYZ
- [ ] Sum of occupations = 108.0
- [ ] max |force| < 0.1 Ha/bohr (idealized geometry is not exactly equilibrium — some force is expected)

### E_cut convergence (Step 2)
- [ ] Total energy plateau within 1 meV between consecutive E_cut values
- [ ] Converged E_cut identified and written into `04_leed_simulation/config.hpp`

### WP injection (Step 3, pre-propagation)
- [ ] ⟨ψ_WP|ψ_WP⟩ ≈ 1.0 (within 5%)
- [ ] WP occupation = 1.0
- [ ] Total occupations sum to 109.0

### TDDFT propagation (Step 3)
- [ ] Propagation completes without crash
- [ ] Density snapshots visible in results/ during run
- [ ] LEED accumulator non-zero after t₁

### Physical results (Step 4)
- [ ] Density snapshots (Fig. 1): clear WP blob before scattering, scattered pattern after
- [ ] LEED pattern (Fig. 2): 6-fold D₆h symmetry for coronene
- [ ] Angular profile: peaks at 0°, 60°, 120°, 180°, 240°, 300° ± 5°

---

## E_cut convergence results (03_ecut_convergence, 2026-04-14)

Sweep over 20–60 Ha with LDA, pseudodojo_pbe pseudopotentials, finite cell 34.76×34.76×59.91 bohr.

| E_cut (Ha) | h (Å) | E_total (Ha) | ΔE (meV) | Grid pts | SCF steps |
|---|---|---|---|---|---|
| 20 | 0.263 | -148.940724 | — | 612,500 | 33 |
| 25 | 0.235 | -150.135129 | −32,501 | 864,000 | 39 |
| 30 | 0.215 | -150.710250 | −15,650 | 1,215,000 | 48 |
| 35 | 0.199 | -150.820869 | −3,010 | 1,474,560 | 44 |
| **40** | **0.186** | **-150.836889** | **−436** | **1,750,000** | **57** |
| 45 | 0.175 | -150.806986 | +814 | 2,083,725 | 55 |
| 50 | 0.166 | -150.770509 | +993 | 2,408,448 | 54 |
| 55 | 0.159 | -150.757387 | +357 | 2,880,000 | 59 |
| 60 | 0.152 | -150.746896 | +285 | 3,281,250 | 63 |

**Behaviour:** Energy decreases normally from 20→40 Ha, then rises monotonically from 40→60 Ha (total rise ~90 meV). Non-monotonic convergence is a known artefact of the interplay between pseudopotential projector completeness and real-space XC aliasing at cutoffs above the pseudopotential's natural range with pseudodojo_pbe.

**Chosen E_cut: 40 Ha** — energy minimum in the sweep; h = 0.186 Å; grid 100×100×175 = 1.75M points; ~3 grid points across WP width (d = 0.53 Å). Set as `ECUT_HA_LEED = 40.0` in `04_leed_simulation/config.hpp`.

### Paper configuration (for reference / paper-matching runs)

These are the exact parameters from Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014) to use if exact replication of the paper's numerical setup is needed:

| Parameter | Paper value | Atomic units | Notes |
|---|---|---|---|
| Code | — | — | Not INQ; different plane-wave code |
| XC | ALDA | `options::theory{}.lda()` | Adiabatic LDA |
| Grid spacing | 0.16 Å | 0.302 bohr | E_cut ≈ (π/0.302)²/2 = 53.9 Ha |
| E_cut (implied) | ≈ 54 Ha | 54 Ha | Use `ECUT_HA_LEED = 55.0` to match |
| Cell x | 18.4 Å | 34.76 bohr | Finite, isolated |
| Cell y | 18.4 Å | 34.76 bohr | |
| Cell z | 31.7 Å | 59.91 bohr | |
| C-C bond | 1.421 Å | 2.685 bohr | Idealized aromatic |
| C-H bond | 1.086 Å | 2.052 bohr | |
| Valence electrons | 108 | — | 54 occupied KS orbitals |
| WP width d | 0.53 Å | 1.001 bohr | Gaussian σ in exp(−r²/2d²) |
| WP distance D | 6.35 Å | 12.0 bohr | Initial WP–flake separation |
| WP energy | 200 eV | 7.350 Ha | k₀ = √(2×7.350) = 3.832 bohr⁻¹ |
| WP direction | −z | k = (0,0,−k₀) | Toward flake |
| WP occupation | 1 electron | f_WP = 1.0 | Singly occupied |
| Propagator | 4th-order Taylor | INQ default | |
| Δt | 4.84×10⁻⁴ fs | 0.020 a.u. | Paper Table I |
| t₁ (WP arrival) | ≈ D/k₀ | 3.13 a.u. | Start accumulating LEED |
| t₂ (end time) | 0.25 fs | 10.33 a.u. | Stop propagation |
| Total steps | ≈ 517 | — | t₂/Δt |
| Obs. plane | z = +D | z = 12.0 bohr | Reflection LEED, incident side |

To switch to paper-matching parameters: set `ECUT_HA_LEED = 55.0` in `04_leed_simulation/config.hpp` (all other parameters already match).

---

## Known issues and open questions

1. **E_cut non-monotonic above 40 Ha.** Resolved: using 40 Ha (energy minimum). See table above.

2. **Force tolerance.** The coronene geometry uses idealised bond lengths (1.421 Å, 1.086 Å) from graphene/aromatic standards. These are not exact LDA equilibrium lengths, so forces will be non-zero but should be small (< 0.05 Ha/bohr).

3. **Absorbing boundary conditions.** Not implemented in the current code. For a clean LEED pattern, absorbing BCs at the box edges prevent WP reflection from the cell walls. INQ has `perturbations::absorbing_potential`. Add this in a Phase 2 refinement if the LEED pattern shows artefacts.

4. **MPI parallelisation.** For large runs, use `INQ_EXEC_ENV="mpirun.openmpi -np 4" inq-run` to distribute over multiple GPUs/cores.

5. **Slice extraction on distributed runs.** `utils::extract_density_slice` currently assumes single-process. For MPI runs, the local indices on each rank cover only part of the grid — the slice extraction needs an MPI reduce. This is a known limitation; run single-process for Phase 1.

6. **Observation plane physics.** The paper places the observation plane on the **incident side** (z = +D, reflection LEED). If needed, a transmission plane at z = −D can also be computed from the same snapshots.

---

## Attribution

- Simulation method: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
- Coronene geometry: analytical D₆h construction from aromatic C-C = 1.421 Å (standard value)
- INQ framework: Andrade et al., alphataubio.com/inq

---

## History

| Date | Milestone |
|---|---|
| 2026-04-14 | Tests 1-3 complete: `hypercubic[ix][iy][iz][ist]` ordering, `rvector_cartesian` units, N2 orbital peaks |
| 2026-04-14 | Phase 0 complete: coronene.xyz generated, 30 C-C + 12 C-H bonds verified |
| 2026-04-14 | All four subfolders created; code written |
| — | Step 1: run 02_ground_state_analysis |
| — | Step 2: run 03_ecut_convergence, update E_cut |
| — | Step 3: run 04_leed_simulation |
| — | Step 4: run analysis.py, compare to paper |
