# Plan: N=162 localised-jellium mass-pair WP runs (σ=1, E=100 eV, CAP η=−1.0)

Status: IN PROGRESS (authoring run machinery). Branch: `overnight-gaussian-classical`.
Owner session started 2026-07-19.

## Goal (user spec, 2026-07-19)

Two matched **quantum wavepacket** projectile runs through a **genuine
162-electron localised jellium slab**, one on each GPU. Identical in every
parameter **except the projectile mass**:

- **Run 1 (GPU 0):** projectile mass **1** (INV_MASS = 1.0).
- **Run 2 (GPU 1):** projectile mass **2** (INV_MASS = 0.5). Bath electrons stay
  mass 1 in BOTH runs (per-orbital `electrons.inverse_mass()[0][wp_idx]` fork —
  the inq-study engine). "Only the projectile mass changes."

User-locked parameters:
- 162 electrons, **genuine** count (NOT density-matched — user decision
  2026-07-19; the existing localised slabs are N=52/82/234). 162 is a magic number.
- Localised jellium slab **z ∈ [−12.5, +12.5]** (25 Bohr thick), smoothening
  (erfc edge_width) = **1.0 Bohr**.
- Two CAPs, **10 Bohr each**, at the two cell boundaries; **η = −1.0**.
- σ_WP = **1.0 Bohr**. E = **100 eV**. Sim time = **100 a.u.**
- **dx = 0.40** Bohr — user states this is the FINEST allowed (memory ceiling).
- "Comfortable Lz so the WP is between the CAP and the jellium, ~4σ from each."
- Full energy decomposition recorded EACH timestep + the extensive observable suite.

## Geometry (Lz sizing, the 4σ rule)

σ_WP = 1 → 4σ = 4 Bohr. Launch the WP on the −z side:

```
 edge   CAP(10)     WP        slab face .. slab .. slab face      CAP(10)   edge
-31.2 |--------| -16.5  | -12.5 =========== +12.5 |          |--------| +31.2
      -31.2→-21.2  4.7σ  4σ                                    +21.2→+31.2
```

- WP launch z = **−16.5** (4.0 Bohr = 4σ from the −z slab face).
- CAP inner face −21.2 → WP standoff 4.7 Bohr = 4.7σ (≥ 4σ ✓).
- **Lz = 62.4 Bohr** (156 pts at dx=0.40; ≥ the 61-Bohr tight-4σ minimum, small margin).
- Tight-4σ minimum would be Lz=61; 62.4 chosen for a slightly more comfortable
  standoff and a cleaner grid count.

## Cell / density / grid (locked)

| Quantity | Value | Notes |
|---|---|---|
| Lx = Ly | 70.4 Bohr (176 pts) | sized so N=162 at r_s≈5.68 in the 25-Bohr slab |
| Lz | 62.4 Bohr (156 pts) | 4σ geometry above |
| dx | 0.40 Bohr | user ceiling (finest allowed) |
| n0 | 162 / (70.4²·25) = 1.3054e-3 a₀⁻³ | neutral to 162 e exactly |
| r_s | 5.68 | matches the canonical localised-jellium density (5.665) to 0.3% |
| N_electrons | 162 | genuine |
| edge_width | 1.0 Bohr | erfc-softened slab faces (GS-study H1) |
| extra_states | 18 | larger DOS than N=82 (50×50); pilot confirms sufficiency |
| grid | 176×176×156 ≈ 4.83M pts | ~101 states incl. WP |

**Memory (planning estimate):** 101 states × 4.83M × 16 B × 2.5 ≈ **19.5 GB**
(the 2.5× work-array factor is conservative; actual peak likely ~14–16 GB). This is
the reason dx cannot go below 0.40 on a 24 GB A30. The pilot smoke confirms the fit;
if it OOMs, trim extra_states (→14) or Lz before production.

## Projectile kinematics

E = 100 eV = 3.6749 Ha. k0 = √(2·m·E_Ha):
- **Run 1 (m=1):** k0 = 2.711 Bohr⁻¹, v = 2.711 a.u.
- **Run 2 (m=2):** k0 = 3.834 Bohr⁻¹, v = 1.917 a.u. (INV_MASS=0.5 → v=k0/m).

Traversal launch(−16.5) → far CAP inner face(+21.2) ≈ 37.7 Bohr:
- m=1: ~13.9 a.u. to reach far CAP; m=2: ~19.7 a.u. So the WP is CAP-absorbed by
  ~15–20 a.u.; the remaining ~80 a.u. of the 100-a.u. window is bath relaxation.
- No chirped focus (`.focus_z` needs dx=0.333; the config itself notes dx=0.40
  fails it). Plain drifting Gaussian → σ_ψ disperses; acceptable since the packet
  is absorbed quickly. This is the honest σ=1 free-Gaussian behaviour.

## Grid-cutoff (aliasing) gate — PASSED (mandatory, pre-launch)

`cutoff_guard.py --kind wp --sigma-wp 1.0 --spacing 0.40`:
- Run 1 (E=100 eV, p0=2.71): PASS, aliased tail 0.00%, k_Nyq=7.85 ≥ p0+3σ_p=4.83.
- Run 2 (m=2 → effective E=200 eV so p0=3.83): PASS, tail 0.00%, k_Nyq=7.85 ≥ 5.96.
- σ_p = 1/(√2·σ_WP) = 0.707 (NOT 1/(2σ)). E_cut(dx=0.40)=839 eV ≫ 100 eV.

## Time integration

- dt = **0.04 a.u.**, N_STEPS = **2500** (= 100 a.u.). Precedent: effmass_sigma1
  validated ETRS at dt=0.04 with k0=5.693 (higher than here). Pilot checks energy
  drift; if |ΔE_total| > 1 mHa, drop to dt=0.02 (→5000 steps).
- Propagator ETRS, γ-only, LDA, inq-study (mass-fork) engine.
- **Cost projection:** grid ~2.5× the 125³/101-state reference (16 s/step) ⇒
  ~40 s/step ⇒ ~28 h/run at 2500 steps. Two GPUs in parallel ⇒ ~28 h wall clock.
  Per `checkpoint-dont-block`: launch at full scope, checkpoint every **500 steps**
  (user directive 2026-07-19; kill loses ≤500 steps) plus a final checkpoint, WARN
  the user with the measured smoke projection; the user owns the kill. Interior
  checkpoints at steps 500/1000/1500/2000; `EM_RESUME=1` extends from the last.

## Observables

- **Per-step full energy decomposition** → `energy_decomp.csv`: total, kinetic,
  hartree, xc, external, nonlocal, ion, ion_kinetic, exact_exchange, nvxc,
  eigenvalues (from `data.energy().*`; the 8 non-diagnostic terms sum to total).
- **Extensive suite at WRITE_EVERY=8** (~313 frames), cloned from `fullsuite_wp`:
  density VTIs total/system/gs_system/wp + complex wavefunction_wp (WF_EVERY=40),
  density_delta (+coarse), observables.csv (energies+current+dipole+L2),
  state_energies.csv + occupations_vs_time.csv (5×WRITE_EVERY), eigenvalues dump,
  momentum_distribution, wp_momentum_stats, wp_real_space_stats, orbital_overlap
  (WP-only every 10 steps) + overlap_full (t=0 and t=final excitations).
- density_total + density_wp at EQUAL cadence (canonical bath n_sys = total − wp;
  min-obs-set rule). Manifest via `minimum_observable_set.hpp` (RunType jellium_wp).

## File layout (ADR 0007 / muon_mass_fork convention)

- Config: `shared/configs/slab_n162_L70x70x62.hpp` (NEW).
- GS builder: `scripts/mass_pair_n162/gs/run.cpp` → converge →
  `shared_gs/slab_n162_L70x70x62_dx0p40/` (PREREQUISITE — no N=162 GS exists).
- WP run (both masses, env-selected): `scripts/mass_pair_n162/wp/run.cpp`.
- Orchestrator: `scripts/mass_pair_n162/orchestrate.py` (Python; GS → pilot →
  two production runs on GPUs 0/1 → analyse → email).
- Per-run analyse: `scripts/mass_pair_n162/analyse.py` (template).
- Run outputs: `mass_pair_n162/{m1,m2}/results/…` (logs gitignored).

## Steps

1. [ ] Config header (SlabN162_L70x70x62).
2. [ ] GS run.cpp; build; converge GS; validate (SCF conv, N=162, energy finite).
3. [ ] WP run.cpp (fullsuite obs + mass fork + checkpoint/resume + per-step energy
   decomp + plain launch + CAP 10 Bohr η−1.0); build once.
4. [ ] Pilot smoke (≤20 steps, m=1): confirm memory fit, WP norm∈[0.95,1.05],
   energy drift, timing → cost projection.
5. [ ] Launch production: m=1 on GPU0, m=2 on GPU1 (INV_MASS 1.0 / 0.5), 2500 steps.
6. [ ] analyse.py both; run notebook + density GIFs; catalogue; email.

## Validation gates

- Pre-launch: cutoff gate PASSED (above). Memory pilot. GS SCF converged (tol 1e-4).
- Post-pilot: norm, energy drift (<1 mHa target), no NaN.
- Post-run: run_completed=true, energy-decomposition closure, WP absorbed.

## Open risks

- Memory ~19.5 GB est. on 24 GB A30 (tight; pilot-gated).
- 162 may not be a clean closed shell for a SLAB (magic numbers are 3D-isotropic);
  T=100 K smearing + extra_states handles partial occupation — not a blocker.
- Cost ~28 h/run; managed by checkpoint-don't-block, not gated.
