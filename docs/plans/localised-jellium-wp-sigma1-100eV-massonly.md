# Plan: σ=1 mass-only WP through r_s≈5.67 slab (100 eV) + vacuum CAP calibration

Status: **in progress** (design locked with user 2026-07-09).
Owner sweep: `ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/sigma1_massonly/`

## Goal & decisions (user, 2026-07-09)

Prior increased-mass σ=2 runs (`effmass_sigma2`, `effmass_12h`) rejected: **CAP
reflection** + **packet too wide**. New run: a **concentrated σ=1** projectile,
**mass-only** (accept dispersion), inside the **1–2 h** budget. Second GPU runs a
**vacuum calibration** to measure CAP reflection + true spreading before trusting
the slab result.

**Sanity-check verdict (grounded):** a σ=1 packet **cannot** be kept rigid across
the 25-Bohr slab by mass tuning — fractional spread `√(1+(d/k₀σ²)²)` only improves
as 1/√m, and m is capped by both aliasing (≤5.5 at dx=0.4) and cost (heavier =
slower = more steps). At the cost-feasible mass the packet is ~3.4× wider at the
slab centre than the classical Gaussian. This run is therefore the **honest σ=1
dispersion-dominated end-member**, not a width-matched comparison. (User chose this
knowingly over chirp-focus / σ≈1.5.)

## Locked config — slab run

| Item | Value | Note |
|---|---|---|
| σ_WP | 1.0 Bohr (density std 0.71) | |
| E | 100 eV | drift KE = ½·m·v² |
| m_eff | 3.45 mₑ (INV_MASS = 0.289855) | cost-capped for ≤2 h; clean aliasing |
| k₀ | 5.0356 Bohr⁻¹ = √(2·E·m) | **fork-corrected** momentum |
| v | 1.460 a.u. | = k₀/m |
| Grid | dx = 0.40, 50×50×64 (125×125×160) | "same grid"; 4σ standoffs shrink L_z 90→64 |
| Slab | 25 Bohr thick, faces ±12.5, N=82, r_s≈5.67 | new GS `slab_n82_L50x50x64_dx0p40` |
| Launch z₀ | −16.5 | 4σ_WP to near face (−12.5) AND to CAP inner edge (−20.5) |
| CAP | η=−0.6, sin² band, peak ±26.25, width 11.5, band [±20.5,±32] | **gentler than σ=2's −1.0**; graded already |
| dt | 0.05 (dt·E_cut=1.5) | stable |
| N_STEPS | 1192 = 3× traversal (launch→far face 29 Bohr) | |
| write_every | 4 (≈298 frames) | 300-frame cadence rule |
| ckpt_every | 300 (3 interior checkpoints + final) | pause/continue |

Aliasing (cutoff_guard, fork-corrected E_eff = E·m = 345 eV): **PASS**, tail 0.00%,
k_Nyq 7.85 ≥ p₀+3σ_p = 7.16.

## Vacuum calibration (second GPU)

Same box/grid/CAP/mass but **no slab** (empty box, 1 spectator + WP, `non_interacting`).
Launch −16.5, k₀/m/CAP identical. Measures:
1. **CAP reflection** — WP norm drainage curve + does a reflected −z lobe appear
   (density VTI). This is the direct test of the user's reflection complaint.
2. **True spreading** vs the free-Gaussian oracle σ_ρ(t)=√(σ_ρ0²+(t/2mσ_ρ0)²).

Gates the slab CAP: if reflection is non-trivial, retune η/width before trusting
the slab stopping. Cheap (2 orbitals vs 82) — minutes.

## Files

- `shared/configs/slab_n82_L50x50x64.hpp` (clone of L90, LZ=64, dx=0.40)
- `scripts/muon_mass_fork/sigma1_massonly/gs/run.cpp` (→ `shared_gs/slab_n82_L50x50x64_dx0p40`)
- `scripts/muon_mass_fork/sigma1_massonly/wp/run.cpp` (slab WP, checkpoint/resume, graded CAP)
- `scripts/muon_mass_fork/sigma1_massonly/vacuum/run.cpp` (CAP+spreading calibration)

## Validation

- cutoff_guard: PASS (above).
- Vacuum calibration BEFORE trusting slab stopping (gates CAP).
- S extraction: initial-drag window (light-projectile rule) — packet decelerates.
- phantom-absorbed-orbital caveat applies once CAP eats the WP.

## Not done / open

- Matched classical rigid-Gaussian twin (σ_pot=0.71, m=3.45) — build AFTER the WP
  config is confirmed by the vacuum calibration.
- Run-notebooks for both once complete.
