# Wavepacket Spreading Investigation

## Physical spreading formula

For a free Gaussian wavepacket (m=1, ħ=1 in atomic units):

    σ(t) = σ₀ √(1 + t²/σ₀⁴)

where t is time in a.u. and σ₀ is the initial Gaussian width in bohr.

This follows from the time-dependent Schrödinger equation for a free particle.
Source: any quantum mechanics textbook, e.g. Griffiths §2.4.

### Timescale for doubling

σ doubles (σ(t) = 2σ₀) when t = σ₀² √3.

| σ₀ | doubling time (a.u.) | doubling time (fs) |
|----|----------------------|--------------------|
| 0.501 bohr (0.265 Å) | 0.435 a.u. | 0.011 fs |
| 1.002 bohr (0.53 Å)  | 1.74 a.u.  | 0.042 fs |
| 3.780 bohr (2.0 Å)   | 24.7 a.u.  | 0.598 fs |

**Key point**: A σ=0.53 Å wavepacket doubles its width in only 0.042 fs. LEED
experiments run orders of magnitude longer than this in classical terms, but
the WP reaches the target molecule at 200 eV before significant spreading.

## LEED paper analysis (Tsubonoya et al. PRB 90, 035416, 2014)

- WP energy: 200 eV → v = √(2·E/m) = √(2 × 200/27.21) = 3.834 bohr/a.u.
- Impact distance: D = 6.35 Å = 12.0 bohr
- Transit time to coronene: t_impact = D/v = 12.0/3.834 = 3.13 a.u. = 0.076 fs
- WP width at impact: σ(3.13) = σ₀ √(1 + 3.13²/σ₀⁴)
  - For σ₀ = 1.002 bohr: σ(3.13) ≈ 1.002 × √(1 + 9.79) ≈ 3.31 bohr (×3.3 spread)

**This is the discrepancy**: naively σ triples before reaching coronene.
However, Tsubonoya et al. saw negligible spreading effects in their LEED patterns.
Possible explanations:
1. The TDDFT propagation is not free (Hartree + XC potentials act on WP).
2. The pattern is recorded time-integrated (spread WP still produces similar pattern).
3. At 200 eV the WP is already very broad relative to lattice spacing.

## Jellium run_01_base observations (2026-04-20)

### What was seen

- **Frame 0 (t=0):** Gaussian WP density clearly visible — compact, localised blob near z=0 face
  of the periodic cell, on top of the (near-uniform) jellium background. This is the manually
  written frame: `density::total()` (40 e⁻ jellium) + `rho_wp` (1 e⁻ WP orbital).
- **Frames 1+ (t>0):** WP is no longer visually distinguishable. The density looks like pure
  uniform jellium background with no localised feature.

Frame 0 is NOT a ground-state image — it is the t=0 state immediately after WP injection, with
the WP orbital superimposed manually. It is the only frame that shows the WP explicitly, because
subsequent frames write `density::total()` which may or may not include the propagated extra state.

### Hypothesis: density::total() excludes extra states during propagation

If `density::total()` in INQ real-time output only sums the first N_occ bands and excludes
extra (unoccupied) states, then frames 1+ show only the 40 jellium electrons. The WP orbital
is propagated internally but never written to the density output. The WP "vanishes" because
it is invisible in the written field, not because it is physically gone.

Evidence for this:
- Jellium N_elec mean = 40.785 across all 101 frames. The systematic +0.785 offset is a
  finite-grid normalisation artifact (constant, not related to WP). If the WP were included
  in frames 1+, we would expect ~41.785 on average. The observed ~40.785 is consistent with
  the WP being excluded from `density::total()` after t=0.
- The t=0 manual addition was introduced precisely because `density::total()` was known to
  exclude extra states at t=0 (handover note). It is likely consistent throughout propagation.

**Implication for all runs:** The density movies and screen patterns will not show the WP
orbital explicitly. Screens accumulate `density::total()` at each step — again WP excluded.
This would mean the screens only record the jellium/coronene response, not the WP passage.

### Things to investigate / TODO

1. **Confirm `density::total()` behaviour during RT propagation.**
   - Option A: read INQ source (`inq/src/observables/density.hpp`) to check if extra states
     are included in the RT density sum.
   - Option B: run a minimal 1-electron extra-state-only system; if N_elec in frame 1 = 0,
     it is confirmed that extra states are excluded.

2. **Fix the screen accumulator if needed.**
   - If `LeedPatternAccumulator::accumulate()` uses `density::total()`, it misses the WP.
   - Should instead accumulate the WP orbital density: `density::orbital(electrons, state_idx)`.
   - This requires passing `state_index` from the injection report into the RT callback.

3. **Fix density writing during propagation.**
   - For all runs, the RT density frame writes should add `rho_wp` just as at t=0, or
     use a different observable that includes the extra state.

4. **Alternative: finite jellium test** (user suggestion).
   - Create a finite-cell jellium (N_elec small, e.g. 4-8 e⁻) with a localised positive
     background confined to a small sphere/box region. Finite boundary conditions prevent
     periodic wrapping. Inject WP far from the jellium blob and propagate.
   - This would isolate WP propagation from the periodic-cell background subtlety, and give
     a clear visual of WP–electron-gas scattering.

5. **Verify coronene runs won't have the same problem.**
   - Coronene has 108 electrons; the WP is 1 extra state. Same issue as jellium.
   - Density movies will show only coronene density fluctuations, not WP passage.
   - If screens accumulate only `density::total()`, the LEED pattern will be from coronene
     response, not from the WP. This may or may not be the intended physical quantity.

---

## To-do

- Use free-propagation run_01_base z-profile to extract σ(t) numerically.
- Compare numerical σ(t) to analytic formula above: verify TDDFT free-particle limit.
- Find (D, E_kin) regime where σ < 2σ₀ on arrival (projectile-like regime).
- Compare coronene LEED patterns at D=3,6,10,15,20 Å to assess practical impact.

## Free propagation expected outcomes

For run_01_base (σ=0.53 Å, 200 eV):
- WP starts at z ≈ Lz − 5σ ≈ 84.9 bohr, moves in −z.
- Hits z=0 boundary at t ≈ Lz/v = 89.9/3.834 = 23.4 a.u. ≈ 0.57 fs.
- At boundary crossing: σ ≈ 1.002 √(1 + 23.4²/1.002⁴) ≈ 23 bohr (×23 spread).
- After boundary reflection: WP propagates back; σ(t) still measurable from z-profile width.
- Expected: energy constant, N_elec = 1.0 ± 0.001, σ(t) matches analytic in early t.
