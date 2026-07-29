# A3 — semi-empirical far-field forensics

Campaign: `docs/campaigns/localised_jellium_parameter_study_2/`, task A3. Generated 2026-07-11
(autonomous phase). Interpretation: **user's**.

Target: the non-zero constant field beyond ~±15 Bohr in
`docs/reports/09-07-2026-meetng-emilio/figures/s1_3_semiempirical_field_potential.png`,
where Gauss's law (enclosed charge ≈ 0) predicts E → 0.

## (i) The analysis chain, step by step (line-by-line recheck)

Source density (both models): planar-mean n_e(z) of the **Lz=160, 3D-periodic, sharp-edge
(edge_width=0)** GS `runs/extend_r160/gs_lz160_p3` (E_GS = −160.9299 Ha), N = 82, spacing 0.5.

1. `load_vti(..., expect_centered_axis='z')` → physical-order n_e(x,y,z); `.mean(axis=(0,1))` → n_e(z). ✓ correct loader, no fftshift.
2. Background n₊(z) = sharp top-hat, n₀ = 1.312e-3, |z| ≤ 12.5. Grid-commensurate: 50 cells × 0.5 Bohr = exactly 25 Bohr ⇒ areal charge exact (n₀·25 = 82/2500). ✓
3. Symmetrise n_e ← ½(n_e + n_e[::-1]) (removes the raw −8.9e-3 dipole; physical slab symmetric). ✓
4. ρ = n₊ − n_e. Neutrality: ∫ρ dz = −1.9e-15 (machine-exact, raw AND symmetrised). ✓
5. φ(z) = −2π Σ ρ(z′)|z−z′| dz′ (1D Poisson sheet stack, O(N²) sum). ✓ standard.
6. Gauge: φ −= ½[φ(z_min+5) + φ(z_max+5)] (constant shift — cannot create a field). ✓
7. Gaussian convolution `gaussian_filter1d(φ, σ_pot/dz)`, σ_pot = σ_WP/√2 = 0.354 (felt by the Gaussian projectile). ✓ convention-correct.
8. E(z) = −∇φ via `np.gradient`. ✓
9. Second curve "analytical uniform-slab net": electrons idealised as top-hat of half-width a_e = √3·σ_z(n_e) = **15.391 Bohr** (NOT grid-commensurate — see below), same areal charge.

No coding error found in the chain. The plateau is faithfully computed from the inputs.

## (ii) Enclosed charge Q(|z| < Z) — the user's Gauss-law premise, tested

| Z (Bohr) | Q enclosed (e, over the 50×50 cell) |
|---:|---:|
| 12.5 | +3.24 |
| 15 | +0.89 |
| 20 | +0.40 |
| 30 | **+0.39** |
| 50 | **+0.39** |
| 79 | +0.02 |

The premise "enclosed charge ≈ 0 beyond ±15" is **not satisfied by the data**: a +0.39 e
deficit persists for every window 20 ≤ Z ≤ 50 and only closes at the box edge.

**Identity check**: E_plateau = 2π·Q_enc/A = 2π·(0.3915/2500) Ha/Bohr = **0.0268 eV/Bohr**
— exactly the measured plateau (0.0268). The plateau IS the enclosed-charge deficit; the
model's electrostatics are self-consistent.

## (iii) Density spill — where the missing 0.39 e sits

| region | electrons | % of 82 |
|---|---:|---:|
| \|z\| > 15 | 1.12 | 1.4% |
| \|z\| > 25 | 0.392 | 0.48% |
| \|z\| > 50 | 0.392 | 0.48% |
| \|z\| > 70 | 0.322 | 0.39% |

n_e never reaches 0: a near-uniform **vacuum floor of ~8.4e-6 e/Bohr³** extends to the box
edges (×2500 Bohr² × ~45 Bohr of vacuum ≈ 0.39 e). The user's suspicion — "non-zero density
flowing out ... up to the boundaries" — is **confirmed as the mechanism**.
*Inference (labelled)*: a genuine surface evanescent tail would be immeasurably small 50 Bohr
from the face, so the flat floor is most plausibly the SCF/eigensolver numerical floor of
the periodic GS density (0.6% of bulk density), not physics. Proposed check (not run): compare
the floor across the h2 Lz-sweep GSs — a numerical floor stays ~constant, a physical tail decays.

## (iv) Plate thickness (dz) sensitivity

| grid | E(+30) eV/Bohr | net Q on grid |
|---|---:|---:|
| native dz = 0.5 | +0.0268 | 0.0000 |
| dz × 2 = 1.0 | +0.0261 | 0.0000 |
| dz × 4 = 2.0 | +0.2489 | +3.28 e (top-hat quadrature broken) |

At the native spacing the plate discretisation is NOT the cause; only a 4× coarsening breaks it.

## (v) The w parameter (erfc edge softening)

| w (Bohr) | E(+14) | E(+16) | E(+30) plateau |
|---:|---:|---:|---:|
| 0 (sharp) | +0.150 | +0.057 | +0.0268 |
| 0.5 | +0.150 | +0.057 | +0.0268 |
| 1.0 | +0.140 | +0.057 | +0.0268 |
| 2.0 | +0.078 | +0.054 | +0.0268 |

w reshapes the **near-face** field only; the far plateau is w-independent (softening moves
charge locally, preserving every enclosed-charge integral for Z outside the edge region).
**w is not the cause.** (Note: the GS behind this figure was itself computed at edge_width = 0,
so model and run backgrounds already match.)

## Ablations on the other suspects

- **Symmetrisation**: raw density gives an asymmetric plateau (+0.0285 / −0.0250); the
  symmetrised ±0.0268 is their reconciliation. Not the cause, only the ± asymmetry.
- **Gaussian convolution**: plateau identical with convolution off. Not the cause.
- **The "analytical uniform-slab" curve's plateau (−0.0395 eV/Bohr, opposite sign)** has a
  *different*, purely numerical origin: its electron top-hat edge a_e = 15.391 Bohr is not
  grid-commensurate, so the grid quadrature loses −0.578 e ⇒ 2π·0.578/2500 = 0.0395. A
  half-cell edge correction (or analytic integration) removes it.

## Ranked causes (evidence-ordered, for the user's verdict)

1. **Vacuum-floor charge** (~0.39 e pooled near the box edges): quantitatively reproduces the
   semi-empirical plateau via E = 2πQ_enc/A. Dominant.
2. **Grid-non-commensurate analytic top-hat** (−0.58 e): fully explains the analytic curve's
   opposite plateau. Independent, fixable artefact.
3. w parameter, plate thickness (at native dz), symmetrisation, convolution, gauge: **excluded**.
