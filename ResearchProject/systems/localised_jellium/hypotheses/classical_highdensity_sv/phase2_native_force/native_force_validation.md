# Phase 2 — native analytic Hellmann–Feynman force validation

Validates `inqkit::dynamics::projectile_force_analytic` — the INQ-native local
HF force on a Gaussian-charge projectile,

    F = − ∫ V_proj(r−R) · ∇n(r) dr,   V_proj = poisson(gaussian_density(σ_pot)),

built byte-for-byte from the density-gradient integrand INQ uses for the local
ionic force (`inq/src/observables/forces_stress.hpp:182–187`). σ_pot = 0.35355
(= σ_WP/√2 for σ_WP = 0.5; the width axis is labelled by σ_WP = 0.5).

Two binaries under `scripts/classical_highdensity_sv/`:
`force_native_analytic/` (Test 1) and `force_vs_native/` (Test 2). Both build
clean (no compile fixes needed — see verdict). All runs finite box (periodicity
0) ⇒ free-space Poisson, GPU (LDA), venv post-processing.

Figure: `native_force_validation.png`.

## Test 1 — closed-form two-Gaussian force

Fixed Gaussian "density" (σ_s = 0.5, ∫ = 1) at z = 0; projectile Gaussian
(σ_pot = 0.35355) swept z_c = 1…14 (dx = 0.5, L_z = 85). Compare `F_ours`
(analytic operator) and `F_fd` (existing finite-difference operator
`projectile_force_z`) to the closed form

    F(d) = erf(a·d)/d² − (2a/√π)·e^(−a²d²)/d,   a = 1/(√2·√(σ_pot²+σ_s²)) = 1.15470.

Sign: both operators return **+F(d)** (repulsive), ratio ≈ **+1** (no flip).

| quantity | ours / analytic | fd / analytic |
|---|---|---|
| median ratio | **1.00084** | 1.00124 |
| median \|ratio−1\| | **0.11 %** | 0.33 % |
| max \|ratio−1\| | **0.60 %** | 0.90 % |
| max shape dev (\|ratio/median−1\|) | **0.51 %** | 1.02 % |

The analytic force is tighter than the FD (as expected: no finite-difference
truncation), matching the closed form to grid accuracy. **Test 1: PASS.**

## Test 2 (DECISIVE) — our analytic force == INQ's native force on a ghost UPF

A **clean** ghost UPF `ghost_sigma0p354.upf` was generated with
`inqview.io.gaussian_psp.generate_gaussian_psp(template, sigma_wp=0.5)` from the
`electron_gaussian_wpsigma0p5.upf` template (NOT the legacy flipped-tail files).
Its `PP_LOCAL` is `V_loc(r) = C·erf(r/(√2·σ_charge))/r` with σ_charge = 0.35355,
C = 2 (template Rydberg tail: V·r → 2 at large r; INQ converts Ry→Ha on read).
Verified positive/repulsive everywhere and equal to the reference C·erf/r to
machine precision on the mesh — **same sign as our V_proj = +poisson(n_proj)**
for a −1 projectile.

GS system: neutral **He** atom (real nontrivial density) at the origin + one
**ghost H** ion at {0,0,z_c}. Because the ghost has z_valence = 0 and 0
projectors, `result.forces[ghost]` is a **pure local HF force** −∫V_loc·∇n —
INQ's own operator (orbital-sum gradient Σ occ·∇|φ|²). We compare it to
`projectile_force_analytic(electrons.density(), …)` (FFT-gradient of the summed
density). Finite box, LDA.

### dx = 0.5 (Bohr)

| z_c | F_INQ.z | F_ours.z | ratio (ours/INQ) |
|---|---|---|---|
| 3.0 | 0.208211 | 0.206112 | 0.98992 |
| 3.5 | 0.156208 | 0.158222 | 1.01289 |
| 4.0 | 0.121560 | 0.120049 | 0.98757 |
| 5.0 | 0.078845 | 0.077694 | 0.98541 |

mean ratio 0.9939, max |ratio−1| = 1.46 %, rms 1.26 %. The x,y components are
0 to 1e-10 (on-axis) for both. The sign of (ratio−1) oscillates (−,+,−,−) ⇒
not a systematic scale/sign error but grid noise.

### dx = 0.4 (Bohr) — convergence check

| z_c | F_INQ.z | F_ours.z | ratio (ours/INQ) |
|---|---|---|---|
| 3.0 | 0.208562 | 0.208535 | 0.99987 |
| 3.5 | 0.156774 | 0.156885 | 1.00071 |
| 4.0 | 0.121657 | 0.121571 | 0.99929 |
| 5.0 | 0.078900 | 0.078946 | 1.00059 |

mean ratio **1.00011**, **max |ratio−1| = 0.071 %**, rms 0.058 %.

Refining dx 0.5 → 0.4 collapses the discrepancy (e.g. z_c = 5: 1.46 % → 0.06 %).
This proves the residual is pure spatial discretization: INQ's orbital-sum
gradient (Σ occ·∇|φ|², sharply peaked near the He core) and the FFT-gradient of
the summed density are mathematically the same field and converge to the SAME
force as dx → 0. **No header edit is required** — the FFT-gradient form is
correct; the two forms agree to 0.07 % at dx = 0.4.

### Verdict — Test 2: PASS

`projectile_force_analytic` reproduces INQ's native local Hellmann–Feynman force
on a real ghost-UPF ion to **≤ 0.07 % at dx = 0.4** (≤ 1.5 %, non-systematic, at
the coarser dx = 0.5). A perturbation projectile of potential V_proj and a
pseudopotential ion of the same V_loc feel the **same** HF force — the operator
is validated for use in the classical stopping-power machinery.

## Compile fixes to `projectile_force.hpp`

**None.** The already-written `projectile_force_analytic` /
`projectile_force_analytic_z` compiled and ran correctly as delivered
(`operations::gradient(density)`, `vproj.linear().cbegin()`,
`gpu::reduce(basis.local_size())`, `cell.to_cartesian(force_cov)`) — no API
adjustments were needed. The gradient form was NOT switched to the orbital-sum
form because the dx-convergence test above shows the FFT-gradient already agrees
with INQ to grid accuracy.

## Provenance

- Binaries: `scripts/classical_highdensity_sv/force_native_analytic/run.cpp`,
  `scripts/classical_highdensity_sv/force_vs_native/run.cpp` (finite box, GPU/LDA).
- Ghost UPF: `force_vs_native/ghost_sigma0p354.upf` (also copied here);
  σ_charge = 0.35355, C = 2 (Rydberg), erf/r verified.
- CSVs here: `force_native_analytic.csv`, `force_vs_native_dx0p5.csv`,
  `force_vs_native_dx0p4.csv`.
- Header under test: `inq-stack/include/inqkit/dynamics/projectile_force.hpp`.
