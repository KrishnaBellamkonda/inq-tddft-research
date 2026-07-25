# Plan — Wide-WP vs classical, full PBC comparison (localised jellium)

**Goal (user, 2026-07-06):** rerun the `wp_per2_E300_long` configuration but under
**full 3D periodic boundary conditions** (periodic in z too), for a **matched pair** —
a wavepacket projectile AND a classical Gaussian projectile at the **same Gaussian
width** (σ_WP = 3.5 Bohr). Then build **two run notebooks** + **one comparison
("phase") notebook** contrasting WP vs classical.

Motivation: the open-z `wp_per2_E300_long` non-plateau was shown to be a phantom
absorbed-orbital bookkeeping artifact ([[reference_phantom_absorbed_wp_orbital_energy]]),
NOT a boundary effect. PBC removes the open-z monopole (E_GS → ~−86 Ha, clean) and is
~2× cheaper (no 2D-Poisson z-doubling). This run also cross-checks that the phantom
artifact persists under PBC (it should) while giving the WP−classical quantum-vs-classical
stopping comparison the campaign wants.

## Fixed configuration (shared by both runs)
- Box 50×50×111 Bohr, dx=0.40, slab half-width 12.5 (axis z), N=82, r_s=5.667.
  Config header `shared/configs/slab_n82_L50x50x111.hpp` (geometry constants reused).
- **Periodicity 3 (full PBC)** — the ONE change vs `wp_per2` (was periodicity 2).
- CAP: two-sided sin², **η=−1.0 Ha, 14 Bohr/side**, inner faces ±41.5 (MID 48.5/111,
  WIDTH 14/111) — matched across WP and classical.
- Projectile: E=300 eV, k0/v0=4.696 a.u., launch_z=−26.5, +z. σ_WP=3.5.
- τ ≈ 70.8 a.u. (same as `wp_per2`).

## The two runs
| | WP | Classical |
|---|---|---|
| projectile | injected Gaussian wavepacket (KS orbital), σ_WP=3.5 | Ehrenfest Gaussian-e ion, UPF `electron_gaussian_wpsigma3p5.upf` (σ_pot=2.475=σ_WP/√2) |
| dt | 0.04 | 0.02 (Ehrenfest stability) |
| N_STEPS | 1773 | 3540 (matched τ) |
| S extraction | pre-absorption momentum-centroid −dKE_WP/ds (primary); re-ledgered ΔE plateau (cross-check) | initial-drag −dKE_ion/ds over early v≥0.85·v0 window (electron_track.csv) |
| GPU | 0 | 1 (concurrent) |

Both label σ = **3.5** (σ_WP convention; σ_pot only in methods footnote).

## Steps
1. **New GS** (PBC-111): `shared_gs/slab_n82_L50x50x111_h0p40_pbc`. Expect E_GS ≈ −86 Ha
   (clean, no monopole). ~40 min, GPU. Gate: SCF converged, n0/r_s correct, 82 e⁻.
2. Edit + rebuild 3 binaries (periodicity, GS dir, classical config+CAP). Compile-gate.
3. Launch WP (GPU0) + classical (GPU1) concurrently, detached, liveness-guarded.
4. On both complete: build 2 run notebooks (run-notebook skill) + 1 comparison notebook
   in `hypotheses/wide_wp/`. Email on completion/failure.

## Validation gates
- Tier A: GS SCF converged; both RT runs start (E finite, no NaN); N drift sane.
- Physics sanity: classical initial-drag S(300 eV) finite; WP centroid drift extractable;
  compare to open-z result and Lindhard high-v tail (S≈0, <1 eV over slab).
- Per [[reference_light_projectile_deceleration]]: gate on clean initial-drag slope, NOT
  v-drift; both are light projectiles that decelerate.

## Files
- Edits: `scripts/wide_wp/{gs,wp,classical}/run.cpp`.
- New: orchestrator `scripts/wide_wp/run_pbc_pair.py`; notebooks in `hypotheses/wide_wp/`.
- No `inq/` edits (immutable); engine = inq-study.

## Open / deferred
- The real code fix for the phantom orbital (norm-weight state-60 kinetic in the energy
  assembly) is a separate inqkit task — this run will exhibit the artifact and we
  re-ledger in post (drop state 60). Not blocking the comparison.
