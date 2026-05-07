# Audit: inqkit multi-kpoint compatibility

Branch: `features/inqkit-multikpoint` (worktree of main repo).
Plan: `docs/plans/inqkit-multikpoint-bloch-viz.md`.

Tag legend:
- `kpoint-summed` — already correct; INQ sums over k-points internally.
- `requires-explicit-kpoint` — accepts a `kpoint_index` argument; multi-k
  callers must pass one. Default of 0 is retained for legacy single-k
  call sites but is **deprecated for multi-k systems**.
- `single-kpoint-assumed-needs-fix` — silently reads only kpoint 0; needs
  a kpoint-aware variant before use on multi-k systems.
- `n/a` — no kpoint logic.

## inqkit/fields/

| Function | Tag | Notes |
|---|---|---|
| `density::total(electrons)` | `kpoint-summed` | Reads `electrons.density()` which is the k-summed density. Single-rank-basis assumed (line 52 throws otherwise). Verified Phase 1 of the Li 54-atom run. |
| `density::orbital(electrons, n, k=0)` | `requires-explicit-kpoint` | Accepts `kpoint_index` with default 0. 15+ legacy single-k callers in `Tutorial/`, `ResearchProject/` rely on the default; do not drop it. **Multi-k callers must pass `k` explicitly.** |
| `density::fft_shift_index` | `n/a` | Pure index helper. |
| `orbital::wavefunction(electrons, n, k=0)` | `requires-explicit-kpoint` | Complex-field analogue of `density::orbital`. Same default-of-0 deal. Throws on multi-rank basis/state, on `spinor_dim != 1`, on out-of-range indices. |
| `RealField3D` / `ComplexField3D` | `n/a` | Pure structs. |

## inqkit/observables/

| Function | Tag | Notes |
|---|---|---|
| `eigenvalue_dump.hpp::dump_eigenvalues` | `kpoint-summed` (writes per k) | Iterates `electrons.kpin()`, writes `kpoint_index, kx, ky, kz, weight, state_index, evalue, occ` long-format. Single-rank-of-kpoints assumed (`kpin_part().start()` offset documented). Multi-k correct. |
| `density_delta.hpp::DensityDelta` | `kpoint-summed` | Operates on a `RealField3D` (k-summed density), no per-k logic. |
| `center_of_density.hpp` | `kpoint-summed` | Reads `electrons.density()`. |
| `momentum_distribution.hpp` | TODO | Not audited yet — out of scope for this plan; revisit when needed for the v=0.450 run's D2/D4 bundle. |
| `state_energy_writer.hpp` | `single-kpoint-assumed-needs-fix` | Comment in `eigenvalue_dump.hpp` notes it ports per-rank logic; behaviour on multi-k unclear. **Out of scope** here; flag in handover. |
| `orbital_overlap.hpp::OrbitalOverlap` | `single-kpoint-assumed-needs-fix` | Constructs reference wavefunctions via `fields::orbital::wavefunction(electrons, i)` with no explicit `k` — silently reads kpoint 0 only. Used by coronene wavepacket runs (single-k). For multi-k metallic propagation we either pass an explicit `k`, or generalise to an N_k × N_states matrix. **Not needed for the immediate Phase-10b production-apply** (no propagation yet); flag for follow-up. |

## inqkit/wavepacket/

| Function | Tag | Notes |
|---|---|---|
| `wavepacket.hpp::WavePacket::inject_into_last_extra_state` | `single-kpoint-assumed-needs-fix` | Writes into `kpin()[0]` slot. Multi-k callers would need an explicit `kpoint_index`. **Out of scope** for the current plan — no wavepacket injection in the v=0.0123 / v=0.450 propagation runs. |
| `injection_report.hpp` | `n/a` | Plain struct. |

## inqkit/real_time/

| Function | Tag | Notes |
|---|---|---|
| `real_time_session.hpp::RealTimeSession` | `kpoint-summed` | Generic step callback — any per-k logic lives in the user-supplied lambda. No fix needed at this layer. |
| `step_context.hpp::StepContext` | `kpoint-summed` | Energies, currents, dipoles fed in are already INQ k-summed scalars. |

## inqkit/screens/

| Function | Tag | Notes |
|---|---|---|
| `plane_screen.hpp::PlaneScreen` | `kpoint-summed` | Reads `electrons.density()` slice. |
| `leed_pattern_accumulator.hpp` | `kpoint-summed` | Same. |

## inqkit/io/

All writers are field-agnostic. No kpoint logic.

## Summary of action items for this plan

1. **A2 (this plan)**: Keep `density::orbital` and `orbital::wavefunction`
   defaults; document that multi-k callers must be explicit. Add an
   inline `static_assert`-style comment near each. **Done**: this
   document is the canonical reference.
2. **B1**: New `dump_orbitals_per_kpoint` driver iterates `(band, k)`
   pairs explicitly — calls `orbital::wavefunction(electrons, n, k)`
   with both arguments specified.
3. **Out of scope, flagged for follow-up**:
   - `state_energy_writer.hpp` audit + multi-k variant.
   - `orbital_overlap.hpp` multi-k generalisation.
   - `wavepacket.hpp` multi-k injection.
   - `momentum_distribution.hpp` audit.

## Single-rank-basis assumption (orthogonal to k-points)

Every helper that walks the real-space grid in serial currently throws
when the basis communicator has size > 1. This is fine for single-GPU
runs (k-parallelism only); a multi-rank basis path is not on the
critical path for this plan.
