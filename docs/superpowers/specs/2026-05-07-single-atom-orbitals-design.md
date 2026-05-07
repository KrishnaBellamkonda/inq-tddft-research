# Single-atom orbitals — design

**Date:** 2026-05-07
**Author:** Hand-off from interactive brainstorm
**Status:** Approved by user, ready for implementation plan

## Goal

Produce ground-state Kohn–Sham orbitals for three isolated atoms (H, Li, Al)
in a vacuum cubic box, suitable for visual inspection of the s/p/d shell
structure in ParaView. Each atom gets its own self-contained subfolder under
`Tutorial/single-atom-orbitals/`, and each run writes the total electron
density plus the complex wavefunction and `|ψ|²` for *every* computed
Kohn–Sham state (occupied + extras).

## Non-goals

- No real-time propagation, no perturbation, no spectroscopy.
- No spin-polarised DFT — restricted KS suffices for orbital shapes.
- No automated ParaView render pipeline; the `results/` VTI series is the
  hand-off point. The user drives ParaView interactively.
- No comparative analysis between atoms; this is a setup-and-write tutorial,
  not a benchmark.

## Requirements

1. **Three runs**, one per atom, each in `Tutorial/single-atom-orbitals/<sym>/run.cpp`
   where `<sym>` ∈ {`h`, `li`, `al`}.
2. **Cubic finite cell, L = 30 bohr**, atom at the centre `(L/2, L/2, L/2)`.
3. **LDA, 60 Ry cutoff**, gamma-point only.
4. **`extra_states(30)`** so every run has at least 30 empty states above the
   highest partially-filled level.
5. **Smearing** via `.temperature(0.001_Ha)` to make SCF robust against the
   odd-electron occupation pattern (1, 1, 3 valence electrons respectively).
6. **Broyden mixing**, `energy_tolerance(1e-6_Ha)`, `max_steps(1000)`.
7. **Output for every state** `i ∈ [0, n_states)`:
   - `results/orbital_density/orbital_<i:04d>_density` — `|ψ_i(r)|²` (real).
   - `results/orbitals/orbital_<i:04d>` — `ψ_i(r)` (complex).
8. **Output once per run:**
   - `results/density/density_total` — total electron density.
   - stdout: GS total energy, density-norm sanity check (∫ρ ≈ N_valence),
     and the eigenvalue ladder (one line per state with the eigenvalue and
     occupation, for shell-label inspection).

## Atom-specific parameters

| Atom | Z (pseudopotential valence) | Expected `n_states` |
|---|---|---|
| H  | 1 | 31 (1 partially filled + 30 extras)  |
| Li | 1 (1s² in core) | 31 |
| Al | 3 (Ne core) | 32 (ceil(3/2)=2 + 30) |

## Architecture

Each `run.cpp` is a stand-alone INQ program of ~80 lines:

```
main()
├── build cell + ions
├── build electrons (cutoff, extra_states, temperature, gamma)
├── ground_state::initial_guess + ground_state::calculate (LDA, Broyden)
├── print energies and per-state eigenvalues + occupations
├── write total density (RealField3DWriter)
└── for i in 0..n_states-1:
        write |ψ_i|² (RealField3DWriter)
        write ψ_i    (ComplexField3DWriter)
```

`n_states` is read at runtime from `electrons.kpin()[0].spinor_set_size()`
(matches the runtime check pattern already used in `inqkit/fields/orbital.hpp`
and `inqkit/fields/density.hpp`).

## File layout

```
Tutorial/single-atom-orbitals/
├── h/run.cpp
├── li/run.cpp
└── al/run.cpp
```

Each `run.cpp` is built and run via `inq-run` from inside its own folder
(produces `build/`, `run`, `results/`). No shared library code is added —
the three programs are intentionally self-contained tutorial files.

## Validation

Per-run, scripted inside `run.cpp`:

- SCF must reach `1e-6 Ha` energy tolerance (`gs.energy.total()` printed).
- ∫ρ(r) d³r ≈ N_valence (1, 1, 3) — printed via `operations::integral`.
- `n_states` printed and asserted ≥ 31 for H/Li and ≥ 32 for Al.

User-driven, post-run:

- ParaView inspection of the orbital VTI series — H should show
  `1s, 2s, 2p×3, 3s, 3p×3, 3d×5, …` ladder; Li similar with the
  pseudopotential's tighter 2s; Al should show 3s, 3p×3 partially filled
  followed by the same 3d/4s/4p ladder.

## Risks / open questions

- **Pseudopotential availability.** INQ's default pseudo-set must contain
  H, Li, and Al. If any is missing, fall back to `pseudo::set::pseudodojo_pbe()`
  (note: that's PBE pseudos used inside an LDA calculation, which is
  acceptable for orbital-shape work but documented in the run.cpp comment).
- **Particle-in-a-box vs. Rydberg states.** The 30-bohr cubic finite box
  hard-walls the wavefunctions; high-lying "Rydberg-like" states at
  i ≳ 15–20 will be box modes, not true atomic Rydberg states. This is
  expected and a feature, not a bug, for the educational goal.
- **Restricted KS on odd-electron atoms.** With temperature smearing, the
  partially-filled state will have an occupation between 0 and 2; the
  orbital shape is still well-defined.
