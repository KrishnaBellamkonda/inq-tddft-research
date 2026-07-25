# Plan: CAP-in-jellium baselines

Spec (with resolved decisions in XML tags): `docs/prompts/cap_in_jellium/baseline_runs.md`.
Glossary: "CAP in jellium baselines" section of `CONTEXT.md`.
Status: design LOCKED via grill 2026-06-17; implementation starting.

> Every absorption / ε number stays **PROVISIONAL** until the inq-study engine
> regression (Task #7) passes. This is the first use of the inq-study CAP in an
> *interacting* (LDA Hartree+XC) bath — the pilot is the interacting-RT sanity gate.

## Goal

Establish credible CAP baselines in r_s=5.69 jellium (L=50, N=162) before the
later energy sweep, by characterising — not minimising — how a two-sided sin² CAP
(20 Bohr total, η=−0.5) perturbs the system in three settings:

- **B0** GS reference (CAP-free, reused `gs_L50_cubic_N162_dx0p40`).
- **B1** CAP on, no projectile → bath-drainage reference.
- **B2** CAP + classical σ=0.5 electron, 100 eV.
- **B3** CAP + σ=0.5 Gaussian WP, 100 eV.

B1–B3 share one **140-a.u. (~7000-step, dt=0.02) ETRS** window so B1 is the exact
subtraction reference for B2–B3.

## Geometry (locked)

Box 50 Bohr, z∈[−25,+25]. CAP slabs [−25,−15]∪[+15,+25] (10 Bohr/side, η=−0.5).
Free region [−15,+15]=30 Bohr. Launch z₀≈−13 (4σ inside the −z CAP edge), +z,
exit through far CAP. Boundary rule now references CAP edges (±15), not box faces.

## New code (each ships a known-case test before the long runs)

| Item | Location | Kind | Test |
|---|---|---|---|
| E-field kernel (idea 1) | `inqview/analysis/efield.py` | derived (numpy) | uniform→0; sinusoid analytic; Gaussian erf |
| E-field viz | `inqview/visualisation/` | viz | theme/geometry only |
| current-density VTI writer | `inqkit/io/` | primary | round-trip + ∫J = global current |
| flux reducer ∮J_z | `inqkit/screens/` | primary | uniform-flow analytic |
| region N(t) (free + slabs) | `inqkit/` reducer | primary | box-integral = N total |
| per-orbital energy ⟨ψᵢ|H|ψᵢ⟩(t) | run.cpp/inqkit | primary | GS → εᵢ |
| total-system n(k)=Σfᵢ|FFT ψᵢ|² | inqkit | primary | GS jellium → shell peaks |

Observable enumeration per baseline + cadences: see the prompt's
`<observable_enumeration>` tag. Plane/flux screens: 9 planes z∈{−20…+20} step 5.

## Build order

0. **This plan.**
1. **E-field kernel** + formula-validation agent + 3 known-cases (pure Python,
   engine-independent → first).
2. **inq-study build sanity**; verify `inq/`↔`inq-study` byte-identical except the
   sanctioned CAP edits (self_consistency complexify, reduce.hpp CUB).
3. New C++ observables (table above), each with a `code-test` known-case +
   `docs/validation/test-catalogue.md` row.
4. `scripts/cap_baselines/run.cpp` (env-driven b1/b2/b3, built once vs inq-study) +
   `dispatch.py` (2-GPU, `cudaMemGetInfo` probe, email per baseline) + **100-step
   B1 pilot**: confirm s/step, energy stays real, total E + norm decrease smoothly
   (absorption signature, no NaN), GS loads, slabs absorb.
5. Long runs B1→B2→B3 across both GPUs; emails per baseline; auto-build
   `hypotheses/cap_baselines/cap_baselines_study.ipynb` (notebook-making skill);
   final email.

## Validation gates

- Pilot (step 4) gates the long runs (Tier-A).
- Each new observable: known-case test locked before use (code-test skill).
- E-field: formula-validation agent + user agreement before "locked".
- Continuity / CAP-sink cross-check (dN_free/dt vs ∮J·dA) as a physics sanity
  check on B1.

## Key references

- inq-study CAP: complexified `vscalar` (self_consistency.hpp:182), GS bit-identical.
- `observables::current_density` (inq/src/observables/current.hpp:26) — native J(r)
  field incl. nonlocal [r,V_nl] term; do NOT reconstruct from orbitals.
- De Giovannini, Larsen & Rubio 2014 (sin² CAP form).
- Vacuum predecessor: `docs/handovers/absorbing-boundary.md` (ETRS-not-CN; η=−0.5/L=20).
