# Plan — Monomial CAP in inq-study + built-in widen-L baseline

Date: 2026-06-15
Status: in progress. Part A (baseline) launched; Part B (implementation) underway.
Parent goal: reflection error ε → low-single-% / ~0 at E≈10 eV, ideally keeping a
THIN (L≈5 Bohr) absorber. Grilled outcome of `cap-thin-absorber-tuning.md`:
the built-in sin² CAP floors at ~10–14 % at L=5 (no shape knob), so ε→0 at fixed
small L needs a different absorber SHAPE.

## Background / evidence

- INQ's built-in CAP (`inq/src/perturbations/absorbing.hpp`) is **sin²-only**,
  hard-coded; constructor `(amplitude, mid_pos, width)` — no shape/order parameter,
  and it is the ONLY CAP class in the engine (verified tree-wide: zero
  `monomial`/`polynomial` hits). So a monomial CAP must be a NEW perturbation.
- Measured width-dependence (cap_real knobs, η=−0.5, E≈22 eV): ε = 14.1 % (L5),
  1.29 % (L10), 3e-3 % (L20) — ~1 order of magnitude per +5 Bohr.
- The built-in sin² is a **hump** (zero at both slab edges) → absorption vanishes
  at the box wall. A **monomial ramp** `η·s^n` (s=0 at inner edge → 1 at wall)
  puts peak absorption at the wall, where reflection happens — expected to beat
  sin² at fixed L. `n` tunes onset smoothness (reflection theory: Riss & Meyer
  1996; CAP family in De Giovannini, Larsen & Rubio 2014, arXiv:1409.1689 §IV).
- inq/ is IMMUTABLE — the new perturbation goes in **inq-study only**. It relies
  on the same scalar-potential complexification already in inq-study (makes any
  imaginary potential propagate), so ε stays PROVISIONAL until **Task #7**.

## Part A — built-in widen-L / (L,η) baseline at E=10 eV  [LAUNCHED]

`scripts/cap_Lopt_E10/dispatch.py`: fix E=10 eV (k0≈0.857), 2D sweep
L ∈ {6,8,10,12,15} × η ∈ {−0.30,−0.50,−1.00} = 15 runs, built-in sin² CAP
(reuses cap_sweep binary, no rebuild). Runs → `cap_Lopt_E10/run_cap_*`, analysis
→ `hypotheses/cap_Lopt_E10/`. Answers: the smallest L (and its η) reaching a target
ε at 10 eV — the sin² reference the monomial must beat.

## Part B — `absorbing_monomial` perturbation in inq-study

**B1. Header** `inq-study/src/perturbations/absorbing_monomial.hpp` (mirrors
absorbing.hpp): class `absorbing_monomial : public perturbations::none`,
constructor `(quantity<energy> amplitude, double mid_pos, double width, int order)`.
Potential over the fractional slab `[mid−w/2, mid+w/2]`:
`vk += i·amplitude·pow(s, order)`, `s=(z−(mid−w/2))/width ∈ [0,1]` (ramp; max at
the outer/wall edge). amplitude<0 absorbs. Comment cites the source.

**B2. Known-case test** (code-test skill): in-header Catch2 unit test — build a
small `basis::real_space`, call `potential()` on a complex field, check the
imaginary part equals `amplitude·s^order` at sample fractional points (and 0
outside the slab). Pure value check, no propagation. Runs under inq-study ctest.
Record in `docs/validation/test-catalogue.md`.

**B3. Run machinery** `scripts/cap_monomial/run.cpp` = copy of cap_sweep/run.cpp
with 3 edits: include the new header; read `CAP_ORDER` (default 2); construct
`perturbations::absorbing_monomial cap(eta*1_Ha, mid_frac, width_frac, order)`.
Identical observables/geometry/ε so it is directly comparable to the sin² runs.
Build ONCE against inq-study (validates the header compiles + engine accepts it).

**B4. Benchmark runs** at fixed L=5, E=10 eV: order n ∈ {1,2,3,4} × a small η set,
→ `cap_monomial/run_*`, analysis `hypotheses/cap_monomial/`. Compare ε_monomial(n)
vs ε_sin² (=0.14–0.20 at L=5) — does any (n,η) reach low-single-%? Then, if
promising, an n×η×E sweep to place the trough at 10 eV.

**B5. (Stretch) transmission-free CAP** (Manolopoulos, J. Chem. Phys. 117, 9552,
2002): the form designed for ε→0 at short L. Same scaffolding; defer until the
monomial benchmark is in.

## Validation gates
- Header is substantive code → ships the B2 known-case test (code-test skill).
- ε from any monomial run is PROVISIONAL until Task #7 (inq-study engine ctest).
- inq/ untouched; verify `diff -q inq/src/perturbations inq-study/src/perturbations`
  shows ONLY the new file (plus any pre-existing inq-study deltas).

## Execution checklist
1. [x] Part A dispatch + launch (built-in baseline).
2. [ ] B1 header `absorbing_monomial.hpp` in inq-study.
3. [ ] B2 in-header known-case test; catalogue row.
4. [ ] B3 `cap_monomial/run.cpp`; build once vs inq-study; smoke 1 run.
5. [ ] B4 benchmark grid at L=5,E=10; build `hypotheses/cap_monomial/` report.
6. [ ] Email comparison (sin² vs monomial vs widen-L) to chiddukanna@gmail.com.
7. [ ] Update handover `docs/handovers/absorbing-boundary.md`.
