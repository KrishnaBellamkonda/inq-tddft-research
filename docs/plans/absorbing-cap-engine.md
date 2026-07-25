# Plan — True integrated CAP via complexified scalar potential (inq-study)

Status: **active** (started 2026-06-14, auto mode)
Parent task: `docs/plans/absorbing-boundary.md` (MFA done); handover
`docs/handovers/absorbing-boundary.md`. ADR-relevant: `inq/` immutable
(`.claude/rules/inq-immutable.md`) — **all engine edits in `inq-study/` only**.

## 1. Goal

Make INQ's in-built `perturbations::absorbing` (a region-restricted sin² CAP,
`V = +i·η·sin²` in a fractional z-slab) actually FUNCTION, by **complexifying the
scalar potential** in `inq-study`, then run WP-in-vacuum depth-sweep CAP
simulations and report them richly.

## 2. Why this is needed (grounded finding, 2026-06-14)

A probe (`vacuum/hypotheses/01_cap_real/tests/cap_probe/run.cpp`) FAILED TO
COMPILE: `absorbing.hpp:45  double += inq::complex`. Root cause:
`self_consistency::update_hamiltonian` builds `vscalar` as
`field<real_space, double>` (real), `self_consistency.hpp:176`, and hands it to
`pert_.potential(time, vscalar)` (`:189`). A real field cannot hold the CAP's
imaginary part → `absorbing` is dead code (its unit test never calls
`potential()`). **Ordering consequence:** the in-built CAP cannot be exercised
*before* the complexify; the complexify is the enabler.

Propagator decision (grounded): **ETRS** (`etrs.hpp` — `exp(−iHdt)` via
`operations::exponential_*`, NO renormalisation → absorption survives). **NOT
Crank–Nicolson** (`crank_nicolson.hpp:139,147,162,165` orthogonalises +
normalises every step → undoes absorption). CN is unusable for any absorber.

## 3. Phase 1 — Complexify the scalar potential (inq-study only)

Make the scalar-potential path that flows `pert → vscalar → vks → ham` carry a
complex value, with real outputs (energy, eigenvalues) taking the real part.

Target files (currently byte-identical to `inq/` — verified):
- `src/hamiltonian/self_consistency.hpp` — `vscalar` → complex; `energy.external`
  uses `real(...)`; copy into `vks` keeps complex.
- `src/hamiltonian/ks_hamiltonian.hpp` — `scalar_potential_` field_set → complex.
- `src/hamiltonian/scalar_potential.hpp` — non-spinor apply already a complex
  multiply; ensure types compile for complex potential.
- Energy/eigenvalue consumers of the scalar potential → `real(...)` where a real
  number is required (CAP energy is artificial; excluded from reported energetics).

Method: change the types, let the compiler enumerate consumers, fix each with
real-part extraction, keep diffs minimal. Ground state (perturbation = `none`)
has zero imaginary part → must remain bit-comparable to baseline.

Build: `INQ_SOURCE=<repo>/inq-study inq-run` (inq-run targets `${INQ_SOURCE}`,
`inq-run:90`). First build is full (~15–20 min, no existing inq-study build dir).

### Phase 1 validation (code-test)
- Compiles against inq-study.
- GS of the free-WP vacuum cell unchanged vs inq/ baseline (imag part = 0).
- Smoke CAP run: WP hits the slab, total WP norm DECREASES (absorption real, not
  no-op); cross-check ε vs `cap_toy.py` at one (E,L).

## 4. Phase 2 — Depth-sweep CAP runs (built-in `perturbations::absorbing`)

WP in vacuum, CAP slab at the FAR end of the box (the end the WP travels toward),
sweep CAP **depth** η. Geometry = MFA geometry (so ε comparable). Each run ships
the **full minimum observable set** (ADR 0006) so the data is analysable:
energetics, current, dipole, WP real-space + momentum stats, density variants,
occupations/eigenvalues at final step.

Depth menu (η, Ha, negative = absorbs): {−0.05, −0.1, −0.25, −0.5, −1.0} at a
fixed mid-range E and L (e.g. E≈22 eV, L=20), plus a couple of (E) points. Exact
list finalised at launch; ETRS; dt=0.01; dx adaptive; NPERP small.
Location: `vacuum/run_cap_<E>_<L>_eta<...>/` (flat, ADR 0007).

## 5. Phase 3 — Report .ipynb

`vacuum/hypotheses/01_cap_real/` (ADR 0007): governing equation + method header,
a sample **density gif**, **energetics** plots, ε-vs-depth, **data-path links**
to each run dir, observations. Figures via canonical theme (ADR 0004); venv
python.

## 6. Phase 4 — Further inq-study runs

Iterate on the changed engine as needed (e.g. depth × width, propagator
robustness), per findings.

## 7. Risks

- Complexify blast radius (energy/mixing assume real V) — mitigate with real-part
  extraction + GS bit-comparison gate.
- ETRS Taylor stability with strong η — dx floor + watch norm trace; back off η.
- First inq-study build is slow / may surface external_libs config gaps.
