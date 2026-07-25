# Plan: two-sided absorber study — sin² CAP vs mask, comfortable-region search

Status: DESIGN LOCKED via `/grill-with-docs` (2026-06-16). Predecessor context:
`docs/plans/cap-thin-absorber-tuning.md`, `docs/plans/cap-monomial-inq-study.md`,
`docs/handovers/absorbing-boundary.md`. Successor (deferred): the σ=0.5 baseline
prompt `docs/prompts/absorbing_boundary/sigma0p5_baselines_with_locked_params.md`.

## Goal

Characterise the **built-in sin² CAP** against the **mask-function absorber (MFA)**
in a **two-sided** geometry (absorber on BOTH z-boundaries) using a clean,
quasi-monochromatic packet, and from the resulting reflectivity maps **let the user
lock a "comfortable region"** — a width L (and CAP depth η) that holds ε low across
the energy band. That locked set feeds the deferred σ=0.5 production-baseline task.

The user owns the comfortable-region verdict: the notebook presents the maps and
curves; it does NOT auto-declare (L*, η*).

## Physics setup (LOCKED)

- **Packet:** energy-scaled `σ = 4√2/k₀` (the previous method) ⇒ constant ~12 %
  relative momentum spread, quasi-monochromatic at every energy. Box scales with
  energy (`Lcell_z = 6σ + L`), so there is no dispersion/overflow wall.
  *(σ = 0.5 Bohr was considered and explicitly deferred — it disperses violently and
  yields survival-fraction, not reflectivity; that is the successor task.)*
- **Two-sided geometry:** total absorber width **L split L/2 at each z-end**
  (e.g. L=20 ⇒ 10 Bohr each end). The forward-moving packet only ever meets the
  **far L/2**; the near absorber catches back-scattered / k<0 flux. ⇒ the L grid
  must run higher than the single-sided sweeps (a two-sided L is ~half as effective
  on the beam as a single-sided L of the same number).
- **Launch:** WP injected a margin inside ONE boundary (negligible density in that
  boundary's CAP region at t=0), moving toward the far boundary.
- **ε (figure of merit):** `ε = inner-region norm(|z| < inner edge) / N₀` at the
  stop time — surviving un-absorbed fraction; lower = better. Same definition as the
  predecessor sweeps (`inqkit::absorbers::inner_region_norm`).
- **t_absorb:** first time the inner-region norm falls below **1 % of N₀** — how
  *fast*, complementing how *well*.

## Implementation (LOCKED)

- **Two-sided CAP — no new engine code.** Compose two built-in slabs via
  `perturbations::sum` / `operator+`:
  `cap = absorbing(η·Ha, mid_hi, w_half) + absorbing(η·Ha, mid_lo, w_half)`,
  with `w_half = (L/2)/Lcell_z` and `mid_hi/mid_lo` the fractional centres of the
  two end-slabs (the `absorbing` ctor already wraps fractional `mid_pos`). Reuses
  the existing (Task #7-gated) `inq-study` `absorbing` class; `inq/` untouched.
  - Test: task-specific **mechanism check** in
    `hypotheses/twosided_cap_vs_mask/tests/` — confirm both ends absorb and a
    symmetric two-sided launch gives left/right-symmetric ε.
- **Two-sided mask — small inqkit addition.** The current
  `inqkit::absorbers::MaskAbsorber` ramps only the high-z orientation; add a
  **two-sided variant** (M=1 for `|z| < z_abs0`, ramps to 0 at both walls), applied
  in the propagate callback after each ETRS step. Wrapper-level, no engine
  dependency.
  - Test: library-generic **feature test** in
    `inq-stack/tests/include/inqkit/absorbers/` (mask shape symmetric, M(inner)=1,
    M(wall)=0).
- **Propagator:** ETRS (absorption requires it; CN renormalises and cancels it).
- **Provisional:** all CAP ε remain PROVISIONAL until the `inq-study` engine
  regression (Task #7) passes; the mask ε do not depend on it.

## Sweep grids (LOCKED)

- **E (eV, log):** 1, 2, 4, 8, 16, 32, 64, 100 *(main 10⁰–10²)* + 300, 1000 *(bonus)*.
- **L total (Bohr, ½ each end):** 10, 16, 20, 26, 30 (per-end 5, 8, 10, 13, 15).
- **η (Ha, CAP only — mask has no depth knob):** −0.3, −0.5, −0.7, −1.0.
- **Anchor E:** 10 eV (ε-vs-L cut and "best" GIF selection).

## Run structure — 3 phases, each emailed on completion (~140 runs, 2 GPUs)

1. **CAP η-sweep** at fixed L=20: E × η  (≈40 runs) → email **ε(E) | η**.
2. **CAP L-sweep** at η=−0.5: E × L  (≈50 runs) → email **ε(E) | L (CAP)**.
3. **Mask L-sweep:** E × L  (≈50 runs) → email **ε(E) | L (mask)** + side-by-side
   + best-CAP/best-mask density GIFs.

Email = key plot(s) per phase via `inqview.email.send_run_email` (threaded). η=−0.5
is the working depth for phases 2–3 (the predecessor `cap_Lopt_E10` found −0.5
optimal); phase 1 shows whether that holds across energy — the user may re-pick.

## Results catalogue (LOCKED — 7 items)

1. ε(E) at different η (CAP), fixed L.
2. ε(E) at different L — CAP and mask (two panels).
3. Side-by-side CAP vs mask (shared axes/colorbar; shared-colorbar rule).
4. Density GIFs of the best CAP and best mask run (lowest ε at anchor E).
5. Time-to-absorb: t_absorb per run and t_absorb(E) for CAP and mask.
6. ε vs L at the anchor energy (CAP at η*, mask) — the "how thin can we go" cut.
7. ε(E, L) heatmap per absorber (CAP, mask; shared clim).

(No auto-locked comfortable-region box — the user reads (L*, η*) off the maps.)

## Notebook

`ResearchProject/systems/vacuum/hypotheses/twosided_cap_vs_mask/` study notebook,
authored per the `notebook-making` skill (context → formulas-with-terms →
reconstructable setup → linked source files → results → takeaway), auto-built by the
dispatcher tail once at end of the batch. Figures `.png`; GIFs for the two best runs.

## Caveats to state in the notebook

- ε is clean reflectivity here (σ ∝ 1/k₀, ~12 % spread) — contrast with the deferred
  σ=0.5 survival-fraction task.
- The near absorber mostly catches back-scatter; with a clean forward packet its
  contribution is small but it completes the two-sided geometry the production runs
  use.
- 1 eV uses a ~130 Bohr box (energy-scaled) — feasible but the largest cell.
- CAP ε PROVISIONAL until Task #7.

## Deferred / not in scope

- σ=0.5 production baselines → successor prompt (already written).
- Monomial CAP — dropped (predecessor showed n=1 ramp ≈ sin² at L+1 Bohr; not
  competitive enough to keep in this comparison).
