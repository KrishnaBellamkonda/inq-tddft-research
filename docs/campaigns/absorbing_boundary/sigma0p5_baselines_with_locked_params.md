---
id: ab-sigma0p5-baselines
area: absorbing_boundary
title: "sigma=0.5 production-packet absorber baselines (locked CAP/MFA)"
status: ready
hypothesis: "The locked (L*, eta*) absorber set yields a usable survival-fraction baseline eps(E) for the sigma=0.5 dispersing production packet."
handover: docs/handovers/absorbing-boundary.md
tasks:
  - { name: "sigma=0.5 baseline sweep (CAP + MFA, E grid)", done: false }
  - { name: "baseline table + study notebook", done: false }
blocked_reason: ""
---

# Production-packet (σ=0.5 Bohr) baselines using the locked CAP/mask parameter set

> **Status: PRELIMINARY.** This prompt is written *before* its predecessor task
> has finished. It consumes the "comfortable region" parameter set that the
> two-sided CAP-vs-mask reflectivity sweep
> (`docs/plans/<twosided-cap-vs-mask>.md`, this task's successor) is going to
> **lock in**. Do not start until that set exists; fill the `LOCKED:` slots below
> from its handover/notebook first.

<identity>
You are a scientific-computing researcher working on first-principles TDDFT
simulations in this repository. You use `inq-study/` for any engine-level CAP work
and never touch `inq/`. You follow the repo rules (file-placement, validation
gates, scientific grounding, commit-message, notebook-making).
</identity>

<description>
The predecessor task characterised the **sin² CAP** and the **mask-function
absorber (MFA)** with a *quasi-monochromatic* packet (σ = 4√2/k₀, ~12 % momentum
spread) and produced clean reflectivity curves ε(E) over a width L and depth η
sweep, then **locked a "comfortable region"** — the (L*, η*) that holds ε low
across the energy band. That locked set is the absorber configuration intended for
production.

This task asks the complementary, production-facing question: **how does that
locked absorber behave for the packet I actually launch in my research runs — a
sharply localised σ = 0.5 Bohr Gaussian — across energies?** The answer is the set
of **baselines** I will quote when running real systems.

Critically, σ = 0.5 is a *different physical regime* from the characterisation
packet, and the predecessor task established why (carry these facts forward; do not
re-derive from scratch, but DO verify them):

- A σ=0.5 packet is **far from monochromatic**: Δk/k₀ = 1/(√2·σ·k₀) ≈ 165 % at
  10 eV, 52 % at 100 eV, ~17 % only by ~1000 eV. The "energy" label is the central
  launch energy E = ½k₀², not a sharp line.
- It **disperses violently**: a free Gaussian spreads as σ(t)=σ₀√(1+(t/2σ₀²)²);
  with σ₀=0.5 the spreading time is ~0.5 a.u., so by the time it drifts to a
  boundary it has fanned out to **tens–hundreds of Bohr** at low/mid energy
  (σ@arrival ≈ travel/k₀). **No box size contains this** — travel scales with the
  box, dispersion scales with travel.
- Therefore ε(E) here is **not** monochromatic beam reflectivity. It is the
  **net surviving (un-absorbed) fraction** of the real dispersing packet — a
  property of *absorber-applied-to-your-workload*. The **two-sided** absorber
  (CAP/MFA on BOTH z-boundaries) is **physically required**: the large k<0
  amplitude at low E sends real flux backwards into the near boundary.
- 1 eV is infeasible (dispersion overflows any practical box); start at ~3 eV.
</description>

## LOCKED inputs (fill from the predecessor task before running)

- `LOCKED: L* (total absorber width, Bohr)   = ____`
- `LOCKED: η* (CAP depth, Ha, near −0.5)      = ____`
- `LOCKED: mask L for MFA (Bohr)              = ____`
- `LOCKED: box / launch geometry rule         = ____`  (two-sided; WP launched a
  margin inside one boundary with negligible density in the CAP region, moving
  toward the far boundary)

## Tasks

<task>
<name>σ=0.5 production-packet baseline sweep</name>

Using the LOCKED absorber set unchanged, run the σ = 0.5 Bohr packet at a grid of
launch energies (suggested 3, 10, 30, 100, 300, 1000 eV — confirm in a grill).
For each energy and for BOTH absorbers (sin² CAP at η*, MFA):

1. Box sized moderately (cost/artifact knob — full containment is impossible);
   two-sided absorbers; τ long enough that ε converges (all flux has reached an
   absorber). Verify convergence of ε vs τ on at least one energy.
2. Record the **survival fraction** ε(E) for CAP and MFA, the **time-to-absorb**
   (t at which inner-region norm falls below a chosen threshold), and a density
   GIF of one representative energy each.
3. The notebook (per the `notebook-making` skill) MUST state up front that the
   y-axis is *survival fraction of a σ=0.5 dispersing packet*, not monochromatic
   reflectivity, and contrast it against the predecessor's clean ε(E) at the same
   (L*, η*) so the gap between "absorber property" and "workload behaviour" is
   explicit. These ε remain PROVISIONAL until the inq-study engine regression
   (Task #7) passes.
</task>

## Deliverables
- A sweep folder + `hypotheses/<sweep>/` study notebook (ADR 0007), auto-built via
  the dispatcher tail (notebook-making skill).
- A short baseline table (E × absorber → survival ε, t-absorb) I can quote in
  production-run planning.
