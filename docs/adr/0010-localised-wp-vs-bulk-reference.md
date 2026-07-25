# ADR 0010 — Localised-slab WP S(E) overlaid on bulk classical + Lindhard references

- **Status:** accepted
- **Date:** 2026-06-26
- **Context scope:** localised jellium quantum-stopping-power sweep (`qsp_phase5`)

## Context

Phase 5 measures the quantum (wavepacket) electronic stopping power S(E) of a
σ_WP=0.5 projectile in the **localised jellium slab** (finite 25-Bohr slab,
two-sided CAP, deposit/L_z method) across v ∈ {1.3, 2.0, 3.0, 4.0, 5.0, 6.0}.

The user wants this plotted against a classical reference at the same width and a
linear-response reference. The only existing matched-width classical S(v) **curve**
is the **bulk** jellium `06_sigma_convergence` sweep (continuous-traversal,
slope-method), and the canonical linear-response reference (`lindhard_elf
.stopping_power_point`) is likewise a **bulk** point-charge result. The localised
classical run-set has exactly **one** reliable point (the phase-4 park-method
v=2.0 = 0.25 eV/Bohr); building a localised classical *curve* would mean re-running
the park-method classical at every swept velocity.

Two further facts shape the decision:

- **σ convention (√2 trap).** The bulk figure's σ is the *charge std* σ_q, with
  σ_WP = √2·σ_q. So the σ_WP=0.5 classical curve is the **σ_q=0.354 (`sigma0p35`)**
  set — *not* the file named `sigma0p5` (which is σ_WP=0.707). See the
  `sigma-wp-convention` rule.
- **Density matches.** Bulk r_s=5.69 ≈ localised slab r_s=5.666, so density is not
  a confound.

## Decision

Overlay the **localised-slab WP S(E) points on the BULK classical (σ_WP=0.5 =
σ_q=0.354) and BULK point-charge Lindhard references**, all three labelled
explicitly as *bulk reference*, and add the single localised park-method classical
point (v=2.0) as a marked **geometry-matched check**. Do **not** re-run a localised
classical curve.

## Consequences

- **Positive.** Reuses validated, already-computed references; gives an immediate
  S(E) estimate across the full velocity range (classical curve to v=3.0, Lindhard
  analytic beyond); keeps the sweep to 5 WP runs (~10 h on 2 GPUs) instead of
  doubling it with classical reruns.
- **Negative / caveat (the reason this ADR exists).** Slab-WP vs bulk-classical is a
  **geometry mismatch**: a finite-slab deposit/L_z number compared against a
  continuous-traversal slope number. It is an **estimate**, not an apples-to-apples
  S(E). Every plot/email/notebook states this; the lone localised park point is the
  only strictly geometry-matched comparison.
- WP points are reported with a convergence flag: high-v converge to true values;
  slow points (and the 54 eV reuse) are **upper bounds** (the WP is not fully
  absorbed by τ).

## Alternatives considered

- **Re-run a localised classical curve** (geometry-matched). Rejected for cost: it
  doubles the run count and adds classical build/dispatch to the autonomy loop, for
  a comparison whose headline conclusion (WP ≫ classical, ≫ Lindhard) the bulk
  reference already shows.
- **Drop the classical overlay, keep only Lindhard.** Rejected: the matched-width
  classical curve is the most physically informative comparison the user wants.
