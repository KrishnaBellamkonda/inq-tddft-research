# ADR 0012 — Agnostic-library governing-PDE discovery, validated by three walls

- **Status:** accepted
- **Date:** 2026-07-03
- **Context scope:** `docs/campaigns/ml-patterns/pattern-finding-in-wp-classical-runs.md`
  — the bulk-jellium PDE-discovery redo (Track B) and any successor that discovers
  a governing differential equation for an induced-density field.
- **Extends:** ADR 0011 (held-out split anti-p-hacking) to the equation-discovery
  setting.

## Context

The bulk-jellium redo re-centres the campaign on **discovering a governing
differential equation** for the induced bath density `n_bath(r,t) − n_bath^GS`,
separately for a classical (point) and a quantum (wavepacket) projectile, then
comparing them. The scientific goal is *"differential equations suggestive of
some physics"* — i.e. the discovered terms should be interpretable (a plasmon
restoring term, a dispersive term, an advective term, a projectile source), not a
black-box surrogate.

Two library strategies were on the table:

- **Physics-priored library.** Seed the SINDy/PDE-FIND operator library with the
  cold-plasma / RPA terms we expect (`−ω_p²n`, `∇²n`, `v·∇n`, moving source).
  Robust and well-grounded, but it **bakes the answer in**: recovering the plasmon
  term is then near-tautological, and a genuinely novel term cannot appear.
- **Broad agnostic library.** A large generic operator set (derivatives to high
  order, products to cubic) with **minimal physics priors**; sparsity decides
  which terms live. Faithful to *discovery*, but with a high spurious-term risk —
  and, inside an autonomous refine-until-it-validates loop, a direct route to
  p-hacking a plausible-looking equation.

The user chose the **broad agnostic library**. That choice is only defensible if
the spurious-term / p-hacking risk is controlled by validation strong enough that
a term surviving it is trustworthy *regardless* of how it was found.

## Decision

Track-B PDE discovery uses a **broad, agnostic operator library with minimal
priors**, and a discovered term is admitted as "physics" **only if it survives all
three independent walls**:

1. **Pinned calibration/held-out cell split** (extends ADR 0011). The library
   scope, sparsity threshold, denoising and subtraction-ladder choices are tuned
   **only on calibration cells**; reported coefficients and the admit/reject
   decision are read from **held-out cells**. The split is pinned in the campaign
   prompt, not chosen at runtime.
2. **Temporal forward-prediction.** On a held-out trajectory, fit on the early
   time-window, **integrate the discovered PDE forward**, and score the prediction
   against the later window. A governing equation must predict dynamics it was not
   fit on — a coefficient that fits `∂ₜn` pointwise but fails forward-integration
   is rejected.
3. **Bootstrap coefficient stability.** The active-term set and coefficient signs
   must persist across resampled calibration subsets. Terms that flicker in/out
   under resampling are spurious and rejected.

**Physical interpretation is assigned post-hoc**, after the surviving sparse set is
fixed — a physics-interpreter names each survivor (does it match `−ω_p²n`? `∇²n`?)
against known plasma physics. Interpretation never feeds back into the library.

The autonomous loop **may** refine the analysis (library scope, sparsity, denoising)
to maximise validation, but only against the **calibration** split and the
**forward-prediction on the fit-window**; the admitted equation and the
classical-vs-WP comparison are always reported from **held-out** data. All refine
attempts are logged. A refute / inconclusive / partial is a valid reported outcome.

## Consequences

- **Positive.** A surviving term is a genuine discovery, not a seeded expectation;
  the recovery of a plasmon or dispersion term (if it happens) is evidence, not
  assumption. The three walls are jointly far harder to game than any one, so the
  agnostic library is safe to use inside an autonomous loop.
- **Negative / cost.** More compute per candidate equation (forward-integration +
  bootstrap over cells), and a real possibility the honest outcome is
  "no stable interpretable PDE survives" — which must be reported as such, not
  tuned away. The broad library can still surface physically meaningless survivors
  that pass numerically; the post-hoc interpreter must flag these as
  "uninterpreted term" rather than forcing a physical name.
- **Reversibility.** Low cost to revert to a physics-priored library later if the
  agnostic route proves too noisy; the walls carry over unchanged.
