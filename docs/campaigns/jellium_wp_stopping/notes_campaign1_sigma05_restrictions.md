# Notes — Campaign 1 (quantum-stopping-power) σ=0.5 restrictions

> Recorded 2026-06-24 during the brainstorming session. These are **known
> limitations**, not blockers (user decision: proceed with Campaign 1 anyway).
> To be referenced from the eventual Campaign-1 run notebook. Sources: the
> `03_cap_stopping` baselines (`qa_jellium_slab_baselines.ipynb`) + the
> sanity-check agent (2026-06-24).

## Why σ=0.5 makes the E_total stopping ledger hard

A small-σ wavepacket is the most *point-like* (Lindhard-relevant) projectile, but
it is also the least numerically extractable with the total-energy method:

1. **Catastrophic spreading.** Free-Gaussian law `σ(t)=σ√(1+(t/σ²)²)` gives a
   **72× width growth** over 18 a.u. at σ_WP=0.5 (the measured 32× in the baseline
   was already slab/CAP-damped). The packet smears across the whole box.
2. **Stalled centroid ⇒ Δz undefined.** The survival-weighted centroid `⟨z⟩`
   decelerates and **stalls at ≈+5 Bohr** (never exits the far face at +12.5), so
   `dE/dx = ΔE/Δz` has **no well-defined Δz** — the packet neither traverses a
   clean path nor fully transmits.
3. **No-wrap and full-absorption are mutually exclusive in a 50-Bohr box.** The
   transmitted front periodic-wraps at t≈14.9 a.u. (v≈2.7), but the surviving
   **slow, low-k_z tail drains sub-linearly** and needs ≫18 a.u. to fully absorb.
   You cannot have both "no wrap" and "fully absorbed" here.
4. **Energy-ledger contamination.** `Formula 2 = E_total(final) − E_GS` equals the
   deposited energy **only at full absorption** (unreachable above). `Formula 1`
   is contaminated by the WP **zero-point KE ≈ 82 eV** (vs 100 eV drift) plus the
   **SIE ≈ 4.5 eV** (one-electron LDA self-interaction, quantified from `p3_wp`).

## Resolutions (to apply in Campaign 1, not blockers)

- **Primary estimator candidate — force/work integral** on the projectile,
  `S = (1/Δz)∫⟨ψ_WP|−∂_z V_ind|ψ_WP⟩ v dt` over the pre-CAP transit. It is *local
  to the projectile*, so it survives spreading and needs neither full absorption
  nor a box-spanning Δz. **Smoke-test on existing p5_wp data before committing.**
- **Or** a much larger z-box (≥80–100 Bohr) + longer run to let the slow tail
  drain — expensive, and Δz stays ambiguous.
- Either way: **subtract the zero-point KE and bound the SIE** (per-σ vacuum-WP
  control) before reporting any "quantum stopping" number.

## σ choice caveat (from this session)

- **σ_WP=0.5 classical ≈ point-Lindhard is validated** (localised baseline
  S=0.706 vs point 0.716 eV/Bohr = 0.99×). Good — σ=0.5 is point-faithful.
- The **analytical finite-σ Lindhard (`stopping_power_sigma`) over-suppresses**
  vs the runs (predicts 0.78× for σ_WP=0.5; the run is 0.99×), so it **cannot**
  be used to judge whether σ_WP=1.0 is "good enough." Deciding σ=1 vs σ=0.5
  needs a dedicated **σ_WP=1.0 classical S(v) run** vs point-Lindhard. See
  `brainstorming-jellium-campaigns.ipynb`.
