# Rule: σ means the wavepacket width (σ_WP)

Apply to: `docs/presentations/`, `docs/reports/`, `ResearchProject/systems/**/hypotheses/`,
`ResearchProject/systems/**/shared/configs/`, run `analyse.py`, plot/figure
scripts (`build_*.py`, `make_*_postproc.py`), notebook builders, and any slide,
caption, axis label, or table that reports a Gaussian width. Always on.

## The one rule

**Whenever a Gaussian width σ is named, labelled, plotted, or tabulated, it is the
wavepacket width σ_WP** — the standard deviation of the projectile *wavepacket*,
not of the classical Gaussian potential.

- The classical Gaussian **potential**'s width is a *derived internal* quantity
  `σ_pot = σ_WP/√2` (so the classical charge std equals the WP density std). It is
  used only when generating the UPF. **Never relabel a run, curve, or axis by
  σ_pot.**
- A classical run "matched" to a WP run is reported at its **σ_WP** (the shared
  label). Surface `σ_pot = σ_WP/√2` only in a methods footnote when the UPF
  generation is being described — never as the headline σ.
- Example: the σ_WP = 0.5 Bohr pair uses `electron_gaussian_sigma0p35.upf`
  (σ_pot ≈ 0.354). Both the WP run and its classical ghost are labelled **σ = 0.5**.

## Why

A single σ axis must mean exactly one thing across the whole thesis and every
presentation, or WP and classical curves are silently shifted by √2 and become
non-comparable. Standardising on σ_WP keeps every S(v)/S(E) overlay, σ-sweep, and
caption on one consistent scale. (User decision, 2026-06-25, for the supervisor
presentation and all downstream work.)

See [[reference_sigma_matching_convention]] for the √2 derivation and UPF
generation mechanics.
