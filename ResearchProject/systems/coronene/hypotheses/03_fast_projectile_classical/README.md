# 03 — Fast-projectile classical limit

## Hypothesis

A fast (E = 800 eV), tightly-localised (σ = 0.33 Bohr ≪ atomic scales) WP
should behave classically: it should propagate ballistically through the
molecule with minimal hybridisation, the WP density should remain spatially
compact for the duration of the run, and the LEED pattern should be
dominated by single-shot diffraction off the molecular potential rather
than by stationary-state hybridisation.

This is the high-velocity, weak-coupling limit of the scattering problem
and is the regime in which a classical Coulomb trajectory should be a
reasonable approximation.

## Runs collated

- `run_E800_s0p33` — E = 800 eV, σ = 0.33 Bohr.
- `run_base`        — paper reference, for contrast.

## Produce comparison artefacts

```
python scripts/coronene_postprocess.py hypothesis \
  --hypothesis-dir hypotheses/03_fast_projectile_classical \
  --runs run_E800_s0p33=$(realpath run_E800_s0p33) \
         run_base=$(realpath run_base)
```
