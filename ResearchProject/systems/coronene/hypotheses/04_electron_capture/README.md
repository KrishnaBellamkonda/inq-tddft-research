# 04 — Electron capture probe

## Hypothesis

A slow (E = 30 eV), spatially delocalised (σ = 3 Bohr) WP can dwell over
the molecule long enough for non-trivial occupancy of the molecule's
unoccupied (or weakly occupied) KS orbitals to develop. If the asymptotic
WP-overlap with one or more bound GS KS orbitals is non-negligible at
t → t_final (i.e. the bar-chart values do not tend back to zero after
the WP has nominally cleared the box), this is consistent with a
transient or persistent electron capture event.

## Runs collated

- `run_E30_s3` — E = 30 eV, σ = 3 Bohr.
- `run_base`   — paper reference, for contrast.

## Key diagnostic

`run_E30_s3/results/analysis/overlap/wp_overlap_with_gs_orbitals.gif`
should be inspected for residual occupancy near the end of the run. The
hypothesis-level comparison plots are useful for an at-a-glance side-by-
side, but the per-run overlap GIF is the primary artefact.

## Produce comparison artefacts

```
python scripts/coronene_postprocess.py hypothesis \
  --hypothesis-dir hypotheses/04_electron_capture \
  --runs run_E30_s3=$(realpath run_E30_s3) \
         run_base=$(realpath run_base)
```
