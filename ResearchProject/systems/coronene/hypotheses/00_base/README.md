# 00 — Base run diagnostic summary

Reference point for every other hypothesis. Single run.

## Hypothesis

The Tsubonoya 2014 paper geometry, replicated with the corrected
coronene.xyz (atoms at z = 0, INQ centred-cell convention `[-L/2, +L/2]`),
yields a 6-fold-symmetric LEED pattern at the paper window. This run is
the comparison reference for all other hypotheses; no separate hypothesis
is being tested here.

## Runs collated

- `run_base` — b = 12 Bohr, σ = 1.0 Bohr, E = 200 eV, cell 35 × 35 × 60 Bohr.

## Produce comparison artefacts

```
python scripts/coronene_postprocess.py hypothesis \
  --hypothesis-dir hypotheses/00_base \
  --runs run_base=$(realpath run_base)
```

This collation is mainly a sanity check — most useful artefacts for the
base run live under `run_base/results/analysis/`.

## Source

Tsubonoya, Hu & Watanabe, *Phys. Rev. B* **90**, 035416 (2014). See
`docs/sources/tsubonoya-2014-coronene-leed.md`.
