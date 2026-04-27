# 05 — Box length and impact distance

## Hypothesis

The simulated LEED pattern should depend on:

1. **Impact distance b** — at fixed cell size, larger b ⇒ "far-field"
   regime where the angular pattern is well-defined; smaller b ⇒ "near
   field" with WP density penetrating the molecular charge cloud.
2. **Box length L_z** — too small a box truncates the WP before it has
   propagated past the molecule; too large a box uses GPU time without
   added physics. There should be a minimum L_z at which the LEED pattern
   stops shifting with further box growth.

Three runs probe both axes: a larger-box pair scans b at fixed cell, and
a smaller-box run shows the truncation effect.

## Runs collated

- `run_b18_35x35x80` — large box (L_z = 80 Bohr), b = 18 Bohr (= 1.5 × base).
- `run_b6_35x35x80`  — large box (L_z = 80 Bohr), b = 6 Bohr (= 0.5 × base).
- `run_35x35x40`     — small box (L_z = 40 Bohr = 2/3 × base), b at base (12 Bohr).
- `run_base`         — paper reference (L_z = 60 Bohr, b = 12 Bohr).

## Produce comparison artefacts

```
python scripts/coronene_postprocess.py hypothesis \
  --hypothesis-dir hypotheses/05_box_length_and_distance \
  --runs run_b18_35x35x80=$(realpath run_b18_35x35x80) \
         run_b6_35x35x80=$(realpath run_b6_35x35x80) \
         run_35x35x40=$(realpath run_35x35x40) \
         run_base=$(realpath run_base)
```
