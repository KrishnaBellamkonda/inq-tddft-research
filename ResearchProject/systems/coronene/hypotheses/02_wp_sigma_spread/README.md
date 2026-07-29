# 02 — Wave-packet σ spread (narrow / medium / wide)

## Hypothesis

Wider real-space σ ⇒ narrower momentum spread ⇒ better-resolved LEED
fringes (more "monochromatic"); narrower σ ⇒ broader momentum spread ⇒
washed-out angular features.

## Runs collated

- `run_s0p33` — σ = 0.33 Bohr (= base/3, narrow).
- `run_base`  — σ = 1.0 Bohr (medium / paper reference).
- `run_s3`    — σ = 3.0 Bohr (= 3 × base, wide).

All other parameters (b, E, cell) are at the Tsubonoya base.

## Produce comparison artefacts

```
python scripts/coronene_postprocess.py hypothesis \
  --hypothesis-dir hypotheses/02_wp_sigma_spread \
  --runs run_s0p33=$(realpath run_s0p33) \
         run_base=$(realpath run_base) \
         run_s3=$(realpath run_s3)
```

Outputs as in `01_wp_energy_spread`.
