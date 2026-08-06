# lz_bulk_sweep — PILOT REPORT (v = 3.0, all four boxes, both halves)

Generated 2026-08-05 15:53 UTC (SLURM job n/a)

## VERDICT: PASS — production released (report-only mode, nothing gated)

### Notes (non-gating)

- s0p5_L15_v3p0: plateau not settled (drift -0.140 eV over the last 10 %)
- cl_s0p5_L15_v3p0: plateau not settled (drift -0.916 eV over the last 10 %)
- s0p5_L35_v3p0: plateau not settled (drift -0.281 eV over the last 10 %)
- cl_s0p5_L35_v3p0: plateau not settled (drift -1.675 eV over the last 10 %)
- cl_s5p0_L15_v3p0: plateau not settled (drift -0.708 eV over the last 10 %)
- cl_s5p0_L35_v3p0: plateau not settled (drift -1.313 eV over the last 10 %)

## Pilot S values (corrected deposit, eV/Bohr)

| run | S_deposit | S (no E_PS cut) | norm_final | settled | steps |
|---|---|---|---|---|---|
| `cl_s0p5_L15_v3p0` | 0.373 | 0.863 | nan | False | 2053/2053 |
| `s0p5_L15_v3p0` | 0.189 | 0.189 | 3.09e-02 | False | 2053/2053 |
| `cl_s0p5_L35_v3p0` | 0.387 | 0.769 | nan | False | 2780/2780 |
| `s0p5_L35_v3p0` | 0.160 | 0.160 | 3.01e-02 | False | 2780/2780 |
| `cl_s5p0_L15_v3p0` | 0.131 | 0.519 | nan | False | 2543/2543 |
| `s5p0_L15_v3p0` | 0.309 | 0.309 | 6.41e-10 | True | 2543/2543 |
| `cl_s5p0_L35_v3p0` | 0.133 | 0.441 | nan | False | 3270/3270 |
| `s5p0_L35_v3p0` | 0.445 | 0.445 | 3.83e-10 | True | 3270/3270 |

## S(L) ordering vs the L = 25 anchors (INFO)

- sigma=0.5 wp: S(15)=0.189  S(25,anchor)=0.167  S(35)=0.160 eV/Bohr  [monotone in L]
- sigma=0.5 classical: S(15)=0.373  S(25,anchor)=0.381  S(35)=0.387 eV/Bohr  [monotone in L]
- sigma=5.0 wp: S(15)=0.309  S(25,anchor)=0.396  S(35)=0.445 eV/Bohr  [monotone in L]
- sigma=5.0 classical: S(15)=0.131  S(25,anchor)=0.133  S(35)=0.133 eV/Bohr  [monotone in L]

## GS interior bulk-likeness (WARN-only)

- s0p5_L15: n(z=0)/n0 - 1 = -3.22 % [WARN], interior peak-to-peak 22.9 % of n0
- s0p5_L35: n(z=0)/n0 - 1 = +0.84 % [ok], interior peak-to-peak 17.7 % of n0
- s5p0_L15: n(z=0)/n0 - 1 = -3.64 % [WARN], interior peak-to-peak 22.6 % of n0
- s5p0_L35: n(z=0)/n0 - 1 = +0.76 % [ok], interior peak-to-peak 20.4 % of n0

## Production cost projection (WARN, never a gate)

- s0p5_L15/wp: 1.75 s/step (measured) -> 7302 production steps = 3.5 h
- s0p5_L15/classical: 1.30 s/step (estimate) -> 7302 production steps = 2.6 h
- s0p5_L35/wp: 4.46 s/step (measured) -> 9886 production steps = 12.3 h
- s0p5_L35/classical: 4.00 s/step (estimate) -> 9886 production steps = 11.0 h
- s5p0_L15/wp: 1.98 s/step (measured) -> 9047 production steps = 5.0 h
- s5p0_L15/classical: 1.70 s/step (estimate) -> 9047 production steps = 4.3 h
- s5p0_L35/wp: 4.89 s/step (measured) -> 11632 production steps = 15.8 h
- s5p0_L35/classical: 4.80 s/step (estimate) -> 11632 production steps = 15.5 h
- TOTAL projected production: 70 GPU-h (+ vacuum baselines; proceeding is the default — kill with scancel if unwanted)

Plan: `docs/plans/jellium-slab-extend-Lz.md` · Handover: `docs/handovers/jellium-slab-extend-Lz.md`

