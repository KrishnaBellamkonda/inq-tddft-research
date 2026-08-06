# Effective width per run — sigma56_sv wavepackets

`sigma_r(t) = sqrt(sigma_x^2+sigma_y^2+sigma_z^2)`; all widths in Bohr.

**⟨σ_r⟩ 1%-window** is the adopted label: the mean of σ_r(t) from
t = 0 until the packet norm has fallen by 1 % (user decision,
2026-08-03). `t(1%)` is where that window closes; compare it with
`t_exit`, the in-slab transit end — the window should cover the
crossing. `⟨σ_r⟩ full` is the superseded full-run mean, which keeps
averaging long after the CAP has destroyed the packet (see
`norm(t_f)`) and so measures a smeared remnant, not the projectile.

| run | v | steps | t_f | σ_r(0) | σ_r@exit | σ_r(t_f) | **⟨σ_r⟩ 1%-window** | t(1%) | t_exit | ⟨σ_r⟩ full | norm(t_f) |
|---|---|---|---|---|---|---|---|---|---|---|
| `s5p0_v2p0` | 2.0 | 1973/4360 | 78.9 | 6.12 | 8.64 | 25.53 | **7.76** | 28.2 | 20.0 | 12.27 | 1.8e-05 |
| `s5p0_v2p5` | 2.5 | 3488/3488 | 139.5 | 6.12 | 7.88 | 29.42 | **7.34** | 23.2 | 16.0 | 21.31 | 1.5e-08 |
| `s5p0_v3p0` | 3.0 | 2907/2907 | 116.3 | 6.12 | 7.41 | 25.17 | **7.05** | 19.8 | 13.3 | 18.54 | 5.3e-10 |
| `s5p0_v3p5` | 3.5 | 2491/2491 | 99.6 | 6.12 | 7.1 | 27.09 | **6.85** | 17.2 | 11.4 | 18.42 | 3.7e-10 |
| `s6p0_v2p0` | 2.0 | 4360/4360 | 174.4 | 7.35 | 9.03 | 30.29 | **8.45** | 28.2 | 20.0 | 22.56 | 8.1e-07 |
| `s6p0_v2p5` | 2.5 | 3488/3488 | 139.5 | 7.35 | 8.52 | 26.86 | **8.15** | 23.1 | 16.0 | 19.71 | 1.1e-08 |
| `s6p0_v3p0` | 3.0 | 2907/2907 | 116.3 | 7.35 | 8.2 | 23.6 | **7.95** | 19.6 | 13.3 | 18.33 | 6.1e-10 |
| `s6p0_v3p5` | 3.5 | 2491/2491 | 99.6 | 7.35 | 7.99 | 26.28 | **7.81** | 17.0 | 11.4 | 18.36 | 3.9e-10 |

## Per-sigma means (what a legend label would read)

| σ_WP | **⟨σ_r⟩ 1%-window** | ⟨σ_r⟩ full (superseded) | n runs |
|---|---|---|---|
| 5 | **7.25** | 17.63 | 4 |
| 6 | **8.09** | 19.74 | 4 |
