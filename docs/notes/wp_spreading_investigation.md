# Wavepacket Spreading Investigation

## Physical spreading formula

For a free Gaussian wavepacket (m=1, ħ=1 in atomic units):

    σ(t) = σ₀ √(1 + t²/σ₀⁴)

where t is time in a.u. and σ₀ is the initial Gaussian width in bohr.

This follows from the time-dependent Schrödinger equation for a free particle.
Source: any quantum mechanics textbook, e.g. Griffiths §2.4.

### Timescale for doubling

σ doubles (σ(t) = 2σ₀) when t = σ₀² √3.

| σ₀ | doubling time (a.u.) | doubling time (fs) |
|----|----------------------|--------------------|
| 0.501 bohr (0.265 Å) | 0.435 a.u. | 0.011 fs |
| 1.002 bohr (0.53 Å)  | 1.74 a.u.  | 0.042 fs |
| 3.780 bohr (2.0 Å)   | 24.7 a.u.  | 0.598 fs |

**Key point**: A σ=0.53 Å wavepacket doubles its width in only 0.042 fs. LEED
experiments run orders of magnitude longer than this in classical terms, but
the WP reaches the target molecule at 200 eV before significant spreading.

## LEED paper analysis (Tsubonoya et al. PRB 90, 035416, 2014)

- WP energy: 200 eV → v = √(2·E/m) = √(2 × 200/27.21) = 3.834 bohr/a.u.
- Impact distance: D = 6.35 Å = 12.0 bohr
- Transit time to coronene: t_impact = D/v = 12.0/3.834 = 3.13 a.u. = 0.076 fs
- WP width at impact: σ(3.13) = σ₀ √(1 + 3.13²/σ₀⁴)
  - For σ₀ = 1.002 bohr: σ(3.13) ≈ 1.002 × √(1 + 9.79) ≈ 3.31 bohr (×3.3 spread)

**This is the discrepancy**: naively σ triples before reaching coronene.
However, Tsubonoya et al. saw negligible spreading effects in their LEED patterns.
Possible explanations:
1. The TDDFT propagation is not free (Hartree + XC potentials act on WP).
2. The pattern is recorded time-integrated (spread WP still produces similar pattern).
3. At 200 eV the WP is already very broad relative to lattice spacing.

## To-do

- Use free-propagation run_01_base z-profile to extract σ(t) numerically.
- Compare numerical σ(t) to analytic formula above: verify TDDFT free-particle limit.
- Find (D, E_kin) regime where σ < 2σ₀ on arrival (projectile-like regime).
- Compare coronene LEED patterns at D=3,6,10,15,20 Å to assess practical impact.

## Free propagation expected outcomes

For run_01_base (σ=0.53 Å, 200 eV):
- WP starts at z ≈ Lz − 5σ ≈ 84.9 bohr, moves in −z.
- Hits z=0 boundary at t ≈ Lz/v = 89.9/3.834 = 23.4 a.u. ≈ 0.57 fs.
- At boundary crossing: σ ≈ 1.002 √(1 + 23.4²/1.002⁴) ≈ 23 bohr (×23 spread).
- After boundary reflection: WP propagates back; σ(t) still measurable from z-profile width.
- Expected: energy constant, N_elec = 1.0 ± 0.001, σ(t) matches analytic in early t.
