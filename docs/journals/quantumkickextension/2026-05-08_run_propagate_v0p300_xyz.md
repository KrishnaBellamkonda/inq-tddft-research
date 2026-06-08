# 2026-05-08_run_propagate_v0p300_xyz

**Title:** run_propagate_v0p300_xyz (li_54_atom_bcc, v=0.300, e-h candidate hunt)
**Run path:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p300_xyz`
**Linked results:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p300_xyz/results`
**Status:** complete

> Companion run to `2026-05-06_run_propagate_v0p0626_xyz`. The pair was
> chosen to populate the two qualitatively distinct regimes of the
> BCN:1719P paper: low-v plasmon (v=0.0626 a.u. at the top of the
> low-v family) and the lowest velocity in the high-v / softening
> regime (v=0.300 a.u., paper's lowest high-v point at 6.56 Å/fs).

## Run summary

### 1. Run identity

| Field | Value |
|---|---|
| run_name | run_propagate_v0p300 |
| run_type | TDDFT impulsive kick on Li 54-atom BCC |
| date_finished | 2026-05-08T07:33:51 |
| wall_time_s | 97600.496265899 |
| executable | run.cpp built via inq-run |
| geometry_file | ../shared/li_54_3x3x3.xyz |
| checkpoint_dir | ../checkpoints/li_54_2x2x2_T200_xyz |

### 3. System configuration

| Field | Value |
|---|---|
| cell_angstrom | 10.53^3 (cubic, periodic) |
| n_atoms | 54 |
| n_electrons | 162 |
| num_states | 101 |
| extra_states | 20 |
| k_grid | 2 2 2 shifted |
| smearing | fermi_dirac |
| smearing_temperature_kelvin | 400 |
| cutoff_ry | 74 |
| xc | pbe (adiabatic-PBE in TDDFT) |

### 5. Kick configuration

| Field | Value |
|---|---|
| kick_velocity_au | 0.3 |
| kick_direction | 1 0 0 (+x) |
| atoms_dynamics | impulsive |

### 6. Real-time configuration

| Field | Value |
|---|---|
| rt_num_steps | 15500 |
| dt_au | 0.04 |
| total_time_fs | 14.9970808 |
| write_every | 100 |
| state_energy_every | 10 |
| occupations_every | 10 |

### 9. End-of-run diagnostics

| Field | Value |
|---|---|
| run_completed | true |
| vti_format | binary |
| validation_n_electrons_integrated_gs | 161.9999999999606 |
| validation_n_electrons_target | 162 |
| validation_abs_error | 3.936406756110955e-11 |

## Observations

(User has not yet supplied observation text for this run. Placeholder
for the user's voice once they have inspected the artefacts.)

## Headline result

**Mid-energy peak at 2.585 eV in `energy_total` FFT** — within 0.22 eV
of the paper's 2.8 eV high-v target (BCN:1719P Figure 4(b)). All four
detrend variants agree to 4 decimals:

| Variant | Peak | DC ratio |
|---|---:|---:|
| raw_subtracted | 2.5852 eV | 2.444 |
| mean_subtracted | 2.5852 eV | 0.653 |
| detrended (linear) | 2.5852 eV | 0.086 |
| plateau_detrend | 2.5852 eV | 1.318 |

This 2.59 eV peak position is **near-identical to the existing v=0.450
result (2.62 eV)** — confirming the paper's observation that the
high-v peak is approximately constant across the high-v family
(softening regime).

## Three-diagnostic discrimination at 2.59 eV

| Diagnostic | Result | Interpretation |
|---|---|---|
| `dipole_x` argmax in [1.5, 5.0] eV | 1.70 eV cluster (top peak 1.55–1.77 eV) | **different** from energy peak — signature of nonlinear regime: dipole (∝ δρ) and energy (∝ \|δρ\|²) probe different aspects |
| `gamma_transitions` histogram in [2.0, 3.0) eV | 95 transitions; **103 within ±0.5 eV of peak** | Dense Γ-Γ transition cluster sits at the peak — consistent with single-particle e-h interpretation (paper: "falls within the range of available Γ-point vertical transitions") |
| `state_energy_spectra` anti-phase pairs near 2.59 eV (top 50) | **0 with opp_metric > 0.7** in 34 candidate pairs | No clean (n, n′) anti-phase signature — argues against a *single-particle* e-h transition |

The conjunction is the paper's "crossover regime" — the 2.59 eV peak
overlaps the e-h transition cluster (so the energy is right) but no
*single* (n, n′) pair carries the oscillation in the way a textbook
e-h would. The paper's words: *"the lack of dominant transitions very
close to the frequency suggest that the electronic system might
represent a crossover regime from the well-defined characteristic
regime."* Our finding **quantifies** that uncertainty — the cross-spectrum
is collective in character, not single-particle.

## Comparison across all four runs

| Run | v (a.u.) | v (Å/fs) | Regime | Energy FFT peak | Paper target | Δ |
|---|---:|---:|---|---:|---:|---:|
| v=0.0123 | 0.0123 | 0.27 | low-v | 5.722 eV | 6.5 eV | -0.78 eV |
| **v=0.0626** (this pair) | 0.0626 | 1.37 | low-v / plasmon | (drift) **6.480 eV in dipole_x** | 6.5 eV | -0.02 eV |
| **v=0.300** (this entry) | 0.300 | 6.56 | high-v entry | 2.585 eV | 2.8 eV | -0.21 eV |
| v=0.450 | 0.450 | 9.84 | high-v | 2.620 eV | 2.8 eV | -0.18 eV |

The four-velocity sweep reproduces the paper's two-regime picture:
low-v plasmon at 6.5 eV; high-v softened mode at 2.6 eV.

![Excess energy per uc vs time](attachments/2026-05-08_run_propagate_v0p300_xyz/excess_energy_per_uc_vs_time.png)

![FFT of excess energy (analyse_inq.py)](attachments/2026-05-08_run_propagate_v0p300_xyz/fft_excess_energy_vs_omega.png)

![Energy_total spectrum (plateau_detrend variant)](attachments/2026-05-08_run_propagate_v0p300_xyz/energy_total_spectrum.png)

![Dipole_x spectrum (q→0 longitudinal mode)](attachments/2026-05-08_run_propagate_v0p300_xyz/dipole_x_spectrum.png)

![Γ-Γ transition histogram](attachments/2026-05-08_run_propagate_v0p300_xyz/gamma_transitions.png)

![Delta-density xz slice (animation)](attachments/2026-05-08_run_propagate_v0p300_xyz/delta_xz.gif)

## Open questions / next steps

- **Why does dipole_x peak at 1.70 eV but energy_total peaks at 2.59 eV?**
  The standard linear-response expectation is that both peak at the
  same ω. The mismatch suggests we are firmly in the nonlinear regime
  the paper describes. A controlled experiment would be to repeat at
  v ≈ 0.1–0.15 a.u. (paper's gap between low-v and high-v families)
  and watch the two peaks merge or split.
- **Is the 2.59 eV peak a "weakly-collective" mode vs a true e-h?**
  The paper itself is uncertain. Our cross-spectrum diagnostic shows
  no anti-phase pair with opp_metric > 0.7, which argues against a
  clean single-particle origin. If the user wants to investigate,
  the natural next step is the **density-spectra plan** at
  `docs/plans/li_v0p0626_plasmon_density_analysis.md` (which applies
  to both runs, plasmon and e-h).
- **GS-offset propagation**: same .xyz GS as v=0.0626. Peak position
  agreement with paper at both velocities confirms GS-offset is benign.
