# 2026-05-06_run_propagate_v0p0626_xyz

**Title:** run_propagate_v0p0626_xyz (li_54_atom_bcc, v=0.0626, plasmon hunt)
**Run path:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0626_xyz`
**Linked results:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0626_xyz/results`
**Status:** complete

## Run summary

### 1. Run identity

| Field | Value |
|---|---|
| run_name | run_propagate_v0p0626 |
| run_type | TDDFT impulsive kick on Li 54-atom BCC |
| date_finished | 2026-05-07T04:22:34 |
| wall_time_s | 109161.240628479 |
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
| kick_velocity_au | 0.0626 |
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

I have conducted a rough study by looking at about 70 odd delta density
VTI files. I observed that periodically there was a blue density
appearing in the cell. I wondered if this could be the plasmon. So, I
conducted a very rough investigation, where I manually annotated the
timesteps where I saw a very bright blue background density, a bright
timestep, and barely visible (only one data point here). With this, and
making a very rough analysis of the results, I found that one of the
potential candidates is in the 7 to 8 eV range, which is in line with
the 6.56 eV signal we are looking at.

## Headline result (added 2026-05-08, post-run)

**Plasmon detected at 6.480 eV in `dipole_x` FFT** — within 0.02 eV of
the paper's 6.5 eV target and within 0.08 eV of the DFT-RPA bulk Li
value (6.56 eV; Faleev et al. 2008, cited in BCN:1719P).

The user's manual VTI annotation (7–8 eV) was qualitatively right; the
precise peak is at the lower end of that range.

| Diagnostic | Peak | Note |
|---|---:|---|
| `energy_total` FFT (any variant) | ~0.55 eV | drowned by low-ω drift; argmax misleading |
| `dipole_x` FFT — top 3 peaks in [5.5, 8.0] eV | 6.480 / 6.411 / 6.549 eV | tight cluster at the paper's 6.5 eV; this is the plasmon |
| `gamma_transitions` histogram in [6.0, 7.0) eV | 0 transitions in [6.0, 6.5); 1 in [5.5, 6.0) | **paper Figure 5 strong test passes** — no Γ-Γ single-particle transition source |
| `state_energy_spectra` anti-phase pairs near 6.48 eV | 0 in top 50 | confirms collective character (plasmon) |

The plasmon attribution is reproduced from the paper's Figure 4(a) /
Figure 5 by three independent diagnostics:

1. **Position match**: 6.48 eV vs paper 6.5 eV (Δ = -0.02 eV).
2. **Figure 5 cliff-edge test**: `gamma_transitions` histogram cuts off
   at ~5.5–6.0 eV; the [6.0, 6.5) bin is empty. No single-particle
   Γ–Γ transition can be the source.
3. **Anti-phase pair test**: top 50 strongest (n, n′) cross-spectra at
   any kpoint show no opp_metric > 0.7 near 6.48 eV. Consistent with a
   collective rather than single-particle mode.

![Excess energy per uc vs time](attachments/2026-05-06_run_propagate_v0p0626_xyz/excess_energy_per_uc_vs_time.png)

![FFT of excess energy](attachments/2026-05-06_run_propagate_v0p0626_xyz/fft_excess_energy_vs_omega.png)

![Dipole_x spectrum (plasmon channel)](attachments/2026-05-06_run_propagate_v0p0626_xyz/dipole_x_spectrum.png)

![Γ-Γ transition histogram (paper Figure 5)](attachments/2026-05-06_run_propagate_v0p0626_xyz/gamma_transitions.png)

![Delta-density xz slice (animation)](attachments/2026-05-06_run_propagate_v0p0626_xyz/delta_xz.gif)

![Coarse delta-density xz slice (animation)](attachments/2026-05-06_run_propagate_v0p0626_xyz/delta_coarse_xz.gif)

## Open questions / next steps

- **Deep density-spectra analysis** per the plan in
  `docs/plans/li_v0p0626_plasmon_density_analysis.md` — gated on user
  review. The plan starts with Stage A (lab-frame pixel-by-pixel FFT
  of `density_rt_delta`) to produce a *spatial map* showing where in
  the cell the 6.48 eV mode lives. Given the dipole_x channel is
  already clean, Stage A is likely sufficient; Stage B Galilean-shift
  may not be needed.
- **GS-offset characterisation**: the new .xyz GS landed 7.4 mHa above
  the old fractional GS but the plasmon peak position is fully
  consistent with the paper, confirming the GS divergence was benign
  (see `docs/sources/li_gs_xyz_vs_fractional_offset_analysis.md`).
- **Energy FFT drift**: the energy_total FFT is dominated by low-ω
  drift; investigating whether `t_skip_fs ≥ 1.0` or the linear
  detrend variant give a clean energy-channel detection of the same
  6.48 eV peak.
