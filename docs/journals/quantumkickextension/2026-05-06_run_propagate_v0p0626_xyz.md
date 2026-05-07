# 2026-05-06_run_propagate_v0p0626_xyz

**Title:** run_propagate_v0p0626_xyz (li_54_atom_bcc, v=0.0626, plasmon hunt)
**Run path:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0626_xyz`
**Linked results:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0626_xyz/results`
**Status:** running

> ⚠️ Entry created while the run is in progress. `run_summary.txt` is
> written by the run.cpp template only at the end of the propagation,
> so the canonical config table is not yet available. The
> "Configured parameters" section below mirrors the source of truth in
> `shared/configs/base_li_54.hpp` (struct `Li_54_v0p0626`). When the run
> completes (~25 h from launch on 2026-05-05 22:03 BST), this entry will
> be updated to (a) replace the configured-parameters table with the
> verbatim `run_summary.txt` two-column table, (b) flip the status, and
> (c) append the post-run delta-density investigation requested by the
> user.

## Configured parameters (placeholder — to be replaced by `run_summary.txt` post-run)

### 1. Run identity

| Field | Value |
|---|---|
| run_name | run_propagate_v0p0626_xyz |
| run_type | TDDFT impulsive kick on Li 54-atom BCC |
| date_started | 2026-05-05 22:03 BST |
| executable | run.cpp built via inq-run |
| geometry_file | ../shared/li_54_3x3x3.xyz |
| checkpoint_dir | ../checkpoints/li_54_2x2x2_T200_xyz |

### 3. System configuration

| Field | Value |
|---|---|
| cell_angstrom | 10.53^3 (cubic, periodic) |
| n_atoms_expected | 54 |
| n_electrons | 162 |
| num_states | 101 (per kpoint × 8 → 808 total) |
| extra_states | 20 |
| k_grid | 2 2 2 shifted MP |
| smearing | fermi_dirac |
| smearing_temperature_kelvin | 400 |
| cutoff_ry | 74 |
| xc | pbe (adiabatic in TDDFT) |

### 5. Kick configuration

| Field | Value |
|---|---|
| kick_velocity_au | 0.0626 |
| kick_direction | 1 0 0 (+x) |
| atoms_dynamics | impulsive |
| qball_velocity_class | low-v family (0.0123–0.0626 a.u.) |
| target_peak_ev | ~6.5 (DFT-RPA bulk Li plasmon at 6.56 eV) |

### 6. Real-time configuration

| Field | Value |
|---|---|
| rt_num_steps | 15500 |
| dt_au | 0.04 |
| total_time_fs | ~15.0 |
| write_every | 100 (density VTI cadence) |
| state_energy_every | 10 |
| occupations_every | 10 |

### 9. End-of-run diagnostics (pending)

| Field | Value |
|---|---|
| run_completed | (pending) |
| wall_time_s | (pending) |
| validation_n_electrons_integrated_gs | 162.000 (verified at GS load) |
| ground_state_offset_vs_old_fractional_ha | +0.0074 (see source note `li_gs_xyz_vs_fractional_offset_analysis.md`) |

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

## Open questions / next steps

- Conduct an in-depth post-run study: examine the delta-density signal
  proper. A very strong delta-density signal would be due to the ions
  themselves moving (rigid translation of the ionic potential drags the
  electron cloud with it). The plasmon signal sits *on top* of that.
- Read the **unfiltered** `density_rt_total/*.vti` frames; compute the
  density-over-time at every grid point.
- Apply a temporal Fourier transform of the (delta) density at each
  grid point. The plasmon should appear as a peak around 6.5 eV in the
  averaged-over-cell spectrum.
- If the temporal FFT of the raw delta density does NOT show a clear
  6.5 eV peak (because the ion-translation component dominates),
  investigate ways to subtract / weave out the rigid ion-motion
  component before the FFT (e.g. Galilean shift to the comoving frame,
  per-pixel polynomial detrend, or projection onto specific Fourier
  modes).
- Once a clean approach is identified, write a well-thought-out,
  scientific plan to identify the plasmon mode in this run. **The plan
  must be reviewed by the user before execution.**
- Cross-check the plasmon attribution against:
  - the existing `dipole_x` FFT (the q→0 longitudinal-density proxy —
    see `docs/sources/dipole_as_q0_density_projection.md`);
  - the `gamma_transitions` histogram to confirm the absence of any
    Γ–Γ vertical transition near 6.5 eV (paper Figure 5);
  - the new `state_energy_spectra` per-state εN(t) FFT and anti-phase
    pair diagnostic to confirm no single-particle (n, n′) pair carries
    the 6.5 eV oscillation.
