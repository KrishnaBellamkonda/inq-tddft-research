# Results Folder Structure Specification for INQ / inqkit / inqview Runs

This document defines the preferred `results/` directory structure for INQ / inqkit / inqview wave-packet real-time TDDFT, LEED, jellium, quantum-kick, and related runs.

The goal is **categorisation**, not rigid templating. Most of the output infrastructure already exists in previous runs. The task is therefore to place existing raw files, analysis products, VTI files, screen data, overlap data, and summaries into a consistent structure that is easy to browse, easy to copy, and easy to process with Python and ParaView.

Only one file is treated as a required textual template:

```text
results/run_summary.txt
```

All other files should follow the naming/category conventions below, but this document does **not** prescribe exact column templates or detailed file internals unless needed for categorisation.

---

# 1. Core principles

## 1.1 Top-level structure must stay simple

The top-level `results/` directory should contain only:

```text
results/
├── run_summary.txt
├── raw/
└── analysis/
```

Do **not** create extra top-level files such as:

```text
manifest.txt
simulation_config.txt
derived_values.txt
command_line.txt
git_info.txt
timing_summary.txt
```

All of that information belongs inside `run_summary.txt`.

---

## 1.2 Use `raw/` for numerical data

Use:

```text
results/raw/
```

for numerical data that is written by the simulation or produced by deterministic post-processing.

Examples:

```text
.raw
.meta.txt
.csv
.dat
.vti
```

This includes:

- raw density fields
- raw orbital fields
- raw screen/LEED arrays
- overlap matrices
- time-domain observable CSVs
- FFT/spectrum CSVs
- VTI files used by ParaView

---

## 1.3 Use `analysis/` for human-facing outputs

Use:

```text
results/analysis/
```

for plots, GIFs, movies, rendered frames, contact sheets, filtered images, and visual summaries.

Examples:

```text
.png
.gif
.mp4
```

Do **not** use:

```text
processed/
derived_data/
```

The word `analysis/` is the standard name for visual and post-processed outputs.

---

## 1.4 Raw CSVs go in `raw/observables/`

All numerical CSV files associated with observables or observable-derived quantities should go under:

```text
results/raw/observables/
```

This includes both direct simulation outputs and post-processing numerical outputs such as FFT spectra.

Examples:

```text
observables.csv
energy_components.csv
current_components.csv
dipole_components.csv
fft_total_energy.csv
fft_current_x.csv
fft_current_y.csv
fft_current_z.csv
dipole_spectrum_x.csv
wp_center_of_mass.csv
wp_width.csv
charge_conservation.csv
norm_conservation.csv
```

Do not create a separate `derived_data/` folder for FFTs or other derived numerical arrays.

---

## 1.5 Observable plots go in `analysis/observables/`

All plots made from `raw/observables/` should go under:

```text
results/analysis/observables/
```

Examples:

```text
total_energy_vs_time.png
all_energy_components_vs_time.png
current_components_vs_time.png
fft_total_energy.png
fft_current_components.png
dipole_spectrum_components.png
```

Keep this folder mostly flat. The filename should make the plot content clear.

---

## 1.6 Time-stepped GIFs must use a fixed colour scale

Any GIF, movie, or frame sequence made from time-stepped data must use a **single colour scale for the entire animation**.

This applies to:

- density slice GIFs
- density volume GIFs
- wavepacket density GIFs
- instantaneous LEED screen GIFs
- overlap heatmap GIFs
- wavepacket-overlap GIFs
- orbital-density GIFs

The colour scale must not auto-rescale independently at each timestep, because that makes the time evolution visually misleading.

Recommended rule:

```text
For a time-series visualisation, compute vmin/vmax or log-scale limits from the entire series first, then render every frame with those fixed limits.
```

For strongly peaked data, use a fixed log scale or fixed percentile-based limits, but the limits must still be fixed through time.

---

# 2. Canonical directory tree

```text
results/
├── run_summary.txt
│
├── raw/
│   ├── ground_state/
│   ├── wavepacket/
│   ├── density/
│   ├── orbitals/
│   ├── screens/
│   ├── overlap/
│   ├── observables/
│   └── vti/
│
└── analysis/
    ├── ground_state/
    ├── observables/
    ├── density/
    ├── screens/
    ├── overlap/
    └── orbitals/
```

---

# 3. `run_summary.txt`

## 3.1 Purpose

`run_summary.txt` is the only top-level summary file. It should contain enough information to understand, reproduce, and diagnose the run without needing several separate metadata files.

It should include:

- run start and end time
- wall time
- host and working directory
- command used to launch the run
- code/version information
- system configuration
- ground-state configuration and outcome
- wavepacket configuration and injection diagnostics
- real-time configuration
- screen configuration
- output map
- end-of-run diagnostics
- warnings or failure status

---

## 3.2 Required template

Use this structure for `results/run_summary.txt`:

```text
RUN SUMMARY
===========

1. Run identity
---------------
run_name:
run_type:
date_started:
date_finished:
wall_time:
host:
working_directory:
executable:
command:
gpu_id:
notes:

2. Code and environment
-----------------------
inq_version:
inqkit_version:
inqview_version:
compiler:
cuda_version:
git_branch:
git_commit:
git_dirty:
python_version:
paraview_version:

3. System configuration
-----------------------
cell_type:
cell_lengths_bohr:
boundary_conditions:
geometry_file:
n_ions:
n_electrons:
n_occupied_states:
extra_states:
wp_state_index:
kpoints:
cutoff_ha:
grid_shape:
grid_spacing_bohr:
temperature:
xc_functional:
pseudopotentials:

4. Ground-state configuration and outcome
----------------------------------------
gs_algorithm:
energy_tolerance:
max_steps:
mixing_type:
mixing_alpha:
converged:
n_scf_iterations:
final_total_energy:
final_kinetic_energy:
final_hartree_energy:
final_xc_energy:
final_external_energy:
final_nonlocal_energy:
final_ion_ion_energy:
final_energy_residual:
final_density_residual:
homo_index:
lumo_index:
homo_energy:
lumo_energy:
homo_lumo_gap:
fermi_energy:
total_charge_integral:
normalisation_check:

5. Wavepacket configuration and injection
-----------------------------------------
wp_enabled:
wp_center_bohr:
wp_sigma_bohr:
wp_k0_bohr_inv:
wp_energy_ev:
wp_direction:
wp_occupation:
wp_state_index:
orthogonalised:
orthogonalisation_tolerance:
norm_before:
norm_after:
max_overlap_before:
max_overlap_after:
passed_tolerance:

6. Real-time configuration
--------------------------
rt_num_steps:
dt_au:
total_time_au:
write_every_density:
write_every_screens:
write_every_overlap:
propagator:
observables_enabled:
density_outputs_enabled:
screen_outputs_enabled:
overlap_outputs_enabled:

7. Screen configuration
-----------------------
n_screens:
screen_orientation:
screen_positions_bohr:
instantaneous_screen_cadence:
time_windowed_screen_windows:
total_screen_accumulation_range:

8. Output map
-------------
raw_ground_state:
raw_wavepacket:
raw_density:
raw_screens:
raw_overlap:
raw_observables:
raw_vti:
analysis_ground_state:
analysis_observables:
analysis_density:
analysis_screens:
analysis_overlap:
analysis_orbitals:

9. End-of-run diagnostics
-------------------------
run_completed:
error_status:
last_completed_step:
final_time_au:
final_total_energy:
energy_drift:
charge_drift:
norm_drift:
max_density:
min_density:
max_wp_density:
final_wp_center:
final_wp_width:
warnings:
```

If the run fails early, still write `run_summary.txt`. In that case:

```text
run_completed: false
last_completed_step: <last completed step>
warnings: <reason or exception if available>
```

---

# 4. `results/raw/`

`raw/` stores all numerical data. It should be organised by scientific category, not by implementation detail.

```text
results/raw/
├── ground_state/
├── wavepacket/
├── density/
├── orbitals/
├── screens/
├── overlap/
├── observables/
└── vti/
```

---

# 5. `results/raw/ground_state/`

## Purpose

Ground-state outputs belong here. This includes data needed to judge whether the SCF calculation converged and whether the initial electronic state is sensible.

## Category layout

```text
results/raw/ground_state/
├── summary.txt
├── scf_history.csv
├── eigenvalues.csv
├── occupations.csv
├── ground_state_observables.csv
├── density_system/
├── density_gs_orbitals/
└── checkpoint/
```

## Notes

- `summary.txt` is ground-state specific. It is not a replacement for `run_summary.txt`.
- `density_system/` stores the total ground-state density before wavepacket injection.
- `density_gs_orbitals/` stores individual GS orbital densities.
- `checkpoint/` stores INQ restart/checkpoint data.
- The exact CSV columns can follow the existing run infrastructure. This document does not require a new CSV template.

---

# 6. `results/raw/wavepacket/`

## Purpose

Wavepacket setup and injection diagnostics belong here.

## Category layout

```text
results/raw/wavepacket/
├── wavepacket_config.txt
├── injection_report.txt
├── orthogonality_report.csv
├── density_wp_initial/
└── wavefunction_wp_initial/
```

## Notes

- Use this folder for the initial packet parameters and injection diagnostics.
- The initial WP density and initial complex WP wavefunction should be separated because they are different objects.
- If the run has no injected wavepacket, this folder may be absent.

---

# 7. `results/raw/density/`

## Purpose

Real-time density fields belong here.

There must be three separate density categories whenever a WP is present:

```text
density_rt_system/
density_rt_wp/
density_rt_total/
```

where:

```text
density_rt_total = density_rt_system + density_rt_wp
```

The distinction is important:

- `density_rt_system/` tracks the original occupied system density.
- `density_rt_wp/` tracks the injected wavepacket density.
- `density_rt_total/` is the full visible density.

## Category layout

```text
results/raw/density/
├── density_rt_total/
├── density_rt_system/
└── density_rt_wp/
```

## Notes

- Large real-space 3D fields should continue to use the existing `.raw + .meta.txt` convention.
- Do not remove `density_rt_total/`; it is required for full-density visualisation.

---

# 8. `results/raw/orbitals/`

## Purpose

Complex real-time orbital wavefunctions belong here.

This is separate from orbital densities and VTI orbital visualisation files.

## Category layout

```text
results/raw/orbitals/
├── homo/
├── lumo/
├── wp/
└── selected_orbitals/
```

## Notes

- Keep complex orbital raw files here.
- Orbital-density VTI files used for ParaView belong in `raw/vti/`, not here.
- Existing complex-field conventions can be reused.

---

# 9. `results/raw/screens/`

## Purpose

Raw LEED/screen detector arrays belong here.

There are three screen categories:

```text
total/
instantaneous/
time_windowed/
```

## Category layout

```text
results/raw/screens/
├── screen_config.csv
├── window_ranges.csv
├── total/
├── instantaneous/
└── time_windowed/
```

---

## 9.1 Total screens

Use:

```text
results/raw/screens/total/
```

for final accumulated or time-averaged screen patterns.

Example:

```text
screen_00.dat
screen_01.dat
screen_02.dat
```

---

## 9.2 Instantaneous screens

Use:

```text
results/raw/screens/instantaneous/
```

for single-timestep screen snapshots.

This folder must be flat. Do not create timestep subfolders.

Use filenames like:

```text
screen_00_t000000.dat
screen_00_t000003.dat
screen_01_t000003.dat
```

---

## 9.3 Time-windowed screens

Use:

```text
results/raw/screens/time_windowed/
```

for screen patterns accumulated over timestep windows.

This folder must be flat. Do not create `window_00/`, `window_01/`, or per-screen subfolders.

Use filenames like:

```text
screen_00_t000000_to_t000300.dat
screen_01_t000000_to_t000300.dat
screen_00_t000300_to_t000600.dat
```

The filename must encode:

```text
screen index + start timestep + end timestep
```

---

## 9.4 Screen coordinate mapping requirement

The Python LEED plotting code must use the correct detector-plane coordinate mapping.

A known failure mode is that the diffraction/screen pattern appears split into four pieces and placed at the plot edges. This usually indicates that the screen array is being interpreted with the wrong coordinate origin, array ordering, transpose convention, or shift convention.

Therefore:

- `screen_config.csv` and/or existing `.dat` headers must provide enough information to map array indices to detector-plane coordinates.
- The plotting code must not blindly assume that array index `(0, 0)` corresponds to the physical lower-left corner unless the raw writer actually uses that convention.
- The plotting code must be explicit about whether the screen axes are centred around zero, cell-centred, or stored in raw index order.
- Do not apply `fftshift`, `ifftshift`, `roll`, quadrant swapping, or transpose operations to LEED images unless the raw data convention requires it and the transformation is documented in the analysis code.
- A direct raw-index plot and a coordinate-mapped plot should be easy to compare during debugging.

Recommended analysis diagnostic outputs:

```text
results/analysis/screens/coordinate_checks/
├── screen_00_raw_index_plot.png
├── screen_00_coordinate_mapped_plot.png
├── screen_00_axis_orientation_check.png
└── ...
```

If the pattern is split across the four plot edges, fix the coordinate mapping in the loader/plotter rather than treating it as a physical LEED feature.

---

# 10. `results/raw/overlap/`

## Purpose

Overlap data belongs here. The overlap matrices track how real-time orbitals project onto the ground-state KS orbital basis.

This is especially important for tracking the injected wavepacket.

## Category layout

```text
results/raw/overlap/
├── index.csv
├── overlap_000000.csv
├── overlap_000001.csv
├── overlap_000002.csv
└── ...
```

The folder should be flat.

## Required wavepacket-overlap coverage

The overlap of the wavepacket state with **all ground-state KS orbitals at all saved timesteps** must be available from this folder.

The preferred implementation is:

- keep the full overlap matrix snapshots as already planned; and
- ensure the WP column/row corresponding to the injected state is included at every timestep.

If a convenience CSV is produced for easier plotting, place it here too:

```text
results/raw/overlap/wp_overlap_with_gs_orbitals.csv
```

This convenience file is optional if the same information is fully recoverable from the matrix snapshots and `index.csv`, but the analysis code must be able to generate the WP-overlap GIF either way.

---

# 11. `results/raw/observables/`

## Purpose

All tabular observables and observable-derived numerical arrays belong here.

## Category layout

```text
results/raw/observables/
├── observables.csv
├── energy_components.csv
├── current_components.csv
├── dipole_components.csv
├── norm_conservation.csv
├── charge_conservation.csv
├── wp_center_of_mass.csv
├── wp_width.csv
├── fft_total_energy.csv
├── fft_current_x.csv
├── fft_current_y.csv
├── fft_current_z.csv
├── dipole_spectrum_x.csv
├── dipole_spectrum_y.csv
└── dipole_spectrum_z.csv
```

## Notes

- `observables.csv` is the main time-domain observable file.
- Energy, current, and dipole component CSVs may be direct outputs or analysis-ready extracted files.
- FFT and spectrum CSVs belong here, not in `analysis/` and not in `derived_data/`.
- Exact column details can follow the existing infrastructure.

---

# 12. `results/raw/vti/`

## Purpose

VTI files used by ParaView belong here.

Although VTI files are derived from raw field arrays, they are data products that are repeatedly reused for visualisation. Therefore they live under:

```text
results/raw/vti/
```

not under `analysis/`.

---

## 12.1 ParaView-friendly flat VTI rule

ParaView automatically detects and selects file series when related `.vti` files are placed in the same folder with sortable filenames.

Therefore, VTI folders must be **flat by series**.

Do not nest each timestep or orbital inside separate subfolders if the files are intended to be selected together in ParaView.

---

## 12.2 Category layout

```text
results/raw/vti/
├── density_gs_system/
├── density_gs_orbitals/
├── density_rt_total/
├── density_rt_system/
├── density_rt_wp/
└── orbitals/
```

---

## 12.3 Real-time density VTI files

All timesteps for a density series must live directly inside the relevant folder.

Use:

```text
results/raw/vti/density_rt_total/
├── density_rt_total_t000000.vti
├── density_rt_total_t000010.vti
├── density_rt_total_t000020.vti
└── ...

results/raw/vti/density_rt_system/
├── density_rt_system_t000000.vti
├── density_rt_system_t000010.vti
├── density_rt_system_t000020.vti
└── ...

results/raw/vti/density_rt_wp/
├── density_rt_wp_t000000.vti
├── density_rt_wp_t000010.vti
├── density_rt_wp_t000020.vti
└── ...
```

Do not use:

```text
results/raw/vti/density_rt_wp/t000000/density.vti
results/raw/vti/density_rt_wp/t000010/density.vti
```

Flat files make it much easier to select the full time series in ParaView.

---

## 12.4 Ground-state orbital VTI files

Ground-state orbital density VTI files must be directly inside:

```text
results/raw/vti/density_gs_orbitals/
```

Use:

```text
results/raw/vti/density_gs_orbitals/
├── orbital_0000.vti
├── orbital_0001.vti
├── orbital_0002.vti
├── orbital_0003.vti
└── ...
```

Do not use:

```text
results/raw/vti/density_gs_orbitals/orbital_0000/density_t000000.vti
results/raw/vti/density_gs_orbitals/orbital_0001/density_t000000.vti
```

The flat form is easier to copy, inspect, and select in ParaView.

---

## 12.5 Real-time orbital VTI files

Real-time orbital VTI files must also be flat enough for ParaView selection.

Use filenames that encode both the orbital index and timestep:

```text
results/raw/vti/orbitals/
├── orbital_0040_t000000.vti
├── orbital_0040_t000010.vti
├── orbital_0040_t000020.vti
├── orbital_0041_t000000.vti
├── orbital_0041_t000010.vti
└── ...
```

If only one orbital is exported over time, the folder can contain just that orbital series.

If many orbitals are exported and ParaView grouping becomes inconvenient, a later extension may split by series, but the first implementation should prefer flat files because copying and manual visualisation are easier.

---

# 13. `results/analysis/`

`analysis/` stores human-facing analysis products.

```text
results/analysis/
├── ground_state/
├── observables/
├── density/
├── screens/
├── overlap/
└── orbitals/
```

---

# 14. `results/analysis/ground_state/`

## Purpose

Ground-state plots and visual summaries belong here.

## Category layout

```text
results/analysis/ground_state/
├── scf_convergence.png
├── eigenvalue_spectrum.png
├── occupations.png
├── density_gs_system.png
├── gs_orbital_gallery.png
└── ground_state_summary.png
```

## Notes

These plots should help answer:

- Did the SCF converge?
- Are the final energies stable?
- Are the occupations sensible?
- Are the HOMO/LUMO and orbital densities plausible?
- Is the ground-state density physically sensible?

---

# 15. `results/analysis/observables/`

## Purpose

Plots made from `raw/observables/` belong here.

## Category layout

```text
results/analysis/observables/
├── observables_summary.png
├── total_energy_vs_time.png
├── kinetic_energy_vs_time.png
├── hartree_energy_vs_time.png
├── xc_energy_vs_time.png
├── all_energy_components_vs_time.png
├── excess_energy_vs_time.png
├── current_x_vs_time.png
├── current_y_vs_time.png
├── current_z_vs_time.png
├── current_components_vs_time.png
├── current_magnitude_vs_time.png
├── dipole_x_vs_time.png
├── dipole_y_vs_time.png
├── dipole_z_vs_time.png
├── dipole_components_vs_time.png
├── fft_total_energy.png
├── fft_energy_components.png
├── fft_current_x.png
├── fft_current_y.png
├── fft_current_z.png
├── fft_current_components.png
├── dipole_spectrum_x.png
├── dipole_spectrum_y.png
├── dipole_spectrum_z.png
└── dipole_spectrum_components.png
```

## Notes

- Keep observable plots in this folder, not in separate nested `energy/`, `current/`, or `spectra/` folders unless the folder becomes too large later.
- Use filenames that clearly describe the plotted quantity.

---

# 16. `results/analysis/density/`

## Purpose

Density animations and rendered density visualisations belong here.

## Category layout

```text
results/analysis/density/
├── total_xy.gif
├── total_xz.gif
├── total_yz.gif
├── system_xy.gif
├── system_xz.gif
├── system_yz.gif
├── wp_xy.gif
├── wp_xz.gif
├── wp_yz.gif
├── total_density_volume.gif
├── system_density_volume.gif
└── wp_density_volume.gif
```

## Fixed colour-scale rule

All density GIFs must use a fixed colour scale over time.

Do not let each frame autoscale independently.

For example:

- `wp_xz.gif` must use one fixed scale for all WP xz frames.
- `total_density_volume.gif` must use one fixed transfer function / colour range over the full time series.

---

# 17. `results/analysis/screens/`

## Purpose

LEED/screen visualisations belong here.

## Category layout

```text
results/analysis/screens/
├── total/
├── instantaneous/
├── time_windowed/
├── filtered/
├── spectra/
└── coordinate_checks/
```

---

## 17.1 Total screen plots

```text
results/analysis/screens/total/
├── all_screens_grid.png
├── all_screens_grid_log.png
├── screen_00.png
├── screen_00_log.png
├── screen_01.png
├── screen_01_log.png
└── ...
```

---

## 17.2 Instantaneous screen animations

```text
results/analysis/screens/instantaneous/
├── screen_00_time_evolution.gif
├── screen_01_time_evolution.gif
├── screen_00_contact_sheet.png
├── screen_01_contact_sheet.png
└── ...
```

All instantaneous screen GIFs must use a fixed colour scale through time.

---

## 17.3 Time-windowed screen plots

This folder must be flat.

```text
results/analysis/screens/time_windowed/
├── screen_00_t000000_to_t000300.png
├── screen_00_t000000_to_t000300_log.png
├── screen_01_t000000_to_t000300.png
├── screen_01_t000000_to_t000300_log.png
├── screen_00_t000300_to_t000600.png
├── screen_00_t000300_to_t000600_log.png
└── ...
```

Do not use:

```text
results/analysis/screens/time_windowed/window_00/
results/analysis/screens/time_windowed/screen_00/
```

---

## 17.4 Filtered screen plots

```text
results/analysis/screens/filtered/
├── screen_00_hard_mask.png
├── screen_00_radial_background_subtracted.png
├── screen_00_fitted_cross_subtracted.png
├── screen_00_template_subtracted.png
├── screen_00_fft_notch_filtered.png
├── screen_00_raw_over_background_ratio.png
└── ...
```

For filtered time-windowed screens, include the timestep window in the filename:

```text
screen_10_t000300_to_t000600_radial_background_subtracted.png
```

---

## 17.5 Screen spectra

```text
results/analysis/screens/spectra/
├── leed_fft_screen_00.png
├── leed_fft_screen_01.png
└── ...
```

For time-windowed spectra:

```text
leed_fft_screen_00_t000000_to_t000300.png
```

---

## 17.6 Coordinate-check plots

Use this folder to debug the LEED plotting coordinate mapping:

```text
results/analysis/screens/coordinate_checks/
├── screen_00_raw_index_plot.png
├── screen_00_coordinate_mapped_plot.png
├── screen_00_axis_orientation_check.png
└── ...
```

This is specifically meant to catch the failure mode where the LEED pattern is split into four pieces and appears at the edges of the image.

---

# 18. `results/analysis/overlap/`

## Purpose

Overlap visualisations belong here.

These outputs should show both:

1. the evolution of the full overlap matrix; and
2. the overlap of the wavepacket state with all ground-state KS orbitals over time.

## Category layout

```text
results/analysis/overlap/
├── overlap_matrix_heatmap.gif
├── wp_overlap_with_gs_orbitals.gif
├── wp_overlap_heatmap.png
├── orbital_survival_heatmap.png
└── dominant_overlap_vs_time.png
```

## Required WP-overlap GIF

The file:

```text
wp_overlap_with_gs_orbitals.gif
```

is required for WP runs with overlap tracking.

It must show how the injected wavepacket overlaps with the original ground-state KS orbitals as a function of time.

This can be visualised as:

- a bar chart animation over GS orbital index;
- a heatmap animation; or
- a scrolling/time-resolved matrix view.

The colour scale must be fixed through time.

---

# 19. `results/analysis/orbitals/`

## Purpose

Orbital visualisation outputs belong here.

## Category layout

```text
results/analysis/orbitals/
├── homo_lumo_density.png
├── selected_orbital_gallery.png
├── selected_orbital_slices.gif
└── wp_orbital_density.gif
```

## Fixed colour-scale rule

Any orbital GIF made from time-stepped orbital data must use a fixed colour scale through time.

---

# 20. Complete preferred tree

```text
results/
├── run_summary.txt
│
├── raw/
│   ├── ground_state/
│   │   ├── summary.txt
│   │   ├── scf_history.csv
│   │   ├── eigenvalues.csv
│   │   ├── occupations.csv
│   │   ├── ground_state_observables.csv
│   │   ├── density_system/
│   │   ├── density_gs_orbitals/
│   │   └── checkpoint/
│   │
│   ├── wavepacket/
│   │   ├── wavepacket_config.txt
│   │   ├── injection_report.txt
│   │   ├── orthogonality_report.csv
│   │   ├── density_wp_initial/
│   │   └── wavefunction_wp_initial/
│   │
│   ├── density/
│   │   ├── density_rt_total/
│   │   ├── density_rt_system/
│   │   └── density_rt_wp/
│   │
│   ├── orbitals/
│   │   ├── homo/
│   │   ├── lumo/
│   │   ├── wp/
│   │   └── selected_orbitals/
│   │
│   ├── screens/
│   │   ├── screen_config.csv
│   │   ├── window_ranges.csv
│   │   ├── total/
│   │   ├── instantaneous/
│   │   └── time_windowed/
│   │
│   ├── overlap/
│   │   ├── index.csv
│   │   ├── overlap_000000.csv
│   │   ├── overlap_000001.csv
│   │   └── ...
│   │
│   ├── observables/
│   │   ├── observables.csv
│   │   ├── energy_components.csv
│   │   ├── current_components.csv
│   │   ├── dipole_components.csv
│   │   ├── norm_conservation.csv
│   │   ├── charge_conservation.csv
│   │   ├── wp_center_of_mass.csv
│   │   ├── wp_width.csv
│   │   ├── fft_total_energy.csv
│   │   ├── fft_current_x.csv
│   │   ├── fft_current_y.csv
│   │   ├── fft_current_z.csv
│   │   ├── dipole_spectrum_x.csv
│   │   ├── dipole_spectrum_y.csv
│   │   └── dipole_spectrum_z.csv
│   │
│   └── vti/
│       ├── density_gs_system/
│       ├── density_gs_orbitals/
│       ├── density_rt_total/
│       ├── density_rt_system/
│       ├── density_rt_wp/
│       └── orbitals/
│
└── analysis/
    ├── ground_state/
    ├── observables/
    ├── density/
    ├── screens/
    │   ├── total/
    │   ├── instantaneous/
    │   ├── time_windowed/
    │   ├── filtered/
    │   ├── spectra/
    │   └── coordinate_checks/
    ├── overlap/
    └── orbitals/
```

---

# 21. Implementation rules for Claude Code

## 21.1 Do not over-template existing files

Most file contents are already defined by previous runs. The priority is to categorise them correctly.

Do not rewrite all writers just to match new exact column templates.

Only `run_summary.txt` needs the explicit template in this document.

---

## 21.2 Do not create forbidden top-level metadata files

Do not create:

```text
manifest.txt
simulation_config.txt
derived_values.txt
command_line.txt
git_info.txt
timing_summary.txt
```

Put that information in:

```text
results/run_summary.txt
```

---

## 21.3 Do not create `processed/` or `derived_data/`

Use:

```text
results/analysis/
results/raw/observables/
```

---

## 21.4 Keep time-windowed screen folders flat

Forbidden:

```text
results/raw/screens/time_windowed/window_00/
results/analysis/screens/time_windowed/window_00/
```

Required style:

```text
results/raw/screens/time_windowed/screen_00_t000000_to_t000300.dat
results/analysis/screens/time_windowed/screen_00_t000000_to_t000300.png
```

---

## 21.5 Keep ParaView VTI series flat

Real-time density VTI series must be selectable directly in ParaView from one folder.

Good:

```text
results/raw/vti/density_rt_wp/density_rt_wp_t000000.vti
results/raw/vti/density_rt_wp/density_rt_wp_t000010.vti
results/raw/vti/density_rt_wp/density_rt_wp_t000020.vti
```

Bad:

```text
results/raw/vti/density_rt_wp/t000000/density.vti
results/raw/vti/density_rt_wp/t000010/density.vti
```

Ground-state orbital VTI files should also be flat:

```text
results/raw/vti/density_gs_orbitals/orbital_0000.vti
results/raw/vti/density_gs_orbitals/orbital_0001.vti
```

---

## 21.6 Preserve system/WP/total density distinction

Always distinguish:

```text
density_rt_system
density_rt_wp
density_rt_total
```

Do not only write system and WP densities.

---

## 21.7 Fix LEED coordinate mapping before interpreting patterns

If a LEED pattern appears split into four quadrants at the image edges, treat it as a coordinate-mapping or plotting issue first.

Do not interpret this as physics until the following are checked:

- array orientation
- transpose convention
- origin convention
- coordinate extent
- whether the detector axes are centred
- whether any shift/roll operation was accidentally applied
- whether the `.dat` loader is reading coordinate columns correctly

Add coordinate-check plots under:

```text
results/analysis/screens/coordinate_checks/
```

---

## 21.8 Add WP overlap with GS KS orbitals at all timesteps

For WP runs with overlap tracking, the analysis must be able to visualise:

```text
wavepacket overlap with all ground-state KS orbitals as a function of time
```

Raw data source:

```text
results/raw/overlap/
```

Required analysis output:

```text
results/analysis/overlap/wp_overlap_with_gs_orbitals.gif
```

This GIF must use a fixed colour scale or fixed y-axis scale through time.

---

## 21.9 Use fixed colour scales for all time-stepped GIFs

For every GIF/movie made from timestep data:

```text
compute global colour limits first
render every frame using those limits
```

This applies to density, screens, overlap, and orbital animations.

---

# 22. Practical shell checks

## Check top-level files

```bash
find results -maxdepth 1 -type f -print
```

Expected:

```text
results/run_summary.txt
```

---

## Check for forbidden top-level files

```bash
find results -maxdepth 1 \( \
  -name 'manifest.txt' -o \
  -name 'simulation_config.txt' -o \
  -name 'derived_values.txt' -o \
  -name 'command_line.txt' -o \
  -name 'git_info.txt' -o \
  -name 'timing_summary.txt' \
\) -print
```

This should print nothing.

---

## Check for forbidden folders

```bash
find results -type d \( -name 'processed' -o -name 'derived_data' \) -print
```

This should print nothing.

---

## Check that time-windowed raw screens are flat

```bash
find results/raw/screens/time_windowed -mindepth 1 -type d -print
```

This should print nothing.

---

## Check that time-windowed analysis screens are flat

```bash
find results/analysis/screens/time_windowed -mindepth 1 -type d -print
```

This should print nothing.

---

## Check for nested VTI files that may be inconvenient for ParaView

```bash
find results/raw/vti -mindepth 3 -name '*.vti' -print
```

This should normally print nothing for the first implementation, because VTI series should be flat within their category folders.

---

# 23. Summary of the most important rules

1. Only one top-level summary file:

   ```text
   results/run_summary.txt
   ```

2. Use:

   ```text
   results/raw/
   results/analysis/
   ```

3. Do not use:

   ```text
   processed/
   derived_data/
   ```

4. Put raw and derived numerical CSVs in:

   ```text
   results/raw/observables/
   ```

5. Put observable plots in:

   ```text
   results/analysis/observables/
   ```

6. Keep time-windowed screen files flat:

   ```text
   screen_00_t000000_to_t000300.dat
   screen_00_t000000_to_t000300.png
   screen_00_t000000_to_t000300_log.png
   ```

7. Keep VTI files flat within their ParaView series folder.

8. Always distinguish:

   ```text
   density_rt_system
   density_rt_wp
   density_rt_total
   ```

9. Treat split-at-the-edges LEED images as a coordinate-mapping bug until proven otherwise.

10. Add WP overlap against all GS KS orbitals over time, and visualise it as:

    ```text
    results/analysis/overlap/wp_overlap_with_gs_orbitals.gif
    ```

11. All time-stepped GIFs and movies must use a fixed colour scale through time.