# Test catalogue — inqkit (C++) + inqview (Python)

A single index of every test in the inq-stack suite: what each test does and
what it asserts against. Companion to `inqkit-tests.md` (designs),
`inqkit-errors.md` / `inqview-findings.md` (bugs + decisions), and the
validation dossiers in this folder.

## How validation occurs (read first)

- **Expected values are derived, not captured.** Every assertion compares code
  output to a value known *before* the test was written — an analytic result,
  a unit conversion, a conservation law, or a documented table. Never "assert ==
  whatever the code prints now" (anti-circularity, IV-M10).
- **Tiers.**
  - inqkit **pure** — C++17 only, no INQ engine; runs anywhere in seconds.
  - inqkit **engine** — links the INQ engine (CUDA/MPI/libxc); GPU; runs on the
    lab box. Includes one real `mpirun -np 2` test.
  - inqview is **entirely pure-Python/numpy** (ADR 0005) — no GPU, no INQ, small
    committed fixtures; `np.allclose` tolerances (never bit-exact).
- **`xfail` = an executable to-do.** A known bug is encoded as a `strict xfail`
  asserting the CORRECT behaviour. It stays "expected-fail" until the fix lands,
  then becomes an unexpected-pass and pytest demands the marker be removed. The
  suite is the bug list.
- **Status snapshot (2026-06-10):** inqkit 26 pure + 19 engine cases green;
  inqview 29 passed + 2 xfailed (the two captured bugs) + `test_lindhard`
  (8 passed / 5 xfail, dynamical-Lindhard sign WIP).

---

# inqkit (C++ / Catch2)

## Pure tier — `inq-stack/tests/cpp/*.cpp` (9 files, 26 cases)

| Test case | What it does | Asserts against |
|---|---|---|
| `fft_shift_index: even-size table (size 6)` | Maps every index of a size-6 axis | The documented even-size shift table |
| `fft_shift_index: odd sizes keep origin at 0` | Shifts an odd-length axis | Physical origin stays at FFT index 0 |
| `fft_shift_index: origin maps to 0 (all parities)` | Index of physical origin | →0 for even and odd sizes |
| `fft_shift_index: bijection` | Shifts all indices | Permutation of `[0,size)` — no collisions |
| `ComplexField3DWriter: _real/_imag round-trip` | Writes then reads a complex field | Real/imag raw arrays byte-preserved |
| `flatten_index: x-slowest z-fastest` | Linearises a 3D index | Documented memory ordering |
| `flatten_index: unique contiguous index` | Flattens every cell | Bijection onto `[0,N)` |
| `step_suffix: 6-digit tag` | Formats a step number | Zero-padded 6-digit string |
| `RealField3DWriter: .raw+.meta round-trip` | Writes then reloads a real field | Field values + metadata preserved |
| `jellium shells: shell table to N=162` | Fills closed shells | Documented magic-number shell table |
| `jellium shells: partial last shell` | Requests a non-magic count | Last shell truncated to `n_states` |
| `DensityDelta: lazy reference + zero` | First snapshot | Captures ref, returns Δn=0 |
| `DensityDelta: t=0 fixed base` | Many snapshots | Δ always vs t=0, not rolling previous |
| `DensityDelta: explicit set_reference` | Manual vs lazy ref | Identical Δ result |
| `DensityDelta: L2 scales with dV` | Computes ‖Δn‖₂ | Sums cells × voxel volume |
| `DensityDelta: grid mismatch throws` | Wrong-shape reference | Throws |
| `center_of_density: uniform centroid` | Uniform 1D field | Centroid at geometric centre |
| `center_of_density: weighted centroid` | Non-uniform field | ∫r·n/∫n in Bohr |
| `center_of_density: zero field guard` | All-zero field | Stays at origin (w>0 guard) |
| `center_of_density: axes not transposed` | Distinct x/y/z weights | Three distinct coordinates |
| `center_of_density: empty field throws` | Empty input | Throws |
| `VTIImageDataWriter: x-fastest reorder (T06)` | Writes a VTI | Axis mapping preserved through reorder |
| `Vec3: default zero` | Default-constructed Vec3 | (0,0,0) |
| `Vec3: dot/norm2/norm` | Vector algebra | Analytic dot/length values |
| `Vec3: arithmetic operators` | +, −, scalar × | Component-wise results |
| `Vec3: project onto unit direction` | k̂ projection | Analytic scalar projection |

## Engine tier — `inq-stack/tests/cpp/engine/*.cpp` (18 files, 20 cases)

| Test case | What it does | Asserts against |
|---|---|---|
| `density::total integrates to N` | He GS density | ∫n dV = 2; grid dims match basis |
| `coord_mapping Test B: EVEN grid (L=10.0)` | Off-centre He, even grid | density argmax within dx of atom (node convention) |
| `coord_mapping Test B: ODD grid (L=10.5)` | Off-centre He, odd grid | Same, odd parity |
| `orbital::wavefunction normalised` | Extracts a KS orbital | Complex field, ∫|ψ|²=1 |
| `density semantics T02/E02` | Inject WP, then propagate | `density()` cache is stale post-inject, refreshed by propagation |
| `orthogonalisation T29/E03` | Strong-overlap WP inject | Single-pass MGS limit; iterated pass needed (E03 fix) |
| `ObservablesWriter header+row` | Writes observables.csv | Header/row match the `ObservableSelection` |
| `PlaneScreen::extract z=0 slice` | Single-rank slice of He | Slice equals the z=0 plane density |
| `plane_screen parallel (E01)` | `mpirun -np 2` slice | Slice agrees across ranks (Allreduce fix) |
| `LeedPatternAccumulator = Σ slice·dt` | Accumulates screen frames | Pattern = time-integrated slice |
| `OrbitalOverlapMatrix t=0 identity` | Overlap at t=0 | Identity block + WP column |
| `OccupationsWriter (mock Viewables)` | He+WP occupations | Occupation row via duck-typed mock |
| `MomentumDistribution peaks at |k0|` | Inject WP at k₀ | Binned n_wp(|k|) peak in the |k₀| bin |
| `wp_momentum two-route (T28/T04)` | ⟨p⟩ real-space vs reciprocal | Both routes agree and equal k₀ |
| `WPMomentumStats::compute() moments` | compute() on injected WP | ⟨p⟩, ⟨p²⟩ of the WP; N>0 |
| `wp_real_space moments` | Injected WP real-space | ⟨r⟩ = centre, Var = σ²/2 |
| `WPRealSpaceStats::compute()` | compute() route | ⟨r⟩ = centre, N≈1 (node convention; caught E04) |
| `StateEnergyWriter (short propagate)` | Per-state ⟨ψ|H|ψ⟩ via RT | n_states finite E_expect rows per step |
| `dump_eigenvalues He GS table` | GS eigenvalue dump | Eigenvalue table written + sane |
| `free WP: non-interacting Gaussian` *(integration)* | Inject free Gaussian WP, propagate non-interacting T=3 a.u. | Analytic free-particle: ⟨z⟩=k₀T, Var=σ²/2+T²/(2σ²), ⟨p⟩=k₀, norm & E_kin conserved, moved+spread (IV-M11) |
| `StateNormWriter: every orbital norm ≈ 1` *(feature)* | compute() per-state ∫|ψ_i|²dV for He+2 extra+WP | one entry per state; every norm ≈ 1; WP slot ≈ 1 (norm-per-state diagnostic) |

---

# inqview (Python / pytest)

All pure tier (ADR 0005). Location: `inq-stack/python/tests/` (+ the existing
`inq-stack/python/inqview/postprocess/test_lindhard.py`).

## `test_fourier.py` — windowed FFT kernel (10 cases, 1 xfail)

| Test | What it does | Asserts against |
|---|---|---|
| `peak_is_at_f0[boxcar/hann/hamming/blackman]` | FFT a known on-bin tone | Peak at `f0` within ½ bin, every window |
| `boxcar_amplitude_is_A` | FFT a tone, boxcar | One-sided amplitude == A (exact) |
| `hann_amplitude_is_A_pending_iv_e03` | FFT a tone, Hann | Amplitude == A — **xfail until IV-E03 `/win.sum()`** |
| `hann_amplitude_currently_reduced…` | Same, current code | Amplitude ≈ A·mean(hann) ≈ 0.5A (pins the gap) |
| `zero_pad_preserves_peak` | zero_pad 1 vs 4 | Peak f0 + amplitude unchanged |
| `dc_offset_hijacks_peak_without_detrend` | Tone + 50 DC, no detrend | Max moves to ω≈0 (todo.txt #2 bug) |
| `detrend_removes_dc_and_recovers_f0` | Tone + offset + drift, detrend | True peak at f0 restored |

## `test_kl_divergence.py` — KL helpers (6 cases)

| Test | What it does | Asserts against |
|---|---|---|
| `normalise_sums_to_one` | Normalise a histogram | Σ=1; `[1,1,2]→[.25,.25,.5]` |
| `normalise_zero_input` | All-zero histogram | Returns zeros (no /0) |
| `kl_self_is_zero` | KL(P‖P) | =0 (Gibbs equality) |
| `kl_known_value_nats` | KL([.5,.5]‖[.25,.75]) | =0.5ln2+0.5ln(2/3) ≈0.14384 nats |
| `kl_nonnegative` | 20 seeded random P,Q | KL ≥ 0 (Gibbs) |
| `kl_is_asymmetric` | KL(P‖Q) vs KL(Q‖P) | Not equal (not a metric) |

## `test_screens_io.py` — LEED loader (4 cases)

| Test | What it does | Asserts against |
|---|---|---|
| `header_fields_parsed` | Load hand-built `.dat` | label/z/time/n_accum/nx/ny/dx/dy |
| `fftshift_moves_origin_peak_to_centre` | Corner peak at FFT (0,0) | Lands at centre `[2,2]`; corners empty; mass conserved |
| `origin_and_extent_are_centred` | Loader origin override | origin=−L/2; extent `(−2,2,−2,2)` |
| `missing_file_raises` | Nonexistent path | `FileNotFoundError` |

## `test_fields.py` — field dataclasses (7 cases)

| Test | What it does | Asserts against |
|---|---|---|
| `meta_derived_geometry` | FieldMeta props | shape, num_points=24, voxel volume, bytes, is_real |
| `meta_bad_dtype_raises` | Invalid dtype | ValueError |
| `realfield_accepts_matching_array` | Valid RealField3D | shape + min/max/mean |
| `realfield_rejects_shape_mismatch` | Wrong shape | ValueError |
| `realfield_rejects_non_float` | int array | ValueError |
| `complexfield_magnitude_phase_array` | real=3, imag=4 | |ψ|=5, phase=atan2(4,3), array=3+4j |
| `complexfield_rejects_shape_mismatch` | Mismatched parts | ValueError |

## `test_wake.py` — shared colour scale (3 cases)

| Test | What it does | Asserts against |
|---|---|---|
| `symmetric_about_zero` | shared_clim over arrays | (−m,m), m = global max|value| |
| `asymmetric_mode_starts_at_zero` | symmetric=False | (0, max) |
| `percentile_clips_lone_spike` | 99×1 + one 1000 | vmax = 95th pct, spike suppressed |

## `test_theme.py` — canonical visualisation theme (6 cases) — FEATURE built 2026-06-10

| Test | What it does | Asserts against |
|---|---|---|
| `cmap_roles_are_the_designed_values` | `cmap_for(role)` | sequential→inferno, diverging→RdBu_r, phase→twilight |
| `unknown_cmap_role_raises` | bad role | ValueError |
| `fixed_dimension_constants` | ONE_COL_IN / TWO_COL_W_IN | (3.5,3.0) / 7.0 |
| `figure_one_col_has_exact_size_and_fixed_axes` | one-col figure | size==(3.5,3.0); axes rect fixed (panels align) |
| `figure_two_col_is_seven_inches_wide` | two-col figure | width==7.0 |
| `apply_theme_installs_designed_rcparams` | apply_theme() | font 10, axes lw 0.8, ticks in, cmap inferno |

## `test_energy_components.py` — functional energy flow (6 cases) — FEATURE built 2026-06-10

| Test | What it does | Asserts against |
|---|---|---|
| `external_is_recovered_as_residual` | E_ext = total−(kin+H+xc) | known external term |
| `components_sum_to_total_exactly` | Σ components | == E_total (exact) |
| `delta_components_sum_to_delta_total` | Σ ΔE | == ΔE_total |
| `conserved_total_gives_zero_drift` | constant E_total | ΔE_total ≈ 0 |
| `redistribution_in_ev` | ΔE per component (eV) | ΔE_kin=0.5 Ha·27.2114; components sum to total |
| `missing_column_raises` | drop a column | ValueError |

## `test_energy_components_render.py` — energy-flow renderers (4 cases) — FEATURE built 2026-06-10

Data-contract tests (not pixels, ADR-0005): assert the rendered artists carry
EXACTLY the dataclass numbers, proving render consumes / never recomputes.

| Test | What it does | Asserts against |
|---|---|---|
| `kernel_decomposition_known` | guard the synthetic inputs | E_ext=[1,2,2]; Σ components == E_total |
| `bars_heights_equal_breakdown` | `render_initial_vs_final_bars` | two BarContainers' heights == `breakdown('initial'/'final')` |
| `flow_lines_ydata_equal_dE_in_eV` | `render_flow_lines` | each line's y-data == `dE_<c> · HA_TO_EV` (kin/H/xc/ext/total) |
| `breakdown_gif_writes_file` | `render_breakdown_gif` to tmp | GIF file exists, non-empty (skips if no Pillow writer) |

## `test_wp_integrity_from_run.py` — from-run assembly (4 cases) — GLUE built 2026-06-10

Builds a synthetic run dir (real CSV layout) and checks `assemble_from_run`
against analytically known values.

| Test | What it does | Asserts against |
|---|---|---|
| `assemble_time_and_sigma` | parse momentum + real-space CSVs | time=[0,0.04]; σ_r=√(Σσ²)=[√1.5, 2] |
| `assemble_kl_known` | KL of per-step n_wp vs initial | kl[0]=0; kl[1]=½ln2=0.346574 |
| `assemble_ipr_is_nan` | ipr unavailable (no WP density VTI) | all NaN (documented) |
| `assemble_reference_previous_first_is_zero` | frame-to-frame mode | kl[0]=0; kl[1]=½ln2 |

## `test_wp_integrity.py` — WP-integrity metrics (8 cases) — FEATURE built 2026-06-10

| Test | What it does | Asserts against |
|---|---|---|
| `kl_self_zero_and_known_value` | momentum KL | 0 for P‖P; 0.14384 nats known case |
| `kl_nonnegative_unnormalised_inputs` | raw histograms | KL ≥ 0 |
| `ipr_localised_greater_than_delocalised` | IPR spike vs uniform | spike > uniform |
| `ipr_uniform_closed_form` | IPR of uniform-N | == 1/N |
| `variance_recovers_gaussian_width` | weighted variance of a Gaussian | ≈ σ² |
| `variance_zero_density_is_zero` | empty density | 0 |
| `wpintegrity_holds_series` | dataclass container | shapes + IPR decay |

## `test_plasmon_spectrum.py` — plasmon peak-locator (4 cases) — FEATURE built 2026-06-10

| Test | What it does | Asserts against |
|---|---|---|
| `phasor_peaks_at_omega_p_for_all_q` | time-FFT of e^{−iω_p t} | peak at ω_p for every q-mode |
| `one_over_q_squared_weighting` | equal-amplitude phasors | loss·q² constant (the 1/q² weight) |
| `complex_fft_separates_plus_and_minus_omega` | complex FFT | power on ONE side of ω (IV-E01 fix) |
| `axial_extraction_picks_the_right_mode` | δn=cos(2πmz/nz) | only axial mode m populated; q_m correct |

## `test_center_of_density.py` — Python COD (5 cases) — FEATURE built 2026-06-10

| Test | What it does | Asserts against |
|---|---|---|
| `single_point_centroid_is_its_node_coordinate` | δ at index (i,j,k) | origin+(i,j,k)·dx (node convention) |
| `gaussian_centroid_recovers_centre` | Gaussian density | its centre (tol absorbs edge truncation) |
| `half_cell_convention_is_offset_by_dx/2` | node vs (i+½)·dx | difference == dx/2 exactly (documents E04) |
| `bath_is_total_minus_wp` | WP/total/bath compare | wp≈+2, bath≈−2, total between (IV-M02) |

## `test_wp_integrity.py` — kl_series additions (built 2026-06-10)
Added `kl_series` (drift-from-launch + frame-to-frame rate): starts at 0; rises
from launch as the WP drifts; 0 frame-to-frame for a steady WP; bad reference raises.

## `test_fourier.py` — subtract= additions (built 2026-06-10)
IV-E03 coherent-gain fix flipped the Hann xfail to a real pass
(`windowed_amplitude_is_A_after_coherent_gain_fix`). Added IV-M12 `subtract=`
cases: initial/mean/detrend remove DC & keep the peak; only detrend suppresses a
genuine-drift low-band leakage; `none` lets a DC offset hijack the peak; bad
mode raises; `FourierResult.subtract` records the choice.

## `test_deps_clean.py` — ADR-0003 invariant (3 cases, ENFORCED 2026-06-10)

The package split landed: the xfail flipped to live parametrized cases.

| Test | What it does | Asserts against |
|---|---|---|
| `analysis_import_is_matplotlib_and_vtk_free[inqview]` | Subprocess `import inqview` | matplotlib/VTK absent (lazy top-level `__init__`) |
| `…[inqview.analysis]` | Subprocess `import inqview.analysis` | matplotlib/VTK absent (clean analysis pkg) |
| `…[inqview.analysis.fourier]` | Subprocess imports the moved kernel | matplotlib/VTK absent |

## `test_lindhard.py` — analytic Lindhard/RPA (existing; 8 pass, 5 xfail)

| Test | What it does | Asserts against |
|---|---|---|
| `constants_rs_5p69` | r_s=5.69 derived constants | kF=0.337, ω_p=0.1276 Ha, E_F |
| `static_limit` | χ⁰(q→0,ω=0) | = −N(E_F) = −kF/π² (Giuliani-Vignale 4.39) |
| `plasmon_zero_q_limit` | ω_pl(q→0) | → ω_p |
| `plasmon_dispersion_positive_slope` | ω_pl(q) | Strictly increasing (Bohm-Gross) |
| `imag_zero_outside_continuum` | Im χ⁰ above e-h continuum | ≈ 0 |
| `imag_nonzero_inside_continuum` | Im χ⁰ inside continuum | < 0 |
| `chi0_array_broadcasting` | Array inputs | Correct broadcast shapes |
| `f_sum_rule[…]`, `high_omega_limit`, `bethe_limit` | Dynamical Lindhard | **xfail** — high-ω sign WIP |

---

## Tests still to add (per the plan + findings)

- **inqkit** — `test_free_wp_engine.cpp`: free Gaussian WP, non-interacting,
  assert analytic σ_r(t)/⟨p⟩=k₀/ballistic centroid/norm/energy (IV-M11).
- **inqview** — once kernels are built in the restructure: `plasmon_spectrum`
  (peak@ω_p, 1/q²), `center_of_density` (E04 dx/2 cross-check), `wp_integrity`
  (free-WP σ_r(t)), `energy_components` (Σ==E_total), `gs_projected` (t=0
  identity), `theme` (figure_one_col size, role→cmap); the free-space-WP
  integration fixture; and the IV-M12 `subtract=` baseline test.
