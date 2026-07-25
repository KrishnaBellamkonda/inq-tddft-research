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
| `test_heuristics: electron_gas_scales_analytic` | HEG scales for r_s=4 | k_F, ω_p=√(3/r_s³), k_TF, λ_F, n₀ analytic identities |
| `test_heuristics: timescales_constant_velocity` | projectile entry/exit/box-edge times | (z_face−z₀)/v closed form (slab END = far face) |
| `test_heuristics: zero_point_sigma_half` | WP zero-point KE 3/(4σ²) | 3.0 Ha (=81.6 eV) for σ=0.5; charge std σ/√2 |
| `test_heuristics: norm_absorption_split` | total/orbital/bath-overflow absorbed | N(0)−N(t) split arithmetic |
| `test_heuristics: spreading_factor` | σ_z(t)/σ_z(0) + max | ratio + max identities |
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
| `T0.4 annulus integrates to n0·π(Rout²−Rin²)·Lz` *(cylindrical-jellium)* | Build annulus n₊ (R_in=3, R_out=7, w=0, tube∥z, L=20) | ∫n₊ = n₀·π(R_out²−R_in²)·L_z = 25.13 (ε=0.08; two curved radial surfaces) — PASS |
| `T0.5 annulus = outer cyl − inner cyl (bore carved)` *(cylindrical-jellium)* | ∫annulus(3,7) + ∫cyl(0,3) vs ∫cyl(0,7); solid cyl = annulus(R_in=0,w=0) | additivity identity holds (ε=0.05): inner complement carves the correct bore, no sign/radius error — PASS |
| `T0.6 erfc-smoothed (w=1) annulus conserves charge` *(cylindrical-jellium)* | Production-config annulus (w≈1, 4-Bohr wall ≫ w) | ∫n₊ ≈ sharp-edge target within ε=0.10 (erfc edges shift ≲ few %; caller rescales n₀ for exact neutrality) — PASS |

> Annular `background_shape` (cylindrical-jellium campaign): the erfc edge profile,
> ½-height crossovers at both radial edges, and axial uniformity are proven
> ANALYTICALLY by the `formula-validation` agent (VERDICT: CONFIRM, 2026-06-28);
> T0.4–T0.6 above validate the grid builder via decomposition-safe integrals.

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

## inqkit full energy-component streaming — FEATURE built 2026-07-07

`ObservablesWriter` extended to stream every INQ energy accessor per step
(`external`, `non_local`, `ion`, `ion_kinetic`, `exact_exchange`, `nvxc`,
`eigenvalues`) alongside the pre-existing `total/kinetic/hartree/xc`. Additive:
new `ObservableSelection` flags default off, so existing runs' CSV schema is
unchanged. Wiring: `step_context.hpp` (fields) → `real_time_session.hpp` (copies
from `data.energy().<accessor>()`) → `observables_writer.hpp` (cols/vals). No
`inq/`/`inq-study/` edit (reads public accessors only).

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| components reconstruct total | classical `buildsmoke` observables.csv (2026-07-07) | integration (real run) | `total == kinetic+external+non_local+hartree+xc+exact_exchange+ion+ion_kinetic`; abs diff 7.8e-13 Ha (machine precision) | PASS |
| schema self-describing | header row of the same CSV | integration | 12 energy columns present in declared order; per-run `sel.energy_*=true` in `{classical,wp}/run.cpp` | PASS |

## wp/run.cpp t=0 density save (`LJ_SAVE_DENSITY`) — FEATURE built 2026-07-08

Env-gated block in `campaign_autorun/wp/run.cpp` writes `density_wp` (=|ψ_WP|² via
`density::orbital`), `density_total`, `density_bath` as VTIs at t=0 for the
screening/WP-potential test. Reads public inqkit accessors only; no `inq/`/`inq-study/` edit.

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| WP charge normalised | `runs/screening_wp/wp_r{4,12}_p2/.../density_wp.vti` | integration (real run) | `∫n_WP·dV = 1.0000` e (both radii), loads via `inqview.load_vti` (physical order) | PASS |
| WP source ≈ ideal Gaussian | same VTI vs analytic | numeric | radial `n_WP` tracks Gaussian(σ_ρ=σ_WP/√2=0.354) at core (ratio ~1.0–1.1, r<1.2 Bohr) | PASS |
| FFT-Poisson validated | notebook WPPOT cell | numeric | Python periodic FFT-Poisson reproduces analytic `erf(r/(√2σ_ρ))/r` for a unit Gaussian to RMS 3.4e-4 Ha (0.4–6 Bohr) before use on n_WP | PASS |
| screening baseline zero | `density_total − density_gs_system` | integration | `n_slab(t=0) − n_GS` bit-identical, max|Δ| = 0.0e+00 (no instantaneous screening); NOTE `density::total` returns 82 e slab-only (WP excluded — `density.hpp` TODO) | PASS |

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

## `test_efield.py` — FFT-Poisson E-field kernel (4 cases) — FEATURE built 2026-06-17

Formula-bearing: independently CONFIRMED by the formula-validation agent (4π
factor, −iGφ̃ gradient sign, ρ=−n direction, G=0 neutralizing background, analytic
single-mode field). For the CAP-in-jellium baselines (idea 1). Native atomic units.

| Test | What it does | Asserts against |
|---|---|---|
| `uniform_density_gives_zero_field` | uniform n | E≡0 (G=0 removed = neutralizing background) |
| `cosine_density_matches_analytic_field` | n=n₀+A cos(G₀z) | E_z = −(4πA/G₀) sin(G₀z) (machine precision) |
| `gaussian_charge_matches_isolated_erf_field` | Gaussian charge σ≪L | isolated erf field f(r)/r² (rel 5%) |
| `si_units_are_a_constant_rescale` | atomic vs SI | one constant factor 5.142e11 V/m |

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

## `test_lindhard_elf.py` — corrected ELF + point-charge reference (16 pass, 2026-06-14)

| Test | What it does | Asserts against |
|---|---|---|
| `thomas_fermi_static_limit` | ε(q→0,0) | 1 + k_TF²/q² |
| `f_sum_rule[q]` | ∫ω Im[−1/ε]dω | (π/2)ω_p² to <1 % at all q (plasmon retained) |
| `imag_eps_only_in_continuum` | Im ε support | 0 outside e-h continuum + plasmon |
| `stopping_qgrid_convergence` | S_LR(v;σ) q-grid | <0.5 % under n_q refinement |
| `stopping_positive_and_monotone_lowv` | low-v friction | 0 < S(0.2) < S(0.4) |
| `stopping_point_converged[v]` | **point-charge reference** | <1 % under qmax-margin + n_q refinement (the one analytical curve) |
| `stopping_point_above_finite_sigma` | point vs σ-suppressed | S_point > S(σ=0.2,0.35,0.5) at fixed v |

---

## Tests still to add (per the plan + findings)

- **inqkit** — `test_free_wp_engine.cpp`: free Gaussian WP, non-interacting,
  assert analytic σ_r(t)/⟨p⟩=k₀/ballistic centroid/norm/energy (IV-M11).
- **inqview** — once kernels are built in the restructure: `plasmon_spectrum`
  (peak@ω_p, 1/q²), `center_of_density` (E04 dx/2 cross-check), `wp_integrity`
  (free-WP σ_r(t)), `energy_components` (Σ==E_total), `gs_projected` (t=0
  identity), `theme` (figure_one_col size, role→cmap); the free-space-WP
  integration fixture; and the IV-M12 `subtract=` baseline test.

## Absorbing boundary — mask absorber (free_wp, 2026-06-13)

| Test | File | Tier | Known-case oracle | Status |
|---|---|---|---|---|
| gate1: mask shape | `ResearchProject/systems/vacuum/tests/gate1_mask_absorber/run.cpp` (T1) | unit (host) | sin² mask: M(z0)=1, M(z0+L/2)=0.5, M(z0+L)=0, monotone — analytic, tol 1e-12 | PASS |
| gate1: ε reducer | same (T2) | unit (engine) | symmetric Gaussian split (mid-grid cut): ε_all=1, ε_none=0, ε_half=0.5 — analytic, tol 1e-3 | PASS |
| gate1: mask fidelity | same (T3) | integration (engine) | M≡1 absorber ⇒ trajectory bit-identical to baseline (|ΔN|,|Δz|<1e-9) | PASS |
| gate1: mask feedthrough | same (T4) | integration (engine) | sin² absorber ⇒ surviving-norm drop ∈ (0.5, base], base≈1 | PASS |
| ε formula | `inq-stack/include/inqkit/absorbers/mask_absorber.hpp` | formula-validation agent | ε=∫_{z<z_abs0}|ψ|²/N0 = paper Eq. 7 — CONFIRM, conditional on σ=4√2/k₀ (enforced) | LOCKED |

## Absorbing boundary — monomial CAP (inq-study, free_wp, 2026-06-15)

| Test | File | Tier | Known-case oracle | Status |
|---|---|---|---|---|
| monomial shape (absorbs) | `ResearchProject/systems/vacuum/hypotheses/cap_monomial/tests/monomial_shape_check/run.cpp` (check 1) | integration (engine, inq-study) | `absorbing_monomial` removes WP norm (ε<1, absorbed>0) ⇒ imaginary potential builds + propagates | PASS |
| monomial order monotonicity | same (check 2) | integration (engine, inq-study) | V=iη·sⁿ, s∈[0,1] ⇒ lower order absorbs more ⇒ ε(n=1)<ε(n=4); measured 0.223<0.548 (sin² hump cannot have an order signature) | PASS |
| in-header unit | `inq-study/src/perturbations/absorbing_monomial.hpp` | unit (engine, ctest) | construct + `has_potential()` true, not a uniform E-field | PASS (build) |

Notes: monomial CAP is an inq-study-ONLY new perturbation (inq/ pristine — diff
shows only the new file). ε results from it stay PROVISIONAL until Task #7 (the
inq-study engine ctest validating the scalar-potential complexification). Source:
De Giovannini–Larsen–Rubio 2014 §IV; Riss & Meyer 1996.

Notes: formula-validation agent FLAG→CONFIRM (σ must be pinned to 4√2/k₀; config
enforces it). test-validation agent caught the `-1` WPRealSpaceStats sentinel bug
and the exactly-on-peak ε=0.5 oracle bug; both fixed before lock.

## Absorbing boundary — two-sided CAP vs mask (2026-06-16)

| Test | File | Tier | Known-case oracle | Status |
|---|---|---|---|---|
| two-sided mask shape | `inq-stack/tests/include/inqkit/absorbers/test_mask_shape.cpp` | unit (pure, ctest) | `sin2_mask_value_twosided`: symmetric M(−z)=M(z), M=1 inner, M=0 both walls, M=½ at per-end midpoint, equals single-sided ramp at \|z\| | PASS |
| two-sided CAP mechanism | `ResearchProject/systems/vacuum/hypotheses/twosided_cap_vs_mask/tests/mechanism_check.py` (cap) | integration (engine, inq-study) | two summed `absorbing` slabs (`perturbations::sum`) in H absorb >90% at anchor (E≈10, L=20): measured absorbed=0.992, ε=0.78% | PASS |
| two-sided mask mechanism | same (mask) | integration (engine) | `TwoSidedMaskAbsorber` in callback absorbs >90% at anchor | (see mechanism_check) |
| kinematic LEED diffraction | `inq-stack/tests/python/inqview/analysis/test_diffraction.py` | unit (pure, pytest) | `diffraction_pattern`: cosine grating period λ → peak at \|k\|=2π/λ; DC removed → specular≈0; kx spacing=2π/(nx·dx); real input → centro-symmetric (odd N) | PASS (5/5) |
| planar Δn(z,t) reducer | `inq-stack/tests/python/inqview/analysis/test_planar_density.py` | unit (pure, pytest) | `planar_delta_map`/`planar_profile`: Σ_xy keeps z; cell_area scales; t0 column =0; known +0.5/+1.0 slab increment tracked; bad shapes raise | PASS (5/5) |
| canonical VTI loader (index↔coord) | `inq-stack/tests/python/inqview/visualisation/test_field_io.py` | unit (vtk, pytest) | `load_vti`: physical order (x[0]=left edge), cell-centred axes, data round-trips with NO fftshift, asymmetric feature stays put; `expect_centered_axis` passes centred slab and FIRES on edge-split | PASS (3/3) |
| FFT-pipeline panel (stage fidelity) | `inq-stack/tests/python/inqview/visualisation/test_fourier_panel.py` | unit (mpl, pytest) | `fft_pipeline_panel`/`fft_stages`: stage-6 amplitude == `FourierTransform.transform()` (no re-derivation), detrended==`_apply_subtract`, windowed==detrended·window, pad len==zero_pad·n, 6 axes, 5 eV tone peak recovered (angular ℏω=2πf convention), transient skip respected; **default subtract='mean'** + FFT axes overlay a **detrend comparison** line == an independent detrend transform (user verdict 2026-06-25) | PASS (4/4) |
| FFT default baseline = mean | `inq-stack/tests/python/inqview/analysis/test_fourier.py` | unit (pure, pytest) | `FourierTransform()` default `subtract=='mean'` (verdict 2026-06-25); legacy `detrend=True/False`→`detrend`/`none`, explicit `subtract=` wins; bare-default transform of (50+tone) strips the offset → peak at f0 not DC | PASS (19/19 file) |
| loss-function locator (BUG-A/B fix) | `inq-stack/tests/python/inqview/pipeline/test_density_fourier_loss.py` | unit (pure, pytest) | `density_fourier.loss_locator`: complex phasor e^{+iω0t} peaks at ω0 in +freq half (BUG-A: `.real` folding quarters the \|·\|² peak, ratio 0.25); `\|n_q\|²/q²`∝1/q² across modes (BUG-B). Verdict 2026-06-25; peak-LOCATOR not −Im[1/ε] | PASS (3/3) |
| per-run tube generator (geometry + S) *(cylindrical-jellium)* | `ResearchProject/systems/cylindrical_jellium/hypotheses/annular_sv/per_run.py` (validated in-session 2026-06-30) | sanity (real run data) | `per_run.py`: GS density peaks at \|x\|=8.8 Bohr (centre of wall band [5,13]), bore hollow → load_vti physical order, NO centre↔edge swap; `stopping_analysis` uses the `stopping-power-extraction` skill kernels (Method A, ΔE_total primary + −dKE_ion cross-check + N-guard, early v≥0.85·v0 window); energy conservation verified ΔE_total +0.045 ≈ −ΔKE_ion +0.044 Ha; ion overlay = cyan marker present in classical+WP GIFs; report 0 errors, 180/180 images resolve | PASS |
| stopping skill kernels (`stopping_power._selftest`) | `.claude/skills/stopping-power-extraction/stopping_power.py` | unit (pure, self-test) | known-slope recovery (incl. `fixed_time_fraction`), slab converged vs not-converged, N-drainage guard — all assertions pass | PASS |
| per-run FFT via fourier skill *(cylindrical-jellium)* | `per_run.fft_panels` + `inqview.visualisation.fourier_panel.fft_pipeline_panel` | sanity (real run data) | every FFT-driven observable (current_z, energy_total) rendered as the audited 6-stage panel (6 axes confirmed, plasmon band ħω_p=√(3/r_s³)·27.211 eV shaded); pipeline raw `fft_*`/`spectra/` figs EXCLUDED from notebook (0 leftover); same exclusion + audited-panel coverage added to `run-notebook` builder for all future notebooks | PASS |
| r_s=2 stopping FLAGGED (open verdict) *(cylindrical-jellium)* | per-run skill analysis, manifest flags | sanity (real run data) | dense short-cell L_z=10 runs: ΔE_total vs −dKE_ion channels diverge (rs2_v0p15 ratio 2.04 r²=0.23; rs2_v0p30 1.18/0.69; rs2_v0p45 1.11) — surfaced, NOT averaged; user owns the verdict (accept ke_ion / accept ΔE_total / rerun longer L_z) | FLAGGED |

Notes: NO new engine code — two-sided CAP is composition of the existing
`inq-study` `absorbing` via `perturbations::sum`; two-sided mask is wrapper-only
(`inq-stack`). `inq/` untouched by this task (but carries a PRE-EXISTING unrelated
`viewables.hpp` `ham()` edit, flagged to the user — not mine). CAP ε PROVISIONAL
until Task #7. Plan: `docs/plans/twosided-cap-vs-mask.md`.

---

# ml-patterns campaign-local kernels (Python, T1 pre-gate, 2026-07-01)

Campaign-local (`docs/campaigns/ml-patterns/kernels/`), NOT promoted to inqview.
Tests: `docs/campaigns/ml-patterns/tests/test_kernels.py` (numpy-only, runs in
seconds). Each formula-bearing kernel also passed an INDEPENDENT
`formula-validation` agent (given only the formula + source, not the test).

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| POD rank-2 recovery | `kernels/pod.py::pod` | unit (pure) | a field built from exactly 2 spatial modes → first 2 POD modes capture >99.9% energy, mode-3 negligible, recovered 2D subspace spans the planted modes (Eckart-Young). formula-validation: CONFIRM (Brunton & Kutz; Halko-Martinsson-Tropp randomized SVD) | PASS |
| POD randomized ≈ deterministic | `kernels/pod.py::_randomized_svd` | unit (pure) | randomized SVD + power iterations matches deterministic leading singular values < 5% | PASS |
| DMD damped sinusoid | `kernels/dmd.py::dmd` | unit (pure) | x(t)=exp(−γt)[cos(ωt)φ1+sin(ωt)φ2] → exact-DMD recovers ω and −γ to <1%. formula-validation: CONFIRM (Tu et al. 2014 exact DMD; Schmid 2010) | PASS |
| DMD windowed | `kernels/dmd.py::dmd` | unit (pure) | frequency recovered to <2% over a sub-window (non-stationary guard) | PASS |
| F_WP analytic | `kernels/formfactor.py::F_WP` | unit (pure) | exp(−q²σ_pot²/2) exact (Jackson form factor) | PASS |
| radial spectrum Gaussian width | `kernels/formfactor.py::radial_power_spectrum` | unit (pure) | FFT magnitude of a real-space Gaussian (std s) recovers s to <5%; no fftshift (physical order) | PASS |
| q_ratio of two Gaussians | `kernels/formfactor.py::q_ratio` | unit (pure) | R(q) follows exp(−q²(s1²−s2²)/2), median rel < 10% | PASS |
| F_ONCV from UPF | `kernels/formfactor.py::F_ONCV_from_upf` | sanity (real UPF) | F_ONCV(q)≈1 at low q (Coulomb-tail-subtracted radial FT of the actual `electron-ONCV-1.2.upf` local potential); |F−1|<5% for q≤1.9 1/Bohr establishes the q-range where the T2 prediction reduces to exp(−q²σ_pot²/2). formula-validation: CONFIRM (Jackson radial FT) | PASS |
| PDE-FIND advection | `kernels/pdefind.py::discover_pde_1d` | unit (pure) | translating Gaussian u(x−ct) → recovers u_t=−c·u_x (c=1.3) to <0.15, spurious terms <0.15. STRidge (Rudy et al. 2017). formula-validation: CONFIRM | PASS |
| PDE-FIND diffusion | `kernels/pdefind.py::discover_pde_1d` | unit (pure) | decaying Fourier modes → recovers u_t=ν·u_xx (ν=0.4) to <0.08, spurious <0.1 | PASS |
| PDE-FIND wave/plasma (2nd order) | `kernels/pdefind.py::discover_pde_1d` | unit (pure) | standing modes, dispersion Ω²=c²k²+ω² → recovers u_tt=c²·u_xx−ω²·u (c²=1, −ω²=−4) with correct term labels; 2nd-diff time deriv. formula-validation: CONFIRM | PASS |
| PDE-FIND forward-predict (Wall 2) | `kernels/pdefind.py::forward_score` | unit (pure) | discovered advection PDE forward-integrates (RK4) to held-out later frames, rel-L2 < 0.2 | PASS |
| PDE-FIND bootstrap (Wall 3) | `kernels/pdefind.py::bootstrap_stability` | unit (pure) | true u_xx term active in >80% of resampled-subset refits; spurious terms not | PASS |
| PDE-FIND noise robustness | `kernels/pdefind.py::discover_pde_1d` | unit (pure) | diffusion recovered under 1% additive noise with Gaussian time+space smoothing, ν to <0.15 | PASS |
| form-factor exact recovery | `kernels/formfactor_residual.py::residual_test` | unit (pure) | n_WP=F(q)·n_cl exactly → σ_fit≈σ_WP (matches_sigma_wp), t_flatness<0.05, r²>0.98, no high-q excess. Linear-response null n_ind=χ·V_ext (Lindhard/RPA) | PASS |
| form-factor high-q excess | `kernels/formfactor_residual.py::residual_test` | unit (pure) | injected high-q excess (nonlinear/quantum fingerprint) → highq_excess_over_noise >10σ vs clean baseline | PASS |
| form-factor t-drift flag | `kernels/formfactor_residual.py::residual_test` | unit (pure) | ratio with a t-drift (trajectory/deceleration mismatch) → t_flatness >5× the flat case (>0.1) | PASS |
| Fork-A collapse selects σ_WP | `kernels/formfactor_residual.py::collapse_fork_a` | unit (pure) | a(σ) with a=0.5σ² → slope 0.5, selects σ_WP (√2-trap resolved empirically) | PASS |
| Fork-A collapse selects σ_pot | `kernels/formfactor_residual.py::collapse_fork_a` | unit (pure) | filter built at σ_pot=σ/√2 → a=0.25σ², slope 0.25, selects σ_pot | PASS |
| Gaussian exponent fit | `kernels/formfactor_residual.py::fit_gaussian_exponent` | unit (pure) | clean exp(−a q²) → a recovered to <1e-6, r²>0.999 (weighted log-linear LS) | PASS |
| radial spectrum localizes k | `kernels/formfactor_residual.py::radial_spectrum` | unit (pure) | single cosine at wavevector k → power peaks in the abs(q)=abs(k) shell (within one bin); full FFT, no fftshift | PASS |

# Muon per-state mass fork (inq-study engine, 2026-07-06)

Fork of the INQ replica `inq-study` adding a tunable per-state mass (muon /
band-structure). Design: `docs/campaigns/muon_projectile/inq_study_engine_notes.md`;
plan + full validation matrix: `docs/plans/muon-mass-fork-implementation.md`.
Test file `inq-study/tests/muon_mass_fork.cpp`, built in `inq-study/build-cpu`,
run via `ctest -R muon_mass_fork`.

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| per-state kinetic factor (T1.1) | `tests/muon_mass_fork.cpp` | inq-study engine (CPU so far) | plane wave e^{ik·r} → `laplacian_states` gives `factor[ist]·(−|k|²)·ψ`, factor=−0.5·inverse_mass (3 electrons −0.5, 1 muon −0.5/206.77); diff < 1e-8 | PASS |
| fourier add variant (ks_ham :235) | `tests/muon_mass_fork.cpp` | inq-study engine | `laplacian_add_states` accumulates a 2nd kinetic term → 2× the analytic value; diff < 5e-7 | PASS |
| electrons unchanged (T1.3) | `tests/muon_mass_fork.cpp` | inq-study engine | electron states' result bit-identical whether or not the muon slot is present (no leakage) | PASS |
| expectation-value + ledger (T1.2/1.5) | ks_hamiltonian orbital_set test | inq-study engine | ⟨T⟩=k²/2m per state; apply and expectation factors identical | TODO |
| bit-for-bit inert-when-off (T0) | electron GS/RT vs pristine inq | inq-study engine | all-mass-1 → energies identical to unforked engine | TODO |
| GPU build + kernels | (nvcc build) | inq-study engine (GPU) | `_states` kernels compile + pass under nvcc | TODO |
| WavePacket `focus_z()` chirp | `scripts/muon_mass_fork/vacuum_focus/run.cpp` | inq-study engine (GPU) | converging launch focuses a σ=1 packet to its waist (density std σ_WP/√2 = 0.707) at the focal point: at **dx=0.333**, σ_z 0.864→**min 0.7071 at t=1.470** (ideal 0.7071 @ τ=1.475); transverse (unchirped) spreads monotonically. dx=0.40 under-focuses (0.775) — chirp's k0+3σ_p tail at Nyquist. | PASS (dx≤0.333) |
| density-GIF linear\|log 2-panel | `inqview/visualisation/density_gifs.py::_save_gif` | inqview (real-data smoke) | every battery GIF renders LINEAR\|LOG side by side — density: linear + `LogNorm`; Δn (delta0/dstep): linear + `SymLogNorm` (linthresh=vmax/100). Both branches produce valid multi-frame GIFs on synthetic + effmass_sigma1 stacks (269 KB / 273 KB). | PASS |
| run-notebook E_total \| N(t) panel | `run_notebook_builder.py::energy_and_number_fig` | inqview (real-data smoke) | side-by-side E_total(t) and N(t)=∫n dV; N from `electron_number.csv` (classical) else ∫`density_total` VTIs (WP). effmass_sigma1: WP **53.00→52.00** (1 e⁻ absorbed), classical **52.00→52.00** (conserved). | PASS |

# Periodicity-3 full-component mirror + quantum self-energy (campaign_autorun, 2026-07-09)

Full-component periodicity-3 remake of the H0 insertion sweep so individual energy
components are offset-free (G=0→0), for the WP-vs-classical self-energy analysis.
Runs `runs/h0_p3/{wp,cl}_r{4..40}_p3` (dispatcher `rerun_h0_p3.py`); shown in Part III
of `hypotheses/campaign_autorun_study/theoretical_slab_model.ipynb`. Handover:
`docs/handovers/campaign-autorun-review-organisation.md` (2026-07-09).

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| p3 decomposition exact | `runs/h0_p3/*/observables.csv` | sim provenance | sum(kinetic+external+hartree+xc+nonlocal+ion+ion_kinetic+exact_exchange) == total for all 12 runs; max |Δ| = 9.9e-13 Ha | PASS |
| p3 physical component signs | notebook Part III P2 waterfall | analysis | G=0→0 gives E_hartree>0 (e–e repulsion), E_external<0 (e–background well); opposite sign to p2's 0.5·rc²-offset components | PASS |
| zero-point KE = 3/(4σ²) | notebook Part III P3/P4 | analysis | d(kinetic)=E_kin(WP)−E_kin(CL) = +3.004 Ha (81.7 eV) at every r == analytic 3/(4·0.5²)=3.000 Ha; 0.1% | PASS |
| self-XC r-independent | notebook Part III P3/P4 | analysis | d(xc) = −0.605 Ha (−16.5 eV), flat across r=4..40 (σ<0.005 Ha) | PASS |
| WP self-Hartree analytic↔numeric | notebook Part III P4 | analysis | closed form 1/(2σ_ρ√π)=0.798 Ha (σ_ρ=0.354) vs isolated FFT-Poisson on saved n_WP = 0.774 Ha; agree ~3% | PASS |
| charged-cell convention caveat | notebook Part III P4 | analysis | raw dHartree(r=40) = −29 eV (p3) vs −274 eV (p2) vs physical +22 eV — convention-dependent, matches neither (net −1 charged cell) | PASS (documented) |

# Energy book-keeping campaign — E_proj_bg + exact d(H+E) decomposition (2026-07-11)

Campaign `docs/campaigns/localised_jellium_parameter_study_2/`; validation note
`docs/validation/e-proj-bg-dual-route.md`; engine analysis
`hypotheses/campaign_autorun_study/b1_decomposition.py`; results notebook
`hypotheses/campaign_autorun_study/energy_book_keeping_campaign.ipynb`.

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| E_proj_bg dual-route | `b1_decomposition.py` + notebook §B1 | analysis | closed-form periodic mean-zero vs independent 3D-FFT grid solve, σ_pot entered separately; max \|A−B\| = 0.20 eV on 80 eV scale (Lz=120); ≤0.23 eV re-validated at Lz=90 | PASS |
| E_proj_bg limiting cases | notebook §B1 | analysis | point-charge limit (σ→0.05) at z_p=−50: 48.97 vs 48.97 eV; slab-centre closed form | PASS |
| B1 known-case gate (exact t=0 identity) | `b1_decomposition.py` | analysis | E_wb+E_selfH+E_bgw−E_ghb reproduces measured d(H+E) of h0_p2 ledger at r={4,12,28,40} to ±4 eV on 40–170 eV terms | PASS (±4 eV documented) |
| ghost UPF parsed by data | notebook §B1 | provenance | PP_LOCAL = +erf(r/0.5)/r Ha exactly on mesh [0,50] Bohr; pure +1/r tail to mesh end; z_valence=0 | PASS |
| ghost truncation ablation | notebook §B1 | analysis | ONLY truncated-at-50-Bohr WITH lateral images reproduces measured energies (gap +3.4 eV); images=0 fails by +96 eV; untruncated by −200…−510 eV | PASS (decisive) |
| gs_ghost SCF far control | `gs_ghost/runs/ghost_r28_p2` | sim | r=28 screening response must vanish toward far field (B2 built-in control): SCF gain −0.4 eV at r=28 vs −2.0 (r=12), −3.0 (r=4) — monotone toward 0 | PASS |
- `projectile_background_energy` (inqkit/jellium) | E_proj_bg = ∫n_proj·v_bg (ideal) + −∫n₊·v_ion (impl) | CPU smoke: n_proj_norm=1.0003, ideal==∫n₊·poisson(ρ_proj) reciprocity 0.0000 eV; formula-validation CONFIRM (Jackson §1.11, gauge cancels in WP−CL for matched single-e pairs) | smoke_eprojbg/run.cpp | 2026-07-12
| absorbing_wrap CAP profile | `hypotheses/cap_fix/tests/test_wrap_profile.py` | unit | wrap cos² bump: peak |η| exactly at periodic boundary, smooth across wrap (slope gap 0), ∫W dz = 15η == two-sided integral, same |z|>25 footprint; two-sided sin² confirmed W=0 at boundary (the topology gap) | PASS |
| cap_fix binary smoke (wrap mode) | `scripts/cap_fix/results/smoke_wrap` | sim | 10-step wrap-mode RT run on GPU 1: builds, propagates (E monotone at floor), charge.csv N=53.000 (52 slab + 1 WP; stale pre-propagator row = 52 documents WP-refresh timing), run_completed=true | PASS |
| twin_decompose engine | `.claude/skills/twin-run-analysis/tests/test_twin_decompose.py` | unit+golden | classical-vs-WP decomposition: synthetic fixture reproduces documented table; on-disk golden pair (proj_perturbation σ=0.5 r=12) → dKin 81.74, dXC −16.47, residual 20.81, SIE 4.34 eV; parity catches mismatch; SIE=R+dXC identity; at-rest drift≈0 | PASS (10/10; +3 dynamic Rung-2) |
| check_twin gate | `.claude/skills/twin-run-generation/check_twin.py` | validator | twin parity + full-decomposition-columns + U_proj_bg availability; PASS on golden pair, FAIL on parity-broken/identical-projectile/missing-U_proj_bg; writes twin_manifest.json | PASS |
| Projectile Ehrenfest integrator | `inq-stack/tests/include/inqkit/dynamics/test_projectile.cpp` | unit (pure) | velocity-Verlet: zero-force→const V, constant-force→V=at/R=½at² exact, a=F/m scaling, harmonic energy bounded (symplectic) | PASS (4 cases/125 assertions) |
| P-dyn dynamic twin gates | `hypotheses/twin_dynamics/pdyn_k1_study.ipynb` | sim+analysis | HF-force energy conservation 0.0003 eV (classical); t=0 collapses to golden 20.81/81.7/−16.47/4.34; residual collapses 20.81→1.03 eV as WP disperses (σ_z 0.35→3.66 = analytic free-dispersion) | PASS |
