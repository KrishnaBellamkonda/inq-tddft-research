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
| `T0.7 filled cylinder carries full n0 ON the axis (w>0)` *(proximity-ladder)* | `cylinder_mask` at d = 0…R with w = 0.5 | mask = 1.0 on the axis (**not** 0.5), uniform to 1e-6 across the interior, ½ on the outer edge, <1e-6 outside. Guards the `annulus(R_in=0, w>0)` trap: the erfc step is centred ON its edge, so `background_mask(0,0,w) = ½` would put n₊ = n₀/2 exactly where the projectile flies — PASS |
| `T0.8 hollow tube still carves its bore (w>0)` *(proximity-ladder)* | `annulus_mask` at the production geometry R_in=10, R_out=14, w=0.5 | axis <1e-12, ½ on the bore edge, 1.0 mid-wall, ½ on the outer edge — the new `cylinder` shape does not leak into the hollow path — PASS |
| `T0.9 cylinder integrates to n0·πR²·Lz (builder branch)` *(proximity-ladder)* | `make_localised_background` with `background_shape::cylinder` (R=7, L=20) | ∫n₊ = n₀πR²L_z = 30.79 (ε=0.08, one curved surface); and result is independent of `inner_radius`, which the shape ignores — PASS |
| `f_bore(0) gate tolerance is grid-aware` *(proximity-ladder)* | `rayleigh_tol_pc` in `proximity_ladder/wp/run.cpp`, checked standalone against the four rungs' MEASURED t=0 deviations (2026-08-03 smoke stage) | tol = 1 % floor ⊕ `100·(df/dR)/f·(dx/2)` gives 1.00/1.10/2.42/7.34 % at r10/r08/r06/r04; measured 0.00/0.087/0.919/1.799 % all pass with 12.7×/2.6×/4.1× margin. **Still rejects a real defect**: a packet mis-sized by 25 % in σ deviates 25.2 % vs 7.34 % tol → CAUGHT (3.4× margin). Replaces a fixed ±1 % that was 120× more demanding at r04 than r10 and had false-aborted a correct run — PASS |

> Annular `background_shape` (cylindrical-jellium campaign): the erfc edge profile,
> ½-height crossovers at both radial edges, and axial uniformity are proven
> ANALYTICALLY by the `formula-validation` agent (VERDICT: CONFIRM, 2026-06-28);
> T0.4–T0.6 above validate the grid builder via decomposition-safe integrals.
>
> `background_shape::cylinder` (proximity-ladder campaign, 2026-08-02): a FILLED
> tube is a separate shape, never `annulus` with R_in = 0. **T0.5 could not have
> caught this** — it probes R_in = 0 at w = 0, where `background_mask(d,0,0) = 0`
> for all physical d ≥ 0 and the composition is accidentally right. The defect
> existed only in the SOFTENED branch, i.e. only in production (w = 0.5), and was
> maximal exactly on the tube axis where a channeling projectile flies. T0.7–T0.9
> close that gap. Ladder ground states additionally gate on it at runtime (gate 3b:
> background n₊/n₀ on the axis must read ~1 filled / ~0 hollow).

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

## Wavepacket orthogonalisation loss (near-launch campaign, 2026-08-01)

New `InjectionReport` fields (`norm_pre_ortho`, `norm_pre_renorm`,
`removed_weight`, `sum_overlap_sq`, `overlap_by_state`,
`ortho_closure_residual()`). They exist because `norm_after` is measured AFTER
renormalisation and is ≈1 by construction, so it cannot express how much of the
packet the Gram–Schmidt projection removed — which is the whole question when a
wavepacket is launched inside a slab's electronic spill-out.
Plan: `docs/plans/effective-sigma-near-launch.md`.

| Test | File | Tier | Known-case oracle | Status |
|---|---|---|---|---|
| `removed_weight` at k₀=0 | `inq-stack/tests/include/inqkit/wavepacket/test_wp_ortho_loss_engine.cpp` | engine | **analytic** 8π^{3/2}σ³e^{−σ²k₀²}/V for a constant occupied state = 0.044546624; got 0.0445464505 (4e-6 rel, tol 2 %) | PASS |
| `removed_weight` at k₀=1 | same | engine | same closed form × e^{−4} = 8.158999e-4; got 8.159026e-4 (3e-6 rel, tol 5 %) — pins the e^{−σ²k₀²} suppression | PASS |
| ortho closure identity | same | engine | Σᵢ\|⟨ψᵢ\|ψ_wp⟩\|² == ‖ψ‖²_pre − ‖ψ‖²_post (KS states mutually orthonormal); residual **0.0 bit-exact**, tol 1e-10·lhs | PASS |
| `max_overlap` consistency | same | engine | with one state below the WP slot, max_overlap == √(sum_overlap_sq): 0.2110603006 both routes | PASS |
| back-compat of `norm_after` | same | engine | ortho ⇒ norm_after==1 (renormalised) while removed_weight>0.01; no-ortho ⇒ removed_weight==0 and norm_pre_renorm==norm_pre_ortho | PASS |

18 assertions / 4 test cases, all PASS (job 32528040, A100).
State 0 is overwritten with the constant 1/√V by hand, so the oracle depends on
no SCF result — it tests inqkit's bookkeeping only.

**Live regression on real data** (job 32528019): the `inject_scan` program
re-injected the far-launch packet at z = −24 and reproduced the completed
campaign's recorded `max_overlap = 3.691564855e-4`
(`wp/results/v2p0/raw/observables/wp_config.txt`) to 12 significant figures.

| Python-side | File | Tier | Known-case oracle | Status |
|---|---|---|---|---|
| `kz_marginal` + `gaussian_fit_quality` | `inq-stack/python/inqview/visualisation/field_io.py` | pure (numpy, on committed run data) | far-launch t=0 packet is an undeformed Gaussian: recovered ⟨k_z⟩ = 1.999791 vs 2.0 (−0.01 %), σ_kz = 1.414473 vs 1/(√2σ)=1.414214 (+0.018 %), skew −0.005, excess kurtosis +0.022, R² vs **analytic** = 1.000000 | PASS |

---

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
| OrbitalKineticStats (extensive kinetic) | `systems/vacuum/hypotheses/cap_norm_investigation/extensive_kinetic/` | sim+identity | vacuum double-sided-CAP pair (dcap_extkin/dcap_baseline, 700 steps): per-step identity Σocc·T_i/norm_i == energies.csv:kinetic EXACT (0.0 Ha, all 701 steps); t=0 kin_bare == analytic ½k₀²+3/(4σ₀²) = 14.777 Ha; E_corr = total−kin+kin_bare tracks E0·norm to 0.00 eV at norm 3.5e-6; post-hoc e_kin_ha·norm agrees to 2.5e-9 eV; cost 0.42 ms/step (0.14%, 1 orbital) | PASS |
| Gaussian projectile self-energy (direct ledger) | `inq-stack/tests/include/inqkit/jellium/test_projectile_self_energy.cpp` | unit (pure) | e_pp = 1/(2σ√π) analytic self-Coulomb == numeric radial quadrature to 1e-4 (σ=0.354/0.5/1.0); σ_pot=0.354 → 0.798 Ha. Guards the compute_coulomb_direct constant | PASS (4 assertions) |
| Direct-potential projectile (no kink / no sheet) | `systems/localised_jellium/hypotheses/classical_highdensity_sv/dyn_direct/run_v4p5_direct.ipynb` | sim+analysis | v=4.5 replica, direct erf/r perturbation+force+ledger vs old charge run: (T3) all curvature maxima IN-SLAB (z≈−8) not at wall (old z=42.2); (T2) e_ps positive & monotone→0 (5.8→3.1→0.59 Ha) vs old sheet drift to −15; e_pp const 0.7979 Ha (std 2e-16); (T5) conservation drift 3.3e-3 eV, across-wall std 1.3e-5 eV; (T4/finding) sheet inflated in-slab KE-loss OLD 5.27 vs NEW 3.97 eV (v=2 pilot pair 23.1 vs 19.2), so old S(v) over-estimated ~20–35%; S(v=4.5)_direct=0.18 (old 0.28) | PASS |

---

# Annular-tube channeling twin — KS-orbital stopping in the bore (2026-08-01)

Plan `docs/plans/cylindrical-channeling-ks-stopping.md`; handover
`docs/handovers/cylindrical-channeling-ks-stopping.md`; engine
`ResearchProject/systems/cylindrical_jellium/hypotheses/channeling_twin/channeling_stopping.py`.

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| `radial_occupancy` shells + moments | `inq-stack/tests/include/inqkit/observables/test_radial_occupancy_engine.cpp` | engine (C++) | RAYLEIGH law for an on-axis Gaussian: f_bore = 1−exp(−R²/2σ_d²), ⟨r⊥⟩ = σ_d√(π/2), ⟨r⊥²⟩ = 2σ_d²; exact 3-shell partition to 1e-12; minimum-image case (tail wraps a transverse face → 0.98, non-periodic impl gives ~0.6); off-axis ⟨r⊥²⟩ = μ²+2σ_d² exact | WRITTEN, not yet compiled |
| minimum-image Hellmann–Feynman force | `inq-stack/tests/include/inqkit/dynamics/test_projectile_force_minimum_image_engine.cpp` | engine (C++) | closed form for φ_drag = cos(2πz/L): E_R = e^{−k²σ²/2}cos(kZ), F_z = e^{−k²σ²/2}k sin(kZ). Deep inside the cell both kernels agree and hit it; straddling the face ONLY the min-image one does (the clipped one is asserted to deviate >5 %); axis-general form == the z wrapper; transverse force = 0 | WRITTEN, not yet compiled |
| channeling stopping engine | `.../hypotheses/channeling_twin/tests/test_channeling_stopping.py` | analysis (unit) | constant-S closed form (dp/dt = −S ⇒ p linear ⇒ trapezoid s₄ EXACT): classical S and all four S_ij recover the input to **1e-8**; s₃ (circular, unwrapped) == s₄ across a face crossing; window follows the MEASURED f_bore (breach at t=10 → window ends at 10, beating the 23.3 formula); freeze detection separates frozen from +50 % growth; (T1−T2)(0) == 3/(4σ²) = 1.2755 eV; all three verdict branches incl. "clean channeling but S still differs → look at E_PP"; resume-segment concat de-duplicates the boundary step | **PASS (8/8)** |
| comparison-notebook cell guards | `.../hypotheses/channeling_twin/tests/test_comparison_notebook_cells.py` | static | every generated code cell compiles; NO non-raw string literal containing a backslash (a valid Python escape silently eats the LaTeX command — `"$\approx$"` → BEL, cost one full notebook execution to find); no mathtext command matplotlib lacks (`\le` vs `\leq`); balanced `$` spans; the RESULT/PREMISE/MECHANISM figure blocks and the mandatory density GIF are still emitted and DISPLAYED | **PASS (5/5)** |
| comparison notebook end-to-end | `.../hypotheses/channeling_twin/build_comparison_notebook.py` | build smoke | executed against a synthetic twin encoding S = 0.20 eV/Bohr: 23 cells, **0 error outputs**, 9 inline figures, verdict returned `AIM MET`, both halves recovered S = 0.200 to 1e-8 (synthetic artefacts deleted afterwards) | **PASS** |
| cutoff/aliasing guard (pre-run, mandatory) | `.claude/skills/tddft-simulations/cutoff_guard.py` | validator | WP (σ_WP = 4, 50 eV, dx = 0.5): aliased tail **0.00 %** (σ_p = 0.1768, k_Nyq = 6.2832); classical: E_cut = 537 eV ≥ 1.10 × 50 eV | **PASS (both halves)** |
| GS gates (in-binary) | `.../scripts/channeling_twin/gs/run.cpp` | sim gate | ∫n dV = N exactly (neutrality; G=0 cancellation requires it); num_states = N/2 + extra (the RT binaries must load the same system); **bore depletion n̄_bore/n̄_wall < 0.5** (if the bore is not electron-poor there is no channel). No E_GS reproduction gate — this is a NEW system with no published reference, and inventing one would be worse than omitting it | NOT YET RUN |
| WP t=0 analytic gates (in-binary) | `.../scripts/channeling_twin/wp/run.cpp` | sim gate | norm = 1; ⟨p_z⟩ = k₀ (1 %); σ_pz² = 1/(2σ²) (5 %); T1 = (k₀²+3σ_p²)/2 (2 %); T1−T2 = 3/(4σ²) (5 %); **CIRCULAR** centroid = launch_z ±0.05 and circular spread = σ_WP/√2 (5 %) — the naive ⟨z⟩ is NOT gated because the packet straddles the face at t=0 by construction; f_bore(0) vs the Rayleigh value (1 %); ⟨r⊥⟩(0) = σ_d√(π/2) (5 %) | NOT YET RUN |
| twin parity gate | `.claude/skills/twin-run-generation/check_twin.py --dynamic` | validator | periodicity/Lz/spacing/N/sigma_wp/launch_z/gs_dir identical; `projectile` field DIFFERS; full energy-decomposition columns in both; projectile.csv carries proj_z/energy_proj_ke/energy_proj_bg_ideal | NOT YET RUN |
| run correctness gates (in-binary) | both `run.cpp` | sim gate | WP: no CAP ⇒ H Hermitian and t-independent ⇒ `energy_total` drift < 1e-3 eV. Classical: `E_electronic + KE_proj + U_proj_bg` drift < 0.05 eV (a drift means force and perturbation disagree). Both: pairwise ledger closes to INQ's own `energy_hartree`/`energy_external` | NOT YET RUN |

---

# Bulk-jellium KS-orbital twin — interaction energies + phase space (2026-08-01)

Plan `docs/plans/bulk-jellium-ks-stopping.md`; handover
`docs/handovers/bulk-jellium-ks-stopping.md`; engine
`ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/ks_stopping.py`.

| Test | Path | Tier | Asserts against | Status |
|---|---|---|---|---|
| `local_stopping` rolling-OLS S(z) | `.../hypotheses/bulk_ks_stopping/tests/test_ks_stopping.py` | analysis (unit) | planted CONSTANT slope recovered to 1e-9 everywhere in the interior; planted LINEAR ramp S(s)=a+bs recovered to 1e-6 at the window centre; edges are FILLED with the nearest interior value, not extrapolated (an extrapolated edge manufactures a fake "S spikes at impact"); all-NaN when fewer than 2·half_width+1 samples | **PASS (4/4)** |
| `Interactions` clipping detector | `.../hypotheses/bulk_ks_stopping/tests/test_ks_stopping.py` | analysis (unit) | onset = start of the **trailing contiguous** clipped run, NOT the global min — launch rows sit ~4e-9 under 1.0 from discretising the Gaussian on the grid, and a global min mislabels t=0 as the onset (the bug actually made on 2026-08-01); `clip_time` = inf when never clipped; a WP half never reports clipping whatever its norm does; `in_window()` truncates a requested window back to the onset | **PASS (5/5)** |
| `PairPhase.divergence` | `.../hypotheses/bulk_ks_stopping/tests/test_ks_stopping.py` | analysis (unit) | two velocity histories separating linearly cross 5 % of v₀ at the analytically known time (±0.02); identical histories return NaN, not a spurious index 0; the gap is ABSOLUTE — a slower wavepacket triggers at the same time as a faster one | **PASS (3/3)** |
| notebook cell guards (ALL 12 builders) | `.../hypotheses/bulk_ks_stopping/tests/test_notebook_cells.py` | static | Covers 8 run + 4 phase notebooks. Every generated code cell compiles (catches `\r` from `\rangle`/`\rm` and `\n` from an f-string annotation — both are line terminators; these broke 5 phase cells AND all 8 run-notebook array tasks, 2026-08-01, the same bug as the channeling-twin builder the same day). NO control character in a cell literal, read from the **AST** (the silent variants `\a \b \f \v` compile but render wrong). **Every `$...$` span parsed by matplotlib's own MathTextParser** — caught `\frac12mv^2`, valid Python AND valid TeX but not valid mathtext, which no blocklist would have had. Balanced `$`; interaction-energy section present in every run notebook; density GIF emitted AND `display()`ed. **6 negative self-tests prove each guard FIRES on the real defect, +1 that correct cells do not trip** | **PASS (79/79)** |
| interaction-energy closure (per run, in-notebook + standalone) | `ResearchProject/systems/jellium/scripts/verify_interactions_closure.py`, and section 7 of each run notebook | sim gate | classical `E_SS == energy_hartree`; WP `E_SS+E_PS+E_PP == energy_hartree` (vs INQ's own scalar — the internal identity is exact by construction and proves nothing alone); `E_SB = E_PB = E_BB = 0` (bulk: uniform background ⇒ φ₊ ≡ 0); classical E_PP constant on unclipped rows | **PASS, 8/8 halves** (max resid 5e-13 Ha) |
| phase-notebook analysis smoke, all 4 pairs | `<scratchpad>/validate_phase_cells.py` | build smoke | every analysis cell run headless against real data for all four families; σ-matching gate `ΔE_PP(0) < 1e-3 eV`; rigid-cloud gate; every fit window CLEAR of its clipping onset; local S cross-checks the established windowed fits (classical 0.363 vs 0.377, WP-drift 0.054 vs `S_24` 0.057) | **PASS (4/4 families, 0 failures)** |

### 100 eV high-density case study (`bulk_ks_stopping_rs4/case_study_100eV`, 2026-08-02)

Run: `venv/bin/python -m pytest ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/tests/test_case_study.py -q` → **21 passed**
(campaign suite total **123**, was 102). Builder:
`.../hypotheses/bulk_ks_stopping_rs4/case_study_100eV/make_case_study.py`.

| Test | Tier | Asserts against | Status |
|---|---|---|---|
| `Meta` provenance parsing (4 tests) | analysis (unit) | cell/N_e/σ/z₀/dt/fit-window parsed from SYNTHETIC `run_summary.txt`+`wp_config.txt` with planted values; r_s recovers an analytically planted 2.000 Bohr to 1e-4; ω_p = √(4πn); a missing run half raises rather than producing mislabelled figures. Every figure title and the results file are built from these, so a parsing slip mislabels all 16 silently | **PASS (4/4)** |
| Ehrenfest centroid integration (2 tests) | analysis (unit) | z₀+∫⟨p_z⟩dt is EXACT (1e-12) for constant p_z and (1e-10) for a linear ramp — trapezoid is exact on both, so any deviation is an indexing/dt bug, not discretisation. This is the position axis the user specifically asked for in place of the density centroid | **PASS (2/2)** |
| T₂ vs T₁ contraction (2 tests) | analysis (unit) | T₂ = ⟨p⟩²/2m uses the MEAN of all three components (planted px,py,pz); with zero spread planted, T₁−T₂ ≡ 0. Confusing ⟨p⟩² with ⟨p²⟩ would not look wrong on any plot | **PASS (2/2)** |
| figure-margin guard (3 tests) | static/layout | `verify_margins()` passes clean art, **FIRES on ink in row 0** (negative self-test — the real defect: a fixed canvas crops an overrunning title and still reports success; it hit 15 of 16 figures in the first cut), and reports a blank figure | **PASS (3/3)** |
| title auto-shrink (2 tests) | static/layout | a 200-char title is shrunk below 9 pt and reported; a short title is left untouched (positive control, so the helper cannot "fix" what is already fine) | **PASS (2/2)** |
| mathtext + control-char guards (3 tests) | static | all 34 `$...$` spans in the builder parsed by matplotlib's own `MathTextParser`; **negative self-test** that `$\frac12 mv^2$` still raises (valid Python AND valid TeX, invalid mathtext); no `\a \b \f \v` control chars in any string literal, read from the AST | **PASS (3/3)** |
| S additivity on the real run | analysis (integration) | S(T₁) = S(T₂) + S(T_var) to 1e-10 — OLS slope is linear in the ordinate, so a failure means the three fits used different windows or path coordinates, which is exactly how a plausible-but-wrong S is produced | **PASS** |
| zero-point spread + beam energy (2 tests) | physics (integration) | T₁−T₂ at t=0 equals 3/(4σ_ψ²) = **5.102 eV** at σ_ψ=2 (ψ-width convention, `.claude/rules/sigma-wp-convention.md`); **T₂(0)**, not T₁(0), equals the nominal 100 eV. Pins the correction to the 2026-08-01 handover, which quoted 2.55 eV from the density-width convention | **PASS (2/2)** |
| σ-matching validated FROM the data | physics (integration) | E_PP(0) agrees between the two halves (0.176996 Ha) because the classical UPF is generated at σ_pot = σ_WP/√2 — a check ON the convention, not an input to it; and the classical E_PP is constant to <1e-9 Ha (a rigid cloud cannot spread) | **PASS (2/2)** |

### Channeling twin — refined analysis (`hypotheses/channeling_twin`, 2026-08-02)

Run: `cd ResearchProject/systems/cylindrical_jellium/hypotheses/channeling_twin && venv/bin/python -m pytest tests/ -q` → **36 passed**.

| Test | File | Pins |
|---|---|---|
| `test_t2_minus_t1_is_exactly_the_variance_term` | `tests/test_refined.py` | `T2 - T1 == var(p)/2m` to 1e-12 — the drift/spread split the study rests on |
| `test_t2_reconstruction_matches_inq_kinetic_energy` | ” | our moment reconstruction == engine `e_kin_ha` |
| `test_label_swap_is_the_users_convention_not_the_engines` | ” | `T1` is the SMALLER (drift) branch; fails if renamed back to `ks_stopping.py`'s convention |
| `test_p_integral_path_recovers_a_known_track` | ” | cumulative trapezoid exact for the track that generated it |
| `test_classical_energy_budget_closes` | ” | `d(E_bath) + d(KE_proj) == 0` |
| `test_interaction_deltas_start_at_zero` | ” | delta columns zero at t=0 for both halves |
| `test_fit_recovers_a_planted_stopping_power` | ” | constant-deceleration fixture ⇒ `dT1/ds = a` exactly; recovers S to 1e-6 |
| `test_fit_returns_nan_not_garbage_on_an_empty_window` | ” | empty window ⇒ NaN, not a spurious slope |
| `test_momentum_slices_keeps_the_whole_k_axis` | ” | **REGRESSION**: `_concat_segments` collapsed the 128-bin k axis to 1 bin |
| `test_momentum_slices_rejects_a_collapsed_k_axis` | ” | the guard raises rather than returning a one-bin frame |
| `test_nearest_slices_snaps_and_does_not_interpolate` | ” | distributions are snapped, never blended |
| `test_production_wp_identity_and_free_reference` | ” | identities hold on real GPU data; `T2-T1(0) == 3/(4σ²)` |
| `test_production_classical_budget_closes` | ” | real integrator: closure < 1e-3 eV, ΔKE = −5.1256 eV |
| `test_production_momentum_distribution_is_a_distribution` | ” | 128 bins/step, peaked at k₀ = 1.917 |
| `test_production_unwrapped_path_starts_at_the_launch_point` | ” | **REGRESSION**: `proj_z_unwrapped` one-step lag |
| 8 static notebook guards | `tests/test_refined_notebook_cells.py` | cells parse; no FORM FEED/CR leaked from `\frac`/`\rangle` in non-raw strings; mathtext subset; balanced `$`; required sections; T1/T2 table present; window defaults to unset |

### 2-D momentum map (`inqview.visualisation.field_io.kz_kperp_map`, 2026-08-02)

Run: `venv/bin/python -m pytest inq-stack/tests/python/inqview/visualisation/test_kz_kperp_map.py -q` -> **7 passed**.

| Test | Pins |
|---|---|
| `test_map_is_normalised_and_shaped_correctly` | sum = 1, non-negative, k_z sorted |
| `test_longitudinal_marginal_recovers_the_drift` | `<k_z>` and `var(k_z)` EXACT (rel 1e-6) — k_z is never binned |
| `test_transverse_marginal_is_rayleigh_not_gaussian` | shell Jacobian present: mode at `sigma_p`, NOT 0; bias sign pinned so a lost Jacobian (which halves `<k_perp^2>`) fails |
| `test_agrees_with_the_existing_1d_kz_marginal` | cross-check vs independently-tested `kz_marginal` |
| `test_drift_moves_weight_along_kz_only` | deceleration visible on k_z, invisible on k_perp — the discrimination the function exists for |
| `test_binning_default_is_one_transverse_grid_spacing` | no bins finer than the grid supports |
| `test_rejects_an_empty_field` | raises rather than returning NaNs |

### Channeling twin refined analysis — additions (2026-08-02)

| Test | Pins |
|---|---|
| `test_kz_asymmetry_is_exactly_neutral_on_a_symmetric_distribution` | CDF interpolation; the naive `kz>mean` count returns 0.454 on an exactly symmetric packet |
| `test_kz_asymmetry_detects_a_planted_skew` | correct SIGN of skew and of `frac_above_mean` |
| `test_kz_asymmetry_accepts_a_2d_map` | marginalises a full map |
| `test_impulse_ratio_is_one_for_identical_twins` | ratio == 1 for identical deceleration |
| `test_impulse_ratio_scales_with_a_weakened_drag` | halved impulse -> 0.5, not 0.25 (guards against ratioing the ENERGY) |
| `test_combined_projectile_coupling_is_the_sum_of_the_two_terms` | `dE_PS + dE_PB`, zero at t=0 |
| `test_production_momentum_map_round_trips_to_recorded_moments` | map vs `wp_momentum_stats.csv`; catches a wrong FFT ordering |
| `test_production_combined_coupling_agrees_far_better_than_its_parts` | sum < 0.35 eV while `E_PS` alone > 2.0 eV |

### WP self-interaction correction (`inqkit::SelfInteractionCorrection`, 2026-08-02)

Plan: `docs/plans/wp-self-interaction-correction.md` (§4). Engine test runs in
the chan-tests SLURM gate; the run-level tiers are jobs in the
`submit-channeling-sic.sh` chain, each `afterok`-gating the next.

| Test | Pins |
|---|---|
| `test_wp_sic_engine` / kick semantics | real multiplicative kick: density and ⟨p_z⟩ invariant at any dt_eff (zero-force, ∫n∇v_H[n]=0); large dt_eff MUST inflate var(p_z) (phase gradient is real momentum); norm exact; no projection in vacuum |
| `test_wp_sic_engine` / Q projection | after an exaggerated leak, ⟨ψ_j\|ψ_wp⟩ ≤ 1e-10 restored, norm 1, `max_overlap_pre`>0 and `norm_removed`>0 reported, bath columns bit-untouched |
| `test_wp_sic_engine` / D1 run-consistency | `u_self == energy_hartree` and `exc_self == energy_xc` (1e-9 rel) for a 1-electron system where n_total = n_wp; polarised-PZ exchange (×2^{1/3}) asserted OUTSIDE tolerance — the spin-consistency defect the review found |
| Tier V vacuum `sic_pzrun` (run-level, HARD) | var(p_z) drift < 0.1 %; σ_dens(t_end) within 0.5 % of √(σ²/2+t²/2σ²); \|⟨p_z⟩−k0\| < 1e-4; E_corrected drift PASS<1e-5 eV / WARN<1e-3; binary exits 4 on failure → afterok blocks Tier B + production |
| Tier V vacuum `sic_h` (run-level) | zero-force + E_corrected ladder; σ under-spread EXPECTED (xc self-binding remains) — reported, discriminates the xc share of the SIE |
| Tier B jellium 200-step (run-level, HARD) | cum_norm_removed < 1e-3; max_overlap_pre < 1e-3; \|ΔE_corrected\| < 0.02 eV; first contact of the projection with a real occupied manifold (vacuum has none) |
| Production (run-level, SOFT) | same three, WARN-and-report (checkpoint-dont-block); E_total drift explicitly NOT a gate under SIC — the conserved quantity is E_corrected (plan §0/D2) |
| vacuum analysis layer, SIC arm (`vacuum/hypotheses/wp_selfinteraction/tests/test_selfinteraction.py`, 3 new, 2026-08-02) | `load_all` skips absent theories but REQUIRES the reference; header-only `sic.csv` ⇒ no SIC data (`u_self is None`); a corrected run planted ON the free solution surfaces its diagnostics, passes the reference's own `numerics_gate` (the Tier V criterion), reports ~0 % excess in `summary_table`, and the hartree/lda xc-difference row survives the extra rows — 13/13 PASS |
| report-2 self-Hartree figure layer (`docs/reports/report2/drafts/draft1/figures/self_hartree/tests/test_selfhartree_data.py`, 19 new, 2026-08-05) | The analytic one-electron self-energy kernel behind `sec:results-self-hartree`, against values known BEFORE it was written: `U = 1/(2a√π)` closed form; `|E_x|/U = 0.67840955` exactly and width-independently over four decades of `a` (derived analytically, cross-checked by radial quadrature of ∫n^{4/3} to 4e-16) — the structural claim the whole subsection rests on; `E_c` deliberately NOT ∝ 1/a (if it ever were, figure SH(b) would be a flat line); reproduction of the cylindrical study's independently measured σ_WP = 4 terms (2.7139 / −1.8412 / −0.6320 eV); and agreement with INQ's OWN `exc_self_ha` from all six vacuum `sic_pzrun` runs to 1e-9 rel (parametrised per σ) — i.e. the plotted curve IS the functional, not a model of it. Plus the run-data invariants the figures depend on: scaled-box `E_PP(0)·σ` constant to 1e-9 Ha and within 0.5 % of the analytic `[1/√(2π) − ξ/36]`; fixed-box `A` recovering the analytic Gaussian coefficient to 0.5 % with `C` in (0.3, 1.0) eV; classical `E_PP` constant (rigid cloud, bound 1e-2 eV vs the observed 4e-9); `E_PP(0)` equal across twin halves (validates σ_pot = σ_WP/√2); `a(0) = σ_WP/√2`; SIC-PZ returning the packet to free evolution (<1e-6); the fixed-time reading REFUSING to extrapolate below σ_WP = 4; and the two clocks running in opposite directions (the paragraph-3 claim, asserted on the data) — 19/19 PASS |
