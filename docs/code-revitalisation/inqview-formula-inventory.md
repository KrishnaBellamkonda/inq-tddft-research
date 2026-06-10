# inqview — formula & numerical-logic inventory

**Purpose.** Companion to `inqkit-review-and-next-steps.md`. This is the
metric/formula table requested in its "Materials required for the Python
inqkit code-review" section: every numerical/formula-bearing function in the
Python `inqview` library, the file it lives in, what it computes, and whether
a test currently exists. Use it as the checklist for deciding which unit and
integration tests to write.

**How it was produced.** Auto-extracted from the Understand-Anything knowledge
graph built over `inq-stack/python/inqview` at commit `1926e0a`. Scope was
restricted to formula/numerical/technical logic — plotting, rendering, and
notification modules (`report1/`, `plots.py`, `paraview.py`, `config.py`,
`defaults.py`, `email.py`, `postprocess/{layout,paraview_3d,compare}.py`,
`scripts/render_density_series.py`) were excluded via `.understandignore`.

> **Note on the "Quantity & formula" column.** These descriptions are
> auto-summarised from the code (docstrings + body) by the graph builder.
> They are a *navigation aid for your review*, not an independent verification
> that the code is correct. Confirming each formula against its source
> (paper/textbook) and pinning it with a known-case test is exactly the review
> work this table is meant to drive.

**Coverage at a glance.** 66 formula/convention units across 24 files. Only
`postprocess/lindhard.py` currently has a test (`test_lindhard.py`, an f-sum-rule
check). Every other row is **untested**.

## Interactive dashboard

A browsable knowledge graph of this same scope is live locally:

- **URL:** http://127.0.0.1:5174/?token=cf89ecc3beaa59bc04bb3fd19bcdc29d
- **Graph file:** `inq-stack/python/inqview/.understand-anything/knowledge-graph.json`
- Filter nodes by the `formula` / `convention` tags to see only the rows below.
- The built-in **Tour** (11 steps) walks the code in correctness-audit order:
  field/IO conventions → FFT-shift & VTI coordinates → orchestration spine →
  Lindhard → loss function → stopping chain → projection/occupation →
  energy ledger → spectra → KL/density.

## Suggested review / test priority

| Tier | Why | Files |
|---|---|---|
| **1. Conventions (everything depends on these)** | A wrong index/FFT-shift/VTI axis order silently corrupts every downstream metric (review-doc point #12). | `data.py` (`_reshape_flat_array`), `vti.py` (`_vtk_point_data_flatten`/`write_vti`), `screens.py` (`load_leed_pattern` fftshift), `fourier.py` (rfft normalisation) |
| **2. Core response formulae** | The physics results hinge on these; benchmark against analytic limits. | `lindhard.py` (extend beyond f-sum rule), `spectral_weight*.py` loss function, `density_fourier.py` |
| **3. Stopping & projection metrics** | Sign conventions & normalisations (see anti-wake sign memo). | `knudsen_ke.py`, `wake.py`, `stopping.py`, `gs_projected_occupations.py`, `gamma_transitions.py` |
| **4. Bookkeeping** | Conservation/closure checks. | `energy_balance.py`, `bath_energy.py`, `kl_divergence.py`, `wp_trajectory.py` |

## Full inventory (grouped by file, audit order)

### `data.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_parse_float_triplet` | Parses a whitespace-separated 3-value float triplet (used for origin and spacing in Bohr), enforcing exactly three components. | **no** |
| `_read_flat_real_array` | Reads a raw binary file into a flat numpy array via np.fromfile using the meta dtype, validating that the on-disk byte count matches the expected element count before returning. | **no** |
| `_reshape_flat_array` | Reshapes a flat binary array into (nx,ny,nz) per layout: 'x_slowest_z_fastest' (flat=ix*ny*nz+iy*nz+iz) reshapes C-order directly; 'x_fastest_z_slowest' (flat=iz*ny*nx+iy*nx+ix) reshapes to (nz,ny,nx) then transposes (2,1,0). | **no** |
| `infer_meta_path` | Derives the .meta sidecar path from a raw data path by stripping known data suffixes and substituting the meta extension. | **no** |
| `load_complex_field` | Loads a complex 3D field from separate real/imag binary dumps, resolving part paths from meta (or inferring the imag path by name substitution), reshaping both per layout, and building a ComplexField3D. | **no** |
| `load_meta` | Parses a .meta sidecar file into a FieldMeta, extracting nx/ny/nz, origin and spacing triplets (Bohr), layout, dtype, units and file references with strict key validation. | **no** |
| `load_real_field` | Loads a real 3D field by inferring/parsing its meta, reading the flat binary value array, reshaping per the layout convention, and wrapping it in a RealField3D. | **no** |

### `fields.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `ComplexField3D` | Dataclass holding separate real and imag 3D arrays with FieldMeta; assembles the complex array and computes magnitude = sqrt(re^2+im^2) and phase = arctan2(im, re). | **no** |
| `FieldMeta` | Grid metadata dataclass for a 3D field: shape (nx,ny,nz), num_points, numpy_dtype mapping, voxel_volume_bohr3 = dx*dy*dz, expected byte counts, and real/complex predicates. | **no** |

### `fourier.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `FourierTransform` | Configurable one-sided FFT engine for TDDFT observable time-series: windows, detrends, zero-pads, applies np.fft.rfft, computes freq = np.fft.rfftfreq(n_pad, d=dt) and amplitude = \|rfft\|/N with interior bins doubled, plus energy/current/dipole convenience transforms. | **no** |
| `transform` | Core FFT routine (FourierTransform.transform): validates and masks the time/value series, optionally detrends, applies the window, zero-pads to n*zero_pad, runs np.fft.rfft, and returns one-sided amplitude=\|rfft\|/N (interior bins x2) with freq from rfftfreq. | **no** |

### `vti.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_vtk_point_data_flatten` | Flattens an (nx,ny,nz) array into VTK point order by transposing (2,1,0) then ravelling, so x is the fastest-varying index as VTK ImageData requires. | **no** |
| `_vtk_xml_scalar_type` | Maps a numpy dtype to the corresponding VTK XML scalar type name (Float32/Float64/Int32/...), raising TypeError for unsupported dtypes. | **no** |
| `write_vti` | Writes a RealField3D as an ASCII VTK XML ImageData (.vti) file with WholeExtent 0..n-1, Origin and Spacing in Bohr; flattens the array so x varies fastest in VTK point order, written as POINT_DATA. | **no** |

### `screens.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `LeedPattern` | Dataclass for a loaded LEED diffraction screen: holds the centred 2D data plus grid geometry, exposing x_axis/y_axis (np.arange*spacing+origin), imshow extent_bohr, and an inverse_fft real-space reconstruction wrapper. | **no** |
| `load_leed_pattern` | Parses a LEED screen text dump (two header comment lines + data rows), applies np.fft.fftshift so the FFT-natural origin moves to the array centre, and overrides origin to (-Lx/2,-Ly/2) so the pattern spans [-L/2,+L/2]. | **no** |

### `postprocess/_common.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `post_ifw_window_au` | Computes the post-instantaneous-field-window time (a.u.) from launch z, box length L, WP sigma, and velocity; marks when the WP has cleared the launch region for valid stopping-power analysis. | **no** |
| `post_ifw_window_from_summary` | Parses run_summary.txt for launch_z, L, sigma and velocity (or k0) and returns the post-IFW analysis window in a.u.; the run-summary-driven wrapper around post_ifw_window_au. | **no** |
| `sigfigs` | Rounds a float to n significant figures (1/2 fig variants), returning the value with leading-digit-based decimal precision; used for compact plot-title number formatting. | **no** |

### `postprocess/_ifft.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `reconstruct_real_space` | Reconstructs a real-space image from a reciprocal-space pattern via 2D inverse FFT. Supports method='complex' (Re of ifft2 of fftshifted data), 'intensity' (sqrt amplitude proxy), and 'amplitude', with optional 2D Hann windowing. | **no** |

### `postprocess/lindhard.py` — has tests

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_F_imag` | Dimensionless imaginary part G(z,u) of the T=0 retarded Lindhard function: pi*u/2 in Region I (0<=u<=\|1-z\|, z<1), pi(1-(z-u)^2)/(8z) in Region II (\|1-z\|<u<=1+z), zero outside the electron-hole continuum. | yes (test_lindhard.py) |
| `_F_real` | Dimensionless real part F(z,u) of the Lindhard function: 0.5 + (1/8z)[(1-(z-u)^2)ln\|(z-u+1)/(z-u-1)\| + (1-(z+u)^2)ln\|(z+u+1)/(z+u-1)\|], with z=q/2kF, u=omega/(q vF); regularised logarithms. Static limit F->1. | yes (test_lindhard.py) |
| `chi0` | Complex Lindhard susceptibility chi0(q,omega)= -(kF/pi^2)(F_real + i F_imag) for the 3D FEG, broadcasting q and omega; prefactor kF/pi^2 = N(E_F), static-limit sign chi0(q->0,0) = -N(E_F). | yes (test_lindhard.py) |
| `epsilon_rpa` | RPA dielectric function epsilon(q,omega) = 1 - v(q) chi0(q,omega), with the Giuliani-Vignale chi0<0 convention so the plasmon sits at Re(epsilon)=0. | yes (test_lindhard.py) |
| `loss_function` | Energy-loss function L(q,omega) = Im[-1/epsilon_RPA(q,omega)], the spectral weight of single-particle and plasmon excitations. | yes (test_lindhard.py) |
| `plasmon_omega` | Plasmon dispersion: 'plasma' returns omega_p=sqrt(4*pi*n) with n=kF^3/(3*pi^2); 'bohm_gross' returns sqrt(omega_p^2 + (3/5) vF^2 q^2). Atomic units. | yes (test_lindhard.py) |
| `stopping_power` | Electronic stopping power S(v)=(2/(pi v^2)) integral over ln q of [integral_0^{qv} omega*Im(-1/eps) domega], unit charge in jellium. Log-spaced q grid (qmax default 2v+kF), trapezoidal omega integration, integrated over ln q to carry the 1/q weight. Returns Ha/Bohr. | yes (test_lindhard.py) |
| `vq` | 3D Coulomb interaction in Fourier space v(q)=4*pi/q^2 (atomic units), guarding q=0 by returning 0 via an infinite denominator. | yes (test_lindhard.py) |

### `postprocess/test_lindhard.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `test_f_sum_rule` | Verifies the f-sum rule integral of omega*Im[-1/eps] domega equals pi*omega_p^2/2 for the RPA loss function, a key conservation check on the dielectric response. | — (is the test) |

### `postprocess/spectral_weight.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_free_wp_density_on_grid` | Computes the analytical free-WP density n_WP(r,t)=(2*pi*s_t^2)^(-3/2) exp(-\|r-r_c(t)\|^2/(2 s_t^2)) with s_t^2=sigma^2+t^2/(4 sigma^2) on a cubic N^3 grid using minimum-image periodic distances, for subtraction of the WP's own density. | **no** |
| `_reference_curves` | Returns plasmon dispersion omega_pl=omega_p sqrt(1+3q^2/(5 kF^2)), particle-hole boundaries omega_+- = \|q^2/2 +- q kF\|, and kinematic line omega=q v0 for the (q,omega) map overlay. | **no** |
| `run` | On-axis spectral-weight driver: loads density VTI series, computes delta_n, 3D-FFTs and Hann-windowed time-FFTs to W(q_z,omega), subtracts the analytic WP density, forms chi=delta_n_resp/V_ext and loss L=-(4*pi/q^2)Im chi, and writes the raw and response (q,omega) maps. | **no** |

### `postprocess/spectral_weight_full.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_free_wp_density_on_grid` | Evaluates the analytical free-wavepacket density n_WP(r,t)=(2*pi*s_t^2)^(-3/2) exp(-\|r-r_c(t)\|^2/(2 s_t^2)) with spreading s_t^2=sigma^2+t^2/(4 sigma^2) on the cubic grid (minimum-image periodic), to be subtracted before response extraction. | **no** |
| `_reference_curves` | Computes overlay reference curves on the (q,omega) map: Bohm-Gross plasmon omega_pl=omega_p sqrt(1+3q^2/(5 kF^2)), particle-hole boundaries omega_+- = \|q^2/2 +- q kF\|, and kinematic line omega=q v0. | **no** |
| `run` | Full spectral-weight driver: loads the density VTI series, forms delta_n, applies 3D spatial FFT and Hann-windowed temporal FFT, builds the raw \|delta_n(q,omega)\|^2 weight and the WP-subtracted response/loss maps, and renders (q,omega) colourmaps with reference curves. | **no** |

### `postprocess/density_fourier.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_bohm_gross_omega` | Bohm-Gross plasmon dispersion omega(q)^2 = omega_p^2 + (3/5) vF^2 q^2 + q^4/4 (atomic units), used to overlay the expected resonance on the density-Fourier spectrum. | **no** |
| `run` | Driver: loads the VTI density series, forms delta_n(r,t)=n(r,t)-n(r,0), 3D-FFTs each frame to extract axial n_q_m(t) for m up to m_max, then time-FFTs each component, writing n_q_vs_time and n_q_spectrum CSV/PNG outputs with the Bohm-Gross overlay. | **no** |

### `postprocess/stopping.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_free_particle_sigma` | Analytic free-particle Gaussian WP width sigma(t)=sigma_0 sqrt(1+(t/(2 m sigma_0^2))^2) (m=1) used as the no-interaction baseline for WP spreading. | **no** |
| `run` | Phase entry: plots dE_kinetic(system) vs projectile z (stopping-power proxy), classical stopping force F_z vs z, and WP sigma_z(t) against the analytic free-particle spreading curve. | **no** |

### `postprocess/knudsen_ke.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_compute_stopping_power` | Computes the one-sided stopping-power slope -d(E_kin)/dz via forward differences in z (eV/Bohr, positive = WP slowing down), plus the shifted-to-zero dE_kin. | **no** |
| `_from_momentum_distribution` | Retroactive estimator: from the \|k\|-binned WP histogram computes <\|k\|^2>(t)=sum_bin k^2 n_wp / sum_bin n_wp and E_kin_ha = <\|k\|^2>/2 per step. | **no** |
| `run` | Phase entry: selects native (wp_momentum_stats) or retroactive (histogram) E_kin source, merges trajectory, computes the stopping power, and writes the knudsen KE CSV and E_kin/stopping-power plots. | **no** |

### `postprocess/wake.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `bath_line_z` | Computes the 1D bath z-profile by summing n_system over x,y and weighting by dx*dy (units e/Bohr); returns (z, line, time_au). | **no** |
| `bath_volume` | Returns the 3D bath density n_system = n_total - n_wp (or n_total for classical runs), subtracting the WP only from an EXACT-step density_wp partner (raises if require_exact_wp and none exists) to avoid a moving-WP dipole residual. Returns (n_system, origin, spacing, time_au, step). | **no** |
| `wp_centroid_z` | Computes the wavepacket centroid z(t) as the first moment sum(z*w_z)/sum(w_z) of the density_wp z-profile (clipped non-negative) at the nearest WP frame; returns None for classical runs. | **no** |

### `postprocess/gamma_transitions.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_build_transitions` | Enumerates all occupied->unoccupied KS single-particle pairs at a k-point (occ>occ_threshold, occ<unocc_threshold) and tabulates their transition energies eps_n' - eps_n. | **no** |
| `_gap_ev` | Computes the HOMO-LUMO gap (eV) at a k-point as min(unoccupied eigenvalue) - max(occupied eigenvalue) using band-fraction thresholds. | **no** |
| `_load_fft_peaks` | Loads the excess-energy FFT spectrum CSV and returns the top-N (energy_ev, amplitude) peaks within an energy ceiling, for overlaying on the transition histogram. | **no** |
| `_select_gamma_kpoint` | Selects the k-point closest to Gamma by minimising \|k\|=\|\|(kx,ky,kz)\|\| across the unique k-points, returning its index and band slice. | **no** |
| `run` | Phase entry: loads eigenvalues+occupations, selects the Gamma k-point, builds the transition table and gap, overlays FFT peaks, and writes the histogram PNG and transitions CSV. | **no** |

### `postprocess/gs_projected_occupations.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_detect_homo_index` | Returns the highest state index with occupation >= threshold (default 0.5), giving the HOMO under INQ's spin-paired f_max=2 storage. | **no** |
| `_project_full` | Computes GS-projected occupations n_i^GS = sum_j f_j(0) O_ij as the matrix-vector product O_sq @ f_init over the evolved states. | **no** |
| `_project_proxy` | Computes GS-projected occupations using shell-averaged proxy overlaps: contrib_s(i) = (g_s/\|P_s\|) sum_{j in P_s} f_j(0) O_ij, summed over degenerate shells, with a column book-keeping check. | **no** |
| `run` | Phase entry: loads overlap snapshots, initial occupations and shell map, projects evolved occupations onto GS orbitals (full or proxy), and plots GS-orbital occupation evolution and depletion. | **no** |

### `postprocess/observables.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_build_variants` | Returns the four signal-preprocessing variants (raw_subtracted, mean_subtracted, linear-detrended, plateau_detrend) ahead of FFT; plateau_detrend subtracts the late-time mean per the QBall recipe. | **no** |
| `_hann_fft` | Applies a Hann window, zero-pads by pad_factor, and rffts a signal, returning freq/omega/energy(eV)/amplitude with amplitude normalised by N (not N_padded) for cross-run comparability. | **no** |

### `postprocess/state_energy_spectra.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_hann_fft` | Plateau-detrend then Hann-window, zero-pad and rfft a per-state energy signal, returning energy(eV)/amplitude/phase and the complex spectrum normalised by N. | **no** |
| `_peak_in_band` | Returns the dominant (omega_peak_eV, amplitude) of a spectrum masked to a [lo,hi] eV band, killing the DC tail and high-frequency noise. | **no** |
| `run` | Phase entry: FFTs each per-state KS energy epsilon_N(t), records per-state dominant peaks, and ranks opposite-phase (n,n') pairs by cross-spectrum amplitude (phase ~ pi) to discriminate electron-hole transitions from plasmons. | **no** |

### `postprocess/eigenvalues_gs.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_fermi_level_ev` | Estimates the Fermi level (eV) as the midpoint between the highest occupied and lowest unoccupied eigenvalue, using band-fraction occupations and falling back to the median when no gap is found. | **no** |
| `_plot_dos` | Computes a Gaussian-broadened density of states sum_i w_i N(eps;eps_i,sigma) on an energy grid, writes the DOS CSV, and plots it with the Fermi level marked. | **no** |

### `postprocess/energy_balance.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `run` | Computes the energy ledger per step: dE_WP (WP single-state), occupation-weighted dE_bath over non-WP states, observed total-energy drift, and the unaccounted residual indicating excitation into initially-empty states; writes CSV and plot. | **no** |

### `postprocess/bath_energy.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `run` | Reads per-orbital state_energies.csv, excludes the WP state index, sums occupation-weighted KS orbital energies per step to give bath energy E_bath(t) in Ha/eV plus its drift from t=0, and plots it. | **no** |

### `postprocess/kl_divergence.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `_kl` | Computes the Kullback-Leibler divergence sum_k p_k log(p_k/q_k) over bins where q_k>0, with an EPS floor on p to keep p log p finite when a bin transiently vanishes. | **no** |
| `run` | Phase entry: pivots the WP momentum histogram to (time x \|k\|), normalises each row, and computes KL(P_t\|\|P_0) against the t=0 reference distribution per step, writing CSV and a KL-vs-time plot. | **no** |

### `postprocess/wp_trajectory.py`

| Function / class | Quantity & formula (from graph summary) | Tested? |
|---|---|---|
| `run` | Phase entry: from observables.csv plots WP centre-of-density <r>(t) per axis, finite-difference velocity v_z(t) (should equal launch k0 at t=0), and the density-fluctuation L2 norm over time. | **no** |
