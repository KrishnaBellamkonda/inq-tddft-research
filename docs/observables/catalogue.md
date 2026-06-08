# Observables Catalogue — TDDFT Simulations

Comprehensive reference for every raw and post-processed observable across
all simulation types. Supersedes the coronene-only `docs/observables_reference.md`.

---

## 1. Raw observables (written by run.cpp)

### 1.1. Universal (all simulation types)

| Observable | File | Format | Cadence | Source |
|------------|------|--------|---------|--------|
| Total energy | `observables.csv` col `energy_total` | CSV | every step | `ObservablesWriter` |
| Kinetic energy | `observables.csv` col `energy_kinetic` | CSV | every step | `ObservablesWriter` |
| Hartree energy | `observables.csv` col `energy_hartree` | CSV | every step | `ObservablesWriter` |
| XC energy | `observables.csv` col `energy_xc` | CSV | every step | `ObservablesWriter` |
| Current (x,y,z) | `observables.csv` cols `current_{x,y,z}` | CSV | every step | `ObservablesWriter` |
| Dipole (x,y,z) | `observables.csv` cols `dipole_{x,y,z}` | CSV | every step | `ObservablesWriter` |
| Density L2 fluctuation | `observables.csv` col `density_l2` | CSV | every WRITE_EVERY | `DensityDelta` |
| GS eigenvalues | `eigenvalues/eigenvalues.csv` | CSV | once (t=0) | `eigenvalue_dump` |
| GS occupations | `eigenvalues/occupations.csv` | CSV | once (t=0) | `eigenvalue_dump` |
| Run summary | `run_summary.txt` | text | stub at start, final at end | direct ofstream |

### 1.2. Density fields

| Observable | File pattern | Format | Cadence | Sim types |
|------------|-------------|--------|---------|-----------|
| GS system density | `vti/density_gs_system/` | VTI (binary) | once (t=0) | all |
| RT system density | `vti/density_total/` or `vti/density_rt_total/` | VTI (binary) | every WRITE_EVERY | all |
| RT bath density | `vti/density_system/` or `vti/density_rt_system/` | VTI (binary) | every WRITE_EVERY | all |
| RT WP density | `vti/density_rt_wp/` | VTI (binary) | every WRITE_EVERY | coronene only |
| Density delta (raw) | `vti/density_delta/` | VTI (binary) | every WRITE_EVERY | jellium |
| Density delta (coarse) | `vti/density_delta_coarse/` | VTI (binary) | every WRITE_EVERY | jellium |

**Note:** Jellium runs write `density_total/` and `density_system/` (without `_rt_` prefix). Symlinks `density_rt_total → density_total` etc. are created by `install_schema_shims()` in analyse.py for pipeline compatibility.

### 1.3. Wave-packet specific (WP runs only)

| Observable | File | Format | Cadence | Source |
|------------|------|--------|---------|--------|
| WP config | `wp_config.txt` | text | once | direct ofstream |
| WP injection report | `wp_injection_report.txt` | text | once | direct ofstream |
| WP momentum stats | `wp_momentum_stats.csv` | CSV | every WRITE_EVERY | `WPMomentumStats` |
| WP real-space stats | `wp_real_space_stats.csv` | CSV | every WRITE_EVERY | `WPRealSpaceStats` (optional) |
| WP initial density | `wavepacket/density_wp_initial.vti` | VTI | once (t=0) | coronene only |
| WP initial wavefunction | `wavepacket/wavefunction_wp_initial.vti` | VTI (complex) | once (t=0) | coronene only |

**WP momentum stats columns:** `step, time_au, px_mean, py_mean, pz_mean, px2_mean, py2_mean, pz2_mean, sigma_px2, sigma_py2, sigma_pz2, e_kin_ha, norm_check`

**WP real-space stats columns:** `step, time_au, x_mean, y_mean, z_mean, x2_mean, y2_mean, z2_mean, sigma_x2, sigma_y2, sigma_z2, norm_check`

### 1.4. State-resolved observables

| Observable | File | Format | Cadence | Sim types |
|------------|------|--------|---------|-----------|
| State energies E_i(t) | `state_energies.csv` | CSV (long) | every 5×WRITE_EVERY | jellium (WP+classical) |
| State variance σ²_E | included in state_energies.csv | CSV | same | jellium |
| Occupations f_i(t) | `occupations_vs_time.csv` | CSV (long) | every 5×WRITE_EVERY | jellium |
| Momentum distribution n(\|k\|,t) | `momentum_distribution.csv` | CSV | every 10×WRITE_EVERY | jellium |
| Gamma transitions | `eigenvalues/gamma_transitions.csv` | CSV | once (from GS eigenvalues) | jellium WP |

### 1.5. Orbital overlap

| Observable | File pattern | Cadence | Sim types |
|------------|-------------|---------|-----------|
| WP-only overlap | `overlap/index.csv` + `overlap_NNNNNN.csv` | every 10 steps | jellium WP, coronene |
| Full O_ij matrix | `overlap_full/index.csv` + snapshots | t=0, t=final (WP); t=0, t=mid, t=final (classical) | jellium |
| Proxy overlap | `overlap_proxies/index.csv` + snapshots + `shells.csv` | every PROXY_STRIDE | jellium |

### 1.6. Classical projectile specific

| Observable | File | Format | Cadence | Source |
|------------|------|--------|---------|--------|
| Projectile track | `electron_track.csv` | CSV | every step | direct ofstream |

**Columns:** `step, time_au, x, y, z, vx, vy, vz, fx, fy, fz`
(Note: `fx,fy,fz` are placeholder zeros — actual force recovered from dv/dt in post-processing)

### 1.7. LEED screens (coronene only)

| Observable | File pattern | Cadence |
|------------|-------------|---------|
| Full-time accumulators | `screens/total/screen_NN.dat` | once (end of run) |
| Per-screen physics windows | `screens/time_windowed/screen_NN_t*_{forward,back}.dat` | once (end of run) |
| Paper-window accumulators | `screens/time_windowed/screen_NN_t*_paper.dat` | once (end of run) |
| Instantaneous snapshots | `screens/instantaneous/screen_NN_tNNNNNN.dat` | every SCREEN_SNAP_EVERY |
| Screen config | `screens/screen_config.csv` | once |
| Window ranges | `screens/window_ranges.csv` | once |

### 1.8. v2 additions (dt=0.01 runs with WP wavefunction saving)

| Observable | File pattern | Format | Cadence | Sim types |
|------------|-------------|--------|---------|-----------|
| WP orbital density | `vti/density_wp/density_t*.vti` | VTI (real) | every WF_WRITE_EVERY | jellium WP v2 |
| WP orbital wavefunction | `vti/wavefunction_wp/wavefunction_t*.vti` | VTI (complex) | every WF_WRITE_EVERY | jellium WP v2 |
| WP orbital wavefunction (coronene) | `vti/wavefunction_wp_rt/wavefunction_t*.vti` | VTI (complex) | every 5×WRITE_EVERY | coronene v2 |
| All-orbital wavefunction dump | `vti/orbitals_tf/orbital_*.vti` | VTI (complex) | once at t_IFW | jellium WP v2 (planned) |

**Note on density::total:** Verified 2026-05-24 that `density::total(electrons)` **includes** the WP orbital
(peak 0.174 vs background 0.0013 e/Bohr³). Both `density_total` and `density_system` VTI series are
byte-identical in jellium runs — no bath/WP separation in raw output. Bath-only density must be
reconstructed in post-processing: n_bath = n_total − |ψ_WP|² (requires saved WP wavefunction).

---

## 2. Post-processed observables (written by analyse.py)

### 2.1. Common base (ALL simulation types must produce these)

These outputs form the **minimum evaluation set** — every completed run
must produce them for consistent cross-run comparison.

| Output | Location | Phase | Description |
|--------|----------|-------|-------------|
| `REPORT.md` | `analysis/REPORT.md` | main | Auto-generated physics summary |
| `observables_summary.png` | `analysis/observables/` | observables | 3-row panel: energy, current, dipole vs time |
| `total_energy_vs_time.png` | `analysis/observables/` | observables | Energy conservation check |
| `all_energies_vs_time.png` | `analysis/observables/` | observables | Per-component energy evolution |
| `current_components_vs_time.png` | `analysis/observables/` | observables | J_x, J_y, J_z vs time |
| `dipole_components_vs_time.png` | `analysis/observables/` | observables | μ_x, μ_y, μ_z vs time |
| `density_fluctuation_l2.png` | `analysis/observables/` | observables | σ²_n(t) fluctuation metric |
| `fft_total_energy.png` | `analysis/observables/` | observables | Energy spectrum |
| `fft_current_{x,y,z}.png` | `analysis/observables/` | observables | Current spectra |
| `dipole_spectrum_{x,y,z}.png` | `analysis/observables/` | observables | Dipole absorption spectra |
| `eigenvalue_levels.png` | `analysis/observables/eigenvalues/` | eigenvalues_gs | KS orbital level diagram |
| `eigenvalues_dos.png` | `analysis/observables/eigenvalues/` | eigenvalues_gs | Density of states |
| `eigenvalue_bars.png` | `analysis/observables/eigenvalues/` | eigenvalues_gs | Coloured bar chart |
| `eigenvalue_table.txt` | `analysis/observables/eigenvalues/` | eigenvalues_gs | Numerical eigenvalue table |
| `density_gs_system_xy.png` | `analysis/ground_state/` | gs | GS density slice (xy plane) |
| `density_gs_z_profile.png` | `analysis/ground_state/` | gs | GS density z-profile |
| `gs_occupations.png` | `analysis/ground_state/` | gs | GS occupation bars |
| `layout_xz.png` | `analysis/layout/` | layout | Cell geometry + WP launch layout |
| Extended spectra | `analysis/observables/spectra/` | observables | 4 variants: raw, mean, detrended, plateau per quantity |

### 2.1b. Expanded mandatory set (added 2026-05-25)

These outputs are now mandatory for every run. They reproduce the report
figures from `docs/reports/14-05-2026-meeting-emilio/figures/` and
`docs/reports/2026-05-21-meeting-emilio/figures/`.

**Energy decomposition (Group B):**

| Output | Location | Description | Reference figure |
|--------|----------|-------------|-----------------|
| `energy_decomposition_classical_vs_wp.png` | `analysis/observables/` | 6-panel ΔE vs Δz (total, kinetic, Hartree, XC, bath sum, WP slot) with IFW shading. Classical vs WP overlay. | fig 05/06 |
| `energy_bookkeeping_bar.png` | `analysis/observables/` | Bar chart: ΔE_kinetic, ΔE_hartree, ΔE_xc at t_IFW | fig 07, fig C1 |

**Density evolution (Group D, beyond animations):**

| Output | Location | Description | Reference figure |
|--------|----------|-------------|-----------------|
| `density_z_profile_evolution.png` | `analysis/observables/` | z-profile heatmap vs time | fig 09 |
| `delta_density_xz_snapshots.png` | `analysis/observables/` | δn(x,y=0,z) at 4 selected times (lab frame) | fig D1 |
| `z_profile_diff_vs_free.png` | `analysis/observables/` | z-profile difference: jellium − free propagation |  |
| `density_diff_vs_free.png` | `analysis/observables/` | 2D δn difference: jellium − free at selected times |  |

**Eigenvalue evolution (Group E, beyond GS):**

| Output | Location | Description | Reference figure |
|--------|----------|-------------|-----------------|
| `ks_eigenenergy_evolution.png` | `analysis/observables/` | KS orbital energies vs time (static, all orbitals) | fig 10 |

**Orbital analysis (Group F):**

| Output | Location | Description | Reference figure |
|--------|----------|-------------|-----------------|
| `gs_basis_decomposition.png` | `analysis/observables/` | Per-orbital Δn_i^GS at t_end: depletion (occ) + excitation (virt) | fig 11 |
| `overlap_heatmap_log_wp.png` | `analysis/overlap/` | \|⟨ψ_i^GS\|ψ_j(t_end)⟩\|² heatmap, log scale, diagonal masked | fig 14 (WP) |
| `overlap_heatmap_log_classical.png` | `analysis/overlap/` | Same for classical companion | fig 14 (classical) |
| `overlap_heatmap_diff_wp_vs_classical.png` | `analysis/overlap/` | Difference: WP − classical overlap |  |

**Momentum and trajectory (Group G):**

| Output | Location | Description | Reference figure |
|--------|----------|-------------|-----------------|
| `momentum_band_free_vs_jellium.png` | `analysis/observables/` | ⟨p_z⟩ ± σ_p vs centroid z, 2-panel free vs jellium | fig M_A |
| `sigma_xyz_vs_time.png` | `analysis/observables/` | σ_x(t), σ_y(t), σ_z(t) from wp_real_space_stats.csv |  |

**Advanced (Group I — v2 runs with WP wavefunction saving):**

| Output | Location | Description |
|--------|----------|-------------|
| `wp_momentum_distribution_before_after.png` | `analysis/observables/` | \|ψ̃_WP(k)\|² at t=0 vs t_IFW. Gaussian → distorted by scattering. |
| `planewave_decomposition.png` | `analysis/observables/` | Evolved KS orbitals projected onto plane-wave basis (E = ℏ²\|G\|²/2m). Shows e-h transitions directly. |
| `spectral_weight_response.png` | `analysis/observables/` | W_resp(q_z, ω) with exact WP subtraction |
| `loss_function.png` | `analysis/observables/` | L(q_z, ω) = −(4π/q²) Im[χ] |
| `secondary_electron_yield.png` | `analysis/observables/` | δ(t) = ∫_vacuum n_SE dr in proxy vacuum region |

**2D momentum-space scattering map (Observable 2, verified 2026-05-25):**

Compute Δ|ψ̃(k_z, k_⊥)|² = P(t_f) − P(t_0) where P is the cylindrically-averaged
momentum distribution of the WP orbital. Reveals scattering channels:
- Blue blob at k₀: beam depletion (always present)
- Red inside elastic ring k²=k₀²: inelastic forward scattering (energy loss)
- Red on ring: elastic deflection (angle change, no energy loss)
- Red outside ring: non-physical (energy gain — numerical error)
The elastic ring k_z² + k_⊥² = k₀² geometrically separates energy-loss (inside) from
energy-conserving (on ring) events.

Implementation: load complex ψ_WP VTIs at t=0 and t_f, 3D FFT with dV normalisation,
cylindrical-average |ψ̃|² into (k_z, k_⊥) bins, subtract. Verified on
`run_wp_n162_L50_E100_sigma1_v2`: 33% beam depletion, inelastic-dominant (red inside ring),
elastic scattering suppressed (no weight on ring — consistent with V_ion(q≠0)=0 in jellium).

Requires: complex WP wavefunction VTIs (v2 runs).

**Plane-wave decomposition (new observable):**

For jellium, the GS orbitals are plane waves ψ_n^GS(r) = e^{iG_n·r}/√V with
energy E_n = |G_n|²/2 (a.u.). At the final IFW step, FFT each evolved orbital
ψ_j(r, t_f) to get its plane-wave amplitudes c_{n,j} = ⟨G_n|ψ_j(t_f)⟩.
Then the plane-wave-basis occupation is:

n_PW(E) = Σ_j f_j Σ_n |c_{n,j}|² δ(E − E_n)

Comparing n_PW(E) at t=0 (step function at E_F) vs t_f reveals:
- Occupied PW states that lost electrons (holes below E_F)
- Empty PW states that gained electrons (excitations above E_F)
- The energy scale of e-h transitions directly

Requires: all-orbital wavefunction dump at t_f (v2 runs only, ~3.2 GB).

### 2.2. Density animations (all types with VTI output)

| Output | Location | Phase |
|--------|----------|-------|
| `density/{total,system,delta}_{xy,xz,yz}.gif` | `analysis/density/` | density |
| `density/{total,system,delta}_z_profile.gif` | `analysis/density/` | density |
| `density/*_log.gif` | `analysis/density/` | density |
| `density/*_coarse_*.gif` | `analysis/density/` | density (if coarse VTI present) |

### 2.3. Jellium WP specific

| Output | Location | Phase |
|--------|----------|-------|
| `wp_position_vs_time.png` | `analysis/observables/` | wp_trajectory |
| `wp_velocity_vs_time.png` | `analysis/observables/` | wp_trajectory |
| `sigma_z_analytic_vs_time.png` | `analysis/observables/` | stopping |
| `gamma_gamma_transitions.png` | `analysis/observables/eigenvalues/` | eigenvalues_gs |
| `ks_energies_absolute.gif` | `analysis/observables/` | state_energies |
| `ks_energies_delta.gif` | `analysis/observables/` | state_energies |
| `occupations_absolute.gif` | `analysis/observables/` | occupations |
| `occupations_delta.gif` | `analysis/observables/` | occupations |
| `momentum_distribution.gif` | `analysis/observables/` | momentum |
| `momentum_heatmap.png` | `analysis/observables/` | momentum |
| `knudsen_ke_vs_t.png` + `.csv` | `analysis/observables/` | knudsen_ke |
| `kl_divergence_vs_t.png` + `.csv` | `analysis/observables/` | kl_divergence |
| `energy_balance.png` + `.csv` | `analysis/observables/` | energy_balance |
| `state_energy_spectra/*.png` | `analysis/observables/state_energy_spectra/` | state_energy_spectra |
| `n_q_m{1..6}.{csv,png}` | `analysis/observables/` | density_fourier (custom) |
| Overlap heatmap | `analysis/overlap/overlap_heatmap_t_end.png` | analyse_extras |
| GS-projected occupations | `analysis/observables/gs_projected_occupations/*.png` | analyse_extras |

### 2.4. Jellium classical specific

| Output | Location | Phase |
|--------|----------|-------|
| `classical_force_fixed.png` | `analysis/observables/` | analyse_extras |
| `delta_E_total_vs_time.png` | `analysis/observables/` | analyse_extras |
| `delta_E_total_vs_z.png` | `analysis/observables/` | analyse_extras |
| `running_slope_vs_z.png` | `analysis/observables/` | analyse_extras |
| `stopping_force_vs_z.png` | `analysis/observables/` | stopping |
| `dE_kinetic_vs_z.png` | `analysis/observables/` | stopping |
| `bath_energy_vs_time.png` + `.csv` | `analysis/observables/` | bath_energy |
| Same overlap heatmap + GS-projected occupations | `analysis/overlap/` + `observables/gs_projected_occupations/` | analyse_extras |
| Same state_energies/occupations/momentum GIFs | as WP above | same phases |

### 2.5. Coronene specific

| Output | Location | Phase |
|--------|----------|-------|
| `gs_orbital_gallery.png` | `analysis/ground_state/` | gs |
| `wp_position_vs_time.png` | `analysis/observables/` | wp_trajectory |
| Screens total grid | `analysis/screens/total/all_screens_grid.png` | screens |
| Per-screen patterns | `analysis/screens/total/screen_NN.png` + `_log.png` | screens |
| Screens IFFT | `analysis/screens/ifft/screen_NN_ifft_amp.png` + `_patterson.png` | screens |
| Screens time-windowed | `analysis/screens/time_windowed/*.png` | screens |
| Screens instantaneous | `analysis/screens/instantaneous/*.gif` | screens |
| Screen coordinate checks | `analysis/screens/coordinate_checks/*.png` | screens |
| WP-GS overlap animation | `analysis/overlap/wp_overlap_with_gs_orbitals.gif` | overlap |

---

## 3. Cross-reference: run inventory and observable coverage

### 3.1. Jellium WP runs (standard density, L=50, r_s≈5.69)

| Run | E (eV) | σ (Bohr) | Tier 1 | Tier 2 | WP mom stats | density_fourier | GS-proj occ | REPORT.md |
|-----|--------|----------|--------|--------|-------------|----------------|------------|-----------|
| run_wp_n162_L50_E20 | 20 | 5.0 | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E25 | 25 | 5.0 | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ |
| run_wp_n162_L50_E50_v2 | 50 | 5.0 | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E100 | 100 | 5.0 | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E300_v2 | 300 | 5.0 | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E600_v2 | 600 | 5.0 | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E100_sigma0p5 | 100 | 0.5 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E100_sigma1 | 100 | 1.0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E100_sigma3 | 100 | 3.0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E100_sigma8 | 100 | 8.0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E20_sigma1 | 20 | 1.0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ |
| run_wp_n162_L50_E25_sigma1 | 25 | 1.0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ |
| run_wp_n162_L50_E50_sigma1 | 50 | 1.0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_wp_n162_L50_E300_sigma1 | 300 | 1.0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

### 3.2. Jellium WP runs (high density, L=30, r_s≈3.41)

| Run | E (eV) | σ (Bohr) | Tier 1 | Tier 2 | REPORT.md |
|-----|--------|----------|--------|--------|-----------|
| run_wp_n162_L30_E100_highdens | 100 | 0.5 | ✓ | ✓ | ✓ |
| run_wp_n162_L30_E100_highdens_sigma1 | 100 | 1.0 | ✓ | ✓ | ✓ |
| run_wp_n162_L30_E50_highdens_sigma1 | 50 | 1.0 | ✓ | ✓ | ✓ |
| run_wp_n162_L30_E200_highdens_sigma1 | 200 | 1.0 | ✓ | ✓ | ✓ |
| run_wp_n162_L30_E300_highdens_sigma1 | 300 | 1.0 | ✓ | ✓ | ✓ |

### 3.3. Classical runs (standard density)

| Run | E (eV) | Tier 1 | electron_track | stopping | bath_energy | REPORT.md |
|-----|--------|--------|---------------|---------|------------|-----------|
| run_classical_n162_L50_E20 | 20 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_n162_L50_E25 | 25 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_n162_L50_E50_attempt2 | 50 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_n162_L50_E100 | 100 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_n162_L50_E600_v2 | 600 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_e1500_L50_cubic | 1500 | ✓ | ✓ | ✓ | ✓ | ✓ |

### 3.4. Classical runs (high density)

| Run | E (eV) | Tier 1 | electron_track | stopping | analyse_extras | REPORT.md |
|-----|--------|--------|---------------|---------|---------------|-----------|
| run_classical_n162_L30_E100_highdens | 100 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_n162_L30_E50_highdens | 50 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_n162_L30_E200_highdens | 200 | ✓ | ✓ | ✓ | ✓ | ✓ |
| run_classical_n162_L30_E300_highdens | 300 | ✓ | ✓ | ✓ | ✓ | ✓ |

### 3.5. Coronene runs

| Run | σ (Å) | E (eV) | Impact | LEED screens | overlap | REPORT.md |
|-----|-------|--------|--------|-------------|---------|-----------|
| run_base | 0.53 | 200 | centre | ✓ (20 screens) | ✓ | ✗ |
| run_cc_bond | 0.53 | 200 | C-C bond (4.028 Bohr x-offset) | ✓ (20 screens) | ✓ | ✓ |
| run_E30 | 0.53 | 30 | centre | ✓ | ✓ | ✗ |
| run_E800 | 0.53 | 800 | centre | ✓ | ✓ | ✗ |
| run_s0p33 | 0.33 | 200 | centre | ✓ | ✓ | ✗ |
| run_s3 | 3.0 | 200 | centre | ✓ | ✓ | ✗ |
| run_E800_s0p33 | 0.33 | 800 | centre | ✓ | ✓ | ✗ |
| run_E30_s3 | 3.0 | 30 | centre | ✓ | ✓ | ✗ |
| run_b6_35x35x80 | 0.53 | 200 | near (b=6 Bohr) | ✓ | ✓ | ✗ |
| run_b18_35x35x80 | 0.53 | 200 | far (b=18 Bohr) | ✓ | ✓ | ✗ |
| run_35x35x40 | 0.53 | 200 | centre (small box) | ✓ | ✓ | ✗ |

### 3.6. Custom/specialised observables

| Observable | Runs where used | Source |
|------------|----------------|--------|
| density_fourier (axial n_q modes) | All jellium with density VTI | `analyse_extras` + `density_fourier.py` |
| classical_force_fixed (F_z from dv/dt) | All classical | `analyse_extras` |
| delta_E_total_vs_z (windowed S) | All classical | `analyse_extras` |
| running_slope_vs_z (box-deficit diagnostic) | All classical | `analyse_extras` |
| GS-projected occupations (n_i^GS vs time) | All jellium with proxy overlap | `analyse_extras` + `gs_projected_occupations` phase |
| LEED IFFT (inverse FFT of screen patterns) | coronene run_cc_bond, run_base | `screens` phase |
| ParaView 3D volumes | coronene run_base only | `paraview` phase (manual) |

---

## 4. Observable gap analysis

### Missing REPORT.md (runs that completed but lack post-processing)

| Run | Reason |
|-----|--------|
| run_wp_n162_L50_E25 | No analyse.py run |
| run_wp_n162_L50_E20_sigma1 | analyse.py not run with full template |
| run_wp_n162_L50_E25_sigma1 | analyse.py not run with full template |
| All coronene runs except run_cc_bond | Pre-date per-run analyse.py pattern |

### Inconsistencies to resolve

1. **WP momentum stats**: Only present in σ=1 and σ-sweep runs (added later in campaign). Absent from early σ=5 energy sweep runs.
2. **density_fourier**: Requires VTI density series. Not produced if density phase skips.
3. **GS orbital VTIs**: Only coronene writes per-orbital GS densities. Jellium skips (expensive per-element loop).
4. **LEED IFFT**: Only run_cc_bond has the dedicated `screens/ifft/` subdirectory. Other coronene runs have IFFT data mixed into `screens/total/`.
