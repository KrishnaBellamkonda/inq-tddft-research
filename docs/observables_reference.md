# Observables Reference — Coronene WP-RT-LEED Framework

Authoritative reference for every numerical artefact written and read by the
coronene replication framework under
`ResearchProject/systems/coronene/`. Lists every output path,
the C++/Python source that populates it, and the schema. Companion to:

* `docs/results_folder_structure_spec.md` — the binding `results/` layout spec.
* `docs/notes/postprocess-algorithms.md` — the math and the file:line for
  every algorithm.
* `docs/notes/coronene-geometry-correction.md` — the z = L/2 → z = 0 fix.
* `docs/handovers/coronene-cumulative.md` — narrative overview.

> The legacy version of this file documented the *jellium* WP-RT runs and is
> obsolete; that history lives in git only.

---

## 1. The C++ run writes (every coronene propagation run)

Sources: `ResearchProject/systems/coronene/shared/cpp/run_template.hpp`,
shared `RealField3DWriter` / `ComplexField3DWriter` / `ObservablesWriter`
under `inq-stack/include/inqkit/io/`, and the screen accumulators under
`inq-stack/include/inqkit/screens/`.

| Observable | C++ source | Output path | Frequency |
|---|---|---|---|
| `run_summary.txt` (final) | run_template.hpp final block | `results/run_summary.txt` | once at end (stub at start) |
| GS density (system, before WP) | `density::total(electrons)` | `results/raw/vti/density_gs_system/density_gs_system.vti` | once |
| GS orbital densities | `density::orbital(electrons, i)` for `i ∈ [0, n_occupied)` | `results/raw/vti/density_gs_orbitals/orbital_NNNN.vti` | once each |
| WP config + injection report | direct ofstream | `results/raw/wavepacket/{wavepacket_config.txt, injection_report.txt}` | once |
| Initial WP density | `density::orbital(electrons, wp_idx)` | `results/raw/vti/.../density_wp_initial.vti` (under `wavepacket/`) | once |
| Initial WP wavefunction (complex) | `fields::orbital::wavefunction(electrons, wp_idx)` (Phase-3 fix: now fftshifted) | `results/raw/wavepacket/wavefunction_wp_initial/wavefunction_wp_initial.vti` | once |
| `density_rt_system` (real, RT) | `density::total(electrons)` | `results/raw/vti/density_rt_system/density_t<step>.vti` | every `WRITE_EVERY` steps |
| `density_rt_wp` (real, RT) | `density::orbital(electrons, wp_idx)` | `results/raw/vti/density_rt_wp/density_t<step>.vti` | every `WRITE_EVERY` steps |
| `density_rt_total` (real, RT) | `density_rt_system + density_rt_wp` (pointwise) | `results/raw/vti/density_rt_total/density_t<step>.vti` | every `WRITE_EVERY` steps |
| Total energy, KE, Hartree, XC | `ObservablesWriter` | `results/raw/observables/observables.csv` | every step |
| Current `J_x, J_y, J_z` | `ObservablesWriter` | same CSV, additional columns | every step |
| Dipole `μ_x, μ_y, μ_z` | `ObservablesWriter` | same CSV, additional columns | every step |
| WP-overlap row | `OrbitalOverlapMatrix::snapshot_wp_only` | `results/raw/overlap/overlap_<step>.csv` + `index.csv` | every step |
| Total LEED screens (full-time accum.) | `LeedPatternAccumulator::save` | `results/raw/screens/total/screen_NN.dat` | end of run (20 screens) |
| Per-screen physics-window accums. | `LeedPatternAccumulator::save` | `results/raw/screens/time_windowed/screen_NN_t<lo>_to_t<hi>_{forward,back}.dat` | end of run |
| Paper-window accums. | `LeedPatternAccumulator::save` | `results/raw/screens/time_windowed/screen_NN_t<lo>_to_t<hi>_paper.dat` | end of run |
| Instantaneous LEED snapshots | `PlaneScreen::extract` + save | `results/raw/screens/instantaneous/screen_NN_t<step>.dat` | every `SCREEN_SNAP_EVERY` steps |
| Screen config + window ranges | direct ofstream | `results/raw/screens/{screen_config.csv, window_ranges.csv}` | once at start |

Density writers all use `emit_raw=false, emit_vti=true,
vti_format=binary` (no `.raw` sidecars). Filename pattern
`density_t<step:06d>.vti` from `RealField3DWriter::write(field, time, step)`.

---

## 2. Three density categories (and how they relate)

INQ's `density::total(electrons)` returns the sum of occupied KS orbital
densities **excluding** the WP extra state, because the WP is injected into
a slot with explicit `occupation = 1.0` outside the originally-occupied
manifold. Therefore:

```
density_rt_system(r, t) = ρ_KS,occupied(r, t)
density_rt_wp(r, t)     = |ψ_wp(r, t)|²  (single-orbital density)
density_rt_total(r, t)  = density_rt_system + density_rt_wp     ← computed at write time
```

The pointwise add lives in `run_template.hpp::add_real_fields`. All three
fields share the same grid, origin (= −L/2), and spacing so the add is
trivial.

---

## 3. WP-only orbital overlap (replaced the full O_ij matrix)

The full `n_ref × (n_ref + 1)` matrix snapshot at every step was prohibitively
expensive. Phase-1 swapped it for the **WP-only row**:

```
O_i,wp(t) = |dV · Σ_r conj(ψ_i^GS(r)) · ψ_wp(r, t)|²,   i ∈ [0, n_occupied)
```

Implementation: `inqkit::observables::OrbitalOverlapMatrix::snapshot_wp_only`
in `inq-stack/include/inqkit/observables/orbital_overlap.hpp:90`.

**CSV format per snapshot** (one row per file):

```
# step=<k> time_au=<t> n_ref=<n_occupied> mode=wp_only
o0,o1,o2,...,o<n_ref-1>
```

`results/raw/overlap/index.csv` columns: `step,time_au,file`.

---

## 4. Per-screen physics windows (Phase-3 corrected)

Source: `coronene::layout::compute_screen_window` in
`ResearchProject/systems/coronene/shared/cpp/leed_screen_layout.hpp`.

Three accumulators per screen:

| Accumulator | Window | Saved as |
|---|---|---|
| `acc_full[k]` | All steps (no gating) | `total/screen_NN.dat` |
| `acc_screen_window[k]` | Per-screen physics window (below) | `time_windowed/screen_NN_t<lo>_to_t<hi>_{forward,back}.dat` |
| `acc_paper[k]` | Global paper window `[T1, T2] = [0.077, 0.25] fs` | `time_windowed/screen_NN_t<lo>_to_t<hi>_paper.dat` |

The per-screen window is designed to **exclude the unscattered Gaussian WP**
crossing the screen plane:

* **Forward screens (z_screen < 0, transmission side)**:
  * `t_start = max(0, (b + σ − z_screen) / |k|)` — when the WP trailing
    edge has cleared the screen on its way down. After this time, density
    at z_screen is the transmitted+diffracted contribution.
  * `t_end = N_steps · dt` — until end of run.

* **Backscattering screens (z_screen ≥ 0, reflection side)**:
  * `t_start = max(0, (b + σ − z_screen) / |k|)` — same expression.
    For z_screen > b + σ this clamps to 0 (screen never under the WP).
  * `t_end = (b + L_z/2 − σ) / |k|` — when the rebound forward leading
    edge reaches the +L_z/2 box face, before periodic-boundary wrap-around.

`screen_config.csv` records `(screen_index, z_bohr, label, kind, t_start_au, t_end_au)`
per screen so the postprocess can reproduce the same naming without
re-deriving the windows.

---

## 5. LEED screen `.dat` file format

Written by `LeedPatternAccumulator::save` in
`inq-stack/include/inqkit/screens/leed_pattern_accumulator.hpp`.

```
# label=<screen_NN> z=<z_bohr> total_time=<T_au> n_accum=<N>
# nx=<NX> ny=<NY> dx=<DX> dy=<DY> origin_x=0.000000 origin_y=0.000000
v_00 v_01 ... v_0,NX-1
v_10 ...
...
v_NY-1,0 ... v_NY-1,NX-1
```

Two important conventions:

* **Array order is FFT-natural**: `(0, 0)` is the physical origin
  `(x = 0, y = 0)` (cell centre). Array indices > N/2 map to negative
  physical coordinates.
* **`origin_x = origin_y = 0` in the header is correct** for that array
  ordering, but downstream consumers must apply `np.fft.fftshift` to
  recover a centred image; the Python loader does this automatically.

`inqview.load_leed_pattern` (`inq-stack/python/inqview/screens.py`) applies
`np.fft.fftshift(data)` on read and **overrides** the metadata origin to
`(-Lx/2, -Ly/2)` so `LeedPattern.extent_bohr` spans
`[-Lx/2, +Lx/2, -Ly/2, +Ly/2]`. Every consumer downstream (the postprocess,
the IFFT helper, plot_leed_pattern) sees a centred-frame field.

---

## 6. Inverse-FFT reconstruction (Phase-3 addition)

LEED is intensity-only — phase is lost. Two reconstruction methods are
provided in `inq-stack/python/inqview/postprocess/_ifft.py`:

| Method | Math | Interpretation |
|---|---|---|
| `patterson` (default) | `IFFT(|F|²)` (Wiener–Khinchin) → `ρ ⋆ ρ` | Pair-correlation peaks at every interatomic separation in the projected density. Physically defensible. |
| `amp_only` | `|IFFT(√|F|²)|²` (zero phase) | Coarse heuristic; not a true reconstruction. |

Public API:

```python
import inqview as iv

pat = iv.load_leed_pattern("results/raw/screens/total/screen_03.dat")
recon = pat.inverse_fft(method="patterson", hann=True)   # shape (ny, nx)
```

Postprocess phase `screens` writes one PNG per method per screen:

```
results/analysis/screens/total/
  screen_NN.png                       (existing — linear)
  screen_NN_log.png                   (existing — log)
  screen_NN_ifft_patterson.png        (new — Patterson)
  screen_NN_ifft_amp.png              (new — amp-only)
```

---

## 7. Coordinate conventions

INQ stores all 3D real-space fields and orbitals in **FFT-natural** order:

* Array index `(0, 0, 0)` = physical position `(0, 0, 0)` (cell centre).
* Indices `> N/2` along any axis wrap to negative physical coordinates.
* `basis.symmetric_range_begin(i) = -N_i / 2`; the centred frame spans
  `[-L_i/2, +L_i/2]`.

Two FFT-ordering bugs were fixed in Phase 3:

1. **`PlaneScreen::iz_nearest`** — converts physical `z_target` to grid
   index. Previously clamped negative indices to 0; now wraps via
   `((iz % Nz) + Nz) % Nz`. Files:
   `inq-stack/include/inqkit/screens/plane_screen.hpp`. Symptom before
   fix: every transmission screen at z < 0 was sampling z = 0 (the
   molecule plane) instead of the requested plane.

2. **`fields::orbital::wavefunction`** — exports complex KS orbitals.
   Previously skipped the `fft_shift_index` mapping that the sibling
   `density::total` / `density::orbital` exporters apply, so the
   `ComplexField3D` was FFT-natural while metadata claimed left-to-right
   physical layout. Now applies the same shift as density.hpp. Files:
   `inq-stack/include/inqkit/fields/orbital.hpp`. Symptom before fix:
   complex orbital VTIs and Python slices showed the wavefunction in
   scrambled positions relative to the origin.

Both bugs were silent — array indices were valid, no exception was
raised — and only surfaced because every "transmission" LEED screen
showed the static coronene cloud (the molecule plane the screens were
mistakenly sampling).

---

## 8. Postprocess derived outputs (Python)

Source: `inq-stack/python/inqview/postprocess/`. Each phase consumes
files under `results/raw/` and writes to `results/analysis/`.

| Phase | Inputs | Outputs |
|---|---|---|
| `summary` | `run_summary.txt` | appends post-processing block |
| `gs` | `raw/vti/density_gs_system/`, `raw/vti/density_gs_orbitals/` | `analysis/ground_state/{density_gs_system_xy.png, gs_orbital_gallery.png}` |
| `layout` | `run_summary.txt` + `raw/screens/screen_config.csv` | `analysis/layout/layout_xz.png` |
| `observables` | `raw/observables/observables.csv` | `analysis/observables/{total_energy_vs_time, all_energies_vs_time, current_components_vs_time, dipole_components_vs_time, observables_summary, fft_*}.png`; numerical FFTs to `raw/observables/fft_*.csv` |
| `density` | `raw/vti/density_rt_{total,system,wp}/` | `analysis/density/{total,system,wp}_{xy,xz,yz}{_log}.{gif,mp4}` |
| `screens` | `raw/screens/{total,instantaneous,time_windowed}/`, `screen_config.csv` | `analysis/screens/{total,instantaneous,time_windowed,coordinate_checks}/...` (linear + log per artefact; IFFT-Patterson + IFFT-amp on each total screen) |
| `overlap` | `raw/overlap/index.csv` + per-step CSVs | `analysis/overlap/wp_overlap_with_gs_orbitals{,_log}.{gif,mp4}` |
| `orbitals` | `raw/vti/orbitals/` (RT orbitals; currently empty) | (skipped) |
| `paraview_3d` | `raw/vti/density_rt_{system,wp}/` | `analysis/density/paraview_3d/volume_overlay_{view_headon,view_3q}.{gif,mp4}` |
| `paraview` | (legacy single-series renderer) | (opt-in via `--with-paraview`) |

Cross-run hypothesis comparisons (`coronene_postprocess.py hypothesis`)
emit:

```
hypotheses/<NN>_*/
  README.md                                 (predefined)
  leed_total_grid.png
  peak_intensity_vs_label.png
  energy_drift_overlay.png
  current_z_overlay.png
  physics/
    current_{x,y,z}_overlay.png
    dipole_{x,y,z}_overlay.png
    energy_total_spectrum_overlay.png
    {current,dipole}_{x,y,z}_spectrum_overlay.png
    wp_overlap_residual_at_t_final.png
```

---

## 9. Screen z-positions (run_base, L_z = 60 Bohr, 20 screens)

Generated by `coronene::layout::screen_z_positions(60.0)`:

```
-29.0  -27.4  -24.4  -20.9  -17.9  -14.5  -11.4   -7.9   -4.9   -1.5
  0.07   3.5    6.5   10.0  13.0  16.4  19.4  22.9  25.9  29.0
```

(Approximate; exact values include the `±0.07 / ±0.13` parity jitter to
avoid grid-plane aliasing.) Screens at z < 0 are forward (transmission);
z ≥ 0 are backscattering. Other runs use the same generator with their
own `L_z`.

---

## 10. Numerical algorithms — quick math index

For full derivations see `docs/notes/postprocess-algorithms.md`.

| # | Quantity | Form |
|---|---|---|
| 1 | DFT of an observable | `X[k] = Σ x[n] · w[n] · exp(-2π i k n / N)` with Hann `w` and pre-detrend |
| 2 | Orbital overlap | `O_ij = |dV · Σ_r conj(ψ_i^GS(r)) · ψ_j(r, t)|²` |
| 3 | MGS at WP injection | `ψ_wp ← ψ_wp − Σ_i ⟨ψ_i^GS|ψ_wp⟩ ψ_i^GS`, then renormalise |
| 4 | WP construction | `ψ_wp(r) = (πσ²)^{-3/4} exp(-|r-b|²/(2σ²)) exp(i k₀·r)` (centred frame) |
| 5 | End-of-box time | `t_end = (b + σ + L_z/2)/|k|`; `N_steps = round(t_end/dt)` |
| 6 | `density_rt_total` | pointwise add of `system + wp` at write time |
| 7 | LEED fftshift | `data ← np.fft.fftshift(data)` on read; origin overridden to `-L/2` |
| 8 | Patterson IFFT | `IFFT(|F|²) = ρ ⋆ ρ` (Wiener–Khinchin) |

---

## 11. Extended preprocessed spectra

For each of `dipole_z`, `current_z`, and `energy_total`, the postprocess
builds three preprocessed signals before FFT-ing:

| Variant | Preprocessing | Purpose |
|---|---|---|
| **raw_subtracted** | `s − s(0)` | Removes the initial value; cheapest detrend; preserves any drift correlated with `t = 0`. |
| **mean_subtracted** | `s − ⟨s⟩` | Removes DC; minimum bias on stationary signals. |
| **detrended** | `scipy.signal.detrend(s, type='linear')` | Removes a linear least-squares fit; isolates oscillatory content from packet drift. Most physical for WP scattering runs. |

All three variants use the same downstream pipeline: Hann window →
zero-pad by `pad_factor = 4` → `np.fft.rfft` → `np.fft.rfftfreq(N_pad, dt_au)`.
Zero-padding **only smooths the visible curve**; the intrinsic spectral
resolution (peak width ≈ `1/(N·dt_au)`) is unchanged.

**Frequency / energy axes**:

```
freq_au   = np.fft.rfftfreq(N_pad, d=dt_au)        # cycles / a.u.-time
omega_au  = 2π · freq_au                           # angular frequency, Ha
energy_ev = 27.21138625 · omega_au                 # photon-energy axis
```

Plots cap the displayed range at 200 eV (a comfortable upper bound for
KS-orbital-difference physics in coronene).

**Output paths** (per run, compartmentalised by quantity):

```
results/analysis/observables/spectra/
  current/   spectrum_current_z_{raw_subtracted,mean_subtracted,detrended,compare}.png
  dipole/    spectrum_dipole_z_{raw_subtracted,mean_subtracted,detrended,compare}.png
  energy/    spectrum_energy_total_{raw_subtracted,mean_subtracted,detrended,compare}.png

results/raw/observables/spectra/
  current/   spectrum_current_z_<variant>.csv     (cols: freq_au, omega_au, energy_ev, amplitude)
  dipole/    spectrum_dipole_z_<variant>.csv
  energy/    spectrum_energy_total_<variant>.csv
```

The `_compare.png` per quantity overlays the three variants on a single
axes for direct comparison — peaks that survive in *detrended* are most
physically meaningful.

**Code**: `inq-stack/python/inqview/postprocess/observables.py`
(`_extended_spectra`, `_hann_fft`, `_build_variants`, `_plot_compare`).

---

## 12. Ground-state orbital eigenvalues + occupations

The KS orbital eigenvalues (`ε_i`) and their fractional occupations (`f_i`)
of the ground state are essential context for every WP-RT run: they fix
the HOMO/LUMO gap, mark the WP injection slot's energy, and seed the
overlap analysis. Every coronene run therefore writes them once, right
after `electrons.load(<checkpoint>)` and before WP injection.

**Writer**: `run_template.hpp::run_propagation` calls a small helper
`coronene::run_template::dump_eigenvalues(electrons, dir)` that flushes
two CSVs.

**Output paths**:

```
results/raw/observables/eigenvalues/
  eigenvalues.csv      # cols: state_index, eigenvalue_ha, eigenvalue_ev
  occupations.csv      # cols: state_index, occupation
```

The state index runs `0..n_states-1` (with `n_states = n_occupied +
extra_states`); the WP slot is `state_index = n_occupied + extra_states - 1`
(populated to occupation 1.0 only after the WP injection block runs).

**Retrofit for existing runs**: `scripts/extract_eigenvalues_from_log.py`
parses each `save_gs/<sig>/run.log` (the GS save's INQ log already prints
all eigenvalues per state at SCF convergence) and copies the resulting
CSVs into every run that loaded the matching checkpoint.

**Postprocess viz** in `analysis/observables/eigenvalues/`:

| File | Content |
|---|---|
| `eigenvalues_levels.png` | horizontal level diagram, eV scale, occupied vs unoccupied colour-coded; HOMO/LUMO marked; WP slot dashed-line annotated. |
| `eigenvalues_dos.png` | density-of-states histogram (Gaussian-broadened, `σ_DOS = 0.1 eV` default), eV axis. |
| `eigenvalue_table.txt` | plain-text dump of state, ε_ha, ε_ev, occ for quick reference. |

These are produced by the `observables` postprocess phase whenever the
input CSVs exist; phase is skipped silently otherwise.
