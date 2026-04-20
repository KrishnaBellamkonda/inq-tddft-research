# Plan: WP Real-Time Propagation — free, jellium, coronene

## Context

The `inqkit::WavePacket` injection API is validated (4-run coronene injection tutorial, all
norm_after=1.0). The next step is TDDFT propagation after injection, across three systems that
bracket increasing complexity: (1) free particle in a box (no ions — measures raw spreading),
(2) jellium (WP + uniform electron gas), (3) coronene (WP-molecule scattering from varying
distances).

A secondary goal is to understand why Tsubonoya et al. (PRB 90, 035416, 2014) saw negligible WP
spreading — the free-propagation runs answer this quantitatively.

---

## New directories

| Directory | System |
|---|---|
| `Tutorial/free-propagation-wp-rt/` | Free particle, no ions, finite cell |
| `ResearchProject/jellium/jellium-wp-rt/` | Jellium electron gas + WP |
| `ResearchProject/systems/coronene/coronene-wp-rt/` | Coronene + WP at varying D |

---

## 1. New Python API: `inq-stack/python/inqview/defaults.py`

**Create** this file. `config.py` already has `PlotDefaults`, `RenderDefaults`, `Theme` — no
changes needed there.

```python
def default_density_movie(
    series: FieldSeries,
    output_dir: Path | str,
    pv_executable: Path | str | None = None,
    fps: int = 12,
    render: VolumeRenderSpec | None = None,
    image_size: tuple[int, int] = (1600, 1200),
    atoms: AtomSpec | None = None,
    frame_stride: int = 1,
) -> dict[str, Path]:
    """FieldSeries → VTI series → ParaView PNG frames → GIF.
    Returns dict with keys 'gif', 'frames_dir', 'vti_dir'."""
```

Pipeline: `convert_real_series_to_vti` → scalar range from `VTISeriesResult.data_min/max`
→ `pv.render_vti_series` → `pv.build_gif`. Also implement `default_wavepacket_movie`
(same signature — for WP orbital density series).

**Update `inq-stack/python/inqview/__init__.py`** to export both functions.

---

## 2. Tutorial/free-propagation-wp-rt/

### Cell and electrons

- Cell: `34.771 × 34.771 × 89.856` bohr, **finite** (LEED paper dimensions)
- No ions
- `options::electrons{}.cutoff(40.0_Ha).extra_states(1)` → 0 occupied + 1 extra = WP at index 0
- GS call trivially empty (call for API consistency; max_steps=10 sufficient)
- WP center: `(Lx/2, Ly/2, Lz − 5σ)` — near top face, moving in −z
- Cell wall reflections occur at high momentum; σ(t) still measurable from z-profile width

### TDDFT parameters

- `dt = 0.02_atomictime`, `n_steps = 10 000` ≈ 4.83 fs
- Write density every 100 steps → 100 frames
- Write observables (energy, current_z) every step via `ObservablesWriter`

### Run table (one variable at a time from base)

| Run | σ (Å) | E_kin (eV) | k direction |
|-----|-------|-----------|-------------|
| `run_01_base` | 0.53 | 200 | −z |
| `run_02_low_momentum` | 0.53 | 50 | −z |
| `run_03_high_momentum` | 0.53 | 800 | −z |
| `run_04_tilted_45` | 0.53 | 200 | 45° in xz plane: kx = kz = −k₀/√2 |
| `run_05_transverse_x` | 0.53 | 200 | +x |
| `run_06_wide_sigma` | 2.0 | 200 | −z |
| `run_07_narrow_sigma` | 0.265 | 200 | −z |

### run.cpp structure (per run — only WP params change)

```cpp
// Cell + empty ions
auto cell = systems::cell::orthorhombic(LX*1.0_b, LY*1.0_b, LZ*1.0_b).finite();
systems::ions ions(cell);

// 0 occupied + 1 extra state
auto electrons = systems::electrons(ions,
    options::electrons{}.cutoff(40.0_Ha).extra_states(1));

// Trivial GS
ground_state::initial_guess(ions, electrons);
ground_state::calculate(ions, electrons, options::theory{}.lda(),
    options::ground_state{}.energy_tolerance(1e-4_Ha).max_steps(10));

// Inject WP (no orthogonalise needed — 0 occupied states)
auto wp = inqkit::WavePacket{}
    .center(WP_CX, WP_CY, WP_CZ).sigma(WP_SIGMA_BOHR).k0(WP_KX, WP_KY, WP_KZ);
auto report = wp.inject_into_last_extra_state(electrons, 1.0);

// t=0 density: density::total() is zero (empty system) + WP orbital
auto rho_occ = inqkit::fields::density::total(electrons);
auto rho_wp  = inqkit::fields::density::orbital(electrons, report.state_index);
for (std::size_t i = 0; i < rho_occ.values.size(); i++) rho_occ.values[i] += rho_wp.values[i];
density_writer.write(rho_occ, 0.0, 0);

// RT propagation
inqkit::RealTimeSession rt(ions, electrons, WRITE_EVERY);
rt.add([&](inqkit::StepContext const& ctx) {
    density_writer.write(inqkit::fields::density::total(*ctx.electrons), ctx.time_au, ctx.step);
});
real_time::propagate(ions, electrons, [&](auto const& data) { rt.step(data); },
    options::theory{}.lda(),
    options::real_time{}.num_steps(N_STEPS).dt(DT_AU*1.0_atomictime).observables_current());
```

### analysis.py structure (per run)

1. Load RT density series from `results/density_rt/`
2. Validate N_elec = 1.0 at each frame
3. Compute z-profile (integrate xy) → extract σ(t) from Gaussian width
4. Plot σ(t) vs analytic `σ₀√(1 + t²/σ₀⁴)` (a.u.)
5. Plot E(t) from observables CSV (expect constant = E_kin)
6. `default_density_movie(series, VISDIR / "density", pv_executable=PV_EXE)` → GIF

---

## 3. ResearchProject/jellium/jellium-wp-rt/

### Cell and electrons

Based on `ResearchProject/jellium/01_ground_state/run.cpp` (N=40, L=40 bohr, r_s=7.26 a₀):

```cpp
systems::ions ions(systems::cell::cubic(40.0 * 1.0_b).periodic());
auto electrons = systems::electrons(ions,
    options::electrons{}
        .spacing(0.50 * 1.0_b)
        .extra_electrons(40)
        .extra_states(3)
        .temperature(0.00862 * 1.0_eV),
    input::kpoints::gamma());
```

- Fermi smearing required (degenerate |n|²=3 shell at N=40)
- Full GS SCF needed (Broyden, 1e-4 Ha, max_steps=300)
- WP center: `(Lx/2, Ly/2, 5*sigma_bohr)` — near z=0 face, moving in +z
- Periodic cell: WP wraps after one traversal (~0.25 fs at 200 eV); 5 fs sees ~20 transits
- Add `.orthogonalise_against_occupied(electrons)` to WavePacket builder

### Run table

| Run | σ (Å) | E_kin (eV) | k direction |
|-----|-------|-----------|-------------|
| `run_01_base` | 0.53 | 200 | +z |
| `run_02_low_energy` | 0.53 | 50 | +z |
| `run_03_high_energy` | 0.53 | 400 | +z |
| `run_04_tilted_45` | 0.53 | 200 | 45° from +z |
| `run_05_wide_sigma` | 2.0 | 200 | +z |
| `run_06_narrow_sigma` | 0.265 | 200 | +z |

### t=0 density

```cpp
auto rho_occ = inqkit::fields::density::total(electrons);   // 40-electron jellium
auto rho_wp  = inqkit::fields::density::orbital(electrons, report.state_index);
for (std::size_t i = 0; i < rho_occ.values.size(); i++) rho_occ.values[i] += rho_wp.values[i];
density_writer.write(rho_occ, 0.0, 0);   // 41 e⁻ total at t=0
```

---

## 4. ResearchProject/systems/coronene/coronene-wp-rt/

### Cell and electrons

Identical to `Tutorial/coronene-wave-packet-with-inqkit/`:
- 34.771 × 34.771 × 89.856 bohr, finite, `coronene_centered.xyz`
- `options::electrons{}.cutoff(40.0_Ha).extra_states(3)`
- LDA, Broyden, mixing=0.1, ndim=8, 1e-4 Ha, max_steps=300

### Run table (vary D; σ=0.53 Å, E=200 eV unless noted)

| Run | D (Å) | σ (Å) | E_kin (eV) | Notes |
|-----|-------|-------|-----------|-------|
| `run_01_d635_base` | 6.35 | 0.53 | 200 | Paper reference |
| `run_02_d3` | 3.0 | 0.53 | 200 | Close approach |
| `run_03_d10` | 10.0 | 0.53 | 200 | Intermediate |
| `run_04_d15` | 15.0 | 0.53 | 200 | Far |
| `run_05_d20` | 20.0 | 0.53 | 200 | Very far — spreading dominates |
| `run_06_projectile` | 6.35 | 0.265 | 800 | Narrow high-E projectile limit |

### t=0 density

```cpp
// After injection (108 occupied + WP):
auto rho_occ = inqkit::fields::density::total(electrons);   // 108 e⁻
auto rho_wp  = inqkit::fields::density::orbital(electrons, report.state_index);
for (std::size_t i = 0; i < rho_occ.values.size(); i++) rho_occ.values[i] += rho_wp.values[i];
density_writer.write(rho_occ, 0.0, 0);   // 109 e⁻ total at t=0
```

---

## 5. docs/notes/wp_spreading_investigation.md

Create with:
- Physical spreading formula: `σ(t) = σ₀ √(1 + t²/σ₀⁴)` (a.u., m=ħ=1)
- For σ₀=1.0 bohr: doubles at t ≈ 1.73 a.u. ≈ 0.042 fs
- LEED paper ran T=0.25 fs; WP at 200 eV (v=3.834 bohr/a.u.) travels 12 bohr in 3.13 a.u.
  = 0.076 fs → arrives at coronene well before σ doubles
- To-do: use free-propagation run_01 z-profile to extract σ(t); compare to analytic;
  find (D, E_kin) regime where σ < 2σ₀ on arrival

---

## 6. Shared C++ utility for field addition (per run.cpp)

Inline function in each run.cpp:

```cpp
static void add_field_inplace(inqkit::fields::RealField3D & a,
                               inqkit::fields::RealField3D const& b) {
    for (std::size_t i = 0; i < a.values.size(); i++) a.values[i] += b.values[i];
}
```

---

## 7. Implementation order

1. `inq-stack/python/inqview/defaults.py` + update `__init__.py`
2. `Tutorial/free-propagation-wp-rt/` — 7 × (run.cpp + analysis.py)
3. `ResearchProject/jellium/jellium-wp-rt/` — 6 × (run.cpp + analysis.py)
4. `ResearchProject/systems/coronene/coronene-wp-rt/` — 6 × (run.cpp + analysis.py)
5. `docs/notes/wp_spreading_investigation.md`
6. Start runs: free-propagation run_01_base first for early validation, then sequentially

---

## 8. Validation targets

| Check | Expected |
|---|---|
| Free prop N_elec = 1.0 per frame | 1.000 ± 0.001 |
| Free prop energy constant | < 0.1% drift over 5 fs |
| Free prop σ(t) matches analytic formula | < 5% deviation |
| Jellium N_elec = 41.0 per frame | 41.0 ± 0.05 |
| Coronene N_elec = 109.0 per frame | 109.0 ± 0.1 |
| GIF produced with visible density motion | visual inspection |

---

## 9. File count

| Component | New files |
|---|---|
| `inqview/defaults.py` + `__init__.py` update | 2 |
| `docs/notes/wp_spreading_investigation.md` | 1 |
| `docs/plans/wp_rt_propagation.md` (this file) | 1 |
| free-propagation: 7 × (run.cpp + analysis.py) | 14 |
| jellium: 6 × (run.cpp + analysis.py) | 12 |
| coronene: 6 × (run.cpp + analysis.py) | 12 |
| **Total** | **42** |

Coronene runs need `coronene_centered.xyz` copied into each sub-directory.
