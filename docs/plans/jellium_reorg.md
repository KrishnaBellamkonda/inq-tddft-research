# Plan: Jellium folder reorganisation (mirror coronene canonical layout)

Status: **draft, not started**
Owner: skcb2
Created: 2026-04-30
Related handover: `docs/handovers/jellium_reorg.md` (to be created on first execution step)

---

## 1. Goal

Restructure `ResearchProject/systems/jellium/` so that:

1. The **top-level layout** mirrors `ResearchProject/systems/coronene/` exactly
   (`shared/`, `scripts/`, `save_gs/`, `checkpoints/`, `configurations/`,
   `hypotheses/`, `run_*/`, `dispatch.log`).
2. Every per-run **`results/` tree** mirrors the coronene `results/` schema
   (`raw/{ground_state,wavepacket,density,vti,observables,overlap,screens}` +
   `analysis/{density,ground_state,layout,observables,overlap,screens}` +
   `run_summary.txt`).
3. Each run computes **all observables that the coronene runs compute**, plus
   one new observable: the **full GS↔evolved KS-orbital overlap matrix**
   `O_ij(t) = |<ψ_i^GS | ψ_j(t)>|²` (already supported by
   `inqkit::observables::OrbitalOverlapMatrix::snapshot()` — only the call
   site needs to switch from `snapshot_wp_only()` to `snapshot()` and the
   post-processor needs new visualisation hooks).
4. Legacy material is moved out of the canonical tree into
   `Tutorial/jellium-legacy/` (parallel to the prior `Tutorial/coronene-legacy/`
   cleanup).
5. New runs are produced **from scratch** (not by reorganising old outputs)
   so the new layout is the source-of-truth from step 0.

---

## 2. Current → target layout

### 2.1 Current jellium/

```
ResearchProject/systems/jellium/
├── 01_ground_state/                  (legacy, GS exploration)
├── 02_ground_state_convergence/      (legacy, convergence study)
├── 03_free_gaussian_wp_propagation/  (legacy, free-WP propagation)
├── jellium-analytical/               (legacy, analytical reference)
└── jellium-wp-rt/                    (production-ish, to be replaced)
    ├── compare_observables.py
    ├── jellium_hypotheses.py
    ├── jellium_spectra.py
    ├── jellium-wp-rt.log
    ├── run_all_wp_rt.sh
    ├── hypotheses/{00_base, 01_wp_energy_spread, 02_wp_sigma_spread,
    │                03_open_vs_closed_shell, 04_tilted_propagation,
    │                05_electron_capture}
    └── run_0{1..7}_*/                (run.cpp + flat results/)
```

### 2.2 Target jellium/

```
ResearchProject/systems/jellium/
├── checkpoints/                      # GS .save dirs, one per cell+N config
│   ├── gs_L40_cubic_N38/             # closed-shell base (variants 1–6)
│   └── gs_L40_cubic_N40/             # open-shell (variant 7)
├── configurations/
│   └── jellium_wp_rt_base/           # frozen description of the canonical config
│       ├── README.md                 # provenance + parameters
│       └── paths.hpp                 # checkpoint path constant (matches coronene pattern)
├── hypotheses/                       # promoted from jellium-wp-rt/hypotheses
│   ├── 01_wp_energy_spread/
│   ├── 02_wp_sigma_spread/
│   ├── 03_open_vs_closed_shell/
│   ├── 04_tilted_propagation/
│   └── 05_electron_capture/
├── save_gs/
│   ├── gs_L40_cubic_N38/{run.cpp, build/, run.log, profile.dat}
│   └── gs_L40_cubic_N40/{run.cpp, build/, run.log, profile.dat}
├── shared/
│   ├── configs/      base.hpp + variant headers
│   └── cpp/          eigenvalues_writer.hpp, leed_screen_layout.hpp,
│                     results_paths.hpp, run_template.hpp
├── scripts/
│   ├── jellium_postprocess.py        # analogue of coronene_postprocess.py
│   ├── dispatch_runs.py
│   ├── repostprocess_all.sh
│   ├── retrofit_eigenvalues.py
│   └── run_queue.txt
├── run_base/
├── run_E50_s0p53/                    # was run_02_low_energy
├── run_E400_s0p53/                   # was run_03_high_energy
├── run_E200_s0p53_tilt45/            # was run_04_tilted_45
├── run_E200_s2p0/                    # was run_05_wide_sigma
├── run_E200_s0p265/                  # was run_06_narrow_sigma
├── run_E200_s0p53_N40/               # was run_07_open_shell
└── dispatch.log
```

Per-run directory contents are identical to coronene:

```
run_*/
├── run.cpp             # thin wrapper that picks a Cfg and calls run_propagation<Cfg>
├── run                 # built executable (gitignored)
├── run.log
├── profile.dat
├── build/              # cmake build dir (gitignored)
└── results/
    ├── run_summary.txt
    ├── raw/
    │   ├── ground_state/
    │   │   ├── density_system/                  (VTI series)
    │   │   ├── density_gs_orbitals/             (VTI series, one per occupied orbital)
    │   │   ├── eigenvalues.csv
    │   │   ├── occupations.csv
    │   │   └── summary.txt
    │   ├── wavepacket/
    │   │   ├── density_wp_initial/              (VTI)
    │   │   ├── wavefunction_wp_initial/         (complex VTI)
    │   │   ├── wavepacket_config.txt
    │   │   ├── injection_report.txt
    │   │   └── orthogonality_report.csv
    │   ├── density/
    │   │   ├── density_rt_total/                (VTI series)
    │   │   ├── density_rt_system/               (VTI series; jellium background)
    │   │   └── density_rt_wp/                   (VTI series)
    │   ├── vti/
    │   │   ├── density_gs_system/
    │   │   ├── density_gs_orbitals/
    │   │   ├── density_rt_total/
    │   │   ├── density_rt_system/
    │   │   ├── density_rt_wp/
    │   │   └── orbitals/
    │   ├── observables/
    │   │   ├── observables.csv                  (step, t, energies, current, dipole)
    │   │   └── eigenvalues/                     (copied from checkpoint)
    │   ├── overlap/
    │   │   ├── index.csv
    │   │   └── overlap_NNNNNN.csv               (full n_ref × (n_ref+1) matrix)
    │   └── screens/
    │       ├── total/                           (full-time accumulator, one .dat per screen)
    │       ├── instantaneous/                   (flat: <label>_tNNNNNN.dat)
    │       ├── time_windowed/                   (per-screen physics + paper window)
    │       ├── screen_config.csv
    │       └── window_ranges.csv
    └── analysis/
        ├── density/                             (slices, GIFs, ParaView 3D GIF)
        ├── ground_state/                        (GS density plots, eigenvalue table)
        ├── layout/                              (sketch of cell + screens + WP)
        ├── observables/                         (energy, current, dipole vs t)
        ├── overlap/                             (heatmap_t<step>.png, diagonal_vs_t.png,
        │                                         wp_column_vs_t.png, orbital_population.csv)
        └── screens/                             (2D + smoothed spectra + ParaView 3D GIF)
```

---

## 3. Mapping table (current → new)

| Current path                                                             | New path                                                                  | Action |
|---|---|---|
| `jellium/01_ground_state/`                                               | `Tutorial/jellium-legacy/01_ground_state/`                                | git mv |
| `jellium/02_ground_state_convergence/`                                   | `Tutorial/jellium-legacy/02_ground_state_convergence/`                    | git mv |
| `jellium/03_free_gaussian_wp_propagation/`                               | `Tutorial/jellium-legacy/03_free_gaussian_wp_propagation/`                | git mv |
| `jellium/jellium-analytical/`                                            | `Tutorial/jellium-legacy/jellium-analytical/`                             | git mv |
| `jellium/jellium-wp-rt/run_0{1..7}_*/`                                   | `Tutorial/jellium-legacy/jellium-wp-rt/run_0{1..7}_*/`                    | git mv |
| `jellium/jellium-wp-rt/hypotheses/`                                      | `jellium/hypotheses/` (promoted up one level, regenerated from new runs) | git mv (then delete stale figures, regenerate) |
| `jellium/jellium-wp-rt/{compare_observables.py, jellium_hypotheses.py, jellium_spectra.py, run_all_wp_rt.sh}` | `Tutorial/jellium-legacy/jellium-wp-rt/`                | git mv (will be superseded by new `scripts/jellium_postprocess.py`) |
| `jellium/jellium-wp-rt/jellium-wp-rt.log`                                | `Tutorial/jellium-legacy/jellium-wp-rt/`                                  | git mv |
| `jellium/01_ground_state/jellium_utils.hpp`                              | starting point only — content copied & adapted into `jellium/shared/cpp/` headers (do NOT keep dual sources) | manual port |

After the moves, `ResearchProject/systems/jellium/` is empty of legacy content
and contains only the new canonical tree.

---

## 4. Run renaming (deviations from base)

Base parameters (frozen from `jellium-wp-rt/run_01_base/run.cpp`):

- Cell: `cubic(40.0 bohr).periodic()`
- Electrons: `extra_electrons(38) extra_states(3) temperature(0.00862 eV) spacing(0.50 bohr)`
- Theory: LDA (ALDA in TDDFT), Γ-only k-points
- WP: σ=0.53 Å (1.0015 bohr), E=200 eV, k₀=3.834 bohr⁻¹, +z, centred at
  `(L/2, L/2, 5σ)` = `(20, 20, 5.008) bohr`
- Propagation: dt=0.020 a.u., N_steps=417, write_every=2, screen_snap_every=3
- Screens: 20 z-positions in [0.5, 39.5] bohr (existing layout retained)

Variants — each is a `struct : Base { … overrides … };`:

| New name                   | Was                | Overrides vs Base                             | N_STEPS |
|---|---|---|---|
| `run_base`                 | `run_01_base`      | (none)                                        | 417 |
| `run_E50_s0p53`            | `run_02_low_energy`| `WP_EKIN_EV = 50`                             | 834 |
| `run_E400_s0p53`           | `run_03_high_energy`| `WP_EKIN_EV = 400`                           | 295 |
| `run_E200_s0p53_tilt45`    | `run_04_tilted_45` | `WP_KX = +k₀/√2, WP_KZ = +k₀/√2`              | 350 |
| `run_E200_s2p0`            | `run_05_wide_sigma`| `WP_SIGMA_ANG = 2.0` ⇒ σ=3.779 bohr, WP_CZ=5σ | 480 |
| `run_E200_s0p265`          | `run_06_narrow_sigma`| `WP_SIGMA_ANG = 0.265` ⇒ σ=0.501 bohr, WP_CZ=5σ | 480 |
| `run_E200_s0p53_N40`       | `run_07_open_shell`| `N_ELECTRONS = 40` (uses `gs_L40_cubic_N40`)  | 417 |

`N_STEPS` numbers are inherited from the legacy runs (already validated for
loop-back margin in the cubic-periodic cell). Future work: derive N_STEPS at
compile time from a `compute_n_steps()` helper analogous to coronene's.

---

## 5. New shared C++ library (jellium/shared/cpp/)

All four headers live in namespace `jellium::*`. They are direct ports of
the coronene equivalents with the listed changes:

### 5.1 `results_paths.hpp` (namespace `jellium::results`)

Identical to `coronene/shared/cpp/results_paths.hpp` modulo the namespace.
Same path strings — the `results/` schema is identical.

### 5.2 `eigenvalues_writer.hpp` (namespace `jellium::eigenvalues`)

Identical to coronene equivalent (the routine just reads
`<checkpoint>/kpin0000000000/eigenvalues.csv` etc. and copies into
`results/raw/observables/eigenvalues/`). No system-specific logic.

### 5.3 `leed_screen_layout.hpp` (namespace `jellium::layout`)

Adapted from coronene's. Differences:

- Cubic-periodic cell with z ∈ [0, L_BOHR] (not centred). Screen z-positions
  carry over from `jellium-wp-rt/run_01_base/run.cpp` (the existing 20-screen
  jittered layout: 0.5 / 2.53 / … / 39.5 bohr).
- `compute_screen_window(z_screen, wp_cz, wp_sigma, k0_z, L, t_total, n_sigmas)`
  re-derived for cubic-periodic geometry with periodic boundary handling
  (back-scattering window models a wrap rather than a hard rebound).
- `screen_label(k)` and `zero_pad6(step)` identical to coronene.
- `N_SCREENS = 20`.

### 5.4 `run_template.hpp` (namespace `jellium::run_template`)

Direct port of coronene's `run_template.hpp` with these substitutions:

| Coronene                                                  | Jellium |
|---|---|
| `systems::cell::orthorhombic(LX,LY,LZ).finite()`          | `systems::cell::cubic(L_BOHR).periodic()` |
| `systems::ions::parse(geometry_xyz_path, cell)`           | `systems::ions(cell)` (empty) |
| `options::electrons{}.cutoff(...).extra_states(...)`      | `options::electrons{}.spacing(0.50 bohr).extra_electrons(N_ELECTRONS).extra_states(EXTRA_STATES).temperature(0.00862 eV)` |
| signature `run_propagation(name, geom_path, gs_path)`     | signature `run_propagation(name, gs_path)` (no geometry path) |
| `overlap_obs.snapshot_wp_only(...)`                       | `overlap_obs.snapshot(...)` ← **the new observable** |
| `density_rt_jellium`                                      | unchanged at I/O level: written to `results/raw/density/density_rt_system/` and `results/raw/vti/density_rt_system/` (same as coronene) |
| stub run_summary `geometry_file`                          | `n_ions = 0; geometry_file = "(none, jellium)"; n_electrons = N_ELECTRONS` |
| paper-window `T1_AU/T2_AU` from Tsubonoya 2014            | jellium has no canonical paper window — set `T1_AU = 0, T2_AU = DT_AU * N_STEPS` (full propagation) and document this in the variant header. Future work: re-derive from a jellium reference. |

The full body otherwise tracks coronene line-for-line (stub run_summary,
GS density write, GS orbital writes, WP injection + reports, three RT
density writers, observables CSV, three screen accumulators per screen, final
run_summary with the same section numbering).

---

## 6. Shared configs (jellium/shared/configs/)

### 6.1 `base.hpp`

```cpp
#pragma once
#include <cmath>
namespace jellium::config {

inline constexpr double ANG_TO_BOHR = 1.8897259886;
inline constexpr double HA_TO_EV    = 27.21138625;
inline constexpr double FS_TO_AU    = 41.341374575751;

inline constexpr double const_sqrt(double x) { /* identical to coronene */ }
inline constexpr double k0_from_ev(double e) { return const_sqrt(2.0 * e / HA_TO_EV); }

struct Base {
    // Cell (cubic, periodic)
    static constexpr double L_BOHR = 40.0;
    static constexpr double LX_BOHR = L_BOHR;
    static constexpr double LY_BOHR = L_BOHR;
    static constexpr double LZ_BOHR = L_BOHR;

    // Electronic structure
    static constexpr int    N_ELECTRONS    = 38;
    static constexpr int    EXTRA_STATES   = 3;
    static constexpr double SPACING_BOHR   = 0.50;
    static constexpr double TEMPERATURE_EV = 0.00862;
    static constexpr double SCF_TOL_HA     = 1.0e-4;
    static constexpr int    SCF_MAX_STEPS  = 300;
    static constexpr int    SCF_MIX_NDIM   = 8;
    static constexpr double SCF_MIX_ALPHA  = 0.1;

    // No CUTOFF_HA in jellium (uses spacing instead). Stub for log compatibility.
    static constexpr double CUTOFF_HA = 0.0;

    // Wave packet
    static constexpr double WP_SIGMA_ANG  = 0.53;
    static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;
    static constexpr double WP_EKIN_EV    = 200.0;
    static constexpr double WP_K0         = k0_from_ev(WP_EKIN_EV);

    static constexpr double WP_CX_BOHR = 0.5 * L_BOHR;
    static constexpr double WP_CY_BOHR = 0.5 * L_BOHR;
    static constexpr double WP_CZ_BOHR = 5.0 * WP_SIGMA_BOHR;

    static constexpr double WP_KX = 0.0;
    static constexpr double WP_KY = 0.0;
    static constexpr double WP_KZ = +WP_K0;             // +z launch (matches legacy)
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // Real-time
    static constexpr double DT_AU             = 0.020;
    static constexpr int    N_STEPS           = 417;
    static constexpr int    WRITE_EVERY       = 2;
    static constexpr int    SCREEN_SNAP_EVERY = 3;

    // LEED window — placeholder (no jellium reference paper currently used).
    static constexpr int    N_SCREENS = 20;
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double WP_ENVELOPE_SIGMAS = 2.0;
};

} // namespace jellium::config
```

### 6.2 Variant headers

Each is a 5–10 line struct. Examples:

```cpp
// shared/configs/E50_s0p53.hpp
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E50_s0p53 : Base {
    static constexpr double WP_EKIN_EV = 50.0;
    static constexpr double WP_K0      = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KZ      = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;
    static constexpr int    N_STEPS    = 834;
};
} // namespace jellium::config
```

```cpp
// shared/configs/E200_s0p53_tilt45.hpp
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E200_s0p53_tilt45 : Base {
    static constexpr double WP_KX = +Base::WP_K0 / 1.41421356237;
    static constexpr double WP_KZ = +Base::WP_K0 / 1.41421356237;
    static constexpr int    N_STEPS = 350;
};
} // namespace jellium::config
```

```cpp
// shared/configs/E200_s2p0.hpp
struct E200_s2p0 : Base {
    static constexpr double WP_SIGMA_ANG  = 2.0;
    static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;
    static constexpr double WP_CZ_BOHR    = 5.0 * WP_SIGMA_BOHR;
    static constexpr int    N_STEPS       = 480;
};
```

```cpp
// shared/configs/E200_s0p53_N40.hpp
struct E200_s0p53_N40 : Base {
    static constexpr int N_ELECTRONS = 40;
    // Same N_STEPS as base (417) — propagation parameters unchanged.
};
```

(Full set in §3.)

---

## 7. save_gs and checkpoints

Closed-shell base (variants 1–6) share one checkpoint:

`save_gs/gs_L40_cubic_N38/run.cpp`

```cpp
#include "../../shared/configs/base.hpp"
using Cfg = jellium::config::Base;

int main() {
    using namespace inq; using namespace inq::magnitude;
    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(Cfg::SCF_TOL_HA * 1.0_Ha)
            .max_steps(Cfg::SCF_MAX_STEPS)
            .broyden_mixing()
            .mixing_ndim(Cfg::SCF_MIX_NDIM)
            .mixing(Cfg::SCF_MIX_ALPHA));
    electrons.save("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N38");
    return 0;
}
```

Open-shell variant gets its own:

`save_gs/gs_L40_cubic_N40/run.cpp` — identical, with
`Cfg = jellium::config::E200_s0p53_N40` and
`electrons.save(".../checkpoints/gs_L40_cubic_N40")`.

Variant `run.cpp` files are then trivial:

```cpp
// run_base/run.cpp
#include "../shared/configs/base.hpp"
#include "../shared/cpp/run_template.hpp"
int main() {
    return jellium::run_template::run_propagation<jellium::config::Base>(
        "run_base",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N38");
}
```

```cpp
// run_E200_s0p53_N40/run.cpp
#include "../shared/configs/E200_s0p53_N40.hpp"
#include "../shared/cpp/run_template.hpp"
int main() {
    return jellium::run_template::run_propagation<jellium::config::E200_s0p53_N40>(
        "run_E200_s0p53_N40",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N40");
}
```

---

## 8. New observable: full GS↔evolved overlap matrix

### 8.1 Definition

`O_ij(t) = |<ψ_i^GS | ψ_j(t)>|²`

- i ∈ [0, n_ref-1] indexes the n_ref = wp_idx GS reference orbitals
  (occupied subspace at t=0, frozen at WP-injection time).
- j ∈ [0, n_ref] indexes evolved orbitals; column `j = n_ref` is the WP.

### 8.2 Already-existing infrastructure

`inqkit::observables::OrbitalOverlapMatrix::snapshot(electrons, t, step)`
computes the full matrix, writes `output_dir/overlap_NNNNNN.csv`, and
appends to `index.csv`. Performance: ~50–500 ms per step at 80³ for n_ref≈19,
documented in the header.

The current jellium `run_01_base.cpp` already calls `snapshot()`. Coronene's
template uses `snapshot_wp_only()` (cheaper). The new jellium template
uses `snapshot()` to satisfy the new requirement.

### 8.3 Validation

- **t=0 identity check**: at step 0, `O_ij ≈ δ_ij` for i,j < n_ref and
  `O_{i, n_ref} ≈ 0` (WP was orthogonalised against occupied states). Print
  the diagonal max-deviation and the WP-column max in the run log; flag
  values larger than 1e-3.

### 8.4 Post-processing

`scripts/jellium_postprocess.py` writes:

| File                                                     | Description |
|---|---|
| `analysis/overlap/heatmap_t<step>.png`                   | n_ref × (n_ref+1) matrix at sampled time slices |
| `analysis/overlap/heatmap_anim.gif`                      | animated heatmap |
| `analysis/overlap/diagonal_vs_t.png`                     | `O_ii(t)` per GS orbital — measures loss of orbital identity under WP perturbation |
| `analysis/overlap/wp_column_vs_t.png`                    | `O_{i, n_ref}(t)` — which GS state the evolved WP "looks like" over time |
| `analysis/overlap/orbital_population_vs_t.csv`           | `Σ_i O_{ij}(t)` for each evolved column j |
| `analysis/overlap/orbital_population_vs_t.png`           | the same as a stacked-area plot |
| `analysis/overlap/identity_residual.csv`                 | `‖O(t) − I_extended‖_F` (Frobenius distance from identity-extended-by-zero) |

Implementation note: matrix is small (≤ 20×21) per step, ≤ 500 steps → trivial
to load all of `overlap_*.csv` into a single (T, n_ref, n_ref+1) numpy tensor.

---

## 9. Post-processing script

`scripts/jellium_postprocess.py` is a near-line-for-line port of
`coronene/scripts/coronene_postprocess.py` with these changes:

- Reads `results/raw/density/density_rt_system/` instead of
  `density_rt_system/` produced by the coronene template (path is identical
  by §5.4 — no code change).
- Adds the §8.4 overlap visualisations.
- LEED screen post-processing identical (the schemas in `screens/` match).
- Layout sketch (`analysis/layout/`) draws a cubic-periodic box with no atoms.
- All outputs as `.png` per `.claude/rules/file-placement.md`.
- Standard 6 outputs for LEED runs preserved (per memory:
  `analysis/screens/{combined_grid_2d.png, …}` + 2D GIF + ParaView 3D GIF).

---

## 10. Validation plan

Per `.claude/rules/testing.md`, propose the following menu before any
expensive run. The user picks which to authorise.

**Cheap (run by default)**

1. Compile-only: `inq-run --build-only` for `save_gs/gs_L40_cubic_N38/run.cpp`
   and `run_base/run.cpp` — confirm the new headers resolve and the templates
   instantiate without compile errors.
2. Header port unit checks: small `Tutorial/` smoke driver that calls
   `jellium::results::*` and confirms the 30+ paths all resolve and create
   directories on demand.
3. `compute_screen_window` unit test: hand-computed cases for
   {forward, backward (periodic wrap)} screens at three z values; check
   `t_start, t_end` to within 1 dt.

**GS sanity (one short SCF)**

4. Run `save_gs/gs_L40_cubic_N38/`: SCF energy must agree with the legacy
   `run_01_base/run.cpp` GS energy (already in
   `Tutorial/jellium-legacy/jellium-wp-rt/run_01_base.log` after the move) to
   within 1e-4 Ha. Also: density at z = L/2 ≈ uniform jellium target
   `n = 38 / L^3 = 5.94e-4 e/bohr³` ± 5 %.

**WP injection sanity (no propagation)**

5. Run `run_base` with `N_STEPS = 0` (override locally): GS load + WP inject
   only. Confirm:
   - `report.norm_after ∈ [0.97, 1.03]`
   - `report.max_overlap ≤ 1e-3`
   - `overlap_000000.csv`: diagonal of n_ref×n_ref block ≥ 1−1e-6, off-diagonal
     ≤ 1e-6, WP column ≤ 1e-3.
   - All 30+ raw paths populated.

**Short propagation**

6. Run `run_base` with `N_STEPS = 20`: 1-frame VTI series in each density
   stream, ~10 overlap CSVs, no NaNs in `observables.csv`, energy drift over
   20 steps < 1e-3 Ha. Check `analysis/` outputs after running the
   post-processor.

**Full base run**

7. Full `run_base` (N_STEPS=417): observables.csv energy_total drift compared
   to the legacy `run_01_base` should agree (same physics) within 1e-4 Ha
   over the run; final overlap-matrix Frobenius residual `‖O(T) − I‖_F`
   recorded as scientific output, not a pass/fail criterion.

**Full suite**

8. Six remaining variants (`E50_s0p53` … `E200_s0p53_N40`). User-authorised
   before kick-off.

---

## 11. Order of execution

| Step | Action                                                              | Validation gate |
|---|---|---|
| 1   | Write this plan (this file)                                          | — |
| 2   | Create `docs/handovers/jellium_reorg.md`                             | — |
| 3   | `git mv` legacy material into `Tutorial/jellium-legacy/` (§3 table). Sanity check: `find ResearchProject/systems/jellium -maxdepth 2 -type d` returns only the new skeleton dirs after step 4. | — |
| 4   | Create empty skeleton: `mkdir -p jellium/{shared/{configs,cpp},scripts,save_gs,checkpoints,configurations/jellium_wp_rt_base,hypotheses}` | — |
| 5   | Port `results_paths.hpp`, `eigenvalues_writer.hpp`,                  | unit checks #1, #2 |
|     | `leed_screen_layout.hpp` (adapted), `run_template.hpp`               | unit check #3 |
|     | into `jellium/shared/cpp/`                                           |   |
| 6   | Write `shared/configs/base.hpp` and `save_gs/gs_L40_cubic_N38/run.cpp` | compile #1, GS sanity #4 |
| 7   | Write `run_base/run.cpp` and the 6 variant headers + run dirs         | compile #1 (all), WP-inject sanity #5 |
| 8   | Run `run_base` short propagation                                      | #6 |
| 9   | Port `coronene_postprocess.py` → `scripts/jellium_postprocess.py`. Run on `run_base/results/`. | manual review of all outputs |
| 10  | Add overlap-matrix visualisations (§8.4)                              | manual review |
| 11  | Run `run_base` full propagation                                       | #7 |
| 12  | Write `save_gs/gs_L40_cubic_N40/run.cpp`, run it                      | analogue of #4 |
| 13  | Run remaining 6 variants (user-authorised)                            | #8 |
| 14  | Regenerate `hypotheses/` from new outputs                             | manual review |
| 15  | Delete the now-superseded legacy files in `Tutorial/jellium-legacy/`? | user decision |

---

## 12. Open decisions (need user input before step 5)

A. **Open-shell N value** (variant `run_E200_s0p53_N40`): keep N=40 from the
   legacy `run_07_open_shell`? Default: yes.
B. **GS checkpoint sharing**: variants 02–06 reuse `gs_L40_cubic_N38`; variant
   07 uses `gs_L40_cubic_N40`. Default: yes.
C. **Tilted variant cell**: keep cubic L=40 box, only change `WP_KX/WP_KZ`?
   Default: yes (matches legacy `run_04_tilted_45`).
D. **Periodic vs finite cell**: keep `.periodic()` (physically required for
   homogeneous electron gas — finite jellium is unphysical and does not
   match the legacy runs). Default: yes.
E. **`hypotheses/`**: regenerate-only from new runs (recommended) vs copy
   stale figures. Default: regenerate-only.

If any default above is wrong, edit this plan before step 5.

---

## 13. Sources / attribution

- Coronene canonical layout: `ResearchProject/systems/coronene/` (in-repo).
- Coronene template: `coronene/shared/cpp/run_template.hpp`,
  `results_paths.hpp`, `leed_screen_layout.hpp`, `eigenvalues_writer.hpp`,
  `shared/configs/tsubonoya_2014_base.hpp`.
- Existing jellium numerics (frozen as the new Base): legacy
  `jellium-wp-rt/run_01_base/run.cpp` (parameters trace to
  earlier exploration; no published reference paper specific to this
  N=38 / L=40 bohr / 200 eV combination is on file — flag as
  *internal-reference, not literature*).
- Tsubonoya, Hu, Watanabe, *Phys. Rev. B* **90**, 035416 (2014) — used by
  coronene as the paper-window reference; **not used** by jellium runs
  (T1_AU/T2_AU set to full-propagation in §6.1).
- `inqkit::observables::OrbitalOverlapMatrix`:
  `inq-stack/include/inqkit/observables/orbital_overlap.hpp` (project-local
  implementation, header-documented).

---

## 14. Notes / known caveats

- **Loop-back margin** in cubic-periodic cell: legacy `N_STEPS` values were
  hand-tuned to keep the WP from wrapping. Keep them initially; in a follow-up
  derive `N_STEPS` at compile time from `compute_n_steps()` adapted for
  periodic geometry.
- **No paper window for jellium**: the `T1_AU/T2_AU` mechanism is preserved
  only for layout/postprocess compatibility with the coronene template. The
  paper accumulator output for jellium is just the full-time accumulator and
  may be removed in a follow-up.
- **VTI metadata for periodic cells**: confirm `inqkit::io::RealField3DWriter`
  emits correct ParaView origin/spacing for `[0, L]` cubic-periodic — coronene
  uses centred `[-L/2, L/2]` so the writer assumption may need a check
  during the §10 step-9 visualisation review.
- **Eigenvalue retrofit**: use `scripts/retrofit_eigenvalues.py` (port from
  coronene) on the GS checkpoint if the SCF run predates the eigenvalue
  writer wiring.
