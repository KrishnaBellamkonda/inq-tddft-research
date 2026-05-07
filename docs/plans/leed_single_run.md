# Plan: Coronene LEED — focused single-run (Lz=31.7 Å)

## Context

Previous 6 coronene runs produced no satisfactory LEED pattern. Root causes: box too large,
WP spread too much before hitting target, screens poorly positioned, no shared colour scale
for comparison. This plan sets up a single carefully designed run with a smaller box, tighter
WP (σ=1 Å), lower energy (100 eV), 22 screens for full spatial coverage, and GIF/VTI
output at every 3 steps. Workflow is gated: GS → user approves → RT → user approves → run.

---

## Computed parameters

### Unit constants

```
ANG_TO_BOHR = 1.8897259886
HA_TO_EV    = 27.21138625
AU_TO_FS    = 0.02418884327
```

### Box (all in bohr in code)

| Quantity | Angstrom | Bohr |
|---|---|---|
| Lx = Ly | 18.48 | 34.9222 |
| Lz | 31.7 | 59.9043 |
| Z_flake = Lz/2 | 15.850 | 29.9522 |

### Wave packet

| Quantity | Value |
|---|---|
| σ | 1.0 Å = 1.8897 bohr |
| E_kin | 100 eV = 3.67493 Ha |
| k0 = sqrt(2·E_Ha) | 2.71107 bohr⁻¹ |
| v = k0 | 2.71107 bohr/a.u. |
| D (WP–flake gap) | 6.35 Å = 11.9997 bohr |
| WP_CZ = Z_flake + D | — | 41.9519 bohr |
| WP_CX = Lx/2 | — | 17.4611 bohr |
| WP direction | −z | |

### 5σ sanity check

- 5σ = 9.449 bohr
- WP to top wall: 59.9043 − 41.9519 = **17.95 bohr = 9.50σ ✓**
- Mirror transmission screen to bottom wall: 17.95 bohr = 9.50σ ✓

### Timing — **NEEDS USER CONFIRMATION**

| Quantity | Atomic units | Femtoseconds |
|---|---|---|
| t_wp (WP reaches flake, D/v) | **4.426 a.u.** | **0.1071 fs** |
| t_final (WP reaches z=0, WP_CZ/v) | **15.474 a.u.** | **0.3743 fs** |
| step_wp (start accumulating screens) | **step 222** | — |
| dt | 0.02 a.u. | — |
| N_steps = ceil(t_final/dt) | **774 steps** | — |

### Screens (22 total)

All screens accumulate from step 222 to step 774 (WP arrival to end).

**Backscatter side (z > Z_flake):**

| Label | z (bohr) | Note |
|---|---|---|
| bs_main | 41.9519 | WP start, main backscatter |
| bs_10 | 57.1814 | equally spaced in (Z_flake, Lz) |
| bs_09 | 54.4585 | |
| bs_08 | 51.7355 | |
| bs_07 | 49.0126 | |
| bs_06 | 46.2897 | |
| bs_05 | 43.5668 | |
| bs_04 | 40.8439 | |
| bs_03 | 38.1209 | |
| bs_02 | 35.3980 | |
| bs_01 | 32.6751 | closest to flake |

Spacing = 29.9521 / 11 = 2.7229 bohr. No screen at Z_flake=29.9522 ✓

**Transmission side (z < Z_flake):**

| Label | z (bohr) | Note |
|---|---|---|
| tr_01 | 27.2292 | closest to flake |
| tr_02 | 24.5063 | |
| tr_03 | 21.7834 | |
| tr_04 | 19.0605 | |
| tr_05 | 16.3375 | |
| tr_06 | 13.6146 | |
| tr_07 | 10.8917 | |
| tr_08 | 8.1688 | |
| tr_09 | 5.4458 | |
| tr_10 | 2.7229 | |
| tr_main | 17.9524 | mirror of bs_main (Z_flake − D) |

Spacing = 29.9522 / 11 = 2.7229 bohr. No screen at Z_flake=29.9522 ✓

---

## New coronene XYZ file

Existing `coronene_centered.xyz` centroid: **(9.200, 9.200, 23.775) Å**
(made for original large box, Lz/2 = 44.928 bohr = 23.775 Å)

New required centroid: **(9.2408, 9.2408, 15.8496) Å** (= Lx/2, Ly/2, Lz/2 in new box)

Shift applied to every atom: **Δx = +0.0408, Δy = +0.0408, Δz = −7.9254 Å**

File: `Tutorial/coronene-leed/run_01/coronene_leed.xyz`

---

## Directory layout

```
Tutorial/coronene-leed/run_01/
  coronene_leed.xyz        ← new centred XYZ (generate first)
  gs.cpp                   ← Phase 1: GS + WP injection + density write
  gs_to_vti.py             ← Phase 1: convert GS density outputs to VTI
  rt.cpp                   ← Phase 2: RT propagation (write after GS approval)
  analysis.py              ← Phase 2: screens, GIFs, LEED plots
  results/
    density_gs/            ← t=0 total density (coronene + WP)
    density_wp_gs/         ← t=0 WP-only density
    vti/                   ← VTI files from gs_to_vti.py
    density_rt/            ← RT total density every 3 steps
    density_wp_rt/         ← RT WP-only density every 3 steps
    observables.csv        ← all observables, every step
    screens/               ← 22 .dat LEED pattern files
    leed_plots/            ← individual + shared-scale PNG plots
    gifs/                  ← animated GIFs from density series
```

---

## Phase 1: Ground state + WP injection (`gs.cpp`)

1. Build cell: `orthorhombic(34.9222_b, 34.9222_b, 59.9043_b).finite()`
2. Parse `coronene_leed.xyz`
3. Compute GS with `extra_states(3)`, cutoff 40 Ha, LDA, Broyden mixing
4. Build WP: center=(WP_CX, WP_CX, WP_CZ), σ=1.8897 bohr, k0=(0,0,−2.71107), `.orthogonalise_against_occupied(electrons)`
5. Inject: `report = wp.inject_into_last_extra_state(electrons, 1.0)`
6. Print `report.max_overlap`. If `report.max_overlap > 1e-3`, re-orthogonalise and re-inject (loop max 3 times), print final overlap.
7. Write `results/density_gs/` with `density::total(electrons)` (already includes WP since occ=1)
8. Write `results/density_wp_gs/` with `density::orbital(electrons, report.state_index)`

**Key note**: `inject_into_last_extra_state` sets `occ[WP]=1.0`, so `density::total()` already includes the WP. No manual `add_field_inplace` needed. This is confirmed from the previous run analysis.

---

## Phase 1: VTI conversion (`gs_to_vti.py`)

Uses `inqview.vti.convert_real_meta_to_vti`:

```python
from inqview import vti, data
sim = data.SimulationData(RUN_DIR)
series_tot = sim.field_series("results/density_gs")
series_wp  = sim.field_series("results/density_wp_gs")
# Convert single t=0 frame each
vti.convert_real_meta_to_vti(list(series_tot.files)[0],
    output_path=RUN_DIR/"results/vti/density_gs_t0.vti", array_name="density")
vti.convert_real_meta_to_vti(list(series_wp.files)[0],
    output_path=RUN_DIR/"results/vti/density_wp_gs_t0.vti", array_name="density_wp")
```

User opens both VTI files in ParaView to confirm: coronene at z=29.95 bohr, WP at z=41.95 bohr.

---

## Phase 2: RT propagation (`rt.cpp`) — write after GS approval

Key parameters:
- `N_STEPS = 774`, `DT_AU = 0.02`, `WRITE_EVERY = 3`
- `STEP_WP_START = 222` (accumulate screens only from this step)

Structure:
1. Same cell/ions/GS computation as `gs.cpp`
2. Same WP injection
3. Two `RealField3DWriter`s: `density_rt/` and `density_wp_rt/`, write every 3 steps
4. `ObservablesWriter` every step (all observables: energy, kinetic, hartree, xc, current xyz, dipole xyz)
5. 22 `LeedPatternAccumulator` objects (one per screen)
6. RT session for density (WRITE_EVERY=3):
   ```cpp
   rt.add([&](StepContext const& ctx) {
       density_writer.write(density::total(*ctx.electrons), ctx.time_au, ctx.step);
       wp_density_writer.write(density::orbital(*ctx.electrons, report.state_index), ...);
   });
   ```
7. RT session for observables + screens (every step):
   ```cpp
   rt_obs.add([&](StepContext const& ctx) {
       obs_writer.append(ctx);
       if (ctx.step >= STEP_WP_START) {
           sc_bs_main.accumulate(..., report.state_index);
           // ... all 22 screens
       }
   });
   ```
8. Save all 22 screen `.dat` files

---

## Phase 2: Analysis (`analysis.py`) — write after GS approval

1. Load `density_rt/` series → `convert_real_series_to_vti` → `results/vti/density_rt/`
2. Load `density_wp_rt/` series → VTI → `results/vti/density_wp_rt/`
3. Build GIFs from VTI series (using `pvbatch` or `imageio`)
4. Load all 22 screen `.dat` files
5. Plot grid of 22 screens in logical order:
   - Row 1 (backscatter, right to left, far to near): bs_main, bs_10…bs_01
   - Row 2 (transmission, near to far): tr_01…tr_10, tr_main
6. Two versions: individual colour scale (per panel) and shared colour scale (all panels same range)
7. Save to `results/leed_plots/leed_individual.png` and `results/leed_plots/leed_shared.png`
8. Plot observables summary (energy, currents, dipole)
9. Print N_elec per frame (expect ≈ 109 ± 0.5 for coronene + 1 WP)

---

## Verification checklist

| Check | When | Pass criterion |
|---|---|---|
| WP visible in `density_wp_gs_t0.vti` | After GS | Gaussian blob at z≈41.95 bohr |
| Coronene visible in `density_gs_t0.vti` | After GS | Ring structure at z≈29.95 bohr |
| `report.max_overlap < 1e-3` | After GS | Printed value ≤ 1e-3 |
| N_elec ≈ 109 every RT frame | After RT | Mean deviation < 0.5 |
| Energy drift < 1e-3 Ha over run | After RT | Check observables.csv |
| All 22 screen files written | After RT | `ls results/screens/*.dat` = 22 files |
| LEED plots generated | After analysis | Both PNG files exist |

---

## Files to create (in order)

1. `Tutorial/coronene-leed/run_01/coronene_leed.xyz` (shift from existing)
2. `Tutorial/coronene-leed/run_01/gs.cpp`
3. `Tutorial/coronene-leed/run_01/gs_to_vti.py`
4. `Tutorial/coronene-leed/run_01/rt.cpp` (after GS approval)
5. `Tutorial/coronene-leed/run_01/analysis.py` (after GS approval)
6. `docs/plans/leed_single_run.md` (copy of this plan, project-tracked)

---

## Timing estimate

- GS: ~5–10 min (coronene with extra_states(3), Lz=59.9 bohr box)
- RT: 774 steps × ~2.5 s/step (estimated from previous runs, Lz≈60 bohr) ≈ **32 min**
- Analysis + GIFs: ~5 min
