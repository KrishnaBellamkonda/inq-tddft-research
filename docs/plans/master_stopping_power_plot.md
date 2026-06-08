# Plan: Master Stopping Power Plot

**Created**: 2026-05-25
**Status**: Ready for implementation
**Priority**: High — this is the central figure of the report

---

## Objective

Build the master stopping power plot S(v) or S(E) from raw simulation data with full traceability. Two versions: standard density (L=50, r_s ≈ 5.69) and high density (L=30, r_s ≈ 3.41). Each version uses two WP stopping power definitions plus the classical Ehrenfest reference.

---

## Design decisions (locked in)

| Decision | Choice |
|---|---|
| Primary WP series | σ=1 Bohr and σ=5 Bohr (both have energy sweeps) |
| Supplementary WP points | σ=0.5, 3, 8 at E=100 eV only |
| Classical line | All available classical runs at matched density |
| x-axis | Projectile energy E (eV), log scale |
| y-axis | Stopping power S (eV/Bohr), log scale |
| Time axis units | femtoseconds throughout |
| t_start | t = 0 (injection time) |
| t_end | Interference-free: when centroid + 4σ_r(t) reaches the periodic boundary at L/2 |
| WP S definition 1 | S₁ = −ΔE_kin / Δz, where E_kin = ⟨p⟩²/2m + σ_p²/2m |
| WP S definition 2 | S₂ = −ΔE_KS / Δz, where E_KS is the KS orbital energy of the WP state |
| Classical S | Linear regression of E_system(z) over the matched time window; uncertainty from regression stderr |
| Matched windows | Classical run uses the same [t_start, t_end] as its paired WP run |
| v1 vs v2 | Use v2 wherever available; flag v1 points with open/hollow markers |
| σ=8 flag | Mark as boundary-compromised (only 3σ rear clearance, no v2 exists) |

---

## Data inventory

### L=50 standard density (r_s ≈ 5.69)

#### WP σ=1 energy sweep

| E (eV) | Run directory | v2? | 4σ clear? | has mom stats | has state_E |
|---|---|---|---|---|---|
| 20 | `run_wp_n162_L50_E20_sigma1_v2` | Y | ✓ | Y | Y |
| 25 | `run_wp_n162_L50_E25_sigma1_v2` | Y | ✓ | Y | Y |
| 50 | `run_wp_n162_L50_E50_sigma1_v2` | Y | ✓ | Y | Y |
| 100 | `run_wp_n162_L50_E100_sigma1_v2` | Y | ✓ | Y | Y |
| 200 | `run_wp_n162_L50_E200_sigma1_v2` | Y | ✓ | Y | Y |
| 300 | `run_wp_n162_L50_E300_sigma1_v2` | Y | ✓ | Y | Y |

#### WP σ=5 energy sweep

| E (eV) | Run directory | v2? | 4σ clear? | has mom stats | has state_E |
|---|---|---|---|---|---|
| 20 | `run_wp_n162_L50_E20` | N (v1) | ✓ (σ=5, launch=-5, clear=20) | Y | Y |
| 25 | `run_wp_n162_L50_E25` | N (v1) | ✓ | Y | Y |
| 50 | `run_wp_n162_L50_E50_v2` | Y | ✓ | Y | Y |
| 100 | `run_wp_n162_L50_E100_v2` | Y | ✓ | Y | Y |
| 300 | `run_wp_n162_L50_E300_v2` | Y | ✓ | Y | Y |
| 600 | `run_wp_n162_L50_E600_v2` | Y | ✓ | Y | Y |

#### WP supplementary at E=100

| σ (Bohr) | Run directory | v2? | 4σ clear? | Note |
|---|---|---|---|---|
| 0.5 | `run_wp_n162_L50_E100_sigma0p5` | N | ✓ | |
| 3 | `run_wp_n162_L50_E100_sigma3` | N | ✓ | |
| 8 | `run_wp_n162_L50_E100_sigma8` | N | ✗ (3σ only) | Flag in plot |

#### Classical L=50

| E (eV) | Run directory | v2? | Note |
|---|---|---|---|
| 20 | `run_classical_n162_L50_E20` | N | v1 only |
| 25 | `run_classical_n162_L50_E25` | N | v1 only |
| 50 | — | — | **MISSING** (incomplete) |
| 100 | `run_classical_n162_L50_E100_v2` | Y | |
| 200 | — | — | **MISSING** (never run) |
| 300 | — | — | **MISSING** (v2 incomplete) |
| 600 | `run_classical_n162_L50_E600_v2` | Y | |

### L=30 high density (r_s ≈ 3.41)

#### WP σ=1

| E (eV) | Run directory | v2? |
|---|---|---|
| 50 | `run_wp_n162_L30_E50_highdens_sigma1_v2` | Y |
| 100 | `run_wp_n162_L30_E100_highdens_sigma1_v2` | Y |
| 200 | `run_wp_n162_L30_E200_highdens_sigma1_v2` | Y |
| 300 | `run_wp_n162_L30_E300_highdens_sigma1_v2` | Y |

#### WP σ=0.5 at E=100

| σ | Run directory | v2? |
|---|---|---|
| 0.5 | `run_wp_n162_L30_E100_highdens` | N |

#### Classical L=30

| E (eV) | Run directory | v2? |
|---|---|---|
| 50 | `run_classical_n162_L30_E50_highdens` | N |
| 100 | `run_classical_n162_L30_E100_highdens` | N |
| 200 | `run_classical_n162_L30_E200_highdens` | N |
| 300 | `run_classical_n162_L30_E300_highdens` | N |

---

## Computation pipeline

### Step 1: Compute interference-free time window per WP run

For each WP run, read `run_summary.txt` to extract:
- `wp_center_bohr` → launch position z₀
- `wp_sigma_bohr` → σ
- `wp_k0_bohr_inv` → k₀ (z-component)
- `cell_bohr` → L
- `dt_au`, `rt_num_steps`

Compute:
```
σ_r(t) = sqrt(σ²/2 + t²/(2σ²))          # density width
z_centroid(t) = z₀ + k₀ × t               # centroid position
z_leading(t) = z_centroid(t) + 4 × σ_r(t)  # 4σ leading edge
t_end: solve z_leading(t_end) = L/2         # boundary hit
```

Use `scipy.optimize.brentq` to solve for t_end.

Output: `{run_name: (t_start=0, t_end, n_steps_end)}` dictionary.

### Step 2: Extract WP stopping power — Definition 1 (momentum)

For each WP run, read `wp_momentum_stats.csv`:
- Columns: `step, time_au, pz_mean, sigma_pz2, e_kin_ha`
- Compute `E_kin(t) = pz_mean² / 2 + sigma_pz2 / 2` (atomic units, m=1)
- Convert to eV: `E_kin_eV = E_kin × 27.211`

Stopping power:
```
ΔE = E_kin(t_end) - E_kin(0)                    # energy change (should be negative)
Δz = z_centroid(t_end) - z_centroid(0)           # distance traversed
S₁ = -ΔE / Δz                                   # eV/Bohr (positive = energy loss)
```

### Step 3: Extract WP stopping power — Definition 2 (KS orbital energy)

For each WP run, read `state_energies.csv`:
- Find the WP state (state_index = `wp_state_index` from run_summary, typically 100)
- Extract eigenvalue at t=0 and t=t_end
- `ΔE_KS = ε_WP(t_end) - ε_WP(0)`
- `S₂ = -ΔE_KS / Δz`

**Caveat**: state_energies.csv format needs verification. Check column names and whether it stores time-resolved eigenvalues or just initial.

### Step 4: Extract classical stopping power

For each classical run, read `observables.csv`:
- Columns include: `step, time_au, energy_total, ...`
- The total energy of the electronic system changes as the classical point charge transfers energy
- Need to identify which column tracks the relevant energy (likely `energy_total` or `energy_kinetic`)

Also need the classical particle position z(t) — this may be in a separate file or computed from the Ehrenfest trajectory.

**Key**: use the same time window [0, t_end] as the paired WP run at the same energy. If no paired WP run exists (e.g. E=600 classical has no σ=1 WP), use the WP σ=5 window.

Stopping power from linear regression:
```
E_system(z) over the window → linear fit → slope = dE/dz = S_classical
stderr of slope → uncertainty bar
```

### Step 5: Verify window matching

For each (WP, classical) pair at the same energy:
- Confirm the classical run's total time ≥ t_end of the WP
- If classical run is shorter, truncate to the classical run's duration and note this

### Step 6: Generate the plot

#### L=50 version (`fig_master_stopping.png`)

Layout: single panel, log-log.

Data series:
1. **Classical Ehrenfest** (L=50): black squares connected by line, error bars from regression. Open squares for v1-only points.
2. **WP σ=1**: red circles, v2 filled, v1 open
3. **WP σ=5**: blue triangles, v2 filled, v1 open
4. **WP supplementary** (σ=0.5, 3, 8 at E=100): distinct markers (diamond, pentagon, hexagon), σ=8 with a warning hatch or cross
5. **Bethe–Bloch asymptote**: grey dashed v⁻² line

Two sub-versions: one using S₁ (momentum definition), one using S₂ (KS orbital definition).

#### L=30 version (`fig_master_stopping_highdens.png`)

Same layout but for L=30 data. Only σ=1 WP sweep + σ=0.5 supplementary + classical line.

---

## Provenance requirements

The plotting script must print, for each data point:
1. Run directory used
2. v1 or v2
3. Time window [t_start, t_end] in both a.u. and fs
4. E_kin(0), E_kin(t_end), ΔE, Δz, S for each definition
5. Whether the 4σ boundary clearance is satisfied at launch

This output goes to stdout when the script runs, providing full traceability.

---

## Files to create

| File | Purpose |
|---|---|
| `inq-stack/python/inqview/report1/stopping_power_data.py` | Computation module: reads CSVs, computes windows, extracts S |
| `inq-stack/python/inqview/report1/fig_master_stopping.py` | Plotting script using the computation module |
| `docs/reports/report1/figures/fig_master_stopping.png` | L=50 plot, S₁ definition |
| `docs/reports/report1/figures/fig_master_stopping_ks.png` | L=50 plot, S₂ definition |
| `docs/reports/report1/figures/fig_master_stopping_highdens.png` | L=30 plot, S₁ definition |
| `docs/reports/report1/figures/fig_master_stopping_highdens_ks.png` | L=30 plot, S₂ definition |

---

## Missing data — runs needed

| Run | Energy | Density | Priority |
|---|---|---|---|
| `run_classical_n162_L50_E50_v2` | 50 eV | L=50 | High — gap in classical line |
| `run_classical_n162_L50_E200` | 200 eV | L=50 | High — no classical at this energy |
| `run_classical_n162_L50_E300_v2` | 300 eV | L=50 | High — v2 incomplete |
| `run_classical_n162_L50_E25_v2` | 25 eV | L=50 | Medium — v1 exists |
| `run_classical_n162_L50_E20_v2` | 20 eV | L=50 | Medium — v1 exists |
| `run_wp_n162_L50_E100_sigma8_v2` | 100 eV | L=50 | Low — supplementary point, needs launch correction |

---

## Dependent plots to rebuild from this pipeline

### `fig_matched_pair` — WP–classical matched-pair comparison

**Current state**: `/local/data/public/skcb2/tddft/docs/reports/report1/figures/fig_matched_pair.png`
**Script**: `inq-stack/python/inqview/report1/fig_matched_pair.py`

**What's wrong with the current version**:
- Panel (a) shows S vs v/v_F with only 3 classical + 2 WP points — will have many more from the pipeline
- Panel (b) shows |ΔE_WP| vs σ — this belongs in `fig:sigma-sweep`, not here
- The plan description specifies panel (b) should be **Hartree energy time series** ΔE_H(t) comparing classical (positive, Coulomb wake) vs WP (negative, depletion anti-wake) at E=100 eV
- Stopping power values were computed with unknown time windows (pre-pipeline)

**Changes needed**:
1. **Panel (a)**: S(v) or S(E) from the pipeline, using the interference-free window S₁ definition. Show classical line (black, all energies), WP σ=1 (red), WP σ=5 (blue). Use v2 runs, flag v1 with open markers. Add Bethe–Bloch asymptote. x-axis as E (eV) rather than v/v_F for consistency with the master plot.
2. **Panel (b)**: Replace σ-sweep with Hartree energy time series. Data source: `observables.csv` column `energy_hartree` (or equivalent) for both the classical run and the σ=1 WP run at E=100 eV. Plot ΔE_H(t) = E_H(t) − E_H(0) vs t (fs). Classical should be positive (the point charge creates a Coulomb wake that increases Hartree energy). WP should be negative (the diffuse WP creates a depletion that decreases Hartree energy). This sign opposition is a key physical result.
3. **Time axis in fs** throughout.
4. All stopping power values from the `stopping_power_data.py` module (shared with master plot).

### `fig_definition_comparison` — Stopping power definition bar chart

**Current state**: `/local/data/public/skcb2/tddft/docs/reports/report1/figures/fig_definition_comparison.png`
**Script**: `inq-stack/python/inqview/report1/fig_definition_comparison.py`

**Changes needed**:
- Rebuild using S₁ (momentum) and S₂ (KS orbital) from the pipeline
- Add the classical reference value with uncertainty bar
- Two sub-panels: (a) σ=1 at E=100 (definitions should converge), (b) σ=5 at E=100 (definitions should diverge)
- All values traceable to the pipeline's provenance output

### `fig_sigma_sweep` — σ-dependence at E=100

**Current state**: existing plot not yet reviewed in detail
**Changes needed**:
- Panel (a): |ΔE_WP| vs σ on log-log (the content currently misplaced in fig_matched_pair panel b)
- Panel (b): KL divergence or other momentum-space coupling metric vs σ
- Use pipeline-computed values with proper interference-free windows

### `fig_sigma_rs` — σ/r_s scaling collapse

**Changes needed**:
- |ΔE_WP|/|ΔE_cl| vs σ/r_s using pipeline values at both densities (L=50 and L=30)
- Both densities should collapse onto a single curve if σ/r_s is the controlling parameter

---

## Validation checks

1. At σ=1 (σ ≪ r_s ≈ 5.69): WP stopping power should approach classical within ~5%
2. At σ=5 (σ ≈ r_s): WP stopping power should be substantially lower than classical
3. S₁ and S₂ definitions should agree at small σ and diverge at large σ
4. Classical S should decrease monotonically with E at high E (Bethe–Bloch regime)
5. All v2 runs at σ=5 should give different S than their v1 counterparts (boundary correction matters)
6. Energy conservation: ΔE_kin(WP) + ΔE_system(bath) ≈ 0 (check this as a sanity test)
