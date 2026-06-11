---
name: tddft-simulations
description: End-to-end workflow for planning, configuring, running, and post-processing INQ TDDFT simulations. Covers jellium WP/classical, coronene LEED, and ionic kick systems. Handles observable selection, GPU scheduling, ground-state validation, pilot evaluation, transient-period exclusion, and automated dispatch with Gmail notifications.
---

# TDDFT Simulation Skill

Use when the user asks to run a TDDFT simulation, set up a propagation run,
launch a parameter sweep, or post-process simulation results in any system
under `ResearchProject/systems/`.

## When to use

- "Run a simulation at E=200 eV on jellium"
- "Set up a coronene LEED run with different impact parameter"
- "Launch the sigma sweep at these energies"
- "Post-process the results from run_wp_..."
- "Make a stopping power comparison plot"
- Any request to create run.cpp, configs, analyse.py, or dispatch scripts

## When NOT to use

- Creating a brand-new material system from scratch (use `build-run` + `simulation-validation`)
- Modifying INQ source code (use standard development workflow)
- Pure literature/theory questions (use `literature-review`)

---

## Phase 0: Identify simulation type

Ask the user:

> What type of TDDFT simulation?
> 1. **Jellium wave-packet** — Gaussian WP projectile in periodic jellium bath
> 2. **Jellium classical** — Ehrenfest point-electron projectile in jellium
> 3. **Coronene LEED** — WP scattering off coronene with LEED screen accumulation
> 4. **Ionic kick** — Delta-kick or continuous perturbation of an ionic system
> 5. **Free-WP validation** — Non-interacting WP propagation (propagator test)

Record the type. All subsequent phases branch on this.

---

## Phase 1: System and ground-state validation

### 1a. Locate system infrastructure

Check that the target system directory exists under `ResearchProject/systems/<system>/`:
- `shared/configs/` — at least one base config `.hpp`
- `shared/cpp/run_template.hpp` — propagation template
- `checkpoints/` or `save_gs/` — at least one GS checkpoint

If missing, stop and tell the user:
> "This system doesn't have the infrastructure for TDDFT runs yet.
> Use the `build-run` skill to set up the ground state first."

### 1b. Validate ground state

Read the GS checkpoint's `results/run_summary.txt`. Verify:

| Check | Criterion | Action if fail |
|-------|-----------|----------------|
| SCF converged | `ground_state_energy_ha` present and non-NaN | **HARD STOP** |
| SCF tolerance | tol ≤ 1e-4 Ha (preferably 1e-6) | **HARD STOP** if > 1e-4 |
| Electron count | Matches expected N for the system | **HARD STOP** |
| Closed-shell (jellium) | N is a magic number (2, 14, 38, 54, 66, 114, 162, 226...) | **WARN** if not |

### 1c. Check for existing GS checkpoints

List all `save_gs/` directories. If multiple exist, present them to the user
with their (L, N, dx) parameters. If the user's requested parameters match
an existing GS, reuse it. If not, flag that a new GS is needed (defer to
`build-run` skill).

---

## Phase 2: Configuration planning

### 2a. Collect run parameters from user

Required parameters (type-dependent):

| Parameter | Jellium WP | Jellium Classical | Coronene | Ionic Kick |
|-----------|-----------|------------------|----------|------------|
| Energy (eV) | ✓ | ✓ | ✓ | — |
| σ (Bohr) | ✓ | — | ✓ | — |
| Cell L (Bohr) | ✓ (default from GS) | ✓ | ✓ | ✓ |
| Grid dx (Bohr) | ✓ (default from GS) | ✓ | ✓ | ✓ |
| N_electrons | ✓ (from GS) | ✓ | from geometry | from geometry |
| dt (a.u.) | default 0.020 | default 0.020 | default 0.020 | default 0.020 |
| Impact offset | — | — | ✓ (WP_CX, CY, CZ) | — |
| Kick direction/strength | — | — | — | ✓ |

### 2b. Nyquist validation (PLANNING PHASE — before any build)

For WP simulations, compute:
```
k₀ = sqrt(2 × E_eV / 27.211)
σ_k = 1 / (σ × sqrt(2))
k_max = k₀ + 3 × σ_k
dx_max = π / k_max
```

**If dx > dx_max: HARD STOP.** Report:
> "Grid spacing dx={dx} exceeds Nyquist limit dx_max={dx_max:.3f} for
> E={E} eV, σ={σ} Bohr. Either reduce dx (requires new GS) or reduce E."

Display the Nyquist table for the user showing all requested energies.

### 2c. Memory estimation (PLANNING PHASE)

Compute estimated GPU memory:
```
grid_points = (L / dx)³
n_states = N_electrons / 2 + extra_states
orbital_memory_GB = n_states × grid_points × 16 / 1e9
total_estimate_GB = orbital_memory_GB × 2.5  (safety factor for work arrays)
```

| Threshold | Action |
|-----------|--------|
| > 20 GB | **HARD STOP** — will not fit on A30 (24 GB) |
| > 16 GB | **WARN** — tight fit, may deadlock at init |
| ≤ 16 GB | OK — proceed |

### 2d. Compute simulation duration

**For WP runs:** Compute self-spread-capped N_STEPS:
```python
# Boundary rule: launch_z = -L/2 + 4σ
# Self-spread cap: solve launch_z + v·t + 3·σ_density(t) = +L/2
# where σ_density(t) = (σ/√2)·√(1 + (t/σ²)²)
# N_STEPS = min(self_spread_cap, boundary_rule_cap)
```

For σ ≤ 2 Bohr, self-spread usually dominates. For σ ≥ 3, boundary rule dominates.

**For classical runs:** Use boundary rule directly:
```
traversal = L - 5σ  (or L - 4σ for relaxed rule at large σ)
t_total = traversal / v
N_STEPS = ceil(t_total / dt)
```

**For coronene:** Use `compute_n_steps(Lz, offset, sigma, k0, dt)` from the base config.

### 2e. WRITE_EVERY and frame cadence

Target ~300 frames for VTI series:
```
WRITE_EVERY = max(1, round(N_STEPS / 300))
```

Report estimated frames: `N_STEPS / WRITE_EVERY`.

---

## Phase 3: Observable selection

<!-- min-obs-set: canonical = inq-stack/include/inqkit/observables/minimum_observable_set.hpp -->
> **Canonical required set (Cluster O, ADR 0006).** The *required* observables
> per run-type are defined ONCE in
> `inq-stack/include/inqkit/observables/minimum_observable_set.hpp`, which the
> run writes as `results/observables_manifest.json` at startup and `validate_run`
> checks post-run. The Tier 1/2 tables below are the **operational view** of that
> set — cadence, writer, and failure-modes — NOT a second source of truth; they
> must stay consistent with the `.hpp` (the Cluster-O drift eval enforces it).
> Tier 3 is run-specific **optional** extras beyond the required set.

### 3a. Auto-enable Tier 1 (always-on, all types)

| Observable | Output | Writer |
|------------|--------|--------|
| observables.csv | energy_total/kinetic/hartree/xc, current_xyz, dipole_xyz | ObservablesWriter |
| eigenvalues/ | GS KS eigenvalues + occupations | eigenvalue_dump |
| density_delta L2 | Integrated density fluctuation σ²ₙ(t) | DensityDelta |
| run_summary.txt | Run metadata (stub at start, final at end) | direct ofstream |

### 3b. Auto-enable Tier 2 (standard for simulation type)

**Jellium WP:**
- density_total + density_system VTI (per WRITE_EVERY, bulk-copy)
- **density_wp VTI (COMPULSORY, per WRITE_EVERY — EQUAL CADENCE)** — the WP
  orbital density |ψ_wp|². MUST be written at the SAME cadence as
  density_total (`wf_write_every == write_every`) for any run intended for
  wake / density-difference analysis. The canonical bath density is
  `n_system = n_total − n_wp` computed in post-processing (the saved
  `density_system` field is NOT run-independent — see §7.0), and that
  subtraction is only exact when density_wp exists at the SAME step as
  density_total. **Failure mode (do not repeat):** the `_v2` jellium runs
  saved density_wp ~10× coarser than density_total (e.g. σ1_v2: 32 wp vs
  317 total frames). Nearest-step subtraction across that gap leaves a
  charge-neutral MOVING-WP DIPOLE residual (~85% of the σ=1 wake signal,
  2026-06-01) and forces wake plots to be sampled only at the sparse
  wp-frame times. Original (pre-template-fix) runs saved 0 wp frames →
  unusable for any wake/difference analysis. Equal cadence is non-negotiable
  when a run is meant to feed `inqview.postprocess.wake`.
- density_delta VTI (raw + coarse, per WRITE_EVERY)
- state_energies.csv (per 5×WRITE_EVERY)
- occupations_vs_time.csv (per 5×WRITE_EVERY)
- momentum_distribution.csv (per 10×WRITE_EVERY)
- orbital_overlap WP-only (per 10 steps)
- orbital_overlap proxies (per 5×WRITE_EVERY)
- orbital_overlap full (t=0 and t=final only)
- wp_momentum_stats.csv (per WRITE_EVERY)

**Jellium Classical:**
- Same density/state_energies/occupations/momentum as WP
- electron_track.csv (every step — projectile pos/vel/force)
- orbital_overlap full (t=0, t=mid, t=final)
- orbital_overlap proxies (per 5×WRITE_EVERY)
- NO wp_momentum_stats, NO wp_real_space_stats

**Coronene LEED:**
- density_rt_system + density_rt_wp + density_rt_total VTI (per WRITE_EVERY)
- orbital_overlap WP-only (every step)
- LEED screen accumulators (full-time + per-screen-window + paper-window)
- LEED instantaneous snapshots (per SCREEN_SNAP_EVERY)
- wp_momentum_distribution (per EVERY 2 timesteps)
- state_energies

**Ionic Kick / Free-WP:**
- density_total VTI (per WRITE_EVERY)
- Minimal: observables.csv + density_delta

### 3c. Ask about Tier 3 (user-requested extras)

Present to the user:

> The following optional observables are available. Select any you need:
> - [ ] WP real-space stats (σ_x, σ_y, σ_z evolution — requires gamma-only k-point)
> - [ ] GS orbital density VTIs (one VTI per occupied orbital at t=0 — large output)
> - [ ] WP initial wavefunction VTI (complex field, ~500 MB for large grids)
> - [ ] density_delta coarse VTI (3 Bohr bins — cheaper than raw, good for visualisation)
> - [ ] Custom observable (describe what you need)

### 3d. Suggest additional observables based on aims

If the user's stated aim involves:
- "stopping power" → ensure electron_track (classical) or wp_momentum_stats (WP)
- "LEED" or "diffraction" → ensure LEED screens enabled
- "plasmon" or "collective" → suggest density_fourier (axial modes)
- "secondary electrons" → flag as not yet implemented; suggest density_delta monitoring
- "convergence test" → suggest extra_states sweep (x20/x40/x80)
- "loss function" or "S(q,omega)" → ensure density VTIs saved at sufficient cadence;
  requires long propagation (T >> 2*pi/omega_p for spectral resolution).
  For r_s=5.69: omega_p=3.47 eV → need T > 10 a.u. minimum.
  For r_s=3.41: omega_p=7.49 eV → need T > 5 a.u. but the L=30 box
  limits trajectory to ~5 a.u. at E=100 eV. Use L=50 runs for spectral analysis.
- "e-h transitions" or "GS decomposition" → ensure orbital_overlap full at t=0
  and t=final (uses GS-projected occupation analysis). Post-processed by
  inqview `overlap` phase producing fig_gs_decomposition.png — a bar chart
  showing per-orbital occupation changes and charge-balance accounting.
  **MANDATORY for any GS-occupations / GS-decomposition graph — save the
  orbitals, not just the overlap scalars:**
  1. **Ground-state KS orbital wavefunctions** (the full reference basis, the
     complex fields, at t=0).
  2. **Final-state KS orbital wavefunctions** (all evolved orbitals ψ_j(t_end),
     the complex fields — NOT only the WP orbital).
  3. **The overlap matrix** itself (`orbital_overlap full` at t=0 and t=final).
  **Why:** `OrbitalOverlapMatrix` computes `|⟨ψ_i^GS|ψ_j(t)⟩|²` in-loop from GPU
  memory and writes ONLY the squared scalars onto the *truncated* `n_ref`
  computed GS basis (default 100). It never persists the orbitals. From the
  squared, basis-truncated overlaps you CANNOT (a) recover the complex orbitals
  (phase is lost), nor (b) re-project the evolved orbitals onto *new / analytic*
  GS orbitals above the cutoff. So the "grey region" (charge excited above
  orbital `n_ref`) is unrecoverable in post — it requires a full re-run.
  Persisting the GS + final KS orbital wavefunctions (or running with `n_ref`
  raised to cover the expected excitation, e.g. ~200 for jellium) is the only
  way to (i) close the charge discrepancy to zero and (ii) extend the
  decomposition onto analytically-known high-energy orbitals afterwards.
  (Discovered 2026-06-01 on `run_wp_n162_L50_E25_sigma1_v2`: only the WP
  orbital was saved in real space, so the +0.497 e bath e-h excitation above
  orbital 100 could not be filled without re-propagating. See
  `docs/handovers/mphil_midterm_presentation.md`.)

### 3e. Post-processing observable catalogue

Key diagnostic plots produced by `analyse.py` and their purposes:

| Observable | Plot | What it shows |
|------------|------|---------------|
| GS decomposition | `fig_gs_decomposition.png` | Per-orbital occupation change at t_end. Bar chart: occupied (red, deplete), unoccupied (blue, gain), with charge-balance sidebar. |
| Energy decomposition | `fig_energy_decomposition.png` | (a) Total-system energy components vs time; (b) WP-specific KE decomposition (<p>^2/2m, sigma_p^2/2m). Critical: panel (a) is total system, panel (b) is WP-only. |
| Loss function | `fig_loss_function.png` | 2D heatmap S(q,omega) with Lindhard boundaries and plasmon dispersion overlaid. Requires density VTI time series. Frequency resolution = 2*pi/T_propagation. |
| DD1 density diff | `fig_DD1_density_diff_grid.png` | Bath response isolated: n_jellium(r,t) - n_free(r,t). Requires matched free-WP companion run. |
| Density z-profiles | `fig_density_profile_comparison.png` | 1D z-profiles comparing WP-in-jellium against free WP at multiple timesteps. |
| Plasmon FFT | `fig_plasmon_fft.png` | Fourier spectrum of density oscillation time series. Identifies plasmon mode at omega_p. |
| Momentum before/after | `fig_momentum_before_after.png` | WP momentum distribution at t=0 and t=end. Shows deceleration and broadening. |

---

## Phase 4: Config and run.cpp creation

### 4a. Create config header

Follow the naming convention:
- **Jellium:** `shared/configs/{base_description}_{energy}_{cell}[_{variant}].hpp`
  - Example: `electron_proj_E50_L50_cubic_sigma1.hpp`
  - Example: `highdens_n162_L30_E200_sigma1.hpp`
- **Coronene:** `shared/configs/{descriptor}_{cell}.hpp`
  - Example: `cc_bond_35x35x60.hpp`

Each config header MUST include:
- Detailed physics comment explaining the run's purpose
- Nyquist calculation in the header comment
- Self-spread or boundary-rule N_STEPS derivation
- For WP: both WP and Classical struct variants (twin Common_-derived structs)

Inherit from the system's base config. Override only changed parameters.

### 4b. Create run directory

Naming convention:
```
run_{type}_{N_electrons}_{cell_spec}_{energy}[_{variant_tags}]/
```

Examples:
- `run_wp_n162_L50_E50_sigma1/`
- `run_classical_n162_L30_E200_highdens/`
- `run_cc_bond/`

### 4c. Create run.cpp

Use the system's existing run.cpp template as the base. Change ONLY:
- Config `#include` path
- `using Cfg = ...` typedef
- `RUN_NAME` string
- `GS_DIR` path
- Header comment

**Do NOT modify callback logic, observable writers, or output structure.**

### 4d. Create analyse.py

Copy from the system's canonical analyse.py template:
- **Jellium WP:** `run_wp_n162_L50_E100_sigma1/analyse.py`
- **Jellium Classical:** `run_classical_n162_L30_E100_highdens/analyse.py`
- **Coronene:** `run_cc_bond/analyse.py`

Change ONLY the docstring/run name. The template includes:
- `install_schema_shims()` — symlink bridge for path schema mismatches
- `analyse_extras` import — windowed stopping power, overlap heatmaps, GS-projected occupations
- `density_fourier` — axial density mode detection
- `physics_block()` — energy budget summary
- `render_markdown()` → `REPORT.md`

---

## Phase 5: Pilot run (if new config)

### 5a. Determine if pilot is needed

Pilot is **mandatory** if:
- First run with this (system, energy range, σ, grid spacing) combination
- First run with a new GS checkpoint
- User explicitly requests it

Pilot is **skipped** if:
- A nearby parameter point has already been validated (e.g. E=100 validated, now running E=200)
- User explicitly says "skip pilot, I trust this config"

### 5b. Run single simulation

Launch one run on a free GPU:
```bash
CUDA_VISIBLE_DEVICES=<free_gpu> inq-run
```

**Always use the venv Python for post-processing:**
```bash
/local/data/public/skcb2/tddft/venv/bin/python3 analyse.py
```

### 5c. Present diagnostic plots to user

After the pilot completes, present these key diagnostics:

1. **Energy drift** — `total_energy_vs_time.png` (conservation check)
2. **All energy components** — `all_energies_vs_time.png`
3. **Density fluctuation** — `density_fluctuation_l2.png` (should be smooth, no spikes)
4. **WP centroid position** — `wp_position_vs_time.png` (should be monotonic for single-pass)
5. **Eigenergies over time** — `ks_energies_absolute.gif` or first frame

If Gmail notification is requested, email these plots to the user via:
```python
from inqview.email import send_run_email
send_run_email(subject, body, attachments=[...], to="chiddukanna@gmail.com")
```

### 5d. Post-run validation checks

| Check | Criterion | Action |
|-------|-----------|--------|
| Run completed | `run_completed = true` in run_summary.txt | **HARD STOP** if false |
| WP norm | `norm_after` ∈ [0.95, 1.05] | **HARD STOP** if outside |
| Energy drift | \|ΔE_total\| < 1 mHa over trajectory | **WARN** if exceeded |
| max_overlap | < 0.01 after orthogonalisation | **WARN** if exceeded |

Wait for user approval before proceeding to batch.

---

## Phase 6: Batch dispatch

### 6a. GPU queue-based scheduling

Poll GPUs via nvidia-smi. Launch next job on first free GPU:

```python
def find_free_gpu() -> int | None:
    """Return GPU ID with < 500 MB used, or None."""
    result = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used",
         "--format=csv,noheader,nounits"],
        capture_output=True, text=True)
    for line in result.stdout.strip().split("\n"):
        idx, mem = line.split(",")
        if int(mem.strip()) < 500:
            return int(idx.strip())
    return None
```

Launch each run with:
```bash
CUDA_VISIBLE_DEVICES=<gpu_id> inq-run
```

After each run completes:
1. Run `analyse.py` with venv Python
2. **Update the run catalogue** (see `tddft-run-catalogue` skill):
   `venv/bin/python3 .claude/skills/tddft-run-catalogue/scan_runs.py --run <run_dir>`
3. Send Gmail notification if requested
4. Check for next queued run

### 6b. Environment setup

Every subprocess MUST have:
```python
env = {
    **os.environ,
    "CUDA_VISIBLE_DEVICES": str(gpu_id),
    "PATH": "/local/data/public/skcb2/tddft/shared/bin:" + os.environ["PATH"],
    "INQ_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share",
    "PSEUDOPOD_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
}
```

### 6c. Post-processing

**Always use the venv Python:**
```
/local/data/public/skcb2/tddft/venv/bin/python3 analyse.py
```

Never use system Python — it lacks VTK/pyvista dependencies.

---


# TODO: This is to be linked with the minimum observable set, and the required analysis can be directly implemented using deterministic hooks
## Phase 7: Post-processing and analysis

### 7.0 Canonical bath density + plot-comparison rules (MANDATORY)

**Canonical bath ("system") density — always use `total − wp`.**
The jellium bath/system density is
`n_system(r,t) = n_total(r,t) − n_wp(r,t)` (full electronic density minus the
injected WP orbital, occupation 1.0); the induced wake is
`Δn_system(r,t) = n_system(r,t) − n_system(r,t0)`.
DO NOT trust the saved `density_system` VTI field for the bath — it is
INCONSISTENT across run generations (verified 2026-06-01 by integration):
old runs save `density_system` = WP-INCLUDED (163 e for N=162), new `_wf`
runs save it bath-only (162 e). Only `total − wp` is run-independent
(N_total=163 = 162 bath + 1 WP, N_wp=1.000, total−wp=162=N_e, conserved at
all t). Classical runs have no WP orbital → `n_system = n_total` directly.
Helper + loaders: `inqview.postprocess.wake` (`bath_volume`, `bath_line_z`,
`bath_slice_xz`, `wp_centroid_z`). WP centroid for overlays comes from the
density_wp first moment (or `wp_real_space_stats.csv`).

**Shared-colorbar rule (directly-compared plots).** Panels meant to be
compared directly (e.g. WP wake vs classical wake) MUST share ONE identical
colour scale — compute it once over all such panels (`wake.shared_clim`).
A difference panel (WP − classical) is an independent quantity → its OWN
colour scale. Provide BOTH linear and log (symlog for signed) views.

#### Recipe: induced-wake difference plots (the canonical 3 steps)

To compare how two systems (e.g. a WP projectile vs a classical electron)
perturb the bath, build the plot in EXACTLY these three steps. Helpers live in
`inqview.postprocess.wake`.

**Step 1 — `n_system` at each timestep (EXACT same-timeframe subtraction).**
For a WP run: `n_system(r,t) = n_total(r,t) − n_wp(r,t)`, with n_total and n_wp
read from the SAME exact step (`wake.bath_volume`, which snaps to a frame that
has an exact density_wp partner and raises otherwise). NEVER subtract a
density_wp from a different (nearest) step: the WP is a moving Gaussian, so a
mismatched-step subtraction leaves a charge-neutral MOVING-WP DIPOLE residual
that looks like the WP is still in the wake (caught 2026-06-01; was ~85% of the
σ=1 signal). When density_wp is saved sparsely (v2 runs: ~10× coarser than
density_total), SAMPLE THE PLOT AT `wake.wp_frame_times(run)` so every frame has
an exact partner. For a CLASSICAL run there is no WP orbital → `n_system = n_total`.
Verify: `∫ n_system dV = N_electrons` (162 here) at every frame.

**Step 2 — induced wake `Δn_system(t) = n_system(t) − n_system(t0)`.**
Subtract the t=0 bath from every frame (`wake.bath_line_z` / `bath_slice_xz`
give the 1D z-profile and 2D xz slab; subtract the t0 result). Known-case:
`Δn_system(t0) = 0` exactly. This isolates the bath's RESPONSE to the projectile.

**Step 3 — cross-system metric difference at a MATCHED timestep.**
To isolate what is unique to one system, difference the two induced wakes at the
SAME physical time t: `D(r,t) = Δn_system^A(r,t) − Δn_system^B(r,t)` (e.g.
A = WP, B = classical). Both A and B must be evaluated at the same t (snap each
to its nearest available frame at that t). Display A and B in panels that SHARE
one colour scale (`shared_clim` over both); the difference D is an independent
quantity → its OWN colour scale. Mark the WP centroid (`wp_centroid_z`) so the
reader can tell co-moving screening from leading/trailing structure. Provide
linear AND symlog versions. Fix the colour/axis scale ONCE over all frames of an
animation — never per-frame.

Reference implementation: `ResearchProject/systems/jellium/scripts/wake_movie_driver.py`.

### 7a. Transient-period exclusion

**Auto-compute the interference-free window (IFW):**

For WP runs:
```
t_IFW = time when 3×σ_density(t) tail reaches far cell face
σ_density(t) = (σ/√2) × √(1 + (t/σ²)²)
```

For classical runs:
```
Δz window = [3, min(28, traversal-3)] Bohr from launch
Equivalent time: t ∈ [3/v, Δz_max/v]
```

**Show the computed window to the user:**
> "Transient exclusion window: t ∈ [{t_start:.2f}, {t_end:.2f}] a.u.
> (Δz ∈ [{z_start:.1f}, {z_end:.1f}] Bohr). Override? [Enter to accept]"

### 7b. Mandatory post-processing outputs (common evaluation base)

Every completed run MUST produce these outputs for consistent cross-run evaluation.
Verify their presence after analyse.py completes; re-run if missing.

**Group A — Time-series diagnostics (PNG):**
- `observables_summary.png` — 3-row panel: energy, current, dipole vs time
- `total_energy_vs_time.png` — energy conservation / drift check
- `all_energies_vs_time.png` — per-component (kinetic, Hartree, XC) evolution
- `current_components_vs_time.png` — J_x, J_y, J_z
- `dipole_components_vs_time.png` — μ_x, μ_y, μ_z
- `density_fluctuation_l2.png` — σ²_n(t) integrated fluctuation

**Group B — Energy decomposition (PNG, per-run + cross-run):**
- `energy_decomposition_classical_vs_wp.png` — 6-panel ΔE vs Δz (kinetic, Hartree, XC, total, bath sum, WP slot) with IFW shading. Ref: fig 06.
- `energy_bookkeeping_bar.png` — bar chart of ΔE_kinetic, ΔE_hartree, ΔE_xc at t_IFW. Ref: fig 07, fig C1.

**Group C — Spectral analysis (PNG):**
- `fft_total_energy.png`, `fft_current_{x,y,z}.png`, `dipole_spectrum_{x,y,z}.png`
- `spectra/{energy,current,dipole}/spectrum_*.png` — 4 detrending variants

**Group D — Density evolution (PNG + GIF):**
- `density/{total,system,wp}_{xy,xz,yz,z_profile}.gif` — density animations
- `density/{total,system}_{xy,xz,yz,z_profile}_log.gif` — log-scale
- `density_z_profile_evolution.png` — z-profile heatmap vs time. Ref: fig 09.
- `delta_density_xz_snapshots.png` — δn(x,y=0,z) at multiple times (lab frame). Ref: fig D1.
- `z_profile_diff_vs_free.png` — z-profile difference between jellium run and matched free-WP propagation
- `density_diff_vs_free.png` — δn between jellium and free propagation at selected times

**Group E — GS and eigenvalue structure (PNG):**
- `eigenvalues/eigenvalue_levels.png` — KS level diagram
- `eigenvalues/eigenvalues_dos.png` — density of states
- `ground_state/density_gs_system_xy.png` — GS density slice
- `ground_state/density_gs_z_profile.png` — GS z-profile
- `ground_state/gs_occupations.png` — occupation bars
- `layout/layout_xz.png` — cell geometry + WP launch layout
- `ks_eigenenergy_evolution.png` — KS orbital energies vs time (static version of GIF). Ref: fig 10.

**Group F — Orbital analysis (PNG + GIF):**
- `gs_basis_decomposition.png` — per-orbital Δn_i^GS at t_end (depletion below Fermi, excitation above). Ref: fig 11.
- `overlap_heatmap_log_wp.png` — |⟨ψ_i^GS|ψ_j(t_end)⟩|² heatmap, log scale, diagonal masked. Ref: fig 14.
- `overlap_heatmap_log_classical.png` — same for classical companion (if available).
- `overlap_heatmap_diff_wp_vs_classical.png` — difference between WP and classical overlap matrices.
- `ks_energies_absolute.gif`, `ks_energies_delta.gif` — eigenvalue time evolution
- `occupations_absolute.gif`, `occupations_delta.gif` — occupation dynamics

**Group G — Momentum and trajectory (PNG):**
- `wp_position_vs_time.png`, `wp_velocity_vs_time.png` — WP centroid tracking. Ref: fig 12.
- `momentum_band_free_vs_jellium.png` — ⟨p_z⟩ ± σ_p vs centroid z, free vs jellium. Ref: fig M_A.
- `momentum_distribution.gif`, `momentum_heatmap.png` — n(|k|, t)
- `sigma_xyz_vs_time.png` — σ_x(t), σ_y(t), σ_z(t) WP width evolution (from wp_real_space_stats)

**Group H — Stopping power and energy balance (PNG + CSV):**
- `knudsen_ke_vs_t.{png,csv}` — Knudsen KE stopping power (WP only)
- `kl_divergence_vs_t.{png,csv}` — momentum distribution divergence
- `energy_balance.{png,csv}` — conservation audit
- `classical_force_fixed.png` — F_z from dv/dt (classical only)
- `delta_E_total_vs_z.png` — windowed S ± SE (classical only)
- `running_slope_vs_z.png` — box-deficit diagnostic (classical only)
- `bath_energy_vs_time.{png,csv}` — bath energy (classical only)

**Group I — Advanced (v2 runs with WP wavefunction saving):**
- `wp_momentum_distribution_before_after.png` — |ψ̃_WP(k)|² before and after collision
- `momentum_difference_map_2d.png` — Δ|ψ̃(k_z, k_⊥)|² scattering map with elastic ring overlay. Reveals inelastic/elastic/backscatter channels. Verified 2026-05-25.
- `planewave_decomposition.png` — evolved KS orbitals decomposed into plane-wave basis (occupations mapped onto PW energies)
- `spectral_weight_response.png` — W_resp(q_z, ω) with exact WP subtraction
- `loss_function.png` — L(q_z, ω) = -(4π/q²) Im[χ]
- `secondary_electron_yield.png` — δ(t) in proxy vacuum region (if applicable)

**Markdown:**
- `REPORT.md` — auto-generated physics summary

### 7c. Type-specific post-processing outputs

**Jellium WP (additional to common base):**
- `wp_position_vs_time.png`, `wp_velocity_vs_time.png` — WP centroid tracking
- `ks_energies_absolute.gif`, `ks_energies_delta.gif` — eigenvalue time evolution
- `occupations_absolute.gif`, `occupations_delta.gif` — occupation dynamics
- `momentum_distribution.gif`, `momentum_heatmap.png` — momentum-space analysis
- `knudsen_ke_vs_t.{png,csv}` — Knudsen KE method stopping power
- `kl_divergence_vs_t.{png,csv}` — momentum distribution divergence
- `energy_balance.{png,csv}` — energy conservation audit
- `density/{delta}_{xy,xz,yz,z_profile}.gif` — density perturbation animations
- `overlap/overlap_heatmap_t_end.png` — WP-GS overlap matrix (analyse_extras)
- `gs_projected_occupations/*.png` — effective GS occupation evolution (analyse_extras)
- `n_q_m{1..6}.{csv,png}` — axial density Fourier modes (density_fourier)

**Jellium Classical (additional to common base):**
- `classical_force_fixed.png` — F_z = m·dv/dt recovered from electron_track
- `delta_E_total_vs_time.png` — bath energy gain with window markers
- `delta_E_total_vs_z.png` — windowed stopping power S ± SE (Δz ∈ [3, 28])
- `running_slope_vs_z.png` — box-deficit diagnostic
- `stopping_force_vs_z.png`, `dE_kinetic_vs_z.png` — stopping curve
- `bath_energy_vs_time.{png,csv}` — bath energy accumulation
- Same eigenvalue/occupation/momentum GIFs and overlap/GS-projected outputs as WP

**Coronene LEED (additional to common base):**
- `gs_orbital_gallery.png` — all orbital shapes at t=0
- `wp_position_vs_time.png` — WP centroid
- `screens/total/all_screens_grid.png` — LEED pattern overview
- `screens/total/screen_NN.png` + `_log.png` — per-screen patterns
- `screens/ifft/screen_NN_ifft_{amp,patterson}.png` — inverse FFT analysis
- `screens/time_windowed/*.png` — per-window patterns
- `screens/instantaneous/*.gif` — screen time evolution
- `overlap/wp_overlap_with_gs_orbitals.gif` — WP-GS overlap animation

### 7d. Standard pipeline phases by type

**Jellium WP:** summary, gs, layout, observables, eigenvalues_gs, gamma_transitions,
wp_trajectory, state_energies, state_energy_spectra, occupations, momentum,
knudsen_ke, kl_divergence, energy_balance, density, overlap, orbitals

**Jellium Classical:** summary, gs, layout, observables, eigenvalues_gs, state_energies,
state_energy_spectra, bath_energy, stopping, gs_projected_occupations,
occupations, momentum, energy_balance, density, overlap, orbitals
(+ analyse_extras: classical_force_fixed, delta_E_total_vs_z with windowed S)

**Coronene:** summary, gs, layout, observables, eigenvalues_gs, wp_trajectory,
state_energies, state_energy_spectra, occupations, momentum, density,
overlap, screens, orbitals

### 7e. Post-processing reference

Full observable catalogue (raw + post-processed, per simulation type):
`docs/observables/catalogue.md`

### 7f. Multi-run comparison scripts

Place in `scripts/` under the system directory. Common patterns:
- `plot_stopping_power_vs_energy.py` — S(v) across energy sweep
- `plot_sigma_sweep.py` — ΔE vs σ at fixed energy
- `run_comparison.py` — overlay observables from multiple runs

Results from multi-run analyses go in `hypotheses/<topic>/physics/`.

---

## Phase 8: File placement rules

| File type | Location |
|-----------|----------|
| Config headers | `ResearchProject/systems/<system>/shared/configs/<name>.hpp` |
| Run directories | `ResearchProject/systems/<system>/run_<type>_<params>/` |
| Run source | `run_<name>/run.cpp` |
| Per-run analysis | `run_<name>/analyse.py` |
| Run results | `run_<name>/results/{raw,analysis}/` |
| Shared C++ templates | `shared/cpp/run_template.hpp` |
| Shared Python analysis | `shared/python/analyse_extras.py` |
| Multi-run scripts | `scripts/<script_name>.py` |
| Hypothesis results | `hypotheses/<topic>/physics/*.png` |
| GS checkpoints | `save_gs/<gs_name>/` and `checkpoints/<gs_name>/` |
| Pseudopotentials | `shared/pseudopotentials/` |

---

## Naming conventions

### Run directories
```
run_{type}_{N_electrons}_{cell_spec}_{energy}[_{variant_tags}]/
```
- `type`: wp, classical, free_wp, plasmon
- `N_electrons`: n162, n138, n128
- `cell_spec`: L50, L30, L60 (cubic implied; add x dimensions for orthorhombic)
- `energy`: E50, E100, E200, E600 (in eV)
- `variant_tags`: sigma1, sigma0p5, highdens, v2, x40, knudsen, tilt45

### Config headers
```
{physics_description}_{energy}_{cell}[_{variant}].hpp
```
- Example: `electron_proj_E100_L50_cubic_sigma1.hpp`
- Example: `highdens_n162_L30_E200_sigma1.hpp`

### GS checkpoints
```
gs_{cell_spec}_{N_electrons}_{grid_spacing}[_{variant}]
```
- Example: `gs_L50_cubic_N162_dx0p40`
- Example: `gs_35x35x60_cut40`

---

## Critical: VTI coordinate mapping

INQ stores grids in **FFT-natural order** (index 0 = physical coordinate 0 = cell centre).
The `inqkit::RealField3DWriter` applies `fft_shift_index()` when writing VTIs, so
VTI files are in **physical order** (left-to-right: −L/2 to +L/2).

**However**, LEED screen `.dat` files are in FFT-natural order. When loading and
plotting screen data, ALWAYS apply `np.fft.fftshift(data)` before plotting.
Without this, diffraction peaks appear at array corners instead of the centre.

Full explanation: `docs/notes/coronene-geometry-correction.md`.

---

## Checklist (use as TodoWrite items)

- [ ] Identify simulation type
- [ ] Locate system infrastructure (configs, GS checkpoints, templates)
- [ ] Validate ground state convergence
- [ ] Collect run parameters from user
- [ ] Nyquist validation (PLANNING PHASE)
- [ ] Memory estimation (PLANNING PHASE)
- [ ] Compute N_STEPS (self-spread cap or boundary rule)
- [ ] Select observables (Tier 1+2 auto, Tier 3 user)
- [ ] Create config header
- [ ] Create run directory + run.cpp
- [ ] Create analyse.py from template
- [ ] Pilot run (if new config)
- [ ] Present diagnostic plots for user evaluation
- [ ] Post-run validation (norm, energy drift, conservation)
- [ ] User approval for batch
- [ ] Dispatch batch (queue-based GPU scheduling)
- [ ] Post-process each run (venv Python, analyse.py)
- [ ] Update run catalogue (`tddft-run-catalogue` skill, `scan_runs.py --run`)
- [ ] Gmail notification per run (if requested)
- [ ] Transient-period exclusion in analysis
- [ ] Multi-run comparison plots (if sweep)
- [ ] Update handover
