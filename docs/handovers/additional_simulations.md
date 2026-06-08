# Handover: Additional Simulations for Report 1

---

## Milestone: 2026-05-25 — skill expanded, v2 dispatch running, analysis plan locked

### Current status

**In progress.** 15 jellium v2 runs dispatching on 2 GPUs (queue-based). Coronene v2 queued (GPUs occupied). Skill updated with comprehensive mandatory observables (Groups A–I). All design decisions locked for post-processing.

### What changed

**LEED coordinate fix verified:**
- `rewindow_leed.py` now applies `np.fft.fftshift()` to screen `.dat` files before plotting
- Physical axes in Bohr added (extent = [-L/2, +L/2])
- Root cause documented: LEED `.dat` files are in FFT-natural order; VTI files are already in physical order (writer applies `fft_shift_index`)
- Reference: `docs/notes/coronene-geometry-correction.md`

**WP target orbital identified:**
- Target: C₁–C₂ radial spoke bond midpoint at (4.028, 0.0) Bohr
- Orbital 25 has highest density at this point (0.0215 e/Bohr³) — a HOMO-region π orbital
- Orbital density plots saved: `run_cc_bond/results/analysis/screens_5sigma/orbital_25_at_bond.png`

**Coronene v2 run created:**
- `ResearchProject/systems/coronene/run_cc_bond_v2/` — same parameters as run_cc_bond
- Coronene `shared/cpp/run_template.hpp` modified: adds `wp_wf_rt_wr` (ComplexField3DWriter) saving ψ_WP every 5×WRITE_EVERY steps
- Output: `results/raw/vti/wavefunction_wp_rt/wavefunction_t*.vti` — complex WP orbital for momentum distribution |ψ̃(k)|²

**Jellium v2 pilot validated:**
- `run_wp_n162_L50_E100_sigma1_v2`: completed, energy drift 91 μeV, WP norm 1.000, 32 WP wavefunction frames, 13 GB disk
- WP orbital extraction overhead: ~14s per call at 125³ grid — acceptable

**15 jellium v2 runs dispatched:**
- Dispatch PID running, queue-based GPU scheduling
- 9 configs created (+ E100 v2 pilot config updated with Classical struct)
- 5 WP σ=1 standard + 4 WP σ=1 HD + 3 classical standard + 3 classical HD

**Skill expanded with mandatory observables (Groups A–I):**

| Group | Name | Description | Data source |
|-------|------|-------------|-------------|
| A | Time-series diagnostics | Energy/current/dipole vs time | observables.csv |
| B | Energy decomposition | 6-panel ΔE vs Δz + bar chart | observables.csv + electron_track |
| C | Spectral analysis | FFT spectra with detrending | observables.csv |
| D | Density evolution | z-profile heatmap, δn snapshots, diff vs free propagation | density VTIs |
| E | GS & eigenvalue structure | Level diagram, DOS, GS slices, KS eigenenergy evolution | eigenvalues + state_energies |
| F | Orbital analysis | GS basis decomposition, overlap heatmaps (WP+classical+diff), occupations | overlap data |
| G | Momentum & trajectory | Centroid, momentum band (free vs jellium), σ_xyz evolution, n(k,t) | wp_momentum_stats + wp_real_space_stats |
| H | Stopping power | Knudsen, KL divergence, energy balance, force, windowed S | Various CSVs |
| I | Advanced (v2 only) | WP momentum before/after, plane-wave decomposition, loss function, SE yield | WP wavefunction VTIs |

**New observable: plane-wave decomposition of evolved KS orbitals.**
Jellium GS orbitals are plane waves with known energies E_n = ℏ²|G_n|²/2m. By projecting the time-evolved orbitals onto this plane-wave basis (FFT of each orbital at t_f), we map occupations onto plane-wave energies. This reveals:
- Which plane-wave states gained occupation (excitations)
- Which lost occupation (holes)
- The electron-hole transition structure directly
Requires all-orbital wavefunction dump at t_f (v2 runs save this).

### Files touched

- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/run_cc_bond/rewindow_leed.py` — fftshift fix + physical axes
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/cpp/run_template.hpp` — added wp_wf_rt_wr for WP wavefunction saving
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/run_cc_bond_v2/` — created (run.cpp + analyse.py + rewindow_leed.py)
- `/home/raid/skcb2/skcb2/tddft/.claude/skills/tddft-simulations/SKILL.md` — expanded mandatory observables (Groups A–I), added VTI coordinate mapping note
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/*_v2.hpp` — 9 new + 1 updated config headers
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_*_v2/` — 15 new run directories
- `/local/data/public/skcb2/tddft/ResearchProject/dispatch_v2_runs.py` — queue-based dispatch for 15 runs

### Tests and validation

- LEED fftshift: verified visually — patterns now centered with molecular structure visible
- v2 pilot: energy drift 91 μeV ✓, WP norm 1.000 ✓, max_overlap 0.00126 ✓
- WP orbital extraction feasibility: 14s overhead per call, acceptable at WF_WRITE_EVERY=10

### Known issues / blockers

- Coronene v2 queued but cannot launch until a GPU frees from jellium dispatch
- `density::total` INCLUDES WP orbital (verified) — all existing density_delta analysis was WP-contaminated. v2 runs save ψ_WP separately for exact subtraction.
- Plane-wave decomposition requires all-orbital wavefunction dump (~3.2 GB, ~40 min extraction at t_f). Only v2 runs will have this data.
- Some extra v2 run directories created by agent (σ=5 variants like run_wp_n162_L50_E100_v2, run_wp_n162_L50_E600_v2) — not in the planned 15, not dispatched

### Assumptions still in play

- WP orbital extraction cost (~14s at 125³) scales linearly with grid size — L=30 runs should be faster
- The plane-wave decomposition assumes GS orbitals are well-approximated as plane waves (exact for uniform jellium at Γ-only k-point)
- Coronene momentum distribution analysis assumes the phase of ψ_WP is needed (density alone insufficient) — verified: density loses phase information

### Exact next steps

1. **Monitor jellium v2 dispatch**: check `dispatch_v2.log` and Gmail for completion (~12-18h total)
2. **Launch coronene v2**: when a GPU frees, run `CUDA_VISIBLE_DEVICES=<free> inq-run` in `run_cc_bond_v2/`
3. **Post-process v2 runs**: run analyse.py (venv) on each completed run
4. **Implement Group B–G observables**: write Python post-processing scripts for:
   - Energy decomposition 6-panel (fig 06 pattern)
   - Energy bookkeeping bar (fig 07/C1 pattern)
   - GS basis decomposition bar (fig 11 pattern)
   - Overlap heatmaps WP + classical + diff (fig 14 pattern)
   - Density z-profile evolution (fig 09 pattern)
   - Momentum band free vs jellium (fig M_A pattern)
   - σ_xyz evolution from wp_real_space_stats
   - z-profile and density diff vs free propagation (fig D1 pattern)
5. **Implement Group I observables** (v2 data required):
   - WP momentum distribution before/after from FFT of ψ_WP VTIs
   - Plane-wave decomposition from all-orbital wavefunction dump
   - Loss function with exact WP subtraction
   - Secondary electron yield in proxy vacuum region
6. **Cross-run comparison plots**: stopping power vs energy (updated with v2 data), energy ledger bar across campaign
7. **Update observables catalogue** (`docs/observables/catalogue.md`) with new mandatory items

---

## Milestone: 2026-05-24 (late) — coronene LEED fixed, v2 pilot launched

### Current status

**In progress.** Coronene LEED 5σ re-windowing complete. Jellium v2 pilot run launched (PID 3571948, GPU 0): `run_wp_n162_L50_E100_sigma1_v2` with dt=0.01 + WP density/wavefunction saving. Awaiting pilot completion to validate the modified C++ template before scaling to 16 runs.

### What changed this sub-milestone

- **Coronene LEED fixed**: `run_cc_bond/rewindow_leed.py` recomputes screen accumulations with 5σ IFW windows. Results in `results/analysis/screens_5sigma/`. Screens 10-11 have empty windows (too close to molecule). 8 backscattering + 10 transmission screens valid.
- **v2 config created**: `shared/configs/electron_proj_E100_L50_cubic_sigma1_v2.hpp` — dt=0.01, N_STEPS=950, WF_WRITE_EVERY=10
- **v2 run.cpp created**: `run_wp_n162_L50_E100_sigma1_v2/run.cpp` — adds `density_wp` VTI + `wavefunction_wp` complex VTI writers at WF_WRITE_EVERY cadence
- **Pilot launched**: PID 3571948 on GPU 0, compiling + propagating

### Previous design decisions (unchanged)

Three workstreams ready:
1. Coronene LEED window fix (post-processing, ready to implement)
2. Jellium v2 re-runs (16 runs, C++ template modifications needed)
3. Secondary electron observable (new, full 4-quantity implementation in periodic box)

### Critical finding: density::total INCLUDES the WP orbital

Verified empirically on `run_wp_n162_L50_E100_sigma1`:
- `density_total` at WP launch position (z=−21): **0.174 e/Bohr³** (134× background n₀=0.0013)
- `density_system` and `density_total` are **byte-identical** (MD5 match)
- The `observables_reference.md` claim that `density::total` excludes the WP was **wrong**

**Consequence:** All existing `density_delta` VTI files contain the WP's own density mixed with the bath response. The spectral weight analysis in `spectral_weight.py` and `spectral_weight_full.py` was dominated by the WP signal, confirming the earlier critical analysis. The free-WP analytical subtraction partially worked but was imprecise because the interacting WP ≠ free WP.

**Solution for v2 runs:** Save ψ_WP(r,t) separately → compute |ψ_WP|² in post-process → exact bath-only δn_bath = n_total − |ψ_WP|² − n₀.

### Decisions made

**Coronene LEED fix (5σ interference-free windows):**
- Backscattering screens (z > 0): t_start = when initial WP centroid is 5σ past screen; t_end = when reflected WP centroid is 5σ before screen
- Transmission screens (z < 0): t_start = 0; t_end = when WP centroid − 5σ reaches −Lz/2
- Close-to-molecule screens may have zero-width windows — acknowledged, skip those
- Uses existing instantaneous screen snapshots from `run_cc_bond` — no re-run needed

**Jellium v2 re-runs (16 runs):**
- σ=1 standard (L=50): E={20, 25, 50, 100, 200, 300} — 6 WP runs
- σ=1 high-density (L=30): E={50, 100, 200, 300} — 4 WP runs
- Classical companions: 6 runs (matched energies, both densities)
- dt = 0.01 a.u. (was 0.02) → 2× more timesteps for better temporal resolution
- WRITE_EVERY tuned for ~300 frames
- New C++ outputs per WRITE_EVERY:
  - `density_wp` (|ψ_WP|²) as real VTI — via `density::orbital(electrons, wp_idx)`
  - `wavefunction_wp` (ψ_WP) as complex VTI — via `fields::orbital::wavefunction(electrons, wp_idx)`
- New C++ output at final IFW timestep only:
  - ALL 81 orbital wavefunctions as complex VTIs (~3.2 GB one-time dump) — for SE energy spectrum
- All existing observables preserved (observables.csv, density_total, density_delta, overlap, state_energies, momentum_distribution, wp_momentum_stats, etc.)

**Secondary electron observable (4 quantities):**
1. **Yield δ(t)**: ∫_vacuum n_SE dr, where n_SE = Σ_i|ϕ_i|² − n₀ and "vacuum" = proxy region |z − z_WP(t)| > R in periodic box
2. **Energy spectrum dN/dE**: FFT each jellium orbital in vacuum region at t_f, histogram by E = ℏ²k²/2m
3. **Emission current j_SE(t)**: flux through counting surface Σ at distance z₀ from WP trajectory
4. **Angular distribution dN/dΩ**: from momentum-space decomposition of jellium orbitals at t_f
- Implementation: full, in periodic cubic box (no slab geometry)
- Per-orbital wavefunctions saved at final IFW step for Quantities 2+4
- Quantities 1+3 computed from n_SE = n_total − |ψ_WP|² − n₀ (no per-orbital extraction needed at every step)

**Loss function L(q,ω):**
- Exact WP subtraction: δn_bath = n_total − |ψ_WP|² − n₀
- Full q_⊥-integrated pipeline from `spectral_weight_full.py`
- Expected to show clean plasmon ridge and P-H continuum with exact subtraction

### Known issues / blockers

- **Per-orbital extraction cost**: `density::orbital(electrons, wp_idx)` uses per-element GPU→host loop (~30s–minutes per call at 125³). Must verify feasibility at dt=0.01 cadence. May need to increase WRITE_EVERY for WP density/wavefunction outputs.
- **All-orbital dump at t_f**: 81 orbitals × 125³ × 16 bytes = 3.2 GB. Single extraction takes ~40 minutes. Must be done only once, at the final IFW step.
- **Secondary electron proxy in periodic box**: No true vacuum region. The proxy (|z − z_WP| > R) is an approximation. Periodic images may contaminate the "vacuum" at late times.
- **Disk space**: 16 runs × ~22 GB (WP wavefunction VTIs) = ~350 GB total. Check available disk.

### Assumptions still in play

- `density::orbital(electrons, wp_idx)` is feasible at WRITE_EVERY cadence with dt=0.01 (needs empirical test on one run)
- The periodic box proxy for "vacuum" captures meaningful secondary electron physics despite no actual surfaces
- The interference-free window for stopping power analysis is correctly computed from the self-spread cap formula
- WP momentum distribution from FFT of ψ_WP is equivalent to the `momentum_distribution.csv` WP column (should be, by construction)

### Exact next steps

1. **Fix coronene LEED windows** (post-processing): recompute screen accumulations with 5σ criteria on existing `run_cc_bond` instantaneous snapshots. Re-visualise.
2. **Modify jellium `run_template.hpp`**: add WP density + wavefunction VTI writers at WRITE_EVERY, add all-orbital dump at configurable final step.
3. **Create new config headers**: 16 configs with dt=0.01 and updated N_STEPS.
4. **Pilot one run**: test the modified template on E=100 σ=1 to verify WP extraction cost and output correctness.
5. **Launch 16 jellium v2 runs** via `/tddft-simulations` skill with GPU queue-based dispatch.
6. **Implement SE post-processing**: Python module for δ(t), dN/dE, j_SE(t), dN/dΩ. Smoke test on pilot run.
7. **Re-run loss function** with exact WP subtraction on v2 data.
8. **Update skill and catalogue** with new observables.

---

## Milestone: 2026-05-24 — 8 simulations complete, spectral weight observable implemented, skill + catalogue created

### Current status

**Complete** — all 8 new TDDFT simulations finished and post-processed. Stopping power comparison plot produced. New `tddft-simulations` skill and observables catalogue created. Spectral weight / loss function observable (Stages 1–4) implemented and tested on `run_wp_n162_L50_E20` with full q_⊥ integration. The loss function shows physically plausible structure but deconvolution is numerically fragile — bath-only density output from C++ is recommended for clean results.

### What changed

**Simulations (8 new runs, all completed):**

| Run | Type | E (eV) | Density | Wall time | S or ΔE_kin |
|-----|------|--------|---------|-----------|-------------|
| run_cc_bond | Coronene WP | 200 | n/a | 52 min | LEED screens captured |
| run_wp_n162_L50_E50_sigma1 | Jellium WP σ=1 | 50 | standard | 98 min | ΔE_kin = 4.31 eV |
| run_wp_n162_L30_E50_highdens_sigma1 | Jellium WP σ=1 | 50 | high | 15 min | ΔE_kin = 4.43 eV |
| run_wp_n162_L30_E200_highdens_sigma1 | Jellium WP σ=1 | 200 | high | 12 min | ΔE_kin = 4.65 eV |
| run_wp_n162_L30_E300_highdens_sigma1 | Jellium WP σ=1 | 300 | high | 11 min | ΔE_kin = 4.52 eV |
| run_classical_n162_L30_E50_highdens | Classical | 50 | high | 23 min | S = 1.75 eV/Bohr |
| run_classical_n162_L30_E200_highdens | Classical | 200 | high | 16 min | S = 0.67 eV/Bohr |
| run_classical_n162_L30_E300_highdens | Classical | 300 | high | 15 min | S = 0.48 eV/Bohr |

**Post-processing fixes:**
- All analyse.py re-run with venv Python (`/local/data/public/skcb2/tddft/venv/bin/python3`) — original dispatch used system Python which lacked VTK
- Classical analyse.py files replaced with proper template from `run_classical_n162_L30_E100_highdens/analyse.py` (original agent-created versions were missing `analyse_extras`, `density_fourier`, `install_schema_shims`, and several pipeline phases including `bath_energy`, `gs_projected_occupations`)
- Dispatch script patched to use venv Python for future runs

**New infrastructure:**
- `tddft-simulations` skill: `.claude/skills/tddft-simulations/SKILL.md` (504 lines, 8 phases)
- Observables catalogue: `docs/observables/catalogue.md` (all raw + post-processed observables per sim type, run inventory with coverage checklist)
- Stopping power plot: `ResearchProject/systems/jellium/stopping_power_vs_energy_all.png` + script

**Spectral weight / loss function observable:**
- `inq-stack/python/inqview/postprocess/spectral_weight.py` — on-axis (q_⊥=0) version (Stages 1–4)
- `inq-stack/python/inqview/postprocess/spectral_weight_full.py` — full q_⊥-integrated version (~2.7 GB buffer, single-pass VTI loading)
- Tested on `run_wp_n162_L50_E20` (σ=5, 344 frames, T=20.6 a.u.)
- W_resp shows bath response structure; L(q,ω) shows absorption (L<0) at low q near plasmon, sign change near Landau onset q_c
- Deconvolution coverage very low (0.01%) because free-WP analytical subtraction leaves large residual — the interacting WP differs from the free WP due to orthogonalisation, stopping, and exchange-correlation

### Files touched

| File | Action |
|------|--------|
| `ResearchProject/systems/coronene/shared/configs/cc_bond_35x35x60.hpp` | Created |
| `ResearchProject/systems/coronene/run_cc_bond/{run.cpp,analyse.py}` | Created |
| `ResearchProject/systems/jellium/shared/configs/electron_proj_E50_L50_cubic_sigma1.hpp` | Created |
| `ResearchProject/systems/jellium/shared/configs/highdens_n162_L30_E{50,200,300}_sigma1.hpp` | Created (3) |
| `ResearchProject/systems/jellium/run_wp_n162_L50_E50_sigma1/{run.cpp,analyse.py}` | Created |
| `ResearchProject/systems/jellium/run_wp_n162_L30_E{50,200,300}_highdens_sigma1/{run.cpp,analyse.py}` | Created (3) |
| `ResearchProject/systems/jellium/run_classical_n162_L30_E{50,200,300}_highdens/{run.cpp,analyse.py}` | Created (3), analyse.py replaced with proper template |
| `ResearchProject/dispatch_additional_sims.py` | Created, patched (venv Python) |
| `ResearchProject/systems/jellium/plot_stopping_power_vs_energy.py` | Created |
| `ResearchProject/systems/jellium/stopping_power_vs_energy_all.png` | Created |
| `inq-stack/python/inqview/postprocess/spectral_weight.py` | Created (on-axis version) |
| `inq-stack/python/inqview/postprocess/spectral_weight_full.py` | Created (full q_⊥ version) |
| `.claude/skills/tddft-simulations/SKILL.md` | Created (504 lines) |
| `docs/observables/catalogue.md` | Created |
| `docs/handovers/additional_simulations.md` | Updated |

### Commands run

```bash
# Dispatch 8 simulations (pair-based, 2 GPUs)
nohup python3 dispatch_additional_sims.py > dispatch_additional_sims.log 2>&1 &

# Re-run analyse.py with venv (after fixing classical templates)
for d in run_cc_bond run_wp_*_sigma1 run_wp_*_highdens_sigma1 run_classical_*_highdens; do
  cd "$d" && /local/data/public/skcb2/tddft/venv/bin/python3 analyse.py > analyse_rerun.log 2>&1
done

# Stopping power plot
/local/data/public/skcb2/tddft/venv/bin/python3 plot_stopping_power_vs_energy.py

# Spectral weight smoke tests
cd run_wp_n162_L50_E20
/local/data/public/skcb2/tddft/venv/bin/python3 -c "from inqview.postprocess.spectral_weight_full import run; run('results', ...)"
```

### Tests and validation

- **All 8 inq-run simulations**: completed (8/8 `run_completed = true`)
- **All 8 analyse.py**: produced REPORT.md (8/8 ✓). Classical runs exit rc=1 only because `knudsen_ke` phase doesn't apply — all physics phases pass.
- **GS convergence**: verified for all 3 checkpoints (coronene, jellium L=50, jellium L=30)
- **Nyquist**: all runs within limits; E=300 σ=1 HD marginal at 98%
- **Stopping power trend**: classical high-density S monotonically decreases with E (S = 1.75, 1.19, 0.67, 0.48 eV/Bohr at E = 50, 100, 200, 300) — consistent with Bethe-Bloch
- **Spectral weight**: validated against known ω_p = 3.47 eV (Bohm-Gross) — feature visible at correct frequency in W_resp

### Trusted sources used

- Tsubonoya, Hu & Watanabe, Phys. Rev. B 90, 035416 (2014) — coronene geometry
- Bohm-Gross plasmon dispersion: ω_pl(q) = ω_p √(1 + 3q²/(5k_F²))
- Lindhard P-H continuum: ω_± = |q²/2 ± q·k_F|
- Free-particle Gaussian spreading: s_t² = σ² + t²/(4σ²)

### Attribution notes

- Spectral weight pipeline (Stages 1–4) follows the methodology described in the user's instructions (derived from standard linear-response / loss-function theory, e.g. Pines & Nozières)
- analyse_extras windowed stopping power (Δz ∈ [3, 28] Bohr) from `shared/python/analyse_extras.py`

### Known issues / blockers

1. **Spectral weight deconvolution is fragile**: free-WP analytical subtraction leaves large residual because the interacting WP is modified by orthogonalisation, stopping force, and XC effects. Only 0.01% of (q,ω) points pass the V_ext threshold. **Solution**: write bath-only density from C++ (see next steps).
2. **E=600 σ=1 deferred**: Nyquist requires dx=0.30 GS (not available)
3. **Classical analyse.py rc=1**: `knudsen_ke` phase fails on classical runs (expects WP momentum data). Non-critical — should be removed from classical PIPELINE_PHASES in template.
4. **Some early σ=1 runs missing REPORT.md**: `run_wp_n162_L50_E20_sigma1` and `E25_sigma1` have analyse.py but were never run with the full template

### Assumptions still in play

- Free-WP spreading formula is an upper bound for the actual WP spread in the bath (conservative for N_STEPS cap)
- The analytical free-WP subtraction in the spectral weight pipeline assumes the WP's ballistic trajectory is close to free propagation — valid for short times but degrades as stopping effects accumulate
- The q_⊥-integrated spectral weight treats all q_⊥ modes equally — no weighting by the WP's momentum-space envelope

### Exact next steps

1. **For clean loss function L(q,ω)**: modify the jellium `run_template.hpp` to write a **bath-only density** VTI series (`density_bath = density_total - density_orbital(wp_idx)`). This eliminates the need for analytical WP subtraction and gives a clean δn_resp directly from the simulation. The `density::orbital(electrons, wp_idx)` call is expensive (per-element GPU→host loop) but only needs to be called every WRITE_EVERY — similar cadence to existing density_delta. Alternatively, write `density_total - density_wp_initial_shifted` using the stored WP from the injection step.

2. **Re-run spectral weight on longer simulations**: the E=20 σ=5 run (T=20.6 a.u.) gives dω = 0.52 eV resolution. Runs with σ=5 at lower energies (E=1.5, plasmon runs) may have even longer T and better spectral resolution.

3. **Complete missing analyse.py runs**: run `analyse.py` on `run_wp_n162_L50_E20_sigma1` and `E25_sigma1` with proper template and venv.

4. **Integrate results into report**: the stopping power plot, coronene C-C bond LEED, and spectral weight maps are ready for Report 1 figures.

5. **Refine tddft-simulations skill**: add `knudsen_ke` to the "skip for classical" list; add spectral_weight as an optional Tier 3 observable.

---

## Milestone: 2026-05-24 (early) — dispatch launched, 8 simulations building

### Current status

Superseded by the milestone above. Original dispatch launched all 8 runs successfully (inq-run 8/8). Initial analyse.py failed due to system Python (no VTK) — fixed by re-running with venv.

### Config headers created (5)
- `cc_bond_35x35x60.hpp` (coronene)
- `electron_proj_E50_L50_cubic_sigma1.hpp` (jellium)
- `highdens_n162_L30_E{50,200,300}_sigma1.hpp` (3 files)

### Validation performed
- GS convergence: 3 checkpoints verified
- Nyquist: all runs within limits
- Self-spread caps: computed for all σ=1 runs
- Memory: all runs well within 24 GB A30 limit
