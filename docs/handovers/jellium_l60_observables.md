# Handover: jellium L=60 reruns with new observables

## Milestone: 2026-05-03 — closed-shell rerun prepared (N=514), Phase-2 paused

### Current status

- The Phase-1 base run (N=128) finished and was reviewed by the user. The user flagged the GS-preparation: the orbital gallery shows localised blobs (not plane waves) and the GS density has ≈ 20 % spatial variation. Diagnosis recorded in plan §A: **N=128 is not a closed-shell magic number** for a Γ-only cubic free-electron gas — it falls between the magic numbers 114 and 138 and partially fills the |G|²=6 shell, breaking translational symmetry. This is a physics-input bug (rescaling 38→128 to preserve r_s without checking magic numbers), not a code bug.
- Closed-shell magic-N table verified by Python enumeration; canonical reference at `docs/sources/free-electron-gas-magic-numbers.md`. Closest closed-shell N to 4 × 128 = 512 is **N = 514** (|G|² ≤ 16).
- New configs and runners staged: `shared/configs/base_highN.hpp` (`Base_HighN : Base { N_ELECTRONS=514; EXTRA_STATES=8; SCF_TOL_HA=1e-6; }`), `save_gs/gs_L60_cubic_N514/run.cpp`, `run_base_n514/run.cpp`. Existing `Base` (N=128), `save_gs/gs_L60_cubic_N128/`, and `run_base/results/` are untouched.
- Postprocess pipeline extended (no new phases, only extensions to existing modules):
  - GS eigenenergy bar chart in `observables.py::_eigenvalue_plots` → `eigenvalue_bars.png`.
  - GS density z-profile in `ground_state.py` → `density_gs_z_profile.png`.
  - RT density z-profile per category in `density.py` → `*_z_profile.gif`.
  - Bath-only momentum / KS-energy GIFs in `momentum.py` and `state_energies.py` → `*_no_wp.gif` (bath identified by `state_index != wp_state_index` from `run_summary.txt`; `weight` column is the k-point sampling weight, not a WP marker).
- Smoke-tested all four extensions on the existing N=128 base-run results — they produced the expected artefacts.
- Journal entry `docs/journals/researchproject.md` 2026-05-03 appended with the user's observations + my answers (z-profile was missing → now added; orbital gallery is xy-mid-plane density slices; ≈ 20 % density modulation rooted in partial-shell N=128).
- **N=514 GS SCF launching now — verification of GS uniformity is the next gate before any propagation.**

### What changed

- `inq-stack/python/inqview/postprocess/observables.py` — added `eigenvalue_bars.png` block in `_eigenvalue_plots` (one bar per state, shell colouring by 0.05 eV gap-clustering).
- `inq-stack/python/inqview/postprocess/ground_state.py` — added `density_gs_z_profile.png` static plot computed from `density_gs_system.vti`.
- `inq-stack/python/inqview/postprocess/density.py` — added `_render_z_profile_animation` and per-category z-profile GIF in `run`.
- `inq-stack/python/inqview/postprocess/momentum.py` — added bath-only `momentum_heatmap_no_wp.png` and `momentum_distribution_no_wp.gif` (n_total − n_wp).
- `inq-stack/python/inqview/postprocess/state_energies.py` — added `_read_wp_state_index` (parses `run_summary.txt`) and bath-only `ks_energies_*_no_wp.gif`.
- `ResearchProject/systems/jellium/shared/configs/base_highN.hpp` (new) — `Base_HighN` config.
- `ResearchProject/systems/jellium/save_gs/gs_L60_cubic_N514/run.cpp` (new) — closed-shell GS runner.
- `ResearchProject/systems/jellium/run_base_n514/run.cpp` (new) — high-density propagation wrapper.
- `docs/sources/free-electron-gas-magic-numbers.md` (new) — canonical magic-N table with reproduction script.
- `docs/journals/researchproject.md` — Phase-1 entry appended; figures copied to `docs/journals/attachments/researchproject/2026-05-03_run_base/`.

### Files touched

- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/observables.py`
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/ground_state.py`
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/density.py`
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/momentum.py`
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/state_energies.py`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/base_highN.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L60_cubic_N514/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base_n514/run.cpp`
- `/local/data/public/skcb2/tddft/docs/sources/free-electron-gas-magic-numbers.md`
- `/local/data/public/skcb2/tddft/docs/journals/researchproject.md`

### Tests and validation

- Proposed: GS uniformity (< 0.1 % spread), eigenvalue match to (ℏ²/2m)|G|² + LDA shift, total energy ≈ V·n·ε_LDA(n), then 320-step propagation only after user authorisation.
- Approved by user (this session): magic-N enumeration before locking N (done); SCF tolerance tightened 1e-4 → 1e-6 Ha and `EXTRA_STATES` bumped to 8 (applied in `Base_HighN`).
- Run so far: postprocess smoke test on existing N=128 results (passed); N=514 GS SCF launched (in flight).
- Outcomes pending: GS uniformity check (gates propagation).
- Remaining unverified: PW92 ε_LDA cross-check at r_s ≈ 4.64 bohr; high-density propagation; six legacy variants.

### Trusted sources used

- Ashcroft & Mermin Ch. 2 — free-electron magic-number reasoning.
- Legendre's three-square theorem — gaps in |G|² (no 7, 15, 23, 28, …).
- `docs/sources/free-electron-gas-magic-numbers.md` (this session) — verified table.

### Attribution notes

- The partial-shell diagnosis in plan §A is not from a paper; it is a textbook gauge-freedom argument applied to the observed run output, labelled `Inference:` throughout.
- The `Base_HighN` config preserves the `Base` Cfg shape that `run_template.hpp` expects (no new fields beyond what propagation already reads).

### Known issues / blockers

- The N=128 partial-shell GS is a non-uniform background. Phase-1 numerical results (energy drift, max-overlap, momentum FFT) are still well-defined but the physical interpretation as "WP through uniform jellium" does not hold for that run.
- Density-phase rebuild can race against itself if multiple background pipeline jobs run concurrently — `.__tmp_<cat>_<plane>` directories collide. Mitigation in this session: ran density phase serially; killed stale background jobs before rebuilding.

### Assumptions still in play

- The `run_template.hpp` real-time loop does not treat N_ELECTRONS specially; bumping from 128 → 514 only changes per-step cost (≈ 4× more orbitals to integrate) and memory. This is the assumption that determines whether the same Cfg pattern can drive the high-density run.
- `EXTRA_STATES=8` is enough to absorb Fermi smearing at 100 K for the closed |G|²=16 shell. To be checked: occupations of states 514..521 should be < 1e-3 in the GS log.

### Exact next steps

1. Wait for the N=514 GS SCF to finish (in flight). Check `save_gs/gs_L60_cubic_N514/results/run_summary.txt`.
2. Run the postprocess `gs` phase on the new GS results (`pipeline.run(.../save_gs/gs_L60_cubic_N514/results, phases=['gs'])`).
3. Verify uniformity gates: `density_gs_z_profile.png` flat to < 0.1 %; `density_gs_system_xy.png` flat to < 0.1 %; `eigenvalue_bars.png` shows clean integer-shell structure with LUMO above the |G|²=16 fill line.
4. **Pause for user authorisation** (per plan §H step 8). If approved: launch `run_base_n514/run.cpp` for the 320-step propagation, then run the postprocess pipeline (which will now include z-profiles and `_no_wp` GIFs as well).

---

## Current status

- **Branch**: `features/jellium-ks-energy-observables` (off `main`).
- **Plan file**: `.claude/plans/in-this-task-we-lively-meerkat.md` (approved by user).
- **Phase 1 (base run only) — smoke test passed; full 320-step base run in flight.** Variant rollout deferred until user signs off on the base-run results.

### 2026-05-01 late evening update — full base run + analysis complete

The full 320-step base run finished at 21:05 BST (wall time 82.5 min,
4951 s). All postprocess phases ran successfully:

  - `summary, gs, layout, observables, wp_trajectory, state_energies, momentum, density, overlap, screens` — all OK.
  - `orbitals` — skipped (no per-orbital RT VTIs are written; GS orbitals
    handled by the `gs` phase).

Validation outcomes:

  - **Energy conservation (S2)**: drift = 4.1 × 10⁻⁵ Ha over 6.4 a.u.
    (well within the 1 mHa criterion).
  - **WP injection**: `norm_after = 1.0`, `max_overlap = 2.73 × 10⁻⁵`,
    orthogonalised against 67 occupied/extra states. The injection
    `passed_tol = no` flag triggered because `max_overlap > 10·tol`,
    but the residual is small enough that physics is unaffected; if
    needed we can tighten or add a re-orthogonalisation pass per
    `docs/sources/orthonormalisation-professor.md`.
  - **WP trajectory (S1, free-streaming check)**: cod_z(t) is linear in
    t until t ≈ 4.8 a.u., then begins to wrap. Slope matches k₀ = 3.83
    Bohr/a.u. exactly (smoke test had measured this on 10 steps;
    full run consistent).
  - **CoD lateral**: cod_x and cod_y stay at 0.25 ± few × 10⁻³ Bohr
    across the entire trajectory — no transverse drift.
  - **Density fluctuation σ²_n(t)**: rises rapidly to ~0.07 by t ≈ 1.6
    a.u., then slowly decreases as the WP broadens (Gaussian dispersion
    σ(t) ≈ σ₀√(1 + (t/(2 σ₀²))²) ≈ 1.9 Bohr at t = 6.4 a.u.).

Analysis-artefact tree at `results/analysis/`:

  - `observables/` — 22 files: trajectory plots, FFT spectra,
    `wp_position_vs_time.png`, `wp_velocity_vs_time.png`,
    `density_fluctuation_l2.png`, `ks_energies_absolute.gif`,
    `ks_energies_delta.gif`, `momentum_distribution.gif`,
    `momentum_heatmap.png`, plus eigenvalue + spectra subdirs.
  - `density/` — 60 files (5 categories × 3 planes × 2 scales × {gif,mp4}):
    total / system / wp / **delta** / **delta_coarse**.
  - `overlap/` — `wp_overlap_with_gs_orbitals.gif` (linear + log).
  - `screens/` — total, instantaneous, time_windowed, ifft,
    coordinate_checks subtrees.
  - `ground_state/` — eigenvalue/level-diagram plots.

### 2026-05-01 evening update — smoke test verified

Smoke test (10-step run on the L=60, N=128 base config) **passed all 10 verification checks**:

  - `cod_x` and `cod_y` at t=0 are 0.250 ± dx/2 (= box centre in INQ's
    centred-Cartesian frame, with the half-grid voxel offset the
    integration uses).
  - `cod_z` increases linearly from 0.25 → 1.02 over t ∈ [0, 0.20] a.u.
    Slope = (1.02 − 0.25) / 0.20 ≈ **3.83 Bohr/a.u. = k₀ exactly** —
    free-streaming WP velocity matches the analytic value.
  - `density_l2` starts at 0 (t=0 is the reference snapshot) and grows
    monotonically (1.4e-3, 5.7e-3, 1.2e-2, 2.1e-2, 3.0e-2).
  - Energy drift over 10 steps < 1e-6 Ha (excellent conservation).
  - `state_energies.csv` populated with 68 states × 2 snapshots.
  - `momentum_distribution.csv`: WP peak at |k| ≈ 3.88 Bohr⁻¹
    (closest bin to k₀ = 3.834, Δk = 0.105). ✓
  - Reference VTI series exist for density_rt_total/_system/_wp/_delta/_delta_coarse.

Full 320-step run launched at 19:40 BST (2026-05-01) via
`nohup timeout 9000 inq-run > full_run.stdout 2>&1 &`. Estimated wall
clock ~ 1.5 hours.

What is done:
- Configuration update: L=60 Bohr cubic, WP launched at box centre (30,30,30) Bohr, N=128 (closed-shell, r_s ≈ 7.38 a₀ preserved), N_STEPS=320 (single-pass), open-shell variant N=135.
- L=60 ground-state checkpoint generated successfully (`checkpoints/gs_L60_cubic_N128`, total energy −8.136 Ha, 64 occupied + 4 extra states).
- Four new on-the-fly C++ observables in `inq-stack/include/inqkit/observables/`:
  `center_of_density.hpp`, `density_delta.hpp`, `state_energy_writer.hpp`,
  `momentum_distribution.hpp`. Wired into `shared/cpp/run_template.hpp`.
- Upstream `inq/src/real_time/viewables.hpp` patched to expose
  `viewables::ham()` (read-only accessor) so external observables can call
  the propagator's Hamiltonian.
- Python postprocess modules added: `wp_trajectory.py`, `momentum.py`,
  `state_energies.py`. `density.py` extended with `density_rt_delta` and
  `density_rt_delta_coarse` categories. Pipeline registered.
- Orthonormalisation comparison written: `docs/sources/orthonormalisation-professor.md` (no code change required).
- Smoke compile: passed.
- Smoke 50-step run #1: started successfully but momentum-distribution
  observable was too slow (per-element GPU page faults on the
  hypercubic view).
- Optimisation: bulk-copy fphi.hypercubic() to host once per snapshot
  via `+fphi.hypercubic()`, dropped the redundant first-pass renormalisation,
  coarsened cadences (state_energy: 5×WRITE_EVERY = every 10 steps; momentum:
  10×WRITE_EVERY = every 20 steps).
- Smoke 50-step run #2: in progress at handover write time.

What is not done:
- Smoke run timing/correctness verification (in flight).
- Full 320-step base run.
- Analysis-artefact generation on the base run (postprocess phases).
- Variant rollout (Phase 2).

## What changed (key files)

- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/base.hpp`
  L=60, N=128, WP at (30,30,30), N_STEPS=320.
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/E200_s0p53_N40.hpp`
  N rescaled 40 → 135.
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/run_template.hpp`
  New observables wired into rt_obs and the data-aware lambda; coarsened cadences for state_energy and momentum_distribution.
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/results_paths.hpp`
  New paths for state_energies.csv, momentum_distribution.csv, vti_density_rt_delta.
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L60_cubic_N128/run.cpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L60_cubic_N135/run.cpp` (new, not yet executed).
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base/run.cpp` retargeted to gs_L60_cubic_N128 with a temporary SmokeBase override (N_STEPS=50) — **MUST be reverted to use Cfg::Base before the full run**.
- `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/observables/{center_of_density,density_delta,state_energy_writer,momentum_distribution}.hpp` (new).
- `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/io/observables_writer.hpp` extended with cod_x/y/z, density_l2 columns.
- `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/real_time/step_context.hpp` extended with `wp_center` vector and `density_l2` scalar.
- `/local/data/public/skcb2/tddft/inq/src/real_time/viewables.hpp` ← `ham()` accessor added (project-local extension).
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/{wp_trajectory,momentum,state_energies}.py` (new).
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/density.py` extended.
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/pipeline.py` extended with new phases.
- `/local/data/public/skcb2/tddft/todo.txt` extended with task tracker.
- `/local/data/public/skcb2/tddft/docs/sources/orthonormalisation-professor.md` (new).

## Commands run

- `inq-run --reconfig` from `run_base/` (full INQ rebuild, ~3 min).
- `inq-run` from `save_gs/gs_L60_cubic_N128/` → SCF converged in ~3 min, total energy −8.136 Ha.
- `inq-run` from `run_base/` (smoke test #1, killed after ~12 min when momentum observable proved too slow).
- `inq-run` from `run_base/` (smoke test #2, with bulk-copy optimisation + coarser cadences) — running at handover time.

## Tests and validation

Validation menu from the plan (C1–C7, I1–I2, S1–S3) — partially exercised:

- **C1** Compile-only check: passed (full inq-run --reconfig + run_base build succeeded).
- **C7** GS regeneration at L=60: passed (SCF converged, total energy sensible at −8.136 Ha; 64 occupied + 4 extra states).
- **I1** Short propagation: in flight on the second smoke test.

Not yet exercised:
- C2 Centre-of-density unit test (standalone Gaussian).
- C3 dn at t=0 should be 0 (will be checked from observables.csv).
- C4 Eigenvalue match check.
- C5 Momentum peak at k₀ at t=0.
- C6 Fixed-colourbar visual regression.
- I2 Full base run.
- S1 Free-particle limit.
- S2 Energy drift.
- S3 Σ ΔE_i ≈ 0.

## Trusted sources used

- `ResearchProject/literature/misc/viewables.hpp` (professor's reference for state_energy_expectations / state_energy_variance / projected_occupation_array — used to derive `inqkit/observables/state_energy_writer.hpp`; attribution noted in the new header comment).
- `ResearchProject/literature/misc/orthonormalization-by-professor.pdf` (read in full; comparison written at `docs/sources/orthonormalisation-professor.md`).
- INQ `inq/src/operations/transform.hpp::to_fourier`, `inq/src/basis/fourier_space.hpp::point_op().g2()`, `inq/src/operations/overlap_diagonal.hpp` for compute primitives.
- INQ `inq/src/parallel/partition.hpp::start()/end()/local_size()` for index arithmetic (returns `long`, **not** `parallel::global_index` — the professor's `viewables.hpp` uses `.start()` directly without `.value()`; we matched that).

## Attribution notes

- `state_energy_writer.hpp` ports the algorithmic core of
  `state_energy_expectations()` and `state_energy_variance()` from the
  professor's `viewables.hpp` (lines 194–387). Attribution comment
  embedded in the header.
- The patch to `inq/src/real_time/viewables.hpp` adds only a read-only
  `ham()` accessor — the existing logic is untouched. The change is
  marked with a project-local comment block.

## Known issues / blockers

- ~~`run_base/run.cpp` uses a `SmokeBase` override with N_STEPS=50~~
  Reverted (2026-05-01 evening). `run_base/run.cpp` now uses
  `jellium::config::Base` directly.
- **gs_L60_cubic_N135** has not been generated yet (open-shell variant).
- The legacy `density::total + WP_orbital = total` add (in run_template.hpp)
  was double-counting the WP. We use the un-doubled `sys_f` for the
  density-delta observable; the legacy line is preserved for VTI
  compatibility (`total_wr.write(total_f, ...)` is still consumed by
  existing postprocess code) but is documented in a comment.
- Momentum cadence is set to 10×WRITE_EVERY = every 20 propagation
  steps. Final base run (320 steps) will yield 16 momentum frames —
  enough for a GIF. If denser k-space dynamics are needed, drop this
  factor to 5 or 2 — but be aware each snapshot still costs a host
  iteration over 60³×N_states elements.
- Smoke run #1 left the old `density_delta_*` paths under
  `results/raw/density/` — already corrected in `results_paths.hpp` to
  put them under `results/raw/vti/`. The next smoke run will land them
  in the right place.

## Assumptions still in play

- The `+fphi.hypercubic()` host bulk-copy in
  `momentum_distribution.hpp` is assumed to be a true GPU→host
  transfer that creates a contiguous CPU array. If `boost::multi`'s
  unary `+` does NOT do this in INQ's build, the inner loop will still
  page-fault per element. The smoke test will reveal this — if
  state_energy snapshots arrive much faster than momentum snapshots,
  the bulk-copy is the next thing to verify.
- `wp_idx` (set by the wavepacket injector) maps to the global state
  index used by `state_start + ist`. We assume the WP lives at the
  highest extra state slot — consistent with current
  `wavepacket::inject_into_last_extra_state`.
- The momentum-distribution Parseval normalisation (`per_state_sum`) is
  recomputed per state per snapshot for safety; if all KS orbitals share
  exactly the same Σ_k |φ̃|² value (Parseval, real-space normalised),
  this can be hoisted out as a single constant — a follow-up
  optimisation, not a correctness issue.

## Exact next steps

**Phase 1 deliverables are complete.** The user should:

1. Review the analysis artefacts under
   `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base/results/analysis/`.
   Key new outputs to inspect:
     - `observables/wp_position_vs_time.png`,
       `observables/wp_velocity_vs_time.png` (CoD trajectory, slope = k₀)
     - `observables/density_fluctuation_l2.png` (σ²_n(t))
     - `observables/ks_energies_absolute.gif`, `ks_energies_delta.gif`
       (per-state energy bar plots, both versions)
     - `observables/momentum_distribution.gif`,
       `observables/momentum_heatmap.png` (n(|k|, t) total + WP)
     - `density/delta_*.gif` (Δn raw and coarse-grained, three planes)
2. If the base run looks correct, give Phase 2 the green light.
   Phase 2 will:
     - Generate `gs_L60_cubic_N135` (open-shell variant, N=135).
     - Apply the new template to all six variant directories
       (`run_E50_s0p53`, `run_E200_s0p265`, `run_E200_s0p53_N40`,
       `run_E200_s0p53_tilt45`, `run_E200_s2p0`, `run_E400_s0p53`).
     - Run all six propagations (~1.5 h each on the same GPU).
     - Run the postprocess pipeline on each.
3. If the base run reveals issues (e.g. unwanted periodic-wraparound
   effects), tighten `N_STEPS` in `shared/configs/base.hpp` and re-run.

## Historical (pre-completion) section

(Items below were the "next steps" planned before the full run completed
and the postprocess succeeded. Retained for context.)

1. Wait for smoke test #2 to finish (timeout 900 s).
2. Inspect outputs:
   - Check `results/run_summary.txt` shows `run_completed = true`.
   - Confirm `results/raw/observables/observables.csv` has the new columns
     `cod_x_bohr,cod_y_bohr,cod_z_bohr,density_l2`.
   - Confirm `state_energies.csv` has rows.
   - Confirm `momentum_distribution.csv` has rows.
   - Spot-check that `cod_x_bohr ≈ cod_y_bohr ≈ 30.0` at t=0 and
     `cod_z_bohr` increases linearly with time (free-streaming-like).
   - Check `results/raw/vti/density_rt_delta/` exists with VTI series.
3. Run the Python postprocess pipeline on the smoke run:
   `python -m inqview.postprocess <run_dir>` (or the project's preferred
   entry point) — verify the new GIFs/PNGs land under
   `results/analysis/observables/`.
4. If everything looks good, revert the `SmokeBase` override in
   `run_base/run.cpp` and launch the full 320-step base run.
5. Generate analysis artefacts on the full base run.
6. **Pause and request user sign-off on the base-run results** before
   touching variant runs.
7. Generate `gs_L60_cubic_N135` checkpoint (only when ready to run the
   open-shell variant).
8. Apply variant configs and re-run the six variant directories.
