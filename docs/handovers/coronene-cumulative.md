# Coronene replication — cumulative handover

Read-once narrative summary of the entire coronene WP-RT-LEED replication
effort. The rolling per-session journal lives in
`docs/handovers/coronene-replication.md`; this file is the canonical
"what was done, why, and what's next" digest.

---

## What this work delivered

**A reproducible coronene WP-RT-LEED framework** under
`ResearchProject/systems/coronene/`, replacing every legacy buggy run.
Specifically:

1. **One canonical geometry** — `shared/geometry/coronene.xyz`, atoms at
   z = 0, reusable across any orthorhombic cell whose half-extent
   contains the molecule plus the WP launch height.
2. **Three GS save runs** — `save_gs/{gs_35x35x60_cut40,
   gs_35x35x80_cut40, gs_35x35x40_cut40}/` — one checkpoint per unique
   `(cell, cutoff)` tuple, reused by every propagation that matches.
3. **Ten propagation runs** — `run_base, run_E30, run_E800, run_s0p33,
   run_s3, run_E800_s0p33, run_E30_s3, run_b18_35x35x80, run_b6_35x35x80,
   run_35x35x40` — each a thin Cfg-templated wrapper around
   `shared/cpp/run_template.hpp`. N_steps per run is **physics-derived**
   from the WP arrival time at the far box face (see
   `compute_n_steps()` in `shared/configs/tsubonoya_2014_base.hpp`); per-screen
   accumulator windows are derived from the WP centroid's arrival /
   departure time at each screen (see `compute_screen_window()` in
   `shared/cpp/leed_screen_layout.hpp`).
4. **Multi-GPU dispatcher** — `scripts/dispatch_runs.py` polls
   `nvidia-smi` for free GPUs, holds at most one of *our* jobs per GPU,
   refills slots as runs finish, leaves other users' jobs alone.
5. **Generalisable Python postprocess pipeline** —
   `inq-stack/python/inqview/postprocess/{pipeline, run_summary,
   ground_state, layout, observables, density, screens, overlap,
   orbitals, compare, paraview_3d}.py` — phase-based, with `--rebuild`,
   `--phases <subset>`, log + linear pairs, MP4 + GIF dual format
   (where ffmpeg is available), and a layout xz diagram per run.
6. **Coronene-specific CLI** —
   `ResearchProject/systems/coronene/scripts/coronene_postprocess.py` —
   thin wrapper exposing `run` and `hypothesis` subcommands.
7. **Six hypothesis comparison folders** under `coronene/hypotheses/`,
   each populated with overlay PNGs from `compare.run_hypothesis()`,
   plus a `physics/` subfolder with current/dipole overlays + their
   FFT spectra + the WP residual-at-t_final bar chart.
8. **ParaView 3D volume-render videos** per run (head-on + 3/4 view),
   overlaying `density_rt_system + density_rt_wp` with log-scale colour
   and density-tied opacity.

---

## Key bugs found and fixed (chronological)

1. **`std::sqrt` not constexpr under CUDA** — replaced with a Newton-iteration
   `const_sqrt(x)` in `shared/configs/tsubonoya_2014_base.hpp` so every Cfg
   variant can compute its own `WP_K0` and `N_STEPS` at compile time.
2. **`*_path()` helpers didn't ensure parent dir** — `wavepacket_config.txt`,
   `injection_report.txt`, `screen_config.csv`, `window_ranges.csv` were
   silently dropped on the first three runs. Fixed by
   `ensure_parent(path)` in `shared/cpp/results_paths.hpp`. The three
   already-launched runs (`run_base, run_E30, run_E800`) miss those four
   metadata files; every datum is duplicated in `run_summary.txt` so the
   loss is recoverable without a re-run.
3. **Postprocess `density` phase glob mismatch** — looked for
   `{cat}_t*.vti` but the C++ writer uses the layout's `field_name`
   (`density_t*.vti`). Fixed in `_common.list_vti_series` to glob
   `*_t*.vti`.
4. **Postprocess `gs` phase no-op when only VTI present** — added a
   VTK-based VTI loader in `inqview/postprocess/ground_state.py`; the
   GS slice + orbital gallery now work whether the run wrote raw+meta
   or VTI-only.
5. **LEED `.dat` four-corner-split** — the C++ writer emits FFT-natural
   index order; the Python loader now applies `np.fft.fftshift` and
   overrides the origin to span `[-L/2, +L/2]`. See
   `docs/notes/coronene-geometry-correction.md` for the pedagogical
   reasoning.

---

## Phase-1 timeline (executed)

| Stage | Runs | Walltime | GPU |
|---|---|---|---|
| Smoke build of the whole INQ tree | (any single run) | ~7 min CMake + nvcc | first run only |
| GS save × 3 | gs_35x35x60_cut40 + gs_35x35x80_cut40 + gs_35x35x40_cut40 | 9 + 9 + 7.5 min | parallel across GPUs as available |
| Smoke `run_base` | 560 steps | 47.5 min | 1 |
| Remaining 9 propagations | 1446 + 280 + 551 + 587 + 276 + 1514 + 769 + 613 + 430 steps | 4 h 30 min wall (two GPUs) | 0 + 1 |
| Postprocess (Phase 1 — original) | 10 runs × ~5 min each | ~50 min CPU | n/a |
| Hypothesis comparisons (Phase 1) | 6 folders × ~30 s | ~3 min | n/a |

Total Phase-1 wall time: **~6 h 15 min** end to end.

## Phase-2 work (complete)

| Item | Status |
|---|---|
| LEED fftshift fix (Python loader) | ✅ done |
| Multi-line animation titles + dual-format (gif + mp4) | ✅ done |
| Linear + log heatmap pairs (density GIFs, WP-overlap GIF) | ✅ done |
| Layout xz diagram per run (new `layout` phase) | ✅ done |
| Total energy / all-energies plot split | ✅ done |
| WP-overlap y-axis range fix | ✅ done |
| Hypothesis `physics/` subfolder | ✅ done |
| ParaView 3D `paraview_3d` phase (overlay, two cameras) | ✅ done |
| Re-postprocess all 10 runs + 6 hypotheses with --rebuild | ✅ done |
| Three pedagogical docs (this file + 2 sibling) | ✅ done |
| Phase-2 commit + merge to `main` | ✅ done |

---

## Phase-3 work (FFT-ordering fixes + extended postprocess + eigenvalues)

Triggered by user observation: "every transmission screen shows the
static coronene electron cloud, not a diffraction pattern". Audit
identified two C++ FFT-ordering bugs and one window-logic error.

### Bugs fixed

1. **`PlaneScreen::iz_nearest` clamped negative z to 0**
   (`inq-stack/include/inqkit/screens/plane_screen.hpp`). Replaced clamp
   with FFT-natural wrap `((iz % Nz) + Nz) % Nz`. Every transmission
   screen at physical z < 0 now lands in the correct upper half of the
   FFT-natural array. **Fix:** commit `c8343c5`.
2. **`fields::orbital::wavefunction` did not apply `fft_shift_index`**
   (`inq-stack/include/inqkit/fields/orbital.hpp`). Mirrored the
   identical pattern from `density.hpp:88-99`. Complex-orbital VTIs and
   downstream slices are now spatially correct. **Fix:** same commit.
3. **`compute_screen_window` used "during transit" window**
   (`ResearchProject/systems/coronene/shared/cpp/leed_screen_layout.hpp`).
   Re-derived: forward (z<0) starts at `max(0, (b+σ−z)/|k|)` and runs to
   `N_STEPS·dt`; backscatter (z≥0) starts at the same expression and
   ends at `(b+L_z/2−σ)/|k|`. The σ scaling is configurable via the new
   `WP_ENVELOPE_SIGMAS` config parameter (default 2.0). **Fix:** commit
   `d827db8`.

### New artefacts shipped in Phase 3

| Item | Location | Status |
|---|---|---|
| IFFT helper (Patterson + amp_only) | `inq-stack/python/inqview/postprocess/_ifft.py` + `screens.LeedPattern.inverse_fft` | ✅ done |
| IFFT outputs in `analysis/screens/ifft/` subfolder (linear + log) | per-run | ✅ done |
| `WP_ENVELOPE_SIGMAS` Cfg parameter (=2.0 sigmas) | `shared/configs/tsubonoya_2014_base.hpp` | ✅ done |
| Dispatcher `--clear-results` + pre-clear before each rerun | `scripts/dispatch_runs.py` | ✅ done |
| Extended preprocessed spectra: raw / mean / detrended × {dipole_z, current_z, energy_total} | `inqview/postprocess/observables.py` | ✅ done |
| Compartmentalised `spectra/{current,dipole,energy}/` per run | `inqview/postprocess/observables.py` | ✅ done |
| Zero-padded smoother spectra (`pad_factor=4`) | `inqview/postprocess/observables.py` | ✅ done |
| Jellium spectrum rollout (uses inqview helpers, flat layout) | `ResearchProject/jellium/jellium-wp-rt/jellium_spectra.py` | ✅ done |
| GS eigenvalues + occupations C++ writer | `shared/cpp/eigenvalues_writer.hpp` (called from `save_gs/*/run.cpp`, `run_template.hpp`) | ✅ done |
| Eigenvalue retrofit script | `scripts/retrofit_eigenvalues.py` | ✅ done |
| Eigenvalue postprocess viz (`eigenvalues_levels.png`, `eigenvalues_dos.png`, `eigenvalue_table.txt`) | `inqview/postprocess/observables.py::_eigenvalue_plots` | ✅ done |
| Refreshed `docs/observables_reference.md` | `docs/` | ✅ done |
| `todo_later.md` entries (FFT bugs, jellium relocation, drift-method question, eigenvalues retrofit, 3 future-useful follow-ups) | `docs/` | ✅ done |

### Branch strategy executed

| Branch | Tip commit | Merge into main |
|---|---|---|
| `fixes/fft-ordering-arrays` | `2aadd71` | merged as `f6065f8` |
| `coronene-fft-fixed-base` (run_base re-run + verify) | (transient) | merged |
| `coronene-fft-fixed-rerun` (9-run rerun + ifft/envelope/spectra/eigenvalues) | `31cfd56` | **pending merge** |

### Branch-3 results (all 9 propagation runs done, exit 0)

| Run | Walltime |
|---|---|
| run_E30 | 7509 s |
| run_E800 | 1862 s |
| run_s0p33 | 3453 s |
| run_s3 | 2943 s |
| run_E800_s0p33 | 2013 s |
| run_E30_s3 | 7118 s |
| run_b18_35x35x80 | 6758 s |
| run_b6_35x35x80 | 4895 s |
| run_35x35x40 | 1772 s |

Each run's auto-postprocess (`scripts/auto_postprocess.sh` + watcher)
produced the full `analysis/` tree minus `paraview_3d/` (skipped by
default; needs `--with-paraview`). Run_base was re-run separately on
Branch 2 with --with-paraview.

### Phase-3 closeout (in progress, 28 Apr)

| Item | Status |
|---|---|
| Re-run save_gs/* to populate eigenvalues.csv in checkpoints | 🔄 in progress (gs_35x35x60 + gs_35x35x80 launched 11:06; gs_35x35x40 queued) |
| Run `scripts/retrofit_eigenvalues.py` over all 10 propagation runs | 📋 pending (blocked on save_gs) |
| Re-postprocess `gs` phase per run (eigenvalue viz) | 📋 pending (blocked on retrofit) |
| 6 hypothesis comparisons | 🔄 in progress (started 11:06; 00_base done 4 s; rest expected ≤ 5 min) |
| ParaView 3D videos for the 9 Branch-3 runs (`--with-paraview`) | 📋 pending (run_base already has them) |
| Commit & merge `coronene-fft-fixed-rerun` into main | 📋 pending |
| Update this handover (Phase-3 section + closeout) | ✅ done (this section) |

---

## Future aims (everything the user has flagged for later)

These items are not blockers for the current framework but should be
addressed before paper-quality results are claimed:

1. **Re-run `03_ecut_convergence` with a non-dojo PSP family**
   (e.g. ONCV-PBE from Schlipf & Gygi). The current sweep on dojo PSPs
   shows non-monotonic E_total above ~40 Ha; non-monotonicity is
   pseudopotential-specific and should be confirmed/refuted on a clean
   PSP set. Documented in `docs/todo_later.md`.
2. **Geometry relaxation of coronene** before scattering. The XYZ used
   has bond lengths from a crystallographic template; max force is
   0.063 Ha/Bohr, well above 0.001 Ha/Bohr convergence. A relaxed
   geometry will shift the molecular charge density and thus the
   scattering potential.
3. **Absorbing boundary conditions** in the cell margins. Currently
   reflected electrons accumulate and interfere with the incoming WP,
   masking the true LEED signal at later times. Investigate INQ's
   `perturbations::absorbing_walls` API and add to `shared/configs/`.
4. **Tier-B validation** (energy conservation < 0.1 % over full run,
   restart consistency, CPU/GPU consistency, charge/norm conservation
   within 1e-3, coordinate-mapped LEED matches raw-index-after-shift
   plots). Spec'd in `docs/validation/coronene-replication.md`,
   skipped for this round.
5. **Tier-C convergence sweeps** (dt ∈ {0.020, 0.010, 0.005} a.u.,
   cutoff sweep on the corrected geometry, Tsubonoya 2014 Fig. 2
   quantitative reproduction). Same.
6. **Far-field b-scan trend monotonicity check** — additional b values
   between `run_b6_35x35x80` and `run_b18_35x35x80` to confirm the
   near-field → far-field transition.
7. **Projectile vs paper comparison** — quantitative side-by-side of
   `run_E800` and `run_base` LEED at the paper window.
8. **MPI-aware screen extraction** — `inqkit::screens::PlaneScreen::extract`
   currently assumes single-rank; add `MPI_Allreduce` over the slice
   array before returning. Without this, multi-rank runs silently
   produce incomplete LEED slices.
9. **Transmission vs reflection planes** — paper Fig. 2 is the
   reflection-side pattern at z = +D; a transmission-side plane at
   z = -D would also be physically meaningful. Add a second
   accumulator at the symmetric position.
10. **Legacy run-directory cleanup** — move every legacy buggy run
    (`04_leed_simulation/`, `coronene-wp-rt/`, the old
    `Tutorial/run_diagnoses/run_01..05` rows) under a single
    `coronene/legacy/` subtree. Some of this is already done on the
    user's local; the rest is in `docs/todo_later.md`.
11. **Backfill the 4 missing metadata files** for `run_base`,
    `run_E30`, `run_E800` from each run's `run_summary.txt`. Decision
    in Phase-2 plan: accept the loss; every datum is duplicated.
12. **VESTA visualisation** of the centred coronene geometry to
    confirm the molecular structure visually matches Fig. 1 of the
    paper.
13. **End-to-end smoke test** running `04_leed_simulation` on CPU
    with a reduced parameter set (low cutoff, ~10 steps) to verify
    the full output pipeline produces sensible numbers.
14. **Restart consistency check** — save GS at step 2, reload at step
    5, confirm the first few RT steps reproduce identical results.
15. **GPU/CPU parity check** — same short input on CPU and GPU,
    compare `leed_pattern.txt` element-wise to 1e-5 tolerance.

These are tracked in `docs/todo_later.md`. Each carries a brief rationale
there so future sessions can pick them up without re-deriving the
context.

---

## Where everything lives

| What | Path |
|---|---|
| This handover | `docs/handovers/coronene-cumulative.md` |
| Per-session journal | `docs/handovers/coronene-replication.md` |
| Phase-1 + Phase-2 plan | `docs/plans/coronene-replication.md` (mirror of the live plan in `~/.claude/plans/the-main-aim-of-mutable-taco.md`) |
| Geometry-correction pedagogy | `docs/notes/coronene-geometry-correction.md` |
| Postprocess algorithms | `docs/notes/postprocess-algorithms.md` |
| Tsubonoya 2014 source note | `docs/sources/tsubonoya-2014-coronene-leed.md` |
| Tier-A validation log | `docs/validation/coronene-replication.md` |
| Future tasks | `docs/todo_later.md` |
| Results spec | `docs/results_folder_structure_spec.md` |
| Visualisation rules | `docs/visualisation-instructions-v1.md` |
| New code | `ResearchProject/systems/coronene/{shared/, save_gs/, run_*/, scripts/, hypotheses/}` |
| Generalisable postprocess | `inq-stack/python/inqview/postprocess/` |

**This handover is at:** `/local/data/public/skcb2/tddft/docs/handovers/coronene-cumulative.md`
