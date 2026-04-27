# Handover: coronene-replication

Rolling handover for the coronene WP-RT-LEED replication framework. See
`docs/plans/coronene-replication.md` for the full plan and the long-form
inventory of the legacy buggy runs.

---

## Milestone: 2026-04-26 21:31 — Queue dispatching across both GPUs

### Current status

10 runs to do (3 GS saves + 9 propagations beyond `run_base`). Status:

| Run | State |
|---|---|
| save_gs/gs_35x35x60_cut40 | ✅ done — GS energy −150.837 Ha |
| save_gs/gs_35x35x80_cut40 | ✅ done — walltime 541 s |
| save_gs/gs_35x35x40_cut40 | ✅ done — walltime 451 s |
| run_base                  | ✅ done — walltime 2851 s; 560 steps; full postprocess produced |
| run_E30                   | 🔄 GPU 1 (since 21:14) — 1446 steps |
| run_E800                  | 🔄 GPU 0 (since 21:21) — 280 steps |
| run_s0p33                 | 📋 queued |
| run_s3                    | 📋 queued |
| run_E800_s0p33            | 📋 queued |
| run_E30_s3                | 📋 queued |
| run_b18_35x35x80          | 📋 queued |
| run_b6_35x35x80           | 📋 queued |
| run_35x35x40              | 📋 queued |

Both GPUs saturated. Dispatcher logs at `scripts/dispatch.log`; auto-
postprocess watcher logs at `scripts/auto_postprocess.log`.

### Bugs found and fixed during smoke test

1. **`std::sqrt` not constexpr under CUDA.** The base config used
   `static constexpr double WP_K0 = std::sqrt(...)` which the nvcc front-end
   rejected. Replaced with a Newton-iteration `const_sqrt(double)` defined
   in `tsubonoya_2014_base.hpp`.
2. **Empty wavepacket / screen-config metadata files.** The `*_path()`
   helpers in `shared/cpp/results_paths.hpp` did not call `ensure_dir` on
   their parent. So `wavepacket_config.txt`, `injection_report.txt`,
   `screen_config.csv`, `window_ranges.csv` were silently dropped. Fixed by
   adding `ensure_parent(path)` to every `*_path()` helper. The 3 already-
   in-flight runs (`run_base`, `run_E30`, `run_E800`) miss those four files
   each, but every datum is duplicated in `run_summary.txt` or recoverable
   from screen filenames; the missing files don't affect the postprocess.
3. **Postprocess `density` phase silent no-op.** Glob expected
   `{cat}_t*.vti` but the C++ writer uses the layout's `field_name`
   (`density_t*.vti`). Switched `_common.list_vti_series()` to glob
   `*_t*.vti`.
4. **Postprocess `gs` phase no-op.** Looked only for raw + meta sidecars;
   the new C++ template emits VTI only. Added VTI fallback in
   `inqview/postprocess/ground_state.py` (loads with vtkXMLImageDataReader,
   plots xy mid-plane and orbital gallery).

### Postprocess artefacts confirmed for run_base

`run_base/results/analysis/`:

* `ground_state/` — `density_gs_system_xy.png`, `gs_orbital_gallery.png`.
* `observables/` — 11 PNGs (energy/current/dipole vs time, FFT spectra).
* `density/` — 9 GIFs (total/system/wp × xy/xz/yz), fixed colour scale.
* `screens/` — 20 total panels (linear+log) + 4×5 grid; 20 instantaneous
  GIFs; 60 time-windowed PNGs (forward/back/paper × log+linear); 2
  coordinate-check plots.
* `overlap/wp_overlap_with_gs_orbitals.gif`.

### Dependencies

`vtk 9.6.1` installed into project venv (`/local/data/public/skcb2/tddft/venv`)
during postprocess validation. `imageio`, `matplotlib`, `pandas`, `numpy`
already present.

### Next steps when queue completes

1. Watcher auto-postprocesses each run as it finishes; check
   `scripts/auto_postprocess.log` for failures and re-run any that failed.
2. Once the last propagation finishes, run the 6 hypothesis comparisons:
   ```
   for h in 00_base 01_wp_energy_spread 02_wp_sigma_spread \
            03_fast_projectile_classical 04_electron_capture \
            05_box_length_and_distance; do
     python3 scripts/coronene_postprocess.py hypothesis \
         --hypothesis-dir hypotheses/$h \
         --runs $(grep '^run_' hypotheses/$h/README.md | ...)
   done
   ```
   The argument list is best constructed by hand from each hypothesis
   README.
3. Update `docs/validation/coronene-replication.md` with Tier-A pass/fail
   per run.

---

## Milestone: 2026-04-26 — Code framework in place; awaiting user review

### Current status

All implementation code is on disk under
`/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/`. **No
simulation has been launched yet.** The dispatcher exists and has a
`--dry-run` mode; the user is reviewing before authorising any GPU runs.

The framework consists of:

* shared/ — geometry (canonical z = 0 coronene.xyz), 11 Cfg headers, and
  three C++ helpers (`results_paths.hpp`, `leed_screen_layout.hpp`,
  `run_template.hpp`).
* save_gs/ — three GS save runs (one per (cell, cutoff) tuple).
* 10 propagation run.cpp files at the top level: `run_base`, `run_E30`,
  `run_E800`, `run_s0p33`, `run_s3`, `run_E800_s0p33`, `run_E30_s3`,
  `run_b18_35x35x80`, `run_b6_35x35x80`, `run_35x35x40`. Each is a thin
  wrapper that templates `coronene::run_template::run_propagation<Cfg>`.
* hypotheses/ — six READMEs (00_base..05_box_length_and_distance) that
  declare which `run_*/results/` trees feed into each comparison.
* scripts/coronene_postprocess.py — argparse CLI with `run` and
  `hypothesis` subcommands.
* scripts/dispatch_runs.py — Python multi-GPU dispatcher with `--dry-run`,
  `--gpus`, `--free-mem-mb` (handles the user's "one GPU may already be
  busy" case via nvidia-smi polling).
* scripts/run_queue.txt — fixed launch order: 3 GS saves, then 10
  propagations grouped by hypothesis.
* inq-stack/python/inqview/postprocess/ — generalisable phase modules
  (pipeline, run_summary, ground_state, observables, density, screens,
  overlap, orbitals, compare).

### What changed

- Added `ResearchProject/systems/coronene/shared/` (configs + cpp helpers + xyz).
- Added 10 `run_*/run.cpp` propagation drivers + 3 `save_gs/<sig>/run.cpp` GS savers.
- Added `inq-stack/python/inqview/postprocess/` subpackage.
- Added `ResearchProject/systems/coronene/scripts/{coronene_postprocess.py,
  dispatch_runs.py, run_queue.txt}`.
- Added `ResearchProject/systems/coronene/hypotheses/00_base..05_box_length_and_distance/README.md`.
- Mirrored plan to `docs/plans/coronene-replication.md`.
- Added `docs/sources/tsubonoya-2014-coronene-leed.md`.
- Appended legacy-cleanup entry to `docs/todo_later.md`.
- Verified `inqkit::observables::OrbitalOverlapMatrix::snapshot_wp_only()`
  already exists at `inq-stack/include/inqkit/observables/orbital_overlap.hpp:90`,
  so no new C++ observable was added (plan task #4 closed as no-op).

### Files touched (absolute paths)

- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/configs/{tsubonoya_2014_base,E30,E800,s0p33,s3,E800_s0p33,E30_s3,cell_35x35x80,b18_35x35x80,b6_35x35x80,cell_35x35x40}.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/cpp/{results_paths.hpp,leed_screen_layout.hpp,run_template.hpp}`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/save_gs/{gs_35x35x60_cut40,gs_35x35x80_cut40,gs_35x35x40_cut40}/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/{run_base,run_E30,run_E800,run_s0p33,run_s3,run_E800_s0p33,run_E30_s3,run_b18_35x35x80,run_b6_35x35x80,run_35x35x40}/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/scripts/{coronene_postprocess.py,dispatch_runs.py,run_queue.txt}`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/hypotheses/00_base..05_box_length_and_distance/README.md`
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/{__init__,pipeline,_common,run_summary,ground_state,observables,density,screens,overlap,orbitals,compare}.py`
- `/local/data/public/skcb2/tddft/docs/plans/coronene-replication.md`
- `/local/data/public/skcb2/tddft/docs/sources/tsubonoya-2014-coronene-leed.md`
- `/local/data/public/skcb2/tddft/docs/todo_later.md` (appended)
- `/local/data/public/skcb2/tddft/docs/validation/coronene-replication.md`

### Commands run

Read-only API exploration only. No builds, no GPU launches.

```bash
# Listing/reading commands; no state-changing commands run.
ls /local/data/public/skcb2/tddft/inq-stack/include/inqkit
python3 -c "import ast; ast.parse(...)"   # syntax-check the new Python modules
```

### Tests and validation

- **Proposed**: Tier-A only, on representative runs (every GS save plus a
  results-tree check on every propagation). Per the plan §5, with the new
  GPU-execution check.
- **Approved**: pending user review of the code.
- **Run**: none yet.
- **Outcomes**: n/a.
- **Remaining gaps**: actual SCF / propagation / postprocess validation;
  the new postprocess modules in inqview/postprocess/ have not been
  exercised against any real `results/` tree.

### Trusted sources used

- Tsubonoya, Hu & Watanabe, *Phys. Rev. B* **90**, 035416 (2014). See
  `docs/sources/tsubonoya-2014-coronene-leed.md`.

### Attribution notes

- `shared/cpp/run_template.hpp` is heavily adapted from
  `ResearchProject/systems/coronene/run_propagate_paper_replica/run.cpp` —
  rewritten to (a) emit the spec-compliant `results/` tree, (b) write three
  density categories instead of two, (c) flat snapshot filenames, (d) use
  `snapshot_wp_only` (per the plan §4.6 item 9).
- `save_gs/*/run.cpp` adapted from
  `ResearchProject/systems/coronene/run_save_gs_paper_replica/run.cpp`.

### Known issues / blockers

- The C++ code has not been compiled yet. Likely a few rounds of small
  fixes will surface on first build (template syntax, header name
  collisions). Plan: build `save_gs/gs_35x35x60_cut40/` first as smoke
  test; iterate fixes there before building the rest.
- The postprocess `density` phase requires the `vtk` Python package for
  loading binary VTI; `imageio` for writing GIFs. Both are expected in the
  project venv but have not been confirmed.
- The dispatcher's "GPU initialisation banner" check searches for `"GPU"`
  / `"CUDA"` substrings in `run.log`; this is loose. If INQ's banner text
  changes the heuristic will produce false negatives. Tighten once the
  exact banner string is observed.

### Assumptions still in play

- **INQ XYZ parser uses absolute coordinates as-is** (verified from
  `inq/src/parse/xyz.hpp:57`); the same z = 0 coronene.xyz works for any
  orthorhombic cell whose half-extent contains the molecule and the WP
  start point.
- **`inq-run` is on PATH after `source ~/.bashrc`** — encoded in the
  dispatcher's pre-flight check.
- **Two GPUs are functionally identical** (user-confirmed) — no
  VRAM-aware pinning needed.
- **Cutoff 40 Ha is converged on dojo PSPs**. Re-confirmation when the
  new `cell_35x35x80` and `cell_35x35x40` GSes are computed (compare GS
  energies against the existing `01_geometry/` results scaled for cell
  difference).
- **The paper's `54 Ha` cutoff is unstable on dojo PSPs**; dropping to
  40 Ha is a deliberate departure recorded in the source note.

### Exact next steps

1. **User reviews the code** under
   `ResearchProject/systems/coronene/{shared,scripts,save_gs,run_*,hypotheses}`
   and `inq-stack/python/inqview/postprocess/`.
2. On user approval, build the framework starting with the cheapest GS
   save run as a smoke test:
   ```
   cd ResearchProject/systems/coronene/save_gs/gs_35x35x60_cut40
   source ~/.bashrc
   inq-run run.cpp
   ```
   Iterate compile fixes until the GS save completes and writes a
   checkpoint.
3. Once the first GS save runs cleanly, build a second GS save
   (`gs_35x35x80_cut40`) and `run_base` end-to-end, exercise
   `coronene_postprocess.py run --results run_base/results --rebuild`.
4. Run a `dispatch_runs.py --dry-run` to confirm the queue maps to the
   intended GPU launches.
5. With user permission, launch the full queue
   (`dispatch_runs.py scripts/run_queue.txt`).
6. After completion, run `coronene_postprocess.py hypothesis` for each of
   the 6 hypotheses folders to populate the comparison artefacts.
7. Update this handover with results / observed walltimes / drift values.
