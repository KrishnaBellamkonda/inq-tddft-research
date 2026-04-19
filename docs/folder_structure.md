# tddft — Folder Structure and Utilities

## Top-level Layout

```
tddft/
├── docs/                    ← guides and reference docs
│   ├── compilation.md       ← full compilation pipeline + inq-run + template comparison
│   ├── folder_structure.md  ← this file
│   ├── inq_tutorial.md      ← INQ C++ API reference
│   ├── inq_source_map.md    ← INQ codebase module map (for extension/modification)
│   ├── jellium_gaussian_wave_packet.md ← validation example
│   ├── plans/               ← active task plans (one per substantive task)
│   ├── handovers/           ← session continuation notes and milestone summaries
│   ├── sources/             ← literature notes, source summaries, attribution
│   ├── validation/          ← test matrices, benchmark definitions, comparison plots
│   ├── reports/             ← report drafts and manuscript fragments
│   ├── notes/               ← temporary working notes
│   └── inq-docs/            ← local mirror of alphataubio.com/inq documentation
├── plans/                   ← legacy root-level notes (new work → docs/plans/)
├── shared/                  ← utilities shared across all projects
│   ├── config.sh            ← INQ paths and CUDA settings (edit once)
│   └── bin/
│       └── inq-run          ← one-command build-and-run for any .cpp
├── inq/                     ← INQ source (gitignored; present on disk; header-only C++17, GPU-capable)
├── inq-stack/               ← project-local library on top of INQ (tracked)
│   ├── include/inqkit/      ← C++ header library
│   │   ├── fields/          ← RealField3D, ComplexField3D, density, orbital
│   │   ├── io/              ← field writers, VTI writer, manifest, observables
│   │   ├── detail/          ← grid_layout, filesystem, text_io, validation
│   │   ├── config/          ← simulation_config
│   │   ├── wavepacket/      ← wavepacket.hpp, injection_report
│   │   ├── screens/         ← plane_screen, leed_pattern_accumulator (stubs)
│   │   ├── core/            ← pipeline, session_context, task (stubs)
│   │   └── jellium/         ← analytics
│   ├── python/inqview/      ← Python post-processing and visualisation package
│   │   ├── fields.py        ← RealField3D, ComplexField3D, FieldMeta
│   │   ├── data.py          ← SimulationData, FieldSeries, loaders
│   │   ├── vti.py           ← VTI series conversion (for ParaView)
│   │   ├── paraview.py      ← ParaViewPipeline batch renderer
│   │   ├── plots.py         ← matplotlib helpers
│   │   ├── fourier.py       ← k-space analysis
│   │   └── config.py        ← Theme, PlotDefaults, RenderDefaults
│   └── pyproject.toml       ← pip install -e inq-stack/
├── Tutorial/                ← learning examples (start here)
│   ├── hello-world/         ← H₂ ground state
│   ├── n2/                  ← N₂ ground state + dipole kick TDDFT
│   ├── n2-with-inqkit/      ← N₂ using inqkit field extraction (integration example)
│   ├── li-bcc/              ← Li BCC metal: SCF + ionic kick TDDFT
│   ├── HF/                  ← HF molecule
│   ├── HF-toy-model/        ← HF laser perturbation toy model
│   ├── ions/                ← ionic dynamics
│   └── minimal-gs/          ← Si minimal ground state
├── ResearchProject/         ← production research runs
│   ├── jellium/             ← jellium HEG studies (tracked)
│   │   ├── 01_ground_state/       ← LDA SCF, all 6 tests PASSED
│   │   ├── 02_ground_state_convergence/ ← grid spacing + shell closure
│   │   └── 03_free_gaussian_wp_propagation/ ← standalone FFTW3, validated
│   ├── systems/
│   │   ├── coronene/        ← coronene C₂₄H₁₂ TDDFT LEED simulation (main active)
│   │   │   ├── 01_geometry/
│   │   │   ├── 02_ground_state_analysis/
│   │   │   ├── 03_ecut_convergence/
│   │   │   ├── 04_leed_simulation/   ← active; runs/ subdirs for run_001…005
│   │   │   └── wp_scattering/
│   │   └── jellium/         ← legacy jellium path (older; ResearchProject/jellium/ is canonical)
│   └── literature/tddft/    ← reference epub
└── QuantumKickExtension/    ← QBall reference calculations (gitignored; present on disk)
    └── codebase/            ← .inp/.sys/.vel files for Li, diamond, Al
```

Each "run" directory is completely self-contained:
- **`run.cpp`** — the INQ configuration (the only file you edit)
- **`build/`** — cmake build tree (created automatically; do not commit)
- **`results/`** — output files (`.dat`, saved ground states, etc.)

---

## Documentation Files

| File | Contents |
|---|---|
| `docs/compilation.md` | End-to-end compilation pipeline, `inq-run` internals, comparison with `inq_template` |
| `docs/inq_tutorial.md` | INQ C++ API reference: cell, ions, electrons, theory, ground state, TDDFT, perturbations |
| `docs/inq_source_map.md` | INQ source code map: all 27 modules, data types, data flows, extension points |
| `docs/folder_structure.md` | This file: directory layout and utility reference |
| `docs/jellium_gaussian_wave_packet.md` | Validation example: jellium + Gaussian wave packet |
| `docs/handovers/coronene_wp_scattering.md` | Active task: coronene LEED — run_004 complete, next steps recorded |
| `docs/handovers/jellium_ground_state.md` | Complete: jellium ground state + convergence studies |
| `docs/plans/coronene_wp_scattering.md` | Full coronene LEED plan with validated parameter table |
| `docs/plans/jellium_research.md` | Jellium research plan (experiments 01–02 done; 02_kick and beyond pending) |

---

## inq-stack library

### inqkit (C++, `inq-stack/include/inqkit/`)

Include alongside `<inq/inq.hpp>`. Provides field extraction and I/O independent of
INQ's internal representation.

```
fields/      → RealField3D, ComplexField3D, density.hpp, orbital.hpp
io/          → RealField3DWriter, ComplexField3DWriter, VTI writer, manifest, observables, text_summary
detail/      → grid_layout, filesystem, text_io, validation
config/      → simulation_config
wavepacket/  → wavepacket.hpp, injection_report.hpp
screens/     → plane_screen, leed_pattern_accumulator (stubs, in development)
core/        → pipeline, session_context, task (stubs, in development)
jellium/     → analytics.hpp
```

Working integration example: `Tutorial/n2-with-inqkit/run.cpp`

### inqview (Python, `inq-stack/python/inqview/`)

Install with: `pip install -e inq-stack/` (into `venv/`).

```
fields.py    → FieldMeta, RealField3D, ComplexField3D
data.py      → SimulationData, FieldSeries, load_real_field, load_complex_field
vti.py       → convert_real_series_to_vti, write_vti
paraview.py  → ParaViewPipeline (batch render via pvbatch)
plots.py     → matplotlib helpers
fourier.py   → k-space analysis
config.py    → Theme, PlotDefaults, RenderDefaults
scripts/     → CLI entry points for batch rendering
```

Active development branch: `features/python-paraview`

---

## Utilities

### `shared/config.sh`

Central configuration file. **Edit this** when paths change.

```bash
export INQ_SOURCE="/local/data/public/skcb2/tddft/inq"
export INQ_SHARE_PATH="$INQ_SOURCE/install/share"
export PSEUDOPOD_SHARE_PATH="$INQ_SOURCE/install/share/pseudopod"

export INQ_CUDA_COMPILER="/lsc/opt/cuda-12.6.2/bin/nvcc"
export INQ_CUDA_ARCH="80"            # A30 = sm_80
export INQ_PYTHON_EXE="/local/data/public/skcb2/tddft/venv/bin/python3"
```

### `shared/bin/inq-run`

One-command build-and-run for any self-contained `.cpp` file.
Add to `PATH` via `~/.bashrc`:

```bash
export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
```

**Usage (run from the directory containing the `.cpp` file):**

```bash
inq-run                  # auto-detects the single .cpp, builds + runs (GPU)
inq-run --cpu            # same but CPU-only build (no CUDA)
inq-run --reconfig       # force cmake reconfiguration (e.g. after path change)
inq-run myfile.cpp       # explicit file name
inq-run --help           # show usage
```

**What it does:**
1. Sources `shared/config.sh` for all INQ/CUDA paths
2. Generates a fresh `build/CMakeLists.txt` pointing to the INQ source tree
3. Runs `cmake .` with CUDA flags (unless `--cpu`) — skipped if already configured
4. Builds with `cmake --build . -j$(nproc)`
5. Runs the binary with `INQ_SHARE_PATH` and `PSEUDOPOD_SHARE_PATH` set

**First run** takes ~2–5 minutes (cmake downloads dependencies). Subsequent runs
only recompile if the `.cpp` changed — typically seconds.

---

## Tutorial Build (Alternative to `inq-run`)

The `Tutorial/` directory has a single shared cmake build tree. This is faster
for iterating on multiple tutorial targets but less self-contained.

```bash
# Initial configure (one time)
cd /local/data/public/skcb2/tddft/Tutorial
mkdir -p build && cd build
cmake .. \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/nvcc \
  -DPython_EXECUTABLE=/local/data/public/skcb2/tddft/venv/bin/python3

# Build a specific target
cmake --build . --target h2 -j$(nproc)    # hello-world
cmake --build . --target n2 -j$(nproc)    # n2 TDDFT
cmake --build . --target li_bcc -j$(nproc) # Li BCC TDDFT

# Run (binary lands next to .cpp)
INQ_SHARE_PATH=... PSEUDOPOD_SHARE_PATH=... ./hello-world/h2
# or simply (if ~/.bashrc exports the paths):
./hello-world/h2
```

---

## Environment Variables (add to `~/.bashrc` once)

```bash
# INQ pseudopotential lookup paths
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod

# inq-run utility
export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
```

---

## Key Fix Applied to INQ Source (do not revert)

**File:** `inq/external_libs/gpurun/include/gpu/reduce.hpp`

CUDA ≥ 12.5 ships CUB 2.4+ which evaluates `invoke_result<TransformOp, T>` for
device lambdas — this fails without an explicit return type. The fix wraps the
kernel with `cuda::proclaim_return_type<Type>(kernel)` and uses `std::plus<Type>{}`
instead of `std::plus<>{}`. Without this fix INQ will not compile on A30/CUDA 12.6.2.
