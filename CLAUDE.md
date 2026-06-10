# CLAUDE.md — Core workflow rules

Keep startup context lean. Do not load more tools, agents, MCP servers, or files than the task requires.
Do not spawn agent teams unless the user explicitly requests or approves them.
Keep this file minimal. Put specialised workflows, literature packs, and benchmark procedures into `.claude/skills/` and `.claude/rules/`.

---

## Workflow rules

- Before any substantive implementation, create or update a plan in `docs/plans/`.
- Before stopping, compacting, clearing context, or declaring completion, create or update a handover in `docs/handovers/`.
- Create files only in designated directories. If no suitable directory exists, propose one before creating many new files.
- Scientific claims, modelling choices, and validation claims must be grounded in trustworthy sources or explicitly labelled as uncertain inference.
- Record important sources and attribution notes in plans, handovers, reports, and code comments where relevant.
- No substantive code is complete without recorded validation status.
- For each substantive change, define component tests, integration tests, and scientific benchmark or sanity checks where applicable.
- Suggest test options to the user before running expensive simulations. The user decides which expensive tests to run.
- Never invent paths, APIs, equations, constants, file contents, or test results.
- When uncertain, say uncertain and verify.
- Prefer updating existing docs over creating duplicates.
- Keep plans, handovers, and notes clear, human-readable, and concise.

Extended rules are in `.claude/rules/`. Skills for literature review, simulation validation, report writing, and handover updates are in `.claude/skills/`.

---

## Repository overview

Three tracked components:
- **`inq-stack/`** — project-local library layer on top of INQ:
  - `inq-stack/include/inqkit/` — C++ header library (field extraction, I/O, screens, wavepacket)
  - `inq-stack/python/inqview/` — Python post-processing and visualisation package
- **`ResearchProject/`** — production research experiments (coronene LEED, jellium)
- **`Tutorial/`** — self-contained INQ learning examples

Not tracked (gitignored but present on disk):
- `inq/` — INQ source (header-only C++17, GPU engine — see build section below)
- `QuantumKickExtension/` — QBall reference calculations (Li, diamond, Al)
- `venv/` — Python virtual environment
- `ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/` — ParaView installation

---

## inq-stack library (inqkit + inqview)

### inqkit (C++ headers, `inq-stack/include/inqkit/`)

INQ-facing extraction layer. Include headers alongside `<inq/inq.hpp>`.

| Module | Contents |
|---|---|
| `fields/` | `RealField3D`, `ComplexField3D`, `density.hpp`, `orbital.hpp` |
| `io/` | `RealField3DWriter`, `ComplexField3DWriter`, VTI writer, manifest, observables |
| `detail/` | `grid_layout`, `filesystem`, `text_io`, `validation` |
| `config/` | `simulation_config` |
| `wavepacket/` | `wavepacket.hpp`, `injection_report.hpp` |
| `screens/` | `plane_screen`, `leed_pattern_accumulator` (stubs, in-progress) |
| `core/` | `pipeline`, `session_context`, `task` (stubs, in-progress) |
| `jellium/` | `analytics.hpp` |

Usage example: `Tutorial/n2-with-inqkit/run.cpp`

### inqview (Python, `inq-stack/python/inqview/`)

Post-processing and visualisation. Install into the venv with `pip install -e inq-stack/`.
Restructured into **four dependency-layered sub-packages** (ADR 0003); the
top-level `inqview` re-exports the public names **lazily** (PEP 562) so importing
the deps-clean layers pulls no matplotlib/VTK.

| Sub-package | Role | May import |
|---|---|---|
| `inqview.io` | loaders + field/format dataclasses (`fields`, `data`, `leed`) | numpy only |
| `inqview.analysis` | numeric kernels → frozen dataclasses (`fourier`, `energy_components`, `wp_integrity`, `plasmon_spectrum`, `center_of_density`, `kl_divergence`) | numpy/scipy/pandas |
| `inqview.visualisation` | all rendering (`plots`, `paraview`, `vti`, `style`, `carpets`, renderers) + the canonical theme (ADR 0004) | matplotlib/VTK |
| `inqview.pipeline` | thin phase orchestration (34 phases + `runner`/`frames`/`cod`) | the above |

`inqview.postprocess` is a **deprecated back-compat shim** forwarding to
`inqview.pipeline` (existing run `analyse.py` import the old path). Tests:
`inq-stack/python/tests/` (portable, numpy-only). Deps-clean invariant enforced by
`tests/test_deps_clean.py`. See `docs/handovers/inqview-restructure.md`.

---

## Build and run (INQ configuration files)

The standard workflow uses the `inq-run` wrapper. From any directory containing a `.cpp` file:

```bash
inq-run              # auto-detect .cpp, build + run (GPU, CUDA sm_80)
inq-run --cpu        # CPU-only
inq-run --reconfig   # force cmake reconfigure
```

Prerequisites (add once to `~/.bashrc`):
```bash
export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
```

Full compilation pipeline: see `docs/compilation.md`.

---

## INQ library build (rebuild from source)

```bash
cd inq
mkdir build && cd build
cmake .. --install-prefix=$HOME \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/nvcc
cmake --build . --parallel
make install
```

Tests:
```bash
cd inq/build
ctest --output-on-failure --timeout 2000
INQ_EXEC_ENV="mpirun.openmpi -np 4" ctest --output-on-failure --timeout 2000
```

**Do not revert the CUB fix** in `inq/external_libs/gpurun/include/gpu/reduce.hpp` — required for CUDA 12.5+.

---

## Key documentation

| File | Contents |
|---|---|
| `docs/compilation.md` | Full `inq-run` pipeline, manual CMake, template comparison |
| `docs/inq_tutorial.md` | INQ C++ API reference (cell, electrons, theory, TDDFT, perturbations) |
| `docs/inq_source_map.md` | INQ source code module map and extension guide |
| `docs/folder_structure.md` | Directory layout (ResearchProject, Tutorial, inq-stack) |
| `docs/handovers/coronene_wp_scattering.md` | Coronene LEED — current active task, run_004 complete |
| `docs/handovers/jellium_ground_state.md` | Jellium ground state — complete |

---

## INQ architecture (quick reference)

INQ is **header-only C++17**. All source is under `inq/src/`.

| Module | Role |
|---|---|
| `systems/` | `ions`, `electrons`, `cell` — primary data containers |
| `hamiltonian/` | Kohn-Sham H, XC, pseudopotentials, PAW |
| `ground_state/` | SCF driver (`calculator.hpp`) |
| `real_time/` | TDDFT propagation (`propagate.hpp`) |
| `operations/` | Field operators: overlap, FFT, gradient, I/O |
| `observables/` | Density, forces, dipole, current |
| `perturbations/` | Kick, laser, electric field, absorbing potential |
| `parallel/` | MPI 3D Cartesian decomposition |
| `interface/` | CLI command parsing |

Data flow: `interface/` → `systems/` → `hamiltonian/` → `ground_state/` or `real_time/` → `observables/`

See `docs/inq_source_map.md` for the full module breakdown.

---

## QuantumKickExtension (not tracked)

`QuantumKickExtension/codebase/` — QBall reference calculations (Li, diamond, Al).
Input files: `.inp`, `.sys`, `.vel`, `.xml`. Gitignored; present on disk only.
Reference: Santervás-Arranz, Stengel, Artacho, Phys. Rev. Research 7, 033292 (2025).
