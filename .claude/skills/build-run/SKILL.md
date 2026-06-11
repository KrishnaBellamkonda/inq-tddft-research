---
name: build-run
description: Use when building or running an INQ simulation (any `.cpp` under `ResearchProject/`, `Tutorial/`, or `inq-stack/` examples). Covers the `inq-run` GPU/CPU wrapper, when to write a manual `CMakeLists.txt` for non-INQ programs, and the bashrc-sourcing requirement.
---

# Build and run INQ simulations

## Environment prerequisites

`inq-run` and the project tooling (`pvpython`, the `venv` Python, `pyright`, `clangd`)
require these env vars / PATH entries:

```
INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
PATH includes: /local/data/public/skcb2/tddft/shared/bin
               /local/data/public/skcb2/pyenv/shims
               /usr/bin
```

The two `*_SHARE_PATH` vars are pinned in the `env` block of `.claude/settings.json`
(tracked), so the Bash tool sees them on every call. The **PATH** entries
(`shared/bin`, `pyenv/shims`) come from `~/.bashrc`, which the Bash tool's shell
sources — they are deliberately NOT set in settings (overriding PATH there would
shadow system tools). If a command still fails with "not found", fall back to
`bash -lc '<cmd>'` (login shell re-sources `~/.bashrc`).

## Standard INQ workflow

From the directory containing the `.cpp` configuration file:

```bash
inq-run              # auto-detect .cpp, build + run (GPU, CUDA sm_80)
inq-run run.cpp      # explicit filename
inq-run --cpu        # CPU-only (only if GPU unavailable or user requests it)
inq-run --reconfig   # force cmake reconfiguration
```

`inq-run` auto-generates `CMakeLists.txt` — do NOT write one manually for INQ projects.

**Default to GPU.** Only use `--cpu` if the user explicitly asks or no GPU is
available (`nvidia-smi` returns no devices).

## Multi-GPU launches

```bash
CUDA_VISIBLE_DEVICES=0 nohup inq-run > run.log 2>&1 &
CUDA_VISIBLE_DEVICES=1 nohup inq-run > run.log 2>&1 &
```

## Non-INQ standalone C++ programs

For programs that only use FFTW3, BLAS, or similar system libraries, write a
`CMakeLists.txt` and build manually:

```bash
mkdir -p build && cd build
cmake .. && make -j$(nproc)
./run_wp           # or whatever the executable is named
```

Find FFTW3 via pkg-config in `CMakeLists.txt`:

```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)
target_link_libraries(my_target PRIVATE ${FFTW3_LIBRARIES} m)
```

## Where INQ configuration files live

| Project area | Location | Tracked by |
|---|---|---|
| Production research runs | `ResearchProject/systems/<material>/<task>/run.cpp` | main git |
| Learning examples | `Tutorial/<name>/run.cpp` | own git (separate repo) |
| QBall references | `QuantumKickExtension/<system>/...` | own git (separate repo) |

Library headers used by these runs:
- C++: `inq-stack/include/inqkit/<module>/*.hpp`
- Python post-processing: `inq-stack/python/inqview/*.py` (installed via `pip install -e inq-stack/`)

## Coordinate convention for ions (project-wide)

**Always specify ion coordinates in an `.xyz` file with absolute Cartesian
positions in the [-L/2, +L/2] convention, then load it with
`systems::ions::parse(<path>, cell)`.** Do *not* use
`insert_fractional(...)`, do *not* hard-code positions inside the
`run.cpp`, and do *not* use any other input form.

Why this matters:

- INQ's exported real-space field has its origin at `-L/2` (see
  `inqkit::fields::density::total` in
  `inq-stack/include/inqkit/fields/density.hpp`). When the ion file is
  written in [-L/2, +L/2] coordinates, the atom positions in the
  rendered VTI line up directly with the positions in the source `.xyz`
  — what you typed is what you see.
- Mixing fractional coordinates with the FFT-shifted output silently
  parks atoms at the cell corners (and wraps them under PBC), which is
  hard to spot in the integrated electron count and only shows up when
  someone tries to slice the orbitals in ParaView. We have already lost
  one debugging cycle to this on a Nitrogen smoke test.
- Reading geometry from an `.xyz` keeps the geometry in version control
  as a separate, human-readable artefact, and lets two simulations
  (e.g. GS + propagation) share *exactly the same* atomic positions
  without code duplication.

Recommended pattern:

```cpp
auto a = 7.02_angstrom;                           // cell side
auto cell = systems::cell::cubic(a).periodic();
auto ions = systems::ions::parse("geometry.xyz", cell);
```

with `geometry.xyz` next to the `run.cpp`:

```
N_atoms
Comment line: cell L = 7.02 Å, atoms in [-L/2, +L/2] Å
Li  -3.510  -3.510  -3.510
Li  -1.755  -1.755  -1.755
... (one line per atom)
```

For supercells, generate the `.xyz` with a tiny Python script in
`scripts/` and check it into git alongside the `run.cpp`. The script
itself documents the geometric construction (BCC basis, lattice
constant, supercell multiplicity) so future readers can re-derive the
positions.

## GPU pinning

Always pin the run to a specific GPU using `CUDA_VISIBLE_DEVICES`. Before
launching a long run, query `nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv`
and pick the GPU with the lowest utilisation **and** lowest committed
memory (a 0% utilisation device with 20 GB committed is somebody else's
parked job — do not steal it). Record the chosen GPU index in the run's
handover and `run_summary.txt`.

```bash
CUDA_VISIBLE_DEVICES=<idx> inq-run > run.log 2>&1 &
```

Never launch a run without the pin.

## Full reference

See `docs/compilation.md` for the complete `inq-run` pipeline, manual CMake
instructions, and a comparison with the upstream `inq_template`.
