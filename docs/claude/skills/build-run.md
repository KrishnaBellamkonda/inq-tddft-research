# Skill: Build and Run INQ Simulations

## Prerequisites

`inq-run` requires environment variables set in `~/.bashrc`. In a Claude Code
session these are **not loaded automatically** — source the file first:

```bash
source ~/.bashrc
```

This adds `inq-run` to `PATH` and sets `INQ_SHARE_PATH` / `PSEUDOPOD_SHARE_PATH`.

## Standard workflow

From the directory containing the `.cpp` configuration file:

```bash
source ~/.bashrc
inq-run              # auto-detect .cpp, build + run (GPU, CUDA sm_80)
inq-run run.cpp      # explicit filename
inq-run --cpu        # CPU-only (only if GPU unavailable or user requests it)
inq-run --reconfig   # force cmake reconfiguration
```

`inq-run` auto-generates `CMakeLists.txt` — do NOT write one manually for INQ projects.

## Non-INQ C++ programs (standalone, no INQ dependency)

For programs that only use FFTW3, BLAS, or similar system libraries, write a
`CMakeLists.txt` and build manually:

```bash
source ~/.bashrc
mkdir -p build && cd build
cmake .. && make -j$(nproc)
./run_wp           # or whatever the executable is named
```

Find FFTW3 via pkg-config in CMakeLists.txt:

```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)
target_link_libraries(my_target PRIVATE ${FFTW3_LIBRARIES} m)
```

## Full reference

See `docs/compilation.md` for the complete `inq-run` pipeline, manual CMake
instructions, and a comparison with the upstream `inq_template`.
