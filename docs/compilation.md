# INQ Compilation Guide

This document explains the complete compilation pipeline for INQ `.cpp` configuration
files, the `inq-run` wrapper in detail, and a comparison with the upstream
[`inq_template`](https://gitlab.com/npneq/inq_template) project.

---

## 1. Background: What Is Being Compiled?

INQ is a **header-only C++17 library** (`#include <inq/inq.hpp>`). A user writes a
`.cpp` file (their "configuration file") that calls the INQ API to define a system,
run a DFT or TDDFT calculation, and save results.

This `.cpp` file becomes the `main()` of a standalone executable.  
The executable links against:
- The INQ header library (compiled inline — all INQ logic ends up in your binary)
- CUDA runtime (cublas, cufft, cudart, cusolver) — GPU kernels
- MPI (OpenMPI) — parallel communication
- FFTW — CPU Fourier transforms (used for non-GPU paths)
- BLAS/LAPACK — dense linear algebra
- External dependencies: pseudopod, libxc, spglib (all fetched by CMake)

The binary is self-contained: no shared INQ library is needed at runtime, only the
environment variables pointing to pseudopotential data.

---

## 2. The `inq-run` Wrapper

**Location:** `/local/data/public/skcb2/tddft/shared/bin/inq-run`

`inq-run` is a bash script that encapsulates the full CMake configure+build+run
pipeline into a single command. It requires no arguments in the typical case.

### 2.1 Prerequisites

Add once to `~/.bashrc`:

```bash
# Make inq-run available on PATH
export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"

# Pseudopotential lookup at runtime
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
```

### 2.2 Usage

Run from the directory containing the `.cpp` file:

```bash
inq-run                  # auto-detect the single .cpp, build + run (GPU)
inq-run run.cpp          # explicit filename
inq-run --cpu            # CPU-only build (no CUDA)
inq-run --reconfig       # force cmake reconfiguration (e.g. after path change)
inq-run --help           # print usage summary
```

### 2.3 What `inq-run` Does — Step by Step

**Step 1 — Source shared config**

```bash
source /local/data/public/skcb2/tddft/shared/config.sh
```

This sets:

| Variable | Value |
|---|---|
| `INQ_SOURCE` | `/local/data/public/skcb2/tddft/inq` |
| `INQ_SHARE_PATH` | `$INQ_SOURCE/install/share` |
| `PSEUDOPOD_SHARE_PATH` | `$INQ_SOURCE/install/share/pseudopod` |
| `INQ_CUDA_COMPILER` | `/lsc/opt/cuda-12.6.2/bin/nvcc` |
| `INQ_CUDA_ARCH` | `80` (Ampere / A30) |
| `INQ_PYTHON_EXE` | `/local/data/public/skcb2/tddft/venv/bin/python3` |

**Step 2 — Auto-detect the `.cpp` file**

If no filename is given, the script scans the current directory for exactly one
`.cpp` file. If zero or more than one is found, it exits with an error.

The target name (binary name) is derived from the `.cpp` filename without extension:
- `run.cpp` → binary named `run`
- `gaussian_wave_packet.cpp` → binary named `gaussian_wave_packet`

**Step 3 — Generate `build/CMakeLists.txt`**

Written fresh on every call (so path changes in `config.sh` are always picked up):

```cmake
cmake_minimum_required(VERSION 3.21)

option(ENABLE_CUDA "Enable CUDA GPU support" OFF)

if(ENABLE_CUDA)
    project(<target> CXX Fortran CUDA)
    set(ENABLE_GPU ON)
    set(GPU_LANGUAGE CUDA)
else()
    project(<target> CXX Fortran)
endif()

# Pull in INQ (headers + linked libraries).
# EXCLUDE_FROM_ALL means only the targets needed here are built.
add_subdirectory(/local/data/public/skcb2/tddft/inq inq EXCLUDE_FROM_ALL)

add_executable(<target> /absolute/path/to/run.cpp)
target_link_libraries(<target> PRIVATE inq)

# Place the binary in the run directory, next to the .cpp file.
set_target_properties(<target> PROPERTIES RUNTIME_OUTPUT_DIRECTORY /run/dir)

if(ENABLE_GPU)
    set_source_files_properties(run.cpp PROPERTIES LANGUAGE CUDA)
endif()
```

Key points:
- `add_subdirectory(INQ_SOURCE)` uses the **pre-built local INQ** — no download needed.
- `EXCLUDE_FROM_ALL` avoids compiling INQ tests and benchmarks.
- The `LANGUAGE CUDA` property tells NVCC to compile the `.cpp` file.
- The binary lands in the **run directory** (next to the `.cpp`), not inside `build/`.

**Step 4 — CMake configure**

Runs only on first call (or when `--reconfig` is passed):

```bash
# GPU build (default)
cmake . \
    -DENABLE_CUDA=ON \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/nvcc \
    -DPython_EXECUTABLE=/local/data/public/skcb2/tddft/venv/bin/python3

# CPU-only build (--cpu flag)
cmake . \
    -DENABLE_CUDA=OFF \
    -DPython_EXECUTABLE=/local/data/public/skcb2/tddft/venv/bin/python3
```

On **first run** this takes ~1–5 minutes: CMake walks the INQ dependency tree and
may download or verify external packages (pybind11, spdlog, catch2, spglib, etc.)
via FetchContent.

On **subsequent runs** CMake finds `CMakeCache.txt` and skips reconfiguration.

**Step 5 — Build**

```bash
cmake --build . --target <target> -j$(nproc)
```

Uses all available CPU cores. Only recompiles if the `.cpp` file changed.
Incremental builds: typically a few seconds.

**Step 6 — Run**

```bash
cd /run/dir
exec env \
    INQ_SHARE_PATH="..." \
    PSEUDOPOD_SHARE_PATH="..." \
    ./<target>
```

The binary is executed in the run directory with the pseudopotential environment
variables set explicitly, so results write to the correct locations.

### 2.4 Directory Layout After First Run

```
/path/to/run/
├── run.cpp             ← your configuration file (only file you edit)
├── run                 ← compiled binary (placed here by CMake)
├── build/
│   ├── CMakeLists.txt  ← auto-generated by inq-run (overwritten each call)
│   ├── CMakeCache.txt  ← CMake cache (persists; skip reconfigure on next run)
│   ├── inq/            ← CMake build subtree for INQ
│   └── ...             ← object files, dependency tracking
└── results/            ← auto-created; your output files go here
    ├── tddft.dat
    ├── gs_save/
    └── ...
```

---

## 3. Tutorial Shared-Build Alternative

The `Tutorial/` directory uses a single shared CMake build tree instead of per-run
`build/` directories. This avoids duplicating the INQ CMake configuration for each
example but is less self-contained.

```bash
# Initial configure (once)
cd /local/data/public/skcb2/tddft/Tutorial
mkdir -p build && cd build
cmake .. \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/nvcc \
  -DPython_EXECUTABLE=/local/data/public/skcb2/tddft/venv/bin/python3

# Build a specific tutorial target
cmake --build . --target h2     -j$(nproc)   # hello-world (H2)
cmake --build . --target n2     -j$(nproc)   # N2 TDDFT
cmake --build . --target li_bcc -j$(nproc)   # Li BCC TDDFT

# Run (binary is placed next to the .cpp)
./hello-world/h2
./n2/n2
./li-bcc/li_bcc
```

Use this approach when iterating on multiple tutorial files and you want to avoid
repeated CMake configuration overhead.

---

## 4. Running Without `inq-run`: Manual CMake

If `inq-run` is unavailable (e.g. you are on a different machine), the equivalent
manual sequence is:

```bash
# 1. Create and enter a build directory
mkdir build && cd build

# 2. Write CMakeLists.txt (or copy from Tutorial/)
cat > CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.21)
project(myrun CXX Fortran CUDA)
set(ENABLE_GPU ON)
set(GPU_LANGUAGE CUDA)

add_subdirectory(/path/to/inq inq EXCLUDE_FROM_ALL)
add_executable(myrun /path/to/myrun.cpp)
target_link_libraries(myrun PRIVATE inq)
set_target_properties(myrun PROPERTIES RUNTIME_OUTPUT_DIRECTORY /path/to/rundir)
set_source_files_properties(/path/to/myrun.cpp PROPERTIES LANGUAGE CUDA)
EOF

# 3. Configure
cmake . \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/nvcc

# 4. Build
cmake --build . --target myrun -j$(nproc)

# 5. Run
cd /path/to/rundir
INQ_SHARE_PATH=... PSEUDOPOD_SHARE_PATH=... ./myrun
```

---

## 5. Comparison: `inq-run` vs. Upstream `inq_template`

The [GitLab `inq_template`](https://gitlab.com/npneq/inq_template) is the upstream
starting point for standalone INQ projects. It differs from the local `inq-run`
workflow in several important ways.

### Template build process

```bash
# Clone the template
git clone https://gitlab.com/npneq/inq_template my_project
cd my_project

# Configure (downloads INQ and all deps automatically)
mkdir build install && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_BUILD_TYPE=Release
make -j4
make install

# Run example
cd runs
./nitrogen
```

The template `CMakeLists.txt` uses **FetchContent** to download INQ from GitLab:

```cmake
FetchContent_Declare(
    inq
    GIT_REPOSITORY https://gitlab.com/npneq/inq.git
    GIT_TAG master
)
FetchContent_MakeAvailable(inq)

add_subdirectory(runs)
```

### Side-by-side comparison

| Aspect | Local `inq-run` | GitLab `inq_template` |
|---|---|---|
| **INQ dependency** | `add_subdirectory()` against local pre-built copy | FetchContent — downloads INQ + all deps from GitLab |
| **First configure time** | ~1–5 min (deps already on disk) | Longer (full network download required) |
| **Subsequent builds** | Seconds (incremental, cache reused) | Same (fast once configured) |
| **GPU: NVIDIA CUDA** | CUDA 12.6.2, sm_80 (Ampere A30); `-DENABLE_CUDA=ON` | Any CUDA 11.5+; user supplies arch with `-DCMAKE_CUDA_ARCHITECTURES` |
| **GPU: AMD HIP/ROCm** | Not configured on this machine | Supported: `-DENABLE_HIP=1 -DCMAKE_HIP_ARCHITECTURES=gfx90a` |
| **Per-project isolation** | Yes — each `.cpp` dir gets its own `build/` | User manually writes and maintains `CMakeLists.txt` |
| **Ease of use** | One command: `inq-run` | Must understand CMake and write project files |
| **Portability** | Tied to `config.sh` paths on this machine | Fully portable; runs anywhere with network access |
| **Ideal use case** | Fast iterative R&D on this machine | Starting a new standalone project from scratch |
| **INQ version** | Whatever is in `inq/` (fixed, known working) | Tracks `master` (may break on new commits) |

### GPU support in the template

The template itself does not enable GPU by default; the user adds GPU flags at
configure time:

```bash
# NVIDIA (Turing, RTX / GTX 16xx):
cmake .. -DENABLE_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES=75 \
         -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc

# NVIDIA (Ampere, A100 / RTX 30xx):
cmake .. -DENABLE_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES=80

# AMD (Frontier, MI250X):
cmake .. -DENABLE_HIP=1 -DCMAKE_HIP_ARCHITECTURES=gfx90a \
         -DCMAKE_CXX_STANDARD=17

# Intel (Aurora):
cmake .. -DENABLE_SYCL=1
```

The local `inq-run` hardcodes sm_80 (A30 GPU on this machine) in `config.sh`.
To target a different architecture, edit `INQ_CUDA_ARCH` in `shared/config.sh`.

---

## 6. Important: CUDA ≥ 12.5 / CUB 2.4+ Fix

**Do not revert this change.**

**File:** `inq/external_libs/gpurun/include/gpu/reduce.hpp`

CUDA 12.5+ ships CUB 2.4+, which evaluates `invoke_result<TransformOp, T>` for
device lambdas at compile time. This fails without an explicit return type.

The fix applied:
- Wraps GPU kernels with `cuda::proclaim_return_type<Type>(kernel)`
- Changes `std::plus<>{}` to `std::plus<Type>{}` (deduction fails for device code)

Without this fix, INQ will not compile on CUDA 12.5+ (including CUDA 12.6.2 used
on this machine).

---

## 7. GPU Architecture Codes Reference

| Code | Architecture | Example GPUs |
|---|---|---|
| 61 | Pascal | GTX 1060, P100 |
| 70 | Volta | V100 |
| 72 | Volta (embedded) | Quadro GV100 |
| 75 | Turing | RTX 20xx, GTX 16xx |
| 80 | Ampere | A100, A30, RTX 30xx |
| 86 | Ampere (consumer) | RTX 3090 |
| 90 | Hopper | H100 |

This machine uses `INQ_CUDA_ARCH=80` (A30, Ampere). Change this in
`shared/config.sh` if targeting a different GPU.
