# Handover: CSD3 device setup (repo migration)

**Rolling file.** Latest milestone at top.
**Repo:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research` (branch `quantum-stopping-power`)
**Machine:** CSD3, login node `login-q-1` (RHEL8), GPU partition `ampere` (A100, sm_80)
**Plan:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/csd3-setup-cuda121-build.md`

---

## 2026-07-29 — engine build unblocked; four migration bugs fixed

### Context

Repo copied from the previous device (`/local/data/public/skcb2/tddft/`). `setup.sh`
failed; the transcript is `/home/skcb2/setup.log` (5499 lines). Task: diagnose and fix.

### Root causes found (all verified, none inferred)

**1. CUB fix incompatible with CUDA 12.1 — 105 compile errors. FIXED.**

`inq-local.patch` used `cuda::proclaim_return_type<Type>(kernel)` in
`gpu/reduce.hpp` — documented in `setup.sh:11` as "the CUB fix for CUDA 12.5+".
The previous device had CUDA 12.6.2.

**CSD3's newest CUDA toolkit is 12.1.** Verified exhaustively: `/usr/local/software/cuda/*`
stops at 12.1; `module avail cuda` under the rhel7 default env AND `rhel8/default-amp`
both stop at `cuda/12.1`; `nvhpc` is 21.x (CUDA 11 era). Nothing >= 12.4 exists.

CUDA 12.1's libcu++ `__detail::__return_type_wrapper`
(`/usr/local/software/cuda/12.1/include/cuda/functional:73-78`) has an
**unconstrained greedy template ctor and no copy ctor**. When thrust's
`transform_input_iterator_t` copy-constructs it, that ctor outranks the implicit
copy ctor and tries to build the *lambda* from a *wrapper*:
`cuda/functional(77): error: no instance of constructor "lambda [](auto)->auto::<unnamed>"`.
libcu++ constrained it in a later release — hence 12.6 OK, 12.1 broken.

**Fix:** replaced `cuda::proclaim_return_type` with a local `fixed_return` functor
(plain struct: implicit copy ctor + fixed return type). This **preserves** the CUB
fix's intent (pinning the reduction value type) on 12.5+ while also working on 12.1,
so `CLAUDE.md`'s "do not revert the CUB fix" constraint holds.

Applied to (user approved editing `inq/`, overriding `.claude/rules/inq-immutable.md`
for this change):
- `inq-local.patch` — **regenerated with `git -C inq diff`**, not hand-edited (a
  hand-edited version had corrupt hunk offsets). Verified: `patch -p1 --dry-run`
  applies cleanly to a pristine `44f73d9527ab` export.
- `inq/external_libs/gpurun/include/gpu/reduce.hpp`
- `inq-study/external_libs/gpurun/include/gpu/reduce.hpp` (byte-identical copy;
  **submodule change is UNCOMMITTED — needs its own commit**)

**Verified** by reproducer (`<scratchpad>/repro.cu`) with INQ's own flags from
`inq/CMakeLists.txt:47` against nvcc 12.1.105:
`proclaim_return_type` → the identical production error; `fixed_return` → **compiles**.

**2. OOM kill on the login node — 16x `cicc died due to signal 9`. FIXED.**

`setup.sh:56` ran `cmake --build --parallel` with no job count → one nvcc per core
(76 cores), while the per-user login cgroup limit is **20 GB**
(`/sys/fs/cgroup/memory/user.slice/user-<uid>.slice/memory.limit_in_bytes`). Each
`cicc` on INQ's template-heavy TUs takes several GB. `login-q-1` also has **no GPU**.

**Fix:** bounded parallelism everywhere —
- `setup.sh` → `--parallel "${BUILD_JOBS:-4}"`
- `shared/bin/inq-run` → `-j"${INQ_BUILD_JOBS:-4}"` (was `-j$(nproc)`)
- `shared/config.sh` → exports `INQ_BUILD_JOBS` (default 4)
- new `shared/bin/build-inq.slurm` → builds on an `ampere` node instead

**3. `shared/config.sh` still held previous-device absolute paths. FIXED.**

`INQ_SOURCE=/local/data/public/skcb2/tddft/inq`,
`INQ_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/nvcc`,
`INQ_PYTHON_EXE=/local/data/public/.../venv/bin/python3` — none exist here, so
`inq-run` could not have worked. Rewritten to derive from the repo root
(`shared/` is one level below it) with env overrides preserved.

**4. `.claude/settings.json` PreToolUse hook pointed at the old device. FIXED.**

`python3 /local/data/public/skcb2/tddft/.claude/hooks/commit_message_check.py` →
every Bash call aborted. Repointed to `$CLAUDE_PROJECT_DIR/.claude/hooks/...`.

### Build status — IN PROGRESS

Two failed attempts, then a working split. **Compute nodes have no outbound
network**, so INQ's `FetchContent` cannot clone Catch2/pybind11 from github.com:
```
Failed to clone repository: 'https://github.com/catchorg/Catch2.git'
FetchContent.cmake:1622 ... external_libs/gpurun/CMakeLists.txt:18
```
Job 32352533 died on exactly this, because its `CLEAN=1` step wiped the build dir
(and with it the previously-fetched deps) before configuring on the compute node.
**Cost: the 614 pre-existing object files (433 MB) from the original device were
discarded.** They were rhel7-stack objects and not reusable against
`rhel8/default-amp` anyway, but the wipe was avoidable.

**Working split — configure on login node (has network), compile on GPU node:**
1. `shared/bin/configure-inq.sh [inq|inq-study]` — login node. **Succeeded**
   (Catch2 v3p3p2 + pybind11 v2.12.0 downloaded, "Generating done").
2. `sbatch shared/bin/build-inq.slurm` — ampere node, compile only; it now hard-errors
   with a pointer to step 1 if `CMakeCache.txt` is missing.

| Job | Result |
|---|---|
| 32352533 | FAILED — configure on compute node, no network (see above) |
| 32352675 | FAILED at 92% — CUDA 11.4/12.1 header clash (cause 5 below) |
| 32353243 | FAILED at link — gcc 8.5 `std::filesystem` (cause 6 below) |
| 32353761 | **RUNNING at handover** — clean toolchain, log `inq-build-32353761.out` |

**5. CUDA 11.4 and 12.1 headers mixed on CPATH. FIXED.**

`module load rhel8/default-amp` pulls `cuda/11.4` in as a **locked dependency**
(it cannot be unloaded). That leaves both toolkits visible:
`CPATH = .../cuda/12.1/include : ... : .../cuda-11.4.../include`, `CUDA_HOME` = 11.4.
nvcc 12.1 implicitly includes its own `crt/sm_80_rt.hpp` *and* picks up 11.4's off
CPATH, so every sm_80 TU dies with:
```
crt/sm_80_rt.hpp(141): error: more than one instance of overloaded function
"__nv_associate_access_property_impl" has "C" linkage
```
Killed job 32352675 at 92% — `inq_executable`, `_pinq`, all tests, all benchmarks
(NOT just test targets).

**Verified** by re-running the exact failing nvcc command for
`inq/src/main/unit_tests_main.cpp` on the login node: fails with 11.4 on CPATH,
compiles cleanly without it.

**Fix:** `shared/bin/csd3-env.sh` strips `cuda-11.4|cuda/11.4` from `CPATH` and
`LD_LIBRARY_PATH` and pins `CUDA_HOME`/`CUDA_PATH` to 12.1, with a guard that
hard-fails if 11.4 reappears.

**6. Host compiler was gcc 8.5 — `std::filesystem` link failures. FIXED.**

`rhel8/default-amp` loads **no compiler module**, so the host compiler defaulted to
system gcc **8.5.0**, where `std::filesystem` lives in a separate `-lstdc++fs` that
INQ's CMake does not link. 6 targets failed with dozens of
`undefined reference to std::filesystem::create_directories(...)` /
`path::_M_split_cmpts()`. The original device used gcc 9.3.0, where `<filesystem>`
is in libstdc++ proper — so this never surfaced before.

**Verified** 2026-07-29 with a `std::filesystem` probe under nvcc:
gcc 8.5.0 → 4 undefined refs; gcc 9.4.0 → 0, links OK.

**Fix:** `csd3-env.sh` loads `gcc/9.4.0/gcc-11.2.0-tfj3hud` — chosen to match the
gcc-9.4.0 ABI of the openmpi 4.1.1 / ucx stack `rhel8/default-amp` already loads
(zen2 spack view, same as those dependencies). Changing the host compiler requires
a **CLEAN reconfigure** (CMake caches `CMAKE_CXX_COMPILER`); done via
`CLEAN=1 ./shared/bin/configure-inq.sh inq`.

Final toolchain, now consistent end to end:
`gcc 9.4.0 + CUDA 12.1 (nvcc 12.1.105) + openmpi 4.1.1`, sm_80.

**CORRECTION to an earlier note in this file:** CSD3 *does* have a newer CUDA —
`cuda/12.8.1/gcc/kdeps6ab`, visible only after `module load rhel8/ampere-env/2025-06-01`
(and it has a clean CPATH out of the box). It is NOT used because that environment
ships no openmpi INQ can link against. Revisit if MPI appears there. The
`fixed_return` fix is version-agnostic, so it works on either toolkit.

Module stack is now `rhel8/default-amp` + `cuda/12.1` for configure, build, and runs
(the original failed build mixed the legacy **rhel7** stack — gcc 9.3.0 +
`Cluster-Apps/openmpi/4.0.4` — on a RHEL8 machine). Note `rhel8/default-amp` supplies
openmpi 4.1.1/gcc-9.4.0 and host gcc is `/usr/bin/gcc` 8.5.0.

**`sbatch`/`squeue` only work from the rhel8 module env** — the default login env has
a rhel7 SLURM whose plugin stack fails on missing `liblua-5.1.so`. Wrap as:
`bash -c '. /etc/profile.d/modules.sh; module purge; module load rhel8/default-amp; sbatch ...'`

SLURM: account `mphil-nikiforakis-skcb2-sl2-gpu`, partition `ampere`.

### NOT done / next steps

1. **Confirm job 32352675 built cleanly** — `tail inq-build-32352675.out`, expect
   "BUILD OK". Nothing is validated until this passes.
2. **Python environment is broken.** `venv/bin/python3` symlinks to `/usr/bin/python3`
   = **3.6.8**; `inq-stack/pyproject.toml` requires **>=3.10**; `inqview` is NOT
   installed (setup.sh aborted at step 3, so step [4/4] `pip install -e inq-stack/`
   never ran). `/usr/bin/python3.11` and `3.12` exist. **Deliberately not touched yet**:
   INQ's configure cached `PythonInterp` = `<repo>/venv/bin/python` 3.6.8 against
   `/usr/lib64/libpython3.6m.so` for the `_pinq` bindings, so swapping the venv
   mid-build would break that link. Do it after the build completes, then decide
   whether `_pinq` is needed (inqview itself is pure Python and does not need it).
3. **`inq-study/` submodule change is uncommitted** — needs its own commit.
4. **Commit the repo changes** (`inq-local.patch`, `setup.sh`, `shared/config.sh`,
   `shared/bin/*`, docs). Not committed yet.
5. **User-requested validation (2026-07-29, mid-session): wavepacket-in-vacuum
   Gaussian-broadening check** across several energies and starting sigma —
   matrix drafted below, **awaiting user approval before launching** (per
   `.claude/rules/validation-gates.md`, expensive runs are the user's call).

### Engine build: DONE (job 32355421, 2026-07-29 18:22)

`inq` built and installed, prefix
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/inq/install`.
**Linkage verified, not just exit status:**
`objdump -p inq/install/bin/inq | grep NEEDED` →
`libcufft.so.11`, `libcudart.so.12`, `libcublas.so.12`, `libcublasLt.so.12`,
`libcusolver.so.11` (all CUDA 12.x), and `ldd` reports no unresolved libraries.
`inq-study` is CONFIGURED against the same 12.1 toolkit (`CUDA_cufft_LIBRARY` =
`/usr/local/software/cuda/12.1/lib64/libcufft.so`); `inq-run` compiles it from
source via `add_subdirectory`, so a full build of it is not needed.

**9. A "successful" build that was silently mis-linked to CUDA 11.4. FIXED.**

**Job 32353761 reported `BUILD OK` with 0 compile and 0 link errors but was
UNUSABLE.** The first sweep exposed it at runtime:
```
./run: error while loading shared libraries: libcufft.so.10
```
`libcufft.so.10` is CUDA **11.x**'s cuFFT soname (12.x ships `libcufft.so.11`).
`objdump -p` showed BOTH `inq/install/bin/inq` and the run binary NEEDing
`libcufft.so.10` / `libcudart.so.11.0` / `libcublas.so.11`.

**Cause:** the cause-5 scrub was INCOMPLETE. It stripped `CPATH` and
`LD_LIBRARY_PATH` but missed:
- **`LIBRARY_PATH`** (4 stale entries) — what gcc/ld searches at LINK time;
- **`CMAKE_PREFIX_PATH`** (1 stale entry) — what `find_package(CUDAToolkit)` searches.

So headers resolved to 12.1 while libraries resolved to 11.4 — a combination that
compiles and links **silently** and only fails on execution.

**Fix:** `csd3-env.sh` now scrubs `CPATH`, `LD_LIBRARY_PATH`, `LIBRARY_PATH`,
`CMAKE_PREFIX_PATH`, `PKG_CONFIG_PATH`, `LD_RUN_PATH`; prepends 12.1's `lib64`;
pins `CUDAToolkit_ROOT`; and the guard now checks all four critical variables.
Both engines were CLEAN-reconfigured and `inq` fully rebuilt afterwards.

**LESSON (do not lose this):** on this cluster a green build is NOT evidence of a
working build. Always verify
`objdump -p <binary> | grep NEEDED | grep -iE 'cufft|cudart|cublas'`
shows 12.x sonames before declaring an engine usable.

**7. Compute nodes have no network — `inq-run` re-triggered FetchContent. FIXED.**

`inq-run` creates a FRESH build dir per run.cpp and `add_subdirectory(INQ_SOURCE)`,
which re-runs FetchContent and dies on a compute node with
`Failed to clone repository: .../Catch2.git`. Fixed generically: `shared/config.sh`
exports `INQ_DEPS_CACHE="$INQ_SOURCE/build/_deps"` and `inq-run` turns every
`<name>-src` there into `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<path>`.
**Verified:** the run's build tree then contains only `*-build` dirs and no local
`*-src`, and `CMakeCache.txt` shows
`Catch2_SOURCE_DIR=.../inq/build/_deps/catch2-src` — reused, not cloned.

**8. `wp_traversal_energy/run.cpp` REQUIRES inq-study, even with `WP_ETA=0`.**

Corrects an earlier assumption in this file ("WP_ETA=0 needs only stock inq" —
WRONG). The run ALWAYS constructs `perturbations::absorbing cap(ETA*1.0_Ha, ...)`
and passes it to `real_time::propagate`, so the CAP path is **instantiated at
compile time** regardless of ETA's runtime value. Stock `inq` keeps the scalar
potential real (`field<real_space, double>`), so the CAP's
`vk[ix][iy][iz] += complex(0.0, ...)` (`src/perturbations/absorbing.hpp:45`) fails
to compile. `absorbing.hpp` is byte-identical in both engines — the difference is
the **complexified scalar potential** in inq-study's `self_consistency.hpp`,
`ks_hamiltonian.hpp`, `ks_hamiltonian`/`propagate`/`electrons` (8 files differ).
`shared/bin/run-dispersion.slurm` therefore sets `INQ_SOURCE=$REPO_ROOT/inq-study`
while keeping the share paths on `inq/install/share`.

### Validation status

| Check | Status |
|---|---|
| `fixed_return` compiles under CUDA 12.1 + INQ's flags | **PASS** (reproducer) |
| `inq-local.patch` applies cleanly to pristine `44f73d95` | **PASS** (`patch --dry-run`) |
| Full `inq` build on ampere | **PASS** (job 32355421) |
| Binaries link CUDA 12.x, not 11.4 | **PASS** (`objdump -p`, `ldd`) |
| WP-in-vacuum Gaussian broadening vs analytic | **PASS** (job 32355881, 9/9) |
| `ctest` on A100 | not started |

### Physics validation: free-Gaussian dispersion — PASS (job 32355881, 43 min)

9/9 runs, sigma_WP in {1,2,3} Bohr x E in {1,10,100} eV, `WP_ETA=0`, built against
`inq-study`. Produced by `scripts/wp_traversal_energy/dispatch_dispersion.py`,
checked by `analyse_dispersion.py`; outputs under
`scripts/wp_traversal_energy/results/dispersion/<run>/`.

| Quantity | Result |
|---|---|
| Worst deviation of sigma_d(t) from `sqrt(sigma0^2/2 + t^2/(2 sigma0^2))` | **0.002 %** |
| Worst k0-dependence of broadening (Galilean invariance) | **0.002 %** |
| Transverse `sigma_x2` vs the same law | **0.000 %** |
| Energy conservation | constant to ~12 significant figures |
| Norm conservation | `norm_check` = 0.999999997 |

**Anti-circularity checks — the agreement is real, not a formula echoing itself:**
- Centroid motion is an INDEPENDENT observable, unused in the width comparison:
  `disp_sig3_E100` moved z_mean −21.2868 → +20.8702 Bohr over t = 15.55 a.u.,
  v = 2.71106 vs k0 = 2.71106 → **0.000 %**. The packet really traversed 42 Bohr.
- `sigma_z2` (17.9341818603) and `sigma_x2` (17.9334678607) DIFFER in the 6th digit —
  the expected numerical asymmetry between moving and stationary axes. An analytic
  shortcut would return identical values.
- `sigma_x2(0)` = 4.5 = sigma0^2/2 for sigma_WP = 3, confirming the sigma convention.

**What this DOES validate:** grid/cell setup, FFT, kinetic operator, ETRS propagator,
wavepacket injection, and — most relevant — the `gpu::reduce` kernels our
`fixed_return` patch rewrote, since `<z>`/`<z^2>` are computed through
`operations::sum`/`integral` built on them.

**What this does NOT validate:** Hartree, XC, pseudopotentials, interacting dynamics,
or the CAP. Near-exact agreement is EXPECTED here — for a kinetic-only Hamiltonian the
FFT propagator advances each Fourier mode by an exact phase `exp(-i k^2 t/2)`, so a
Gaussian evolves essentially exactly. Strong machinery check, NOT evidence the
interacting engine is correct. `ctest` on an A100 remains the outstanding broad check.

### WP-in-vacuum validation — design (AS EXECUTED, job 32355881)

Superseded-by-results note: this section was written as a *proposal for approval*;
it is retained as the **design rationale + analytic reference** for the run that
has now completed. Outcome is the PASS section immediately above. Kept because the
analytic law and the box-sizing constraints are what any future re-run must honour.

Run script: `ResearchProject/systems/vacuum/scripts/wp_traversal_energy/run.cpp`.
Set `WP_ETA=0` → no CAP → pure periodic vacuum, non-interacting, single electron.

CORRECTION (was wrong in the original proposal): this does **NOT** work against
stock `inq`. `WP_ETA=0` disables the absorber *numerically*, but
`perturbations::absorbing` is still instantiated at COMPILE time, so the
translation unit needs the CAP-enabled engine. Build against **`inq-study`**
(`INQ_SOURCE=$REPO_ROOT/inq-study`, as `shared/bin/run-dispersion.slurm` sets).
Building against stock `inq` fails — that was job 32354396.

Analytic reference, stated in that file's header (lines 33-37) and consistent with
the standard free-Gaussian result — for wavefunction amplitude
`psi ~ exp(-(z-z0)^2 / (2 sigma0^2))` the **density** width is

```
sigma_dens(t) = sqrt( sigma0^2/2 + t^2/(2 sigma0^2) )      [atomic units, m_e = 1]
R(t) = sigma_dens(t)/sigma_dens(0) = sqrt(1 + (t/sigma0^2)^2),  spreading time tau = sigma0^2
```

(Consistency check: `Delta_x(0) = sigma0/sqrt2`, `Delta_p = 1/(sqrt2 sigma0)`,
`Delta_x(t)^2 = Delta_x(0)^2 + (Delta_p t/m)^2` reproduces the formula exactly.)

Measured from `results/<out>/raw/observables/wp_real_space_stats.csv`, column
`sigma_z2` (= `<z^2> - <z>^2`), so `sigma_dens = sqrt(sigma_z2)`.

Three independent things this tests:
1. `sigma_dens(t)` tracks the analytic curve for each `sigma0`.
2. **Broadening is independent of k0** (Galilean invariance) — the sharpest check,
   and free across an energy sweep.
3. Centroid `z_mean(t) = z0 + k0 t`, and `E_total` conserved (= kinetic only).

Per `.claude/rules/sigma-wp-convention.md`, runs are labelled by the **wavepacket**
sigma (`WP_SIGMA`); the density width is `sigma_WP/sqrt2` at t=0 and is a derived
quantity, not the label.

Draft matrix (`WP_SIGMA` x energy), `WP_ETA=0`, `WP_DT=0.01`, `WP_WF_EVERY=5`:

| sigma_WP (Bohr) | E (eV) | k0 (a.u.) | notes |
|---|---|---|---|
| 1.0, 2.0, 3.0 | 1, 10, 100 | 0.271, 0.857, 2.711 | 9 runs, full grid |

Box/steps must be sized per (sigma0, k0) so the packet neither wraps the periodic
box nor lets its tails reach the edge (`sigma_z2` is meaningless once density wraps):
requires `launch_z - 5 sigma_dens(0) > -LZ/2` and
`launch_z + k0 T + 5 sigma_dens(T) < LZ/2`. That dispatcher was written:
`ResearchProject/systems/vacuum/scripts/wp_traversal_energy/dispatch_dispersion.py`
(computes `LZ`, `LPERP`, `NSTEPS`, `launch_z` from `(sigma0, k0)`; aborts if the
k-space margin `(pi/h)/(k0 + 4 dk) <= 1`).

Tiering — **Tier 2 was executed** (3 sigma x 3 E = 9 runs, the k0-independence
test). Tier 1 was skipped as redundant once Tier 2 was approved. **Tier 3 (dt and
grid-spacing `h` refinement at one grid point, to separate discretisation error
from physics) was NOT run and remains outstanding** — note the worst deviation was
already 0.002%, so Tier 3 would be tightening an already-tight bound.

---

## 2026-07-30 — measured SLURM throughput; Python env rebuilt (3.6 -> 3.11)

### Measured GPU throughput on `ampere` (all from `sacct`, not estimated)

Hardware: 1x NVIDIA A100-SXM4-80GB, node `gpu-q-21`, `AllocTRES =
billing=32,cpu=32,gres/gpu=1,mem=250G,node=1`. QOS `gpu1`.

Queue wait = `Start - Submit`, from `sacct`:

| Job | What | Submit | Wait | Elapsed | State |
|---|---|---|---|---|---|
| 32352533 | inq-build | 17:16:07 | 26 s | 00:01:00 | FAILED (FetchContent, no net on node) |
| 32352675 | inq-build | 17:20:48 | 2 m 04 s | 00:01:43 | FAILED (CUDA 11.4/12.1 header clash) |
| 32353243 | inq-build | 17:31:28 | 14 s | 00:05:41 | FAILED (gcc 8.5 std::filesystem) |
| 32353761 | inq-build | 17:42:26 | 1 m 15 s | 00:07:03 | "COMPLETED" but MIS-LINKED to CUDA 11.4 |
| 32354260 | wp-dispersion | 17:53:08 | 29 s | 00:00:44 | FAILED (re-triggered FetchContent) |
| 32354396 | wp-dispersion | 17:57:29 | 19 s | 00:02:30 | FAILED (needs inq-study, not inq) |
| 32355421 | inq-build | 18:15:01 | **5 s** | 00:07:00 | **COMPLETED, linkage verified 12.x** |
| 32355881 | wp-dispersion | 18:23:18 | **6 s** | **00:43:43** | **COMPLETED, 9/9 runs** |

Queue wait: min 5 s, median ~27 s, max 2 m 04 s, mean ~37 s (n=8).
**Caveat:** all 8 were submitted in one ~70-minute window on 2026-07-29 evening.
Checked again 2026-07-30: `ampere` had **602 PENDING / 157 RUNNING** with ~18 nodes
in `maint`. These waits are NOT a guaranteed steady state — re-measure before
assuming fast turnaround.

Engine build cost: **7 min** on the GPU node (configure on login node separately,
~1 min). Rebuild is incremental unless `CLEAN=1`.

### Per-run propagation cost (job 32355881, parsed from step timestamps)

Total job 43:43, of which **37.1 min propagation** and ~6.6 min build + startup.
`s/step` is mean; `med` is median (mean is inflated by periodic I/O steps).

| Run | sigma_WP | E (eV) | grid (Bohr) | Mpts | steps | prop time | s/step | med s/step |
|---|---|---|---|---|---|---|---|---|
| disp_sig1_E1   | 1 | 1   | 15x15x15 | 0.05 | 174  | 3.6 s   | 0.021 | 0.011 |
| disp_sig1_E10  | 1 | 10  | 15x15x20 | 0.07 | 174  | 4.3 s   | 0.025 | 0.013 |
| disp_sig1_E100 | 1 | 100 | 15x15x20 | 0.07 | 174  | 4.3 s   | 0.025 | 0.013 |
| disp_sig2_E1   | 2 | 1   | 30x30x35 | 0.49 | 693  | 74.7 s  | 0.108 | 0.059 |
| disp_sig2_E10  | 2 | 10  | 30x30x35 | 0.49 | 693  | 74.0 s  | 0.107 | 0.059 |
| disp_sig2_E100 | 2 | 100 | 30x30x50 | 0.70 | 693  | 99.5 s  | 0.144 | 0.082 |
| disp_sig3_E1   | 3 | 1   | 45x45x50 | 1.58 | 1559 | 506.7 s | 0.325 | 0.213 |
| disp_sig3_E10  | 3 | 10  | 45x45x60 | 1.90 | 1559 | 598.4 s | 0.384 | 0.256 |
| disp_sig3_E100 | 3 | 100 | 45x45x85 | 2.69 | 1559 | 859.0 s | 0.551 | 0.365 |

All at `dt=0.01` a.u., spacing `h=0.4` Bohr, 1 electron, non-interacting, `WP_ETA=0`,
`wf_every=5`, `mom_every=5`.

**Scaling rule of thumb (single non-interacting orbital, this grid range):**
cost/step is roughly linear in grid points — ~0.13 s per Mpt (median, from the
sig3 rows) down to ~0.2 s/Mpt at the smallest grids where launch overhead
dominates. A 2.7 Mpt box costs ~0.37 s/step. **Do NOT extrapolate to interacting
runs** — Hartree/XC + many orbitals change the constant AND the scaling; the
9 runs here carry exactly one occupied state and no SCF.

### Budget and priority tiers (this is the real throughput constraint)

`mybalance` + `sacctmgr show qos`, 2026-07-30:

| Account | Partition QOS | QOS priority | Max wall | Max GPU/user | Hours available |
|---|---|---|---|---|---|
| `mphil-nikiforakis-skcb2-sl2-gpu` | `gpu1` | 5000 | 1-12:00:00 | 64 | **364** (of 365, 1 used) |
| `nikiforakis-sl3-gpu` | `gpu2` | 1000 | 12:00:00 | 32 | ~3,000 |
| `peng-sl3-gpu` | `gpu2` | 1000 | 12:00:00 | 32 | ~3,000 |
| `blakely-sl3-gpu` | `gpu2` | 1000 | 12:00:00 | 32 | ~2,998 |
| `*-sl4-gpu` | `gpu3` | **0** | 12:00:00 | (none set) | scavenger tier |

All 8 jobs above ran on the **SL2** account (`gpu1`, priority 5000) — that is why
waits were seconds. The SL2 GPU pot is only **365 hours**, so it should be spent on
work that needs fast turnaround. Bulk/overnight sweeps belong on an SL3 account
(`gpu2`, priority 1000, ~9,000 h total across three) accepting longer queueing and a
**12 h wall cap** (vs 36 h on `gpu1`) — the 12 h cap matters, hence the
final-timestep-checkpoint rule.

Not verified: actual observed wait on `gpu2`/`gpu3` — no job has been submitted on
those accounts from this device. Priority numbers are the QOS config, not a measured
wait.

### Python environment rebuilt — DONE

Was: `venv/bin/python` = **3.6.8**, `inq-stack/pyproject.toml` requires `>=3.10`,
`inqview` NOT installed (setup.sh died at step 3, so step [4/4] never ran).

Done 2026-07-30:
- Old venv preserved as `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/venv.py36.bak`
  (11 MB) — delete once confident.
- New venv from **`/usr/bin/python3.11`** (3.11.13). Note: the module tree tops out
  at `python/3.8.11` (gcc) and `python/3.11.0-icl` (Intel build); the system RPM
  python3.11/3.12 are newer and need no module, so no module is loaded for Python.
- `pip install -e "inq-stack/[analysis,viz,test]"` → numpy 2.4.6, scipy 1.17.1,
  pandas 3.0.5, matplotlib 3.11.1, vtk 9.6.2, imageio 2.37.4, pytest 9.1.1.
  Login node has PyPI access (verified HTTP 200).
- **Safe to rebuild:** `ENABLE_PYTHON:BOOL=OFF` in both `inq/build/CMakeCache.txt`
  and `inq-study/build/CMakeCache.txt` — the `_pinq` bindings are NOT built, so the
  stale cached `PYTHON_EXECUTABLE=.../venv/bin/python` and
  `PYBIND11_PYTHON_VERSION=3.6.8` are inert. Earlier caution about disturbing them
  was unfounded.

Test suite: `cd inq-stack && ../venv/bin/python -m pytest -q` →
**163 passed, 8 skipped, 5 xfailed, 1 xpassed**. Includes `test_deps_clean.py`
(ADR 0003 layering invariant). The 8 skips are `test_gaussian_psp.py`
("electron-ONCV template not present") — a missing fixture, pre-existing.

### NumPy 2 migration (consequence of the venv rebuild)

numpy 2.4.6 REMOVED two things this repo used. Both were **latent** — masked by
`xfail`/`skip` markers, so a green suite did not prove them working.

1. `float(np.array([x]))` on a 1-element array → `TypeError` (deprecated 1.25,
   removed 2.0). One site: `test_lindhard.py::test_plasmon_zero_q_limit`. The
   library was correct (`plasmon_omega` is documented `-> NDArray`); the TEST was
   wrong. Fixed by indexing `[0]`.
2. **`np.trapz` REMOVED** (renamed `np.trapezoid`). Verified:
   `hasattr(np,'trapz') == False` on 2.4.6. 12 call sites across 8 files:
   - `inqview/pipeline/lindhard.py` (3) → **`scipy.integrate.trapezoid`**, not
     `np.trapezoid`, because `pyproject.toml` declares `numpy>=1.24` where
     `np.trapezoid` does not exist yet; `scipy>=1.10` is pinned so
     `integrate.trapezoid` is valid on every supported numpy. Module already
     imported `scipy.integrate`.
   - 2 test files + 5 run-analysis scripts (`ResearchProject/systems/jellium/...`,
     `.../localised_jellium/hypotheses/campaign_autorun_study/...`,
     `docs/campaigns/ml-patterns/kernels/formfactor.py`) → `np.trapezoid`
     (these run only in this venv).

   The `lindhard.py` sites sit behind xfailing tests, so they were never executed.
   Verified fixed by calling `L.stopping_power` directly: v=1/3/10 →
   S=0.131/0.0364/0.00327 a.u., path executes.

   Side effect: `test_f_sum_rule[0.4]` now **XPASSes**. Its `xfail` had been
   absorbing an `AttributeError` from `np.trapz`, i.e. the documented "dynamical
   Lindhard high-omega sign error" was not actually being exercised at all under
   numpy 2. The other 3 parametrisations still xfail. **The Lindhard sign issue is
   unchanged and still open** — this only means the xfails now test physics rather
   than dying on a missing attribute.

Exhaustive re-sweep for other numpy-2 removals (`np.float_`, `np.NaN`, `np.Inf`,
`np.alltrue`, `np.product`, `np.row_stack`, ...) across all `*.py`/`*.ipynb`
excluding `venv/`: **clean**.

### Still outstanding

- **`ctest` on an A100** — the broad engine check. Not started. Needs user approval
  (expensive). `INQ_EXEC_ENV="mpirun.openmpi -np 4" ctest` for the MPI pass.
- **NOTHING IS COMMITTED.** Working tree carries all of the above. `inq-study/`
  submodule change (`external_libs/gpurun/include/gpu/reduce.hpp`) needs its OWN
  commit inside the submodule first. Per `.claude/rules/commit-messages.md`, split
  production vs infra into separate commits.
- Tier 3 dt/h refinement (see design section above) — not run.
- Old `rerun_*.sh` in `scripts/wp_traversal_energy/` still hold old-device absolute
  paths. Not fixed; `dispatch_dispersion.py` supersedes them.
- `venv.py36.bak` can be deleted.
- Plotting of the sigma(t) curves is now UNBLOCKED (matplotlib present) but not done;
  `analyse_dispersion.py` is still pure-stdlib and does not plot.

---

## 2026-07-30 — PreToolUse hook was silently disabled; 3 further bugs found

### The reported error: "PreToolUse:Bash hook failed with non-blocking status code"

**Fault #10.** `.claude/settings.json` invokes the hook as bare `python3`. On CSD3/RHEL8
that is `/usr/bin/python3` = **3.6.8**, and `commit_message_check.py` had
`from __future__ import annotations` — a hard **SyntaxError** on 3.6 ("future feature
annotations is not defined") → exit **1**.

PreToolUse exit-code semantics: `0` allow, `2` block (stderr fed to the model),
**anything else = non-blocking error** (warning shown, tool still runs). So exit 1 meant
the commit-message guard was **completely non-functional** — not merely noisy. Every
commit went unchecked. The earlier fix in this handover repaired the hook's PATH but not
its interpreter, so it never actually ran on this device.

Not a venv issue: hooks use PATH `python3`, never `venv/bin/python`.

**Fix:** removed the future import from both `.claude/hooks/*.py` — every annotation in
them is a plain name (`str`/`bool`/`Result`), so it bought nothing. Both now carry a
comment forbidding PEP 585 (`list[str]`) / PEP 604 (`X | Y`) annotations, which would
reintroduce the incompatibility. `settings.json` keeps bare `python3` (portable) rather
than hardcoding an interpreter version.

### Fault #11 — file_placement_check.py had a stale old-device REPO

`REPO = "/local/data/public/skcb2/tddft"`, hardcoded. After migration nothing matched, so
`_to_rel()` returned `None` for every in-repo file: a **silent always-allow no-op**. Now
derived from `CLAUDE_PROJECT_DIR`, else three dirnames up from the file. (Not registered
as a hook — only imported by its eval — so it was not producing the reported error, but
it was equally broken.)

### Fault #12 — settings.json had no env block; ~/.bashrc sets nothing

`build-run-env.md` (LOCKED 2026-06-11) specifies the two additive vars pinned in
`settings.json` `env`, PATH deliberately left profile-driven (pinning PATH there would
shadow system git/python3/nvcc). The `env` block was **absent entirely**, and `~/.bashrc`
sets neither the share vars nor `shared/bin` on PATH — the prerequisites CLAUDE.md
documents were never applied on this device. The vacuum sweep worked only because
`run-dispersion.slurm` exports them explicitly.

Added (dirs verified to exist): `INQ_SHARE_PATH` and `PSEUDOPOD_SHARE_PATH` pointing at
`<repo>/inq/install/share[/pseudopod]`. **Takes effect on Claude Code restart**, not in
the session that wrote it.

`run_build_run_env_eval.py` asserted the previous device's literal paths; it now DERIVES
them from its own location. Deliberate asymmetry: `settings.json` keeps an absolute path
(machine config), the eval derives — so a future repo move fails this eval LOUDLY instead
of leaving a stale value unnoticed.

**STILL OUTSTANDING (user's shell, not edited):** `~/.bashrc` has no `shared/bin` on PATH,
so bare `inq-run` is unavailable interactively. Not changed unilaterally — user's profile.
Line to add:
`export PATH="/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/shared/bin:$PATH"`

### Fault #13 — REAL BUG: forbidden word at end-of-message was ALLOWED

Found by spot-checking the block path after the interpreter fix. `_first_forbidden`
exempts path/identifier context by inspecting the neighbouring character. When the match
sits at the very start/end there is no neighbour and the value is `""` — and **the empty
string is a substring of every string in Python**, so `"" in "-/_"` is True and the match
was exempted. `chore(repo): made by Claude` was **silently ALLOWED**.

Fixed with `(before and before in ...) or (after and after in ...)`: absence of a
neighbour is prose context, not path context.

**Why the eval missed it at 22/22:** all 9 REJECT cases placed the forbidden word before
another character (`add Claude trailer`, `Anthropic SDK`, `Claude <x@y.z>`, `Claude Code`).
None ended with it. Green proved nothing about the boundary — the same "green != working"
pattern as the CUDA 11.4 mis-link and the numpy-2 sites hidden behind xfails. Three for
three this session.

### Fault #14 — FALSE POSITIVE: doc heredocs read as commit messages

Discovered when the now-working hook **blocked this very handover write**.
`_extract_commit_message` treated any heredoc in a command merely CONTAINING the substring
`git commit` as the commit message. So
`cat >> handover.md <<'EOF' ... every git commit was unchecked ... EOF`
had its entire body checked and blocked on the first forbidden word. In a repo that
documents the commit rule at length this fires constantly.

Fixed: extraction now anchors on the invocation (`\bgit\s+commit\b`), considers only text
AFTER it, and requires `-F`/`--file` on that invocation's own first line for the heredoc
branch — which every real heredoc commit has (`git commit -F - <<'EOF'`). Verified real
commits still extract (`-m`, multi-`-m`, `git   commit` spacing, `-F -` heredoc) while doc
heredocs return `None`. This handover appending successfully IS the regression test.

### Verification — all 5 programmatic evals PASS under bare python3 (3.6.8)

All five `run_*.py` had the same future-import defect (their own docstrings say to run them
with `python3`), so none could execute. Removed from all five after confirming no PEP
585/604 annotations were present.

| Eval | Result |
|---|---|
| `run_build_run_env_eval.py` | PASS 3/3 (was FAIL: stale paths + missing env block) |
| `run_cluster_o_drift_eval.py` | PASS (19 canonical observables covered) |
| `run_cluster_r_eval.py` | PASS (figure standard) |
| `run_commit_hook_eval.py` | **PASS 32/32** (was 22/22 with two coverage gaps) |
| `run_file_placement_eval.py` | PASS 15/15 |

Commit-hook eval grew 22 -> 32: +3 REJECT (forbidden word last token of subject/body),
+2 ACCEPT (path context that also ends the message stays exempt — proves no over-block),
+1 defensive index-0 sentinel, +4 EXTRACT (spacing, `--amend`, two doc-heredoc
false-positive guards). Regression-proved: reverting the boundary guard in a scratch copy
makes the two new end-of-message cases return `ok=True`.

The index-0 case is labelled a DEFENSIVE sentinel, not a discriminating test: index 0 is
the action-word slot, so the subject-format rule rejects it either way (confirmed the
pre-fix code also rejected it). The `before` guard is unreachable-by-verdict today.

### Hook behaviour matrix (verified on 3.6.8)

| Input | Exit | Meaning |
|---|---|---|
| non-commit command (`ls -la`) | 0 | allow |
| valid message | 0 | allow |
| forbidden word (mid-message or END) | 2 | BLOCK |
| bad subject format / action word | 2 | BLOCK |
| `.claude/` or `docs/claude-*` path context | 0 | allow (exempt) |
| `git commit --amend` (no message) | 0 | allow (cannot judge) |
| doc heredoc mentioning `git commit` | 0 | allow (not a commit) |

### Files changed this milestone

- `.claude/hooks/commit_message_check.py` — future import removed; **boundary bug fixed**;
  **extraction false positive fixed**
- `.claude/hooks/file_placement_check.py` — future import removed; REPO derived
- `.claude/settings.json` — `env` block added
- `.claude/evals/programmatic/run_*.py` (5) — future imports removed
- `.claude/evals/programmatic/run_build_run_env_eval.py` — expectations derived from REPO
- `.claude/evals/programmatic/run_commit_hook_eval.py` — 10 new cases (22 -> 32)
- `.claude/evals/programmatic/build-run-env.md`, `commit-hook.md` — specs updated

Note `.claude/settings.json` still pins `"model": "fable"`, which overrides the session
default on restart.

---

## 2026-07-30 — inqkit test suite: 32/32 PASS on A100 (jobs 32393114 + 32393669)

### Result

| Tier | Tests | Pass | Fail | Skip | Build | ctest |
|---|---|---|---|---|---|---|
| `pure` (CPU, no INQ link) | 12 | **12** | 0 | 0 | 14 s | 0.8 s |
| `engine` (GPU, links INQ, sm_80) | 20 | **20** | 0 | 0 | 438 s | 66 s |
| **Total** | **32** | **32** | **0** | **0** | 7.5 min | 67 s |

Job 32393114, `gpu-q-13`, A100-SXM4-80GB, elapsed 00:08:47. Built against **stock
`inq`** (no engine test needs the CAP; `test_mask_shape` is pure-tier). ~1,356 Catch2
assertions total. Engine timings 1.04–8.18 s each; slowest `density_total` (8.18),
`plane_screen_parallel` (7.75, MPI×2), `density_semantics` (6.11), `plane_screen` (5.41).

Reproduce:
```
# configure on the LOGIN node (compute nodes have no network — fault #3)
cmake -S inq-stack/tests/include -B inq-stack/tests/include/build \
      -DCATCH2_SOURCE_DIR=$PWD/inq/build/_deps/catch2-src
cmake -S inq-stack/tests/include/engine -B inq-stack/tests/include/engine/build \
      -DINQ_SOURCE=$PWD/inq -DINQKIT_INCLUDE=$PWD/inq-stack/include \
      -DENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=80 \
      -DCMAKE_CUDA_COMPILER=$INQ_NVCC -DPython_EXECUTABLE=$PWD/venv/bin/python \
      -DFETCHCONTENT_SOURCE_DIR_{CATCH2,PYBIND11,SPDLOG,SPGLIB}=...
# then, FROM THE REPO ROOT:
sbatch shared/bin/run-inqkit-tests.slurm
python3 shared/bin/collate-inqkit-tests.py
```

Both CMakeLists carry previous-device defaults for `INQ_SOURCE` / `INQKIT_INCLUDE`
(and a stale nvcc/venv path in an engine-CMakeLists comment). They are `CACHE PATH`,
so they are OVERRIDDEN with `-D` above rather than edited — but a future session that
configures without those flags will silently point at `/local/data/public/...`.

### New tooling

- `shared/bin/run-inqkit-tests.slurm` — builds + ctests both tiers on `ampere`, emits
  JUnit XML per tier. `set -uo pipefail` (NOT `-e`): a failing test must not abort the
  job before results are collected.
- `shared/bin/collate-inqkit-tests.py` — parses the JUnit XML into a table. Parsing XML
  rather than scraping ctest console text is deliberate: the console prints only counts,
  and a tier that never BUILT is invisible there. A tier with no XML is reported as
  MISSING, never as zero-and-passing. Self-checked pre-run (both tiers MISSING, exit 1).

### Fault #15 — REPO_ROOT from BASH_SOURCE fails under SLURM

First submission (32392850) died in **3 s**: SLURM copies the batch script to
`/var/spool/slurm/slurmd/job<ID>/slurm_script`, so `${BASH_SOURCE[0]}` resolved
`REPO_ROOT=/var/spool/slurm` and `csd3-env.sh` was not found. `build-inq.slurm` and
`run-dispersion.slurm` already used `${SLURM_SUBMIT_DIR:-$(pwd)}`; the new script now
matches, plus a marker-file guard that fails with an explicit "submit from the repo
root" message (verified by invoking it from `/tmp`). **The job must be submitted from
the repo root.**

### Anti-trivial-pass verification (job 32393669, `gpu-q-1`, 2:46)

An all-pass ctest log is NOT evidence the tests did anything — `--output-on-failure`
captures output only for FAILING tests, so a green log contains no GPU banner by
construction. Checked separately:

- **Assertion counts**: 19/20 report `All tests passed (N assertions)` with N ≥ 2; none
  zero. Largest: `plane_screen` 811, `leed` 402, `free_wp` 22, `eigenvalue_dump` 21.
- **`test_plane_screen_parallel_engine` reports 1 assertion** — this is CORRECT, not a
  red flag: the source contains exactly one `CHECK` (cross-rank variance of the
  plane_screen slice signature, the E01 all_reduce invariant). Its missing
  "All tests passed" prefix is stdout interleaving from 2 MPI ranks; ctest and mpirun
  both returned success, which cannot happen on a failed assertion.
- **GPU execution**: INQ logged `process 0 has gpu id wE3O/g3iLvRBXiI0YMOdrg`.
  Binaries confirmed `NEEDED libcufft.so.11 / libcudart.so.12` (correct 12.x toolkit,
  not the 11.4 that caused fault #6) and carry **sm_80** cubins (`cuobjdump --list-elf`).
- **INCONCLUSIVE, not negative**: an `nvidia-smi --query-compute-apps` poll during a run
  returned nothing. The polling loop was too coarse for a ~8 s test. The GPU-id banner +
  sm_80 cubins + 12.x linkage are the evidence; the nvidia-smi check proves nothing
  either way and should not be cited.
- **Invalid attempt, recorded so it is not repeated**: assertion counts were first
  probed by running the engine binaries ON THE LOGIN NODE, which has no GPU — they
  produced no output. That was not evidence of anything.

### SCOPE — what 32/32 does and does NOT cover

DOES: inqkit's own logic (grid/FFT-shift indexing, field + VTI writers, jellium shells,
absorber mask geometry, projectile kinematics) and its integration with the engine —
15 of 20 engine tests run a real `ground_state::calculate`, 3 also `real_time::propagate`,
1 exercises the MPI decomposition.

Does NOT:
- **No coronene test exists.** See below.
- Hartree/XC/pseudopotential correctness at large — INQ's OWN `ctest` suite is still
  UNRUN and remains the outstanding broad check on the `fixed_return` CUDA patch.
- Nothing here is bit-identical (see below).

### Coronene "bit-to-bit replication" — DOES NOT EXIST (user expectation mismatch)

The user asked to run coronene bit-to-bit replication tests. Searched the whole tree:

- **No coronene test source anywhere.** The only mention in the test tree is
  `test_minimum_observable_set.cpp`, asserting the observable NAME LIST for
  `RunType::coronene` contains `wp_momentum_stats`. No simulation runs.
- The coronene work is a **paper replication** of Tsubonoya, Hu & Watanabe,
  Phys. Rev. B **90**, 035416 (2014), via `inqkit/config/tsubonoya_2014_coronene.hpp`,
  consumed by production runs in `ResearchProject/systems/coronene/` — not by any test
  target.
- `docs/validation/coronene-replication.md` is a Tier A/B/C checklist with **all 13 runs
  × 10 checks unfilled (`—`)**; Tiers B and C marked deferred. There is no recorded
  coronene result on any device, hence no reference to compare against.
- **Nothing in the suite is bit-identical.** The closest, `wp_real_space_replica_engine`,
  is an independent in-test recomputation compared with `margin(0.05)`/`margin(0.15)`.
  Bit equality would be the WRONG target for GPU reductions, whose summation order is
  not guaranteed stable run-to-run — such a test would be flaky by construction.

To get real coronene regression coverage: pin a reference run, store its observables,
assert against a physically-justified tolerance. The tolerance is a deliberate physics
decision, not a bit comparison. NOT started.

---

## 2026-07-30 — INQ's OWN suite: 249/249 (247 in one sweep + 2 after _pinq fix)

### Headline

`ctest` on INQ itself — the broad correctness check on the local CUDA patch — is DONE.

Job 32394792, `gpu-q-22`, A100-80GB, QOS `gpu1`, elapsed **00:29:44**, test time 1777 s,
serial (no `ctest -j`: tests share one GPU; concurrency risks false failures).
**247/249 passed (99%).** The 2 failures were `inq::tests::gaas.py` and `na2+.py`,
BOTH caused by this session's venv rebuild — no C++/CUDA/physics failure anywhere.
After the `_pinq` fix (job 32396491) those 2 pass, so **all 249 pass**.

NOTE ON RIGOUR: 247 passed in a single sweep; the remaining 2 were verified by
targeted rerun (`ctest -R '\.py'`), NOT by a second full sweep. A confirmatory
all-249 run has NOT been done. Cheap (~30 min) if a single clean sweep is wanted
for the record.

**No build was needed** — job 32355421 had already produced `src/all_unit_tests`
(98 MB; runs all 231 unit tests via Catch2 tag filters — they are NOT 231 separate
targets) plus the 11 executables in `tests/`. An earlier claim in this session that
"zero test executables exist" was WRONG: the check globbed `unit_tests*`, which does
not match `all_unit_tests`. A corresponding estimate of "~1.5 h to build 249 CUDA TUs"
was also wrong; runtime was the only cost.

### Results by group (from `inq/build/inq-ctest.log`)

| Group | Pass | Fail | Time (s) |
|---|---|---|---|
| `inq::tests` (33; 19 shell-driven SCF/TDDFT) | 31 | 2 → **0** | 1298.1 |
| `unit_tests::hamiltonian` | 15 | 0 | 41.3 |
| `unit_tests::operations` | 19 | 0 | 34.4 |
| `unit_tests::perturbations` | 12 | 0 | 33.9 |
| `unit_tests::interface` | 25 | 0 | 28.0 |
| `unit_tests::systems` | 3 | 0 | 16.0 |
| `unit_tests::parallel` (MPI) | 12 | 0 | 15.4 |
| `unit_tests::basis` | 10 | 0 | 14.3 |
| `unit_tests::solvers` | 9 | 0 | 12.8 |
| `unit_tests::observables` | 7 | 0 | 12.2 |
| `unit_tests::ground_state` | 5 | 0 | 11.7 |
| `unit_tests::utils` | 11 | 0 | 11.0 |
| `unit_tests::matrix` | 7 | 0 | 10.0 |
| **`gpurun::unit_tests::gpu`** | **7** | **0** | 9.4 |
| `unit_tests::eigensolvers` | 1 | 0 | 8.4 |
| `unit_tests::ionic` | 6 | 0 | 7.8 |
| `unit_tests::real_time` | 5 | 0 | 6.2 |
| `pseudopod::unit_tests::*` | 21 | 0 | 43.1 |
| `libpaw` / `run_paw_lib` / `spglibtest` | 3 | 0 | 2.4 |
| other unit groups (magnitude, math, states, options, mixers, parse, input, inq, bomd, config, physics) | 31 | 0 | ~31 |
| **TOTAL** | **247** | **2 → 0** | 1663 |

**`gpurun::unit_tests::gpu` 7/7 is the key row** — those are the tests for the module
the `fixed_return` patch rewrote (`external_libs/gpurun/include/gpu/reduce.hpp`).
Everything layered above also passes: Hamiltonian, operations, solvers, eigensolvers,
ground_state, real_time, perturbations, PAW, pseudopod, MPI parallel.

Heaviest physics, all passing: `electron_gas_stress` 202.8 s, `electron_gas.sh` 116.8,
`h2+_absorption.sh` 100.4, `diamond_hybrid.sh` 88.0, `al4h1.sh` 80.1,
`nitrogen_non_diagonal.sh` 75.8, `hydrogen_local.sh` 64.5, `silicon.sh` 59.7,
`silicon_hartree_fock.sh` 50.7, `oxygen.sh` 48.5. Covers hybrid functionals,
Hartree-Fock, stress, absorption — none touched by inqkit tests or the vacuum sweep.

### Fault #16 — venv rebuild broke `_pinq` (self-inflicted, and I dismissed the risk)

`inq/build/python/_pinq.cpython-36m-x86_64-linux-gnu.so` was built 2026-07-29 18:19
against the THEN-current venv (3.6.8). The venv was rebuilt to **3.11.13** on
2026-07-30, and a `cpython-36m` extension cannot be imported by 3.11 →
`ModuleNotFoundError: No module named '_pinq'`.

**I had explicitly reassured the user this was safe**, reasoning that
`ENABLE_PYTHON:BOOL=OFF` meant the bindings were not built. Two errors:
1. The artefact was ON DISK; I read one cache variable and never checked the tree.
2. `ENABLE_PYTHON` gates INSTALLING the python API (cache docstring: "Install Python
   API interface"); the `_pinq` target builds regardless. So the flag never implied
   what I claimed.
My original caution (recorded earlier in this handover) was correct and I talked
myself out of it.

### Fixing it — attempt 1 FAILED and made it worse (job 32396176)

Attempt 1 `-U`'d the stale entries and set the MODERN `Python_EXECUTABLE` /
`Python3_EXECUTABLE`. But this build has **`PYBIND11_FINDPYTHON:BOOL=OFF`**, so
pybind11 uses the CLASSIC `FindPythonLibs`, which does not read those — both landed in
the cache as `UNINITIALIZED` and were ignored. With the 3.6 pin removed and nothing
replacing it, the classic search fell back to the OLDEST interpreter present:
```
PYTHON_EXECUTABLE   = /usr/bin/python2.7
PYTHON_INCLUDE_DIRS = /usr/include/python2.7      (PY_VERSION "2.7.18")
PYBIND11_PYTHON_VERSION = 2.7
```
→ `pybind11/detail/common.h:277: error "PYTHON < 3.6 IS UNSUPPORTED."`
i.e. 3.6-wrong became 2.7-rejected. The `.so` was never replaced, so nothing was lost.

### Attempt 2 — SUCCESS (job 32396491, 4:10)

Pin the CLASSIC trio that `FindPythonLibs` actually reads (all three verified present
before submitting; the venv is built on system 3.11 so its headers ARE the system ones,
confirmed via `sysconfig.get_paths()['include']`, not assumed):
```
-DPYTHON_EXECUTABLE=<repo>/venv/bin/python
-DPYTHON_INCLUDE_DIR=/usr/include/python3.11      # PY_VERSION "3.11.13"
-DPYTHON_LIBRARY=/usr/lib64/libpython3.11.so
-DPYBIND11_PYTHON_VERSION=3.11
```
Guard confirmed resolution (3.11.13) BEFORE compiling; built in 174 s; produced
`_pinq.cpython-311-x86_64-linux-gnu.so` (50 MB, 12:38). Import smoke test:
`_pinq OK`, `pinq OK`.

Rerun of all 4 python cases — **4/4 pass**:

| Test | Before | After | Note |
|---|---|---|---|
| `diamond.py` | Passed 0.22 s | Passed 0.10 s | **still SKIPS — no ASE** |
| `gaas.py` | **Failed 0.08 s** | **Passed 23.13 s** | real calculation |
| `h2_ase.py` | Passed 0.16 s | Passed 0.06 s | **still SKIPS — no ASE** |
| `na2+.py` | **Failed 0.10 s** | **Passed 8.03 s** | real calculation |

The runtime jump (0.08 s → 23.13 s; 0.10 s → 8.03 s) is the evidence the fix is real:
those seconds are actual DFT calculations, not a trivial pass.

Script hardening added to `shared/bin/rebuild-pinq.slurm` after attempt 1:
- **Post-configure guard** reads back `PYTHON_VERSION`/`PYTHON_INCLUDE_DIRS` from the
  CMake CACHE (the authority on what nvcc receives) and aborts if < 3.6, rather than
  spending 4 minutes compiling toward a knowable failure. Checking intent instead of
  the cache was the original mistake.
- **Honest failure reporting**: attempt 1 printed only `gmake: *** [_pinq] Error 2`
  because the grep matched `'Error '` first and `head -20` truncated before the real
  pybind11 `#error`. Now surfaces compiler diagnostics first plus the log tail.

### Still-silent tests (NOT fixed, by design decision pending)

`diamond.py` and `h2_ase.py` report **Passed** without exercising the bindings — both
open with
```python
ase_spec = importlib.util.find_spec("ase")
if(ase_spec is None): print("Cannot find ASE, this test will be skipped"); exit()
```
ASE is not installed in the new venv, so they `exit()` at line 5 with status 0 and
ctest counts a pass. So the python tier reads 4/4 while only **2 of 4** genuinely test
anything. `pip install ase` would make all four real; `h2_ase.py` also imports
`ase.calculators.nwchem`, so it may need more than ASE alone. NOT done — awaiting the
user's decision, since it changes which tests actually execute.

### Fault #17 — sbatch needs the module env (REPEAT of a documented fault)

A chained "wait for job, then submit" background command failed with
`sbatch: error: plugin_load_from_file: dlopen(lua.so): liblua-5.1.so`. `sbatch` must
run inside the `rhel8` module environment; the waiter only put slurm on `PATH` without
sourcing `csd3-env.sh`. **This fault was already recorded earlier in this handover and
was repeated anyway.** Any automation that shells out to `sbatch` MUST source
`shared/bin/csd3-env.sh` first.

### Leftover

`inq/build/python/_pinq.cpython-36m-x86_64-linux-gnu.so` (95 MB, stale 3.6 build) is
still present alongside the new 3.11 module. Harmless — Python selects by ABI tag —
but deletable to save space and avoid confusion.

### Cumulative validation state (this device)

| Suite | Count | Result |
|---|---|---|
| INQ own `ctest` (serial) | 249 | **249 pass** (247 sweep + 2 post-fix) |
| inqkit `pure` | 12 | 12 pass |
| inqkit `engine` (GPU) | 20 | 20 pass |
| `inqview` python | 163 | 163 pass, 8 skip, 5 xfail, 1 xpass |
| programmatic evals | 5 suites | all pass (32/32 commit-hook) |
| vacuum free-Gaussian dispersion | 9 runs | PASS, worst dev 0.002 % |

NOT done: INQ's MPI ctest pass (`INQ_EXEC_ENV="mpirun.openmpi -np 4" ctest`) — the
serial pass covered `unit_tests::parallel` 12/12 in-process, but not the 4-rank
launcher path. Nothing is committed.
