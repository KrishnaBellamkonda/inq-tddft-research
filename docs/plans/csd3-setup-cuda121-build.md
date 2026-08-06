# Plan: CSD3 setup — make the engine build under CUDA 12.1

**Created:** 2026-07-29
**Trigger:** `setup.sh` failed on CSD3 (`/home/skcb2/setup.log`, 5499 lines) while
bootstrapping this repo on a new device.
**Branch:** `quantum-stopping-power`

---

## Diagnosis (verified, not inferred)

`setup.sh` ran on CSD3 login node `login-q-1`. Steps [1/4] (submodule) and [2/4]
(inq present, patch already applied) succeeded. Step [3/4] (build inq) failed with
**two independent** errors.

### Cause 1 — CUB fix incompatible with CUDA 12.1 (105 compile errors)

`inq-local.patch` rewrites `gpu/reduce.hpp:61` to

```cpp
thrust::transform_reduce(..., cuda::proclaim_return_type<Type>(kernel), init, std::plus<Type>{});
```

`setup.sh:11` documents this as "the CUB fix **for CUDA 12.5+**"; the previous
device ran CUDA 12.6.2 (`/lsc/opt/cuda-12.6.2/bin/nvcc`, per `CLAUDE.md`).

**This cluster's newest CUDA toolkit is 12.1.** Verified: `/usr/local/software/cuda/`
tops out at 12.1; `module avail cuda` under both the rhel7 default env and
`rhel8/default-amp` tops out at `cuda/12.1`; `nvhpc` is 21.x (CUDA 11 era). No
toolkit >= 12.4 exists on CSD3.

CUDA 12.1's libcu++ `cuda::__detail::__return_type_wrapper`
(`/usr/local/software/cuda/12.1/include/cuda/functional:73-78`) declares an
**unconstrained greedy template constructor** and **no copy constructor**:

```cpp
template <class _Fn>
explicit __return_type_wrapper(_Fn&& __fn) noexcept : __fn_(forward<_Fn>(__fn)) {}
```

When thrust's `transform_input_iterator_t` copy-constructs the wrapper, that
template ctor is a better match than the implicit copy ctor (non-const lvalue),
so it attempts to construct the *lambda* from a *wrapper*:

```
cuda/functional(77): error: no instance of constructor "lambda [](auto)->auto::<unnamed>"
                     matches the argument list
  argument types are: (cuda::__4::__detail::__return_type_wrapper<double, lambda [](auto)->auto>)
```

libcu++ constrained this constructor in a later CUDA release, which is why the
same patch builds on 12.6 and fails on 12.1.

**Verified by reproducer** (`scratchpad/repro.cu`) with INQ's own flags from
`inq/CMakeLists.txt:47` (`-std=c++17 --extended-lambda --expt-relaxed-constexpr
-arch=sm_80`) against `/usr/local/software/cuda/12.1/bin/nvcc`:

| Variant | Result |
|---|---|
| `cuda::proclaim_return_type<Type>(kernel)` (current patch) | FAILS — identical `cuda/functional(77)` error |
| plain functor with fixed return type + implicit copy ctor | **COMPILED OK** |

### Cause 2 — OOM kill on the login node (16x `cicc died due to signal 9`)

`nvcc error : 'cicc' died due to signal 9 (Kill signal)` = SIGKILL from the cgroup
OOM killer, surfacing as `gmake ... Error 9`.

- `setup.sh:56` calls `cmake --build "$ROOT/inq/build" --parallel` with **no job
  count** -> as many jobs as cores.
- `login-q-1` has 76 cores; each `cicc` on INQ's template-heavy TUs needs several GB.
- Per-user login-node cgroup limit measured at
  `/sys/fs/cgroup/memory/user.slice/user-<uid>.slice/memory.limit_in_bytes` = **20 GB**.
- `login-q-1` also has **no GPU** (`nvidia-smi` fails), so `ctest` cannot run there
  regardless.

### Side issue (fixed before diagnosis, unblocking)

`.claude/settings.json` PreToolUse hook pointed at the old device path
`/local/data/public/skcb2/tddft/.claude/hooks/commit_message_check.py`, which does
not exist here — every Bash call aborted. Repointed to
`$CLAUDE_PROJECT_DIR/.claude/hooks/commit_message_check.py`.

---

## Decisions (user, 2026-07-29)

1. Apply the portability fix to `inq-local.patch`, `inq/`, **and** `inq-study/`.
   (Explicitly overrides `.claude/rules/inq-immutable.md` for this change;
   `inq-local.patch` is the sanctioned delta mechanism per `setup.sh:44`.)
2. Build via SLURM on the `ampere` partition, account
   `mphil-nikiforakis-skcb2-sl2-gpu`, using the `rhel8/default-amp` module env.

---

## Steps

1. [x] Fix the broken `.claude/settings.json` hook path.
2. [ ] Rewrite the `reduce.hpp` hunk in `inq-local.patch` to use a local
       `fixed_return` functor instead of `cuda::proclaim_return_type`.
       This **preserves** the CUB fix's intent (pinning the reduction value type so
       thrust/CUB deduces a concrete type rather than `auto`) on CUDA 12.5+ while
       also working on 12.1 — it is NOT a revert, so `CLAUDE.md`'s
       "do not revert the CUB fix" constraint holds.
3. [ ] Apply the same edit to `inq/external_libs/gpurun/include/gpu/reduce.hpp`
       and `inq-study/external_libs/gpurun/include/gpu/reduce.hpp`.
4. [ ] Harden `setup.sh`: bounded `--parallel` (env-overridable), and document the
       CSD3 module env.
5. [ ] Add a SLURM build script and submit it to `ampere`.
6. [ ] Record validation status; write the handover.

## Portability fix (exact form)

```cpp
// Pin the reduction value type without cuda::proclaim_return_type, whose libcu++
// implementation before CUDA 12.4 has an unconstrained ctor that breaks when
// thrust copies the functor. Plain struct => implicit copy ctor, fixed return type.
template <typename Type, typename KernelType>
struct fixed_return {
	KernelType kernel;
	template <typename... Args>
	GPU_FUNCTION Type operator()(Args... args) const { return kernel(args...); }
};
```

`GPU_FUNCTION` is already in scope in `reduce.hpp` (used by `array_access`, line 42).

---

## Validation

| Tier | Check | Status |
|---|---|---|
| A | Reproducer: functor form compiles under CUDA 12.1 + INQ flags | **PASS** (done) |
| B | Full `inq` build completes on `ampere` | pending |
| C | `ctest` on an A100 node | pending — user approval before the expensive sweep |

Nothing is claimed correct until Tier B passes. Tier C is the expensive tier and
needs user sign-off per `.claude/rules/validation-gates.md`.
