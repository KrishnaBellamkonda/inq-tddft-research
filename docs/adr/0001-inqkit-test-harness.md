# 1. inqkit C++ tests live outside the headers, in tests/cpp with a pure/engine ctest split

Date: 2026-06-08
Status: accepted

## Context

inqkit is a header-only C++17 library on top of INQ. INQ's own convention is
to embed Catch2 `TEST_CASE` blocks *inside* each source header, compiled by a
per-file unit-test driver. We need a test harness for inqkit's 38 headers.

Two facts constrain the choice:
- The task workspace statement says a **new folder of tests** is added to
  inq-stack, and project rule "never edit a file not required by the plan"
  discourages mutating all 38 headers just to host tests.
- ~5 headers are **pure** (compile with a bare C++17 compiler: schemas,
  jellium analytics, text parsing, validation predicates); the other ~33 are
  **engine-coupled** (take INQ `basis`/`field` objects, so a test must link
  the full INQ engine — CUDA/MPI/libxc/fftw).

## Decision

- Tests live in `inq-stack/tests/cpp/`, one `test_<header>.cpp` per header,
  mirroring the `include/inqkit/` tree. Headers are left untouched.
- Framework: **Catch2** (already vendored via INQ; same `Approx`/`_a` idiom).
- Build: a `CMakeLists.txt` registers **one ctest test per file**, tagged with
  a label — `pure` (links only Catch2, C++17) or `engine` (links INQ via the
  same discovery `inq-run` uses).
- Selection: `ctest -L pure` runs anywhere; `ctest -L engine` runs only where
  INQ is built (CPU or GPU box).

## Consequences

- Cheap pure tests can run on a GPU-less CI runner; the heavy engine tests are
  isolated to where INQ exists. This split is the foundation of the CI design.
- Per-file ctest granularity supports the "one component at a time, reviewed
  before the next" rule.
- inqkit headers stay free of test code, but tests are physically separated
  from the code they cover (mitigated by the mirrored directory layout).
- Rejected: embedding `TEST_CASE` in each header (edits all 38, fights the
  never-edit rule); a single monolithic test binary (forces the INQ link into
  pure tests, coarse granularity).
