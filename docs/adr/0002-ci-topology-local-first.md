# 2. CI is local-first (git hook on the GPU box); hosted GitHub Actions is a deferred pure-tier backstop

Date: 2026-06-08
Status: accepted

## Context

The inqkit/inqview test suite spans tiers with very different runtime needs:
- **pure** tests: C++17-only or numpy-only; run anywhere in seconds.
- **engine** / integration tests: link the INQ engine (CUDA/MPI/libxc) or
  import CUDA/VTK Python paths; require INQ to be built. GPU tests *must* run
  on the lab box.

GitHub-hosted runners have no CUDA GPU and no INQ build, so they cannot run the
engine/integration tier. The project rule also says "use GPU whenever
available." The repo has a remote (`KrishnaBellamkonda/inq-tddft-research`).

## Decision

- **Primary CI is local**, on the GPU box: a git **pre-push hook** plus
  `make`/`ctest`/`pytest` targets run the **full** suite (pure + engine +
  integration, CPU or GPU) on every change before it leaves the machine.
- **Hosted GitHub Actions is designed-for but deferred.** When added, it runs
  only the **pure** tier (`ctest -L pure`, `pytest -m "not engine"`) as a fast
  green/red backstop on push/PR. GPU/engine tests never run in cloud CI.
- Tests are tagged from day one (ctest `pure`/`engine` labels; pytest
  `engine` marker) so the hosted job is a drop-in later with no rework.
- No self-hosted runner: it would execute arbitrary PR code on the compute
  node, needs IT sign-off, and couples green-CI to box uptime.

## Consequences

- "Runs automatically on every change" is satisfied locally for all tiers;
  the cloud only ever sees what it can actually run.
- The pure/engine tagging is mandatory, not optional — it is what keeps the
  future hosted job honest and cheap.
- Trade-off: until the hosted backstop exists, there is no off-box CI signal;
  a contributor without the GPU box cannot get green CI. Accepted for a
  single-box research repo.
