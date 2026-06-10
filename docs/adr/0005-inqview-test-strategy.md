# 5. inqview test strategy: portable pure-Python suite, analytic + characterization, no renderer tests

Date: 2026-06-10
Status: accepted

## Context

inqview is about to be released to the scientific community, so its test
suite must be **runnable on an arbitrary user's machine** — no GPU, no INQ
engine build, no multi-GB data, and no failures caused by differing
hardware/BLAS/FFT backends. At the same time the suite must guard a large
behaviour-preserving restructure (ADR 0003/0004) and anchor the science.

Two naive options each fail:
- pure bit-exact characterization (`assert == captured golden`) breaks
  across machines because different BLAS/FFT backends reorder float ops;
- pure analytic tests can't cover real file formats or the orchestration.

## Decision

A **combination** strategy, all tiers pure-Python/numpy on committed
fixtures (the INQ engine and GPU are used *once by us* to generate
fixtures, never re-run in the suite):

- **I/O parsing tests** — read small committed data files and compare to
  golden. Parsing is deterministic and identical across machines, so exact
  comparison is allowed here.
- **Analysis-kernel tests** — assert against **analytic** expected values
  from a *simple reduced system whose response is known* (e.g. FFT of a
  pure sinusoid → known peak bin; KL of two known distributions → known
  nats; a reduced linear-response system → known loss function). Compared
  with **physical tolerances** (`np.allclose`), never bit-exact.
- **Integration test** — a **free-space wave-packet propagation** whose
  physics is analytically known (Gaussian spreading σ(t), ⟨p⟩ conserved,
  centroid = k₀·t/m). The run is executed once to produce a **small
  committed output fixture** + analytic golden; the test post-processes that
  committed output in pure numpy and matches the analytic expectation.
- **Visualisation renderers** — **no tests** (you cannot meaningfully
  unit-test a GIF/figure). The only visualisation-layer test is a tiny
  **numeric theme-config test**: `figure_one_col()==(3.5,3.0)`,
  `cmap_for('diverging')=='RdBu_r'`, key rcParams after `apply_theme()`.
  This guards the designed proportions (ADR 0004) without rendering.

**Portability is a hard acceptance criterion.** A test that can only pass on
the lab box (tight tolerance, GPU, local path, huge fixture) is a defect.

## Consequences

- External users get a green suite from `pip install` + `pytest` alone; no
  INQ, no CUDA, no data download.
- The free-space-WP integration fixture doubles as a **physics anchor**:
  golden is paired with an analytic law, so a regression that changes
  physics fails even though the capture would have matched.
- Layering cheap **physical invariants** (bath integral ≈ N_e, KL ≥ 0, WP
  norm ≈ 1) on top of golden turns the fixtures into bug-finders rather than
  bug-freezers.
- Fixtures must be small (target < 5 MB total): real file formats/layout but
  grid-downsampled and few-frame. Committed plainly (no git-lfs) if budget
  holds.
- The suite is **pure-tier** in the project's tier scheme: it never enters
  the engine tier, unlike inqkit. CI placement is trivial (any runner).
- Trade-off: golden values are tolerance-compared, so a test cannot catch a
  sub-tolerance numerical drift. Accepted — physical correctness, not
  bit-reproducibility, is the goal for a portable scientific library.
