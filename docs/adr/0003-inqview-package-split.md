# 3. inqview splits into io / analysis / visualisation / pipeline, with a deps-clean analysis layer

Date: 2026-06-10
Status: accepted

## Context

`inqview` is the Python post-processing + visualisation package
(`inq-stack/python/inqview/`). It grew organically into three de-facto
layers that were never named as such:

- **primitives** — top-level modules (`fields`, `data`, `fourier`,
  `screens`/`LeedPattern` loader, `overlap`, `vti`, `plots`, `paraview`).
- **pipeline phases** — `inqview/postprocess/*.py`, ~23 modules each
  bundling *three* jobs: numeric compute, plotting, and artefact writing.
- **applications** — `report1/**`, `scripts/**` (one-off figure/CLI code).

The bundling inside the middle layer is the core problem the user flagged
("postprocess and visualisation have to be separated"): e.g. `wake.py`
mixes `bath_volume()` (compute) with `shared_clim()` (a colour-scale viz
helper), and every phase imports `matplotlib`/`vtk` at module top. The
library is about to be **released to the scientific community**, so the
public import surface and dependency footprint now matter.

A concrete pain point: a user on a headless cluster node who only wants
the *numbers* (a loss function, a KL-divergence time series, a stopping
power) must today import modules that pull in matplotlib and VTK.

## Decision

Reorganise `inqview` into four dependency-layered sub-packages:

- **`inqview.io`** — loaders and field/format dataclasses. numpy only.
  (e.g. `data`, `fields`, `vti` read side, the `LeedPattern` loader.)
- **`inqview.analysis`** — numeric post-processing kernels. Imports **only**
  numpy/scipy; returns plain dataclasses / ndarrays. Each kernel exposes a
  `compute(...) -> Result` function. (e.g. `fourier`, `wake` bath math,
  `density_fourier` loss function, `kl_divergence`, `overlap`, `stopping`,
  `lindhard`.)
- **`inqview.visualisation`** — all rendering: matplotlib, VTK/paraview,
  GIF assembly, colour-scale helpers. The **only** layer permitted to
  import plotting/VTK libraries. Each renderer takes an analysis `Result`.
- **`inqview.pipeline`** — thin phase orchestration: `compute → plot →
  write artefact`. Holds no maths and no bespoke plotting; wires
  `analysis` to `visualisation` and manages files/logging.

**Invariant (the decisive property):** `import inqview.analysis` (and
`inqview.io`) must pull in **no** matplotlib and **no** VTK. This is
enforced, not conventional — a test asserts the import graph stays clean.

`report1/**` and `scripts/**` remain out of the public API; they are
applications that consume the four packages and are relocated, not tested.

## Consequences

- **Public API for the community** becomes legible: numbers from
  `inqview.io` + `inqview.analysis`, pictures from `inqview.visualisation`,
  turnkey runs from `inqview.pipeline`.
- **Headless/cluster use** works with a numpy/scipy-only install; plotting
  deps become optional extras.
- The deps-clean invariant is **testable** — a pure-tier import-graph test
  fails if any `analysis`/`io` module imports matplotlib or vtk. This is
  one of the first inqview tests to lock.
- **Migration cost**: ~23 phase modules must be split into a `compute()`
  kernel (→ analysis) and a renderer (→ visualisation), with the phase
  reduced to orchestration. Done incrementally; the suite's
  characterization tests guard behaviour through each move.
- The `screens.py` / `overlap.py` **name collisions** (loader vs phase)
  are resolved by the new package boundaries (loader → `io`, phase →
  `pipeline`), not by deletion — they were never true duplicates.
- Trade-off: more packages and a stricter import discipline than a flat
  module list. Accepted because the release audience and the
  compute-without-plotting use case both demand it.
