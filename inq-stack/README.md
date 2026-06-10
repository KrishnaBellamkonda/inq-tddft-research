# inq-stack

A project-local library layer on top of [INQ](https://gitlab.com/npneq/inq) for
real-time TDDFT wave-packet / scattering studies. Two components:

- **inqkit** (`include/inqkit/`) — header-only C++17 extraction layer over INQ:
  field/density/orbital extraction, observable writers, wave-packet injection,
  LEED screens. Include alongside `<inq/inq.hpp>`.
- **inqview** (`python/inqview/`) — Python post-processing & visualisation,
  installable with `pip install -e .` (extras: `[analysis]`, `[viz]`, `[test]`).

## inqview public API (4 layers, ADR 0003)

Dependency-layered so a headless node can compute observables without plotting
deps. `import inqview.analysis` / `inqview.io` pull in **no matplotlib and no VTK**.

| Layer | Import | Contents |
|---|---|---|
| `inqview.io` | numpy only | `RealField3D`, `ComplexField3D`, `SimulationData`, `load_leed_pattern` |
| `inqview.analysis` | numpy/scipy/pandas | `FourierTransform`, `energy_components`, `wp_integrity`, `plasmon_spectrum`, `center_of_density`, `kl_divergence` — each `compute(...) -> frozen dataclass` |
| `inqview.visualisation` | matplotlib/VTK | renderers (`render_*`), the canonical theme (`style`, ADR 0004), VTI/paraview |
| `inqview.pipeline` | the above | thin per-phase orchestration over a run's `results/` |

Public names are also re-exported lazily from the top level
(`from inqview import RealField3D`); `inqview.postprocess` is a **deprecated**
shim for `inqview.pipeline`.

```python
import inqview.analysis as A
ec = A.energy_components.compute("results/raw/observables/observables.csv")
print(ec.redistribution_ev())          # where did the energy go? (eV)
```

## Build & test

```bash
# C++ tests (Catch2; pure tier needs no INQ, engine tier links INQ + GPU)
cmake -S tests/cpp -B tests/cpp/build -DENABLE_CUDA=ON && ctest --test-dir tests/cpp/build

# Python suite (portable, numpy-only)
pip install -e ".[analysis,viz,test]"
pytest python/tests
```

Runs use the `inq-run` wrapper (see `../docs/compilation.md`). Architecture,
decisions, and the observable contract: `../docs/adr/`, `../CONTEXT.md`,
`../docs/observables/`.
