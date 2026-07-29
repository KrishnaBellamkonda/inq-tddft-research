# INQ-TDDFT — electronic stopping-power research codebase

Real-time TDDFT / Ehrenfest study of **electronic stopping power** in jellium and
localised-jellium systems, built on the [INQ](https://gitlab.com/npneq/inq) GPU
engine. This branch is a **self-contained snapshot of the codebase** — the
library, the research run machinery, and the engine deltas needed to build it.
Run outputs, ground states, and notebooks are intentionally excluded (they are
regenerated, not shipped).

## Layout

| Path | What it is |
|---|---|
| `inq-stack/include/inqkit/` | **inqkit** — C++ header layer over INQ: field/orbital extraction, VTI/observables I/O, wavepacket injection, projectile dynamics, mask absorber, jellium energetics |
| `inq-stack/python/inqview/` | **inqview** — Python post-processing: loaders, analysis kernels (Fourier, Lindhard, energy decomposition), density-GIF + figure rendering |
| `inq-stack/tests/` | Catch2 (C++) and pytest (Python) suites |
| `ResearchProject/systems/<name>/` | Per-system run machinery: `scripts/` (`run.cpp`, `analyse.py`, dispatchers), `shared/` configs, pseudopotentials, `hypotheses/` analysis |
| `inq-study/` | **Submodule** — project-modified INQ engine with complex-absorbing-potential (CAP) support (stock upstream cannot compile a CAP run) |
| `setup.sh` + `inq-local.patch` | Engine bootstrap: clone pinned upstream INQ, apply the two sanctioned deltas (CUB fix + `ham()` accessor), build |

## Setup

```bash
git clone --recurse-submodules <this-repo>
cd <repo>
bash setup.sh          # clones+patches+builds inq, inits inq-study, installs inqview
```

Then add to your shell profile (paths printed by `setup.sh`):

```bash
export PATH="$PWD/shared/bin:$PATH"
export INQ_SHARE_PATH=$PWD/inq/install/share
export PSEUDOPOD_SHARE_PATH=$PWD/inq/install/share/pseudopod
```

## Running

Builds use the `inq-run` wrapper (GPU by default). Runs that load a converged
ground state expect it under the system's `shared_gs/` — regenerate these on your
machine with the `ResearchProject/systems/**/save_gs/*/run.cpp` builders before
launching (ground states are large and not shipped). Post-process with a run's
`analyse.py`, which drives the `inqview` pipeline.
