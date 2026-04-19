# Handover: Session Restart / Repository Overview

*Updated: 2026-04-18. Use this as the first file to read in a fresh session.*

---

## Current status

Two active workstreams:

1. **Coronene LEED simulation** — `ResearchProject/systems/coronene/04_leed_simulation/`
   - run_004 (d=0.53 Å, Lz×1.5, 1561 steps) is **COMPLETE** as of 2026-04-16.
   - Analysis done: all 10 figures generated. Key findings in handover.
   - **Next**: interpret LEED k-space pattern; decide between periodic crystal array (Option B) or wider screen spacing for follow-up run.
   - Full details: `docs/handovers/coronene_wp_scattering.md`

2. **inqview VTI/ParaView pipeline** — `inq-stack/` on branch `features/python-paraview`
   - VTI writer (`vti.py`), ParaView batch renderer (`paraview.py`), and new fields layer (`fields.py`) are in active development.
   - Tutorial integration: `Tutorial/n2-with-inqkit/run.cpp` shows inqkit usage.
   - Status: code written; validation not yet recorded (needs test against known input before declaring done).

---

## What this repository is for

GPU-accelerated TDDFT simulations using the INQ C++ library, with a focus on:
- Electron wavepacket scattering from molecules (coronene LEED)
- Jellium HEG ground state and convergence benchmarks
- Developing `inqkit` (C++ extraction layer) and `inqview` (Python post-processing) on top of INQ

---

## Main active components

| Component | Path | Status |
|---|---|---|
| Coronene LEED runs | `ResearchProject/systems/coronene/04_leed_simulation/runs/` | run_004 done |
| inqkit C++ library | `inq-stack/include/inqkit/` | fields/io complete; core/screens stubs |
| inqview Python package | `inq-stack/python/inqview/` | VTI+ParaView in progress |
| Jellium ground state | `ResearchProject/jellium/01_ground_state/` | complete |
| Jellium convergence | `ResearchProject/jellium/02_ground_state_convergence/` | complete |
| Free WP propagation | `ResearchProject/jellium/03_free_gaussian_wp_propagation/` | complete |

---

## Repository structure (practical)

```
inq-stack/
  include/inqkit/    C++ headers — use alongside <inq/inq.hpp>
  python/inqview/    Python package — pip install -e inq-stack/

ResearchProject/
  jellium/           01–03 experiments (all complete)
  systems/coronene/  04_leed_simulation/ active; runs/run_001…run_005

Tutorial/
  n2-with-inqkit/    inqkit integration example (run.cpp)
  hello-world/ n2/ li-bcc/ HF/ ions/ minimal-gs/

docs/
  plans/     coronene_wp_scattering.md (active), jellium_research.md, jellium_gs_visualisations.md
  handovers/ coronene_wp_scattering.md (latest), jellium_ground_state.md
```

Not tracked (gitignored but present on disk): `inq/`, `venv/`, `QuantumKickExtension/`, `ParaView-6.1.0-*/`

---

## Most important docs to read first

1. `docs/handovers/coronene_wp_scattering.md` — run_004 results, next steps
2. `docs/folder_structure.md` — full directory map
3. `docs/compilation.md` — build workflow (`inq-run`)
4. `docs/inq_tutorial.md` — INQ C++ API
5. `docs/plans/coronene_wp_scattering.md` — parameters, validation checklist

---

## Build workflow

```bash
source ~/.bashrc          # sets PATH, INQ_SHARE_PATH, PSEUDOPOD_SHARE_PATH
cd <dir-with-run.cpp>
inq-run                   # GPU build + run (auto-detects .cpp)
inq-run --cpu             # CPU-only fallback
```

**Never activate the `quantum-wave-packet` pyenv when running `inq-run`.**  
For Python visualisations: use `venv/` (`source venv/bin/activate`) or specify `venv/bin/python3`.

---

## Key conventions

- `inq-run` auto-generates `CMakeLists.txt` — do NOT write one manually for INQ projects.
- INQ cell origin is (0,0,0) corner — all atom coordinates must be positive within [0,L].
- Figures saved as `.png` only (never `.pdf` or `.svg` unless explicitly requested).
- GPU execution preferred over CPU whenever available.
- Do NOT revert the CUB fix in `inq/external_libs/gpurun/include/gpu/reduce.hpp`.

---

## Current branch

`features/python-paraview` — VTI/ParaView pipeline for inqview.  
Main branch is `main`.

---

## Known issues / open tasks

1. `inq-stack` core/screens/real_time modules are stubs — not yet implemented.
2. inqview VTI/ParaView pipeline needs a known-case validation test before being declared complete.
3. `ResearchProject/systems/jellium/` is a legacy path (older experiments); `ResearchProject/jellium/` is canonical.
4. `docs/plans/jellium_gs_visualisations.md` has pending visualisation tasks (orbital cube files, XC offset plot, finite-size analysis) awaiting user confirmation.
5. Coronene LEED next run not yet decided — options documented in `docs/handovers/coronene_wp_scattering.md`.
