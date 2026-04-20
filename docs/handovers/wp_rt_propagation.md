# Handover: WP Real-Time Propagation — free, jellium, coronene

## Current status

**Implementation complete on branch `features/leed-screen`.** All C++ and Python code
written and committed (commit 132d2b5). No runs have been started yet.

---

## What was done

### C++ library (inqkit)
- `inq-stack/include/inqkit/screens/plane_screen.hpp` — ported from `04_leed_simulation/utils.hpp`;
  extracts 2D ρ(x,y,z_screen) using GPU_SYNC() + hypercubic() pattern.
- `inq-stack/include/inqkit/screens/leed_pattern_accumulator.hpp` — time-integrates slice;
  saves 2-line header .dat file matching inqview.screens.py format.

### Python library (inqview)
- `inq-stack/python/inqview/screens.py` — `LeedPattern` dataclass + `load_leed_pattern()`
- `inq-stack/python/inqview/plots.py` — `plot_leed_pattern()` added
- `inq-stack/python/inqview/defaults.py` — `default_density_movie()` + `default_wavepacket_movie()`
- `inq-stack/python/inqview/__init__.py` — exports all new symbols
- All imports verified: `python -c "import inqview; ..."` passes.

### Simulation run files
All 19 runs written, each with `run.cpp` + `analysis.py`:

| System | Runs | Location |
|---|---|---|
| Free propagation | 7 | `Tutorial/free-propagation-wp-rt/run_{01..07}_*/` |
| Jellium N=40 | 6 | `ResearchProject/jellium/jellium-wp-rt/run_{01..06}_*/` |
| Coronene | 6 | `ResearchProject/systems/coronene/coronene-wp-rt/run_{01..06}_*/` |

Coronene runs have `coronene_centered.xyz` already copied in.

### Notes
- `docs/notes/wp_spreading_investigation.md` — WP spreading formula + LEED paper timing analysis

---

## Files touched

- `inq-stack/include/inqkit/screens/plane_screen.hpp`
- `inq-stack/include/inqkit/screens/leed_pattern_accumulator.hpp`
- `inq-stack/python/inqview/screens.py` (new)
- `inq-stack/python/inqview/defaults.py`
- `inq-stack/python/inqview/plots.py`
- `inq-stack/python/inqview/__init__.py`
- `Tutorial/free-propagation-wp-rt/run_0{1-7}_*/run.cpp` + `analysis.py`
- `ResearchProject/jellium/jellium-wp-rt/run_0{1-6}_*/run.cpp` + `analysis.py`
- `ResearchProject/systems/coronene/coronene-wp-rt/run_0{1-6}_*/run.cpp` + `analysis.py`
- `docs/notes/wp_spreading_investigation.md`

---

## Key implementation notes

### Two RT sessions pattern
Each run.cpp uses two `RealTimeSession` instances:
- `rt` — fires every `WRITE_EVERY=100` steps, writes density frames
- `rt_obs` — fires every step (write_every=1), writes observables + screen accumulation

The INQ propagate callback calls both: `[&](auto const& d){ rt.step(d); rt_obs.step(d); }`.

### t=0 density
`density::total(electrons)` does NOT include the WP orbital (extra state with occ=1 but not
counted as occupied in the GS density). Must explicitly add `density::orbital(electrons, ist_wp)`
element-wise before writing step 0.

### Free-prop: extra_states(1), no orthogonalisation
0 occupied states → no Gram-Schmidt needed. state_index = 0.

### Jellium: extra_states(3), full Broyden SCF
N=40 degenerate shell requires extra buffer states for Broyden convergence.
Temperature smearing 0.00862 eV (0.000317 Ha). WP orthogonalised against occupied states.

### Coronene: extra_states(3), 54 occupied (state_index=56)
Same setup as injection tutorial run_01. WP orthogonalised. 4 screens per run.

---

## Tests and validation

| Check | Status |
|---|---|
| inqview imports cleanly | ✓ verified |
| screens.py parser format matches leed_pattern_accumulator.hpp output | ✓ manual check |
| C++ files compile | NOT YET — awaiting build |
| RT propagation produces correct N_elec | NOT YET — awaiting run |
| Free-prop σ(t) matches analytic | NOT YET — awaiting run |

---

## Known issues / blockers

- IDE shows false-positive clang errors in all run.cpp files (`inq/inq.hpp` not found).
  This is expected — INQ cmake paths not in IDE config. Same in wavepacket.hpp. Not a real error.
- No runs have been built or started yet.

---

## Exact next steps

1. **Build and start free-propagation run_01_base first** (fastest validation path):
   ```bash
   cd Tutorial/free-propagation-wp-rt/run_01_base
   inq-run
   ```
2. If build fails, check INQ API: main risk is `LeedPatternAccumulator(PlaneScreen{...})` brace
   initialization — may need explicit constructor call `LeedPatternAccumulator(PlaneScreen(z, label))`.
3. After run_01 completes, run `python analysis.py` to verify N_elec=1 and σ(t) vs analytic.
4. If all pass, start remaining free-prop runs sequentially.
5. Then jellium runs (need full Broyden SCF ~300 steps before RT — expect ~2-4 h each).
6. Finally coronene runs (GS SCF ~100 iterations + 10000 RT steps — expect ~8-12 h each).
