# Handover: inqkit/inqview I/O Validation — Tutorial Tutorials

*Updated: 2026-04-18*

---

## Current status

Both tutorial directories created and run. I/O pipeline validated on HF (clean) and coronene (norms pass, but SCF physics caveat — see below).

| Tutorial | SCF | N_elec | HOMO ρ norm | HOMO |ψ|² norm | PNG slices | VTI | ParaView |
|---|---|---|---|---|---|---|---|
| `hf-gs-with-inqkit` | ✓ Converged (28 iters, E=-24.528 Ha) | 8.000 ✓ | 1.0000 ✓ | 1.0000 ✓ | ✓ | ✓ | ✓ |
| `coronene-gs-with-inqkit` | ⚠ Converged to wrong minimum (178 iters, E=+302.29 Ha) | 108.000 ✓ | 1.0000 ✓ | 1.0000 ✓ | ✓ | ✓ | ✓ |

---

## What changed

- Created `Tutorial/coronene-gs-with-inqkit/` with `run.cpp`, `coronene_centered.xyz` (copied from run_004), `analysis.py`
- Created `Tutorial/hf-gs-with-inqkit/` with `run.cpp`, `analysis.py`
- Both built with `inq-run` (GPU, CUDA sm_80); results written via inqkit; loaded and validated with inqview

---

## Files touched

```
Tutorial/coronene-gs-with-inqkit/run.cpp
Tutorial/coronene-gs-with-inqkit/coronene_centered.xyz
Tutorial/coronene-gs-with-inqkit/analysis.py
Tutorial/coronene-gs-with-inqkit/results/density/density_total.{raw,meta.txt}
Tutorial/coronene-gs-with-inqkit/results/orbital_density/orbital_0053_density.{raw,meta.txt}
Tutorial/coronene-gs-with-inqkit/results/orbitals/orbital_0053_{real,imag}.raw + .meta.txt
Tutorial/coronene-gs-with-inqkit/results/visualisation/*.png  (3 slice PNGs)
Tutorial/coronene-gs-with-inqkit/results/visualisation/frames/{density,orbital_density}/frame_000000.png
Tutorial/coronene-gs-with-inqkit/results/visualisation/vti/*.vti

Tutorial/hf-gs-with-inqkit/run.cpp
Tutorial/hf-gs-with-inqkit/analysis.py
Tutorial/hf-gs-with-inqkit/results/  (same structure, orbital index 3)
```

---

## Commands run

```bash
# Build + run C++ GS simulations
cd Tutorial/hf-gs-with-inqkit && inq-run
cd Tutorial/coronene-gs-with-inqkit && inq-run

# Python analysis (inqview)
PYENV=/local/data/public/skcb2/pyenv/versions/3.10.19/envs/quantum-wave-packet/bin/python3
$PYENV Tutorial/hf-gs-with-inqkit/analysis.py
$PYENV Tutorial/coronene-gs-with-inqkit/analysis.py
```

---

## Tests and validation

### HF (all pass)
- N_electrons = 8.000 (expect 8) ✓
- HOMO density norm = 1.0000 ✓
- HOMO |ψ|² norm = 1.0000 ✓
- matplotlib slices, VTI, ParaView frames all written ✓

### Coronene (I/O pass; SCF physics caveat)
- N_electrons = 108.000 (expect 108) ✓
- HOMO density norm = 1.0000 ✓
- HOMO |ψ|² norm = 1.0000 ✓
- matplotlib slices, VTI, ParaView frames all written ✓
- **SCF caveat**: GS converged to E=+302.29 Ha vs run_004's +288.98 Ha (same cell, same settings). The Broyden mixer stalled in a higher-energy density minimum — dn still ~0.03 at termination despite de<1e-4. The density slices may appear diffuse rather than molecule-localised. The I/O pipeline is validated; the coronene physical result is not.

---

## Known issues / blockers

1. **Coronene SCF local minimum**: Tutorial run converges to E=+302.29 Ha instead of run_004's +288.98 Ha. Root cause unclear — same cell, same ions, same mixing parameters, same extra_states(3). Possible causes: (a) different random seed in initial_guess, (b) no restart from saved state. For a physically correct tutorial, start from a saved run_004 ground state, or investigate the mixing initialisation.

2. **Coronene density physically wrong**: The density slices will not show a clean ring-system pattern because the SCF didn't reach the true minimum. The orbital norms are correct (INQ normalises), but the shapes are unreliable.

3. **`pyenv activate` doesn't work in non-login bash shells**: Use the direct python path: `/local/data/public/skcb2/pyenv/versions/3.10.19/envs/quantum-wave-packet/bin/python3`

---

## Assumptions still in play

- HF uses LDA (user-confirmed) with 30 Ry cutoff, bond length 0.917 Å, atoms centred at z=0 in 8 bohr cubic finite cell.
- Coronene uses run_004 cell (34.771×34.771×89.856 bohr), LDA, 40 Ha cutoff, extra_states(3).
- HOMO = orbital index 53 (coronene, 0-based: 54 occupied) and index 3 (HF, 4 occupied).
- inqview installed in `quantum-wave-packet` pyenv.
- ParaView at `/local/data/public/skcb2/tddft/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch`.

---

## CPK sphere rendering (added 2026-04-19)

`AtomSpec` dataclass + CPK colours + VDW radii added to `inq-stack/python/inqview/paraview.py`.
`render_density_from_meta_series` and `render_vti_series` both accept `atoms: AtomSpec | None`.
`_paraview_batch_script()` creates `Sphere` sources per atom in the pvbatch script.

- Fix required: `ColorBy(d, None)` fails in ParaView 6.1. Use `d.ColorArrayName = ['POINTS', '']` instead.
- Background: VisRTX overrides `render_view.Background`; grey background is unavoidable in current PV 6.1 + VisRTX setup.
- Test: HF density render confirmed — green F sphere and white H sphere visible in `results/visualisation/frames/density/frame_000000.png`.

---

## Exact next steps

1. **Real-time writing**: next task is `features/real-time-orbital-density-writing` branch — see `docs/plans/inqkit_realtime_density_writing.md`.
2. **Publication 2D plots**: `plots.py` still a stub — see `docs/plans/inqview_publication_plots.md`.
3. **Coronene SCF**: physics still at wrong minimum; defer until production use requires it.
