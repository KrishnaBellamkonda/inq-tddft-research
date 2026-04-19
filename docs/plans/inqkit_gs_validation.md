# Plan: inqkit/inqview GS I/O validation — coronene and HF tutorials

*Created: 2026-04-17 | Status: COMPLETE*

---

## Context

The inqkit (C++) + inqview (Python) library stack was prototyped in
`Tutorial/n2-with-inqkit/` (single orbital + density, GS only). Before using inqkit in
production runs, the full pipeline (write → read → validate → visualise) needed to be
validated on physically richer systems. Two tutorial directories were created for this
purpose.

---

## Systems and directories

| System | Directory | Key parameters |
|---|---|---|
| HF | `Tutorial/hf-gs-with-inqkit/` | L=16 bohr cubic finite, E_cut=15 Ha, 4 occ. orbitals, HOMO=3 |
| Coronene | `Tutorial/coronene-gs-with-inqkit/` | LX=LY=34.771, LZ=89.856 bohr, E_cut=40 Ha, 54 occ. orbitals, HOMO=53 |

---

## What was built and tested

### C++ `run.cpp` (both tutorials)
- Ground state only (no real-time)
- Writes: total density, HOMO orbital density, HOMO wavefunction (complex)
- Uses `inqkit::fields::density::total()`, `density::orbital()`, `orbital::wavefunction()`
- Uses `inqkit::io::RealField3DWriter` and `ComplexField3DWriter`

### Python `analysis.py` (both tutorials)
- Loads fields with `load_real_field` / `load_complex_field`
- Validates norms: N_electrons, HOMO density norm, HOMO |ψ|² norm
- Matplotlib slice plots (3 planes × 3 fields)
- VTI conversion via `convert_real_meta_to_vti` and `write_vti`
- ParaView volume renders via `ParaViewPipeline`

### Publication-quality ParaView rendering (CPK spheres)
- Added `AtomSpec` dataclass to `paraview.py`
- CPK-coloured `Sphere` sources rendered alongside volume density
- HF atoms (H=white, F=green) confirmed visible in output PNGs
- See `docs/plans/inqview_publication_plots.md` for full spec

---

## Validation outcomes

| Check | HF | Coronene |
|---|---|---|
| SCF converged | ✓ (28 iter, E=-24.528 Ha) | ⚠ wrong minimum (178 iter, E=+302.29 Ha) |
| N_electrons | 8.000 ✓ | 108.000 ✓ |
| HOMO density norm | 1.0000 ✓ | 1.0000 ✓ |
| HOMO \|ψ\|² norm | 1.0000 ✓ | 1.0000 ✓ |
| Matplotlib slices | ✓ | ✓ |
| VTI files | ✓ | ✓ |
| ParaView frames | ✓ | ✓ |
| CPK spheres in render | ✓ (HF tested) | not tested |

---

## Key bugs discovered and fixed

1. **INQ corner-origin convention**: INQ uses (0,0,0) as the cell corner. Atoms inserted
   with negative coordinates or x=y=0 sit at the cell edge and may wrap. Documented in
   `docs/inq_tutorial.md` (ions chapter). Fix: always place atoms near (L/2, L/2, L/2).

2. **HF geometry (original)**: H atom at z=-0.459 bohr → wrapped to far side of cell.
   Fix: L=16 bohr, atoms centered at (8, 8, 8±0.866 bohr).

3. **Coronene SCF minimum**: Extra states(3) required for Broyden stability. Even with
   this, tutorial run converged to E=+302 Ha vs run_004's +288 Ha (different random seed
   in initial_guess). I/O pipeline validated; physics not representative.

4. **`pyenv activate` in non-login shells**: Does not initialise pyenv. Use direct python
   path: `/local/data/public/skcb2/pyenv/versions/3.10.19/envs/quantum-wave-packet/bin/python3`.

5. **`ColorBy(d, None)` in ParaView 6.1**: Fails with "invalid association string NONE".
   Fix: `d.ColorArrayName = ['POINTS', '']` for solid-color sphere rendering.

---

## Open items (not addressed in this task)

- Coronene SCF: reproducing run_004's +288 Ha ground state in the tutorial
- Publication-quality 2D matplotlib slices (`plots.py` still a stub) — see `inqview_publication_plots.md`
- ParaView background colour: VisRTX forces grey background regardless of `render_view.Background`
