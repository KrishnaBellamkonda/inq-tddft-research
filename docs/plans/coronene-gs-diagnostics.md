# Plan: Coronene Ground State Diagnostic Runs

**Branch:** `fixes/coronene-gs` (branched from `main` after merging `features/jellium-inital-exploration`)
**Directory:** `Tutorial/coronene-leed/run_diagnoses/`

---

## Context

A cross-like artifact (central bright spot + horizontal/vertical arms) appears in the coronene
LEED pattern across all 7 TDDFT runs, regardless of wavepacket energy, sigma, cell size, or
coronene orientation. Orbital density visualisation from `Tutorial/coronene-leed/run_07/`
(ground-state only, all 62 KS states written) shows unphysical density shapes, confirming
the error is in the SCF ground state — not in the wavepacket injection or propagation.

**XYZ parser verified clean:** `inq/src/parse/xyz.hpp:57` reads coordinates in Angstroms and
converts to Bohr via `atom_position * unit.in_atomic_units()`. No double-conversion risk in the
`ions::parse(filename, cell)` path. Parser is not the bug.

**Primary suspect:** SCF energy tolerance of `1e-4_Ha` is loose for a system with 54 occupied
states and D6h molecular symmetry. Coronene has doubly-degenerate HOMO and HOMO-1 (by D6h
symmetry). The steepest-descent eigensolver may resolve these into unphysical linear
combinations before the density reaches its symmetry-correct form — the total energy converges
before the individual orbital shapes do.

---

## Diagnostic Logic

| Run | What it isolates | Cross disappears implies |
|-----|-----------------|--------------------------|
| 01 tight_scf | SCF tolerance (1e-4 -> 1e-8 Ha) | Loose convergence was the bug |
| 02 benzene | Different D6h aromatic, far fewer states | Coronene-specific degeneracy issue |
| 03 coronene_hardcoded | Coords via C++ ions.insert, no xyz file | XYZ parser was the bug |
| 04 graphene | Extended honeycomb, ~31 C, no H passivation | Molecule identity matters |
| 05 quarter_coronene | 7C+3H fragment, breaks D6h entirely | Symmetry / degeneracy matters |

---

## SCF Parameters (uniform across all 5 runs)

All runs use tighter convergence than run_07:
- `energy_tolerance(1e-8_Ha)` (was `1e-4_Ha`)
- `max_steps(1000)` (was `300`)
- `broyden_mixing()`, `mixing_ndim(8)`, `mixing(0.1)` — unchanged from run_07
- `cutoff(54.0_Ha)`, `.pbe()`, `.finite()` cell — unchanged from run_07

Only `energy_tolerance` and `max_steps` change vs run_07, keeping the diagnostic clean.

---

## Cell (all 5 runs)

```cpp
static constexpr double LX_BOHR = 34.9222;   // 18.48 Ang
static constexpr double LY_BOHR = 34.9222;
static constexpr double LZ_BOHR = 59.9043;   // 31.7 Ang

auto cell = systems::cell::orthorhombic(
    LX_BOHR * 1.0_b, LY_BOHR * 1.0_b, LZ_BOHR * 1.0_b).finite();
```

Molecule centred at `(LX/2, LY/2, LZ/2)` = `(17.461, 17.461, 29.952)` bohr in all runs.

---

## Outputs (all 5 runs)

Every run writes:
- `results/density/` — total ground-state density (1 frame, via `inqkit::fields::density::total`)
- `results/orbital_density/orbital_XXXX/` — per-KS-state orbital density (all states)
- `results/orbital_density/orbital_index_map.csv`
- `results/checkpoint/` — INQ native checkpoint for later restart
- `results/ground_state_summary.txt` — human-readable parameter record

---

## Run 01 — `run_01_tight_scf/`

**Hypothesis tested:** SCF tolerance is the primary cause.

Geometry: copy `coronene_leed.xyz` from `Tutorial/coronene-leed/run_07/`. No new xyz needed.
System: C24H12, 108 electrons, 54 occupied + 8 extra = 62 states.

Only changes from run_07: SCF tolerance, max_steps, mixing_ndim (see table above).

---

## Run 02 — `run_02_benzene/`

**Hypothesis tested:** Cross is coronene-specific (due to its many near-degenerate occupied states).

Benzene C6H6 in the same finite cell. D6h symmetry like coronene but only 30 electrons, 15
occupied states. HOMO-LUMO gap ~6 eV (LDA) — much less risk of degenerate-orbital artefact.

**`benzene.xyz`:**
```
12
Benzene C6H6 D6h, centred at (9.2408,9.2408,15.8496) Ang  [Lx=Ly=18.48, Lz=31.7 Ang]
C   10.6378   9.2408  15.8496
C    9.9393  10.4506  15.8496
C    8.5423  10.4506  15.8496
C    7.8438   9.2408  15.8496
C    8.5423   8.0310  15.8496
C    9.9393   8.0310  15.8496
H   11.7236   9.2408  15.8496
H   10.4822  11.3908  15.8496
H    7.9994  11.3908  15.8496
H    6.7580   9.2408  15.8496
H    7.9994   7.0908  15.8496
H   10.4822   7.0908  15.8496
```
C-C = 1.397 Ang, C-H = 1.086 Ang. `extra_states(8)` -> 23 total states.

**User must visually approve `benzene.xyz` before build.**

---

## Run 03 — `run_03_coronene_hardcoded/`

**Hypothesis tested:** XYZ file parsing introduces a coordinate offset or misassignment.

Identical geometry to run_01 but all 36 atom positions inserted directly in C++ via
`ions.insert("C", {x_angstrom, y_angstrom, z_angstrom})` using the `_angstrom` magnitude
suffix. INQ's `add_atom()` calls `in_atomic_units()` on the quantity type, which converts
Angstroms to Bohr correctly (`inq/src/systems/ions.hpp:48`). No `.xyz` file at runtime.

If run_03 matches run_01 exactly: parser was innocent. If it differs: parser was introducing
an error.

No xyz to approve. 108 electrons, 54 occupied, `extra_states(8)`.

---

## Run 04 — `run_04_graphene/`

**Hypothesis tested:** Cross appears in any in-plane extended honeycomb system, not just coronene.

Small graphene nanoflake, no hydrogen passivation:
- Graphene lattice constant a = 2.46 Ang, C-C = 1.42 Ang
- Primitive vectors: a1 = (2.46, 0) Ang, a2 = (1.23, 2.131) Ang
- Basis: A at (0,0), B at (1.23, 0.711) per unit cell
- All sites within r < 5.0 Ang from origin selected (~31 C atoms total)
- Patch translated to cell centre (9.2408, 9.2408, 15.8496) Ang
- Stored in `graphene_nanoflake.xyz`

Bare graphene edges produce metallic-like edge states. Use `extra_states(12)` to
accommodate additional near-Fermi unoccupied states.

**User must visually approve `graphene_nanoflake.xyz` before build.**

---

## Run 05 — `run_05_quarter_coronene/`

**Hypothesis tested:** Cross depends on full D6h symmetry of the complete molecule.

Top-left quadrant of coronene (viewed from -z), atoms with x <= 9.2408 Ang AND y >= 9.2408 Ang.

**`quarter_coronene.xyz`:**
```
10
Coronene top-left quarter (7C+3H), centred fragment  [same cell as run_07]
C    8.530340  10.471462  15.849600
C    7.819840  11.702084  15.849600
C    7.819840   9.240840  15.849600
C    6.398840   9.240840  15.849600
C    8.530340  12.932706  15.849600
C    6.398840  11.702084  15.849600
C    5.688340  10.471462  15.849600
H    8.325105  13.999137  15.849600
H    5.577901  12.413038  15.849600
H    4.662166  10.826939  15.849600
```
Fragment has many dangling bonds — intentional for this diagnostic.
`extra_states(8)`, non-spin-polarized LDA.

**User must visually approve `quarter_coronene.xyz` before build.**

---

## File Layout

```
Tutorial/coronene-leed/run_diagnoses/
  run_all_diagnoses.sh           <- GPU-aware sequential/parallel launcher
  run_01_tight_scf/
    run.cpp
    coronene_leed.xyz            <- copy from run_07/
    analysis.py                  <- copy from run_07/
  run_02_benzene/
    run.cpp
    benzene.xyz
    analysis.py
  run_03_coronene_hardcoded/
    run.cpp                      <- no xyz file; all 36 atoms hardcoded
    analysis.py
  run_04_graphene/
    run.cpp
    graphene_nanoflake.xyz
    analysis.py
  run_05_quarter_coronene/
    run.cpp
    quarter_coronene.xyz
    analysis.py
```

---

## Runner Script — `run_all_diagnoses.sh`

GPU-aware launcher with at most 2 simultaneous simulations:

```bash
RUNS=(run_01_tight_scf run_02_benzene run_03_coronene_hardcoded
      run_04_graphene run_05_quarter_coronene)
MAX_PARALLEL=2

# At startup: query nvidia-smi for free GPU count; warn if < MAX_PARALLEL.
# Loop: while runs remain, count running background jobs; if < MAX_PARALLEL, launch next.
# Each run: cd <dir> && inq-run > run.log 2>&1 && python analysis.py >> run.log 2>&1
# Use `wait -n` (bash 4.3+) to wait for any one background job to free a slot.
```

Uses `nvidia-smi` to determine available GPUs at startup. Logs each run to `<dir>/run.log`.
analysis.py runs immediately after each simulation to produce VTI orbital densities.

---

## Phased Implementation

### Phase 1 (write all files)
1. Commit + merge `features/jellium-inital-exploration` to `main`
2. Create `fixes/coronene-gs` branch
3. Write all 5 `run.cpp` + all xyz files + `analysis.py` copies + runner script
4. User reads all files and visually inspects xyz geometries in their viewer

### Phase 2 (after user approval of xyz files)
1. User approves benzene.xyz, graphene_nanoflake.xyz, quarter_coronene.xyz
2. Execute `run_all_diagnoses.sh` (up to 2 runs in parallel on available GPUs)
3. User inspects orbital density VTI files in ParaView

---

## Verification

1. SCF convergence: each run must converge to `1e-8_Ha` within 1000 steps — check `run.log`
2. Orbital count: `results/orbital_density/` must have one directory per KS state
3. ParaView inspection: run_01 orbital shapes should show D6h symmetry if hypothesis is correct
4. Cross-check: run_03 densities should be numerically identical to run_01 (parser validation)
5. Update `docs/handovers/coronene_wp_scattering.md` after each phase

---

## Commits

**Commit 1** (on `features/jellium-inital-exploration`):
- Files: analysis.py (7 runs), inq-stack/python/inqview/screens.py
- Message: `Update jellium wp-rt analysis scripts and inqview screens module`

**Commit 2** (on `fixes/coronene-gs`):
- All of `Tutorial/coronene-leed/run_diagnoses/`
- Message: `Add 5 ground-state diagnostic runs to identify coronene LEED cross artifact`
