# Plan: Diagnose the cross-shaped density artefact in INQ ground states

## Context

INQ ground-state densities of planar π-systems (coronene, benzene, graphene,
quarter-coronene) appear as a cross-like feature streaming along ±x, ±y, ±z
beyond the molecular footprint, while the same systems run in Qball give clean,
physical densities. The user has previously **succeeded** at running coronene
with atoms centred on the z = 0 plane (matching INQ's [-L/2, +L/2] convention)
and the diffraction pattern came out correct; however the **rendered density
still looked split into quarters at the corners of the box**, and so did the
LEED pattern. This points away from a simulation bug and toward a
**density-writer / VTI-conversion indexing bug**: the array values may be
FFT-natural-ordered (origin at array index 0), while the VTI metadata claims
the array origin is at -L/2, so the visualisation shifts the molecule to the
edges of the box.

The user's request now is concrete:

> "I want to test by having a coronene.xyz file along the z = 0 plane. Then,
> I want to run ground-state convergence as a test and visualise the density
> of the orbitals and the total electronic density as .vti files for me to
> examine carefully."

So the immediate goal is a clean, well-instrumented diagnostic run that lets
the user inspect both the raw binary density and the converted .vti, and lets
us pinpoint whether the writer is producing an FFT-natural array under a
[-L/2, +L/2] origin metadata.

## Phase-1/2 findings relevant to the new scope

### Cell convention — INQ uses [-L/2, +L/2]

- `inq/src/systems/cell.hpp:212` — `contains(point)` requires fractional
  coords ∈ `[-0.5, 0.5)`.
- `inq/src/systems/cell.hpp:472–475` — unit test confirms `cell.contains(-5, 5, 5)`
  for cell `(28.62, 90.14, 12.31)`.
- `inq/src/systems/cell.hpp:219–225` — `position_in_cell` wraps to `[-0.5, 0.5)`.

### Internal array index — INQ uses FFT-natural ordering (this is the smoking gun)

- `inq/src/basis/grid.hpp:78–84` — `to_symmetric_range`: array index `i` maps
  to symmetric index via `if(i >= (sizes+1)/2) i -= sizes`. So **array index 0
  corresponds to symmetric index 0**, i.e. **the cell centre, not `-L/2`**.
- `inq/src/basis/grid.hpp:90–95` — `from_symmetric_range`: the inverse, adds
  `sizes` if the symmetric index is negative. So negative half of the
  symmetric range lives in the upper half of the array.
- `inq/src/basis/grid.hpp:109–115` — `symmetric_range_begin = -sizes/2`,
  `symmetric_range_end = sizes/2 + sizes%2`.

### Density writer mis-aligns array layout with VTI metadata

`inq-stack/include/inqkit/fields/density.hpp:43-86` (and the orbital twin at
:155–198) does:

```cpp
field.origin_x_bohr = basis.symmetric_range_begin(0) * spacing[0];   // = -L/2
...
for (int ix = 0; ix < nx; ++ix)
  for (int iy = 0; iy < ny; ++iy)
    for (int iz = 0; iz < nz; ++iz)
      field.values[flatten(ix,iy,iz)] = hc[ix][iy][iz];               // FFT-natural
```

The **origin** field is set to `-L/2` (matching the rendered VTI
`Origin="-17.4611 …"`), but the values are written in FFT-natural array order:
`hc[0][0][0]` is the value at the **cell centre**, not at `-L/2`. The .vti
reader (in `inq-stack/python/inqview/vti.py` and ParaView) then maps array
index 0 to position `Origin + 0*Spacing = -L/2`, placing the cell-centre
values at the cell corner. For an even-sized grid the FFT-natural layout has
`hc[i]` correspond to physical position `(i if i < (size+1)/2 else i-size) *
spacing`, so the array needs an `fftshift` of `(size+1)/2` along each axis
before the metadata-implied origin makes physical sense:

```
correct VTI value at index ix  =  hc[(ix + (size+1)/2) % size][...][...]
```

Right now the writer omits this shift, which would explain *exactly* the
reported behaviour: a centred molecule visualised as four quarters at the
corners of the (xy-slice) view, the same pattern in the LEED image, and so on.
The simulation itself can be entirely correct.

The `inq/src/observables/density.hpp` and `examples/h2.cpp` style code in INQ
that writes cube/xsf files internally should be examined for parity — they
may already do the fftshift correctly, which would prove that the bug is in
*inqkit*'s extractor, not in INQ.

## Hypotheses, post-redirect

| # | Hypothesis | Likelihood | Notes |
|---|---|---|---|
| **H1** | **Density writer mis-aligns FFT-natural array with [-L/2, +L/2] origin metadata** — the real visualisation bug. | **Very high** | Direct code inspection; matches user's "split into quarters at corners" exactly. |
| H2 | The SCF itself is wrong even with centred coordinates | Low | User reports earlier centred run gave a "right diffraction pattern". We must still verify that GS energy is sane (negative, ~-150 Ha for valence PBE). |
| H3 | Both H1 *and* a (smaller) physics issue (e.g. PSP mismatch, Broyden vs. Anderson, smearing) | Possible | Becomes relevant only if H1's fix doesn't completely clean the picture. |

H1 dominates. The plan therefore focuses on a single well-instrumented
diagnostic run that proves H1 is the cause, without yet modifying the
writer. The writer fix (and re-render of existing data) is queued as a
follow-up that requires an explicit user go-ahead, because it is a code
change to `inq-stack/include/inqkit/fields/density.hpp`.

---

# Recommended approach

## Step 1 — Single well-instrumented diagnostic run with centred coronene

Create:

- `Tutorial/coronene-leed/run_diagnoses/run_06_centred_writer_check/`
  - `coronene_centred.xyz`
  - `run.cpp`
  - `analysis.py`

### `coronene_centred.xyz`

Use the **exact** Qball geometry from `coronene-qball/coronene.sys`
(lines 14–50). The 36 atoms sit in the xy-plane at `z = 0`, with the
molecule centre at the origin. Header: `Coronene C24H12 D6h, atoms in INQ
[-L/2,+L/2] convention (Qball-parity), centred at z=0 plane`.

### `run.cpp`

Start from `run_01_tight_scf/run.cpp` and modify:

1. Geometry file → `coronene_centred.xyz`.
2. **Pre-SCF defensive assertion**: after `parse(...)`, iterate over atoms
   and abort with a clear error if any atom is outside `[-Lx/2, +Lx/2]`,
   `[-Ly/2, +Ly/2]`, `[-Lz/2, +Lz/2]`. Print the fractional position of the
   atom that's most distant from the origin. This catches any future xyz
   file written under the wrong convention.
3. **Post-SCF diagnostic dump** (the key new instrumentation): after the
   ground state converges, write to `results/grid_diagnostics.txt`:
   - `nx, ny, nz`
   - `dx, dy, dz` (in Bohr)
   - `Lx = nx*dx`, `Ly`, `Lz`
   - `symmetric_range_begin = (-Lx/2, -Ly/2, -Lz/2)`
   - For each axis, the physical position of the four corners of the array:
     index 0, index size/2 - 1, index size/2, index size - 1. This makes
     the FFT-natural index→position mapping explicit on the run.
   - Cubic (linear) probe of the total density along the x-axis at
     `y = 0, z = 0` for **every** array index, dumped as
     `(ix, hc[ix][0][0], physical_x_under_fft_natural,
     physical_x_under_naive)`. The user can read this file and see at a
     glance whether the density peak sits at array index 0 (FFT-natural,
     molecule centred at origin) or at array index nx/2 (naive, molecule
     centred at L/2).
   - The same probe along `y` and `z`.
   - The integrated total density: `sum(hc) * dx*dy*dz` should equal the
     number of electrons (108 for coronene). This is independent of the
     index ordering and is a clean SCF sanity check.
4. Keep the existing per-orbital write loop, but **also** write the **raw
   binary** density (which `RealField3DWriter` already does in `.raw +
   .meta.txt` form). The user can inspect both the raw binary (where the
   FFT-natural ordering is preserved) and the rendered `.vti` (where the
   ordering is reinterpreted under `Origin = -L/2`).
5. SCF options unchanged from run_01: PBE, cutoff 54 Ha, extra states 8,
   `energy_tolerance(1e-6_Ha)`, `max_steps(1000)`,
   `broyden_mixing().mixing_ndim(8).mixing(0.1)`. We are *not* changing the
   physics — we want to isolate the writer.

### `analysis.py`

Reuse `run_01/analysis.py` to convert `.raw + .meta.txt` to `.vti`. **Do
not modify** the converter yet — we want to see the buggy output that
visualises the molecule at the corners.

### How the user verifies

After the run:

1. **`results/ground_state_summary.txt`** — total energy: should be
   **negative**, of the order of -150 Ha (valence PBE for coronene with NC
   PSPs). Iteration count should be modest (≲100).
2. **`results/grid_diagnostics.txt`** — the x-axis probe of `hc[ix][0][0]`
   should show **density peaks at `ix` = 0** and at the nearby indices
   *and* at `ix` near `nx-1` (because the molecule wraps in the FFT-natural
   array). Equivalently, the centre of mass of `hc` along x is **near
   array index 0**, not near `nx/2`. If so, H1 is fully confirmed.
3. **`results/vti/density/density_t000000.vti`** — open in ParaView. If H1
   is correct: **the molecule will appear in four quarters at the corners
   of the box**, exactly as the user has previously reported. The xy-slice
   at the metadata centre (z = 0 in metadata = `nz/2` in array) will look
   like vacuum, and the slice at the metadata corner will show the
   molecule. **If this is what the user sees**, we have isolated the bug
   to the writer.
4. **`results/vti/orbital_density/orbital_NNNN.vti`** — same expectation
   for individual KS orbitals. The HOMO/HOMO-1 pair should show the same
   "split-into-quarters" behaviour, confirming the bug is per-field, not
   physics-specific.

## Step 2 — Once H1 is confirmed: propose the writer fix (separate plan)

If Step 1 confirms H1, the fix is small and surgical:

```cpp
// in inqkit::fields::density::total / orbital, where the loop is
auto fft_shift = [](int idx, int size) {
  return (idx + (size + 1) / 2) % size;     // moves -L/2 to array index 0
};
for (int ix = 0; ix < nx; ++ix) {
  int sx = fft_shift(ix, nx);
  for (int iy = 0; iy < ny; ++iy) {
    int sy = fft_shift(iy, ny);
    for (int iz = 0; iz < nz; ++iz) {
      int sz = fft_shift(iz, nz);
      auto flat = inqkit::detail::grid_layout::flatten_index(ix, iy, iz, ny, nz);
      field.values[flat] = hc[sx][sy][sz];
    }
  }
}
```

But this is a **separate code change** to `inq-stack/include/inqkit/fields/density.hpp`
and we should not make it before the user has examined Step 1's outputs.
Once the user gives the go-ahead, we open a new task that:

1. Adds the shift to `density::total` and `density::orbital`.
2. Adds a unit-style smoke test (synthetic Gaussian centred at origin →
   converted to .vti → reread → peak at index nx/2, ny/2, nz/2).
3. Re-renders existing diagnostic runs (or at minimum re-renders run_06's
   density) to confirm the molecule is now at the metadata centre.
4. Notes the fix in `docs/handovers/coronene_wp_scattering.md`.
5. Does NOT touch INQ source.

## Step 3 — Only if the SCF is *also* wrong (Step 1 H2 evidence)

If Step 1 shows a positive or wildly off total energy, then in addition to
the writer bug there is a physics issue. We then escalate through:

1. Pseudopotentials: pass Qball's `C_ONCV_PBE-1.2.xml` and
   `H_ONCV_PBE-1.0.xml` via `ionic::species("C").pseudo_file(...)`.
2. Mixing: try `linear_mixing().mixing(0.1)` (INQ has only LINEAR and
   BROYDEN — no Anderson).
3. Smearing: add `options::electrons{}.temperature(10.0_K)`.
4. Cutoff: bump to 100 Ha (most expensive — needs explicit user approval).

These are deferred until the writer hypothesis is settled.

## Critical files to be modified or created

- **NEW**: `Tutorial/coronene-leed/run_diagnoses/run_06_centred_writer_check/coronene_centred.xyz`
- **NEW**: `Tutorial/coronene-leed/run_diagnoses/run_06_centred_writer_check/run.cpp`
- **NEW**: `Tutorial/coronene-leed/run_diagnoses/run_06_centred_writer_check/analysis.py`
- **NO modifications** to INQ source.
- **NO modifications** to `inq-stack/` headers in this plan. The writer
  fix is queued for a follow-up plan once H1 is confirmed.

## Existing files / utilities to reuse

- `Tutorial/coronene-leed/run_diagnoses/run_01_tight_scf/run.cpp` — template
- `Tutorial/coronene-leed/run_diagnoses/run_01_tight_scf/analysis.py` — template
- `Tutorial/coronene-leed/run_diagnoses/coronene-qball/coronene.sys` — source
  of the centred geometry (lines 14–50)
- `Tutorial/coronene-gs-with-inqkit/coronene_centered.xyz` — possibly the
  same centred file from the user's earlier successful diffraction run; we
  should diff it against the Qball coronene.sys before adopting it.
- `inqkit::io::RealField3DWriter`, `inqkit::fields::density::total/orbital`
  — used as-is. We are diagnosing them, not patching them yet.
- `inq-run` wrapper — standard build/run path.

## Validation menu (per `.claude/rules/testing.md`)

- **Tier A** (always run on Step 1):
  - SCF convergence reached at 1e-6 Ha tolerance
  - No NaN/Inf in energies
  - Total energy negative and of physical order (−100 to −200 Ha)
  - `∫ ρ dV ≈ N_e` (= 108 for coronene)
- **Tier B** (visual checks the user does):
  - Examine `results/grid_diagnostics.txt` x/y/z probes — peak at array
    index 0 ⇒ FFT-natural ⇒ writer bug confirmed
  - Open `results/vti/density/density_t000000.vti` in ParaView — molecule
    at corners ⇒ writer bug confirmed
  - Open a frontier orbital `.vti` — same pattern ⇒ writer bug confirmed
- **Tier C** (deferred): full convergence study; spectrum comparison with
  Qball after a TDDFT kick. Not part of this plan.

## Assumptions and open questions

- **Assumption**: the user's earlier successful diffraction-pattern run
  (with centred coords) had a sane total energy. We will verify this in
  Step 1 by inspecting `results/ground_state_summary.txt` rather than
  relying on memory.
- **Assumption**: the FFT-natural index convention is consistent across
  even-sized grids (nx, ny, nz = 120, 120, 200 for our cell at 54 Ha
  cutoff). For odd sizes the shift formula `(idx + (size+1)/2) % size`
  remains correct, but we will only generate even-sized grids in this
  plan.
- **Question for the user (raised by ExitPlanMode)**: do you want me to
  also include in run_06's `analysis.py` an experimental
  fftshift-on-the-Python-side variant that reads the `.raw + .meta.txt`
  and writes a *second* `.vti` with the values pre-shifted? This would
  let you compare buggy vs. fixed renderings *without* changing the
  inqkit C++ writer in this plan. (Recommended: yes, it makes Step 1
  fully diagnostic.)

## What this plan does NOT do

- Modify INQ source.
- Modify the inqkit C++ density writer (Step 2 is a separate plan).
- Re-run the LEED TDDFT pipeline.
- Test pseudopotential, mixing, smearing, or cutoff hypotheses unless
  Step 1 reveals a physics issue beyond the writer indexing.
