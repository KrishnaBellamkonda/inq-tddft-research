# Reflection: coronene density-writer FFT-shift bug

I retrospectively think the diagnostic ladder for the cross-shaped artefact
in INQ density renders should have included, very early, a **coordinate-import
test** that shifts the coronene system rigidly along z by a known offset and
checks whether the rendered density tracks that offset. That single test
would have isolated the writer indexing bug from any physics question
within an afternoon.

## What the missing test would have looked like

1. Run `run_06`-style ground state with coronene atoms centred at z = 0.
2. Repeat the same SCF with all atoms displaced by a known z-shift, e.g.
   `+5 Bohr`, `+10 Bohr`, `+L_z/4`, and (as a corner case) `+L_z/2 - dz`.
3. Render the total density to .vti and read off where the molecular plane
   appears.
4. Decision rule:
   - If the rendered plane tracks the displacement linearly and stays
     centred on the metadata coordinate `z = +5, +10, +L/4, …`, the writer
     is healthy → the cross artefact must be a physics issue.
   - If the rendered plane stays at a fixed visual position (e.g. always
     at the corners) regardless of the shift, OR moves in a predictable
     but *wrong* direction (e.g. wraps when the shift exceeds half the
     cell), the bug is in how the writer maps array indices to Cartesian
     coordinates — and the fix is purely mechanical.

This is the standard input/output sanity check we already require for any
new utility under `.claude/rules/development-feedback-loop.md` (point 2,
*"For coordinate computation: print one point; compare against
hand-calculated grid coordinate"*). I should have applied it at the very
first sign of the artefact rather than after several runs of escalating
SCF tolerance, mixing-history depth, and pseudopotential variations.

## Why this would have shortened the loop

- The FFT-natural / metadata-origin mismatch in
  `inqkit::fields::density::total/orbital` is a writer-side bug. None of
  the tested SCF parameters could ever have removed it. Five diagnostic
  runs (`run_01..run_05`, days of compute) produced essentially the same
  rendered artefact because they all routed the (correct) array through
  the same (incorrect) writer.
- A z-shift sweep would have separated *position computation* (whose
  ground truth I can hand-calculate from the .xyz file) from *SCF physics*
  (whose ground truth requires Qball or a literature value). The former
  is a five-minute test; the latter is a multi-hour run.
- The same test, applied to any new field-extraction code I write
  (orbital wavefunctions, current density, dipole density, screen-plane
  slices), guards against the same class of bug being reintroduced
  silently.

## Going forward

For every new field writer or grid extractor I add to `inq-stack/`:

1. Write a synthetic-input test that injects a known peak at a known
   physical position (e.g. a Gaussian centred at `(+5, 0, 0)` Bohr).
2. Round-trip through the writer and reader.
3. Assert that the rendered peak position matches the input position to
   within one grid spacing.

This is now codified in
`inq-stack/tests/python/test_density_fft_shift_logic.py`, but the broader
lesson is the *coordinate-import sanity check* belongs in the diagnostic
ladder for any future "rendered output looks weird" report — before any
SCF-tolerance, mixing-scheme, or pseudopotential exploration.
