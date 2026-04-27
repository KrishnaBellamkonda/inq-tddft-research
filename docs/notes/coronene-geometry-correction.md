# Coronene geometry correction — `z = L/2` → `z = 0`

## Why this note exists

Every full-RT-LEED coronene scattering simulation produced before the
replication framework used a `coronene.xyz` whose carbon plane sat at
`z = L_z/2` instead of at `z = 0`. The convention conflict between the
`xyz` file and INQ's centred-cell convention caused the molecule to be
silently mis-placed, and the visible LEED pattern to render as a
four-corner-split ghost. This note explains the bug at the level a
new reader can reconstruct the failure mode and the fix from scratch.

## INQ's cell convention

INQ stores orthorhombic cells as **centred** intervals on every axis:

```
x ∈ [-L_x / 2, +L_x / 2]
y ∈ [-L_y / 2, +L_y / 2]
z ∈ [-L_z / 2, +L_z / 2]
```

Source: `inq/src/systems/cell.hpp:212` —
`bool contains(point) { return p[i] >= -0.5 && p[i] < 0.5 ... }` in
fractional coordinates, i.e. the cell occupies `[-L/2, +L/2]` in real
space when reduced to a single image. This is what the
`tsubonoya_2014_coronene.hpp:15-17` comment records:

> "INQ uses [-L/2, +L/2] for orthorhombic cells. The flake therefore
>  sits at z=0, matching coronene_centred.xyz (Qball-parity geometry)."

So an atom at z = 0 sits at the **mid-plane** of the cell. An atom at
z = +L/2 sits at the **+z face** — the cell boundary.

## What the buggy xyz files did

The legacy coronene xyz files (e.g. those used by
`04_leed_simulation/run_001..005`, `coronene-wp-rt/run_01..06`) placed
the entire molecule at z = `L_z/2`. Read directly:

```
36
coronene C24H12
C    1.421000   0.000000   15.849600
C    2.842000   0.000000   15.849600
...                          ^^^^^^^
```

Combined with `cell.orthorhombic(L_x, L_y, L_z).finite()` and
`ions::parse(xyz_file, cell)`, INQ's parser reads positions verbatim
(`inq/src/parse/xyz.hpp:57`) and places them into the cell as-is.
Atom z = +L_z/2 = +15.85 Å is *just inside* the +z face. The molecule
is at the cell boundary, not the cell centre.

That alone is mostly cosmetic — but it interacts with the INQ FFT
grid in a way that breaks the LEED visualisation, as the next section
shows.

## Pedagogical example: 1D, 4 grid points

Take a 1D cell of length `L = 4` with 4 grid points. INQ stores the
grid in **FFT-natural** order:

```
array index :  0      1      2      3
physical z  :  0     +1     -2     -1     (cell centre at index 0)
```

That is, INQ's `to_symmetric_range` maps array index `i` to the physical
coordinate `i*dz` if `i ≤ N/2`, else `(i - N) * dz`. So **array index 0
maps to physical z = 0**, the **cell centre**, and the wrap to the
negative half happens above index `N/2`.

Now place a δ-function at z = +2 = +L/2. In array-index space that's
index 2 — the cell-edge slot. After an FFT, this δ produces the
phase factor `exp(i k · L/2) = exp(i π n)` for the n-th harmonic, i.e.
alternating signs. When the inverse FFT is plotted in raw array order,
the peak appears split between index `N/2` (physical −L/2) and
indices `0` (physical 0) — the famous four-corner ghost.

Now place the same δ at z = 0. It lands at array index 0. No phase
oscillation, no wrap, no ghost. The peak appears cleanly at the cell
centre when plotted in physical coordinates.

## What changed in the corrected `coronene.xyz`

Every atomic z-coordinate was shifted by `−L_z/2` so the molecule sits
at z = 0 in INQ's centred frame. The corrected file is at:

```
ResearchProject/systems/coronene/shared/geometry/coronene.xyz
```

First rows:

```
36
Coronene C24H12, D6h symmetry, C-C=1.421 Ang, C-H=1.086 Ang
C    1.421000   0.000000   0.000000
C    2.842000   0.000000   0.000000
...
```

(Compare to the buggy version: same `x, y`; the `z = 15.8496` column is
now `z = 0`.) Nothing else in the xyz format changed.

This single change makes the file cell-size-agnostic: any orthorhombic
cell whose half-extent contains the molecule's footprint plus the WP
launch height plus a few σ — `L_z ≥ 2(b + ~3σ) ≈ 30 Bohr` minimum —
loads the same xyz unchanged. That is why the new framework reuses
one canonical xyz across cells of `60 Bohr`, `80 Bohr`, and `40 Bohr`.

## The LEED four-corner-split symptom

Even with the corrected geometry, the LEED pattern still appeared
split into four corners when plotted with the legacy Python loader.
This was a **separate** bug, in the visualisation layer, with the same
root cause as above: `LeedPatternAccumulator` writes its `.dat` files
in INQ's FFT-natural order (array index `(0, 0)` = physical origin,
`x = 0, y = 0`), but `inqview.load_leed_pattern` plotted the data
without an `np.fft.fftshift`, so the diffraction peak (at physical
origin) landed at array `(0, 0)` — a corner of the matplotlib image —
with its 4-way symmetric tail distributed to the other three corners.

The fix, ported from
`ResearchProject/systems/coronene/run_propagate_paper_replica/analysis.py`
(`_load_screen_centred`), is a single line in
`inq-stack/python/inqview/screens.py`:

```python
data = np.fft.fftshift(data)
origin_x_bohr = -0.5 * nx * dx_bohr
origin_y_bohr = -0.5 * ny * dy_bohr
```

After shift, array index `(0, 0)` maps to physical `(-L_x/2, -L_y/2)`
and the diffraction peak naturally lands near the image centre.
`LeedPattern.extent_bohr` then auto-spans `[-L_x/2, +L_x/2,
-L_y/2, +L_y/2]`. The companion `coordinate_checks/` plots in
`results/analysis/screens/coordinate_checks/` show both the raw
FFT-natural view and the centred view side by side, exactly as the
spec §17.6 requires.

## Files touched by the geometry correction

- New canonical xyz:
  `ResearchProject/systems/coronene/shared/geometry/coronene.xyz`
- C++ runs that load it: every `run_*/run.cpp` and every
  `save_gs/<sig>/run.cpp` under
  `ResearchProject/systems/coronene/`.
- Variant headers under `shared/configs/` set `WP_CZ_BOHR = +b`
  (positive Bohr above the molecule plane) on a centred cell.
- `inq-stack/python/inqview/screens.py` — fftshift + extent override
  (the LEED visualisation fix).
- `inq-stack/python/inqview/postprocess/screens.py` — coordinate-check
  raw plot now `np.fft.ifftshift`s the (already-shifted) data so the
  raw-vs-mapped comparison is meaningful again.

The legacy buggy run trees (`04_leed_simulation/`, `coronene-wp-rt/`,
`run_propagate_paper_replica/` predecessors) were moved to a
`legacy/` subfolder on disk per the cleanup TODO and are not part of
the new framework.
