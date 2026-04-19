# Plan: inqview publication-quality 2D slice plots

*Created: 2026-04-19*

---

## Goal

Replace the current rough `imshow` slices in `analysis.py` with paper-ready 2D figures
that have physical axes, correct colormaps, contour lines, and optional atom-position
overlays. Each figure conveys exactly one quantity (density, orbital density, or
wavefunction real part) and is ready to drop into a manuscript without post-processing.

The work is entirely in `inq-stack/python/inqview/plots.py` (currently a stub) and
the `analysis.py` scripts that call it. No new modules; no new dependencies beyond
numpy and matplotlib.

---

## Visual specification

### All slice figures

| Property | Value |
|---|---|
| Font family | serif (matches most journal templates) |
| Font size (labels/ticks) | 9 pt |
| Font size (colorbar label) | 9 pt |
| Figure width | 3.5 in (single column) or 7.0 in (double column) |
| Figure height | set automatically from physical aspect ratio |
| DPI | 300 (matches `PlotDefaults.dpi` override for paper) |
| Axis units | Å (ångström) — conversion from bohr applied internally |
| Axis labels | "x (Å)", "y (Å)", "z (Å)" as appropriate for the slice |
| Colorbar | right side, same height as axes, label includes units |
| Colorbar units | e bohr⁻³ (density) · bohr⁻³ (orbital density) · bohr⁻³/² (wavefunction) |
| Contour lines | overlaid in semi-transparent black, see per-type rules below |
| Atom markers | open circles sized to covalent radius, coloured by element (C=grey, H=white edge) |
| Atom tolerance | atom shown on slice if its perpendicular distance < 1 bohr from the plane |
| Output format | PNG at 300 dpi (per project rule) |

### Per-field-type colourmap and normalisation

| Field type | Colourmap | Normalisation | Contour rule |
|---|---|---|---|
| Total density `n(r)` | `cividis` (from `PlotDefaults.scalar_cmap`) | `LogNorm(vmin=1e-4 × peak, vmax=peak)` | 6 log-spaced contours, α=0.4 |
| Orbital density `\|ψ\|²` | `cividis` | `LogNorm(vmin=1e-4 × peak, vmax=peak)` | 6 log-spaced contours, α=0.4 |
| Wavefunction Re(ψ) | `coolwarm` (from `PlotDefaults.signed_cmap`) | `TwoSlopeNorm(vcenter=0)` symmetric | single zero contour in black, lw=0.8 |

`LogNorm` makes low-density tails visible alongside the density peak — essential for
showing both the molecular core density and the long-range tails in one figure.

---

## API — `inq-stack/python/inqview/plots.py`

### Public functions

```python
def set_publication_rc() -> None:
    """Set matplotlib rcParams for single-column publication figures."""
    # sets: font.family=serif, font.size=9, axes.linewidth=0.8,
    #       xtick.major.width=0.8, savefig.dpi=300, figure.dpi=300

def density_slice(
    field: RealField3D,
    plane: Literal["xy", "xz", "yz"],
    index: int | None = None,          # grid index of slice; None = midplane
    *,
    atom_positions: np.ndarray | None = None,   # shape (N, 3), bohr
    atom_symbols: list[str] | None = None,
    units: str = "angstrom",           # "angstrom" or "bohr"
    figwidth_in: float = 3.5,
    log_norm: bool = True,
    contours: int = 6,
    title: str | None = None,
) -> tuple[Figure, Axes]:
    """
    Single publication-quality 2D slice of a real scalar field (density or orbital density).
    Returns (fig, ax) for further annotation before saving.
    """

def orbital_slice(
    field: ComplexField3D,
    plane: Literal["xy", "xz", "yz"],
    index: int | None = None,
    *,
    component: Literal["real", "imag", "magnitude"] = "real",
    atom_positions: np.ndarray | None = None,
    atom_symbols: list[str] | None = None,
    units: str = "angstrom",
    figwidth_in: float = 3.5,
    title: str | None = None,
) -> tuple[Figure, Axes]:
    """
    Single publication-quality 2D slice of a complex orbital field.
    For 'real' and 'imag': diverging coolwarm map, zero contour.
    For 'magnitude': treated identically to density_slice (LogNorm, cividis).
    Returns (fig, ax).
    """
```

### Internal helpers (private, not exported)

```python
_bohr_to_ang = 0.529177210903

def _slice_extent(meta, plane, index)
    # returns (horiz_coords, vert_coords, data_2d) in requested units

def _atom_overlay(ax, atom_positions, atom_symbols, plane, index, plane_coord, tol_bohr, units)
    # projects atoms within tol onto the slice; draws open circles

def _make_colorbar(ax, im, label)
    # attaches a neat colorbar with correct label
```

---

## Changes to `analysis.py` in both tutorials

Replace the current `_save_slices` helper with calls to `plots.density_slice` and
`plots.orbital_slice`. The atom positions are loaded from the xyz file for coronene;
for HF they are computed from the ion coordinates in `run.cpp`.

```python
# Load atom positions from xyz (coronene) or inline (HF)
# atom_positions: np.ndarray shape (N, 3) in bohr
# atom_symbols:   list[str]

set_publication_rc()

# Total density — slice at z = coronene plane (or z = midplane for HF)
fig, ax = density_slice(density, plane="xy", index=None,
                        atom_positions=atom_pos, atom_symbols=atom_sym,
                        title="Total electron density")
fig.savefig(VISDIR / "density_xy.png", bbox_inches="tight")

# HOMO orbital density
fig, ax = density_slice(rho_homo, plane="xy", index=None,
                        atom_positions=atom_pos, atom_symbols=atom_sym,
                        title="HOMO orbital density")
fig.savefig(VISDIR / "homo_density_xy.png", bbox_inches="tight")

# HOMO wavefunction — real part
fig, ax = orbital_slice(psi_homo, plane="xy", index=None,
                        component="real",
                        atom_positions=atom_pos, atom_symbols=atom_sym,
                        title=r"HOMO $\mathrm{Re}(\psi)$")
fig.savefig(VISDIR / "homo_psi_re_xy.png", bbox_inches="tight")
```

Each call produces one figure with one quantity — total density, orbital density, or
wavefunction — with axes in Å, a labelled colorbar, contour lines, and atom positions.

---

## Atom position loading

For `coronene-gs-with-inqkit/analysis.py`: parse `coronene_centered.xyz` directly in
Python (simple line-by-line reader — no external dependency needed).

For `hf-gs-with-inqkit/analysis.py`: hardcode from `run.cpp` values (L=16 bohr,
bond=0.917 Å centred at L/2):

```python
L_bohr = 16.0
half_bond_bohr = 0.917 * 1.8897259886 / 2
cx = L_bohr / 2
atom_positions = np.array([[cx, cx, cx - half_bond_bohr],   # H
                            [cx, cx, cx + half_bond_bohr]])  # F
atom_symbols = ["H", "F"]
```

---

## Element colour and size table (internal to plots.py)

| Symbol | Fill colour | Edge colour | Marker radius (pts) |
|---|---|---|---|
| H | white | black | 4 |
| C | #404040 | #404040 | 5 |
| N | #3050F8 | #3050F8 | 5 |
| O | #FF0D0D | #FF0D0D | 5 |
| F | #90E050 | #90E050 | 5 |
| default | #AAAAAA | #AAAAAA | 4 |

---

## 3D ParaView: atom spheres (`paraview.py`)

The ParaView pipeline is extended to render atoms as CPK-coloured spheres alongside
the volume density. A new `AtomSpec` dataclass carries the atom data; it is optional
on all render calls so existing code is unaffected.

### `AtomSpec`

```python
@dataclass
class AtomSpec:
    positions: list[list[float]]   # Cartesian positions in bohr (same frame as VTI)
    symbols:   list[str]           # element symbols, length N
    radius_scale: float = 0.4     # fraction of VDW radius used for sphere size
    opacity:      float = 1.0
    specular:     float = 0.3
    specular_power: float = 20.0
```

### CPK colours and VDW radii (built into `paraview.py`)

Standard CPK colours (H=white, C=grey, N=blue, O=red, F=green, …) and VDW radii in Å
are stored as module-level dicts and passed to pvbatch via the JSON config.

### pvbatch sphere rendering

After the volume pipeline, for each atom a `Sphere` source is created with
`Center = position_bohr`, `Radius = vdw_radius_bohr × radius_scale`,
`ThetaResolution = PhiResolution = 24`. `DiffuseColor` is set to the CPK colour.
Atoms are rendered once before the frame loop (ionic positions are fixed in GS).

### API change

```python
render_density_from_meta_series(..., atoms: AtomSpec | None = None) -> list[Path]
render_vti_series(...,             atoms: AtomSpec | None = None) -> list[Path]
```

## 2D matplotlib: linear normalisation flag

`density_slice` and `orbital_slice` accept `log_norm: bool = True`. When `False`,
a simple `Normalize(vmin=0, vmax=peak)` is used instead of `LogNorm`.

## Out of scope

- Interactive figures (plotly, ipywidgets)
- Any changes to the C++ inqkit side

---

## Validation

Before declaring complete, the following checks must pass for both tutorials:

1. `density_slice` output: colorbar range spans at least 3 orders of magnitude; atom
   markers appear at physically correct positions in the coronene x-y slice.
2. `orbital_slice(component="real")` output: coolwarm map is symmetric about zero;
   zero contour is visible as a clean line.
3. Both functions run without error when `atom_positions=None` (no overlay).
4. All output PNGs are ≥ 300 dpi and saved to `results/visualisation/`.
5. HF: density clearly concentrated on the F atom with a smaller lobe on H.
6. Coronene: density shows D₆h ring symmetry in the x-y midplane.

---

## Execution order

1. Implement `set_publication_rc`, `density_slice`, `orbital_slice` in `plots.py`
2. Test against HF (small, fast to re-run analysis.py)
3. Update `analysis.py` for both tutorials
4. Run both analysis scripts, visually inspect output
5. Record validation results in handover
