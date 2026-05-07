# Plotting and Visualisation Rules

These rules apply to all Python analysis, plotting, GIF, and movie-generation scripts for the project.

The aim is not to rigidly template every plot, but to ensure that all visual outputs are scientifically readable, consistent across runs, and suitable for direct use in reports, slides, and diagnosis.

---

## 1. Numerical formatting

All displayed numerical values must be rounded to **3 significant figures**.

This applies to:

- plot titles
- legends
- colour-bar labels
- annotations
- timestep/time labels in GIFs
- printed summary values
- diagnostic text written onto figures

Examples:

```text
time = 2.36 fs
energy = 0.0184 Ha
frequency = 0.315 a.u.
sigma = 1.89 bohr
```

Do not over-print long floating point values such as:

```text
time = 2.356728193847 fs
```

---

## 2. General line, scatter, and time-series plots

These rules apply to line plots, scatter plots, spectra, current plots, energy plots, dipole plots, conservation checks, and similar figures.

### Axis labels

Every plot must have clearly labelled x and y axes.

Use appropriate physical units:

| Quantity | Preferred unit |
|---|---|
| distance / position | bohr |
| energy | Ha, or eV when more interpretable |
| time | fs |
| frequency | atomic units unless explicitly converted |
| current | atomic units |
| dipole | atomic units |
| overlap | dimensionless |
| density | bohr^-3 |
| LEED/screen intensity | arbitrary units, unless otherwise normalised |

Examples:

```text
Time / fs
Energy / Ha
Current J_z / a.u.
Frequency / a.u.
Position z / bohr
Overlap |<psi_i^GS | psi_j(t)>|^2
```

### Titles

Every plot title must include either:

1. the **run name**, or
2. the **hypothesis name**,

alongside the specific plot description.

Examples:

```text
run_05_d20: Total energy vs time
Coronene LEED cross hypothesis: screen_10 log intensity
Jellium WP spreading: sigma(t)
```

### Legends

Include a legend whenever more than one curve, marker set, component, or dataset is shown.

Legends are required for:

- energy component comparisons
- current component comparisons
- dipole component comparisons
- multiple runs on one plot
- raw vs filtered comparisons
- fitted background vs measured data
- FFT component comparisons

A legend is not required for a single unambiguous curve.

---

## 3. Density plots and density GIFs

These rules apply to total density, system density, wavepacket density, orbital density, and any density-derived slices or movies.

### Axes

All density plots must label both axes and include units.

Examples:

```text
x / bohr
y / bohr
z / bohr
```

For 2D slices, the axis labels must match the slice plane:

| Slice | x-axis | y-axis |
|---|---|---|
| xy | x / bohr | y / bohr |
| xz | x / bohr | z / bohr |
| yz | y / bohr | z / bohr |

### Titles for timestepped density plots

Every frame in a density GIF or movie must include:

- run name or hypothesis name
- density type
- slice plane
- timestep in the format `step 004/600`
- physical time, rounded to 3 significant figures

Example:

```text
run_03_d10: WP density, xz slice, step 004/600, t = 0.968 fs
```

### Fixed colour scale through time

All GIFs or movies made from timestepped density data must use the **same colour scale for every frame**.

This is essential. Do not autoscale each frame independently.

For a given animation, determine the global colour range before rendering frames:

```text
vmin = global minimum over all frames used in the GIF
vmax = global maximum over all frames used in the GIF
```

or, if using percentile clipping for readability, use the same percentile range for every frame:

```text
vmin = global 1st percentile over all frames
vmax = global 99th percentile over all frames
```

The chosen colour scaling policy should be consistent across frames and, where possible, mentioned in the plot or analysis notes.

---

## 4. LEED and screen plots

These rules apply to raw LEED screens, instantaneous screens, time-windowed screens, filtered screens, log plots, and LEED GIFs.

### Coordinate mapping must be physically correct

The Python LEED visualisation must correctly map screen coordinates to image axes.

A known failure mode is producing a diffraction pattern split into four pieces and placed at the image edges. This usually indicates an incorrect coordinate convention, wrap-around mistake, or wrong interpretation of the screen grid ordering.

Before trusting LEED plots, check:

- the screen coordinate origin
- whether coordinates are centred or indexed from zero
- whether `fftshift`/periodic wrapping has been applied incorrectly
- whether x/y axes are transposed
- whether the raw `.dat` file stores coordinates explicitly or only stores values
- whether the plotting code is using the same ordering as the writer

If a LEED pattern appears split into four edge/corner pieces, the visualisation should be treated as incorrect until the coordinate mapping is fixed.

### Axes and labels

All LEED plots must label screen axes with units, usually:

```text
x / bohr
y / bohr
```

or, if reciprocal-space plotting is explicitly used:

```text
k_x / bohr^-1
k_y / bohr^-1
```

### Titles

LEED titles must include:

- run name or hypothesis name
- screen index
- screen type: total, instantaneous, or time-windowed
- timestep or timestep range where relevant
- time or time range where relevant

Examples:

```text
run_05_d20: total LEED screen_10
run_05_d20: instantaneous screen_10, step 120/600, t = 2.9 fs
run_05_d20: windowed screen_10, steps 100-180, t = 2.42-4.36 fs
```

### Fixed colour scale through time

All LEED GIFs or movies made from timestepped screen data must use the same colour scale through time.

This applies to:

- instantaneous LEED GIFs
- time-windowed LEED GIFs
- filtered LEED GIFs
- log-scale LEED GIFs

Do not autoscale frame-by-frame.

---

## 5. Overlap visualisation

Overlap visualisation should represent evolved KS orbitals in terms of the ground-state KS orbital basis.

The central quantity is:

```text
O_ij(t) = |<psi_i^GS | psi_j(t)>|^2
```

where:

- `i` indexes ground-state reference orbitals
- `j` indexes evolved orbitals at time `t`
- the wavepacket state is treated as one of the evolved states when present

### KS orbital overlap plots

Overlap plots for the evolved KS orbital manifold must be shown as **heatmaps**.

Use matrix-style heatmaps where:

```text
rows    = ground-state KS orbitals
columns = evolved KS orbitals
colour  = overlap magnitude
```

Recommended axis labels:

```text
Ground-state KS orbital index i
Evolved KS orbital index j
Overlap |<psi_i^GS | psi_j(t)>|^2
```

For time-dependent overlap movies/GIFs, each frame must include:

- run name or hypothesis name
- timestep in the format `step 004/600`
- time rounded to 3 significant figures
- fixed colour scale through time

The colour scale should usually be fixed over `[0, 1]` for overlap heatmaps unless there is a strong reason to zoom in.

### Wavepacket overlap plots

For the wavepacket specifically, the preferred visualisation is an animated **bar graph** over time.

Use:

```text
x-axis = ground-state KS orbital index i
y-axis = |<psi_i^GS | psi_WP(t)>|^2
```

The title should include:

- run name or hypothesis name
- `WP overlap with GS KS orbitals`
- timestep as `step 004/600`
- physical time rounded to 3 significant figures

Example:

```text
run_03_d10: WP overlap with GS KS orbitals, step 004/600, t = 0.968 fs
```

The y-axis scale must be fixed through the animation, preferably:

```text
0 <= overlap <= 1
```

unless a zoomed-in diagnostic plot is explicitly created and labelled as such.

---

## 6. GIF and movie rules

These rules apply to every animated visualisation.

### Required frame title information

Every frame must include:

- run name or hypothesis name
- data type being visualised
- timestep as `step current/total`
- physical time rounded to 3 significant figures

### Fixed colour or y-axis scaling

All GIFs made from timestepped data must use consistent scaling through time.

This includes:

- density GIFs
- LEED GIFs
- filtered LEED GIFs
- overlap heatmap GIFs
- wavepacket-overlap bar-chart GIFs
- orbital-density GIFs
- current/energy animated diagnostics if any are made

For heatmaps and density plots, fix the colour scale.

For animated bar charts, fix the y-axis range.

For animated line plots, fix both x and y limits unless the animation is explicitly designed to show progressive reveal.

---

## 7. File categorisation expectations

Raw numerical data should stay under:

```text
results/raw/
```

Analysis figures, plots, GIFs, and videos should stay under:

```text
results/analysis/
```

For observables:

```text
results/raw/observables/
```

contains CSV or numerical output.

```text
results/analysis/observables/
```

contains PNGs, GIFs, or other visual outputs.

Do not create a separate `derived_data/` folder unless explicitly requested. Derived numerical CSVs, such as FFT spectra, should be categorised with the observable they came from under `raw/observables/`.

---

## 8. Final checklist before saving a figure

Before saving any plot, GIF, or movie, check:

- [ ] Are numerical values rounded to 3 significant figures?
- [ ] Are axes labelled?
- [ ] Are units included?
- [ ] Does the title include the run name or hypothesis name?
- [ ] Is there a legend where multiple datasets/components are plotted?
- [ ] For timestepped data, does the title include `step current/total`?
- [ ] For timestepped data, does the title include physical time?
- [ ] For GIFs/movies, is the colour scale or y-axis range fixed through time?
- [ ] For LEED plots, has the coordinate mapping been checked?
- [ ] For overlap plots, are KS overlaps shown as heatmaps and WP overlaps as animated bar graphs?