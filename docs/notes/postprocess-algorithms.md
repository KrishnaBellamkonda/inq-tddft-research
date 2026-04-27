# Postprocess algorithms — equations and implementations

Reference for every numerical and physical algorithm used by the coronene
RT-LEED framework's postprocess pipeline. Each section gives the
mathematical definition, the discrete form actually computed, and the
file:line where the code lives.

---

## 1. Discrete Fourier transform of an observable time series

**Definition.** Given a signal `x(t)` sampled on a uniform grid
`t_n = n·dt` for `n = 0, …, N − 1`, its discrete Fourier transform is

```
X[k] = Σ_{n=0}^{N-1}  x[n] · w[n] · exp(-2πi · k · n / N),
                k = 0, …, N − 1.
```

The frequency axis (in atomic units of inverse time) is

```
ω_k = 2π · k / (N · dt),    k = 0, …, N/2.
```

We take `|X[k]|` as the amplitude spectrum. A window `w[n]` reduces
spectral leakage; the default is **Hann**:

```
w[n] = 0.5 · (1 − cos(2π · n / (N − 1))).
```

The signal is detrended before windowing: `x[n] ← x[n] − mean(x)`,
so DC is removed and the energy spectrum's zero-frequency bin
doesn't dominate.

**Implementation.** `inq-stack/python/inqview/fourier.py`,
`FourierTransform.transform_column(df, col)`:

```python
y = (df[col] - df[col].mean()).to_numpy()
y *= self.window(N)  # Hann by default
F = np.fft.rfft(y)
freq = np.fft.rfftfreq(N, d=dt) * 2*np.pi  # angular frequency in 1/a.u.
return FourierResult(freq, np.abs(F) / N, ...)
```

The dipole power spectrum is the squared magnitude (∝ optical
absorption strength via Fermi's golden rule).

---

## 2. Orbital overlap integral

**Definition.** Given the ground-state KS orbital `ψ_i^{GS}(r)` and a
real-time evolved orbital `ψ_j(r, t)`, the time-dependent overlap is

```
S_ij(t) = ⟨ψ_i^{GS} | ψ_j(t)⟩ = ∫ ψ_i^{GS}(r)* · ψ_j(r, t) d³r.
```

The squared modulus,

```
O_ij(t) = |S_ij(t)|²,
```

is the probability that the evolved orbital `j` is found in the
GS-orbital subspace `i` at time `t`. For a wave-packet state injected
into the last extra slot, the row `O_{i, wp}(t)` for `i ∈ [0, n_occupied)`
quantifies how much of the WP has projected onto the molecule's bound
KS manifold.

**Discrete form.** With orbitals sampled on a uniform 3D grid with
volume element `dV = dx · dy · dz`:

```
S_ij(t) ≈ dV · Σ_r  conj(ψ_i^{GS}[r]) · ψ_j[r, t]
```

where `r = (i_x, i_y, i_z)` runs over all grid points.

**WP-only optimisation.** Computing the full `n_ref × n_evolved` matrix
costs `O(n_states · n_occupied · n_grid)` per step. The WP-only variant
computes only the row containing the WP state — `O(n_occupied · n_grid)` —
which for the coronene base run (n_states = 62, n_occupied = 54) is
about 62× cheaper per step. The plan §4.6 item 9 mandates the WP-only
form for production runs.

**Implementation.** `inq-stack/include/inqkit/observables/orbital_overlap.hpp:90`,
`OrbitalOverlapMatrix::snapshot_wp_only(electrons, time_au, step)`:

```cpp
auto wp = fields::orbital::wavefunction(electrons, n_ref_);
std::vector<double> O(n_ref_, 0.0);
for (int i = 0; i < n_ref_; ++i) {
    std::complex<double> inner(0.0, 0.0);
    auto const& ri = ref_wfns_[i].values;
    auto const& ej = wp.values;
    for (std::size_t r = 0; r < n_pts; ++r)
        inner += std::conj(ri[r]) * ej[r];
    inner *= dv_;
    O[i] = std::norm(inner);  // |inner|^2
}
```

`ref_wfns_` are extracted once from the post-injection electrons object
and held in CPU memory; only the WP wavefunction is pulled fresh from
the device per step.

---

## 3. Modified Gram–Schmidt re-orthogonalisation at WP injection

**Definition.** After raw injection of the Gaussian WP into the last
extra-state slot, the WP can have a small projection onto each occupied
orbital `i`. Modified Gram–Schmidt removes it iteratively:

```
ψ_wp ← ψ_wp − Σ_{i = 0}^{n_occupied − 1} ⟨ψ_i^{GS} | ψ_wp⟩ · ψ_i^{GS}
```

Note that the projections are recomputed after each subtraction (the
"modified" in MGS), avoiding the loss of orthogonality the classical
form suffers when overlaps are large. After the loop, the WP is
renormalised:

```
ψ_wp ← ψ_wp / ‖ψ_wp‖.
```

**Implementation.** `inq-stack/include/inqkit/wavepacket/wavepacket.hpp`,
the `if (do_ortho_) { ... }` block inside `inject_into_last_extra_state`:

```cpp
for (int i = 0; i < ist_wp; ++i) {
    // ⟨ψ_i | ψ_wp⟩ — real and imag parts via two GPU reductions
    auto res_re = gpu::run(1, gpu::reduce(n_pts), 0.0, ...);
    auto res_im = gpu::run(1, gpu::reduce(n_pts), 0.0, ...);
    // ψ_wp ← ψ_wp − (re + i·im) · ψ_i
    gpu::run(..., [=] GPU_LAMBDA(...) {
        phicub_[ix][iy][iz][ist_w] -= ...;
    });
}
// Renormalise
auto norm_after = sqrt(reduce(|ψ_wp|²·dV));
phicub_ ← phicub_ / norm_after;
```

The `InjectionReport` records `norm_before`, `max_overlap`,
`norm_after`, and a `passed_tolerance` boolean (max overlap below
the user-supplied tolerance × 10).

---

## 4. Wave-packet construction

**Definition.** The injected wave packet is a 3D Gaussian with momentum
phase factor:

```
ψ_wp(r) = (π σ²)^{-3/4} · exp(-|r - b|² / (2 σ²)) · exp(i k₀ · r)
```

where `b = (b_x, b_y, b_z)` is the WP centroid and `k₀ = (k_x, k_y, k_z)`
is the mean wave vector. The amplitude `(π σ²)^{-3/4}` makes the
real-space Gaussian normalised on its own.

**Coordinate sampling.** INQ stores the real-space grid in FFT-natural
order; coordinates are obtained from
`basis.point_op().rvector_cartesian(ix, iy, iz)`, which internally calls
`to_symmetric_range` so that grid index `ix` maps to physical
`(ix < N/2) ? ix*dr : (ix - N) * dr` — i.e. the centred frame
`r ∈ [-L/2, +L/2]`. Using `ix * dr` directly would silently corrupt
the phase across half the cell.

**Implementation.** `inq-stack/include/inqkit/wavepacket/wavepacket.hpp`,
the GPU kernel inside `inject_into_last_extra_state`:

```cpp
gpu::run(basis.local_sizes()[2], basis.local_sizes()[1],
         basis.local_sizes()[0],
         [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
    auto rvec = point_op_.rvector_cartesian(ix, iy, iz);  // centred frame
    double r2 = (rx-bx)² + (ry-by)² + (rz-bz)²;
    double amp = (π·σ²)^{-3/4} · exp(-r²/(2σ²));
    double ph  = k_x·rx + k_y·ry + k_z·rz;
    phicub_[ix][iy][iz][ist_w] = complex(amp·cos(ph), amp·sin(ph));
});
```

---

## 5. End-of-box arrival time and per-screen physics windows

**End-of-box time.** The propagation is run for the time it takes the
trailing edge of the WP — the `centroid + 1σ` point — to reach the far
end of the box at the WP's mean speed `|k₀|`:

```
t_end = (b + σ + L_z / 2) / |k₀|
N_steps = round(t_end / dt)
```

For the Tsubonoya base config (`b = 12 Bohr, σ = 1.0 Bohr,
L_z = 60 Bohr, |k₀| = 3.834 Bohr⁻¹, dt = 0.020 a.u.`), this gives
`t_end ≈ 11.20 a.u.` ⇒ `N_steps = 560`.

**Per-screen windows.** A screen at z-position `z_s` integrates only
during the interval the WP is physically present at it.

* **Forward screens** (`z_s < b`):
  ```
  t_start = (b − z_s) / |k₀|       # centroid arrives at the screen
  t_end   = (b + σ − z_s) / |k₀|   # trailing edge clears
  ```
* **Backscattering screens** (`z_s ≥ b`): a model that assumes a
  rebound off the molecule at z = 0:
  ```
  t_start = (b + z_s) / |k₀|       # rebound centroid arrives
  t_end   = (b + σ + z_s) / |k₀|   # rebound trailing edge clears
  ```

Each screen has its own accumulator gated on `[t_start, t_end]`. A
single global "paper window" `[T1_AU, T2_AU] = [3.18, 10.34] a.u.`
also runs in parallel on every screen for direct paper-figure
comparison.

**Implementation.** `ResearchProject/systems/coronene/shared/configs/tsubonoya_2014_base.hpp:44-58`
(`compute_n_steps`) and
`ResearchProject/systems/coronene/shared/cpp/leed_screen_layout.hpp`
(`compute_screen_window`).

---

## 6. Three density categories at write time

At each saved RT step, three real-valued density fields are emitted:

```
ρ_system(r, t) = density::total(electrons)               # occupied orbitals only
ρ_wp(r, t)     = |ψ_wp(r, t)|²                            # WP orbital alone
ρ_total(r, t)  = ρ_system(r, t) + ρ_wp(r, t)              # full visible density
```

INQ's `density::total(electrons)` does **not** include the WP extra
state (per `docs/observables_reference.md:27`); the C++ run template
adds the two fields pointwise to obtain `ρ_total`.

**Implementation.** `ResearchProject/systems/coronene/shared/cpp/run_template.hpp`,
the `add_real_fields(a, b)` helper:

```cpp
inqkit::fields::RealField3D out = a;
for (std::size_t i = 0; i < a.values.size(); ++i)
    out.values[i] = a.values[i] + b.values[i];
return out;
```

The two source fields share the same grid layout because they are
produced from the same `electrons` object on the same iteration.

---

## 7. LEED screen-data → image-extent fftshift convention

LEED `.dat` files are written by `LeedPatternAccumulator::save()` in
**FFT-natural** order: array index `(0, 0)` corresponds to physical
origin `(x = 0, y = 0)`, with positive coordinates first then wrapped
negative. To plot in centred coordinates, the loader applies a 2D
`np.fft.fftshift` and overrides the origin so the
`(extent_x_min, extent_x_max, extent_y_min, extent_y_max)` spans
`[-L_x/2, +L_x/2, -L_y/2, +L_y/2]`.

**Implementation.** `inq-stack/python/inqview/screens.py`,
`load_leed_pattern`:

```python
data = np.fft.fftshift(data)
origin_x_bohr = -0.5 * nx * dx_bohr
origin_y_bohr = -0.5 * ny * dy_bohr
```

The companion `coordinate_checks/` raw-index plot uses
`np.fft.ifftshift(pat.data)` so the un-shifted FFT-natural layout is
also visible side-by-side — this is the diagnostic `run_06` originally
identified.

---

## 8. ParaView 3D volume rendering with log-scale + density-tied opacity

For each density series, ParaView's volume mapper builds:

* a **colour transfer function** `LUT(ρ)` with `UseLogScale = 1`,
  rescaled to `[ρ_min, ρ_p99]` so the dim WP tail spans most of the
  colour range while the bright peak still saturates, and
* an **opacity transfer function** `pwf(ρ)` that is piecewise linear
  with three control points (no log scaling on opacity itself):

  ```
  pwf(0)    = 0.00
  pwf(p50)  = 0.05
  pwf(p99)  = 0.60
  ```

  i.e. the bulk volume (below the median) is nearly transparent, the
  median region sees just enough opacity to define a soft envelope,
  and the peak region (above the 99th percentile) is mostly opaque.

The system density uses a "blue" colour preset; the WP density uses
"orange". The two volumes share the camera and the viewport so the
overlay reads at a glance.

**Implementation.** `inq-stack/python/inqview/postprocess/paraview_3d.py`,
`make_volume()` inside `_pvbatch_script()`:

```python
lut = GetColorTransferFunction(array_name)
lut.RescaleTransferFunction(*scalar_range)
lut.UseLogScale = 1
lut.ApplyPreset(color_preset, True)
pwf = GetOpacityTransferFunction(array_name)
pwf.Points = [0,0,0.5,0,  p50,0.05,0.5,0,  p99,0.6,0.5,0]
```

Two cameras are configured by explicit `(position, focal_point,
view_up)` rather than azimuth/elevation chains, so the views are
reproducible across pvbatch invocations.

---

## 9. Multi-line plot-title convention for animations

Per the visualisation rules (`docs/visualisation-instructions-v1.md`),
every animation frame's title is two lines:

```
<run_name>: <data type / scale tag>
step k/N, t = X.XX fs    ← rounded to 3 sig figs
```

**Implementation.** `inq-stack/python/inqview/postprocess/_common.py`,
`title(run_name, what, step=k, total_steps=N, time_au=t_au,
multiline=True)` returns a `\n`-joined two-line string. All animated
phases (density, screens, overlap) call this helper.
