# New-features backlog — inqkit + inqview

The feature **ideas** surfaced in TODO comments across both libraries, promoted
to concrete next-phase work items. Bugs and the restructure live elsewhere
(`inqview-findings.md`, `inqkit-errors.md`, `inqview-restructure-and-tests.md`);
this file is *new capability* only.

**Rule (user-mandated):** every new feature ships with a test that proves it
works (known/analytic input → expected output) before it is used downstream.

Status key: `decided` (locked in findings, ready to build) · `new` (idea
captured this round, needs a quick design pass) · `future` (parked).

## Build progress (2026-06-10) — user selected 10 features
**BUILT + tested ✓ (8):**
- canonical theme → `inqview/visualisation/style.py` (test_theme, 6).
- fourier `subtract=` + coherent-gain fix → `inqview/fourier.py` (IV-M12 + IV-E03;
  the IV-E03 xfail flipped to a real PASS; test_fourier 17).
- `energy_components` → `inqview/analysis/energy_components.py` (test 6; renderer
  bars/lines/GIF still to add in visualisation).
- `wp_integrity` (ipr/momentum_kl/variance/kl_series + dataclass) →
  `inqview/analysis/wp_integrity.py` (test 11; from-run assembly follows).
- `plasmon_spectrum` (peak-locator, complex-FFT, |n_q|²/q², axial extraction) →
  `inqview/analysis/plasmon_spectrum.py` (test 4; fixes IV-E01/E02).
- KL drift-rate (`kl_series` initial/previous) + `(k,t)` carpet renderer →
  `inqview/visualisation/carpets.py` (carpet untested per IV-M10).
- Python `center_of_density` (node convention + WP/total/bath compare, E04
  dx/2 cross-check) → `inqview/analysis/center_of_density.py` (test 5).
- **norm-per-state** (inqkit) → `observables/state_norm_writer.hpp` +
  `test_state_norm_engine.cpp` (every orbital ∫|ψ|²≈1). **BUILT + PASSED**
  (ctest 2.67 s). inqkit engine tier now 19 files / 21 cases.

inqview suite after this batch: **68 passed, 1 xfailed** (only deps-clean).

**REMAINING (2 inqkit) — specs below:**

### current+dipole as Vec3 — spec
- **Where:** `io/observables_writer.hpp` (+ `detail/vec3.hpp`).
- **What:** store current `(Jx,Jy,Jz)` and dipole `(dx,dy,dz)` internally as
  `inqkit::detail::Vec3` units (consistency with `center_of_density`); CSV column
  names UNCHANGED (current_x/y/z, dipole_x/y/z) for backward compat.
- **Test (engine):** run a short propagation with current+dipole selected; read
  back the CSV and assert the row equals the Vec3 components (round-trip), and a
  known-vector unit test of the Vec3 → CSV formatting. Pure Vec3 parts can also
  go in the existing `test_vec3.cpp` (pure tier).

### N-dim plane screen — spec
- **Where:** `screens/plane_screen.hpp` (retires the z-only limitation at line 4).
- **What:** add an `axis` (0=x,1=y,2=z) / arbitrary-normal parameter to
  `PlaneScreen`; `iz_nearest`-style index logic generalised per axis with the
  same FFT-natural wrap.
- **Test (engine):** build a toy field with a known plane; extract the slice
  along x, y, and z; assert each equals the analytic plane values. Add a
  time-averaged-screen variant (`Σ_t ρ·dt / T` of constant frames == the frame).

---

## inqview — new capabilities

| Feature | Source | What to build | Proving test | Status |
|---|---|---|---|---|
| `energy_components` | energy_balance.py:31; todo_later §3; todo.txt | `analysis` kernel: E_kin/H/xc/E_ext(residual) flow from observables.csv → `EnergyComponents` dataclass; viz: initial-vs-final bars + ΔE(t) lines + GIF | `Σ components == E_total`; `Σ ΔE == ΔE_total`; free-WP energy conserved | decided (IV-M07) |
| `wp_integrity` | kl_divergence.py:49 | `analysis`: `WPIntegrity(kl_mom, σ_r, ipr)` | free-WP → kl_mom≈0, σ_r(t) analytic, IPR decays | decided (IV-M05) |
| `plasmon_spectrum` | density_fourier.py:28,33 | rename loss-fn → peak-locator; `axial` + `3d_binned` modes; `PlasmonSpectrum` dataclass | undamped-plasmon phasor → δ-peak at ω_p + 1/q² scaling | decided (IV-M01/M04) |
| `center_of_density` (Py) | wake.py:81 | `analysis`: COD from VTI (node convention) + `CODComparison(wp, total, bath)` | python-COD vs inqkit-CSV differ by exactly dx/2 (documents E04) | decided (IV-M02) |
| KL drift-rate + `(k,t)` carpet | kl_divergence.py:46,53 (Runfeng) | frame-to-frame `KL(P_t‖P_{t−Δ})`; contour of the existing `(time × k-bin)` grid | KL≥0; carpet shape matches grid | decided (IV-M05) |
| canonical theme | _common.py:36; plots.py:40; overlap.py:74; observables.py:111 | `visualisation.style`: semantic cmap roles + fixed-dim figure factory | `figure_one_col()==(3.5,3.0)`; `cmap_for('diverging')=='RdBu_r'` | decided (ADR-0004) |
| fourier `subtract=` | todo.txt #1,#2 | `subtract={'initial','mean','detrend','none'}`; canonical per column | tone+offset+drift → peak unshifted, DC suppressed | decided (IV-M12) |
| band structure (ε vs k-path) | orbitals_per_kpoint.py:30 | multi-k dispersion plot | parabola/ε(k) recovered for a known multi-k case | future (IV-M09) |
| orbital×component energy | energy_balance.py; IV-M07 | per-orbital kinetic split (Hartree/xc non-attributable) | needs saved ψ; kinetic split sums correctly | future |

## inqkit — new capabilities

| Feature | Source | What to build | Proving test | Status |
|---|---|---|---|---|
| N-dim plane screen | screens/plane_screen.hpp:4 | generalise `PlaneScreen` beyond the z axis (x/y/arbitrary normal) | extract a known slice along each axis of a toy field | new |
| time-averaged screen | screens/plane_screen.hpp:5 | accumulate a time-averaged total-density screen | Σ_t ρ·dt / T equals the mean of constant frames | new |
| any-k-point support | wavepacket.hpp:180; observables (gamma-only guards) | lift the gamma-only restriction in WP inject + stats | multi-k toy: per-k moments correct (enables QKE band structure) | new |
| momentum-space Gram-Schmidt | wavepacket.hpp:67 | reciprocal-space MGS variant; compare to real-space | both routes give same orthonormal set within tol | new (experiment) |
| projected-occupation snapshot | todo.txt #7 | re-enable `OrbitalOverlapMatrix::snapshot()` (full, not WP-only) at coarse cadence | t=0 identity `n_i^GS(0)=f_i(0)` (pairs with inqview gs_projected) | new |
| current+dipole as Vec3 | observables_writer.hpp:60 | track current/dipole as a `Vec3` unit (consistency with COD) | round-trip a known vector through the writer | new |
| Gaussian projectile potential | todo.txt #9 | optional Gaussian (vs radial) external potential for scattering | potential matches analytic Gaussian at sample points | new (sim) |
| norm-per-state diagnostic | todo_later §"Other" 2 | `norm_per_state.csv` (Σ|ψ_i|²dV per evolved orbital) | injected WP norm ≈ 1; catches >1e-3 leakage | new |

## Notes
- The inqview `decided` rows are already designed (findings IV-M*) — they get
  built during/after the package restructure, each behind its proving test.
- The inqkit `new` rows need a quick design pass (and, where formula-bearing,
  the formula-validation agent) before implementation — same discipline as the
  inqkit characterization round.
- `screens/plane_screen.hpp` carries TWO ideas (N-dim + time-averaged); the
  N-dim generalisation also retires the z-only limitation noted in
  `todo_later.md` "MPI-aware slice extraction" lineage.
- Stub files (`text_io`, `manifest_writer`, `filesystem`, `text_summary_writer`,
  `validation`) are "write-if-necessary" placeholders — evaluate need before
  writing; not features per se.
