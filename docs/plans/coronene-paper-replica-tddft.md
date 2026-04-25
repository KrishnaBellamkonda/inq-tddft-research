# Plan: TDDFT replica of Tsubonoya 2014 with full observable suite

## Context

Goal: re-validate the `fixes/coronene-gs` branch with a full TDDFT
simulation that mirrors the parameters of the Tsubonoya, Hu, Watanabe paper
(*Time-dependent density-functional theory simulation of electron
wave-packet scattering with nanoflakes*, Phys. Rev. B 90, 035416 (2014)),
and record every observable in `docs/observables_reference.md`. Once the
ground state and the observables come out cleanly, the writer-fix branch
is ready to merge.

The diagnostic run `run_06_centred_writer_check` proved that:

- INQ's SCF is healthy when atoms are truly centred at the origin
  (E_total = -151.42 Ha for coronene PBE, integrated charge = exact 108
  electrons, HOMO/HOMO-1 doubly degenerate).
- The cross artefact in earlier renders was a writer-side index mismatch,
  fixed in `inqkit::fields::density::total/orbital` by
  `fft_shift_index(idx, size) = (idx + (size+1)/2) % size`.

This plan covers the next, decisive validation: reproduce the paper's
coronene + 200 eV WP run end-to-end and confirm:

1. The ground-state density renders correctly (no cross).
2. The TDDFT propagation conserves N_electrons, energy is reasonable.
3. Observables (overlap matrix, dipole, current, screens) match the
   physical expectations from the paper (LEED hexagonal pattern,
   π-plasmon dipole peak around 5 eV, π+σ peak around 17 eV).

## Paper parameters (translated to atomic units)

Reference: Tsubonoya 2014, Eqs. (1)–(9), §III.

| Quantity | Paper value | Atomic units |
|---|---|---|
| Cell | 18.4 × 18.4 × 31.7 Å³ | 34.7710 × 34.7710 × 59.9043 Bohr |
| Boundary | finite (open) | `cell.finite()` |
| XC | ALDA | `options::theory{}.lda()` |
| Pseudopotentials | Troullier–Martins (Ref [24]) + Kobayashi (Ref [25]) | INQ default NC PSPs (or pass Qball ONCV `.xml` to override — TBD) |
| C–C bond | 1.42 Å (fixed) | from the centred xyz; ions clamped |
| Spin | unpolarised | INQ default |
| Grid spacing | 0.16 Å | 0.302 Bohr (≈ matched at INQ cutoff 54 Ha → Δr ≈ 0.291 Bohr) |
| Δt (LEED window) | 4.84 × 10⁻⁴ fs | 0.020 a.u. |
| WP standard deviation `d` | 0.53 Å | 1.0015 Bohr |
| WP centre `b` (offset from flake) | 6.35 Å, perpendicular to flake plane | 12.0 Bohr along z |
| WP wave vector `k` | E_kin = 200 eV | \|k\| = √(2 × 200/27.2114) = 3.834 Bohr⁻¹ along −z (toward flake) |
| LEED window | t₁ = 0.077 fs to t₂ = 0.25 fs | t₁ = 3.18 a.u. → t₂ = 10.34 a.u. |
| Total LEED propagation | 10000 steps | 200 a.u. = 4.84 fs |
| Observation plane | 8 Å from flake | 15.12 Bohr |

Notes
- The paper uses **ALDA**; our existing `coronene-wp-rt` runs used LDA
  (matches). We propose to keep ALDA for fidelity to the paper, *not* PBE
  (which the diagnostic runs used). If you'd prefer PBE for parity with
  `run_06`, that's a one-line change.
- The paper's pseudopotentials are not directly available; INQ's default
  NC set is the natural substitute. Passing the Qball ONCV `.xml` files
  is a one-line override if you want exact parity with the qball
  reference computation we already have.

## Coordinate convention (post-fix branch)

INQ uses `[-L/2, +L/2]` for an orthorhombic cell. The flake therefore sits
at `z = 0` with atoms in `coronene_centred.xyz` (taken verbatim from
`coronene-qball/coronene.sys` lines 14–50, the same file used in
`run_06_centred_writer_check`). The WP launch position becomes:

```
WP centre b = (0, 0, +D)  with  D = 6.35 Å = 12.0 Bohr
WP wave vector k = (0, 0, -|k|)  with  |k| = 3.834 Bohr⁻¹
```

so the WP is placed at `z = +12 Bohr` and travels toward `z = 0` (the
flake plane) and onward to `z < 0`.

The pre-SCF defensive cell-bounds check from `run_06` will be reused; it
aborts if any atom falls outside `[-L/2, +L/2]`.

## Cell sizing

The paper uses Lz = 31.7 Å = 59.90 Bohr. The WP at +12 Bohr leaves room
for `Lz/2 - 12 = 17.95 Bohr ≈ 9.5 Å` of vacuum behind the WP, and similar
in front of the flake. Adequate for the LEED window. We can keep this
exactly. *(Existing `coronene-wp-rt/run_01_d635_base` widened Lz to
89.856 Bohr and offset the flake; we don't need that — at the post-fix
convention the flake at z=0 with paper Lz suffices.)*

## Observable suite (mirror of `docs/observables_reference.md`)

| # | Observable | C++ source | Output | Frequency |
|---|---|---|---|---|
| 1 | GS total density (no WP) | `density::total(electrons)` before injection | `results/density_gs/` | once |
| 2 | GS orbital densities | `density::orbital(electrons, i)` for `i = 0 … wp_idx-1` | `results/density_gs_orbitals/orbital_XXXX/` | once each |
| 3 | RT total density (target + WP) | `density::total + density::orbital(wp_idx)` | `results/density_rt_total/` | every WRITE_EVERY |
| 4 | RT target-only density | `density::total(electrons)` | `results/density_rt_jellium/` *(name kept for analysis-script reuse; coronene plays the "target" role)* | every WRITE_EVERY |
| 5 | RT WP orbital density | `density::orbital(electrons, wp_idx)` | `results/density_rt_wp/` | every WRITE_EVERY |
| 6 | KS overlap matrix `O_ij(t)` | `OrbitalOverlapMatrix::snapshot()` | `results/overlap/overlap_XXXXXX.csv` + `index.csv` | every step |
| 7 | E_total, E_kin, E_Hartree, E_xc(t) | `ObservablesWriter` | `results/observables.csv` | every step |
| 8 | Current J_x, J_y, J_z(t) | `ObservablesWriter` | `results/observables.csv` | every step |
| 9 | Dipole μ_x, μ_y, μ_z(t) | `ObservablesWriter` | `results/observables.csv` | every step |
| 10 | Time-averaged LEED screens (20 z-positions) | `LeedPatternAccumulator` | `results/screens/screen_NN.dat` | end of run |
| 11 | Instantaneous LEED snapshots | `PlaneScreen::extract` | `results/screens_snapshots/step_XXXXXX/screen_NN.dat` | every 3 steps |

## Header / config split

Per your request, paper parameters live in a header so any future run.cpp
can `#include` the same configuration:

- **NEW**: `inq-stack/include/inqkit/config/tsubonoya_2014_coronene.hpp`
  - `static constexpr double LX_BOHR = 34.7710;`
  - `LY_BOHR`, `LZ_BOHR`
  - `WP_SIGMA_BOHR = 0.53 Å × ANG_TO_BOHR`
  - `WP_EKIN_HA = 200 eV / HA_TO_EV`
  - `WP_K0 = sqrt(2 * WP_EKIN_HA)`
  - `WP_OFFSET_BOHR = 6.35 Å × ANG_TO_BOHR` (= +D along z, flake at z=0)
  - `OBSERVATION_PLANE_BOHR = 8 Å × ANG_TO_BOHR` (= +15.12 Bohr; not all
    physically meaningful for our 20-screen layout — kept for reference)
  - `DT_AU = 0.02`
  - `N_STEPS_LEED = 10000`
  - `WRITE_EVERY = 100` (frequent enough for a 200 a.u. window → 100
    density frames; matches `coronene-wp-rt/run_01_d635_base`)
  - `SCREEN_SNAP_EVERY = 30` (jellium uses 3, but at 10000 steps that is
    >3000 snapshot directories per screen; we coarse-grain to keep disk
    use under control — confirm this is acceptable in approval phase)

- **NEW**: `Tutorial/coronene-leed/run_diagnoses/run_07_paper_replica/`
  - `coronene_centred.xyz` (copy of `run_06_centred_writer_check`'s xyz)
  - `run.cpp` includes `inqkit/config/tsubonoya_2014_coronene.hpp` and
    runs the SCF → WP injection → propagation → write all 11 observables.
    Structurally a near-clone of `jellium-wp-rt/run_01_base/run.cpp` with
    coronene system + paper parameters and the post-fix density extractor.
  - `analysis.py` adapted from `jellium-wp-rt/run_01_base/analysis.py`
    (reuses the full diagnostic ladder: N_e conservation, density
    consistency, energy spectrum, current spectrum, density-slice GIFs,
    LEED panel, overlap GIFs).

## Overlap matrix GIF code: investigation step

You mentioned an error in the overlap-matrix GIF production code. I read
the existing implementation in
`ResearchProject/jellium/jellium-wp-rt/run_01_base/analysis.py:255-302`
and the C++ writer at `inq-stack/include/inqkit/observables/orbital_overlap.hpp`.
The C++ side is consistent with the documented schema (rows = GS index `i`,
columns = evolved index `j`, columns 0..n_ref-1 the occupied block,
column n_ref the WP). The Python loader skips the leading `# step= …`
header line via the `try/except ValueError` in `_load_overlap_matrix`,
which on the existing jellium data produces the identity matrix at t=0
and decay thereafter — i.e. *looks* correct on the surface.

What I need from you to land a confident fix:

- Which behaviour did you observe? Candidates I can imagine:
  - **(a)** Off-diagonal GIFs look flat / invisible because `ymax` is
    pinned by the diagonal-1 spike when iterating across all frames →
    all subsequent off-diagonal columns get rescaled by the *same*
    diagonal-frame ymax that was computed for that particular j.
    *(Fix: per-j ymax already, so this isn't it. Re-checking…)*
  - **(b)** Mislabelled axes (the header says `|O_{ij}|^2`, but the
    written value is already a squared magnitude — so for an unmoved
    occupied orbital we expect 1, not 1²; that is fine).
  - **(c)** GIFs out of memory or slow due to per-frame matplotlib state.
  - **(d)** Wrong ordering of rows vs columns (transpose).
  - **(e)** A frame skipped because of the `if m.shape == (n_ref, n_evolved):`
    filter when the matrix file is partially written.
  - **Other** — please describe.

The plan is to:

1. Reproduce the symptom you saw using the existing jellium overlap data.
2. Add a regression test under `inq-stack/tests/python/`
   (`test_overlap_loader.py`) that builds a synthetic identity-at-t0 +
   decaying-diagonal series, runs the loader, and asserts the loaded
   matrices have the expected shape, ordering and values.
3. Fix the code, re-run on jellium data, verify the GIFs render
   correctly, then port the fix into the paper-replica run's analysis.

## Critical files / new files

### Created
- `Tutorial/coronene-leed/run_diagnoses/run_07_paper_replica/run.cpp`
- `Tutorial/coronene-leed/run_diagnoses/run_07_paper_replica/coronene_centred.xyz` (verbatim copy from run_06)
- `Tutorial/coronene-leed/run_diagnoses/run_07_paper_replica/analysis.py`
- `inq-stack/include/inqkit/config/tsubonoya_2014_coronene.hpp`
- `inq-stack/tests/python/test_overlap_loader.py` (overlap loader regression test)

### Modified
- Whatever specific function in
  `ResearchProject/jellium/jellium-wp-rt/run_01_base/analysis.py` proves
  to be at fault, factored into a shared helper under
  `inq-stack/python/inqview/overlap.py` (then imported back into the
  per-run analysis scripts).

### Reused (no edits)
- `inqkit::WavePacket` (`inq-stack/include/inqkit/wavepacket/wavepacket.hpp`)
- `inqkit::observables::OrbitalOverlapMatrix`
- `inqkit::screens::PlaneScreen`, `LeedPatternAccumulator`
- `inqkit::io::ObservablesWriter`, `RealField3DWriter`
- `inqkit::RealTimeSession`, `StepContext`
- `inqkit::fields::density::total/orbital` (post-fix; correct rendering)
- `inqkit::fields::orbital::wavefunction` for the overlap snapshot

## Validation menu (per `.claude/rules/testing.md`)

### Tier A (always run, fast)
- SCF converged at 1e-4 Ha (paper-style tolerance; 1e-6 if cheap)
- E_total negative and within 1 % of `run_06`'s -151.42 Ha
- ∫ ρ_gs dV = 108 (electrons)
- WP injection report: `norm_after ≈ 1.0`, `max_overlap < 1e-6`
- Total density GIF visually shows the molecule at the metadata centre
  (smoke test of the writer fix end-to-end)

### Tier B (run, takes hours)
- N-electron conservation across the full TDDFT run < 0.1 % drift
- `total = jellium + wp` density consistency to < 1e-6
- Energy total drifts < 0.1 % over the LEED window
- Overlap matrix at t=0 is approximately identity on the n_ref block,
  and `O[i, wp_idx] < 1e-3`
- LEED screen `screen_10` (z = 21.07 Bohr midpoint) shows hexagonal
  symmetry (paper Fig. 2(a) reference)

### Tier C (deferred — requires explicit approval)
- Extend propagation by 25000 steps with WP removed, Δt = 0.04 a.u., for
  the dipole excitation σ_x^d(ω). This is the second half of the paper
  protocol. Total compute roughly 3× the LEED window. Skip unless the
  Tier-A/B results are clean.

## Resource estimate

`run_06` SCF on coronene at 54 Ha cutoff took ≈ 18 min wall on one A30
plus ≈ 7 min for the orbital writes. The TDDFT propagation cost is
roughly `N_steps × per-step-cost`. At 10000 steps × 0.5 s/step on a
single A30 (rough estimate from the existing coronene-wp-rt logs), the
LEED window is ≈ 90 min. Plus orbital snapshots, screen accumulation, and
overlap CSVs add ~20–30 % overhead. Expected wall-time: **2.5–3 h on one
GPU**. We have two A30s; we can pin this to GPU 0.

## Decision points awaiting your input

1. **XC**: ALDA (paper) or PBE (run_06 baseline)?
2. **Pseudopotentials**: INQ default or Qball ONCV `.xml`? Preference is
   to keep INQ default for portability and only switch if Tier-A energies
   are noticeably off. Confirm.
3. **Tolerance**: `energy_tolerance(1e-4_Ha)` (paper-style) or `1e-6_Ha`
   (`run_06` style, costs ≈ 3× more iterations)?
4. **`SCREEN_SNAP_EVERY`**: 3 (jellium parity, ~3300 dirs / screen,
   ~66000 dirs total) or 30 (10× coarser, ~6600 dirs total)? Confirm.
5. **`WRITE_EVERY` for densities**: 100 (matches `coronene-wp-rt/run_01_d635_base`,
   100 frames over 10000 steps) — confirm.
6. **Overlap GIF symptom**: which of (a)–(e) did you observe? Or other?

## What this plan does NOT include

- The 25 000-step dipole-excitation phase (Tier C, deferred).
- Re-running `run_01..run_05` raw outputs through the corrected
  pipeline (separate small task, no SCF needed).
- Merging `fixes/coronene-gs` into main — happens after this run is
  validated.

## Approval gate

I will not run the simulation until you reply with:
- selections for items 1–6 above (or "defaults" → I'll proceed with
  ALDA / INQ default PSPs / 1e-4 Ha / SCREEN_SNAP_EVERY=30 /
  WRITE_EVERY=100 / overlap symptom = unspecified, I'll diagnose).
- explicit "go" to launch on GPU 0.
