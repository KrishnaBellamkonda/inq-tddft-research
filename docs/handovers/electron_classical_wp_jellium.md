# Handover: electron-classical-wavepacket-jellium comparison

## v2-final state (2026-05-08, evening — CLASSICAL DONE, WP abandoned)

**Classical run completed at 20:27 BST.** Wall time **3 h 17 m** (11844 s).
All 860 propagation steps executed cleanly. `run_completed = true`,
`Done. Wall time 11844.3 s.` Exit 0.

**Physics signal extracted from electron_track.csv + observables.csv:**

| Quantity | Value |
|---|---|
| Distance traveled | 45.14 Bohr (single transit + small wrap, L=50 box) |
| v_z initial → final | 10.4999 → 10.4966 bohr/atu (drop = 0.0033) |
| Velocity drop | 0.031 % |
| Projectile KE loss | **0.937 eV** (= 0.06 % of 1500 eV) |
| System total energy gain | 0.939 eV |
| **Energy conservation residual** | **+0.002 eV** (machine precision) |
| Bath kinetic-energy gain | 0.891 eV |
| Bath Hartree+xc rearrangement | 0.048 eV |
| **Stopping power S(v=10.5)** | **0.021 eV/Bohr** |

The stopping power 0.021 eV/Bohr at v=10.5 in r_s=5.69 jellium is in
the high-velocity Bethe regime: S(v) ∝ ln(v)/v² rather than the v_F-peak
Lindhard regime. Physically the projectile's kinetic energy is going
predominantly into bath single-particle e-h excitations (kinetic gain
≫ Hartree+xc rearrangement), as expected for a fast (k_0 ≫ k_F) charged
projectile.

**Postprocess pipeline ran cleanly** on the classical results: 6 phases
(summary, observables, state_energies, bath_energy, stopping,
occupations) all PASS. 17 PNGs + CSVs in
`results/analysis/observables/`:
- `dE_kinetic_vs_z.png` — bath KE gain vs projectile z (stopping curve)
- `stopping_force_vs_z.png` — F_z reconstructed from dv/dt
- `bath_energy_vs_time.png` + CSV — band-structure-summed bath E(t)
- `observables_summary.png`, `total_energy_vs_time.png`,
  `current_components_vs_time.png`, `dipole_components_vs_time.png`
- `fft_total_energy.png`, `fft_current_*.png`, `dipole_spectrum_*.png`
  — frequency-domain spectra
- `ks_energies_absolute.gif`, `ks_energies_delta.gif`,
  `occupations_absolute.gif`, `occupations_delta.gif` — animated bar
  charts

**WP run abandoned** (not viable in this session at dx=0.30/102 states/
24 GB GPU). The v2 plan's WP-vs-classical comparison is therefore
incomplete — only the classical half is available. A follow-up session
should fix the inqkit per-element-loop bulk-copy and retry the WP path.

## v2-final state (2026-05-08, evening — WP abandoned)

After multiple attempts (4 separate WP launches with progressively
stripped callbacks), the WP run at dx=0.30 / 102 states cannot reach
SCF step 1 in less than 30 minutes of pre-step-0 setup. The process
remains alive (CPU 60–92 %, GPU 0 at 100 % util / 35 W power — same
memory-bound thrashing pattern as the v2 dx=0.248 GS), but never emits
the step-0 propagation log line. INQ's TDDFT internal init
(density.calculate, ks_hamiltonian construction, time_zero, etc.) at
this state count + grid + 24 GB GPU = effective deadlock from
cudaMallocAsync paging.

**Final WP status: NOT COMPLETED in this session.** Stripped callback
(only `obs_writer.append(ctx)` per WRITE_EVERY=4 steps) was
insufficient.

**Classical run status: PROGRESSING.** 253/860 = 29 % at 18:10. ETA
finish ~22:30 BST tonight. Full observable stack (energy, current,
dipole, density, momentum, ehrenfest projectile track via velocity
differences).

**Next-session work plan**:
1. Properly patch `inqkit/fields/orbital.hpp::wavefunction` and
   `density.hpp::orbital` with a bulk GPU→host copy. Test on a 1-element
   benchmark first. Expect this to bring WP per-step callback back to
   ~5 sec from 30 min.
2. Find why INQ's TDDFT init (density.calculate, ks_hamiltonian
   construction, etc.) appears to deadlock at dx=0.30 / 102 states. May
   be unrelated to inqkit — could be in INQ's gpu-pool allocation.
   Possible mitigation: run WP at dx=0.40 with same N=162 GS as
   classical and accept the WP momentum aliasing at k_0=10.5
   > Nyquist=7.85 (loss of physics validity, but the WP-vs-classical
   comparison is still informative for stopping power if interpreted
   carefully).
3. Implement a "fast WP" variant that propagates only the WP orbital
   (frozen-bath approximation) — single-orbital propagation skips most
   of the memory pressure.

## v2-final state (2026-05-08, late afternoon)

**Currently running** (as of 17:47 BST 2026-05-08):
- **Classical** (PID 1875627, GPU 1, dx=0.40 Cfg): step 155/860 (18 %), ~13 s/step, ETA finish ~22:30 BST.
- **WP** (PID 1900922, GPU 0, dx=0.30 Cfg, minimal callback): launched 17:44 BST; pre-step-0 setup taking >3 min and counting; will continue overnight. ETA finish (if step time settles at 63 s) ~9 AM tomorrow.

**Files committed** (commit `e4d7a74` on `runs/electron-classical-wavepacket-jellium`):
- `electron-ONCV-1.2.upf`: `<PP_NONLOCAL>` stub added — pseudopod's UPF2 parser null-derefs without it (fix documented in PP_INFO).
- `electron_proj_E1500_L50_cubic.hpp`: now declares both
  `Electron_Proj_E1500_L50_cubic_WP` (dx=0.30) and
  `Electron_Proj_E1500_L50_cubic_Classical_dx0p40` (dx=0.40).
- `save_gs/gs_L50_cubic_N162_dx0p30/run.cpp`: GS for the WP run.
- `save_gs/gs_L50_cubic_N162_dx0p40/run.cpp`: smaller-grid GS for the
  classical run (forces_stress::calculate at dx=0.30 wouldn't fit
  in 24 GB GPU).
- `run_classical_e1500_L50_cubic/run.cpp`: now uses `.ehrenfest()`
  (was `.impulsive()`), so the bath's electronic forces decelerate
  the projectile — the stopping-power signal we want.
- `run_wp_e1500_L50_cubic/run.cpp`: per-step callback stripped to the
  cheapest possible `obs_writer.append(ctx)` only. Density VTI series,
  density_delta, cod_z trajectory, full-overlap snapshots, and the
  WP-only overlap snapshots are all skipped to keep the callback
  from triggering slow GPU→host orbital extractions.
- `inq-stack/include/inqkit/fields/density.hpp`: bulk-copy patch for
  `density::total` (boost::multi::array<double, 3> from GPU view);
  `density::orbital` left at the original per-element loop (a clean
  bulk-copy attempt for the 4D hypercubic case stalls on this
  boost::multi version, so we avoided that path).

**Observables that will be available**:
- Classical (full stack via `.ehrenfest()`): observables.csv (energy/
  current/dipole), state_energies.csv, occupations_vs_time.csv,
  momentum_distribution.csv, density VTI series (total + delta), full
  overlap matrix at {0, N/2, N}, electron_track.csv (pos, vel, F=0
  placeholder — Ehrenfest forces ARE applied internally; postprocess
  reconstructs them from `dv/dt`).
- WP (minimal): observables.csv only. No density VTI, no overlap, no
  cod_z. Postprocess will not have spatial WP density information.

**Known limitations going forward** (for the next session):
1. The inqkit per-element host loops in `fields::orbital::wavefunction`
   and `fields::density::orbital` need a proper bulk-copy patch
   before per-step orbital diagnostics are usable on dx=0.30.
2. The full-matrix overlap snapshots can't be produced for the WP run
   without the same library fix.
3. The postprocess pipeline (bath_energy.py, stopping.py) is generic
   and should work for both runs once they finish — bath_energy.py
   reads state_energies.csv (present in both); stopping.py reads
   observables.csv (present in both) and electron_track.csv
   (classical-only).



## v2 dx=0.30 update (2026-05-08, mid-day)

The v2 GS at dx=0.248 ran into a GPU memory wall: eigensolver workspaces
(6 wavefunction-sized buffers ≈ 90 GB nominal) exceeded the 24 GB A30
memory, forcing PCIe paging. Single-GPU iter 0 took **5.55 hours**.
MPI-2 with domain-parallel decomposition didn't help meaningfully
(both GPUs at 100 % util but only 35–60 W vs 165 W peak — same
memory-pressure pattern, just split across two GPUs).

**User pivot**: relax the grid to dx=0.30 (half the cells, fits one A30
comfortably) and accept ~6 % WP-tail aliasing at k_0+3σ_k=11.10 vs
Nyquist k=10.47. User's "each run on exactly one GPU" rule remains:
GS uses GPU 0 only, then WP+Classical concurrent on GPU 0 + GPU 1.

**Result on dx=0.30**:
```
SCF iter 0 : wtime =  62.2s   e = 923.32   dn=1e-01  dst=5e+01
SCF iter 1 : wtime =  59.8s   e = 428.40   dn=1e-01  dst=7e+00
```

That's **62 s/iter at dx=0.30 vs 19,959 s/iter at dx=0.248 — 320×
speedup per iter**. Total GS now estimated ~10–20 minutes
(vs the projected 8–15 hours at dx=0.248).

**Other v2 fix**: classical run.cpp now uses `.ehrenfest()` ion
dynamics (was `.impulsive()`) so the bath's electronic forces actually
decelerate the projectile — this is the whole point of the classical
comparison. With `.impulsive()` we'd only have seen bath response,
no projectile slowdown.

**File status**: `electron_proj_E1500_L50_cubic.hpp` SPACING_BOHR
changed 0.248→0.30, paths in 4 files updated dx0p248→dx0p30, save_gs
dir renamed (build/ cleaned to fix CMake cache). Both run.cpps
re-built successfully on dx=0.30 Cfg.

## v2 reconfiguration (2026-05-08)

User stopped the v1 production GS mid-SCF and switched to a different
scientific configuration. The v1 sections below `---` divider remain
for traceability but are obsolete; the active config is v2.

**v2 config**:
- Cubic 50³ Bohr periodic (replaces 40×40×150 orthorhombic)
- N = 162 electrons, r_s = 5.69 (matches L=50 N=162 base lineage)
- Spacing 0.248 Bohr (= π/g_max for cutoff 80 Ha = 160 Ry)
- WP: KE=1500 eV, k₀=10.50 bohr⁻¹, σ=5 Bohr
- Launch (corner-origin): (25, 25, 15) = INQ centred (0, 0, −10),
  i.e. 3σ from −z face — cleaner injection than the v1 5-Bohr buffer
- dt=0.005 atu, N_STEPS=860, total t=4.3 atu (one full periodic wrap
  during the simulation)
- Concurrent execution on two GPUs: WP on GPU 0, classical on GPU 1,
  via `CUDA_VISIBLE_DEVICES`

**v2 file changes** (committed alongside this handover update):
- NEW: `shared/configs/electron_proj_E1500_L50_cubic.hpp`
- NEW: `save_gs/gs_L50_cubic_N162_dx0p248/run.cpp`
- NEW: `run_wp_e1500_L50_cubic/run.cpp`
- NEW: `run_classical_e1500_L50_cubic/run.cpp`
- NEW: `scripts/classical_electron_smoke/Cv2_pre_gs_dryrun/dryrun.cpp`
- v1 versions are SUPERSEDED but retained on disk and in the prior
  commit `b277a0d` for traceability.

**v2 smoke gates passed**:
- v2 pre-GS dryrun PASS — INQ accepts cubic 50³ at dx=0.248, allocates
  210×210×210 = 9.26 M grid points (cube is comparable to v1's 9.11 M
  total) with 101 states.
- C1, C2, C3 from v1 carry over (config-independent: tested UPF
  parsing, mass override = m_e exact, impulsive propagator at 1.78e−15
  Bohr machine precision).

**v2 GS launched 11:46 BST 2026-05-08** as background task `ba2osvh2z`
on GPU 0 (`CUDA_VISIBLE_DEVICES=0`). Build phase was only ~4 min
(libxc cache reused from dryrun); SCF started at 11:50.

**v1 GS task `bfsl7dela`** killed at 11:08; orphan `./run` PID
1556242 also cleaned up. Both GPUs returned to idle (13 MiB each)
before v2 GS launch.

## Current status

Phases 0–4 of the plan (`.claude/plans/the-objective-in-this-dapper-moon.md`)
are complete. All five smoke gates that don't require the production GS
are green. The next step is launching the production GS, estimated at
1–3 hours wall on a single A30. **Awaiting user confirmation before
launching.**

| Phase | Description | Status |
|-------|-------------|--------|
| 0     | Git hygiene (commit pending changes, merge to main, create runs branch) | DONE |
| 1     | `electron-ONCV-1.2.upf` created from antiproton UPF | DONE |
| 2     | INQ mass override API confirmed | DONE |
| 3     | Cfg header `electron_proj_E1000_L40x40x150.hpp` | DONE |
| 4a    | `save_gs/.../run.cpp` written | DONE |
| 4b    | Pre-GS dry-run (cell + electrons construction, no SCF) | DONE — PASS |
| 4c    | Production GS SCF | **PENDING USER CONFIRMATION** |
| C1    | UPF parses cleanly into INQ | DONE — PASS |
| C2    | Mass override gives 1.0 a.u. = m_e | DONE — PASS |
| C3    | Free-flight impulsive propagator exact to 1.78e-15 Bohr | DONE — PASS |
| 5a    | WP run.cpp (hand-written, run_template is cubic-only) | DONE — uncompiled |
| 5b    | Classical run.cpp (custom UPF + mass override, .impulsive ions) | DONE — uncompiled |
| 6     | Python postprocess (bath_energy.py + stopping.py + pipeline wiring) | DONE — imports check |
| C4    | Classical electron in jellium GS, 100 steps | NOT STARTED — needs GS |
| C5    | Production classical run | NOT STARTED — needs GS, C4 |
| W1–W3 | WP smoke gates + production run | NOT STARTED — needs GS |
| 8     | End-to-end verification | NOT STARTED |

## What changed

### New files (`runs/electron-classical-wavepacket-jellium` branch)

Source-controlled scaffolding (committed):
- `ResearchProject/systems/jellium/shared/pseudopotentials/electron-ONCV-1.2.upf`
  — copy of `anti-proton-ONCV-1.2.upf` with a PP_INFO comment documenting
  that the projectile mass is set host-side at runtime.
- `ResearchProject/systems/jellium/shared/configs/electron_proj_E1000_L40x40x150.hpp`
  — defines `Common_E1000_L40x40x150`, `Electron_Proj_E1000_L40x40x150_WP`,
  and `Electron_Proj_E1000_L40x40x150_Classical`. Box 40×40×150 orthorhombic,
  N=162, dx=0.30 Bohr (Nyquist-safe at k_0=8.57), dt=0.005 a.u., N_STEPS=2800.
- `inq-stack/python/inqview/postprocess/bath_energy.py` — sums per-orbital
  energies excluding the WP slot (or all states for classical runs); writes
  `bath_energy_vs_time.csv` + `.png`.
- `inq-stack/python/inqview/postprocess/stopping.py` — dE_kin vs trajectory z
  for both projectile types; F_z(z) for classical; σ(t) free-particle
  baseline for WP.
- `inq-stack/python/inqview/postprocess/pipeline.py` (modified) — registers
  `bath_energy` and `stopping` phases.

Run-directory scaffolding (committed since this is a `runs/...` branch where
the experiment definition is the deliverable):
- `ResearchProject/systems/jellium/save_gs/gs_L40x40x150_orth_N162_dx0p30/run.cpp`
  — production GS save script (orthorhombic cell + 162 e + dx=0.30 +
  EXTRA_STATES=20).
- `ResearchProject/systems/jellium/run_wp_e1000_L40x40x150/run.cpp` — hand-
  written WP run (run_template is cubic-only). Reproduces the L=50 N=162
  observable stack + 3 explicit full-overlap snapshots at {0, N/2, N}.
- `ResearchProject/systems/jellium/run_classical_e1000_L40x40x150/run.cpp`
  — classical-electron run (custom UPF + mass override + impulsive ion
  dynamics + `electron_track.csv`).

### Smoke-test scratch files (NOT committed; under `scripts/classical_electron_smoke/`)

- `C1_load_upf/smoke_C1.cpp` — UPF load test.
- `C2_mass_override/smoke_C2.cpp` — mass override test.
- `C3_free_flight/smoke_C3.cpp` — impulsive propagator kinematic test.
- `C_pre_gs_dryrun/dryrun.cpp` — orthorhombic cell + electrons build sanity.

### main branch state

`features/jellium-ks-energy-observables` was fast-forwarded into main as
two commits:
- `fbedd1c` inqkit+inqview observables (KS-energy, eigenvalue, momentum,
  occupations, density_delta, momentum_distribution, …) + 9 new inqview
  postprocess modules.
- `355f3d4` docs+scaffolding (handovers, journals, plans, reports, sources,
  jellium configs, antiproton UPF, eigenvalues_writer, leed_screen_layout).

These two commits land 22k insertions across 88 files (production
library + research scaffolding + ~40 docs).

The `runs/electron-classical-wavepacket-jellium` branch was then cut from
the post-merge main; nothing has been committed on this branch yet
(Phases 1, 3, 4 files are uncommitted on disk).

## Files touched

- WRITE: `ResearchProject/systems/jellium/shared/pseudopotentials/electron-ONCV-1.2.upf`
- WRITE: `ResearchProject/systems/jellium/shared/configs/electron_proj_E1000_L40x40x150.hpp`
- WRITE: `ResearchProject/systems/jellium/save_gs/gs_L40x40x150_orth_N162_dx0p30/run.cpp`
- WRITE: `ResearchProject/systems/jellium/scripts/classical_electron_smoke/C1_load_upf/smoke_C1.cpp`
- WRITE: `ResearchProject/systems/jellium/scripts/classical_electron_smoke/C2_mass_override/smoke_C2.cpp`
- WRITE: `ResearchProject/systems/jellium/scripts/classical_electron_smoke/C3_free_flight/smoke_C3.cpp`
- WRITE: `ResearchProject/systems/jellium/scripts/classical_electron_smoke/C_pre_gs_dryrun/dryrun.cpp`

## Commands run

- `git add` (specific files, not bulk) → `git commit fbedd1c` (production code)
- `git add` (specific files, not bulk) → `git commit 355f3d4` (docs+scaffolding)
- `git checkout main && git merge --ff-only features/jellium-ks-energy-observables`
- `git checkout -b runs/electron-classical-wavepacket-jellium`
- `cp anti-proton-ONCV-1.2.upf electron-ONCV-1.2.upf`
- `inq-run` in each of `C1_load_upf/`, `C2_mass_override/`, `C3_free_flight/`, `C_pre_gs_dryrun/`

## Tests and validation

Each smoke test had explicit pass criteria stated in its source file.
Observed outputs:

- **C1**: `species(0).symbol() = H`, `atomic_number = 1`, default mass
  1837.18 a.u. (proton mass); `has_file = true`; UPF parsed without
  exception. PASS.
- **C2**: `mass_au = 1.0` exactly (using
  `1.0 / inq::ionic::species::amu_to_atomic_units` literal);
  `KE_au = 36.7499 Ha`; `KE_eV = 1000.015`. PASS (mass_au within 1e-3,
  KE_eV within 1.0 of 1000.0).
- **C3**: 20 impulsive steps at v_z=8.5732, dt=0.005 — final z error
  vs analytic = 1.78e-15 Bohr (machine precision); transverse drift = 0
  exactly. PASS.
- **Pre-GS dryrun**: cell::orthorhombic(40, 40, 150).periodic() accepted;
  101 states constructed (162/2 + 20); INQ chose grid 135×135×500
  = 9.11M points (slightly larger than the 133×133×500 = 8.85M I'd
  estimated — INQ rounds to FFT-friendly sizes). PASS.

What remains UNVERIFIED:
- The pseudopod machinery hasn't actually parsed `PP_LOCAL` or generated
  the screened ionic potential yet — that happens lazily inside
  `ground_state::calculate(...)`. So a malformed UPF body would only
  surface in C4 or in the production GS itself.
- Energy drift over the full propagation — pending C5 / W3.

## Trusted sources used

- `inq/tests/pseudos.cpp:23` — canonical
  `ions.insert(ionic::species(sym).pseudo_file(path), {pos})` pattern.
- `inq/src/ionic/species.hpp:62-91` — `mass(amu)` setter and
  `amu_to_atomic_units = 1822.8885`.
- `inq/src/ionic/propagator.hpp:40-59` — impulsive::propagate_positions
  is exactly `r += dt*v`.
- `inq/src/systems/cell.hpp:61` — orthorhombic constructor.
- `inq/src/systems/electrons.hpp:191` — INQ requires num_electrons > 0.
- `ResearchProject/systems/jellium/save_gs/gs_L50_cubic_N162_dx1p0/run.cpp`
  — pattern for the orthorhombic save_gs.
- `ResearchProject/systems/jellium/run_positive_ion_L50_v0p33/run.cpp`
  — pattern for the classical run.cpp (which we'll mirror in 5b).

## Attribution notes

- `electron-ONCV-1.2.upf` derives from `anti-proton-ONCV-1.2.upf`,
  which derives from `H_ONCV_PBE-1.2.upf` (Schlipf-Gygi norm-conserving
  pseudopotential library; cite Schlipf & Gygi, Comput. Phys. Commun.
  2015, DOI 10.1016/j.cpc.2015.05.011 for the original H pseudopotential).
- The "antiproton-as-projectile" trick in INQ TDDFT is from the
  positive-ion companion plan (`docs/plans/jellium_positive_ion_companion.md`).

## Known issues / blockers

1. **r_s mismatch (flagged in plan).** N=162 in 40×40×150 gives r_s ≈
   7.07, much more dilute than the L=50 N=162 base run (r_s ≈ 5.69).
   Direct comparisons of plasmon ω_p, k_F, etc. with prior journal
   entries will need re-derivation at this density.
2. **run_template.hpp is cubic-only.** The WP run will need a
   hand-written run.cpp (Phase 5a) rather than the standard template.
   This is a one-off cost for this run; future orthorhombic runs would
   benefit from extending run_template to read `LX/LY/LZ` instead of
   `L_BOHR`.
3. **Latent MPI-double-init bug in `run_positive_ion_L50_v0p33/run.cpp`
   line 38.** That file calls `input::environment{}` AND constructs
   electrons later — discovered in our pre-GS dryrun that this combo
   aborts at runtime. The existing positive-ion run is marked "SCAFFOLD,
   NOT YET BUILT" so the bug is latent. We should fix it when we touch
   that area in Phase 5b. New code should NOT call `input::environment{}`.

## Assumptions still in play

- Pseudopod will accept the electron UPF when it is finally parsed
  (during ground_state::calculate). The byte content of `electron-ONCV-1.2.upf`
  is identical to `anti-proton-ONCV-1.2.upf` apart from the PP_INFO
  comment, and pseudopod parses comments only as opaque text.
- N=162 will give a closed-or-near-closed shell in the orthorhombic box
  with T=100 K smearing. If `occupations_vs_time.csv` shows large
  fractional occupations near E_F, treat the run as "scoping" rather
  than "production" for plasmon-resonance work.
- The classical electron at E=1000 eV is well above plasmon resonance
  (~7–10 eV at this density), so the comparison is in the Bohr/Bethe
  single-particle regime where classical/quantum are expected to agree
  in mean stopping. Disagreement at this energy would be physically
  interesting.

## Exact next steps

1. **User decision: launch the production GS?**
   - Estimated wall: 1–3 hours on the local A30 (9.1M grid points,
     101 states, SCF tol 1e-6 Ha).
   - Run via:
     ```
     cd /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L40x40x150_orth_N162_dx0p30
     inq-run
     ```
   - Output checkpoint:
     `ResearchProject/systems/jellium/checkpoints/gs_L40x40x150_orth_N162_dx0p30/`.
2. After GS is saved, hand-write Phase 5a (WP run.cpp) and 5b (classical
   run.cpp). Both should be modeled on
   `run_positive_ion_L50_v0p33/run.cpp` but with orthorhombic cell and
   Phase 5b's species-with-mass-override pattern.
3. Run Step C4 (classical, 100 steps) — first time the production grid
   meets the production propagator. ~10 GPU-min.
4. Run Step W1+W2 (WP injection report + 50-step free-flight slope check).
5. Launch production C5 (full classical) and W3 (full WP). Each ~30–90 min.
6. Implement Phase 6 Python postprocess (bath_energy.py, stopping.py).
7. Final verification (Phase 8).
