# Handover: electron-classical-wavepacket-jellium comparison

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
