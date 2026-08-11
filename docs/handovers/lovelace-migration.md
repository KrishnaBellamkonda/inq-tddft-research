# Handover — Lovelace GPU migration + graphene Phase 3 bring-up

---

## Milestone: 2026-08-11 — lovelace_test WP+classical runs completed; Phase 3 confirmed done

### Current status

Both `lovelace_test` test runs completed successfully on the local Lovelace A30 GPUs.
Phase 3 (twodef_sv WP 9-run sweep) confirmed all-done (finished 2026-08-06 12:46). The
`lovelace_test` campaign validates the new Lz=90 / dx=0.5 / CAP-10-Bohr bilayer geometry
end-to-end before any S(v) production sweep.

### What changed

- `ResearchProject/systems/graphene/scripts/lovelace_test/classical/run.cpp` — **NEW** (created
  this session): full classical twin, direct erf/r force, Ehrenfest lattice, pairwise
  interactions, checkpoint/resume. Adapted from `twodef_sv/classical/run.cpp`; uses
  `spacing(DX*1.0_b)` instead of `cutoff`; adds `GR_GEOM` env override (local path
  fallback via `TDDFT_ROOT`).
- `shared/bin/run-lovelace-test.sh` — **NEW**: three-phase launcher (GS → WP+classical,
  setsid-detached parallel), idempotent GS skip, sequential binary builds to avoid CMake
  race, `INQ_SOURCE` switched to `inq-study` for dynamics.

### Files touched

- `/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/lovelace_test/gs/run.cpp` — GS builder (created previous session, unchanged)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/lovelace_test/wp/run.cpp` — WP dynamics (created previous session, unchanged)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/lovelace_test/classical/run.cpp` — classical dynamics (created this session)
- `/local/data/public/skcb2/tddft/shared/bin/run-lovelace-test.sh` — campaign launcher (created this session)

### Commands run

```bash
# Launcher (setsid-detached, survived session)
LOGFILE=.../lovelace_test/lovelace_test.log
setsid nohup bash shared/bin/run-lovelace-test.sh > "$LOGFILE" 2>&1 &
# PID: 3282970
```

### Run results

| Run | Energy | σ | dt | Steps | Wall | Completed |
|---|---|---|---|---|---|---|
| WP 150 eV, GPU 0 | k0=3.3204 | 4 Bohr | 0.04 | 500 (20 a.u.) | 637 s | 2026-08-11 06:54 |
| Classical 25 eV, GPU 1 | k0=1.3555 | 4 Bohr (σ_pot=2.828) | 0.04 | 500 (20 a.u.) | 460 s | 2026-08-11 06:51 |

WP output: `.../lovelace_test/wp/results/bi_E150_sigma4_dx0p5/`
Classical output: `.../lovelace_test/classical/results/bi_E25_sigma4_dx0p5/`

Per run: observables.csv, interactions.csv, ions_track.csv (WP adds wp_momentum_stats.csv,
wp_real_space_stats.csv; classical adds electron_track.csv). VTI density frames (51 per
type at SAVE_EVERY=10; 274 MB WP, 170 MB classical).

### Tests and validation

- Proposed / run: GS end-to-end (bilayer, 48 C, 192 e⁻, Lz=90, dx=0.5); WP 150 eV +
  classical 25 eV propagation (500 steps / 20 a.u.). Both `run_completed = true`.
- WP interactions closure step 0: e_ss + e_pp + e_ps = -2011.502 == energy_hartree ✓
- Classical interactions step 0: e_hartree_check matches e_ss (correct for classical) ✓
- WP momentum stats: pz_mean=3.3204 (=k0 ✓), e_kin_ha=5.619 Ha = 152.9 eV, drift KE=150 eV ✓
- Classical electron_track step 0: vz=1.3555, ke=0.9187 Ha = 25.0 eV ✓; forces balanced
  (Fel=-0.48, Fcore=+0.49 → net +0.01 Ha/Bohr small net attraction at z=-19 as expected).
- Remaining gaps: no post-processing / S(v) extraction done yet for this run pair.

### Known issues / blockers

- Classical `e_external_check = 0` in interactions.csv at all steps (should equal E_SB+E_PS).
  Likely because `moving_gaussian_projectile_potential` contributes to a different INQ energy
  ledger term than `energy_external`. The pairwise e_ps term itself passes the Hartree
  closure, so the physical decomposition is intact; only the external closure check is
  anomalous. Needs investigation before using e_ps as the stopping-power proxy.
- Disk at 99% full (74 GB free on /local/data). VTI frames are 445 MB per campaign run;
  must monitor before any production S(v) sweep.
- Phase 3 notebooks watcher (`run-gr2d-p3-notebooks.sh`, PID 1747791): not confirmed
  still running; all 9 runs are done so the watcher should have finished and emailed.
  Verify `p3_notebooks.log` for completion / email sent.

### Assumptions still in play

- `sigma_WP/√2` → `sigma_pot` convention (σ_pot=2.828 for σ_WP=4) is baked into classical
  run.cpp at compile time.
- Impact parameter (2.3244, 1.342) is a hollow/near-C site in the bilayer. Not explicitly
  verified to be a specific crystallographic site.
- inq-study required for both WP and classical dynamics (both use `perturbations::absorbing`
  CAP).

### Exact next steps

1. Check `p3_notebooks.log` to confirm Phase-3 notebook watcher finished; if email wasn't
   sent, relaunch `run-gr2d-p3-notebooks.sh` (idempotent via `.nb_state/` markers).
2. Investigate `e_external_check = 0` anomaly in classical `interactions.csv`: compare
   `energy_external` from `observables.csv` vs `e_ps` from `interactions.csv`; check sign
   convention in `compute_coulomb_direct` (the jellium reference uses +1 charge proton).
3. Write a quick Python analysis (`hypotheses/lovelace_test/build_lovelace_test.py`) to:
   - Extract S(v₀) from classical `electron_track.csv` (early window where vz ≥ 0.85·v₀)
   - Plot WP norm(t) from `interactions.csv` norm_wp column (CAP absorption curve)
   - Compare pairwise e_ps time-series for classical vs WP (stopping mechanism)
4. Design production S(v) sweep: pick 4–6 energies (e.g. 25/50/100/150/200/300 eV) for
   both classical and WP variants; adapt `run-lovelace-test.sh` into a multi-energy
   dispatcher before launching.
5. Commit new artefacts: `lovelace_test/classical/run.cpp`, `run-lovelace-test.sh`,
   `lovelace-migration.md` update (forbidden-word check required).

---

Rolling handover for resuming GPU work on the **lovelace** workstation after the
CSD3 MPhil allocation was exhausted. Repo: `/local/data/public/skcb2/tddft`,
branch `quantum-stopping-power`. Ground states + analysis CSVs were rsynced from
CSD3; nothing else transfers. This file was created 2026-08-06 (the brief
referenced `docs/handovers/lovelace-migration.md`, which did not previously
exist — the pre-existing related handover is
`docs/handovers/device-migration-branch-packaging.md`, about git branch
packaging, not runs).

## Machine gates — ALL VERIFIED 2026-08-06 (host: lovelace)

| Gate | Result | Consequence |
|---|---|---|
| `nvcc --version` | **CUDA 12.9** (≥ 12.4) | **No** CUB `fixed_return` port needed. |
| GPU arch (CUDA probe, `cudaGetDeviceProperties`) | **2× NVIDIA A30, cc 8.0 = sm_80**, 23.5 GB each | `INQ_CUDA_ARCH=80` (A100 default) is correct — unchanged. Same Ampere GA100-class arch as A100. |
| NVML | `nvidia-smi` fails "Driver/library version mismatch" (drv 535.309) | Monitoring only; **CUDA compute works** (probe returned `device_count=2`). Do NOT fall back to CPU. Use `fuser /dev/nvidia{0,1}` + `ps` to check GPU occupancy, not nvidia-smi. |
| GPUs free | `fuser /dev/nvidia0/1` empty; only tl666 CPU `mpirun` runs | Both A30 available. |
| `/local/data` | 147 GB free (98% full) | Graphene fits locally. lz_bulk_sweep (~450 GB) does NOT — symlink to `/data/phy-damysus` (5.9 TB free) before launching it. |
| Ground states | `graphene/shared_gs/gs_3x2_{mono,bi}_cut50_Lz80` present | Phase 3 GS ready. |

## Env to build/run here (bypasses CSD3 Slurm + Environment Modules)

The `shared/bin/run-*.slurm` files are the **authoritative record of each run's
env vars**, NOT executable here (they `module load` + `#SBATCH`). Use `inq-run`
(header-only INQ → compiles `run.cpp` against `INQ_SOURCE` headers directly).
`shared/config.sh` is already device-portable (uses `$(command -v nvcc)` first).
Key exports for the graphene (CAP) runs:

```
INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study      # CAP runs need inq-study; stock inq can't compile CAP
INQ_DEPS_CACHE=/local/data/public/skcb2/tddft/inq/build/_deps   # inq-study/build/_deps absent → reuse upstream (same dep versions) to avoid network clone
INQ_SHARE_PATH=.../inq/install/share
PSEUDOPOD_SHARE_PATH=.../inq/install/share/pseudopod
```

## Priority 1 — graphene Phase 3 (twodef_sv WP), 9 runs, 1 run/GPU × 2 GPUs

- Run dir: `ResearchProject/systems/graphene/scripts/twodef_sv/wp/` (`run.cpp`,
  builds to `./run`). Matrix from `shared/bin/run-gr2d-p3.slurm`:
  - mono E = 15/25/50/100/200/300 eV  (K0 = 1.0500/1.3555/1.9170/2.7110/3.8340/4.6955)
  - bi   E = 15/50/200 eV             (K0 = 1.0500/1.9170/3.8340)
  - bi E25/E100/E300 done in Phase 2 (results dirs `E25/E100/E300`, old naming) — excluded.
- Fixed env per run: `GR_VARIANT`, `LJ_SIGMA=2.0`, `LJ_K0`, `LJ_LAUNCH_Z=-12`,
  `LJ_DT=0.025`, `LJ_CAP_ETA=-1.0`, `LJ_CAP_L=20`, `LJ_OUT=<V>_E<EV>`,
  `LJ_GS_DIR=shared_gs/gs_3x2_<V>_cut50_Lz80`. N_STEPS auto =
  `ceil(1.5·(|LAUNCH_Z|+LZ/2)/K0/DT)` = `ceil(3120/K0)` (LZ=80). mono=72 states,
  bi=144 states (bi ~2× cost/step).
- **Dispatcher: `shared/bin/run-gr2d-p3-local.sh`** (local, non-Slurm analog of
  run-gr2d-p3.slurm). flock'd job queue (LPT order) + 2 workers pinned to GPU
  0/1; calls `./run` DIRECTLY (not inq-run) so the GPUs don't race on `build/`.
  Idempotent (skips `run_completed = true`) + resumable (`LJ_RESUME=1`). Launch:
  `setsid nohup bash shared/bin/run-gr2d-p3-local.sh > <log> 2>&1 &`.

## Priority 2 — lz_bulk_sweep, 30 runs, ~39 GPU-h (NOT STARTED)

Output ~450 GB > local 147 GB free. Before launching, symlink
`lz_bulk_sweep/{wp,classical,vac}/results` → `/data/phy-damysus`. Env records in
`shared/bin/run-lzb-*.slurm`. Graphene stays local.

## MILESTONE — 2026-08-06 06:30: 9-run sweep LAUNCHED on both GPUs

- Migration blocker fixed: `twodef_gs.hpp` baked a CSD3 `/rds/...` geometry path.
  Added `GR_GEOM` env override in `wp/run.cpp:99` (portable; mirrors `LJ_GS_DIR`);
  dispatcher passes local `shared/geometry/graphene_3x2{,_bilayer}.xyz`. Rebuilt
  binary (06:25). Classical + gs run.cpp still carry the baked path — wire `GR_GEOM`
  there on next touch.
- Smoke timing (20-step): mono 2.97 s/step, bi 5.44 s/step (setup-inflated). LIVE
  per-step is lower: bi ~4.0 s/step → revised total ~12 GPU-h, ~6 h wall on 2 A30.
- Grid 45×54×256 = 622,080 pts; ~1.2 GB frames/run → ~11 GB for 9 runs (fits 146 GB
  local; disk guard still recommended).
- Dispatcher `run-gr2d-p3-local.sh` running (setsid-detached, host lovelace),
  LPT order: GPU1=bi_E15, GPU0=bi_E50 first; 7 queued. Log:
  `.../twodef_sv/wp/p3_dispatch.log`. Both jobs healthy (finite steady energy).
- User decision (2026-08-06): AUTO-CHAIN notebooks — per-run notebook (density
  GIFs + ledger) as each run finishes, then Phase-3 study notebook + ledger + ratio
  after all 9, then email. Notebook watcher = CPU-side, must NOT steal GPU.

## STATE (2026-08-06 06:39) — Priority 1 FULLY AUTONOMOUS, running

Two detached daemons on lovelace (survive session end):

1. **Physics** — `run-gr2d-p3-local.sh` (pid 1740819). 2 A30, 1 run/GPU, LPT
   queue. At snapshot: GPU1=bi_E15 (step ~124), GPU0=bi_E50 (step ~110), 7
   queued. Log `.../wp/p3_dispatch.log`. ~6 h wall projected.
2. **Notebooks** — `run-gr2d-p3-notebooks.sh` (pid 1747791), CPU-only. Polls
   every 90 s; builds each `<OUT>_run_notebook.ipynb` (density GIFs + energetics)
   as its run completes (markers in `hypotheses/twodef_sv/.nb_state/`); after all
   9, runs `build_phase3.py` (→ `phase3_ledger.csv`, `phase3_sv.png`,
   `phase3_study.ipynb`) and emails the 4-part summary with the S(E) plot. Log
   `.../wp/p3_notebooks.log`. Disk guard (<15 G → warn+email) fixed.

VALIDATED (not compile-only): geometry fix ran end-to-end (`GEOM_FIX_OK`);
`build_phase3.py` built a sane partial ledger from the 3 existing Phase-2 bi runs
(ratio 25→300 eV: −0.46/0.58/2.81, executed nbclient 0-err); `run_notebook_builder.py`
built E100 notebook (type=wp, 29 cells, 0 errors) with the exact watcher command.

New tracked artefacts: `shared/bin/run-gr2d-p3-local.sh`,
`shared/bin/run-gr2d-p3-notebooks.sh`,
`hypotheses/twodef_sv/build_phase3.py`, `wp/run.cpp` GR_GEOM override.

### Resume / control
- Kill physics: `kill 1740819` (+ its `./run`). Resume/extend: re-run dispatcher
  with `LJ_RESUME=1` (idempotent; skips completed, resumes partial from ckpt).
- Kill notebooks: `pkill -f '[r]un-gr2d-p3-notebooks'`. Relaunch: same command
  (idempotent via `.nb_state/` markers). Do NOT `pkill -f run-gr2d-p3-notebooks.sh`
  from a shell whose own cmdline contains that string — it self-kills (exit 144).

### NOT DONE
- lz_bulk_sweep (Priority 2): 30 runs ~39 GPU-h; symlink
  `lz_bulk_sweep/{wp,classical,vac}/results` → `/data/phy-damysus` (5.9 TB free)
  BEFORE launching (~450 GB > local). Env in `shared/bin/run-lzb-*.slurm`.
- Classical + gs graphene run.cpp still carry the baked `/rds/...` geometry path —
  add `GR_GEOM` there before any classical Phase-3 twin run.
- Commit the new artefacts (commit-messages rule; forbidden words checked).

## Gotchas carried from memory/rules

- Light free-Ehrenfest WP decelerates — but run.cpp sizes N_STEPS by traversal
  (1.5·distance/v), which is correct here; the light-projectile rule governs S
  *extraction*, not this run sizing.
- Frames fill sda1 → at 100% full the Bun runtime SIGBUS-crashes Claude. Watch
  `df -h /local/data` during the 9 runs (each saves ~50 frames).
- Background runs must be `setsid`-detached to outlive the session (harness
  SIGKILLs `run_in_background` on session end).
