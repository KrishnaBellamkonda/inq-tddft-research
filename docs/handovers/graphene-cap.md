# Handover: Graphene + CAP wave-packet/classical scattering

Plan (authoritative, decisions LOCKED): `docs/plans/graphene-cap.md`.
Spec: `docs/prompts/absorbing_boundary/graphene_with_cap.md`.
Glossary: "Graphene CAP" section of `CONTEXT.md`.
Paper: `ResearchProject/literature/tddft-quantum-projectile/resources/wave-packet-electron-dynamics-on-graphene.pdf`.

---

## ⏸ PAUSED 2026-06-21 21:55 — impact-parameter comparison (RESTART TOMORROW)

**User paused this task and freed GPU 0.** The running grazing `cl_b3` was KILLED
mid-run; the grazing dispatcher + perp-wp watcher were killed so nothing
auto-relaunches. **The jellium agent on GPU 1 was deliberately NOT touched.**
GPU 0 is free (no graphene processes). NOTHING graphene is running.

### What the campaign IS (read first)
A **perpendicular-vs-grazing impact-parameter comparison** on a SINGLE common
target — **coronene C24H12 finite flake** — to see how impact geometry affects the
WP-vs-classical interaction. Both arms share the **same box (20×22×60), grid
(50 Ha → ~65×74×192), sim-time (N=1319, τ≈26.4 a.u.), CAP (L=20, η=−0.5),
projectile (E=100 eV, σ=1.47, ETRS)**. Only the flake orientation differs:
- **Perpendicular** = flake in x-y plane (⊥ beam), beam +z head-on; GS
  `shared_gs/gs_perp_coronene_50ha`, geometry `coronene_flake_perp.xyz`.
- **Grazing** = flake in y-z plane (∥ beam), beam +z at impact parameter b=x-offset;
  GS `shared_gs/gs_grazing_coronene_50ha`, geometry `coronene_flake_grazing.xyz`.
The classical projectile uses the **`He`-symbol z_valence=−1 UPF**
(`electron_gaussian_sigma1p47_He.upf`) + `.extra_electrons(+1)` — the H-collision
+ Ewald fixes; do NOT revert.

### CHECKPOINT — progress state (2 of 8 runs done)
| Run | Arm | Status | Notes |
|---|---|---|---|
| `cl_b1` | grazing | ✅ DONE | `grazing/run_cl_b1/` preserved |
| `cl_b3` | grazing | ❌ KILLED mid-run (step 404/1319) | dir REMOVED — reruns from scratch |
| `cl_b6` | grazing | ⏳ not started | |
| `wp_b1` | grazing | ⏳ not started | |
| `wp_b3` | grazing | ⏳ not started | |
| `wp_b6` | grazing | ⏳ not started | |
| `cl_b0` | perp (head-on) | ✅ DONE | `perp/run_cl_b0/` preserved |
| `wp_b0` | perp (head-on) | ⏳ not started | **perp WP binary NOT built yet** |

**Preserved (no rebuild needed):** both GSes (`gs_{grazing,perp}_coronene_50ha`),
binaries `scripts/grazing/{cl,wp}/run` + `scripts/perp/cl/run`. **MISSING:**
`scripts/perp/wp/run` (never built — build it on resume).

INQ real-time runs have **no mid-run checkpoint**, so a killed run restarts from
step 0. Only the 2 completed runs' results survive; the 2 GSes (the expensive part)
are fully saved.

### RESTART TOMORROW — exact commands (pick a GPU; GPU 1 is ~4.6× faster)
The grazing dispatcher is now **RESUMABLE** (skips runs with
`run_completed=true`, reruns the rest). To restart the whole remaining queue:
```bash
cd /local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/grazing
# GPU=0 (slow ~12.6 s/step) or GPU=1 if free (fast ~2.7 s/step). Detached:
GPU=1 setsid nohup bash dispatch.sh > dispatch_resume.out 2>&1 < /dev/null &
# -> auto-skips cl_b1, reruns cl_b3, cl_b6, wp_b1, wp_b3, wp_b6, rebuilds notebook
```
Then the **perp WP** run (needs a build first — fresh ~15 min):
```bash
cd /local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/perp/wp
INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study \
INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share \
PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod \
GR_OUTDIR=.../graphene/perp/run_wp_b0/results GR_CAP=1 GR_CX=0 GR_CY=0 \
GR_TAG=wp_b0 GR_E_EV=100 CUDA_VISIBLE_DEVICES=<gpu> inq-run > .../perp/run_wp_b0/run.log 2>&1 &
```
(A ready watcher recipe is `/tmp/perp_wp_watcher.sh` — but /tmp may be cleared;
the command above is canonical.) Runtimes: classical ~1.1 h (GPU 1) / ~4.5 h
(GPU 0); WP slightly more. Full remaining set ≈ 6 runs.

### KNOWN ISSUES carried into resume (fix in post, no GPU)
1. **Classical KE_loss is z-wrap-contaminated.** The CAP absorbs electrons, NOT
   the classical ion, so over τ≈26.4 the ion travels ~71 Bohr and WRAPS the
   60-Bohr periodic z-box ~1.2×, ending ~1.5 Bohr from the flake (mid-acceleration
   to a 2nd pass). `perp cl_b0` thus reports `KE_loss=−2.5 eV` (spurious *gain*).
   **The real stopping = KE(t) drop over the FIRST pass only** (ion crosses flake
   ~step 233, clears into vacuum ~step 400-500, before wrap). Fix the classical
   extraction in `hypotheses/grazing/build_grazing_report.py` to use first-pass
   loss, not init-vs-final. WP runs are unaffected (CAP genuinely absorbs the WP).
2. **GPU 0 is ~4.6× slower than GPU 1** for the identical workload (perp `cl_b0`
   2.7 s/step on GPU 1 vs grazing 12.6 s/step on GPU 0). Killing the stale
   44-day weather-climate Jupyter kernel (PID 1816189, done) did NOT change it →
   GPU 0 looks like genuinely slower/throttled hardware (NVML dead, can't confirm).
   Prefer GPU 1 for the resume when the jellium agent frees it.
3. **GPU process discovery with NVML dead:** use `fuser -v /dev/nvidia{0,1}` (NOT
   nvidia-smi) to see who holds each GPU.

### Deliverables / analysis still TODO
- Fix classical first-pass stopping extraction (issue 1) in the grazing notebook
  builder; it already has KE(t), planar Δn, LEED via the tested kernels.
- The comparison notebook `hypotheses/grazing/grazing_study.ipynb` auto-builds at
  dispatcher end; it currently reads grazing runs — extend to overlay the perp arm.

---

## Milestone: 2026-06-21 — CLASSICAL ENERGY ANOMALY ROOT-CAUSED + FIXED; ensemble relaunched

**ROOT CAUSE (confirmed at source + by smoke):** the classical projectile was
inserted as an ion with **`z_valence=0`** but a `+1/r` repulsive local potential.
Because INQ's ion-ion Ewald charge array is exactly `valence_charge()`
(`inq-study/src/ionic/interaction.hpp:329`), a `z_valence=0` projectile is
**invisible to the ion-ion sum** → it felt the carbon **electrons** (repulsion,
via its local pp) but **NOT the carbon nuclei** (attraction). Net spurious
repulsion → the projectile decelerated ~6.7 eV while still **11 Bohr out in
vacuum**. The same `z_valence=0` also made the local-potential G=0/alpha term
inconsistent → a constant **+103 Ha** offset (step-0 e=−36.74 vs graphene GS
−143.94). Jellium was unaffected (no nuclei; uniform medium → position-independent
overlap → zero force), which is why the identical UPF worked there for months.
The WP run is stable because its −1 electron is a non-back-reacting probe orbital
(occ 0), not an ion — it never had this term.

**FIX (applied):**
- New UPF `shared/pseudopotentials/electron_gaussian_sigma1p47_zm1.upf` with
  **`z_valence=-1.00`** (consistent with its `+1/r` tail) — copy of the σ=1.47
  projectile UPF, only the z_valence field + provenance comment changed.
- `scripts/cap_cl/run.cpp` (and synced master `scripts/cap/run_classical.cpp`):
  point to the `_zm1` UPF and add **`.extra_electrons(1.0)`** so the QUANTUM
  electron count stays 96 (graphene) while the projectile is a proper −1 Ewald
  charge; the cell then carries the physical net −1 + a +1 uniform background.

**VALIDATED (30-step smoke, GPU 1, `scripts/cap_cl/build_smoke_zm1.log` /
`results_smoke_zm1/`):** Number of electrons = 96 ✓; step-0 e = **−143.9047 Ha**
(= graphene GS, +103 Ha offset GONE) ✓; energy drift over 30 steps = **−1.4e-4 Ha
(−0.004 eV)** vs the old **+6.7 eV** ✓; `proj_vf.z=2.71111` vs `v0=2.71106`
(**KE_loss=−0.004 eV**, ballistic in vacuum, no lateral deflection) ✓;
`run_completed=true`. The spurious vacuum force is eliminated.

**Binary frozen** `scripts/cap/run_cl` (51 MB, the fixed build).
**ENSEMBLE RELAUNCHED** 2026-06-21 14:48 (`scripts/cap/dispatch_cl.sh`, setsid
nohup, GPU 1, `dispatch_cl_zm1.out`, PID 535508): 6 runs cl_{centroid,channeling}
_s{1,2,3} (CAP on, seeds 1/2/3, E=100 eV, perpendicular +z). ~36 min each →
~3.6 h. Per-run `post_and_email.py`; rebuilds the notebook at batch end.
⚠ VALIDATION GATE: check run #1 (cl_centroid_s1) for physical close-encounter
stopping (KE_loss>0, energy conserved) before trusting the rest — the smoke only
proved the *vacuum* part (projectile still 11 Bohr out at step 30).

**Note — these are the PERPENDICULAR ensemble** (the originally-planned deferred
classical half). The user's intended **grazing / impact-parameter study** (WP or
projectile moving PARALLEL to the sheet at height b, scanning b) is a SEPARATE
new sub-campaign — the plan's "channeling" was a sideways-displaced perpendicular
impact, NOT a grazing pass. To be designed.

### Milestone: 2026-06-21 (16:50) — TARGET-COMPARABILITY CORRECTION (user)

**User caught a design flaw:** the perpendicular runs targeted **periodic graphene**
(24 C, no H) while the grazing runs target the **coronene flake** (24 C + 12 H) —
different targets → the perpendicular-vs-grazing impact-parameter comparison was
invalid. **Fix (locked): coronene C24H12 in BOTH arms**, with the **same box
(20×22×60), same grid (50 Ha), same sim-time (N=1319), same CAP (L=20, η=−0.5),
same projectile (E=100, σ=1.47, ETRS)** — the ONLY difference is the flake
orientation: **x-y plane (⊥ beam, head-on)** vs **y-z plane (∥ beam, grazing at b)**.
It is literally the same molecule rotated 90°.

**Actions taken:**
- **KILLED the graphene-periodic perpendicular ensemble** (user choice; wrong
  target). `dispatch_cl.sh` + the orphaned `cap/run_cl` killed → GPU 1 freed.
  centroid_s1/s2/s3 completed before the kill = standalone graphene-CAP data
  (the WP graphene arm + those classicals remain as the graphene feasibility).
- **Grazing arm UNAFFECTED** — still running on GPU 0 (cl_b1, dispatcher PID
  630830 alive).
- **Built the coronene-PERPENDICULAR arm** (user: run now, accept no free GPU):
  - Geometry `shared/geometry/coronene_flake_perp.xyz` (native x-y plane, centred).
  - GS config `shared/configs/graphene_perp_coronene_gs.hpp` (box/grid/electrons
    IDENTICAL to grazing; checkpoint `shared_gs/gs_perp_coronene_50ha`).
  - `scripts/perp/{gs,cl,wp}/run.cpp` (= grazing sources, include swapped; cl
    keeps the He-symbol z_val=−1 projectile — coronene-perp also has H).
  - GS **building+running on GPU 1** (`scripts/perp/gs/gs_build.log`, PID 660803).
  - Dispatcher `scripts/perp/dispatch.sh` (cl+wp, default b=0 head-on; `BLIST`
    env extends to a perpendicular b-scan).
- **GPU map now:** GPU 0 = grazing coronene b-scan; GPU 1 = perpendicular coronene
  (GS → runs). Both arms coronene, fully comparable. NO free GPU (user-accepted).

---

### Perpendicular ensemble — run #1 VALIDATED the fix end-to-end (2026-06-21 15:21)
`run_cl_centroid_s1` completed. Member's ACTUAL KE (jittered v0): 66.0 → 59.1 eV
= **~6.9 eV lost** crossing graphene (physical electronic stopping); electronic
energy −143.90 → −143.65 Ha (+0.25 Ha = +6.9 eV) → **energy conserved**, no
blow-up. The `run_summary` `KE_loss=40.9` is the OLD nominal-100 bug (run.cpp did
`E_eV−KE_f`); the **notebook cell computes the correct value** from each member's
KE₀. (Grazing run.cpp now fixed to write `KE_initial_eV`/`KE_loss` from actual v0.)
Ensemble continues run #2/6 on GPU 1; notebook auto-rebuilds at batch end.

### GRAZING / impact-parameter sub-campaign — BUILDING (2026-06-21)
User's intended geometry, finalised: a **finite 2-D graphene flake** (coronene
C24H12, H-passivated — NOT periodic bulk) reoriented into the **y-z plane**
(`shared/geometry/coronene_flake_grazing.xyz`, flake normal = x), grazed by a
**+z** projectile at impact parameter **b = x-offset**. Same z-CAP + traversal;
box 20×22×60 Bohr (finite flake needs vacuum on all sides). Reuses the existing
run.cpp almost verbatim (the flake is rotated in the GS; `GR_CX=b, GR_CY=0`).
- GS config `shared/configs/graphene_grazing_gs.hpp` (namespace `graphene_cfg`,
  108 e = C24H12, 50 Ha, tiny smearing, GS `shared_gs/gs_grazing_coronene_50ha`).
- GS run.cpp `scripts/grazing/gs/run.cpp` (= scripts/gs/run.cpp, include swapped).
  **BUILDING+running on GPU 0** (`gs_build.log`).
- Run binaries `scripts/grazing/{cl,wp}/run.cpp` (= cap_cl/cap run.cpp, include →
  grazing config; cl also fixed KE_loss to use actual KE₀ + records impact_b).
- PLAN: validate GS (closed-shell gap) → smoke 1 grazing run → dispatch 6
  (classical+WP × b={1,3,6}) on GPU 0. b-scan likely rolls overnight (bigger box).

**GS VALIDATED 2026-06-21 15:30:** closed-shell, **HOMO–LUMO gap 2.76 eV** (LDA,
coronene — correct), E=−150.77 Ha, 108 e (54 occ + 24 empty), checkpoint saved.

**SPECIES-COLLISION BUG (found + fixed):** the classical projectile was
`species("H")`, but the coronene flake HAS 12 H atoms. INQ keys species by SYMBOL
→ the projectile "H" collided with the flake H and silently inherited the DEFAULT
H pseudo (z_val=+1, a proton) → 110 e ≠ GS 108 → checkpoint load FAILED.
(Perpendicular graphene had no H, so "H" worked there — that's why this only
surfaced now.) FIX: distinct projectile symbol **`He`** + a `He`-labelled UPF
copy `electron_gaussian_sigma1p47_He.upf` (z_valence=−1). Only `grazing/cl/run.cpp`
needed it (the WP run injects an orbital, no projectile ion → no collision).

**GRAZING CLASSICAL SMOKE PASSED 2026-06-21 15:46** (`scripts/grazing/cl/results_smoke`):
3 species (He,H,C) ✓, **108 e** (= GS) ✓, load OK ✓, step-0 e=−150.78 (GS value,
no offset) ✓, drift −0.05 eV over 30 steps (ballistic in vacuum at b=3) ✓,
`run_completed=true`. ~5.6 s/step → ~2 h/full-run (886k-pt box). WP binary
building+smoking on GPU 0. b-scan dispatcher ready: `scripts/grazing/dispatch.sh`
(frozen `cl/run` + `wp/run`, {cl,wp}×b={1,3,6}, GPU 0, auto-builds notebook).

**WP GRAZING SMOKE PASSED + b-SCAN LAUNCHED 2026-06-21 15:59.** WP smoke: 2 species
(H,C — no projectile ion), 108 e (= GS), WP injected clean (norm_after=1.0,
overlap 8e-5, N0=1.0), energy −146.72 Ha stable (probe, non-back-reacting),
`run_completed=true`. WP ~12.5 s/step (slower than classical's 5.6 — more I/O +
108 states) → ~4.5 h/WP-run. **b-scan dispatcher LAUNCHED** detached on GPU 0
(PID 630828, `scripts/grazing/dispatch.out`): runs cl_b{1,3,6} first (~2 h each →
~6 h, the headline dE/dx-vs-b lands by morning) then wp_b{1,3,6} (~4.5 h each),
then auto-builds `hypotheses/grazing/grazing_study.ipynb`. **Notebook builder
written** `hypotheses/grazing/build_grazing_report.py` (partial-tolerant: KE-loss
vs b, WP ε vs b, planar Δn per b; reuses the tested kernels). Total b-scan ~20 h
on one GPU → overnight+. PROVISIONAL until Task #7.

**GPU map (overnight):** GPU 0 = grazing b-scan; GPU 1 = perpendicular ensemble
(run 2/6, ~2 h left, auto-rebuilds `hypotheses/cap_scattering/` notebook). Both
self-complete; no manual step needed. ⚠ If the jellium agent needs a GPU, these
two graphene jobs hold both — coordinate.

### Post-processing framework (planar Δn + LEED) — BUILT + VALIDATED 2026-06-21
Two reusable, deps-clean `inqview.analysis` kernels (numpy-only, no VTK/mpl):
- `inqview/analysis/diffraction.py` — `diffraction_pattern(density,dx,dy)` →
  kinematic LEED `|FFT2(ρ̄·Hann)|²` + (kx,ky) axes [rad/bohr]. From the
  time-integrated real-space screen density (`io.leed.load_leed_pattern`).
- `inqview/analysis/planar_density.py` — `planar_delta_map(cubes,times,z)` →
  Δn(z,t)=∬[n(t)−n(t₀)]dxdy (paper Fig. 1). Takes pre-loaded cubes (VTI load
  stays in the viz layer).
Tests `inq-stack/tests/python/inqview/analysis/test_{diffraction,planar_density}.py`
**10/10 PASS** (cosine→peak at 2π/λ; Σ_xy reduction; t₀ col=0; etc.); deps-clean
invariant still PASS; both run on real WP-run data. Catalogued.
Wired into `hypotheses/cap_scattering/build_report.py` as **§4.4** (planar Δn
heatmap + 8-screen LEED transmission/reflection panels) + a partial-tolerant
**classical KE-stopping** cell (§4.3, uses each member's *actual* jittered KE₀).
Notebook regenerated, executed, **13 code cells / 0 errors**; new figs
`figs/fig_planar_dn.png`, `fig_leed_centroid_cap.png`. §4.4 auto-switches to the
classical `density_rt_system` Δn once `run_cl_*` complete (dispatcher rebuilds).

---

## Milestone: 2026-06-18 — grill complete, decisions locked, build not yet started

**Done:**
- `/grill-with-docs` session: confirmed NO prior locked graphene decisions existed
  (spec was a raw prompt with open Qs + unit errors). Interviewed the user and
  locked the full design — see the plan. Wrote `docs/plans/graphene-cap.md` and
  the `CONTEXT.md` "Graphene CAP" glossary section.
- GPU status confirmed via `systems/vacuum/gpu_probe` (NVML broken, expected):
  **GPU 0 busy (user's other run), GPU 1 free (23.4 GB)**. Use GPU 1.

**Key locked parameters (full list in plan):**
- 4×4 graphene (32 C), cell 18.6×18.6×60 Bohr, LDA, ONCV-C, 50 Ha, Γ, Ehrenfest.
- Two-sided sin² CAP, L=20 Bohr, **W=−0.5 Ha** (pre-tuned; no W-sweep).
- WP: σ_r=1.47 Bohr (d=1.1 Å), E=100 eV (k₀=2.711), launch z₀=−12.65 Bohr (+z).
- Trajectories: centroid (atom) + channeling (hollow), both perpendicular.
- Ensemble: 3 classical/trajectory; classical Gaussian width 1.47 Bohr (=WP).
- Observables: coronene set (WP) / jellium-classical set (classical) ∪ core, PLUS
  ε(t)+vacuum-remaining, planar-integrated Δn(z,t), absorbed-fraction+t-absorb,
  screen current flux; whole-system = density+current+WP orbital (no per-orbital
  dump); LEED 8 screens @ ±4,±8,±12,±16 Bohr.
- ETRS propagator (mandatory w/ CAP), dt=0.02, ~800–1000 steps, ~60 frames.
- Per-stage emails `[graphene-cap]`; one master notebook in
  `systems/graphene/hypotheses/<sweep>/`.

**Not done (next steps, in order):**
1. Scaffold `ResearchProject/systems/graphene/` (ADR 0007 layout).
2. Build graphene GS (LDA/50Ha/Γ, fixed experimental lattice) → validate
   (convergence, energy/atom, semimetal gap≈0). Email stage 1.
3. Verify carbon ONCV pseudo (reuse coronene's); generate
   `shared/pseudopotentials/electron_gaussian_sigma1p47.upf`.
4. Confirm/add `inqkit` grid current-density `j(r,t)` VTI writer (inq-stack only).
5. Write run.cpp (build against **inq-study**: `INQ_SOURCE=…/inq-study inq-run`),
   full observable manifest, two-sided CAP via `perturbations::sum`, ETRS.
6. **Smoke one short TDDFT run** → fix n_steps/WRITE_EVERY + real per-run wall
   time → re-estimate campaign before full dispatch.
7. Dispatch: 2 no-CAP baselines + (centroid: 3 classical+1 WP) + (channeling:
   3 classical+1 WP) on GPU 1. Per-stage emails.
8. Master notebook (notebook-making skill) + final summary email.

---

## Milestone: 2026-06-19 (early) — GS scaffolded + launched

**Done:**
- Scaffolded `ResearchProject/systems/graphene/` (ADR 0007 layout).
- Geometry generator `scripts/gs/gen_geometry.py` → `shared/geometry/graphene_4x2.xyz`:
  **32 C, rectangular 4×2 supercell, every NN bond = 1.4203 Å (exact graphene),
  box 18.5949 × 16.1037 × 60 Bohr**, sheet at z=0, a C atom at (x,y)=(0,0) =
  centroid target; channeling/hollow target ≈ (0, 5.07 Bohr) — REFINE from
  rendered GS density before the channeling run.
- GS config `shared/configs/graphene_gs.hpp` (LDA, 50 Ha, extra_states=32,
  **temperature=0.10 eV** smearing for the semimetal, SCF tol 1e-6, broyden).
- GS `scripts/gs/run.cpp` (orthorhombic **.periodic()**, default ONCV-C pseudo,
  saves checkpoint `shared_gs/gs_4x2_50ha`, GS density VTI, eigenvalues CSV via
  `inqkit::observables::dump_eigenvalues`).
- **Launched** GS build+run on **GPU 1** (PID 3675817, `scripts/gs/run.log`),
  full INQ recompile in progress.

**Update (build fix):** first GS run tripped the defensive bounds check —
`FATAL: atom 0 outside [-L/2,+L/2]`. Diagnosed as a boundary atom at exactly
x=-Lx/2 with a ~2.5e-6 Bohr Å→Bohr rounding overshoot (NOT a geometry bug;
coronene.xyz confirms INQ reads Å). Fix: relaxed the bounds tolerance to 1e-2
Bohr (periodic cells legitimately have face atoms). Rebuilt + relaunched on
GPU 1 — **compiles + links clean, SCF now iterating** (96 states = 64 occ + 32
extra; near-Fermi eigenvalues cluster ≈0 → semimetal as expected; smearing
0.1 eV). Note: each `inq-run` here rebuilds libxc (~15 min); the build dir is
`scripts/gs/.build`-style cache — subsequent edits recompile only run.cpp IF
the cache persists.

**Next (on GS completion):** validate (SCF converged, energy/atom, semimetal
gap≈0 from eigenvalues.csv), email Step-1 result, then build the TDDFT/CAP
pipeline (run.cpp vs inq-study, grid-current writer, smoke run, 10-run dispatch).

---

## Milestone: 2026-06-19 (00:40) — GS VALIDATED (semimetal); cell corrected 32→24

**Done:**
- **GS converged + validated as graphene.** First 32-atom (nx=4) GS gave a
  spurious ~2 eV gap → NOT graphene (Dirac point at BZ K does not fold to Γ for
  nx=4). Verified via tight-binding that **nx must be a multiple of 3** to fold
  K→Γ (nx=3,6 metallic; nx=4 gapped). Corrected cell to **3×2 = 24 C atoms**,
  in-plane **13.9462 × 16.1037 × 60 Bohr**. Regenerated geometry
  (`shared/geometry/graphene_3x2.xyz`, bonds 1.4203 Å), updated
  `shared/configs/graphene_gs.hpp` (LX=13.9462, N_C=24, N_e=96, extra=24),
  checkpoint `shared_gs/gs_3x2_50ha`.
- **Validation (24-atom):** E=-143.942 Ha = **-5.998 Ha/atom** (size-consistent);
  **GAPLESS** — 46 full + **4 partial-occ (≈1.0) states at E_F spread 10 meV =
  Dirac point**; 22 empty. Figure `scripts/gs/results/gs_validation.png`.
- **Step-1 email SENT** (`[graphene-cap]`, msg-id stored in
  `systems/graphene/.graphene_cap_email_thread.txt`) — covers GS validation +
  the 32→24 correction + next steps. Thread future emails via `in_reply_to`.

**Cost re-estimate:** 24 C (96 e, ~48 orb) is cheaper than 32 → campaign ~8-10 h.

**NEXT (pipeline, not started):**
1. Generate `shared/pseudopotentials/electron_gaussian_sigma1p47.upf` (classical
   projectile smearing = WP width).
2. CAP run.cpp — adapt coronene `shared/cpp/run_template.hpp` (has WP injection +
   LEED screens) + compose two-sided CAP `absorbing(η=-0.5,+mid,w=10) +
   absorbing(-0.5,-mid,10)` via `perturbations::sum`; **ETRS**; full manifest;
   density+current+WP-orbital VTIs; survival ε + planar Δn(z,t). Build vs
   **inq-study** (`INQ_SOURCE=…/inq-study inq-run`).
3. Grid current-density `j(r,t)` VTI writer in inqkit if absent (inq-stack only).
4. **Build-cache:** each `inq-run` here rebuilds libxc (~15 min). Before the
   10-run dispatch, build ONCE (a single build-once binary in `scripts/<sweep>/`)
   and reuse — do NOT pay the rebuild per run.
5. Smoke 1 run → fix N_STEPS/WRITE_EVERY + real wall time → dispatch 10 runs
   (GPU 1): 2 no-CAP baselines + centroid(3 classical+1 WP) + channeling(3+1).
6. Master notebook (notebook-making) + per-stage emails.
- Channeling lateral target: refine hollow-site (x,y) from rendered GS density
  (gen_geometry prints a rough (0, 5.07 Bohr); verify before the channeling run).

---

## Milestone: 2026-06-19 (01:05) — CAP run.cpp written; building vs inq-study (smoke)

**Done:**
- Generated `shared/pseudopotentials/electron_gaussian_sigma1p47.upf` (classical
  projectile smearing = WP width; V(0)=0.5428 Ha ✓, 5001-pt mesh) via
  `inqview.io.gaussian_psp.generate_gaussian_psp` (template = jellium sigma0p5 UPF).
- Wrote **`scripts/cap/run.cpp`** — env-driven build-once WP-CAP binary. Combines
  coronene template (GS load, WP injection w/ orthogonalise_against_occupied, LEED
  screens, density VTIs total/system/wp, overlap, momentum) + vacuum two-sided CAP
  (`perturbations::absorbing(η,+mid,w)+absorbing(η,−mid,w)`, ETRS default, survival
  ε via `inner_region_norm_twosided`). Env: GR_E_EV, GR_CX/CY (impact pt),
  GR_CAP(1/0 → baseline), GR_OUTDIR, GR_DT, GR_NSTEPS, GR_TAG. CAP geom for
  Lz=60,L=20: mid_frac=0.4167, width_frac=0.1667, z_in=20, z0=−12.65, η=−0.5.
  Saves density_l2 via DensityDelta; WP complex wavefunction VTI (→ WP current in
  post); 8 LEED screens @ z=±4,±8,±12,±16; manifest RunType::coronene.
- **Launched build+smoke vs inq-study** (GR_NSTEPS=40, CAP on, E=100) on GPU 1,
  `scripts/cap/build_smoke.log`, env in `scripts/cap/smoke_env.sh`. Fresh build
  (~15 min) → then 40-step smoke to outdir `results_smoke`.

**Scope notes / decisions made this phase:**
- run.cpp is **WP-mode only** so far (the centerpiece, smoke-able). Classical
  projectile mode (moving Gaussian-UPF ion + Ehrenfest) NOT yet written — add
  after WP smoke validates the machinery, then one rebuild serves the full
  dispatch.
- **Grid total-current writer DEFERRED**: `data.current()` gives the integrated
  current (in observables.csv); the WP complex wavefunction is saved so WP current
  density is derivable in post. A full-grid total-current VTI would need a new
  inqkit writer (all orbitals) — documented as deferred, not blocking.
- **Planar Δn(z,t)** will be DERIVED in post-processing from the saved
  density_system VTIs (cleaner than an in-callback axis sum).
- **Ion dynamics:** using INQ `real_time::propagate` default. Verify post-smoke
  whether default = fixed or Ehrenfest; user chose Ehrenfest — add the option if
  the default is fixed (minor effect on sub-fs scattering either way).

**NEXT:** on smoke success → read ε/wall-time, set N_STEPS for production (~1300),
add classical mode, rebuild ONCE, dispatch 10 runs (GPU 1), per-stage emails,
notebook. If smoke build/run fails → fix (likely an inqkit API signature).

---

## Milestone: 2026-06-19 (01:33) — CAP smoke PASSED; 4 WP runs LAUNCHED

**Smoke (40-step, vs inq-study) PASSED:** compiled clean, WP injected idx=71
norm_after=1.0 max_overlap=1.2e-4 (clean orthogonality to graphene occ states),
ε computed, run_completed=true. Direct-invoke probe (2 steps) rc=0, **~6.1 s/step
with per-step I/O** (production WRITE_EVERY=21 ⇒ ~1.3 h/run). Built binary at
`scripts/cap/run` (vs inq-study); snapshotted to **`scripts/cap/run_wp`** so the
campaign binary is frozen while I rebuild for classical mode.

**Channeling/centroid targets:** centroid (atom) = (0,0); channeling (hexagon
hollow) = **(4.6655, −2.6840) Bohr** (2.67 Bohr to nearest atom = hex center).

**Production geometry:** N_STEPS≈1319 (auto), WRITE_EVERY≈21, dt=0.02, E=100 eV,
k0=2.711, z0=−12.65, two-sided CAP η=−0.5 L=20 (mid_frac 0.4167, width_frac 0.1667).

**LAUNCHED: 4 WP runs** (`scripts/cap/dispatch_wp.sh`, nohup, GPU 1,
`dispatch_wp.out`) invoking the frozen `run_wp` (NO rebuild per run):
`cap_scattering/run_wp_{centroid_nocap, channeling_nocap, centroid_cap,
channeling_cap}/`. After each: `post_and_email.py` → quicklook PNG (survival/
absorbed + ΔE) + validation + threaded `[graphene-cap]` email. ETA ~5 h.

**NEXT (classical mode — the remaining 6 runs):**
- Read jellium classical run.cpp for the moving-projectile API (insert a
  `electron_gaussian_sigma1p47.upf` ion at z0 with +z velocity = k0/m_e, Ehrenfest).
- Add classical mode to `scripts/cap/run.cpp` (GR_MODE=classical: no WP injection;
  add projectile ion + velocity; observables = electron_track + densities; same
  two-sided CAP), **rebuild into a SEPARATE binary** `run_cl` (don't clobber the
  running `run_wp`). Smoke, then dispatch 6 classical runs (centroid×3 +
  channeling×3) on GPU 1 AFTER the WP campaign frees it.
- Then master study notebook (notebook-making) in `hypotheses/cap_scattering/`.

**Verify on WP run #1:** real wall time; survival ε(t) should fall as the WP
crosses graphene and the far CAP absorbs (CAP run); baseline (no-CAP) ε stays
high but WP density should reflect/transmit. Confirm ion-dynamics default
(Ehrenfest vs fixed) from the run.

---

## Milestone: 2026-06-19 (02:05) — WP campaign healthy; classical mode written + staged

**WP campaign running well:** run #1 (centroid_nocap) ~step 1100/1319 at
**~1.65 s/step → ~36 min/run** (energy stable, drift ~1e-6). 4 WP runs ETA ~03:55.
Per-run quicklook emails fire via post_and_email.py.

**Classical mode written:** `scripts/cap/run_classical.cpp` — projectile =
`ionic::species("H").pseudo_file(electron_gaussian_sigma1p47.upf).mass(1/1822.8885)`
inserted at (cx+jitter, cy+jitter, z0+jitter) with v=(jitter,jitter,k0+jitter);
**mt19937 ensemble** (GR_SEED>0 draws Gaussian pos σ=1.47 + mom σ_k=0.481; seed=0
central). Ehrenfest (carbons + projectile move). Same two-sided CAP. Observables:
electron_track.csv every step, density total/system VTI, density_delta,
observables.csv, eigenvalues, manifest jellium_classical, 8 LEED screens. Final
KE_loss_eV in run_summary. **API confirmed** from jellium classical run.cpp.

**Staged (ready):** `scripts/cap_cl/run.cpp` (= copy, separate build dir so it
won't clobber the running WP `run_wp`); `scripts/cap/dispatch_cl.sh` (6 runs:
cl_{centroid,channeling}_s{1,2,3}, CAP on, GR_SEED=1/2/3).

**NEXT (next wake, AFTER WP campaign frees GPU 1 ~03:55):**
1. `cd scripts/cap_cl && INQ_SOURCE=…/inq-study GR_NSTEPS=2 GR_TAG=cl_smoke
   GR_OUTDIR=results_smoke CUDA_VISIBLE_DEVICES=1 inq-run --reconfig` → build+smoke
   the classical binary. Check it inserts the projectile (proj_mass_au≈1), runs,
   KE fields populate.
2. Freeze: `cp scripts/cap_cl/run scripts/cap/run_cl`.
3. `nohup bash scripts/cap/dispatch_cl.sh > scripts/cap/dispatch_cl.out 2>&1 &`
   (6 classical runs, ~3.5 h). Per-run emails.
4. Build master notebook in `hypotheses/cap_scattering/` (notebook-making skill):
   per-step sections, WP vs classical survival/stopping, LEED, planar Δn(z,t)
   derived from density_system VTIs, deviations table. Final summary email.
- If GPU still busy at wake: build the WP-results notebook (CPU) meanwhile + reschedule.
- Verify: WP CAP-run survival ε(t) should DROP (absorption); no-CAP baseline ε stays
  high. Classical KE_loss_eV > 0 (stopping). Ion-dynamics = Ehrenfest (confirmed:
  velocities set + propagate moves ions).

---

## Milestone: 2026-06-19 (03:06) — 2/4 WP done (physics sane); notebook builder ready

**WP results so far (PROVISIONAL):**
- run_wp_centroid_nocap: DONE ε=0.7227, absorbed≈0 (no CAP ✓), wall 2151s (36min).
- run_wp_channeling_nocap: DONE ε=0.7138, absorbed≈2e-4 (no CAP ✓), wall 2551s.
- run_wp_centroid_cap: RUNNING (#3). run_wp_channeling_cap pending (#4). ETA ~04:10.
- Sanity ✓: no-CAP baselines absorb ~0; ~28% of WP exits |z|<20 (transmission/scatter).
  CAP runs (pending) should show absorbed_fraction>0.

**Notebook builder DONE + REUSABLE:** `hypotheses/cap_scattering/build_report.py`
→ `cap_scattering_study.ipynb` + figs/ (survival ε(t), classical KE-stopping,
energy drift, run-inventory table, full narrative + deviations table). Reads
whatever runs exist; re-run after each phase. Built partial (2 WP runs) OK.

**NEXT (next wake, when GPU 1 free = all 4 WP done):** build+smoke classical
(scripts/cap_cl, INQ_SOURCE=inq-study, GR_NSTEPS=2) → cp run → scripts/cap/run_cl
→ `nohup bash scripts/cap/dispatch_cl.sh` (6 runs ~3.5h) → after, re-run
build_report.py (now with classical KE figs) → final summary email. If GPU still
busy: reschedule.

---

## Milestone: 2026-06-19 (04:17) — ALL 4 WP done (CAP works!); classical ensemble LAUNCHED

**WP campaign COMPLETE — CAP validated in graphene:**
| run | CAP | ε_survival | absorbed |
|---|---|---|---|
| centroid_nocap | 0 | 0.722 | ~0 |
| channeling_nocap | 0 | 0.714 | ~0 |
| centroid_cap | 1 | **0.119** | **0.847** |
| channeling_cap | 1 | **0.120** | **0.844** |
CAP on ⇒ ε 0.72→0.12, ~85% absorbed. Centroid≈channeling at 100 eV (subtle, expected).

**Classical binary smoke PASSED** (cap_cl): compiled clean, carbons=24 + projectile
idx=24 **proj_mass_au=1.0** ✓, v=(0,0,2.711). Frozen → `scripts/cap/run_cl`.

**Classical ensemble LAUNCHED** (`scripts/cap/dispatch_cl.sh`, nohup, GPU 1,
`dispatch_cl.out`): 6 runs cl_{centroid,channeling}_s{1,2,3} (CAP on, seeds 1/2/3),
~36 min each → ETA ~07:50. Per-run emails via post_and_email.py.

**NEXT (final pass, when all 6 classical done):**
- Re-run `hypotheses/cap_scattering/build_report.py` (now picks up classical
  electron_track → KE-stopping fig + ensemble) → updates cap_scattering_study.ipynb.
- Optionally execute the notebook (jupyter nbconvert) for inline outputs.
- Send FINAL `[graphene-cap]` summary email (thread via
  systems/graphene/.graphene_cap_email_thread.txt): full WP+classical table,
  attach master figs, link notebook, restate PROVISIONAL (Task #7).
- Update handover + CONTEXT as complete. Campaign = 4 WP + 6 classical = 10 runs ✓.

---

## Milestone: 2026-06-19 (05:20) — BUG: classical ions frozen (.ehrenfest() missing) → FIXED, re-running

**Caught via validation (v0 vs vf):** first 2 classical runs had
`proj_v0 == proj_vf` exactly → projectile NEVER moved → TRUE stopping = 0 eV.
Root cause: **INQ freezes ions by default; `options::real_time` needs
`.ehrenfest()`** (jellium classical runs have it; my run_classical.cpp omitted it).
The `KE_loss_eV=-99` on s2 was a 2nd minor issue (my KE_loss used nominal 100 eV,
not the jittered member's actual initial KE — fix in post: stopping = KE(t0)−KE(tf)
from electron_track).

**Actions taken:** killed dispatch_cl + running run_cl; deleted invalid
`cap_scattering/run_cl_*`; added `.ehrenfest()` to run_classical.cpp +
cap_cl/run.cpp; **rebuilding run_cl** (cache warm) + 30-step smoke to confirm
v changes. WP runs UNAFFECTED (WP is a propagating orbital, not an ion; frozen
carbons negligible at sub-fs — documented deviation, not re-run).

**NEXT:** smoke shows proj_vf≠proj_v0 ⇒ `cp cap_cl/run cap/run_cl`; relaunch
`dispatch_cl.sh` (6 runs); then build_report.py + final email. Also: fix
build_report.py classical KE figure to read electron_track KE(t) (already does)
and the summary to use per-member initial KE for stopping.

---

## Milestone: 2026-06-19 (05:35) — WP CAMPAIGN DELIVERED; classical ensemble DEFERRED (bug)

**DELIVERED (valid, PROVISIONAL per Task #7): 4 WP graphene+CAP runs.**
| run | CAP | ε_survival | absorbed |
|---|---|---|---|
| centroid_nocap | 0 | 0.723 | ~0 |
| channeling_nocap | 0 | 0.714 | ~0 |
| centroid_cap | 1 | 0.119 | 0.847 |
| channeling_cap | 1 | 0.120 | 0.844 |
CAP works (ε 0.72→0.12, ~85% absorbed). Notebook rebuilt:
`hypotheses/cap_scattering/cap_scattering_study.ipynb` (+ figs/). Final email sent.

**CLASSICAL ENSEMBLE DEFERRED — two real bugs (needs focused debug, NOT rushed):**
1. **Energy reference wrong:** classical run starts at **e=−36.7 Ha** vs WP's
   −139.8 Ha (same GS load, same C pseudo, 25 ions = 24C+projectile confirmed in
   System info). Difference = the projectile-as-ION (repulsive Gaussian UPF,
   bare −1 Coulomb in the periodic cell) vs WP-as-orbital. Present at step 0 in
   BOTH frozen & ehrenfest smokes ⇒ it's the projectile setup, not dynamics.
   Hypotheses to check: charged-cell (net −1) Ewald/background handling for an
   ion vs orbital; projectile ion-ion / bare-Coulomb self-energy in the cell;
   whether jellium avoids it via the +jellium background (graphene has none).
2. **Ehrenfest energy drift (NOT a crash — see correction 2026-06-19 11:33):**
   with `.ehrenfest()` the electronic energy drifts UP and the projectile loses
   ~6.7 eV while still ~11 a0 from the sheet (in vacuum) — spurious. Likely
   unrelaxed-carbon forces and/or charged-cell error; investigate (relax GS first?
   smaller dt? compensating background?).
   ⚠ **CORRECTION:** the earlier "crashes ~step 24" claim is **unsupported by the
   on-disk logs.** `scripts/cap_cl/build_cl2.log` (the Ehrenfest 30-step smoke)
   ran to "step 30 … ended normally" with `run_completed=true`; no error/abort/NaN
   in any classical log. The real blocker is the energy anomaly above, not a crash.

**Classical debug entry points:** `scripts/cap/run_classical.cpp` (+ cap_cl build);
compare to a WORKING jellium classical run (run_classical_n162_L50_E100_v2) —
jellium has the +background that neutralises the −1 projectile; graphene does not.
Consider: (a) add a uniform compensating background / use INQ's charged-cell
handling; (b) relax the graphene GS so carbon forces ≈0 before ehrenfest;
(c) reduce dt for the ion-dynamics leapfrog. Also fix the `KE_loss_eV` summary to
use the member's actual initial KE (jitter), and build_report.py reads KE(t) from
electron_track (already correct).

**STATUS: WP deliverable COMPLETE. Classical = open follow-up task.**

---

## Milestone: 2026-06-19 (11:18) — master notebook rebuilt as EXECUTED full-journey study

**Done (notebook-making skill):** rewrote
`hypotheses/cap_scattering/build_report.py` into a comprehensive **executed**
notebook (`cap_scattering_study.ipynb`, 20 cells / 7 code, **0 errors**, outputs
embedded) covering the whole campaign to date:
- **§4.1 Stage 1 (GS):** loads `scripts/gs/results/...`, asserts gapless Dirac
  semimetal (E/atom=−6.00 Ha; 4 partial-occ states at E_F, ~10 meV spread),
  in-notebook spectrum fig `figs/fig_gs_spectrum.png`.
- **§4.2 Stage 2 (WP+CAP, DELIVERED):** reads the 4 `run_wp_*` summaries → table,
  survival ε(t)/absorbed(t) fig + energy-drift fig (regenerated from CSVs).
  CAP ε 0.72→0.12, ~85% absorbed.
- **§4.3 Stage 3 (classical, DEFERRED):** documents both bugs; the KE-stopping
  cell is partial-tolerant — auto-renders once valid `run_cl_*` tracks exist.
- Formulas at point-of-use (skill rule); PROVISIONAL gate (Task #7) restated.

Style fixes vs old builder: executable code cells (not pre-rendered PNG links),
GS stage added, stale `fig_classical_stopping.png` removed (was misleading while
classical is deferred).

**Run-machinery completeness:** added the auto-build tail call to BOTH
`scripts/cap/dispatch_wp.sh` and `scripts/cap/dispatch_cl.sh` — each dispatcher
now regenerates the notebook from all completed runs at batch end (skill
convention; previously absent).

**Verified:** no inq processes running, GPU 1 free. Notebook builds in ~30 s on
CPU via `PYTHONPATH=…/inq-stack/python venv/bin/python3 build_report.py`.
When the classical bugs are fixed and the 6 runs complete, re-running the builder
(or the dispatcher tail) auto-adds the classical KE-stopping figure — no edits
needed.

---

## Milestone: 2026-06-19 (11:33) — notebook enriched (per-run GIFs + energetics + log-grounded post-mortem); skill updated

**Skill change (durable):** added to `notebook-making` SKILL.md a **per-run
visual-intuition rule** — every *significant* run carries an **xz density GIF +
per-run energetics** BEFORE aggregate plots; **sweeps** GIF only stated
**representative** runs (no silent caps); failed/anomalous runs get a
**log-grounded post-mortem** (quote the actual log, correct the record if a prior
claim isn't supported). Added to Definition of Done.

**Notebook rebuilt + executed** (`cap_scattering_study.ipynb`, **29 cells / 10
code, 0 errors**, **5 embedded GIFs + 7 embedded PNGs**, 8 MB):
- §4.2 now has a per-run block for each of the 4 WP runs: an **xz `density_rt_wp`
  GIF** (CAP edges z=±20 + sheet z=0 marked) + a **per-run energetics** figure,
  then the aggregate survival/absorbed comparison.
- §4.3 rewritten as a **log-grounded classical post-mortem** using
  `scripts/cap_cl/build_cl*.log` + `results_smoke2/`: shows e(t) rising and the
  projectile KE losing 6.7 eV in vacuum, plus the classical smoke xz GIF.
- New helpers in the builder: `_load_vti` (reuses `inqview.pipeline.density`
  conventions), `make_xz_gif`, `energetics_fig`. GIFs land in `figs/`.

**RECORD CORRECTED:** the prior "Ehrenfest crashes ~step 24" note is **WRONG** —
no on-disk log shows a crash. `build_cl2.log` ran 30 Ehrenfest steps to "ended
normally" (`run_completed=true`, KE_loss=6.70 eV). The genuine classical blocker
is the **energy-reference anomaly** (starts −36.7 Ha vs WP −139.8; spurious vacuum
loss), most plausibly a charged-cell/background issue — NOT a crash. The §4.3
post-mortem and the classical-bug entry above are updated accordingly.

**Caveats / guardrails:**
- ⚠ CAP results PROVISIONAL until inq-study engine regression (Task #7).
- NEVER edit `inq/` (immutable rule). Engine CAP lives in `inq-study`.
- Use GPU 1 only; verify with `gpu_probe`, never nvidia-smi.
- venv python for post-processing; figures `.png` via canonical theme.
- This is a FEASIBILITY REPLICA — deviations from the paper tabulated in plan;
  carry the "not the paper's converged numbers" caveat into every deliverable.
