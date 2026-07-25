# Handover: N=162 localised-jellium mass-pair WP runs (σ=1, E=100 eV, CAP η=−1.0)

Rolling handover. Plan: `docs/plans/mass-pair-n162-sigma1-cap.md`.
Branch: `overnight-gaussian-classical`. Started 2026-07-19.

## Update: 2026-07-21 — runs KILLED to free GPUs (supersedes the pause below)

The suspend below was escalated to a full teardown at the user's request to
**free GPU memory**. `SIGKILL`ed (in order) the supervisors then the runs+orted:
PIDs 2616425, 2687871 (supervisors), 2678810, 2684101 (runs), 2678878, 2684166
(orted). **Both A30s verified free: 23815/24062 MB each** (`fuser` shows no
holders; NVML-independent `gpu_probe`). The run PIDs in the table below are now
DEAD — do **not** `SIGCONT` them.

**Resume is now checkpoint-based only** (the in-memory steps past the checkpoints
are gone): re-launch each run with `EM_RESUME=1` + the desired `EM_N_STEPS` —
**m1 resumes from step 1500, m2 from step 1000**. m1 loses the 3 in-flight steps
(1500→1503), m2 loses 71 (1000→1071); both recompute on resume. On-disk
checkpoints intact and verified. The two notebook builders (3527359, 3527358)
remain SIGSTOP'd and hold NO GPU memory (CPU-only); resume or kill at will.

## Milestone: 2026-07-21 — BOTH RUNS PAUSED (user-directed), suspended mid-flight

### Current status
**Campaign PAUSED, not stopped or completed.** Both WP runs and both their
supervisors are suspended (`SIGSTOP`, process state `T`) as of 2026-07-21 ~20:07.
Zero data loss — the runs are frozen in memory at their live step; `SIGCONT`
resumes exactly where they left off. Reason: building the two run-notebooks on
the same host saturated the CPU (matplotlib density-GIF battery, ~900 VTI frames
× 9 kinds/run) and starved the runs' host threads, stalling both (m1 dropped to
~75 s/step then 0; m2 absorbed a 946 s single step). Deprioritising the builders
(`nice 19` + `ionice -c3`) was insufficient, so the builders were paused; the
user then directed pausing the runs themselves.

Paused positions (in-memory, live):
- **m1** (mass 1, GPU0): **step 1503 / 2500**, t = 60.12 a.u. Healthy — energy
  −53.1344 Ha, ~46 s/step before pause.
- **m2** (mass 2, GPU1): **step 1071 / 2500**, t = 42.84 a.u. Healthy — energy
  −52.90 Ha, ~50 s/step before pause. (m2 trails: it crashed once at launch
  2026-07-19 15:51 rc=−6, was resumed by `m2_watch.py`.)

On-disk RT checkpoints (the fallback if the suspended procs are ever killed):
- **m1**: `last_step=1500` (written 2026-07-21 20:05:21; 78 orbital state files +
  spin_density verified). Confirmed live: watched `rt_state.txt` flip 1000→1500
  as m1 crossed step 1500, and the run continued cleanly to 1503.
- **m2**: `last_step=1000` (written 2026-07-21 15:24:28). Next disk checkpoint due
  at step 1500 (~429 steps away when paused).
- Checkpoint interval **500 steps** (`EM_CKPT_EVERY=500`) — adherence verified.

### Suspended process inventory (all state `T`)
| PID | what | paused at |
|---|---|---|
| 2678810 | m1 `wp/run` (GPU0) | step 1503 |
| 2684101 | m2 `wp/run` (GPU1) | step 1071 |
| 2616425 | `orchestrate.py` (m1 supervisor) | — |
| 2687871 | `m2_watch.py` (m2 supervisor) | — |
| 3527359 | m1 `run_notebook_builder.py` | mid GIF-battery |
| 3527358 | m2 `run_notebook_builder.py` | mid GIF-battery |

Supervisors were frozen FIRST so neither can misread a suspended run as crashed
and relaunch a competitor on the same GPU (both relaunch only on process *exit*,
but freezing them makes the pause unambiguous). GPUs remain allocated (SIGSTOP
holds device memory); they are NOT freed for other users while suspended.

### To RESUME the campaign exactly (zero loss)
```bash
# order: supervisors first, then runs (so a watcher sees a live, advancing run)
kill -CONT 2616425 2687871      # orchestrate.py, m2_watch.py
kill -CONT 2678810 2684101      # m1, m2 runs  -> continue from steps 1503 / 1071
# (leave the notebook builders 3527359/3527358 paused until the runs finish)
```
If the suspended run processes are ever KILLED (reboot / manual kill) instead of
continued, resume from the last DISK checkpoint via the run's resume branch
(`EM_RESUME=1`, larger `EM_N_STEPS`) — m1 restarts from step 1500 (loses the 3
in-flight steps), m2 from step 1000 (loses 71 in-flight steps). Reference:
`.../scripts/mass_pair_n162/orchestrate.py` `launch(..., resume=True)`.

### Notebooks (deferred, paused mid-build)
Two run-notebooks were being (re)built on the current partial data at
`.../hypotheses/mass_pair_n162/{m1,m2}_run_notebook.ipynb` with full config flags
(rs=5.686, cap-inner=21.2, proj-sigma=0.70711, launch-z=−16.5, v0=2.711/1.917,
e-gs-ha=−53.8388, l-slab=25). They are paused mid GIF-render. Recommended: leave
paused until the runs finish, then `SIGCONT` (or rebuild) so notebooks reflect
near-final data and do not contend with the runs.

### Exact next steps
1. Decide notebooks-now (accept run slowdown, `SIGCONT` builders) vs runs-first
   (keep builders paused, resume runs).
2. To continue the science: `SIGCONT` the four run/supervisor PIDs above.
3. When runs complete (m1 ~46 s/step; m2 ~50 s/step), finish the notebooks on
   the idle host.

## Milestone: 2026-07-19 — autonomous pipeline launched, GS building

### What the user asked for
Two matched **quantum wavepacket** projectile runs through a **genuine
162-electron localised jellium slab** (z=±12.5, smoothening edge_width=1.0),
one per GPU, identical except the projectile mass:
- **m1** — projectile mass 1 (INV_MASS=1.0), GPU 0.
- **m2** — projectile mass 2 (INV_MASS=0.5), GPU 1. **Bath electrons stay mass 1
  in both** (per-orbital `electrons.inverse_mass()[0][wp_idx]` — inq-study fork).
σ_WP=1, E=100 eV, two CAPs 10 Bohr/side η=−1.0, sim time 100 a.u., dx=0.40
(user: finest allowed). Full energy decomposition each timestep + extensive
observable suite. **Checkpoint every 500 steps.** User stepped away 2026-07-19
02:33; "paramount that these runs execute" → full autonomy required.

### Design decisions (all user-confirmed)
- **Genuine 162 electrons** (NOT density-matched) → NEW GS required (none existed).
  Box **70.4×70.4×62.4 Bohr, dx=0.40 → 176×176×156 grid**, n0=1.305e-3, r_s=5.68.
- Lz=62.4 from the 4σ rule: WP launch z=−16.5 (4σ from slab face), CAP inner
  face ±21.2 (4.7σ from WP), CAP region ±[21.2, 31.2] (10 Bohr), outer face=Lz/2.
- k0: m1=2.711, m2=3.834 (=√(2·m·E)). dt=0.04 → 2500 steps (precedent:
  effmass_sigma1). Memory est ~19.5 GB (conservative); both A30s 23.8 GB free.
- Aliasing gate PASSED for both masses at dx=0.40 (0.00% tail).
- Plain (unchirped) launch — the chirp needs dx=0.333; dx=0.40 forbids it.

### Files created (all NEW, on branch)
- `ResearchProject/systems/localised_jellium/shared/configs/slab_n162_L70x70x62.hpp`
- `.../scripts/mass_pair_n162/gs/run.cpp` — GS builder
- `.../scripts/mass_pair_n162/wp/run.cpp` — WP run (fullsuite obs + mass fork +
  checkpoint/resume + per-step energy_decomp.csv + plain launch + 10-Bohr CAP)
- `.../scripts/mass_pair_n162/orchestrate.py` — autonomous orchestrator
- `docs/plans/mass-pair-n162-sigma1-cap.md`

### What is RUNNING now
1. **GS build+run** (I launched it directly): background bash task, GPU 0,
   `INQ_SOURCE=inq-study`. Log:
   `.../scripts/mass_pair_n162/gs/gs_build_run.log`. Compiling inq-study (~75% at
   02:33). Writes checkpoint `shared_gs/slab_n162_L70x70x62_dx0p40/`.
2. **Orchestrator** PID 2390495 (detached, `nohup`), log
   `.../scripts/mass_pair_n162/orch.log`. Autonomously: waits GS → builds WP
   binary (inq-study) → 24-step pilot gate → launches m1(GPU0)+m2(GPU1) →
   auto-resumes on crash (≤3 tries, EM_RESUME=1 from last 500-step ckpt) →
   best-effort analyse → emails chiddukanna@gmail.com at every milestone/failure.

### Autonomy gates (what could still stop it — being watched)
- GS SCF must converge (energy finite). Orchestrator checks; times out at 4 h.
- **WP binary must compile** (inq-study; ~15 min first build). Highest residual
  risk — new merge of fullsuite+effmass+per-step energy writer. On fail →
  "WP BUILD FAILED" email, STOP.
- **Pilot gate** (correctness only, per checkpoint-dont-block): blocks ONLY on
  crash/OOM/NaN/WP-norm∉[0.9,1.1]; WARNS (proceeds) on energy drift. On hard
  fail → "PILOT FAILED" email, production NOT launched.

### How to check status (for the user / next session)
```bash
cat .../scripts/mass_pair_n162/orch.log                    # orchestrator timeline
tail .../scripts/mass_pair_n162/gs/gs_build_run.log        # GS build/SCF
ls   .../shared_gs/slab_n162_L70x70x62_dx0p40/             # GS checkpoint present?
tail .../scripts/mass_pair_n162/wp/rt_m1.log  rt_m2.log    # production progress
cat  .../scripts/mass_pair_n162/wp/results/{m1,m2}/rt_ckpt/rt_state.txt  # last_step
grep run_completed .../wp/results/{m1,m2}/run_summary.txt  # done?
ps -p 2390495                                              # orchestrator alive?
```
Manual resume of a killed run (from cwd `.../mass_pair_n162/wp`, INQ_SOURCE=inq-study):
`EM_OUT=m1 EM_INV_MASS=1.0 EM_K0=2.711 EM_RESUME=1 CUDA_VISIBLE_DEVICES=0 ./run`

### Status: done / partial / not done
- DONE: config, gs/run.cpp, wp/run.cpp, orchestrator, plan; GPUs verified free;
  aliasing gate passed; email verified working; orchestrator launched.
- PARTIAL (in flight): GS compiling; WP not yet built; pilot not yet run.
- NOT DONE: production runs; analyse.py (optional, orchestrator skips if absent —
  28 h runway to add it); catalogue/notebook/journal; run-notebook density GIFs.

### Verified vs unverified
- VERIFIED: geometry/4σ math; aliasing gate; API names (InjectionReport,
  energy().*, jellium::eigenvalues, writer ctors) against compiling templates;
  inq-study build path; email import; GPU free (23.8 GB each).
- UNVERIFIED: WP binary COMPILES (biggest risk); memory actually fits at runtime;
  energy drift at dt=0.04 on this grid; 162 closed-shell for a slab (smearing handles).

## Milestone: 2026-07-19 06:14 — pilot CRASHED (step-0 GPU illegal access); diagnosing

GS converged fine (E=−73.106 Ha). WP binary compiled (0 errors). WP injects OK
(norm=1.000000, max_overlap=2.7e-4, wp_idx=98). Cost measured: **~20.5 s/step ⇒
~14 h per 2500-step run** (better than the 28 h estimate). BUT the 24-step pilot
hit `CUDA ERROR: an illegal memory access` at the FIRST propagation step (t=0
energy row written, then abort). Orchestrator correctly did PILOT-FAILED email +
STOP (no production launched). Orchestrator process has exited.

**Hypothesis (ranked):** step-0 GPU-memory exhaustion from a fullsuite-only
observable that was only ever validated on smaller CUBIC boxes. #1 =
`momentum_distribution` (its `accumulate` does `to_fourier` on ALL 99 states →
~7.85 GB transient on top of ~16–20 GB base → >24 GB). effmass (mass fork + CAP +
density_delta, orthorhombic) ran fine but used NONE of the extras.

**Action taken:** added env toggles to `wp/run.cpp` (default ON) to bisect —
`EM_OBS_MOM`, `EM_OBS_WF`, `EM_OBS_OVL` (unique_ptr-guarded), `EM_OBS_STATE`,
`EM_OBS_DIPCUR` (leave DIPCUR on — coupled to RealTimeSession). Recompiled OK.
Running 6-step bisect: GPU0 = MOM off; GPU1 = all four extras off (core baseline).
Logs `wp/rt_dbg_momoff.log`, `wp/rt_dbg_core.log`.

**Likely fix:** if MOM-off runs, drop `momentum_distribution` from production
(keep `wp_momentum_stats`, which FFTs only the single WP orbital → cheap, gives
⟨p⟩/σ_p) and, if needed, also disable the complex WF VTI. Then re-pilot and
relaunch the orchestrator. The per-step energy_decomp + density VTIs + wp stats +
overlap remain — still an extensive suite.

## Milestone: 2026-07-19 ~07:00 — ROOT CAUSE found + fixed

Bisect verdict (6-step tests): the crash is **GPU memory exhaustion from enabling
`.observables_current().observables_dipole()`** in `rt_opts`. INQ computes the
current density (gradients over all 99 states) each step → another full-field
temporary on top of the ~16–20 GB ETRS peak → >24 GB → step-0 illegal access.
- `dbg_nodipcur` (current/dipole OFF, all extras off): **RUNS clean**, energy
  conserved (E_total −68.61970→−68.61968 over 2 steps, drift ~1e-5 Ha).
- `dbg_nocap` (CAP off, dipcur ON): still crashes → CAP is NOT the cause.
- `dbg_momoff` and `dbg_core` (both dipcur ON): crash → extras are NOT the cause.
`momentum_distribution` is ALSO disabled for production (its accumulate FFTs all
99 states → ~7.65 GB transient — the other memory bomb).

**FIX applied:** `orchestrate.py` BASE now sets `EM_OBS_DIPCUR=0 EM_OBS_MOM=0`
(WF/OVL/STATE stay on). Kept: per-step full energy decomposition (energy_decomp.csv,
the user's key ask — needs no current/dipole), density VTIs (total/system/wp/delta),
wp_momentum_stats, wp_real_space_stats, state_energies, occupations, overlaps.
DROPPED: J(t)/μ(t) current-dipole time series + the all-states momentum histogram.
dt=0.04 verified stable. Cost ~20 s/step ⇒ ~14 h/run.

**Toggles added to `wp/run.cpp`** (env, default ON): EM_OBS_MOM/WF/OVL/STATE/DIPCUR.

## Milestone: 2026-07-19 ~08:10 — refined: base propagation is memory-bound

Further bisect showed it is NOT just current/dipole. `dbg_lean2` (DIPCUR=0 MOM=0
OVL=0, but WF=1 STATE=1) STILL crashes at step 0, while `dbg_nodipcur` (EVERYTHING
off) runs. Conclusion: **the base 99-state ETRS propagation sits near the 24 GB
A30 ceiling**, so EVERY all-states observable independently OOMs — current/dipole
(current field), momentum (FFT all states), overlaps (GS reference copy), AND
state_energies (H applied to all 99 states, ~7.65 GB). `dbg_prod` (OVL on)
confirmed GPU free=0 MB (fully saturated). Only the always-on core fits.

Also learned: with OVL off, setup is FAST (~10 s, no 99×99 t=0 overlap snapshot) —
so a step-0 abort ~80 s after launch is still the first propagation step.

Also hit a **CUDA teardown race**: launching a new run on a GPU right after
`kill -9`-ing the previous one starts with the old context's ~18 GB still held →
false OOM. ALWAYS confirm `gpu_probe` shows ~23800 free before launching.

**PRODUCTION observable set (locked in orchestrate.py BASE):** DIPCUR=0 MOM=0
OVL=0 WF=0 STATE=0 — the confirmed-working minimal set. Keeps (always-on):
per-step energy_decomp (11 components), density VTIs total/system/gs/wp/delta/
delta_coarse, observables.csv (energies+L2), wp_momentum_stats, wp_real_space_stats,
electron_number. Dropped extras are recoverable in post from saved densities.
Open option if extras are wanted later: regenerate GS with fewer extra_states
(→ headroom) — needs a new ~3 h GS.

## Milestone: 2026-07-19 ~08:20 — PIVOT to N=120 (memory confirmed)

Continuous GPU probe during a MINIMAL-set N=162 run: **min GPU free = 0 MB** —
the base 99-state ETRS propagation saturates the 24 GB A30. Intermittent step-0
crashes = borderline OOM (sometimes fits, usually not). Definitive.

**User directive: "let's try 120 electrons or so."** New system built:
- Config: `shared/configs/slab_n120_L60x60x62.hpp` — N=120, Lx=Ly=60.8 (152 pts),
  Lz=62.4 (156), dx=0.40, n0=1.2985e-3, r_s=5.686, extra_states=18 → 78 states.
  ALL z-geometry identical to n162 (slab z=±12.5, CAP ±[21.2,31.2] η=−1.0, WP
  launch −16.5, edge_width 1.0). Only the transverse box shrank.
- Memory: 78 states × 152²×156 ≈ 59% of n162 → ~14 GB used, ~10 GB free. Reliable.
- `gs/run.cpp` + `wp/run.cpp` swapped to SlabN120 (config include, Cfg, GS_DIR,
  cell strings). Orchestrator GSDIR → `shared_gs/slab_n120_L60x60x62_dx0p40`.
- Obs set stays MINIMAL in BASE (DIPCUR/MOM/OVL/WF/STATE all 0) for reliability;
  N=120 headroom MAY allow re-enabling STATE/WF later (pilot will show).

**Dir-name caveat:** campaign dir is still `scripts/mass_pair_n162/` but now runs
N=120 (kept for fast incremental rebuild — reuses compiled inq-study). Config +
run_summary record N=120 correctly. RENAME dir → `mass_pair` at final cleanup.

**Watch out:** CUDA teardown race — after `kill -9`, wait until `gpu_probe` shows
~23800 free before launching a new run on that GPU (else false OOM).

### Status (this milestone)
- IN FLIGHT: N=120 GS build+SCF on GPU0 (`gs/gs_n120_build_run.log`) → checkpoint
  `shared_gs/slab_n120_L60x60x62_dx0p40`. WP binary rebuild (CPU) `wp/wp_n120_build.log`.
- Stale n162 dbg_*/pilot dirs deleted. Old n162 GS checkpoint still on disk (unused).

## Milestone: 2026-07-19 ~08:34 — N=120 GS converging; orchestrator armed

- **N=120 GS SCF RUNNING** on GPU0 (`gs/gs_n120_build_run.log`), ~44 s/iter,
  energy converging (iter4 e=18.8 → toward ~−73). ~1–1.5 h to finish. Writes
  checkpoint `shared_gs/slab_n120_L60x60x62_dx0p40` + `gs/results/run_summary.txt`.
  (Deleted the STALE n162 `gs/results/run_summary.txt` — do not trust an old one.)
- **WP binary built** for N=120 (`wp/run`, 53.7 MB, 0 errors).
- **Orchestrator ARMED**: PID launched 08:33, log `orch_n120.log`, "waiting for GS
  checkpoint". Hardened this session:
  * `gs_ready()` now needs checkpoint dir non-empty AND `run_completed = true` in
    the summary (avoids acting on a stale/partial GS).
  * `wait_gpu_free(gpu, 20 GB)` before the pilot and every production launch
    (kills the CUDA teardown race — the repeated false-OOM trap).
  * BASE obs = minimal (DIPCUR/MOM/OVL/WF/STATE=0) for reliability.
- **Gotcha (cost me time):** `pgrep -f orchestrate.py` / `pkill -f mass_pair...`
  MATCH THE GREP COMMAND ITSELF → false "running" + self-kill (exit 144). Kill
  orchestrators by PID; check "running" by inspecting `/proc/<pid>/cmdline` and
  excluding your own shell.

### Next actions (resume here)
1. Wait for the orchestrator to detect the N=120 GS, run the pilot (minimal set,
   24 steps) — MUST clear step 0 (the N=162 failure point) with GPU free>few GB,
   and emit a REAL per-step cost projection. Watch `orch_n120.log`.
2. If pilot passes → m1(GPU0)+m2(GPU1) launch automatically (checkpoint every 500).
   If pilot fails → read `wp/rt_pilot.log`; if still OOM at N=120, drop extra_states
   in the config and regen GS, or reduce N further.
2. Watch the pilot clear step 0 (the old failure point) with GPU free > few GB.
3. If pilot shows generous headroom, consider re-enabling STATE/WF in BASE for a
   richer suite (optional).
4. Cleanup: rename dir → `mass_pair`; remove old n162 GS + config if abandoning it.
1. Watch GS converge → orchestrator builds WP. If WP compile fails, fix
   `wp/run.cpp`, delete any stale `wp/run`, relaunch orchestrator.
2. Confirm pilot passes (memory/norm/drift + s/step projection email).
3. Confirm m1+m2 launched; then hands-off (emails cover the rest).
4. Write `scripts/mass_pair_n162/analyse.py` (per-run inqview pipeline + density
   GIFs per the notebook-density-gif rule) during the ~28 h runs.
5. On completion: catalogue upsert, run notebooks, journal, commit (2-commit
   hygiene: infra/config vs run provenance; no forbidden words).

---

## Milestone 2026-07-19 13:30 — GS was KILLED, pipeline relaunched

**What stopped and why (ROOT CAUSE — proven from the prior-session transcript).**
The m1/m2 production runs never started; they were gated on the N=120 GS, which
was **killed mid-SCF at 08:36:20 local** (iter 9 of ~20+, e +410→−1.13 Ha toward
≈−73). NOT a crash, NOT OOM, NOT a stray `pkill` — **the HARNESS killed it on
session teardown.** Forensics:
- The GS build+run was launched at 07:24:25Z as a **`run_in_background=True` Bash
  TOOL task** (harness-tracked, task id `bocq44cc0`). The transcript records its
  final status verbatim as **`<status>killed</status>`**.
- The prior session's LAST action was launching the email watcher at 07:36:16Z
  (08:36:19 local). Immediately after, the session **ran out of context / ended**.
- On session end the harness tore down its tracked background tasks → SIGKILL to
  `bocq44cc0` at 08:36:20 local (1 s after the last foreground call returned).
  Hence the clean mid-eigenvalue truncation, instantly-freed GPUs, no error trace.
- The email watcher (PID 2538907) SURVIVED because it was started with
  `nohup … & disown` inside a *foreground* (completed) tool call → orphaned to
  init, NOT a harness task. That asymmetry is the tell.

With no GS checkpoint ever written, the orchestrator (PID 2536720) waited its full
4 h and timed out at 12:33:54 ("GS TIMEOUT" email), then exited.

**THE LESSON (supersedes the earlier `pkill` guess): never launch a long
autonomous run via the Bash tool's `run_in_background` flag — it is bound to the
session and dies when the session ends/compacts. Always detach with
`setsid`/`nohup … & disown` so the process is orphaned and session-independent.**
The relaunch below does exactly this.

**Corrective action taken (13:28–13:30).**
- Killed stale email watcher 2538907 by exact PID.
- Deleted empty `gs/results/run_summary.txt`.
- Relaunched GS **fully detached** (`setsid env CUDA_VISIBLE_DEVICES=0 ./run`,
  own session — immune to shell teardown / pattern-kill): **PID 2615069**,
  log `gs/gs_n120_relaunch.log`. SCF iter 0 e=409.80 (bit-identical → reproducible),
  ~48 s/iter, ETA ~20 min to −73 Ha. GS binary intact (no rebuild).
- Re-armed orchestrator **PID 2616425** (`orch_n120.log`) — waits GS → 24-step
  pilot gate → launch m1(GPU0)+m2(GPU1), ckpt every 500.
- Re-armed "both propagating" email watcher **PID 2616426** (`email_watch.log`).

**Lesson reinforced:** launch every long autonomous job with `setsid` (detached,
session-independent) — NOT the Bash-tool `run_in_background` flag (harness-tracked,
dies on session end; that is what killed the GS). Secondary: kill only by exact
PID, never `pkill -f mass_pair...`.

### Next actions (resume here)
1. Watch `gs/gs_n120_relaunch.log` converge (~20 min) → orchestrator auto-pilots.
2. Pilot MUST clear step 0 with GPU headroom (the N=162 failure point). If it
   OOMs even at N=120, drop `EXTRA_STATES` in `slab_n120_L60x60x62.hpp` + regen GS.
3. Then m1+m2 launch automatically; emails cover the rest.

---

## Milestone 2026-07-19 16:02 — both production runs propagating (m2 crashed once, retried)

**GS converged** (E=-53.8388 Ha, iter 96/300) → **pilot PASSED** (24/24 steps,
dE=+0.0002 Ha, wall 3528 s ≈ **147 s/step ⇒ ~102 h per 2500-step run**) → both
launched 15:49.

**m2 crashed at step 0** on GPU1 (`CUDA ERROR: illegal memory access`, SIGABRT
rc=-6); no checkpoint ⇒ orchestrator marked it permanent-fail after 1 try, emailed
"m2 RUN FAILED". **Root cause: memory still at the edge.** INQ rounds the grid UP
to **160×160×160 = 4.10M pts** (NOT the config's assumed 152×152×156 = 3.60M),
because 152=8·**19** has a bad FFT prime factor so INQ bumps it to 160=2^5·5. The
78-state ETRS peak at 4.10M pts sits right at the 24 GB A30 ceiling → pilot ran at
**0 MB free** (step 24 took 643 s of thrashing) → clearing step-0 alloc is a
coin-flip. Pilot + m1 won it; m2 lost. NOT mass-specific (m1/m2 identical bar mass).

**Fix applied:** relaunched m2 on the now-clean GPU1 (setsid-detached, exact
orchestrator env) — **cleared init on the 1st retry**. BOTH now propagating:
- **m1** pid 2678810, GPU0 — monitored + auto-resumed by orchestrate.py (pid 2616425).
- **m2** pid 2684101, GPU1 — monitored by **m2_watch.py** (pid 2687871, adopts via
  `wp/results/m2/m2.pid`, resumes EM_RESUME=1 on crash, ≤4 tries). Written because
  the orchestrator only tracks m1.

**Residual risk (honest):** both run at 0 MB free. Init cleared (main hurdle) and
the pilot ran 24 steps clean, so per-step crashes should be rare; checkpoint-every-500
+ both watchers cover mid-run deaths. The vulnerable window is steps 0–500 (no
checkpoint yet → a crash there = fresh restart, handled up to N tries). **Fallback
if crashes recur:** reduce EXTRA_STATES 18→~6 (78→~66 states) for real headroom,
regen GS (~1.3 h), relaunch both — NOT done now since both cleared init and are
progressing.

**Grids/params locked:** N=120, dx=0.40 (grid 160^3), 78 states, 2500 steps, dt=0.04,
CAP eta=-1.0 ±[21.2,31.2], sigma_WP=1.0, launch_z=-16.5, k0: m1=2.711 / m2=3.834.
