# Handover — Wide low-spread wavepacket (localised jellium)

Campaign: `docs/campaigns/localised_jellium/wide_wavepacket_lowspread.md`
(id `wide-wavepacket-lowspread`). Goal: fire a **wide, near-rigid WP** through the
localised jellium slab, matched to a same-width classical projectile, so any
WP−classical stopping difference is **purely quantum** (Pauli + interference),
not dispersion. Two-phase: Phase-0 human-in-the-loop design → Phase-1 autonomous
E×σ sweep.

---

## 2026-06-30 — design deliberation (Stage 2/4 via /campaigns); Phase-0 defined

**Status:** `draft`, 1/5. WP parameters + geometry + CAP + observables + Phase-0
all LOCKED in discussion. NOT yet launched. No new code/runs yet.

### Locked parameters
- **σ_WP = 3.5 Bohr, E = 300 eV** (v=4.70 a.u.). Analytical spread (free-Gaussian
  law): **2.6 % at the slab** (matched-σ essentially perfect), 17.6 % at far CAP.
- **Box 50×50×101 Bohr**, **dx=0.40** (refined from phase-5's 0.50 by the cutoff
  guard — see Phase-1). LZ=101 set so the equidistant launch gives
  **exactly 4σ** clearance (user: CAP↔WP gap = 4σ). Launch z₀=−26.5; clearance
  14.0 Bohr to BOTH slab face and CAP inner face. σ=3.5 sits on the 4σ rule.
- **Slab unchanged:** |z|<12.5, N=82, r_s≈5.67 (matched to production).
- **CAP = same as phase-5** (user decision): two-sided sin², **η=−0.7 Ha,
  10 Bohr/side**, region [±40.5,±50.5] in the LZ=101 box.
- **Observable suite = the EXTENSIVE phase-5 suite** (verbatim from
  `scripts/qsp_phase5/wp/run.cpp`): observables.csv (E components, current,
  dipole, density_l2), state_energies, occupations, eigenvalues, momentum_dist
  n(k,t), wp_momentum_stats, wp_real_space_stats (σ(t) — validates spreading),
  overlap/overlap_full, electron_number (N guard), density VTIs
  (total/system/gs_system/wp/delta/delta_coarse), wavefunction_wp. Classical run
  drops WP-specific, adds electron_track.csv (z,vz,F → S=−dKE/dz).

### Key findings (planning notebook)
- **Spreading solved at the slab for σ₀≥3** — it no longer discriminates σ₀; the
  driver moved to CAP load + cost.
- **CAP completeness — real INQ vs toy (§4).** Real INQ two-sided-CAP data
  (`systems/vacuum/hypotheses/twosided_cap_vs_mask/twosided_combined.csv`, has
  E=300 points) shows the phase-5 CAP reflects only **~0.2 %** at η=0.7 for the
  *narrow* benchmark packet (σ=4√2/k₀≈1.2). The 1D `cap_toy` forward model is
  **~5× pessimistic**. CAVEAT: real runs use the narrow packet, NOT our σ≈4 — so
  **wide-σ CAP adequacy is the one open risk**, closed by a Phase-0 wide-σ run
  (fallbacks: η→1.0, or Lhalf→15 / LZ≈111).
- Naming: dataset `L10`=L_total=10 (Lhalf=5/side); phase-5 "10/side"=Lhalf=10=
  dataset `L20`. Toy `L` = per-side = Lhalf.

### Phase-0 (defined; see campaign `<phase0>`)
- P0.1 operating point ✅ done.
- P0.2: matched classical UPF σ_pot=2.475 (verify by V(r)); new LZ=101 slab GS;
  **wide-σ CAP smoke** — gate: residual WP norm <2 %, N drift <2 %.
- P0.3: first matched WP+classical run; extract quantum S (energy method,
  E_ref=E_total(0)−⟨T_WP⟩−E_SIE) cross-checked vs n(k,t) centroid; classical S
  from electron_track. Done when user signs off S extracts cleanly.

### Artifacts (all in `docs/campaigns/localised_jellium/`)
- `wide_wavepacket_planning.ipynb` (executed) + `build_wide_wp_planning_notebook.py`
  (rebuildable: `venv/bin/python3 build_wide_wp_planning_notebook.py`).
- `wide_wp_planning_figs/`: `spreading_vs_sigma.png`,
  `cap_reflectivity_real_vs_toy.png`, `operating_point_table.csv`.

### Phase-1 sweep grid (LOCKED 2026-06-30; campaign `<phase1>`)
- One width σ_WP=3.5; **dx=0.40** set by the cutoff/aliasing guard
  (`k_max=π/dx ≥ k0+4σ_p`, σ_p=1/(√2·σ_WP)=0.202; phase-5's dx=0.50 fails at E≥500).
- **E ∈ {200,280,360,440,520,600} eV** (6 energies, 2-GPU friendly), cutoff margin
  6–20σ_p, spread@slab ≤3.9%. Each E = WP + matched classical → 12 runs + 1 vacuum SIE control;
  σ=0.5 phase-5 reused as reference; Lindhard analytical. Wall ~3–4 h/run
  (0.118 h/au), ~35 h compute / ~17 h on 2 GPUs.
- All high-v tail (v/v_F≈11–20; Bragg peak ~2–6 eV unreachable) → expect falling S(E).
- Verified in the planning notebook §7 (cutoff-guard assert passes, 0 errors).

## 2026-07-01 — Phase-0 execution (S2 done; P0b pair launched)

**Built (scripts/wide_wp/):** `gs/run.cpp`, `wp/run.cpp`, `classical/run.cpp` (all
env-driven, inq-study engine). UPF `shared/pseudopotentials/electron_gaussian_wpsigma3p5.upf`.
Config `shared/configs/slab_n82_L50x50x101.hpp` (dx=0.40, **edge_width=1.0**).

**S2 GS — DONE + validated** (`hypotheses/wide_wp/gs_validation.ipynb`):
checkpoint `shared_gs/slab_n82_L50x50x101_h0p40`. E_GS=**−86.04 Ha**, n0=1.312e-3 ✓,
r_s=5.667 ✓, 82 e⁻ ✓, SCF converged (iter 78).

**Findings/incidents this session (all resolved):**
1. **GS energy −86 Ha vs production −70.2** — finer dx=0.40 resolves the localised
   background's high-k content; NOT a bug (n0/r_s correct). Each S uses its own
   dx-consistent ref so it cancels; flag for the σ=0.5 overlay.
2. **edge_width 0→1.0** — sharp edge (|z|=12.5 not a grid node at dx=0.40) caused
   density sloshing / slow SCF; softening (per GS-study H1) fixed convergence
   (energy unchanged ~−86, so edge wasn't the −86; dx is).
3. **GPU contention** — a failed `pkill` (matches cmdline `./run`, not a path) left
   two GS runs thrashing GPU0 (89 s/iter). Kill by PID; pin CUDA_VISIBLE_DEVICES.
4. **classical include bug** — copied from `fullsuite_classical/` (2 levels deep)
   into `wide_wp/classical/` (3 levels) → includes off by one `../`. Fixed
   (config ../../→../../../, jellium ../../../→../../../../). Re-launched.

**P0b matched pair — RUNNING (2026-07-01 ~00:51):** WP (GPU0, `p0b_wp`, 750 steps
dt0.04) + classical (GPU1, `p0b_classical`, 1500 steps dt0.02) — both τ=30 a.u.,
E=300. **P0a folded in**: the WP run's residual-norm + N(t) is the CAP-completeness
check. ~3.5 h. **To check on completion:** WP residual inner norm <2%, N drift <2%,
E_total plateau, quantum S (energy method, E_ref=E_total(0)−⟨T_WP⟩−E_SIE) vs n(k,t)
centroid; classical S=−dKE/dz; spread@slab≈2.6% from wp_real_space_stats σ(t).

**AUTONOMOUS ORCHESTRATOR LAUNCHED (2026-07-01 02:38, PID nohup):**
`scripts/wide_wp/orchestrate.py` (log `orchestrate.log`). User switched to full
overnight autonomy: it waits for the P0b pair → numeric gate (both complete, WP
E finite, N not collapsed) → runs the 6-energy P1 sweep (WP GPU0 + classical GPU1
per energy, resumable, WRITE_EVERY=20 for speed, per-energy N_STEPS) → per-run S →
`hypotheses/wide_wp/se_wide_wp.csv` + `se_wide_wp.png` → email. Emails at the P0b
gate and at sweep completion. **ETA ~12–16 h (2 GPUs) → done tomorrow afternoon,
not by morning** (each run ~2–3 h at dx=0.40/LZ=101; honestly flagged to user).

**METHOD FINDING (light-projectile rule):** p0b classical track shows STRONG
deceleration (v/v0 → 0.86 already at slab entry, partly in the vacuum approach —
possible pre-contact / PBC self-image drag, cf. classical-projectile-fix). So:
- **Classical S = INITIAL DRAG** −dKE/ds over the early v≥0.85·v0 window (rule),
  NOT slab-average. Test on p0b partial: S(v0=300 eV) ≈ **5.6 eV/Bohr** provisional.
- **WP S = energy method** (phase-5) is a TRAVERSAL AVERAGE for a decelerating
  packet; rule-compliant WP S(v0) = n(k,t) centroid initial drag — data saved,
  extract on review. Flagged in the completion email.

### Next (on review, tomorrow)
1. Review the P0b gate email + the sweep S(E) email; check per-run gates.
2. WP centroid initial-drag extraction (rule-compliant) from saved n(k,t).
3. Assess the vacuum-approach drag (real pre-contact stopping vs PBC self-image) —
   may need the scientific panel / ties to classical-projectile-fix (A-vs-Z).
4. Bound SIE via a vacuum-WP control; then finalise S(E) + quantum component.
5. Autonomy-readiness checklist → flip `draft`→`ready` (or `done`).

## 2026-07-03 — overnight run FAILED; stall diagnosed; gate-review notebook built

**Status of the 2026-07-01 autonomy: PRODUCED NOTHING.** The orchestrator
(`scripts/wide_wp/orchestrate.py`) timed out and halted — `orchestrate.log`:
`[2026-07-01 10:38:51] P0b timeout — halting`. The P1 sweep never ran (no
per-energy dirs, no `hypotheses/wide_wp/se_wide_wp.{csv,png}`).

**Root-cause diagnosis (post-mortem; kernel/GPU logs root-restricted so not 100%
provable, but evidence converges):** both P0b runs died together at ~03:00 on
2026-07-01 (WP 03:00:06 @ step 415/750; classical 02:59:55 @ step 494/1500),
**healthy** to the last log line (no NaN/abort). Ruled out: reboot (uptime 157 d),
NVIDIA driver upgrade (none in apt/dpkg — NVML mismatch is the known cosmetic one),
the orchestrator (its log is 3 lines, silent 02:38→10:38, launched/killed nothing).
Both were launched by hand via **`inq-run` in the FOREGROUND** (no nohup/setsid).
**Most likely cause: SIGHUP when their controlling session closed** — two sibling
foreground jobs die together, 11 s apart = one signal. Two design flaws turned that
into a wasted night: the orchestrator (a) *waits* for a P0b pair it doesn't launch
and can't restart, and (b) `done()` only checks the completion sentinel — no
liveness check — so a dead run leaves it polling the full 8 h.

**Gate-review notebook BUILT (this session):**
`hypotheses/wide_wp/p0b_gate_review.ipynb` (builder `build_p0b_gate_review.py`,
figs `p0b_gate_review_figs/`) — ONE consolidated notebook with the 6 gate criteria,
each: plot + measured number + threshold + **blank "Verdict (you): ____"** (user owns
the verdict). Built on the PARTIAL data. Also built the two full-battery run-notebooks
`p0b_wp_run_notebook.ipynb` + `p0b_classical_run_notebook.ipynb`.

**Findings from the partial data (corrected — first pass had nan/window bugs):**
- C1 spreading: WP ballistic (v=4.70=k₀), σ_z tracks the FREE law within +2.8% at
  slab centre; abs spread +4% (entry face) / +13% (centre) on the *density* std
  (√2 faster than the σ_WP figure) — no anomalous spreading.
- C2 CAP: centroid reached z≈37 (CAP inner 40.5), absorption started — **INCONCLUSIVE**
  (killed before residual→0). Earlier "norm collapse" alarm was a MISREAD of the
  corrupted final row; N_total drift is only −0.21%.
- C3 bath: 83.00→82.82 (−0.21%) ✓. C4 energy: still evolving — INCONCLUSIVE.
- C5 quantum S: energy-method partial meaningless mid-run; n(k,t) centroid 4.70→4.59.
- **C6 classical: initial-drag S=5.6 eV/Bohr BUT net face-to-face ΔKE≈+3.3 eV →
  S_net≈+0.13 eV/Bohr (~40× smaller). KE largely RECOVERS across the slab (v dips to
  0.79·v₀ at centre, recovers on exit) = conservative well dominates.** Whether the
  physical classical S is ~0.13 or 5.6 is a USER interpretation call (ties to
  `stopping-power-extraction` + light-projectile rules) — NOT decided by me.

**User signed off the gate (2026-07-03)** and authorised restarting the failed runs
+ proceeding to the sweep — THEN asked to hold the trigger ("keep everything ready,
start only when I ask"). So everything is STAGED, NOT RUNNING.

**Hardened orchestrator READY (STAGED, not started):** `scripts/wide_wp/orchestrate.py`
was patched to fix the root cause:
- `launch()` now uses `start_new_session=True` (SIGHUP-immune — the 2026-07-01 death).
- New `run_p0b_to_completion()` / `_launch_p0b()`: the orchestrator now **launches the
  P0b pair itself** (detached, exact original config: WP 750/dt0.04/WE4, classical
  1500/dt0.02/WE10, k0=4.6957, launch_z=−26.5, CAP on) with a **liveness guard**
  (kill+retry once if a log is silent >25 min with no sentinel), replacing the broken
  passive 8 h wait. It renames any partial dir aside (non-destructive), then on
  completion rebuilds the gate notebook and flows into the existing per-run-completion-
  checked P1 sweep.
- Syntax-checked (`py_compile` OK).
- A standalone `relaunch_p0b.py` (P0b-only, holds at gate) also exists as an
  alternative if only the pilot is wanted.

**TO START (when the user asks):**
`cd scripts/wide_wp && setsid nohup venv/bin/python3 orchestrate.py > orchestrate.log 2>&1 &`
GPUs verified FREE (ml2218's 48 `lmp` are CPU-only; no nvidia device held).

**Incident 2026-07-03:** the orchestrator was briefly launched then immediately killed
(user said hold); partial p0b data was moved aside and restored intact (WP step 412 /
classical step 490 rows verified). No production data lost. Nothing is running now.

**Next (on the user's go):** start orchestrator → P0b completes+verified → gate notebook
auto-rebuilt on full data → user re-checks criteria 1–6 → sweep runs.

### Verified vs unverified
- Geometry, spreading numbers, gap-fit, sim-time extrapolation, real-vs-toy CAP
  reflectivity: **computed in the executed notebook** (0 errors; toy cross-checked
  bit-for-bit vs cap_toy).
- Wide-σ CAP completeness in a real run: **UNVERIFIED** — the Phase-0 gate.
- Sister campaign `classical-projectile-fix` brainstorm is PAUSED (decisions
  captured in its file's resolved_decisions).

## 2026-07-03 (session 2) — grilling → open-z ENLARGED-BOX long run STAGED

**Trigger:** user's concern that the WP "was not entirely captured by the CAP
before it started looping around". Grilled the fix (`grill-with-docs`) to a shared
design; then built + staged a single long exploratory run. **Not the sweep** —
this run is to WATCH the energy plateau, then the user reviews + launches the sweep.

### Decisions locked (8-question grill)
- **Box z: 101 → 111 Bohr** (+10 CAP runway). z ∈ [−55.5, +55.5].
- **CAP: η −0.7 → −1.0 Ha, 10 → 14 Bohr/side**, region [±41.5, ±55.5], inner ±41.5.
- **Boundary: periodicity 2 (open-z)** — USER CHOICE, overriding the GS study's
  "use PBC" verdict; open-z's net-charge G=0 monopole is knowingly accepted, to
  **debug later** (Q(t) logged in `electron_number.csv` for post-hoc subtraction).
- **Run: WP only, dt=0.04**, no classical pair.
- **Length: τ = 3 × rigid end-to-end traversal at mean k0** = 3·111/4.696 ≈ **71 a.u.
  → N_STEPS = 1775**. WRITE_EVERY=6 (~296 VTI frames, cadence rule), WF_EVERY=40.
- **S method:** energy long-time PLATEAU (primary) + n(k,t) centroid cross-check.

### Physics established (verified against inq source)
- `periodicity(2)` changes ONLY the Poisson kernel (`poisson.hpp:190`, 2D-truncated
  Coulomb); the kinetic propagator is ALWAYS FFT (`laplacian.hpp:34`) → open-z does
  **NOT** stop WP grid-wraparound (box+CAP do). It's an electrostatics choice.
- **Monopole caveat (the one real risk of open-z here):** net-charged cell (jellium
  is an external well, 82 e⁻ + WP, no compensating charge). Under periodicity 2 the
  Hartree G=0 term = `0.5·rc²` (`poisson.hpp:49`), so E_total carries `∝ L_z²·Q²`.
  Constant while Q fixed (cancels during traversal), but STEPS when the CAP drains
  the WP (Q:83→82) — exactly when the long-time plateau forms. So the plateau LEVEL
  is monopole-shifted (GS study measured the analogous ΔQ=1 shift: +4.4→−2.1 eV);
  judge plateau SHAPE now, correct the level later. Under PBC this term is 0 (why the
  GS study said use PBC) — user chose open-z anyway.

### Built this session (compile-gated)
- NEW config `shared/configs/slab_n82_L50x50x111.hpp` (LZ=111, open-z, monopole
  caveat documented). Struct `SlabN82_L50x50x111`.
- `scripts/wide_wp/gs/run.cpp` + `wp/run.cpp`: → new header, `.periodicity(2)`,
  CAP η−1.0/48.5·111/14·111, GS_DIR → `shared_gs/slab_n82_L50x50x111_h0p40_per2`
  (env-overridable). Binaries rebuilt via `cmake --build build --target run`.
- NEW hardened orchestrator `scripts/wide_wp/run_long_wp_per2.py` (py_compile OK):
  GS(if absent)→WP→run-notebook→email. SIGHUP-immune (`start_new_session`),
  liveness guard (kill+retry once on 25-min silent death), MAX_HOURS=14, failure
  alerts. On WP completion: builds the **run-notebook** (skill-local builder,
  full single-run battery incl. the LOCKED Fourier `fft_pipeline_panel`) into
  `hypotheses/wide_wp/wp_per2_E300_long_run_notebook.ipynb`, generates an
  E_total(t) plateau + N(t) figure, and sends a 4-part result email.

### TO START (detached, survives logout)
```
cd .../scripts/wide_wp
setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 run_long_wp_per2.py \
    > run_long_wp_per2.log 2>&1 &
```
GPU: both free (2×25 GB; ml2218's `lmp` are CPU-only), CUDA probe OK. GPU 0 used.

### LAUNCHED 2026-07-03 23:29 — GS running (GPU 0)
- Both binaries compiled clean (compile-gate passed). Orchestrator launched detached
  (`setsid nohup ... run_long_wp_per2.py`, pid 3116218), log `run_long_wp_per2.log`.
- **GS running**: log confirms `Periodicity = 2d (slab)`, cell 50×50×111, N=82,
  r_s=5.667; GPU0 fully allocated, ~33 s/SCF-iter. Initial energy e≈+422 Ha (vs
  ~−86 under PBC) = the open-z G=0 monopole offset active (expected, not a bug).
- **Sentinel bug found + fixed on first launch attempt:** the stale July-1
  `gs/results/run_summary.txt` (cell 50×50×101) fooled `run_to_completion` into
  "GS already complete" → WP launched against a missing per2 checkpoint and FATAL'd
  cleanly (no GPU wasted). Fix: orchestrator now `unlink(missing_ok)`s the stale GS
  summary before launching GS, and the GS-present gate uses `GS_CKPT.iterdir()`.
  Stale summary + partial WP dir cleaned; relaunched OK.
- ETA: GS ~40–45 min → WP ~6 h → run-notebook + 4-part email. Orchestrator emails
  on completion AND on failure; liveness guard (kill+retry once on 25-min silence).

### 2026-07-04 — GS done; WP running; MAX_HOURS timing bug fixed via waiter
- **GS DONE 00:05** (36 min). E_GS = **+65.0 Ha** (open-z monopole-inflated; PBC would
  be ~−86). Checkpoint `shared_gs/slab_n82_L50x50x111_h0p40_per2`, run_completed=true.
- **WP running** from 00:05 (pid 3139197, GPU0). Energy stable (~76 Ha, tiny de/step).
- **Timing bug:** realised **~38 s/step** (open-z 2D-Poisson `enlarge({1,1,2})` doubles
  the z-grid + WRITE_EVERY=6 VTI writes), NOT the ~13 estimated → full 1773 steps
  ≈ **18.7 h (WP done ~18:45 on 07-04)**. The orchestrator's `MAX_HOURS=14` would have
  KILLED the WP at ~13:29 (step ~1269, incomplete → false failure email, no notebook).
- **Fix (no data lost):** retired the orchestrator (killed pid 3116218; the detached WP
  survives), refactored its post-WP logic into `run_long_wp_per2.finish()`, and launched
  a standalone **`finish_wp.py` waiter** (pid 3195981, 24 h deadline) that WATCHES the
  existing WP (never relaunches), with a liveness/death guard, and calls `finish()` on
  the real completion sentinel → run-notebook + 4-part email. Revised WP ETA ~18:45.

### 2026-07-04/05 — WP COMPLETE; KEY RESULT: CAP works, but E_total does NOT plateau
- **WP finished normally 2026-07-04 17:24** (wall 62,395 s ≈ 17.3 h, all 1773 steps,
  `run_completed = true`). Data at **`wp/results/results/wp_per2_E300_long/`**.
- **OUTPUT-PATH BUG (fixed for future):** `wp/run.cpp` prepends `"results/"` to `LJ_OUT`,
  but the orchestrator passed `LJ_OUT="results/wp_per2_E300_long"` → data landed in
  `results/results/...`. The waiter watched `results/...`, so at 17:27 it fired a FALSE
  "PROCESS DIED" email and aborted. **The run was fine.** Fixed `WP_OUT` → bare name +
  `wp_res = WPDIR/"results"/WP_OUT` in `run_long_wp_per2.py` + `finish_wp.py`.
- **RESULT (recovered manually; honest email sent 2026-07-05 03:04):**
  1. **CAP fully absorbs the WP — cleanly.** WP orbital `norm_check` decays 1 → **5e-8**
     AND **N_total 83.000 → 82.000** in lockstep (verified full-resolution
     `electron_number.csv`: first 83.00000, last 81.99997, drop = 1.00003). NO wraparound;
     the box+CAP fix (η−1.0, 14 Bohr/side, LZ=111) WORKS. The original "not captured
     before looping" concern is RESOLVED. **NB (2026-07-06):** a mid-session claim that
     "N stays ≈83" was a self-inflicted indexing bug (indexed the 1774-row N array with
     296-row observables-grid indices → only read t<12 a.u., pre-absorption). There is
     **NO norm-vs-count paradox** — the panel's E4 caught this by reading the files.
  2. **E_total does NOT plateau** — after absorption it OSCILLATES ±20–35 eV with
     period ~24 a.u. ≈ **2·ω_p** (plasmon ringing: ω_p=0.128 Ha, T_p=48.9 a.u.), through
     t=70.9. dE_total(end) = −16.5 eV is a point on the oscillation, not a plateau.
  → **The plateau-based energy S method does not cleanly apply as-is.** Open question
     for the user: time-average over the ringing / longer run to test damping / revisit
     method. (This is exactly the look-first result wanted before the sweep.)
- Fig `hypotheses/wide_wp/wp_per2_E300_long_plateau.png` (E_total(t)+N(t)) sent.
  Run-notebook `hypotheses/wide_wp/wp_per2_E300_long_run_notebook.ipynb` BUILDING
  (bg pid 3601976, log `scripts/wide_wp/notebook_build.log`) — its collective-response
  FFT should confirm/deny the 2·ω_p ringing attribution.
- Whether the ±30 eV E_total swing is genuine plasmon ringing vs an open-z / energy-
  conservation artifact is UNVERIFIED — needs the notebook + user interpretation (a
  closed unitary system after the projectile leaves should conserve E_total exactly;
  the CAP is far from the bath, so the swing is notable). Candidate for `scientific-panel`.

### 2026-07-06 — user switched to FULL PBC; ringing PERSISTS (not an open-z artifact)
- **User edited gs+wp run.cpp → full 3D PBC** (`.periodic()`, checkpoint `_pbc`) and
  launched the run manually. PBC GS done (**E_GS=−99.3 Ha**, physical — vs +65 open-z).
- **PBC WP run (user's, pid 4031791, GPU0)** healthy: `Periodicity=3d`, WP norm=1,
  E_total sane ~−88 Ha. At handover: step ~1073/1773, N:83→82.0002 (CAP absorbing).
  Output `wp/results/wp_pbc_E300` (bare LJ_OUT → no double-results). ETA ~20:00 07-06.
- **KEY FINDING:** E_total STILL oscillates ±~25 eV under PBC (−87.8…−89.7 Ha), same
  as open-z → the swing is **NOT the G=0 monopole**. Leading hypothesis: the CAP
  drains energy as density sloshes through the absorber (a CAP makes E_total a
  non-conserved sink), and/or genuine bath dynamics. So **no clean plateau under
  either BC** so far → the energy-plateau S method needs a rethink (time-average,
  CAP-free diagnostic, or a different estimator). Strong `scientific-panel` candidate.
- **Set up (non-interfering):** `finish_pbc.py` WATCH-ONLY finisher (pid 135962) —
  watches `wp_pbc_E300`, builds run-notebook + emails on completion, NO relaunch,
  liveness/death guard. Did NOT launch anything on the GPUs (user owns the run).
- Open-z run-notebook built OK earlier: `hypotheses/wide_wp/wp_per2_E300_long_run_notebook.ipynb`
  (43 cells; density-GIF battery skipped on an index error — minor, revisit).
- Output-path bug fixed in `run_long_wp_per2.py`/`finish_wp.py` (bare LJ_OUT); the
  user's manual PBC launch already used the bare name so its data path is clean.

### 2026-07-07 00:04 — PBC run COMPLETE; no-plateau CONFIRMED (definitive)
- **`wp_pbc_E300` DONE** (1773 steps, wall 78,316 s ≈ 21.8 h, run_completed=true).
  `finish_pbc.py` built the notebook (`hypotheses/wide_wp/wp_pbc_E300_run_notebook.ipynb`,
  45 cells) + plateau fig + **emailed the user 00:04 07-07**. All processes idle.
- **VERDICT (full trajectory):** E_total oscillates **−87.6…−89.8 Ha (~59 eV p-p)**
  through the whole run — **NO plateau**, identical in magnitude to open-z. N:83→82.00
  (CAP fully absorbs WP). So the ring is **not** the monopole and **not** a BC artifact
  → it is the **always-on CAP acting as a non-Hermitian energy sink** (E_total not
  conserved as density crosses the absorber) and/or genuine bath dynamics.
- **CONCLUSION for the method:** the plateau-based energy S estimator cannot work with
  an always-on CAP under EITHER BC. Robust alternatives: (a) **n(k,t) momentum-centroid**
  S (immune to the CAP energy bookkeeping — recommended); (b) a **CAP-free / CAP-η-varied**
  diagnostic to confirm the sink hypothesis; (c) time-average over the ring. Strong
  `scientific-panel` candidate. This is the decision point before any S(E) sweep.
- Minor: density-GIF battery still skips on an index error in the run-notebook builder
  (both open-z + PBC notebooks) — revisit `run_notebook_builder.py` GIF battery.

### 2026-07-06 — run-notebook delivered; ΔE_total plot added; sweep HELD; panel launched
- **Run-notebook DONE + valid:** `hypotheses/wide_wp/wp_per2_E300_long_run_notebook.ipynb`
  (43 cells now 44, includes the locked Fourier/momentum panel, 0 errors). Markdown
  report referencing pre-rendered PNGs in `..._run_notebook_figs/`.
- **ΔE_total(t) plot added to the notebook (user request):** new cell at the top of
  the Energetics section + fig `..._run_notebook_figs/delta_total_energy_vs_time.png`
  (canonical theme). ΔE_total = (E_total(t)−E_total(0))·27.2114 eV. **Does NOT plateau:**
  late-time (t>30) mean −4.5 eV, **std 21 eV, ptp 62 eV**, end −16.5 eV; dominant
  ω≈0.31 Ha (T≈20 a.u.) ≈ 2.4·ω_p. NOT yet folded into the run-notebook BUILDER
  (`.claude/skills/run-notebook/run_notebook_builder.py`) — TODO if it should persist
  across rebuilds.
- **DECISION (AskUserQuestion 2026-07-06):** sweep geometry = **new open-z (111, η−1.0)**;
  S-extraction = **HOLD — resolve the no-plateau method question first** before committing
  the ~200 GPU-hr sweep (open-z is ~38 s/step ⇒ ~17 h/run × 12 runs). Sweep NOT launched.
- **Scientific panel LAUNCHED** (background Workflow, run `wf_1ed49c5c-f43`, task
  `w8ulfpgbm`): 4 experts (TDDFT methodologist / condensed-matter / jellium-stopping /
  data-custodian) + judge, opus-high. Question: is the ±21 eV ΔE_total ring genuine
  plasmon/collective physics or a method artifact (CAP energy/charge bookkeeping, open-z
  monopole, renormalisation, grid/dt), and can S be extracted anyway? Brief carries the
  real numbers + the N-vs-norm inconsistency as EVIDENCE (my hypothesis kept out per the
  anchoring guard). Verdict pending → relay to user, user owns verdict.

### 2026-07-06 — PANEL VERDICT + verification: the non-plateau is a POST-PROCESSING ARTIFACT
- **Panel done** (run `wf_1ed49c5c-f43`, 9 agents, opus-high). Verdict: the ±21 eV
  ΔE_total ring is a **method artifact (Fork A: phantom absorbed-WP orbital)**, NOT bath
  physics and NOT open-z electrostatics. S(300 eV) ≈ 0 (≪1 eV), consistent with the
  Lindhard high-v tail (v/v_F≈14). ΔE_total/L_z is unusable here.
- **VERIFIED by direct re-ledger of the CSVs (decisive test):**
  - `corr(E_kin_total, e_kin_ha_WP) = +1.000`; implied bath kinetic
    `E_kin_total − e_kin_ha_WP` is **flat, std 0.016 Ha**.
  - **Drop state-60 kinetic → ΔE_total plateaus: late-time std 0.00 eV, deposit +1.7 eV.**
  - Mechanism: CAP removes the electron from the DENSITY (N 83→82) but energy_total still
    counts KS state 60 at occupation 1 with its NORMALIZED ⟨T⟩≈11 Ha (rings 9.8↔11.8 Ha,
    ω≈0.31 Ha) as the CAP reshapes the norm-1e-8 residual. Pure bookkeeping.
  - The `wp_momentum_stats.csv` `norm_check` column is NOT the norm (6.9e7→3.2, blows up
    the naive norm-weight); use `wp_real_space_stats.csv` `norm_check` (→5e-8) or just drop
    state 60 for the re-ledger.
- **Notebook cell corrected** (`3433a2bc`) + fig `delta_total_energy_vs_time.png` now shows
  as-logged (rings) vs phantom-removed (plateaus). Earlier cell text (wrong N + "does not
  plateau") replaced.
- **Fixable in post-processing** — no rerun needed to salvage the energy method: norm-weight
  or drop the absorbed orbital's kinetic from the ledger. The panel recommends the campaign
  adopt the **pre-absorption momentum-centroid −dKE_WP/ds** channel as the primary light-
  projectile S estimator and formally drop ΔE_total/L_z.

### Next (pending user decision)
- **User to ratify the verdict** (they own it). Open questions the panel raised:
  1. Confirm in the inqkit energy-assembly whether state 60 enters energy_total at occ 1
     (crux of the fix — where to norm-weight/drop it). → this is the real code fix to make.
  2. Resolve dipole_z amplitude conflict (±11 vs ~0.03 a.u.) — sets the true bath-mode size.
  3. **Campaign-level:** S(300 eV)≈0 sits at the numerical floor. Lower v toward the Bragg
     peak (v/v_F≈1–2) where S is large & cleanly extractable? And switch production geometry
     to charge-compensated PBC to retire the open-z ambiguity?
  4. Adopt pre-absorption −dKE_WP/ds as primary S estimator, drop ΔE_total/L_z?
- **STILL PENDING (unrelated to the artifact):** why the concentrated-WP runs equilibrated
  but this one appeared not to — likely the same phantom term was smaller/shorter there;
  cheap check = re-ledger a past concentrated run's e_kin_ha. NOT yet done.
- Sweep remains HELD until the S-method (initial-drag vs re-ledgered plateau) + geometry
  (open-z vs PBC) are chosen.
- OPEN: open-z net-charge G=0 monopole correction (deferred, user "debug later").

### 2026-07-06 — PBC WP-vs-classical matched pair LAUNCHED (user directive)
- **User directive:** rerun the `wp_per2` config under FULL 3D PBC (periodic in z too),
  as a matched pair — a WP AND a classical projectile at the same Gaussian width — then
  build 2 run notebooks + 1 comparison ("phase") notebook. (User: "analogous to quantum
  kick"; "do not use workflows".) Plan: `docs/plans/wide-wp-pbc-classical-comparison.md`.
- **Config (the ONLY change vs wp_per2 = periodicity 2→3):** box 50×50×111, dx=0.40,
  **periodicity 3 (full PBC)**, CAP η−1.0/14-Bohr matched across both, E=300 eV, σ=3.5.
  WP dt0.04 N1773; classical dt0.02 N3540 (matched τ≈70.8), UPF
  `electron_gaussian_wpsigma3p5.upf` (σ_pot=2.475). Both label σ=3.5.
- **Edits (rebuilt, compile-gate PASSED all 3):** `gs/run.cpp` + `wp/run.cpp`
  `.periodicity(2)`→`.periodic()` + GS dir `_per2`→`_pbc`; `classical/run.cpp` retargeted
  101→111 config + CAP η−0.7/10 → η−1.0/14 + LJ_GS_DIR env→pbc. run_summary boundary
  strings updated.
- **GS (PBC-111) running** GPU0 (pid 3994910, log `gs/gs_pbc.log`); confirmed
  `Periodicity = 3d (fully periodic)`, cell 50×50×111, r_s=5.667. Expect E_GS≈−86 Ha
  (clean, no monopole). Stale `gs/results/run_summary.txt` (old 111 GS) removed to avoid
  the sentinel trap. Checkpoint → `shared_gs/slab_n82_L50x50x111_h0p40_pbc`.
- **GPU check (NVML down, used CUDA driver API cuMemGetInfo):** dev0 0% free (GS), dev1
  99% free. Other user's python is NOT on a GPU. Both GPUs available for the pair.
- **Orchestrator LAUNCHED** detached: `scripts/wide_wp/run_pbc_pair.py` (pid 4004450, log
  `run_pbc_pair.log`). Waits for GS sentinel+checkpoint → launches WP (GPU0) + classical
  (GPU1) concurrently, per-job liveness guard (kill+retry once, STALL 30min, MAX_HOURS=24)
  → builds 2 run notebooks (run-notebook skill) + comparison notebook
  (`hypotheses/wide_wp/build_wp_vs_classical_pbc.py`, smoke-tested on open-z WP data:
  loaders + initial-drag S + de-ledgered ΔE all OK) → 4-part email. Output paths: bare
  LJ_OUT (`wp_pbc_E300`, `classical_pbc_E300`) → data at `<dir>/results/<OUT>`.
- **ETA:** GS ~30-40min → pair concurrent ~15h (classical dt0.02×3540 is the long pole) →
  done ~tomorrow afternoon. Emails at completion AND failure.
- **S extraction (both):** INITIAL DRAG (light-projectile rule) — classical −dKE_ion/ds over
  early v≥0.85·v0 (electron_track.csv); WP momentum-centroid −d(½pz²)/ds. WP ΔE_total
  re-ledgered (drop absorbed-orbital kinetic) as cross-check. This PBC run also tests that
  the phantom-orbital artifact persists under PBC (predicted: yes — not a boundary effect).

#### 2026-07-06 03:10 — classical run length CORRECTED; orchestrator→waiter
- **Timing bug caught at launch:** WP ~30 s/step (no moving ion) but **classical ~58 s/step**
  (moving Gaussian-ion pseudopotential re-projected each step). At the originally-set 3540
  steps that is **~59 h**, over the 24 h guard → would be killed mid-run. Classical physics
  is done by t≈15 a.u. anyway (projectile absorbed), so τ=70.8 was overkill.
- **Fix:** stopped `run_pbc_pair.py` orchestrator + the 3540-step classical; **left WP running**
  (pid 4031791, untouched); **relaunched classical at N_STEPS=1500 (τ≈30 a.u., dt0.02)** on
  GPU1 (compute pid 4057426 under MPI orted). τ=30 covers launch→slab→full CAP absorption
  + margin (light-projectile rule sizing). `run_pbc_pair.py` CL_STEPS updated 3540→1500.
- **Waiter replaces orchestrator:** `finish_pbc_pair.py` (pid 4063744) WATCHES both running
  runs (never relaunches) → `run_pbc_pair.finish()` on both sentinels. Liveness = **log
  mtime staleness > 30 min** (robust vs the MPI/pid-naming quirk that fooled the first waiter
  attempt — the compute proc is `run` under `orted`, not the full path).
- **IGNORE** the false "[wide-wp PBC] CLASSICAL DIED" email (03:05): first waiter used a
  pgrep pattern + wrong pid (captured the exited launcher, not the MPI compute child). Both
  runs verified stepping after: WP step 71, classical step 1 (e=−94.5 Ha, finite).
- **Revised ETA:** WP ~15 h (done ~18:00), classical ~24 h (done ~03:00 tomorrow) → pair
  completes ~tomorrow morning; waiter builds notebooks + comparison, emails.
- **GS energy note:** PBC-111 E_GS = **−99.33 Ha** (vs box-101 PBC −86; box-dependent Ewald/
  neutralising-background term for the net-charged cell — cancels in each run's own ΔE ref).

#### 2026-07-06 07:12 — CLASSICAL PBC self-image drag DIAGNOSED + FIXED (→ open-z)
- **Problem:** the classical PBC projectile lost **5.6 eV/Bohr in pure VACUUM** (z −26.5→−12.5,
  before any electrons; KE 300→225 eV; constant force). Spurious.
- **Discriminating test (decisive, same PBC/box/E/σ):** the WP-PBC projectile centroid shows
  vacuum drag **−0.05 eV/Bohr** (pz 4.696→4.70 flat, S≈0) — CLEAN. So the artefact is specific
  to the **classical point charge**, not a general PBC field. = the documented PBC z self-image /
  charged-point-in-PBC (Makov-Payne-class) artefact; the delocalised WP is immune. (The track
  NaN at slab entry was a benign ofstream-buffer flush artifact — electronic e evolved smoothly.)
- **User decision (AskUserQuestion):** fix = **classical → open-z (periodicity 2)** (kills the
  z self-image; documented fix). WP stays PBC (BC-insensitive: S≈0 in both open-z and PBC), so
  the comparison stays valid: classical(open-z) vs WP(PBC).
- **Deployed:** stopped the artefact classical-PBC + waiter (WP-PBC untouched, still running);
  `classical/run.cpp` `.periodic()`→`.periodicity(2)`, GS default→per2, boundary string; rebuilt
  (compile-gate OK); relaunched **classical open-z** GPU1: LJ_OUT=`classical_openz_E300`,
  GS=`..._per2`, N_STEPS=**1250** (τ≈25 a.u., ~23 h at ~58 s/step), dt0.02. Loads the open-z per2
  GS. Step 1 e=+69.85 Ha (open-z monopole-inflated; classical S uses the ion KE track, NOT the
  electronic energy, so the offset is irrelevant). No NaN.
- `run_pbc_pair.py` updated (CL_OUT/CL_STEPS/GS_PER2); `finish_pbc_pair.py` CL_LOG→classical_openz.log;
  waiter relaunched (pid parked, watching `classical_openz_E300` + `wp_pbc_E300`).
- **TO CONFIRM as the run progresses:** classical open-z vacuum drag should collapse 5.6→~0
  eV/Bohr (check ~step 30, z≈−23.6: vz should stay ≈4.696, vs PBC's 4.568). The comparison
  notebook will show it. Comparison is now WP(PBC) vs classical(open-z) — matched geometry,
  BC differs only in the projectile-electrostatics knob (WP verified BC-insensitive).

#### 2026-07-07 06:33 — PAIR COMPLETE; open-z fix FAILED — classical S is still the artefact
- **Both runs done.** WP `wp_pbc_E300` 1773/1773 (E_total −89 Ha, PBC/physical); classical
  `classical_openz_E300` 1250/1250 (E_total +70.8 Ha, open-z monopole). Both waiters exited
  cleanly. Notebooks built: `wp_pbc_E300_run_notebook.ipynb` (45 cells),
  `classical_pbc_E300_run_notebook.ipynb` (41 cells), `wp_vs_classical_pbc_comparison.ipynb`
  + `wp_vs_classical_pbc_S.png`. Two emails sent (WP-only 00:04, pair 06:33). Known-harmless
  warnings: density-GIF battery skipped ("list index out of range"); ParaView phases skipped.
- **Comparison reported S_classical=5.56, S_WP=0.046 eV/Bohr (~120×).** The WP number is real
  (diffuse σ=3.5 WP barely couples; vacuum drag ≈0). **The classical number is NOT.**
- **DECISIVE CHECK (contradicts the "TO CONFIRM" above):** over **13.5 Bohr of pure vacuum**
  (electron_track.csv steps 0..154, z −26.5→−13.0, i.e. BEFORE the slab edge at −12.5, no
  electrons present) the classical projectile decelerates at a **constant** force:
  vz 4.696→4.066, KE 11.02→8.27 Ha ⇒ **vacuum drag = 5.564 eV/Bohr**, constant to 0.07%
  (dvz/step mean −4.13e-3, stdev 2.8e-6). S_classical=5.564 = this artefact **exactly**.
- **=> The open-z switch did NOT collapse the drag (5.6 → still 5.6).** So the artefact is
  **not** the z self-image (open-z frees z). It is a constant vacuum self-force specific to the
  sharp classical charge (diffuse WP immune in BOTH BCs) on an x,y-periodic cell ⇒ leading
  suspect = **in-plane (x,y) periodic self-interaction / compensating-background gradient on the
  charged Gaussian**, untouched by the z boundary. The "120× gap" is this artefact, not physics.
- **STATUS:** WP result trustworthy; **classical S unusable until the vacuum self-force is found
  and removed.** Strong `scientific-panel` / `diagnose` candidate. Do NOT report 5.56 as stopping.
- **Next-test ideas (unstarted, need user OK):** (a) classical projectile in a FULLY empty cell
  (no jellium background) at several box widths L_xy → is the drag ∝ 1/L_xy³ (dipole-image) or
  independent (background gradient)? (b) widen L_xy to test in-plane image convergence.
  (c) subtract a measured vacuum-drag baseline from the in-slab track as an interim S estimate.

#### 2026-07-09 14:05 — wide-WP slides built for the Emilio deck (two case studies)
- **Deck:** `docs/reports/09-07-2026-meetng-emilio/emilio_deck_draft1.pptx` (25 slides).
  User asked (2026-07-09) to build the Section-3 wide-WP slides and add the right
  plots. User decision: make **TWO case studies** treated identically —
  `wp_pbc_E300` (full 3D PBC) and `wp_per2_E300_long` (open-z, periodicity 2);
  same σ_WP=3.5, E=300 eV, 1773 steps (~1.71 fs), differ only in boundary.
- **New builder:** `assets/make_s3_wide_wp.py` writes into `figures/`:
  - triptych GIFs `s3_2_wwp_{pbc,openz}_total_{total,dfirst,dprev}.gif` — via
    `hypotheses/_density_views.render_decomposition_views` TOTAL-ONLY path
    (`wp_dir=None`, 296 `density_total` frames; density_wp only 15 so decomposition
    avoided). Rules-compliant: fixed shared clim across frames, linear+log panels.
  - energy PNGs `s3_energy_{pbc,openz}.png` — ΔE_total(t) physical (absorbed
    WP-orbital kinetic `e_kin_ha` removed, per reference_phantom_absorbed_wp_orbital_energy)
    vs as-logged; + N(t) CAP drain (83→82). Time axis in fs.
- **Result (the campaign's plateau question — answered YES):** physical de-ledgered
  ΔE plateaus DEAD-FLAT (late-20%-window std = 0.00 eV) at **+0.6 eV (pbc)** /
  **+1.7 eV (openz)** — near-zero deposited energy, consistent with the wide WP's
  negligible stopping (S_WP≈0.05). As-logged trace rings ±40 eV (phantom orbital KE).
- **Deck wiring** (`build_emilio_deck.py`): added `placeholder()` helper; S3.1
  concentrated-WP left as a BLANK titled placeholder (user fills — wrong example
  previously). Slides 14 (blank), 15/16 (pbc density+energy), 17/18 (openz
  density+energy). Learnings left blank per house rule; substance in captions.
- **Verified:** deck assembles (25 slides, correct pic counts 3/1/3/1); GIFs
  regenerated cleanly (exit 0, ~10 min). NOT verified by eye (figure-work rule: user
  previews). Rebuild: `make_s3_wide_wp.py` then `build_emilio_deck.py`.

#### 2026-07-09 14:45 — energy figs corrected + delta-density GIFs added to run folders
- **Energy slides remade (user: "something is off").** Removed the blue de-ledgered
  "absorbed-orbital-KE-removed" trace — it manufactured a flat plateau the RAW data
  does not have. `s3_energy_{pbc,openz}.png` now show a SINGLE original-data trace
  ΔE_total(t) (from observables.csv only), shaded band relabelled **"in slab"**
  (rigid packet centre between slab faces), + N(t) from electron_number.csv.
  - Raw ΔE_total does NOT plateau: late-20%-window std ≈ 17 eV, peak-to-peak ≈ 52 eV,
    ends −18.7 eV (pbc) / −16.5 eV (openz). This is the known ±40 eV energy ringing
    / no-plateau (CAP-as-energy-sink hypothesis). Deck slide titles changed from
    "…reaches a plateau" → "…total energy and electron number"; captions state it
    oscillates and does not settle. Builder: `assets/make_s3_wide_wp.py::energy_fig`.
- **Delta-density GIFs added to BOTH run notebook-figs folders** (user request, via
  report-figures skill). New generator `hypotheses/wide_wp/make_delta_density_gifs.py`
  → into `wp_pbc_E300_run_notebook_figs/` and `wp_per2_E300_long_run_notebook_figs/`:
  - `density_delta.gif` (Δn=n(t)−n(0), clim ±2.1e-3)
  - `density_delta_instantaneous.gif` (Δn=n(t)−n(t−Δt), clim ±7.3e-5)
  Rules-compliant (linear+log two-panel, fixed shared clim, load_vti physical order).
  Clims match across the two runs → case studies stay comparable. 296 frames each.
