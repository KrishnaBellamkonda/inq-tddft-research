# Handover — Nazarov–Gross mass ladder (dense slab, wide wavepacket)

**Rolling file. Latest milestone at top.**
**Repo:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research` (branch `quantum-stopping-power`)
**Machine:** CSD3, GPU partition `ampere` (A100, sm_80), account `mphil-nikiforakis-skcb2-sl2-gpu`
**Plan:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/nazarov-gross-slab-mass-ladder.md`
**Source note:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/sources/nazarov-2025-quantum-projectile-stopping.md`

---

## 2026-08-05 (2nd) — RELAUNCHED after fixing a walltime bug that would have crippled scheduling

**Symptom.** The first launch's three compile jobs were given an estimated start
of **10:29, a 7.5-hour wait**, on an `ampere` partition holding 563 pending jobs.

**Cause — a real bug, not just a bad default.** `submit_rt` passed the per-run
walltime through the `SBATCH_TIMELIMIT` environment variable, but an explicit
`#SBATCH --time` inside the job script **overrides** that env var. So every job
took the script's ceiling: 45-minute compiles asked for 6–12 hours, and every
production run would have asked 12 hours regardless of its real 3.6 h cost.
Slurm cannot backfill a job into a gap shorter than its request, so the whole
chain was queued as if each step were a maximal job.

**Fix.** `sbatch()` now takes `hours=` and emits a COMMAND-LINE `--time`, which
overrides everything. Requests are now: builds 45 min, GS 3 h, vacuum controls
1.5 h, production runs their own `_rt` estimate (6 h for 2560 steps, 12 h for
5120), notebooks 2 h.

**Action taken.** Killed the (idle, polling) orchestrator, cancelled the three
over-requested build jobs, cleared `state.json` (P0 had not completed, so nothing
was lost), patched, relaunched. New jobs 32877570/1/2 at 45 min each.

**Email is knowingly disabled** (user instruction): no Gmail credentials on this
device. `orch.log` is the record; every phase still logs what it would have sent.

---

## 2026-08-05 — DESIGNED, BUILT, TESTED, LAUNCHED autonomously

### What this campaign is

Test the central claim of Nazarov & Gross (arXiv:2510.26222, 2025): at the same
charge and the same velocity, projectiles of different **mass** feel different
friction — and mass acts through exactly one channel, the projectile's spatial
**width**. Measured with a **bath-deposit** stopping power in a dense jellium
slab, which is the only observable that means the same thing for a classical
Gaussian perturbation and for a wavepacket of any mass.

### The design, and the three numbers that forced it

Read the plan for the full derivation. The three constraints that fixed every
parameter, each calibrated against this repo's own prior runs:

1. **Aliasing** `π/h ≥ M·v + 3/(2σ_WP)` — reproduces the `nazarov_gross` guard
   exactly (M = 2.2, v = 2.711, h = 0.35, σ = 0.5 → 8.964 vs π/0.35 = 8.976).
   At σ_WP = 4 and v = 1.07 the coarsest allowed h is ~2, far above the 0.50 used,
   so **the bath sets the grid, not the projectile** — aliasing is a non-issue.
2. **Timestep** `dt ≤ 0.08·min(M,1)·h²`. The `min` is the trap: one dt advances
   all 124 orbitals and the 103 bath states have m = 1, so a heavy projectile buys
   NO speed-up. Cost is flat for M ≥ 1 and 1/M below it. (I had this wrong as
   `0.08·M·h²` in the first draft of the plan; corrected before any run.)
3. **Traversal.** Below the Bragg peak a σ_WP = 0.5 packet at M = 1 stops in a
   few Bohr at every density. The **wide packet is the escape hatch**: at
   σ_WP = 4 the form factor cuts the coupling ~50×, taking deposit/KE from 6.3 to
   0.16. Density then supplies the KE budget — the existing r_s = 5.665 slab gives
   deposit/KE = 1.63 (stops); **r_s = 2.5 gives 0.30 (crosses)**.

### Locked system

r_s = 2.5011 slab, 15 Bohr thick, 30 × 30 × 120 Bohr box, periodicity(2), h = 0.50,
N = 206 → 123 states. v₀ = 1.0742685 a.u. = **1.40 v_F = 0.875 of the Bragg peak**
(below it, as required). σ_WP = 4.0, launch z = −25, two-sided CAP η = −1.0 Ha over
±[45, 60]. Ladder: classical M→∞, classical M = 1, and WP at M = 3, 1.2, 1.0, 0.5.
σ sweep at M = 1 over σ_WP ∈ {2, 3, 6}. ≈ 41 GPU-hours, ≈ 24 GB.

### Files created (all absolute)

| Path | What |
|---|---|
| `.../shared/configs/slab_n206_L30x30x120_rs2p5.hpp` | source-of-truth config; self-checks against closed forms |
| `.../scripts/ng_mass_ladder/gs/run.cpp` | ground state + hard gates (electron count, r_s, finite E) |
| `.../scripts/ng_mass_ladder/wp/run.cpp` | WP half: mass fork, CAP, interactions.csv, checkpoint/resume |
| `.../scripts/ng_mass_ladder/classical/run.cpp` | classical half as a **moving Gaussian perturbation, not a UPF** |
| `.../scripts/ng_mass_ladder/orchestrate.py` | the autonomous chain (P0…P9) |
| `shared/bin/run-ng-gs.slurm`, `run-ng-rt.slurm`, `run-ng-notebooks.slurm` | job scripts |
| `.../hypotheses/ng_mass_ladder/ng_analysis.py` | S extraction, width, verdicts |
| `.../hypotheses/ng_mass_ladder/make_figures.py` | the 7 NG validation figures |
| `.../hypotheses/ng_mass_ladder/build_notebooks.py` | run notebooks + phase notebooks |
| `.../hypotheses/ng_mass_ladder/tests/test_ng_analysis.py` | **10 tests, all passing** |

Also modified: `ResearchProject/systems/vacuum/scripts/wp_selfinteraction/run.cpp`
gained `WP_INV_MASS` (defaults to 1.0, so every previously published run in that
directory is bit-unchanged) so the Phase-1 control can ask whether the LDA
self-interaction error is *mass-dependent*.

### Verified vs unverified

**VERIFIED**
- Config header compiles and every derived constant (r_s, k_F, ω_p, σ_pot, CAP
  region, aliasing headroom, dt ceilings) self-checks against its closed form.
- All four Python modules import; all three SLURM scripts pass `bash -n`.
- **10/10 analysis tests pass**, each against an analytically known answer:
  `extract_S` recovers a planted slope to 0.2 %; the WP bath deposit correctly
  excludes what the projectile owns; the fit window is genuinely restricted to
  the slab (a 10× outside-slope does not leak in); a projectile that never
  crosses returns NaN rather than a plausible wrong number; segments concatenate;
  the pilot gate fails on NaN energy and passes a healthy pair.
- Pre-flight: `sbatch` available, `inq-study` present, 220 GB free on /rds.

**UNVERIFIED — this is what the campaign measures**
- Nothing has been compiled against INQ yet; P0 is the first thing the
  orchestrator does and a compile failure stops the chain with an email.
- The ~5 s/step cost estimate is scaled from p3's points × states, not measured.
- The dt ceiling is an inference from two runs; the binaries enforce it as a hard
  refusal and the pilot's energy drift confirms it.
- Whether the deposit separates between masses at all — that is the result.

### Known risks recorded in the plan §7

Traversal numbers use the packet's *initial* width and the packet spreads (real
deposit is smaller, traversal easier); `stopping_power_sigma` over-suppresses so
all sizing is a bracket; the σ = 6 sweep point is containment-marginal
(4σ = 24 Bohr in a 30 Bohr cell) and is flagged automatically; LDA SIE may impart
a mass-dependent width error, which Phase 1 measures.

### Defect found in passing (NOT fixed, outside scope)

`ResearchProject/systems/localised_jellium/hypotheses/classical_slab_stopping/analyse_classical_baseline.py:98-101`
calls `stopping_power_point(V0, RS)` / `stopping_power_sigma(V0, RS, …)` where the
signature is `(v, kF, …)` and `RS` (line 34) is r_s. At r_s = 4, v = 1.3 the
Lindhard reference comes out **13.17 eV/Bohr instead of 1.95 — 6.75× too high**.
Fix is `L.kF_from_rs(RS)`.

### How to drive it

```bash
cd /rds/user/skcb2/hpc-work/tddft/inq-tddft-research
setsid nohup venv/bin/python3 \
  ResearchProject/systems/localised_jellium/scripts/ng_mass_ladder/orchestrate.py \
  >> ResearchProject/systems/localised_jellium/scripts/ng_mass_ladder/orch.log 2>&1 &
```

Idempotent: `state.json` records completed phases and every run is skipped if its
`run_summary.txt` says `run_completed = true`. Killing and restarting resumes;
a killed RUN resumes from its checkpoint (the orchestrator does this automatically
on a SLURM `TIMEOUT`). Only a build failure, a bad ground state, or a pilot that
cannot resolve a deposit stops the chain — cost overruns are warnings, per
`.claude/rules/checkpoint-dont-block.md`.

Monitor: `tail -f .../ng_mass_ladder/orch.log`, `squeue -u skcb2`, and one email
per phase to chiddukanna@gmail.com.

### Not done

- Nothing committed to git.
- No run has executed yet — the orchestrator was launched at the end of this
  session and P0 (compile) is the first gate.
- `docs/validation/test-catalogue.md` row for the 10 new tests.
- Run-catalogue rows (`tddft-run-catalogue`) after the runs complete.
