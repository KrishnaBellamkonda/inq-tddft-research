# Handover — Nazarov–Gross validation (fixed-velocity projectile-mass sweeps)

Campaign: `/local/data/public/skcb2/tddft/docs/campaigns/nazarov_gross_comparison/nazarov-gross-comparison.md`
(status `ready` → orchestrator launched; frontmatter `done` flags are flipped by
the orchestrator itself as stages complete).

## 2026-07-12 evening — checkpointing added; LAUNCH ON HOLD (user needs GPUs)

**Policy change (user):** never self-block on projected budget overruns —
checkpoint instead, launch full scope, WARN by email; the user kills/resumes.
Codified in `.claude/rules/checkpoint-dont-block.md` + memory
`feedback_checkpoint_dont_block`.

**Code done:** `wp/run.cpp` now checkpoints (NG_CKPT_EVERY=200 default →
`results/<run>/rt_ckpt/` + `rt_state.txt`; `NG_RESUME=1` resumes bit-faithfully,
segment-suffixed CSVs, run_summary.txt always final-state) — cloned from the
proven `sigma1_massonly/wp/run.cpp` pattern. `orchestrate.py`: budget gate →
WARN-and-proceed; both pilots at full 1400-step target; syntax-checked; binary
rebuilt against inq-study.

**STATE: staged, NOT launched — the user needs the GPUs for another task.**
GS35/GS40 + h=0.40 smoke are done and will be skipped on relaunch (idempotent).
To start when the GPUs are free:
```
cd /local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/nazarov_gross \
  && setsid nohup ../../../../../venv/bin/python3 orchestrate.py >> orch.log 2>&1 &
```
Sequence on launch: r1 m0.5∥m0.71 (~7.4 h) → r2 m1.41∥pilot m10 (1400 steps) →
pilot m1 (1400) — total ~19 h across 2 GPUs, checkpoint every 200 steps, emails
per stage. Kill any run with `kill <pid>` (pids in orch.log); resume by rerunning
orchestrate.py (completed runs skipped; a killed run needs NG_RESUME=1 — see the
rule file for the recipe).

## 2026-07-12 03:00 — h=0.35 cost blow-up; RE-PLAN to h=0.40 (user-locked)

**What happened overnight.** GS35 + GS40 both converged (1.8 h;
E₃₅ = −160.90 Ha, E₄₀ = −161.01 Ha, Phase 1a done + email). Guards passed as
predicted. The h=0.35 m=2.2 smoke then EXECUTED but at ~260 s/step effective
(61→87 s/step propagation, rising, PLUS ~30-min stalls on every 10th-step
observable write; init alone 40 min) vs the ~16 s/step extrapolation ⇒ an
880-step null run ≈ 60+ h, 9× over budget. Energy conservation was FINE
(drift in the 7th decimal) — cost, not physics, failed. Evidence:
`scripts/nazarov_gross/wp/smoke.log` (kept). Diagnosis (inference): 24 GB GPU
memory oversubscription at 143³ × 137 states (managed-memory thrash — p3 wrote
the same suite in seconds at h=0.5) compounded by the 143 = 11×13 FFT radix.

**Intervention.** Orchestrator killed at 03:01 BEFORE it could launch the null
rounds (its original design only budget-gated the pilots — gap now fixed).
Partial smoke results removed; the log retained.

**User re-plan decision (2026-07-12):** null branch at **h = 0.40** (GS40
already converged) with only the **3 surviving rungs {0.5, 0.71, 1.41}** —
m=2.2 dropped (aliasing BLOCK at 0.40; user chose NOT to substitute m=1.8).
Pilots unchanged. Fresh 14 h clock at relaunch.

**Orchestrator changes:** smoke now h=0.40/m=1.41 with a HARD budget gate on
the null rounds (measured 2×880-step estimate must fit remaining wall, else
stop + email); round 2 = m1.41 (GPU0) ∥ pilot m=10 at its full 1400-step target
(GPU1, launched only if it fits); pilot m=1 is the budget-sized tail (floor
800 steps, else deferred + email). Campaign file updated to match
(hypothesis mass list, tasks 1b/1c, locked_parameters spacing + re-plan
rationale, guard_rails, budget).

## 2026-07-11 — Phase 1 designed (grill session) + autonomous launch

### What this is
Test of Nazarov & Gross 2025 (arXiv:2510.26222, source note
`/local/data/public/skcb2/tddft/docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`):
same-charge projectiles of different mass at the SAME velocity feel different
friction — but only in the slow (sub-v_F) regime; at high v mass drops out.

### Locked design (user-locked in the 2026-07-11 grill, decision by decision)
- **Baseline:** the p3 run (`scripts/fullsuite_wp/results/p3_wp`, 2026-06-23;
  ledger notebook `hypotheses/qsp_phase2/quantum_stopping_ledger_p3_26-6-26.ipynb`).
  Slab n234, 50³ Bohr, r_s = 3.996 ⇒ v_F ≈ 0.48 a.u., E_F ≈ 3.1 eV. σ_WP = 0.5,
  E = 100 eV ⇒ v = 2.7110633401 = 5.6·v_F, launch −23, dt 0.02, 880 steps, no CAP.
- **Sweep invariant = VELOCITY** (user pick over fixed-k₀/fixed-E).
- **Spacing bound h ∈ [0.35, 0.40]** (user); all null runs at **h = 0.35**
  (single-grid comparison; aliasing mass ceiling m ≤ 2.2 at the 2%-tail tier).
- **Null branch (4 runs, v = 2.711 = HIGH-v regime ⇒ NG predicts FLAT S(m)):**
  m ∈ {0.5, 0.71, 1.41, 2.2}; **m=1 NOT rerun** (user) — p3 anchors it
  (cross-grid caveat h=0.5 vs 0.35 accepted; flatness judged primarily among the
  4 same-grid runs).
- **Slow pilots (2 runs, v = 0.25 = 0.52·v_F, the DISCRIMINATING regime):**
  m ∈ {1, 10} at h = 0.40 (own GS), launch −13.5, steps budget-sized
  (cap 1500 / floor 800, else deferred). They de-risk the Phase-2 slow
  production ladder.
- **Budget:** 14 h wall × 2 GPUs, orchestrator-enforced from its start time;
  pilots are the sacrificial tail.
- **S extraction:** null = the p3 retained-energy ledger method (comparability);
  pilots = initial-drag window (light-projectile rule; NO v-drift gates).
- **Key regime fact (answers the user's question):** v = 2.711 is NOT the NG
  discriminating regime — the null branch is the control half; only the slow
  branch can positively validate the theorem.

### Infrastructure built this session (all new)
- `ResearchProject/systems/localised_jellium/scripts/nazarov_gross/gs/run.cpp`
  — env-spaced slab GS (NG_SPACING, NG_CKPT, NG_COMPILE_PROBE), engine inq-study.
- `.../nazarov_gross/wp/run.cpp` — p3 clone (fullsuite_wp) + fixed-velocity mass
  envs (NG_MASS, NG_V, NG_SPACING, NG_GS_DIR, NG_OUT, NG_N_STEPS, NG_LAUNCH_Z);
  mass via `electrons.inverse_mass()[0][wp_idx]` (the inq-study fork, same call
  as effmass_pair/quantum/run.cpp:81; fork trusted via muon-mass-fork Phases 1–3
  incl. bit-for-bit regression). No CAP.
- `.../nazarov_gross/orchestrate.py` — Python orchestrator: compile gates → GS35
  (GPU0) ∥ GS40 (GPU1) → cutoff guards (verified pre-launch: 5 pass, m=2.2 warn
  at ~1.7% tail ≤ 2%) → 60-step m=2.2 smoke (energy-drift gate < 1e-3 Ha,
  measures s/step) → null rounds (m0.5∥m0.71 then m1.41∥m2.2) with one-shot
  retry → budget-gated pilots → emails per stage + frontmatter flips.
  Idempotent (skips `run_completed = true`).
- New GS checkpoints will land at `shared_gs/slab_n234_L50_h0p35` and `_h0p40`.
- CONTEXT.md: new glossary section "Nazarov–Gross mass sweep (2026-07-11)"
  (null branch, slow branch, fixed-velocity sweep, aliasing mass ceiling,
  spreading systematic).

### Verified vs unverified
- VERIFIED: cutoff-guard verdicts for all 6 rungs (module call, output above);
  email module imports; both GPUs free (~23.5 GB via .gpuprobe); GS/WP binaries
  — see the launch section below for build outcome.
- UNVERIFIED / extrapolated: per-step wall cost at h=0.35 (~16 s/step estimated
  from p3's 5.6 s/step × grid ratio) — the smoke MEASURES it; dt=0.02 energy
  conservation at h=0.35 — smoke gates it; the m=1 slow pilot's extractability
  (that is the pilot's question).
- Expected timeline: GS ~1–1.5 h → smoke ~0.3 h → null rounds ~2×4 h →
  pilots ~4 h ⇒ ~13–14 h total.

### Not done (next session)
- Task 1e analysis: retained-energy ledgers per null run, S(m) figure,
  spreading-systematic check (`wp_real_space_stats`), pilot initial-drag
  verdicts; notebook under `hypotheses/nazarov_gross/`.
- Phase 2 (slow production ladder) — design gated on the pilot verdicts.
- Run-catalogue rows (tddft-run-catalogue) after runs complete.
- Nothing committed to git yet this session (campaign + scripts + CONTEXT.md +
  this handover are new/modified working-tree files).

### Failure modes + where to look
- Orchestrator log: `scripts/nazarov_gross/orch.log` (timestamps + elapsed h).
- Per-run logs: `scripts/nazarov_gross/wp/<name>.log`, GS logs in `gs/`.
- Any stage failure → email to chiddukanna@gmail.com with the log pointer;
  crash → full-traceback email.
- On resume: rerun `orchestrate.py` — it skips completed runs.
