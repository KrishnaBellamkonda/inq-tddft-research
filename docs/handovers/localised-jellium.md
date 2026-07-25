# Handover — Localised jellium implementation + scattering campaign

Rolling handover. Plan: `docs/plans/localised-jellium.md`. Prompt:
`docs/prompts/localised_jellium/localised_jellium_campaign.md`. Theory:
`docs/notes/localised-jellium-theory.md`. Glossary: `CONTEXT.md` → "Localised
jellium".

## 2026-06-27 (COMPLETE + aliasing finding + diagnostics) — Phase 5

**Sweep COMPLETE 11:42:52** (POSTPROC_DONE), ~14 h, fully autonomous: all 6 points
measured, 6 threaded emails sent, 5 run notebooks (81 cells) + study notebook built.

**KEY SCIENTIFIC FINDING — high-v grid aliasing (user-confirmed via E_cut).**
Grid h=0.5 ⇒ k_Nyq=π/h=6.28, E_cut=½(π/h)²=19.74 Ha=**537 eV**. The WP is not
monochromatic (σ_WP=0.5 ⇒ σ_p=1.0), so its momentum TAIL crosses Nyquist even when
the drift E<537. Aliased fraction 1−Φ((k_Nyq−k₀)/σ_p): v3 0.05%, v4 1.1%, v5 10%,
v6 39% — tracks exactly where S blows up + the late energy slope flips POSITIVE
(energy created, impossible under a CAP). Momentum data confirms: v5/v6 n_wp PEAK
jams the top bin 6.23 (not k₀); v3 peaks correctly at 3.4.

**Corrected results (se_state.csv):** trustworthy **E≤122 eV (v≤3): S≈2.4–2.6
eV/Bohr, FLAT**, all upper bounds (23→2.37, 54→2.39, 122→2.57). 218 eV (v4) =4.50
**borderline** (~1% alias). **340 eV→12.95 & 490 eV→18.90 = ALIASED, EXCLUDED**
(lower-bound flag, +slope). Clean criterion k₀+3σ_p<k_Nyq ⇒ E≲146 eV. High-v fix:
finer grid h≤0.35 (ideally 0.25 → E_cut 2150 eV) + new GS.

**Diagnostics added (2026-06-27, user request):** `build_phase5_comparisons.py` →
figs/cmp_{energy,norm,momentum_kz,centroid_sigma}.png + 6 side-by-side
cmp_density_<E>eV.gif. Study notebook `qsp_phase5_study.ipynb` rebuilt (14 cells, 0
err) with the aliasing section (E_cut=537) + cross-run diagnostics + GIF grid +
corrected takeaway.

**S(E) plot now auto-flags aliasing (2026-06-27):** `build_se_plot.py` renders any
point with `bound=="lower"` (positive late energy slope = grid-aliased) as a gray
✗ "aliased (excluded, off-scale)" marker, and sets the y-limit from the PHYSICAL
points only — so 490 eV no longer distorts the plot.

**v5 (340 eV) RE-RUN at finer grid LAUNCHED 2026-06-27 14:58 (PID 1203248).** Fixes
the v5 aliasing on h=0.35 (k_Nyq=8.98, E_cut=1096 eV; v5 tail now ~4σ inside Nyquist).
Infra `scripts/qsp_phase5/rerun_v5_h035/`: `gs/run.cpp` (GS at 0.35 → new checkpoint
`shared_gs/slab_n82_L50x50x90_h0p35`; env LJ_SPACING/LJ_GS_DIR), `wp/run.cpp` (phase-5
WP + LJ_SPACING/LJ_GS_DIR, deeper include paths), `run_rerun.sh` (GS → 20-step
stability smoke → full WP k0=5 τ=40 n=1000 → chain). **analyse_phase5 now reads
`P5_EGS` env** for the grid-dependent GS anchor (the 0.35 GS energy ≠ −70.2257).
Chain overwrites the se_state `p5_wp_v5p0` row (clean replaces aliased), emails the
corrected S(E), builds `p5_wp_v5p0_h035_run_notebook.ipynb`. ETA ~11 h (GS ~1 h +
WP ~10 h; 2.9× grid). Monitor `scripts/qsp_phase5/rerun_v5_h035/rerun.log`.

**RERUN COST CAUGHT + FIXED (2026-06-27 17:33):** first launch (PID 1203248) GS OK
(−71.857 Ha; +44 eV vs 0.5 — anchor matters!) + dt=0.04 SMOKE PASSED, but the full
run was **~56 h**: at h=0.35 compute is ~60 s/step (5×, not 3×) AND each VTI write
is ~430 s (60 MB/frame × ~7 fields every 3 steps). Killed it (incl. the orted/MPI
child), reused the GS, relaunched **LEAN (PID 1270065, rerun2.log, SKIP_SMOKE=1)**:
n=700 (τ=28), LJ_WRITE_EVERY=25, LJ_WF_EVERY=100000 (wavefunction/wp VTIs off — none
needed for the energy-method S). **Revised ETA ~13–14 h** (compute-bound). 28 obs
points + norm + momentum + density_total frames retained. Config now in run_rerun.sh.

**STILL OPEN:** v6 (490 eV) still aliased/excluded (would need same 0.35 treatment +
the same lean config); mark campaign tasks done in the YAML.

## 2026-06-26 (LAUNCHED) — Phase 5: WP quantum stopping power S(E) velocity sweep (AUTONOMOUS)

**Last phase. Designed via grill-with-docs, fully implemented, validated, LAUNCHED
~21:30 (PID 3832333).** Autonomous overnight dispatcher — no user/Claude in loop.

**What it does:** sweeps the σ_WP=0.5 WP drift energy and measures the quantum
stopping power S=[E_total(t_f)−E_GS]/L_z at each, building one **S(E)** curve
overlaid on BULK classical (σ_WP=0.5=σ_q=0.354 `sigma0p35`) + BULK Lindhard +
the lone localised park point (v=2.0, 0.25). Plot rebuilt + **emailed after each
run** (threaded `[lj-wp-se-sweep]`, root sent 21:30).

**Grid (S vs drift E=½k₀²·27.211):** {23, 54, 122, 218, 340, 490} eV ↔
v∈{1.3, 2.0, 3.0, 4.0, 5.0, 6.0}. **54 eV (v=2.0) REUSED from phase 4** (seeded
into se_state.csv as p4_wp_v2p0, S=2.39, upper bound). **5 NEW runs:**
v∈{1.3,3.0,4.0,5.0,6.0}. τ≈200/v (cap 200). ETA ≈10 h / 2 GPUs (value-first:
GPU0 v6→v5→v4→v3 clean first; GPU1 v1.3 alone). 1-GPU fallback 18.6 h still <24 h.

**System reused unchanged from phase 4:** GS `shared_gs/slab_n82_L50x50x90`, box
50×50×90, slab L_z=25, CAP η=−0.7 faces ±35, launch −23.75, dt=0.04, spacing 0.50
(>0.40 ⇒ no WP-init deadlock). Only k₀ (`LJ_K0`) and τ (`LJ_N_STEPS`) vary — single
env-driven binary `scripts/qsp_phase5/wp/run.cpp` (built `INQ_SOURCE=inq-study`).

**Files (all under systems/localised_jellium/):**
- `scripts/qsp_phase5/run_sweep.sh` — dispatcher (smoke gate→2-GPU→per-run chain).
  Launch: `nohup bash run_sweep.sh > sweep.log 2>&1 &`. Smoke fail ⇒ email abort.
- `scripts/qsp_phase5/wp/run.cpp` — env-driven WP binary (`LJ_K0`).
- `hypotheses/qsp_phase5/analyse_phase5.py` — per-run QSP (E_GS anchor, conv gate,
  N-guard) → `results_<tag>.json` + upsert `se_state.csv`.
- `hypotheses/qsp_phase5/build_se_plot.py` — S(E) overlay + threaded email.
  Caches: `classical_sigma0p5_bulk.csv` (6 pts), `lindhard_ref.npz`.
- `hypotheses/qsp_phase5/build_phase5_notebook.py` — study notebook.
- `.claude/skills/run-notebook/run_notebook_builder.py` — +WP energy-method QSP
  section (`--e-gs-ha`, `--l-slab`); softened the loss-fn "indicative only" note.

**Docs:** plan `docs/plans/localised-jellium-qsp-phase5.md`; campaign
`docs/campaigns/localised_jellium/qsp_phase5_velocity_sweep.md` (id
locjel-qsp-phase5-se-sweep, 0/8, INDEX rebuilt); decision **ADR 0010**
(slab-WP-vs-bulk geometry estimate); CONTEXT.md "quantum S(E) sweep" note.

**VALIDATED before launch (on phase-4 WP data):** analyse_phase5 reproduces
2.391 eV/Bohr exactly; classical σ_q=0.354 extraction = 6 pts (peak 0.98 @ v=0.8,
0.25 @ v=3.0; v=2.0 bulk 0.51 vs localised park 0.25 = the geometry gap);
Lindhard cache built; S(E) PNG built; **seed email SENT 21:30** (thread root
`<178250582265...@inqview.gmail>`). All python py_compile-clean; run_sweep.sh
bash -n clean.

**STATUS AT HANDOVER:** smoke gate building/running (monitor `sweep.log`). On
smoke pass → production auto-starts. Expected emails: seed(1pt)→v6.0(2)→v5.0(3)→
v4.0(4)→v3.0(5)→v1.3(6), then study notebook + `POSTPROC_DONE`. Convergence:
high-v (218/340/490) converge to TRUE values; v=1.3 & 54 are UPPER bounds.

**If resuming:** check `sweep.log` + `se_state.csv`. Failed individual runs are
SKIPPED (sweep continues, partial S(E) still emailed). To re-analyse a finished
run: `analyse_phase5.py <results_dir> <tag> <v>` then `build_se_plot.py --email`.

## 2026-06-26 (later) — WP quantum stopping power: independent skill-grounded re-extraction

User asked to "calculate the quantum stopping power for the WP projectile" (phase 4).
Re-derived independently with the `stopping-power-extraction` skill kernels on
`scripts/qsp_phase4/wp/results/p4_wp/raw/observables/observables.csv` —
**reproduces `analyse_phase4` exactly: S_WP = 2.39 eV/Bohr (= 2.4, 2 s.f.).**
- Method = geometry-correct **Method B (slab deposit/thickness)**, BUT with the
  **WP-specific anchor `E_GS`, not `E_total(t₀)`** — the WP's 136 eV kinetic lives
  inside `E_total(0)` (E_total(0)−E_GS=+140.7 eV), so the kernel's default t₀ anchor
  gives a meaningless −0.12 eV/Bohr. The `E_GS` anchor is the right WP adaptation
  (`energy_method()` already did this).
- **Guard PASS:** N_total drain 1.1% (<2% tol); bath-only overflow 0.0023 e (negligible)
  ⇒ the deposit signal is bath excitation, not CAP-on-the-bath artefact.
- **NOT converged → UPPER bound:** norm_f=0.094 (WP 91% absorbed, not 100%); late
  dE/dt=−0.046 eV/au still draining ⇒ residual WP energy will leave ⇒ converged
  S < 2.39. A converged value needs the WP fully absorbed (longer run / stronger CAP).
- Refs: point-Lindhard(54)=0.448 ⇒ **S_WP/LR=5.3×, S_WP/classical=9.6×**. Unchanged.

## 2026-06-26 (DONE 11:44) — Phase-4 S(v=2.0 / 54 eV) COMPLETE — classical redesign worked

Autonomous pipeline ran end-to-end (production 05:46→11:03, analysis+notebooks→11:44),
fully hands-off. **Park FIRED** — the "may not reach ±35" caveat was UNFOUNDED: the 54 eV
ion reflects SHALLOW and reaches the CAP fast (t=25), unlike the deep-penetrating 100 eV.

**Results (results.json, TAG=p4):**
- **WP energy method: S = 2.39 eV/Bohr** (deposited E_total(t_f)−E_GS = 59.8 eV / 25).
  NOT converged (norm_f=0.094 vs 0.02 gate; slower 54 eV WP less absorbed by τ=100;
  worse than p3's 0.046) ⇒ BOUND. Late slope −0.046 eV/au (energy plateaued). t=0 OK:
  E_system(0)−E_GS=+0.24 eV; ⟨T_WP⟩ analytic 5.0 Ha (136 eV) == run 4.996. N_final 82.09.
  wall 5.07 h. **S_WP(54)≈S_WP(100)=2.38 — WP stopping FLAT 54–100 eV.**
- **Classical energy method (NEW park+remove): S = 0.249 eV/Bohr** (deposited 6.23 eV),
  **CONVERGED (ok)**, park t=25.04 z=−36.94, N-drain 0.24%. wall 4.35 h. Old transit
  method reads 0.0 (parked KE=0) — exactly the p3 failure this fixed.
- Point-Lindhard(54 eV)=0.448. **WP/classical=9.6×; WP/LR=5.3×; classical/LR=0.56×
  (classical reflects, below LR).** Classical CAP transmit 0.85 reflect 0.061; ω_p 3.5 eV.

**Notebooks (executed 0 err):** `hypotheses/qsp_phase4/{qsp_phase4_study,p4wp_run_notebook,
p4cl_run_notebook}.ipynb`. `POSTPROC_DONE` 11:44.
**MINOR (non-fatal, pre-existing):** p4cl pipeline `stopping` phase fails `'fz' not in
index` (track has no force col) — notebook built (60 cells); headline classical S is
from `classical_energy_method` (clean). **S(v) overlay NOT yet rebuilt** with the v=2.0
points. NOT committed.

## 2026-06-26 (earlier) — Phase-4: energy-matched S(v=2.0 / 54 eV) point — RUNS COMPLETE, STUDY NOTEBOOK BUILT

Goal: add ONE on-grid S(v) point for the σ_WP=0.5 WP-in-localised-jellium at
**v=2.0 a.u. = 54.42 eV** (an exact classical S(v) grid energy from
`jellium/hypotheses/06_sigma_convergence/sigma_sweep_report.py` →
`sv_convergence_energy.png`; grid E={0.54,4.9,8.7,13.6,23,54,122} eV). Chosen as the
clean on-grid energy with k₀≥2σ_p (σ_p=1/(2σ)=1.0). GS REUSED
(`shared_gs/slab_n82_L50x50x90`, E_GS=−70.22568 Ha). Plan:
`docs/plans/qsp-phase4-energy-sv-points.md`.

**Decisions (grill 2026-06-26):** one energy (54.42 eV), WP + classical, WP
mechanism UNCHANGED. **Classical REDESIGN** (fixes p3 reflection anomaly):
Ehrenfest dynamical ion → when **|z_ion|≥35** (CAP inner face) the ion is PARKED
and **`ions.remove(0)`** deletes it so its Gaussian-Coulomb potential is EXACTLY
zero → projectile-free segment 2 to τ. Classical S via **energy-deposited method
(WP-like) using the stopping-power-extraction skill** (`slab_stopping_power`,
gated). NOTE: at 54 eV the ion likely reflects and may need τ>100 to reach ±35 —
classical τ sized from smoke (default CL_TAU=160).

**Implementation feasibility (verified by source read):** INQ RT restart with
moving ions is unsupported (`propagate.hpp:15`) → two-segment uses FRESH
`propagate(start_step=0)` chunks; state persists in ions/electrons. Callback fires
at iter=0 (carried-over dup) then 1..n (`propagate.hpp:81,123`) → global step `g`,
SKIP iter==0 on chunks after the first, stamp time `g·dt` (not chunk-local
`data.time()`). `ions::remove(long)` exists (`ions.hpp:236`). CAP/background
`zero_step` are no-ops → re-issuing chunks is safe.

**Files BUILT:** `shared/configs/slab_n82_L50x50x90_E54.hpp`
(struct `SlabN82_L50x50x90_E54`, WP_K0=2.0); `scripts/qsp_phase4/wp/run.cpp`
(phase-3 WP, config swapped); `scripts/qsp_phase4/classical/run.cpp` (two-segment
park+remove; env `LJ_CHUNK`=25, `LJ_PARK_Z`=35, `LJ_PARK_STEP`=smoke hook);
`scripts/qsp_phase4/run_production.sh` (WP τ=100 GPU0 + classical CL_TAU GPU1).
**PENDING:** `post_process_phase4.sh`, `analyse_phase4.py`, `build_phase4_notebook.py`
(build during production).

**Autonomy policy (user, hands-off):** validation-gated pipeline, nothing proceeds
on a red gate. Gate 2 (smoke) = compile + mechanism: WP clean injection/no-NaN;
classical ion advances → forced park (`LJ_PARK_STEP=100`) → `remove` fires →
seg-2 clean, CSV/VTI numbering intact. **If classical smoke FAILS: halt classical,
launch WP production autonomously, report** (don't silently change the physics).
Non-converged S = reported as a bound, NOT a failure. On smoke pass: auto-size
CL_TAU, launch both (2 GPUs), run the phase-3-style watcher → study + 2 run-notebooks
→ proactive headline report.

**SMOKE PASSED (2026-06-26 05:40):** compile OK; chunk continuity EXACT (chunk-1
end e=−64.4357217 == chunk-2 start, no seam jump); park fired (`[park] |z|=17.4
step 101`), `ions.remove(0)` dropped E_total −63.51→−70.207 (≈E_GS ⇒ potential
zeroed); segment-2 ran projectile-free; track global-step numbering clean (parked
rows vz=ke=0). `analyse_phase4.classical_energy_method` validated standalone
(stopping-skill, no exception). Note: at 54 eV the ion decelerates HARD (v 2.0→1.17,
KE 54→19 eV) over 6 Bohr while still outside the slab ⇒ reflects shallow.

**PRODUCTION LAUNCHED (2026-06-26 05:46):** `run_production.sh 100 160 0.04` —
WP p4_wp GPU0 (n=2500, τ=100), classical p4_classical GPU1 (n=4000, τ=160,
chunk=25, park|z|≥35). Classical propagating (step 10 e=−64.665, == smoke ⇒
deterministic), WP building. Watcher ARMED (task blingwy3x,
`post_process_phase4.sh` with build-precheck + tightened `done. parked` marker).
ETA ~8 h (classical long pole) + ~1 h post-proc. On done → study + 2 run-notebooks
→ proactive headline report. **PENDING USER NOTE:** if park doesn't fire by τ=160
(possible — phase-3 100 eV took >100 au to return), `classical_energy_method`
flags it; may need longer τ. Phase-3 notebooks rebuilt OK (3 parser bugs fixed).
NOT committed.

**RUNS COMPLETE + STUDY NOTEBOOK BUILT (2026-06-26 11:15).** Both `./run` processes
finished (`run_completed=true`): WP `p4_wp` t=99.84 (wall ~5.1 h), classical
`p4_classical` ran 4001 steps to t≈160 (wall ~4.4 h). A fresh watcher
(`bn0kkr98n`) replaced the stale `blingwy3x` and woke on completion. Then:
`analyse_phase4.py` (TAG=p4) ran clean → **24 figures (12 PNG + 12 GIF) +
`results.json`**; `build_phase4_notebook.py` → **`qsp_phase4_study.ipynb` executed
with 0 errors**. Artefacts in `hypotheses/qsp_phase4/{figs/,results.json,
qsp_phase4_study.ipynb}`.

**Results (all in `results.json`):**
- **WP energy method: deposited E_total(t_f)−E_GS = 59.8 eV ⇒ S_WP = 2.39 eV/Bohr,
  but NOT converged** — norm_f=**0.094** (gate <0.02 FAILS), late dE/dt=−0.05 eV/au,
  plateau 19.5 au. So S_WP is a **deposited-energy upper bound**, not stopping.
- **t=0 bookkeeping CLEAN:** ⟨T_WP⟩ analytic 136.1 eV vs run 135.9 eV;
  E_system(0)−E_GS = **+0.24 eV ≈ 0**.
- **Classical REFLECTED at 54 eV:** park FIRED (t=25.0, `ions.remove` dropped
  E_total to ≈E_GS ⇒ potential zeroed), ion ended z=−36.9, v→0 — never traversed
  the slab, so face-window S=**0.000**; energy-method deposited 6.2 eV (S=0.25, not
  a clean traversal value). Lower energy than the 100 eV p3 ⇒ reflection more likely.
- **CAP (new two-sided, 90-box):** transmit(+z)=0.849, reflect(−z)=0.061 (~6–7%).
- **Cross-phase verdict (the science):** σ=0.5 WP energy-ledger S = {p2 τ40: 2.73,
  p3 τ100/100eV: 2.38, p4 τ100/54eV: 2.39} — **flat, none converged, all ~5× the
  point-Lindhard 0.448 eV/Bohr**. The ledger contains the packet **zero-point KE
  3/(4σ²)=81.6 eV > drift KE ½v²=54.4 eV**, which never leaves the system ⇒
  measures *deposited* energy, not stopping. **Extending τ cannot fix it** (p3 norm
  0.046 → p4 0.094 actually worse at lower v); need larger σ (muon-like) or a
  first-moment ⟨p⟩ observable. §9 takeaway was CORRECTED from the misleading
  "extend τ" line to this grounded explanation.

**Fixes made:** `analyse_phase4.py` cosmetic phase-3→phase-4 relabels (docstring;
Lindhard log line `S(100eV…)`→`S(54eV…)` — numbers were already computed at
V0=2.0, only the printed label was a phase-3 relic).

**PENDING (not done; await user):** (1) per-run deep-dive notebooks
`p4wp/p4cl_run_notebook.ipynb` (parity with p3 via the inqview pipeline runner;
p3's had a known `[fail] layout` on the `50x50x90` cell string). (2) Add the p4
54 eV WP point to the campaign S(v)/S(E) plots (`qsp_phase2/build_sv_plots.py`).
(3) Classical reflection at 54 eV is a real physical result, not a bug — but the
classical S(v) point is not a clean traversal value. (4) NOT committed.

## 2026-06-26 (later) — Phase-3 big-box run COMPLETE: energy ledger + S(v) placement

The τ=100 big-box production twin (50×50×90) finished (t_final=99.84, 314 frames).
Added `p3` to the `RUNS` config; built `quantum_stopping_ledger_p3_26-6-26.ipynb`
(exec 0 err) + updated the S(v)/S(E) plots to show BOTH WP points.

**WP (clean, the headline):** E_GS(90-box)=−70.22568 Ha (GS run_summary). E_total
−63.406→−68.042. **Full ledger ΔE = E_total(t_f)−E_GS = 59 eV → S = 2.38 eV/Bohr.**
E_jellium(0)−E_GS = +0.36 eV ✓. Gate NEAR-MET (norm remaining 0.044 vs 0.02; slope
−0.086 eV/a.u.) — far closer than p2, so 2.38 is close to converged.
**KEY RESULT:** despite 2.5× longer time + 96% absorption, S only fell 2.73→2.38 vs
p2 — the ~8× offset above LR (0.28) is the **zero-point, NOT a convergence artifact.**
Confirms the validated prediction. Drift-credit again negative (−0.87), rejected.

**Classical p3 — ANOMALY (flagged, NOT a stopping number):** ion launched z=−23.75
toward slab, **stalled at z=−3.1 (KE→0) and reversed**, never crossing the far face
(z_final=−14.7); lost ΔKE_ion=82 eV (E_total method dE=78 eV) — **~10× the ~8 eV a
100 eV point charge should lose over 25 Bohr** (bulk classical S≈0.34). Trapped/
reflected; the equal-potential-face method is invalid (no traversal). The compute
guards this (`traversed` flag) and the notebook §5 shows the anomaly box. **TODO:
investigate the p3 classical run setup** (charge sign / potential depth / launch).

**S(v) plots updated** (`sv_plots_26-6-26_figs/`): both r_s=5.67 WP upper bounds now
shown — ▼ p2 (τ=40, 2.7) + ▲ p3 (τ=100, 2.4), barely apart, both ~8× the LR curve
(0.28). p3 classical NOT plotted (anomalous). NOT committed.

## 2026-06-26 — Energy-ledger physics RESOLVED + p5 notebook (full-ledger upper bound)

Deep grill on which energy-bookkeeping argument extracts σ=0.5 WP stopping.
**Independently validated by a fresh scientific subagent** (full agreement). Verdict:

- **The 82 eV is the WP zero-point (confinement) KE** `3/(4σ²)` (uncertainty-principle,
  ∝1/σ²) — NOT a Hartree WP–bath term (that's the tiny ~0.3–0.5 eV `E_jellium(0)−E_GS`).
- **Master eq (CAP = only sink):** `E_total(0)=E_total(t_f)+E_CAP_removed`. At full
  absorption the WP's 185 eV splits: CAP carries out, bath keeps `E_total(t_f)−E_GS`.
- **Full ledger** `S=(E_total(t_f)−E_GS)/L_z` is the ONLY internally-consistent number:
  **p2 = 2.7, p5 = 4.7 eV/Bohr — UPPER BOUNDS** (zero-point-inflated + not fully absorbed).
- **Drift-credit / "zero-point-persists" REJECTED:** subtracting only the 100 eV drift
  gives p5 = +1.5 (≈ user's 32 eV) BUT **p2 = −0.5 eV/Bohr (impossible: bath below its
  own GS)**. The sign-flip with absorption fraction proves the zero-point does NOT
  persist in the final jellium (it leaves with the absorbed packet). Decisive.
- **First-moment vs second-moment insight:** stopping lives in ⟨p⟩ (drift), zero-point
  in Var(p); the energy ledger sums both → contamination baked in. zero-point/drift ∝
  1/(m²σ²v²) = 0.82 here (worst case). Clean routes: **large-σ / muon** (zp→0) or a
  drift-momentum/force observable — but force/trajectory is UNAVAILABLE (WP+bath are one
  inseparable KS system, no Exact Factorization in INQ), so energy ledger is the route.
- **User decision (2026-06-26): full ledger, as upper bound.** S(v) WP point stays 2.7
  (p2, r_s=5.67) on the existing plots (already a plain marker after annotation removal).

**p5 ("14–20 a.u." = N=234, r_s=4, τ=18) notebook built** — same analysis as the p2
ledger, parametrized. `quantum_stopping_ledger.py` now has a RUNS config (p2/p5);
`build_quantum_stopping_notebook.py [p2|p5]`. Outputs: `quantum_stopping_ledger_26-6-26.ipynb`
(p2, rebuilt) + `quantum_stopping_ledger_p5_26-6-26.ipynb` (NEW), both exec 0 err, both
show the full-ledger headline + the rejected drift-credit row for transparency.
p5 E_GS=−160.992 Ha (static-run log), SIE=4.5 eV, E_jellium(0)−E_GS=+0.52 eV ✓.
NOT committed.

## 2026-06-25 (later 5) — 26-6-26 meeting: quantum-stopping energy ledger + S(v) placement

Grill-with-docs session → built the careful quantum-stopping deliverable for the
26 Jun supervisor meeting. Plan: `docs/plans/quantum-stopping-ledger-26-6-26.md`.
All artefacts in `ResearchProject/systems/localised_jellium/hypotheses/qsp_phase2/`.

**Locked method (user's retained-energy definition):**
`S_WP = [E_total(t_f) − E_jellium(0)]/L_z`, `L_z=25 Bohr`;
`E_jellium(0) = E_total(0) − ⟨T_WP⟩ − E_SIE` (compare to E_GS as a consistency check,
NOT subtracted in S). Data = `p2_wp`/`p2_classical` (τ=40, reused; Phase-3 incomplete).

**Computed ledger (p2_wp, all verified live in the notebook):**
- E_total(0)=−38.94 Ha; ⟨T_WP⟩(0)=6.645 Ha (run-measured `e_kin_ha`, =181 eV: drift
  100 + zero-point 81); E_SIE=4.40 eV=0.162 Ha.
- **E_jellium(0)=−45.749 Ha vs E_GS=−45.759 Ha → +0.26 eV ✓ (self-consistent).**
- E_total(t_f=40)=−43.253 Ha → ΔE=+68 eV → **S_WP=2.7 eV/Bohr — UPPER BOUND ONLY.**

**KEY FINDING (why 2.7 is NOT physical, user flagged it should be comparable to LR):**
the τ=40 WP is unconverged — residual WP orbital STILL carries ⟨T_WP⟩(t_f)=68 eV ≈ ΔE,
gate fails (norm rem. 0.13≫0.02; E_total slope −1.1 eV/a.u.≠0). The bath's +68 eV is
dominated by zero-point dumping (82 eV) + un-absorbed remnant, NOT velocity-stopping —
the documented "σ=0.5 WP stopping unmeasurable via E_total" (qa_v, notes_campaign1).
For a comparable (≈LR) number E_total must drain to ≈−45.2 Ha (≈τ=90 → the Phase-3 run).
NO residual correction reaches the comparable zone from τ=40 data.

**Classical slab (p2_classical, user: "don't expect it to be right"):** ion NOT absorbed
(transits, wraps, re-enters by t=40); E_total method inapplicable. ΔKE two ways:
slab-centre KE-min ("lowest energy point", user choice) S≈3.0 eV/Bohr (over-counts via
conservative well); equal-potential faces (|z|=12.5, 67→55 eV) S≈0.53 eV/Bohr (true).

**Artefacts:** `quantum_stopping_ledger.py` (computation + 2 diagnostic figs),
`build_quantum_stopping_notebook.py` → `quantum_stopping_ledger_26-6-26.ipynb`
(executed 0 err, live tables), `build_sv_plots.py` → 4 figs (S(v)/S(E) × Plot A/B).
Plot A = LR + classical-bulk-σ0.5 + WP-100eV (drawn as upper-bound down-arrow + guide
band to the comparable zone); Plot B adds the classical-slab point (3.0, "not trusted").
New rule `.claude/rules/number-rounding.md` (2 s.f. default / 3 s.f. max).

**LR reference (user correctness check):** electron count verified — 82 background e,
intact at t_f (N_total 83.0→82.13 = 82 bath + 0.13 residual WP). The slab IS at a
formally different density (n0=1.312e-3, r_s=5.667, kF=0.3387) vs the bulk S(v) runs
(1.296e-3, r_s=5.690, kF=0.3373) — but only **0.41% in kF** (matched by design per
slab_n82_L50x50x70.hpp). LR rebuilt at the SLAB kF (0.3387); shifts LR@v=2.711 only
0.279→0.282 — does NOT explain the 10× WP excess. Source: Lindhard 1954 / Lindhard &
Winther 1964 (docs/sources/stopping-power-formulae.md); finite-slab applicability
Quijada 2007 (velocity-gated, holds at v=2.711).

**STATUS:** notebook + ledger + all 4 S(v)/S(E) plots DONE & verified (exit 0), LR now
at slab density. Classical-bulk σ0.5 S(v): 0.37/0.84/0.97/0.94/0.51/0.25 at v=0.2–3.0
(Bragg peak v≈0.8); at v=2.711: bulk-classical ≈0.34, slab-LR=0.282, WP upper-bound=2.72.
NOT committed. **OPEN:** user reviewing the 2.7-upper-bound WP representation — may want
it repositioned to a converged/comparable estimate (needs Phase-3 to full absorption,
or a zero-point-free estimator; neither clean from the τ=40 data).

## 2026-06-25 (later 4) — Phase 3 big-box run: seam CAP engine-fix + setup (IN PROGRESS)

Designing/launching the next-phase run (big-box / long-time σ=0.5 WP+classical pair,
energy-method stopping). Grilled extensively; two background validators (energy-method
+ CAP/plasmon) returned (findings folded in below). **Decisions locked:**

- **Geometry:** box **50×50×90** (z∈[−45,45]); slab |z|<12.5 unchanged; launch z=**−23.75**
  (equidistant — 11.25 Bohr to both the slab face and the CAP inner face → kills the
  t=0 init-absorption); τ=**100 a.u.**; spacing 0.50; ETRS; LDA; inq-study.
- **CAP = TWO-SIDED** (FINAL decision 2026-06-25). User initially chose a single
  seam-centred CAP ("A"); I implemented + Python-validated an `inq-study/absorbing.hpp`
  wrap-fix for it (min-image distance; NEW==OLD for non-seam Δ≤1.4e-15; seam peaks at ±45,
  continuous) — then the user **reverted to the two-sided CAP** ("work with a known devil")
  because it carries the 1.3%-reflection benchmark. So: `inq-study/absorbing.hpp` **restored
  byte-identical to `inq/`** (wrap-fix removed; `inq/` was never touched), and the run.cpp
  uses the original two bumps `cap_lo(η,-40/90,10/90)+cap_hi(η,+40/90,10/90)` → region
  **[±35,±45]**, η=−0.7, 10 Bohr/side (same η + per-side width as the benchmarked 70-box
  CAP, just moved to the bigger box's edges). Reflectivity benchmark carries over. NOTE:
  `rvector` returns reduced coords [−0.5,0.5) (constructor clamps mid_pos there); the qsp
  build's `inq_SOURCE_DIR=inq-study` (CMakeCache-confirmed).
- **Energy method (user's, corrected):** S = [E_total(t_f)−E_GS]/L_z. INQ `energy_total`
  is the clean real ⟨H₀⟩ (CAP does NOT contaminate it — validator verified). Subtract the
  WP **kinetic** ⟨T_WP⟩≈6.675 Ha (NOT eigenvalue ε_WP → would double-count SIE) + SIE(4.40
  eV) at t=0 only; compare to E_GS (expect E_GS + few-eV cross-Hartree, NOT exact; existing
  τ=40 gives E_GS−0.5 eV ✓). **User DEFINES the jellium system = density REMAINING in box**;
  CAP-absorbed = transmitted/reflected WP + secondaries (ledgered as a DIAGNOSTIC, NOT added
  back). ⟨p_z⟩-loss (Method A) REJECTED by user (KS orbital ≠ physical WP). Convergence gate:
  WP norm<0.02 AND energy plateau.
- **Reflectivity caveat (honest):** the 1.3% reflection + cap_sweep curves were for the OLD
  two-bump / one-sided geometries, NOT the seam-centred split. Seam profile is strictly
  gentler (same 20-Bohr absorber, peak-at-seam, 10-Bohr monotonic onset) → expect ≤1.3%,
  to be MEASURED in-run (WP reaches +z CAP at ≈21.7 a.u.).

**Files created:** `shared/configs/slab_n82_L50x50x90.hpp` (LZ 70→90, rest identical);
`scripts/qsp_phase3/{gs,wp,classical}/run.cpp` (gs = copy of phase1; wp/classical = copies
of phase2 with seam CAP `absorbing(η,0.5,20/90)`, launch −23.75, GS dir →
`shared_gs/slab_n82_L50x50x90`, fixed provenance string). Schematics:
`docs/notes/qsp_bigbox_run_schematic.png`, `qsp_cap_topology_comparison.png`.

**90-box GS DONE:** `shared_gs/slab_n82_L50x50x90`, **E_GS=−70.22568 Ha**, 82 e, 61 states,
r_s=5.667. The −24.5 Ha shift vs the 70-box GS (−45.759) is PURELY electrostatic (hartree
+external; kinetic Δ−0.006, xc Δ+0.007 → identical electronic structure) — the charged-slab
-in-PBC finite-size constant that **cancels exactly in E_total−E_GS** (same cell). Benign;
never compare absolute E_GS across box sizes.

**dt=0.04 smoke PASSED** (150 steps, two-sided CAP, both runs): energy smooth/monotonic,
no NaN, ended normally through slab entry; WP norm 83.000→82.995 (**init-absorption FIXED**
by the equidistant launch); E_total(0)−E_GS=+185.5 eV (matches 70-box) and run ⟨T_WP⟩≈6.64
Ha (≈analytic 6.675). ~6 s/step on the 90-box.

**Config validated** by an agent: READY TO RUN, all params match design, **propagator = ETRS**
(not CN — verified `real_time.hpp:183` default), CAP two-sided/engine pristine. Campaign
updated: frontmatter task **P3.1** + Phase-3 LOCKED decisions block. Stale config-header CAP
comment fixed (two-sided). (Minor: run.cpp top-of-file comment blocks still carry phase-2 copy
text — slab_n234_L50, old launch defaults — cosmetic, code correct.)

**STATUS — PRODUCTION RUNNING + AUTONOMOUS POST-PROC ARMED** (launched ~21:20, task bv2rosqz3,
`run_production.sh 0.04`): WP (GPU0, results/p3_wp) + classical (GPU1, results/p3_classical),
dt=0.04, 2500 steps = τ 100 a.u., write_every=8. WP is the slower (~8.5 s/step) ⇒ **ETA ~6 h**.

**Analysis pipeline BUILT + tested** (on the 150-step smoke via P3_TAG=p3v):
- `hypotheses/qsp_phase3/analyse_phase3.py` — 90-box geometry + NEW energy-method diagnostics:
  S=[E_total(t_f)−E_GS]/L_z, the **convergence triple** (norm<0.02 AND |dE/dt| AND plateau),
  the **t=0 cross-term sanity** (smoke: ⟨T_WP⟩ run 180.8 ≈ analytic 181.6 eV; E_system(0)−E_GS=−0.47
  eV ✓), the **CAP-absorbed-energy ledger** (WP/bath split). Two pre-existing P2 bugs FIXED:
  overlap_full parse (`comment="#",header=None` → KS off-diag 0.000→0.215) and efield API
  (`electric_field(...).ez`). `P3_TAG` env switches smoke(p3v)/production(p3).
- `build_phase3_notebook.py` — study notebook assembler (tested, 0 errors).

**AUTONOMOUS WATCHER LAUNCHED** (task bhn7qgxwe, `scripts/qsp_phase3/post_process_phase3.sh`,
log `post_process.log`): polls for both runs' "done.", then FULL analyse_phase3.py (with GIF
battery) → shrink GIFs (470px/64-colour) → build_phase3_notebook.py → both run-notebooks
(p3wp/p3cl, cap-inner 35, launch −23.75) → writes `POSTPROC_DONE`. Fully hands-off. (Fixed a
`grep -c "||echo 0"` double-count bug that would have blocked completion detection, before relaunch.)

**Still pending (NOT auto):** corrected energy diagnostics on the existing τ=40 runs; the P2.3
(large-σ) + P2.4 (classical removal) follow-on runs. Schematic `docs/notes/qsp_bigbox_run_schematic.png`
= two-sided CAP.

## 2026-06-25 (later 3) — notebook presentation fix: GIF display size + titles

Presentation-only (no physics numbers / `results.json` touched). Two problems the user
hit reviewing `qsp_phase2_study.ipynb`:

- **Figures rendered too large.** Cause: both embed helpers emitted Markdown `![]()`
  images → native pixel size. Fix: `embed(path, caption, width=None)` in
  `hypotheses/_nbreport.py` now emits an HTML `<img src width>` when `width` is set;
  the run-notebook builder's own `img()` (`.claude/skills/run-notebook/run_notebook_builder.py`)
  auto-sizes by extension (GIF 360 px, PNG 520 px). Study-notebook widths: GIF 360,
  multi-panel PNG 600, single PNG 520 (constants `W_GIF/W_PNG2/W_PNG1` in
  `build_phase2_notebook.py`).
- **Cryptic GIF titles** ("wp · total · n"). Cause: baked-in title in
  `make_density_gif_battery` used terse codes. Fix (shared, in
  `inq-stack/python/inqview/visualisation/density_gifs.py`): spelled-out `CAT_TTL`/
  `KIND_TTL` ("Total system", "density n(x,z,t)", …) + new `run_title` param with a
  prettifier (`wp→"Wavepacket run"`, `classical→"Classical run"`); `run_label` still
  keys filenames. Run-notebook builder passes `run_title` from `rtype`.

Regenerated GIFs (titles are baked into pixels), re-applied the PIL downscale pass
(target width 470 px, 64-colour) since regen re-inflated them, and rebuilt all three
notebooks: **study 46 cells / 21 `<img>`**, **p2wp 79 / 56**, **p2cl 60 / 42**, all
**0 exec errors, 0 missing refs, 0 leftover native-size images**. GIF totals:
study `figs/` 7.1 MB (max 1.2), p2wp_figs 9.1 MB (max 1.6), p2cl_figs 3.4 MB (max 1.2).
Verified `density_gifs` imports with `run_title`. `inq/` untouched; nothing committed.

## 2026-06-25 (later 2) — classical KE(z) dip-recovery + plots_examples additions

- **Classical transport plot** (`figs/classical_transport.png`, new `classical_transport()` in
  analyse_phase2): **z(t)** + **KE(z)**. Confirms (from track) the projectile KE is NOT
  monotonic — **dips to 30.3 eV at the slab CENTRE (z=+1.3), recovers to 54.8 eV at +12.5** =
  a **conservative mean-field well** (energy borrowed/returned, NOT stopping). Only the net
  loss between **equal-potential** points (symmetric faces ±12.5) is true stopping ⇒ **S=0.507
  eV/Bohr**; centre/asymmetric windows would over-count. **Window choice dominates** — record
  it. Added to study notebook **§6a** + wired into the run-notebook builder (classical branch).
  Added `SKIP_GIFS=1` env to analyse_phase2 for fast PNG-only re-render (reuse existing GIFs).
- Study notebook rebuilt (46 cells, 0 err). Both run-notebooks rebuilt earlier (WP 77, cl 56,
  0 err); classical re-rebuilding now to pick up the transport plot.
- **plots_examples.md NEW items → SKILL.md:** (a) **KS eigen-energy bar GIFs** incl/excl WP +
  Δ(E_i(t)−E_i(0)) — generator = `inqview.pipeline.state_energies`/`state_energy_spectra`
  (example `…/run_wp_n162_L50_E100_sigma1/.../ks_energies_delta_no_wp.gif`); NOT yet emitted for
  qsp_phase2 (phase needs enabling/wiring — FOLLOW-UP). (b) **S(v)/S(E)** curves = sweep-level
  → study notebooks, not single-run.
- **Particle removal — user now prefers PARKING over deletion** (asked why deletion disrupts).
  Rationale recorded: deletion's problems are (1) sudden quench — the screening cloud rebounds
  when the ion's potential vanishes abruptly → spurious collective oscillation/energy kick;
  (2) **periodic charge-neutrality** — removing the projectile's charge unbalances the cell ⇒
  ill-defined G=0 Hartree (implicit compensating jellium) [VERIFY whether the Gaussian pseudo-ion
  carries net charge in this setup]; (3) stock INQ can't delete an ion mid-propagate anyway.
  Parking (freeze v at +35, keep charge) — BUT see corrections below.

  **2026-06-25 (later 3) — CORRECTIONS after reading the projectile UPF + engine:**
  - **Neutrality is NOT a problem (I was wrong to flag it).** The projectile UPF
    `electron_gaussian_sigma0p35.upf` has **z_valence = 0** — it is a pure external Gaussian
    potential mimicking a −1 charge's field, contributing NO valence charge to the electron
    count / Hartree G=0 balance. So zeroing the projectile's charge does NOT unbalance the cell.
    The user is right: **zero-charge is clean here.**
  - **Parking-with-charge is NOT clean (correcting my "negligible at 22.5 Bohr").** The
    projectile potential is a long-range **Coulomb tail** (erf(r/√2σ)/r), NOT Gaussian-
    suppressed: at r=22.5 Bohr it is **~1.2 eV** at the slab (plus periodic images). So a
    frozen-but-charged parked ion still perturbs the slab by ~1 eV ⇒ **zeroing the charge is
    physically preferable to parking** (matches the user's instinct).
  - **Feasibility:** velocity-parking is trivially run.cpp-only (`real_time::propagate` takes
    `ions&` by reference — run.cpp's step lambda mutates the real object). Cleanly **nulling the
    ionic POTENTIAL** mid-run is harder (it's the ionic pseudopotential, rebuilt each Ehrenfest
    step from ion position) — needs a probe to find a hook (zero/scale the species, or an
    inqkit potential-gate). **NEXT: probe scouts (a) velocity mutation honoured? (b) can the
    ion's Gaussian potential be nulled past z≥35?** Then implement zero-charge if feasible, else
    park + quantify/accept the ~1.2 eV tail. τ=40 both.

## 2026-06-25 (later) — P2.1 notebooks upgraded: GIF battery, heuristics, looping dx

Big review-driven upgrade of the P2.1 deliverables + reusable library + skill.

**Physics question answered (data-grounded): is the projectile looping back?**
- **Classical YES** — the CAP absorbs *wavefunctions*, not the classical point charge, so
  the Ehrenfest ion is never absorbed: track z −22→**+62.7** (box [−35,35] ⇒ wrapped to
  physical −7.3) → re-approaches slab → E_total **rises after t≈30** (expected/physical).
- **WP NO** — E_total **monotonic down** (0→−117 eV); N_total 83→82.13; WP orbital norm
  1.000→**0.136** (86% absorbed at +z CAP); centroid advances to +8.66 then the survival-
  weighted mean drifts back (artefact, not reversal). The only WP quantity that *rises* is
  the system kinetic energy (0→+6.6 eV in first ~8 a.u. = bath excitation), then falls.

**New reusable library (inqview), with passing test:**
- `inqview/visualisation/density_gifs.py::make_density_gif_battery` — the
  **{n, Δn=n(t)−n(0), Δn=n(t+dt)−n(t)} × {total, wp, bath}** GIF matrix (bath = n_total −
  n_wp). WP→9 GIFs, classical→3 (Δn = induced wake). Density GIFs share a log scale across
  total/bath (low densities visible; WP total now uses the classical-total scale per
  request); Δ kinds symmetric diverging. Loads via `load_vti` (physical order, NO fftshift).
- `inqview/analysis/heuristics.py::compute_heuristics` (+ building blocks) — groups A–I:
  HEG scales (k_F=0.339, v_F=0.339, E_F=1.56 eV, λ_F=9.28, **ω_p=3.49 eV, T_plasmon=48.9
  a.u. > τ=40 ⇒ plasmon under-resolved**, k_TF=0.657), timescales (**t_exit slab END =
  12.73 a.u.**, t_enter 3.50, box-edge 21.0), zero-point KE 81.6 eV, spreading ×41.3,
  norm/absorption (WP 0.865 e: orbital 0.864 + bath overflow 0.001; classical from track),
  Lindhard refs. Test: `inq-stack/tests/python/inqview/analysis/test_heuristics.py` (5 pass).

**analyse_phase2.py upgraded:** 12-GIF battery (shared classical→WP density scale),
heuristics block → results.json, **slab-exit-time markers** (t=12.73) on conv/classical/
energetics, energetics CAP-param header + "monotonic drain vs periodic re-entry" subtitles,
new **norm_absorption.png** (N_total + WP-orbital norm vs time). Reran (~15 min; GIFs are
large, 5–27 MB each, path-referenced).

**Study notebook** `qsp_phase2_study.ipynb` rebuilt: §4 now the full GIF battery (4a–4d) +
§4e energetics + **§4f looping-back explanation** + §4g norm/absorption; new **§8.5
physical anchors/heuristics** tables; §5 N(t) framing corrected (starts at 83, not 82).

**Run-notebook SKILL + builder made standard (user-authorized skill edit):**
- `run_notebook_builder.py`: auto-renders the density-GIF battery (`make_density_gif_battery`)
  + a **Physical anchors & heuristics** section (`compute_heuristics`) for every run-notebook.
- `SKILL.md`: battery table now includes the 9-GIF density battery, KL metric, momentum-
  difference, heuristics, energy-decomp CAP/slab annotations; references
  `docs/notes/plots_examples.md` as the plot spec (LEED = coronene-only; 2D loss Fourier-gated).
- Run-notebook **CAP display fixed**: added verified `cap_eta_ha=-0.7`, `cap_width_bohr=10`,
  `rs`, `v0_au`, `launch_z_bohr` to both `run_summary.txt`; corrected the stale "η−0.5" string.
- Both run-notebooks **rebuilding in background** (battery+heuristics) — verify via
  `p2wp_runnb.log`/`p2cl_runnb.log` (0 errors expected).

**STILL PENDING (user-gated, NOT started):**
- **Classical-particle removal on wrap-around — DESIGN ACCEPTED 2026-06-25, NOT yet built.**
  Locked spec: **remove the classical ion at z ≥ +35 (box edge)** by **mutating the (const)
  ions** (`StepContext.ions` is `const*` — needs const_cast in the run.cpp step lambda, or a
  non-const ion hook in `inqkit::real_time_session`; NEVER `inq/`). **Keep BOTH runs at the
  same τ = 40 a.u.** (do NOT shorten the classical run — user explicit). Goal: classical
  late-time state = relaxing bath with no projectile, comparable to the WP's absorbed state.
  IMPLEMENTATION PLAN (next): (1) **5-step probe** — does INQ `real_time::propagate` honour a
  mutation of run.cpp's `ions` mid-step, or iterate an internal copy? (2) if honoured →
  run.cpp-only: in the step lambda, once z_ion≥35 set velocity=0 each step (+ zero force) — at
  +35 the ion & its periodic image are ~22.5 Bohr from the slab, charge overlap ~e^(−(22.5/0.35)²)≈0
  ⇒ effectively removed without the Hartree step that true charge-deletion would cause; if NOT
  honoured → add a non-const ion hook to `inqkit::real_time_session`. (3) re-run classical P2.1
  with the gate, re-analyse, confirm late-time energy no longer rises. True ion deletion / charge-zeroing
  mid-run is NOT supported by stock INQ (Hamiltonian built around the fixed ions list) — the
  park-at-edge is the faithful equivalent. code-test + simulation-validation gates apply (new run).
- plots_examples.md report-style figures (system-setup, GS-excitation-decomp, 2D loss) not
  yet added as generators — documented in SKILL.md as the spec.

## 2026-06-25 — Phase-2 P2.1 COMPLETE: both runs + 3 notebooks built

Phase-2 P2.1 (WP+classical convergence/CAP test, 2000 steps / 40 a.u.) **finished
overnight** (`P2_1_DONE` 2026-06-24 20:29, both exit 0; WP wall 7217 s, classical 6836 s).
Autonomous dispatcher ran build→run→analyse cleanly; only the notebook step failed —
`build_phase2_notebook.py` didn't exist yet. **Now built (this session):**

- **Study notebook:** `hypotheses/qsp_phase2/qsp_phase2_study.ipynb` (28 cells, executed 0
  errors, all fig refs resolve) via new `build_phase2_notebook.py` (uses `_nbreport`,
  path-referenced figs/GIFs, harvest-before-rebuild).
- **Run notebooks:** `hypotheses/qsp_phase2/p2wp_run_notebook.ipynb` (63 cells, type=wp) +
  `p2cl_run_notebook.ipynb` (50 cells, type=classical), via the run-notebook skill builder
  (`--cap-inner 25 --rs 5.666 --launch-z -22 --v0 2.711`; classical `--measured-s 0.018632`).
  Figs travel beside in `*_figs/`.
- **Dispatcher wired:** `scripts/qsp_phase2/dispatch.sh` tail now also builds both run-notebooks
  (auto-build convention); added `PYTHONPATH`.

**P2.1 RESULTS (all PROVISIONAL — 40 a.u. test, from `results.json`):**
- **Q1 WP convergence = NO** — late slope **−1.03 eV/a.u.** at t_f, N_final 82.135; deposited
  68.2 eV ⇒ S_WP **2.73 eV/Bohr is an UPPER BOUND only**. P2.2 must lengthen τ.
- **Q2 classical steady state = YES** — KE loss across slab ⇒ S_cl ≈ **0.507 eV/Bohr**
  (vs point-Lindhard 0.282). PROVISIONAL pending `stopping-power-extraction` skill.
- **Q3 CAP = GOOD** — reflect 0.013 (~1%), transmit 0.865; WP init clean (room to orthog).
  **Keep η=−0.7, 10 Bohr.**
- **Dominant surprise: spreading ×41** at σ=0.5 ⇒ comparison dispersion-dominated →
  motivates the large-σ no-spread run.

**Provenance/analysis bugs to fix before P2.2 (logged in study §9):**
1. `run_summary.txt` CAP label stale ("η−0.5"); **compiled value is η−0.7** (run.cpp:69,124).
   Reformatted the two summaries to one-pair-per-line (kept `.orig`) so the parser reads them;
   fix the run.cpp summary writer (one-per-line + correct η) for P2.2.
2. **KS off-diagonal weight = 0** (overlap_full parse/index bug, not physics).
3. **E-field block skipped** (`inqview.analysis.efield` API mismatch; density VTIs saved → recoverable).
4. Classical KE mass/velocity-unit convention — re-extract via stopping-power skill.

NOTE: tried to harden `run_notebook_builder.py::parse_summary` for multi-pair lines — **denied**
(skill self-modification guard); fixed the data instead. **NEXT:** user reviews the 3 notebooks →
author P2.2 production (longer τ + large-σ run) + fix the 4 bugs. NOT committed.

## 2026-06-24 (later 6) — Phase-1 COMPLETE: GS validated + SIE = 4.40 eV

Both Phase-1 GPU runs done (exit 0). Results:
- **GS (P1.1):** `shared_gs/slab_n82_L50x50x70` — converged, **E_GS = −45.75885 Ha**, 82 e, r_s=5.667,
  orthorhombic 50×50×70; interior density mean 0.00136 (≈ n0), **5.6% flat** (Friedel+surface).
  Closed-shell check deferred (non-blocking). Density fig `hypotheses/qsp_phase1/gs_density_xz.png`.
- **SIE short-RT (P1.2):** WP injected far at z_mean(0)=**−32.000**, norm=1.000, σ_z(0)=0.353 (=σ_WP/√2);
  **E_total(0) = −38.95203 Ha**; **KE_WP = 6.64498 Ha = 180.82 eV** (drift 94.2 + zero-point 86.6).
- **SIE both ways (P1.3) — diagnostic validated:** **SIE_a (E_GS+100eV) = +85.22 eV**;
  **SIE_b (clean, −measured KE) = +4.40 eV = THE SIE**; difference +80.82 eV = KE_WP−100 = the
  zero-point the "+100 eV" omits. SIE≈4.40 eV matches old r_s=4 (~4.5 eV) — weakly density-dependent. ✓
- **Notebook:** `hypotheses/qsp_phase1/phase1_gs_sie.ipynb` (5 cells: GS validation + SIE both-ways).
  Scripts: `analyse_sie.py`, `make_gs_density_fig.py`, `build_phase1_notebook.py`. Runs:
  `scripts/qsp_phase1/{gs,sie}/`.
- **Campaign status → paused, 3/4** (P1.1–P1.3 done; P2 open). INDEX refreshed.
- **NEXT:** user analyses Phase-1 notebook → author Phase 2 (CAP 10/η−0.7, sim-time, steady-state
  for BOTH runs, stopping). The 4.40 eV SIE floor carries into the WP−classical comparison. NOT committed.

## 2026-06-24 (later 5) — Phase-1 EXECUTION underway (GS building/running)

Running Phase 1 of `quantum-stopping-power`. Created: `shared/configs/slab_n82_L50x50x70.hpp`
(orthorhombic 50×50×70, N=82, r_s≈5.665, WP launch −32, CAP-off for SIE);
`scripts/qsp_phase1/gs/run.cpp` (GS, orthorhombic cell via `systems::cell::orthorhombic`,
checkpoint → `shared_gs/slab_n82_L50x50x70`); `scripts/qsp_phase1/sie/run.cpp` (WP-far short RT,
10 steps, CAP off, writes observables.csv energy_total + wp_momentum_stats e_kin_ha).
GPUs free (2×25 GB, cudaMemGetInfo). GS build+run launched on GPU0 (bg log
`scripts/qsp_phase1/gs/gs_build_run.log`); compiled clean to 100%, SCF running. NEXT: validate GS
(energy, closed-shell, num_electrons=82), then build+run SIE on GPU, then compute SIE_a/SIE_b +
zero-point cross-check, then Phase-1 notebook. NOT committed.

## 2026-06-24 (later 4) — Phase-1 campaign authored: `quantum-stopping-power` (GS + SIE)

Ran `/campaigns` on Campaign 1. Authored **`docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md`**
(id `jwps-quantum-stopping-power`, status draft, 0/4 tasks; indexed). Locked decisions:
- **Method asymmetry:** classical stopping = Ehrenfest ΔKE_ion via the user's
  **`stopping-power-extraction` skill** (localised-slab branch = ΔE_total/L_z + convergence
  gate); WP stopping = total-energy balance (needs E_total convergence — the make-or-break).
- **Density r_s≈5.67** (matches n162-in-50³ for S(v) comparability); **82 background e**
  (even; nearest closed-shell at GS validation); slab **|z|<12.5**; box **50×50×70** (x,y=50
  preserves in-plane density; z→70 vacuum gives CAP + far-launch room). NOTE: 162 e would FILL
  the box — a localised slab is a sub-volume, so ~81–82 e, NOT 162. Region layout (z): slab
  [±12.5] · free [±12.5,±25] (12.5 Bohr each) · CAP [±25,±35] (10 Bohr each).
- **SIE diagnostic BOTH ways:** SIE_a = E_tot(0)−(E_GS+100eV) [=SIE+zero-point], SIE_b =
  E_tot(0)−E_GS−KE_WP(⟨p²⟩/2) [=SIE]; cross-check SIE_a−SIE_b ≈ zero-point 3/(4σ²)=81.6 eV.
  Explained to user WHY "+100 eV" over-counts (omits the 82 eV zero-point KE). WP+slab energy
  needs RT injection (moving WP ≠ GS eigenstate); few steps, CAP off.
- **Phase-1 tasks (locked):** P1.1 GS slab; P1.2 WP+slab-far short RT (E_tot(0)+KE_WP); P1.3
  SIE both ways + cross-check. **Phase 2 DEFERRED** to next brainstorm (CAP/box/sim-time/
  steady-state-for-both/stopping).
- **Phase-2 CAP (resolved):** **10 Bohr each** ([±25,±35], inner edge ±25), **η=−0.7**, in
  the 50×50×70 box → 12.5 Bohr free region each side, ample far-launch room. (User revised from
  the earlier 16.5-Bohr/60-box that left no launch room. Current baseline CAP was 7.5 Bohr each;
  user had misremembered as 12.5 = the slab half-width.) WP Phase-1 launch z≈−32.
- Phase-2 still to specify next prompt: sim-time, steady-state for both runs, stopping.
- NOT committed.

## 2026-06-24 (later 3) — Campaign decisions + brainstorming notebook + 3-campaign rough draft

New campaign workspace: **`docs/campaigns/jellium_wp_stopping/`**. Sanity-check agent
(general-purpose) validated the design; key outcomes:
- **Muon quantum WP infeasible in stock INQ** — KS-orbital kinetic mass is m_e-hardwired
  (`inq/src/hamiltonian/ks_hamiltonian.hpp:202`); needs an `inq-study` per-orbital-mass fork
  (well-scoped). Classical-muon trivial (ionic mass tunable). Muon = Campaign 3 (FUTURE).
- **Campaign 1 (σ_WP=0.5) E_total ledger compromised** (72× free-spread, stalled centroid →
  Δz undefined, no-wrap vs full-absorption incompatible in 50 Bohr, zero-point 82 eV + SIE
  4.5 eV). User: NOT blockers, proceed; recorded in `notes_campaign1_sigma05_restrictions.md`.
  Resolution: smoke-test a force/work-integral estimator; subtract zero-point + bound SIE.
- **Campaign 2 (large rigid σ) sound** — agent recommends σ_WP≈4 @ E≥300 eV (spread ≤12%,
  ZPKE 1.3 eV, SIE 0.6 eV, fits 50-Bohr box). σ_WP=3 UPF ready; σ_WP=4 needs generation.
  Caveats: matched-σ does NOT fully isolate quantum-ness (WP has Pauli+SIE, ghost neither);
  CAP retune for faster packet. Specifics deferred to Campaign-2 brainstorm.
- **σ decision finding (validated, flips the approach).** Built
  `make_sigma_lindhard_comparison.py` → `sigma_lindhard_comparison.png`. **Point-Lindhard is
  trustworthy** (0.717 eV/Bohr @ r_s=4,v=2.71 ≈ baseline 0.719); **σ_WP=0.5 classical = 0.706
  = 0.98× point** (premise correct). BUT the **analytical finite-σ Lindhard
  `stopping_power_sigma` OVER-suppresses** (predicts 0.77× for σ_WP=0.5 vs the 0.98× run) →
  **cannot judge σ_WP=1.0**. ⟹ σ=1 vs σ=0.5 **undecidable from existing data; needs one cheap
  σ_WP=1.0 classical S(v) run** vs point-Lindhard. √2 trap noted: bulk campaign labels σ_pot,
  localised labels σ_WP.
- **Artifacts (all in `docs/campaigns/jellium_wp_stopping/`):**
  `brainstorming-jellium-campaigns.ipynb` (5 cells, the session record + plot),
  `draft_campaigns.md` (rough draft, 3 campaigns), `notes_campaign1_sigma05_restrictions.md`,
  `make_sigma_lindhard_comparison.py` + png/csv, `build_brainstorming_notebook.py`.
- **Next:** run `/campaigns` skill on **Campaign 1 = "quantum-stopping-power campaign"**. NOT committed.

## 2026-06-24 (later 2) — Brainstorm: literature (Nazarov-Gross EF), SIE quantified (4.5 eV), spread-tradeoff figure

- **Literature (web).** Found **Nazarov & Gross 2025, arXiv:2510.26222 "Stopping power of electron
  liquid for slow quantum projectiles"** — state-of-the-art on EXACTLY this campaign. Source note:
  `docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`. Key: (1) treats projectile
  fully quantum; (2) solves the energy-partition problem via **Exact Factorization** (Ψ=χ(R)Φ(r))
  — the rigorous answer to our §5 "WP+bath inseparable" (INQ does NOT implement EF → our
  total-energy ledger is the approximate route); (3) **mass/width-dependent friction is a purely
  quantum effect** — validates the WP-vs-classical premise AND matching WP σ to classical-potential
  σ; (4) classical (M→∞) limit recovers Lindhard; (5) their projectile is an EIGENSTATE width (does
  not spread) vs our free Gaussian (spreads) — "rigid" for us = large-σ, not truly dispersionless.
  Field also openly states the total-energy method "no longer applies" for quantum projectiles.
- **SIE quantified — NO new run needed.** Realised `p3_wp` (LJ_CAP=0 defaults launch_z=−23) IS the
  user's intended "WP far from slab" run. `SIE = E_total(0)[far] − E_GS_slab − KE_WP`, with **KE_WP
  = measured ⟨p²⟩/2 = 6.645 Ha = 180.8 eV, NOT 100 eV** (the user's "+100 eV" omits ~81 eV
  zero-point+transverse → would overcount SIE by ~81 eV — flagged + corrected). Result: far (−23)
  excess **+4.55 eV**, near (−15.5) +5.02 eV ⟹ **SIE ≈ 4.5 eV**, WP–slab Hartree only ~0.47 eV over
  7.5 Bohr. SIE ~1/s so larger-σ WP → smaller SIE (σ=3 → <1 eV). Recorded as notebook **§5.1**.
- **Added the user's spreading-tradeoff figure** (`fig_spread_tradeoff.png` from
  jellium/hypotheses/cap_baselines → copied as `ref_spread_tradeoff.png`) to §7 with the spreading
  law + `E_min(σ,f)` inversion — the "central difficulty" reminder. Notebook → **35 cells**.
- **Open decision (next):** the σ–energy strategy (Q3) is now informed by literature + SIE. User
  proposed two run families: (1) bounded-spread rigid σ (≤20% spread over transit, σ from E) to
  isolate pure-QM vs classical at MATCHED σ; (2) small-σ + matched classical (same initial σ) to
  isolate spreading effects — user trusts the small-σ stopping more vs Lindhard. Both to be
  analysed/decided, plus density (r_s) and higher-energy choices. NOT committed.

## 2026-06-24 (later) — Next-campaign grill begins: centroid-stall finding + momentum-anchor scheme

New `grill-with-docs` session to brainstorm the **next jellium-slab campaign** (aim: WP vs
classical stopping power through the slab; "ruthlessly empirical"). Phase 1 = extract remaining
learnings from the p5_wp / p5_classical baselines before designing.

- **Reflection diagnostic — data inventory.** Established what the existing runs *can* answer:
  `wp_momentum_stats.csv` already stores **per-axis momentum second moments** (`px2_mean,
  py2_mean, pz2_mean`) → ⟨k_z²⟩ vs ⟨k_⊥²⟩ ARE on disk (pz_mean 2.63→1.71, ⟨k_⊥²⟩ 2.00→2.15
  transverse heating, e_kin 6.65→4.31 Ha). `momentum_distribution.csv` is **1-D |k| only (no
  sign)**. `wavefunction_wp/` has **91 complex VTI frames** (Δt=0.2) → signed `J_z` and full
  3-D `ψ(k)` are computable. Reflection is answerable from existing data **only via the signed
  real-space current** `J_z`, not from sign-blind moments.
- **KEY FINDING — the WP is NON-RIGID; the centroid stalls.** From `wp_real_space_stats`:
  `⟨z⟩(t)` (survival-weighted) rises from −15.5, **stalls at ≈ +5.1** (max +5.6) from t≈9, and
  **never reaches the far face +12.5**; `d⟨z⟩/dt` collapses 2.64→−1.28 (drifts backward late).
  `σ_z` balloons **0.37→11.7 (≈32×)**, so the ±3σ envelope spans the whole box and `⟨z⟩−3σ_z`
  runs backward to −30. ⟹ the user's proposed **"centroid ± 3σ crosses each face" frame-timing
  scheme is mostly undefined** (only the leading edge `⟨z⟩+3σ_z` cleanly crosses: near −12.5 at
  t≈0.4, far +12.5 at t≈4.0). The stall is a **survival-weighted absorption signature** (fast
  forward components removed by +z CAP → slow/back remainder pulls ⟨z⟩ back), NOT physical
  stopping. Classical projectile marches ballistically to +30 (decel 2.71→2.35 = real signal).
- **Built `qa_viii_centroid_trajectory.py`** → `qa_viii_centroid_trajectory.png` (2 panels:
  ⟨z⟩±σ_z trajectory + centroid velocity, WP vs classical). Added **§4.6** to the review
  notebook with the non-rigidity observation, and recorded a **run-independent momentum-snapshot
  anchor scheme (A0–A8, event-anchored not time-anchored)** in §7 — because the next runs change
  jellium density + WP energy, so fixed times don't transfer. Rebuilt notebook → **32 cells**
  (verified: §4.6, anchor scheme, centroid PNG, cell-17 verbatim all present).
- **User steer:** don't over-fit specific snapshot times now (next runs differ); record the
  anchor *scheme*. The full 3-D `ψ(k)` momentum-scattering analysis (sign-resolved n(k_z),
  n(k_⊥), 2-D map, survival ratio testing "high-|k_z| captured first") is **queued for the next
  runs** per the recorded scheme — not yet computed on the baselines.
- **Reflection finding (qa_ix, signed J_z across all 4 internal planes) → §4.7.** In the
  wrap-free window the backward-side planes show **no significant backward current**: −z CAP edge
  (−17.5) integrates to **+0.033 e forward** (net forward, NOT backward), largest instantaneous
  backward excursion anywhere **−0.008**. Forward budget: near face +1.01, far face +0.81, +z CAP
  +0.74 (matches qa_ii). ⟹ left-free refill (§4.4) reads as **spreading tail, not a reflected
  wave**; **near-face reflection negligible**. ⚠️ But **far-face reflection is unmeasured** (would
  reach −z CAP only at t≈22, after the 18-a.u. run) → strongest argument for longer box/run +
  sign-resolved n(k_z) next. Verdict left to the user (presented neutrally). Built
  `qa_ix_cumulative_current_regions.py` (+png/csv); notebook → **34 cells**.
- **Open / next:** continue Phase-1 grill (cumulative current per region across all internal
  planes was offered; reflection-window test; equilibration/sim-time; observable suite incl.
  E-field + loss function; Fourier + stopping-power training as prerequisites), then brainstorm
  every sim parameter, smoke-test, write the `/campaigns` prompt. NOT committed.

## 2026-06-24 — Round-3 audit: KS-orbital framing + robust total-density absorbed count

User re-reviewed `qa_jellium_slab_baselines.ipynb`, re-pasted the original 10-todo prompt,
and asked: (a) cross-check all todos delivered; (b) ensure their in-notebook comments are
preserved verbatim; (c) **stop framing the 0.62 as "0.62 of the wavepacket"** — in the WP run
slab and WP are not cleanly separable, so 0.62 must be stated as the **WP KS-orbital norm**;
(d) verify the **total absorbed = n(t)−n(0) of the total density**.

- **Audit:** all 10 todos confirmed present (table in the chat). Verbatim §8 quotes match the
  original prompt. **Found one user-authored cell (live #17)** — their observation that the
  classical/WP **initial energies are not comparable** (Hartree term; initialise at a distance
  where the electron–slab Hartree term matches; expect SIE to push `E_WP(0)−E_cl(0)` past
  100 eV) — that was **NOT in the build script** and would have been lost on rebuild.
- **Verified the robust measure** directly from `observables/electron_number.csv`:
  `N_total(0)` is integer-clean (**235** WP / **234** classical); absorbed = `N_total(0)−
  N_total(final)` = **0.8309 (WP)** / **0.2156 (classical)**. Decomposition reconciles to
  machine precision: WP-orbital loss 0.62193 + bath 0.20897 = 0.83090. Classical total-density
  absorbed (0.2156) **is** the bath overflow (no WP orbital to separate).
- **Edits to `build_qa_jellium_slab_baselines.py` (then rebuilt → 30 cells):**
  (1) folded user cell #17 in verbatim (after §5 prose) + a response; (2) §4.1 rewritten to
  **lead with the decomposition-free total-density count** (0.831/0.216) then the KS-orbital
  decomposition behind a ⚠️ "KS orbital is not a clean physical separator" caveat; (3) every
  "0.62 of the wavepacket / WP-self / the WP electron itself" → "**WP KS-orbital norm**"
  (§4.1 caption, §5, §7, §8 Todo 2/3); (4) §7 req #6 extended with **matched-Hartree launch**.
  Verified post-rebuild: cell-17 verbatim preserved, §8 quotes intact, dangerous framing gone,
  all 9 figures present.
- **Open:** await user re-review of the updated notebook, then grill Q9 (CAP L/η + box) or
  brainstorm the next campaign from §7. NOT committed (user hasn't asked).

## 2026-06-23 (later 2) — Q1 absorption/transmission/reflection analysis + review notebook

Grill-with-docs Q&A (user's 9-question list). Built the Q1 deep-dive on the p5_wp /
p5_classical CAP runs to inform the NEXT simulation. All in
`ResearchProject/systems/localised_jellium/hypotheses/03_cap_stopping/`.

### Artifacts (all reproducible, path-referenced)
- `qa_iii_absorbed_norm.{py,png,csv}` — cumulative absorbed norm. **WP combined 0.831 =
  WP-self 0.622 + bath overflow 0.209**; classical bath overflow 0.216, projectile 0
  (transmits). Sum-check machine-precision. Resolves the handover's conflicting numbers
  (0.83 = N_total absorbed; 62% = the WP itself).
- `qa_i_region_densities.{py,png,csv}` — projectile charge in 5 z-bands vs t (WP from
  density_wp; classical reconstructed 1-unit Gaussian at the tracked ion z, periodic-wrapped).
- `qa_ii_per_cap_flux.{py,png,csv}` — per-CAP flux from the complex `wavefunction_wp`
  current. **T+R=0.698 vs bookkeeping 0.700 (validated).** Transmission ~0.743.
- `qa_iv_region_detail.{py,png,csv}` — norms overview (twin-axis), bath Δ-per-band,
  reflection diagnostic (free-region). t=0 baseline: **9.7 e bath spill-out** between slab
  and CAPs + 1.0 WP launched in left-free.
- `qa_v_stopping_energy.{py,png}` — Q3/Q8 stopping ledger. **Two baselines diverge by
  85.9 eV ≈ the WP's analytic zero-point KE 81.6 eV** (3/(4σ²), σ=0.5): Formula 1
  (final−init+100eV) = 32.2 eV (S 1.29); Formula 2 (final−E_GS) = 118.1 eV (S 4.72).
  **σ=0.5 WP stopping is UNMEASURABLE via E_total** — zero-point (82 eV) ≈ drift (100 eV),
  WP only 62% absorbed. Classical 0.706 eV/Bohr stands. Q8 naive check fails (−612 eV;
  structural differences). → next campaign needs **larger σ (≈3, zero-point ~2 eV) + full
  absorption**.
- `build_qa_jellium_slab_baselines.py` → **`qa_jellium_slab_baselines.ipynb`** (**29 cells,
  9 figures**, valid JSON, refs resolve) — the user's decision-support review notebook,
  comprehensive: §1 geometry (Q2) · §2 classical projectile (Q4) · §3 WP injection (Q5) ·
  §4 absorption/transmission/reflection (Q1,Q6) · §5 stopping energy ledger (Q3,Q8) · §6
  loss function (Q7) · §7 next-sim requirements · **§8 round-2 reader observations
  (verbatim + validations)**. Renamed from qa_reflection_absorption.
- **Round-2 additions (2026-06-24, user review):** `qa_vi_total_norm_compare.{py,png}`
  (both runs start at 235 — classical electrons + reconstructed projectile; WP loses 0.831,
  classical 0.216 bath-only); `qa_vii_xz_logdensity_final.{py,png}` (final-frame xz log
  density, region lines — spread WP vs compact classical); `qa_v` extended with the **+100 eV
  drift-credited curve**. §5 gained the **energies table** (GS −160.99 Ha/−4380.8 eV; WP t=0
  −154.16/−4194.3, +185.9 above GS; classical t=0 −131.65/−3582.1, +798.3 above GS) and a
  **corrected Formula-2 explanation**: F2=118 eV is dominated by the **residual un-absorbed
  WP (~70 eV)**, not deposited — it equals stopping only at full absorption (user flagged the
  prior "zero-point inflation" wording as wrong). §8 records todos 2/3/4/5/6 verbatim with
  answers (62% = WP-orbital norm over whole box; t0 = first RT frame; reflection needs
  k_z-sign momentum; "WP transmission" = WP orbital current only).

### Key physics findings
- **Bath overflow is run-independent** (~0.21 both runs) — a CAP property. The entire
  WP-vs-classical absorbed-norm difference is the WP electron itself (0.62).
- **Transmission dominates** (~0.74 forward flux, validated).
- **Reflection UNRESOLVED.** Left-free refills to ~0.092 (peak t≈9) = candidate reflection,
  but **indistinguishable from the trailing edge of the 23×-spread WP** with 1D density; the
  −z-edge flux is contaminated by periodic wrap-around (WP reaches +z box edge at t≈14.9).
  Do NOT read reflection as either negligible or large.

### Methodology findings (carry to next campaign)
- **WP current MUST use a spectral (FFT) z-derivative** — finite-difference undercounts the
  k₀=2.71 current by ~28% at dx=0.5 (T+R 0.42→0.70).
- **1D density cannot separate reflection from spreading** → need momentum-resolved (k_z
  sign / 2D k_z–k_⊥).
- **Periodic wrap-around contaminates reflection** → need a longer z-box or earlier stop.

### Next-simulation requirements (the four, for Q9)
(1) longer z-box / stop < t≈15; (2) momentum-resolved (k_z) observable; (3) fine-cadence
flux/current (spectral-deriv or CAP-edge flux screen); (4) CAP L/η tuning from the ε(E,L) maps.

### Grilling status (9-question list)
DONE: Q2 (geometry — periodic-3D bulk vs the localised 2D-periodic slab, verified periodic
in x,y), Q3+Q8 (stopping via E_total referenced to slab GS −160.99 Ha, sign-fixed; run both
formulas, gap = WP self-energy), Q7 (loss function — both deferred until Fourier training),
Q1 (this section), Q6 (answered by qa_iii). PENDING: Q4 (classical projectile / extra_electrons
construction), Q5 (WP init-absorption — norm data can answer), Q9 (CAP L/η for next campaign).
Then brainstorm next campaign → `campaigns` skill.

## 2026-06-23 (later) — √2 width validation + plot/notebook honesty + exact-match mandate

Grill-with-docs follow-up. The σ-label √2 trap was independently validated and the
plots/notebooks relabelled to the **actual** projectile widths.

### Validated (independent fresh-context agent, code+math, file:line cited)
- **WP density std = σ_WP/√2.** `wavepacket.hpp:234,254` builds ψ∝exp(−r²/2σ²) ⇒
  |ψ|²∝exp(−r²/σ²) ⇒ density std σ/√2 (numerically 0.5→0.35355). The generator
  `inqview/io/gaussian_psp.py:126` encodes the same: `sigma_charge = sigma_wp/√2`.
- **The slab pair is NOT exactly matched — ~1 %.** The WP run injects σ_WP=0.5
  ⇒ density std **0.35355**. The classical run loads the **legacy**
  `electron_gaussian_sigma0p35.upf`, whose charge std fits empirically to **0.350**
  (filename label accurate; it is the old charge-std convention, a rounded stand-in
  for 0.354). So WP 0.354 vs classical 0.350 — physically negligible for stopping,
  but the clouds are not identical. Only the WP equals 0.354.

### Done (option iii: relabel existing + exact-match for the future)
- **Matched UPF generated:** `ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf`
  via `generate_gaussian_psp(sigma_wp=0.5)` → charge std **0.353553** (validated: fitted
  charge_std 0.353553, written V matches C·erf(r/(σ_c√2))/r to 5e-11). For FUTURE runs;
  existing slab data NOT re-run (user decision).
- **Plots relabelled / Gaussian-matched curve dropped** (user: compare to point-charge
  Lindhard ONLY — σ chosen to sit in linear-response regime, measured 0.706 ≈ point 0.719):
  - `.claude/skills/run-notebook/run_notebook_builder.py`: new `lindhard_mode`
    (`both` default preserves generality; `point` draws only point-charge + point-only
    narrative). Arg `--lindhard {both,point}`.
  - `06_sigma_convergence/sigma_sweep_report.py`: legend label σ→**σ_q** + legend title
    "σ_q = charge std (σ_WP=√2 σ_q)" + SIGMAS comment. Figures regenerated.
  - `03_cap_stopping/build_cap_stopping_report.py`: setup table → WP density std 0.354 /
    classical charge std 0.350; √2 paragraph rewritten with the 1 % provenance caveat +
    exact-match mandate. Study notebook rebuilt (0 errors).
- **Run notebooks patched** (point-only Lindhard, honest widths): `p5wp_run_notebook.ipynb`
  (cell ae86c2b4, WP density std 0.354) + `p5cl_run_notebook.ipynb` (cell 62708807,
  classical charge std 0.350, 1 % note). `lindhard_stopping.png` regenerated point-only
  for both (S_point 0.719 eV/Bohr). All 3 notebooks valid JSON; no stale bracket text.

### MANDATE (carry to next campaign)
**Every future quantum-vs-classical run uses an EXACT-matched UPF** generated with
`generate_gaussian_psp(sigma_wp=…)` (unified convention), e.g.
`electron_gaussian_wpsigma0p5.upf` — never the rounded legacy `sigmaXpY` files for a
matched pair. Classical charge std must equal WP density std = σ_WP/√2 exactly.

### Rebuild recipe (run notebooks, to stay point-only on a future full rebuild)
`PYTHONPATH=…/inq-stack/python venv/bin/python3 run_notebook_builder.py <results_dir> <out.ipynb> --rs 3.995 --proj-sigma <0.354|0.350> --measured-s 0.025945 --measured-v 2.711 --lindhard point …` (plus the original `--cap-inner --decomp-prefix --launch-z --v0 --run-cpp` args).

## 2026-06-23 — Notebook refresh: decomposition GIFs, CAP lines, Lindhard, exit-time

All notebooks refreshed to the full-suite data + new analyses (user requests this
session). Key artefacts:
- **Decomposition renderer** `hypotheses/_density_views.py`:
  `render_decomposition_views` (WP runs → total/bath/wp × n/Δfirst/Δprev = **9
  GIFs**; bath = total − |ψ_WP|², exact) and `render_classical_views` (classical →
  total e-density + **projectile Gaussian** built from the ion track, σ_pot=0.35 →
  **6 GIFs**). Both take `cap_inner` → dashed **lime** CAP-boundary lines at
  z=±17.5; slab faces cyan dotted at ±12.5. GIFs pre-rendered into
  `hypotheses/<sweep>/<run>_decomp_*.gif` (p3wp 9 no-CAP, p5wp 9 +CAP, p5cl 6 +CAP).
- **`run-notebook` builder extended** (`.claude/skills/run-notebook/
  run_notebook_builder.py`): `--decomp-prefix` (embeds the 9/6 decomposition GIFs in
  a "Three-way density decomposition" section), `--cap-inner` (CAP dashed lines on
  carpets + lead GIF via `_zlines`), `--rs/--proj-sigma/--measured-s/--measured-v`
  (analytical **Lindhard stopping** panel via `inqview.analysis.lindhard_elf.
  stopping_power_sigma`/`_point`), `--launch-z/--v0` (**WP ballistic exit-time**
  plot from `wp_real_space_stats.z_mean`).
- **User cell-13 comments on p5wp ADDRESSED** (not lost — extracted then answered):
  (1) WP exit time — t_exit=(12.5−launch_z)/v0 = **10.33 a.u.** for p5wp; plotted vs
  the measured centroid. (2) Stopping vs analytical Lindhard for the **slab density
  r_s=3.995**: ω_p=5.90 eV, v_F=0.480; the curves BRACKET the measurement —
  Gaussian(σ_q=0.354) S=0.542, point S=0.719, measured(classical ΔKE/x) **S=0.706
  eV/Bohr**. Measured sits ~30% above matched-Gaussian, near point-charge.
  **GAP UNRESOLVED** — candidate causes (none confirmed, flagged as inference in the
  notebook): on-grid projectile harder than σ=0.35 at dx=0.5; nonlinear Z=−1 beyond
  RPA; finite-box/CAP/transient. Verify via dx 0.5→0.25 + Z-scaling. NOT a
  Fourier-from-run-data analysis → loss-function gate NOT triggered (analytical
  theory curve). Source: `docs/sources/stopping-power-formulae.md` (Lindhard–Winther
  1964; Correa 2018).
- **Study notebooks refreshed** (02/03): total-only GIFs → 9/6 decomposition; "pending
  re-run" markers DROPPED; loss-gate + Δω≈9 eV caveat KEPT; deep-dive links added.
  02 = 30 MB, 03 = 41 MB (heavy — all decomposition GIFs embedded), 0 errors.
- **Campaign master refreshed**: notebook map (6 notebooks), full-suite headline S
  from `electron_track.csv`, work-energy corroboration note. 235 KB.
- **Run notebooks** p3wp/p5wp/p5cl rebuilt with decomposition + CAP + Lindhard +
  exit-time (p3wp/p5wp 72-cells; p5cl 57). Loss gate respected (no `spectral_weight`
  in runner PHASES). NOTE: p5wp decomp briefly missed on one rebuild (transient glob
  race) — re-run fixes it; verify decomp_section present after any rebuild.
- **PATH-REFERENCED notebooks (2026-06-23).** User couldn't view the run notebooks
  — they were 25–41 MB (base64-embedded GIFs) and choked viewers. Switched ALL
  notebooks to RELATIVE-PATH markdown image refs: builder `img()` uses
  `os.path.relpath(fig, out_ipynb.parent)`; `_nbreport.embed()` + `set_outdir()`
  same. Run notebooks converted in place (image code-cells → `![](rel)`), study/
  master rebuilt. Now ALL 7 notebooks are **10–22 KB**, 0 errors, 0 missing refs.
  CAVEAT: figures must travel beside the notebooks — decomposition GIFs in
  `hypotheses/<sweep>/`, carpets/Lindhard/exit PNGs in `<stem>_figs/`, pipeline
  figures referenced via `../../scripts/.../results/.../analysis/`. In-place: all
  resolve. The 7: `localised_jellium_campaign_study.ipynb`, `01/slab_validation_
  study`, `02/projectile_slab_study`, `03/cap_stopping_study`, `02/p3wp_run_
  notebook`, `03/p5wp_run_notebook`, `03/p5cl_run_notebook`.
- **LINEAR|LOG 2-panel density plots + readable colorbars (2026-06-23).** User:
  every density plot = linear + log side by side; colorbar gradations fixed across
  the GIF; tick labels readable (2 s.f. + `×10ⁿ` offset — the edges were clipping
  Δn numbers). Implemented:
  - `_density_views.py::_gif` → 2-panel (linear left / log right) single figure;
    positive fields LogNorm(floor vmax/1e3), signed Δn SymLogNorm(linthresh=max/100,
    RdBu_r); norm computed ONCE (fixed across frames); `_readable_cbar` = ScalarFmt
    useMathText + powerlimits((0,0)) + MaxNLocator(5). All 24 decomposition GIFs
    regenerated 2-panel (notebooks path-referenced → auto-pick-up, no rebuild).
  - `run_notebook_builder.py::carpet` + `lead_gif` → 2-panel linear|log +
    `_readable_cbar`; run notebooks rebuilt to refresh those.
  - `report-figures` SKILL.md: **rule 8** (readable colorbar, no edge-clip) +
    **rule 9** (every density map ships linear AND log) added; rule 7 (fixed
    gradations across frames) now actually enforced via single fixed `norm`.
- **ANNOTATION PRESERVATION — harvest-before-rebuild (2026-06-23, plan-approved).**
  User will annotate notebooks directly; rebuilds regenerate from scratch and wiped
  their cells. Implemented in BOTH builders (`run_notebook_builder.py` and
  `hypotheses/_nbreport.py`): every builder cell stamped `metadata.gen="builder"`
  (+ `anchor` slug for headed markdown) via `tag_builder()`; `build()` calls
  `harvest_user_cells(out)` BEFORE regenerating (collects untagged markdown) then
  `reinject()` splices each back at its anchor (tagged `gen="user"`, round-trips
  forever); orphaned annotations → "📌 Carried-over reader annotations" section
  (never dropped). Transition guard: a pre-tagging notebook (no gen tags) harvests
  as nothing → must rebuild once to tag before annotating. Sidecar `<stem>.notes.md`
  pin added to `_nbreport` too (parity). Skills documented (run-notebook +
  notebook-making SKILL.md). **VERIFIED end-to-end** (master study NB: reader cell
  survived rebuild, count=1 no-dup, tagged user, at anchor) + unit-tested
  harvest/reinject + orphan carry-over. Plan:
  `.claude/plans/streamed-sprouting-wolf.md`. Baseline rebuild: study 01/02/03+master
  DONE (tagged, 0 spurious); run notebooks p3wp/p5wp/p5cl rebuilding (bg).
- **3 NEW user analysis notes captured** in `hypotheses/03_cap_stopping/
  p5wp_run_notebook.notes.md` (sidecar, future work — see plan appendix): (1) CAP
  reflection vs absorption %/boundary (needs plane screens at z=±17.5/±25 + flux;
  `PlaneScreen` exists, unused); (2) equilibrated-system quantum stopping
  S=(E_final−E_initial)/dx — needs a LONGER run (≥30 a.u.; p5wp 18 a.u. NOT
  equilibrated, WP 62% absorbed); (3) non-spreading WP run design (σ=1.0+E=300 eV →
  ~5× spread vs 23× now). NOT YET DONE.
- **Two run-machinery TODOs for future runs:** classical `run.cpp` should log `fz`
  (force) directly (currently post-added as m·dv/dt); WP overlap should write to
  `raw/overlap/` (currently symlinked from `raw/observables/overlap/`).

## 2026-06-22 — Grill: observable suite, mapping fix, full-suite re-run

Grill-with-docs session. Decisions:
- **D1 (observable set = default).** The full analysis suite — density VTIs
  `total`/`system`/**`wp`** + `gs_system`, `density_delta` (induced) + a NEW
  Δ-vs-previous-step view, `momentum_distribution`, `overlap`/`overlap_full` +
  `gs_projected_occupations`, `eigenvalues`/`occupations_vs_time`/
  `state_energy_spectra`, and total/WP/bath energy decomposition — is the
  **default-required** set for every run, OPT-OUT only. To be baked into a
  SLIMMED `tddft-simulations` skill + supporting hooks/rules. Extends the
  already-approved 2026-06-15 set in `docs/observables/minimum-set-spec.md`
  (SoT `inqkit/observables/minimum_observable_set.hpp`, ADR 0006).
- **D2 (re-run APPROVED).** Re-run localised-jellium WP (Ph3) + both CAP runs
  (Ph5) on the canonical FULL-SUITE template
  `ResearchProject/systems/jellium/scripts/cap_baselines/run.cpp` (emits
  density_total/system/wp/delta, overlap, momentum, state energies). My bespoke
  minimal runs wrote only `density::total` → cannot reconstruct WP/bath split.
- **D3 (analyse current data NOW).** Proceed with current-data analysis (the
  reconstructable subset: Δn on TOTAL density vs first/prev step, total-density
  momentum FFT) — it unblocks the user. MUST fix the GIF bug first (see below).
- **D4 (standing task).** After the re-runs land, RE-RUN THE ANALYSIS CHAIN on
  them with the full suite (three-way decomposition GIFs, overlap excitations,
  per-component energies, loss function).
- **D5 (loss-function GATE).** `L(q,ω)` IS wanted from the scattering re-runs
  (NOT deferred to a dedicated long run), BUT: (1) the frequency-resolution limit
  (Δω ≈ 2π/T ≈ 9 eV at T≈18 a.u. — cannot resolve the ~6 eV r_s≈4 plasmon) must be
  stated **in bold, clearly, in the notebook**; (2) the loss-function build step is
  HARD-GATED behind a future task — *the user training Claude to conduct Fourier
  analysis correctly* — and may not be implemented before that. See memory
  `feedback-fourier-loss-function-gate`; relates to `reference_loss_function_method`.

### Mapping bug DIAGNOSED (root cause of the recurring index↔coordinate issue)
- inqkit VTIs are written in **physical order** (`Origin="-25"`,
  `RealField3DWriter` applies `fft_shift_index()` at write — `tddft-simulations
  SKILL.md:790-792`). So **VTIs must NOT be fftshifted**; only LEED screen
  `.dat` files are FFT-natural and need `np.fft.fftshift`.
- My `make_*_postproc.py` GIF scripts call `np.fft.fftshift` on VTI data →
  swaps centre↔edge → user saw the slab at the edges, vacuum in the middle.
  **Visualisation bug, NOT a setup bug.** `check_t1_interior.py` is correct
  (reads origin, no fold → `n_of_z.png` is right). Canonical
  `inqview/pipeline/density.py:60,94` already reads `GetOrigin()` — blessed path.
- **Q1 FIX (locked):** new `inqview.io.load_vti` returning physical-ordered
  array + coordinate axes, with a HARD ASSERT self-check; forbid hand-rolled
  fftshift on VTI; hoist the VTI-vs-screen rule into an ALWAYS-ON path-scoped
  `.claude/rules/` file. ADR to follow once Thread A closes.

### EXECUTION (2026-06-22, after grill)
- **Thread A DONE.** `inqview/visualisation/field_io.py::load_vti` — canonical VTI
  loader, physical order, NO fftshift, returns data+axes, HARD-ASSERTS axis/dim
  invariants + optional `expect_centered_axis` (inner-vs-outer-half |n| mass,
  robust to Friedel; fires on a centre↔edge swap). Lazy export `inqview.load_vti`
  (io stays VTK-free, ADR 0003). Always-on rule
  `.claude/rules/vti-coordinate-mapping.md`. Test
  `inq-stack/tests/python/inqview/visualisation/test_field_io.py` (3/3 PASS) +
  catalogue row. Verified on real GS VTI: centred-slab assert PASSES; a
  deliberately fftshifted copy FAILS (bug caught).
- **Current-data unblock (D3) DONE.** All GIFs re-rendered via `load_vti` (no
  fftshift) using shared helper `hypotheses/_density_views.py::render_total_views`
  → per run: `<prefix>_total/_dfirst/_dprev.gif` (n(t), induced n(t)−n(0), flux
  n(t)−n(t−Δt)). Static run passed `centered_assert=True` (machine-proof the slab
  is centred — answers the cell-21 QUESTION: it was the fftshift viz bug, not the
  setup). Buggy `*_xz_density.gif` deleted. `n_of_z.png` y-axis fixed to 10⁻³
  units (cell-23 TODO). Notebook builders updated (diagnosis cells, 3-view embeds,
  pending-re-run markers, BOLD loss-function Δω≈9 eV caveat + Fourier gate) and
  rebuilt.
- **Only the TOTAL row is available from current data** (runs wrote `density::total`
  only); WP/bath rows, momentum, overlaps, per-component energy, L(q,ω) all marked
  "pending full-suite re-run" in the notebooks — NOT faked.
- **D2 re-run — BUILT, VALIDATED, LAUNCHED (2026-06-22 ~13:57).**
  - Two clean binaries (mirror proven templates + localised background + env CAP):
    `scripts/fullsuite_wp/run.cpp` (gold `run_wp` template; modes via LJ_CAP:
    p3wp no-CAP / p5wp CAP) and `scripts/fullsuite_classical/run.cpp`
    (`cap_baselines` b2 + background). Env: LJ_OUT, LJ_CAP, LJ_N_STEPS,
    LJ_WRITE_EVERY, LJ_WF_EVERY, LJ_LAUNCH_Z, LJ_DT.
  - Built `fullsuite_wp` against inq-study (INQ_SOURCE override + stock-inq share
    pins) — **compiled clean**. 60-step validation run COMPLETED; pipeline
    validator `fullsuite_wp/validate_pipeline.py` PASS on all channels
    (density total/system/wp/gs_system + wavefunction_wp + delta/coarse, momentum,
    wp_momentum/real_space stats, overlap/overlap_full, observables/state_energies/
    occupations) and the `load_vti` total/WP/**bath** decomposition (bath≈234 e⁻,
    WP≈1 e⁻, axes physical). ~3.78 s/step on GPU.
  - **ONE known gap (non-blocking):** GS `eigenvalues.csv`/`occupations.csv` not in
    the `shared_gs/slab_n234_L50` checkpoint (gs_slab never wrote them), so
    `copy_from_checkpoint` no-ops → empty `eigenvalues/` dir. Dynamic
    `state_energies.csv` + `overlap_full` cover excitations; only the GS
    band-structure *reference* plot is unavailable. Retrofit later if wanted
    (needs a tiny GS-load dump; no separate retrofit script exists yet).
  - **PRODUCTION LAUNCHED** via `scripts/fullsuite_dispatch.sh` (background):
    p3wp (GPU0, no CAP, 880) + p5wp (GPU1, CAP, 900) concurrent → then p5cl
    (classical, GPU0, builds inq-study first). Both WP runs confirmed propagating
    (GS load OK). ETA ~2.2 h. Logs: `fullsuite_wp/{p3wp,p5wp}_run.log`,
    `fullsuite_classical/p5cl_build_run.log`, `scripts/fullsuite_dispatch.log`.
- **PRODUCTION COMPLETE (2026-06-22 ~17:33, all exit 0, run_completed=true):**
  p3wp (880 steps, 89 frames total/wp/delta), p5wp (900, 91 frames), p5cl (900,
  91 frames total/delta; wp=0 as expected — classical has no WP orbital;
  ion vz 2.711→2.352, decelerated). Outputs in
  `scripts/fullsuite_{wp,classical}/results/<mode>/raw/{vti,observables}/`.
- **Classical KE logging bug FOUND+FIXED.** `ke_ion_ha` was logged with the *amu*
  mass (1/1822.9) → under-scaled ~1822×. TRAJECTORY (z,vz) is CORRECT (INQ uses
  amu internally). Corrected `electron_track.csv` in post (ke=0.5·v², m_e=1 a.u.)
  → KE 3.675→2.767 Ha, ΔKE_full=0.908 Ha (matches earlier minimal run). Source
  lambda fixed for future runs.
- **D4 analysis = `run-notebook` skill** (user added it). Builds deep per-run
  notebooks (carpets total/Δ0/Δstep/wake + 1 lead GIF + energetics + momentum +
  KS excitation + stopping) over `inqview.pipeline`. **Loss gate RESPECTED**:
  `spectral_weight` is NOT in the runner PHASES, so no L(q,ω) is built; the skill
  already prints a bold low-resolution τ-note. Canonical notebooks (2026-06-23,
  the FIXED set — superseded gapped `p*_run.ipynb` deleted):
  `hypotheses/02_projectile_slab/p3wp_run_notebook.ipynb`,
  `hypotheses/03_cap_stopping/{p5wp,p5cl}_run_notebook.ipynb`.
- **D4 DONE — 3 deep-dive run-notebooks built (0 errors each).** Battery:
  density carpets (total/Δ0/Δstep) + lead GIF; momentum (WP + `_no_wp` bath);
  bath_energy + all-energy breakdown; KS-excitation (occupations + ks_energies
  abs/delta, `_no_wp`); dipole/current spectra; KL/knudsen. p5cl drops WP-orbital
  panels. Builder fixes: run_summary one-key-per-line (regex-reformatted the 3) +
  classical KE rescale. Loss NOT built (gate honoured); τ-note auto-included.
- **2026-06-23 TWO PHASE FIXES (rebuilt all 3):** (1) `overlap` phase looked in
  `raw/overlap/` but runs wrote `raw/observables/overlap/` → symlinked; WP-overlap-
  with-GS-orbital excitation GIFs now land for p3wp/p5wp. (2) `stopping.py` hard-
  requires an `fz` force column the classical track lacked → added physical
  Ehrenfest force `fz = m_e·dvz/dt` (m_e=1 a.u.) to `electron_track.csv`; VALIDATED
  by work-energy theorem `−∫fz dz = 0.908 Ha = ΔKE_full`. Classical headline
  stopping figure (ΔKE→S + force-vs-z) now present. For FUTURE runs: log `fz` in
  the classical `run.cpp` and write overlap to `raw/overlap/` (run-machinery TODO).
- **HEADLINE REPRODUCED from full suite:** p3wp dE=−1.13 mHa (conserved); p5wp WP
  absorbed 0.831, bath intact; **p5cl S = 0.02596 Ha/Bohr = 0.706 eV/Bohr**
  (identical to earlier minimal → reproducible). Still PROVISIONAL (single
  trajectory, dx=0.5).
- **REMAINING polish:** run-SET `*_study.ipynb` still carry stale "pending re-run"
  markers (update to cite full-suite run-notebooks); optional GS eigenvalue
  retrofit; Thread B (skill split + observable infra); loss function gated on
  Fourier-training.
- **Thread B (skill split + observable infra) — deferred to its own focused pass.**

## 2026-06-21 — Design locked, Phase-1 implementation started

### Done
- **Design fully locked** via grill-with-docs (8 decisions). See prompt
  `<locked_decisions>`. Headline: localised positive background = a **static
  `inqkit` perturbation** adding `v_bg = −poisson(n₊)` (no `inq/` or `inq-study/`
  edits); **all phases on `inq-study`**; **built-in sin² CAP** two-sided
  (eta=−0.5, 7.5 bohr/side); **slab N=234 @ r_s≈4 (Na)** in a **50 bohr** box;
  projectile **σ=0.5, 100 eV**; classical twin via existing
  `electron_gaussian_sigma0p35.upf`; **20 feature-aligned screens**; stopping
  power `S = ΔE_bath/x`, x=25 bohr, measured AFTER projectile gone, for BOTH
  classical and WP.
- **All 4 validation tiers (T0+T1+T2+T3) approved.**
- Wrote: prompt, plan, this handover, and CONTEXT glossary section.
- **Core code — T0 VERIFIED (compiles + 3/3 assertions pass on GPU):**
  - `inqkit/jellium/localised_background.hpp` — builds n₊ for slab/sphere/box,
    sharp or erfc edge.
  - `inqkit/jellium/background_perturbation.hpp` — static perturbation: caches
    `φ=poisson(n₊)`, adds `−φ` via explicit `gpu::run` loop.
  - `inqkit/jellium/analytics.hpp` — n₀(r_s), k_F, E_F, E_self helpers.
  - **Bug found+fixed by T0 compile:** `inq::gpu::run` → `gpu::run` (gpu is a
    top-level namespace, not under inq::). Both headers.
- **T0 test PASSED:**
  `inq-stack/tests/include/inqkit/jellium/test_localised_background_engine.cpp`
  (wired into engine CMake). Asserts slab ∫n₊=N, sphere ∫n₊=N, perturbation well
  attractive (∫v_bg·n₊<0). All 3 pass. → catalogue row still TODO.
- **GS LAUNCHED** (background, GPU 0): slab N=234 SCF via
  `systems/localised_jellium/scripts/01_slab_validation/gs_slab/run.cpp` +
  config `shared/configs/slab_n234_L50.hpp`. dx=0.5, LDA, checkpoint →
  `shared_gs/slab_n234_L50`. Engine = stock inq.

- **GS CONVERGED + T1 PASSED.** GS energy −160.99 Ha; external(e–bg) −285 Ha
  (attractive ✓); kinetic 0.067 Ha/e ≈ HEG t(r_s=4)=0.069 ✓. T1 interior:
  ⟨n⟩=100.7% n0, max dev 2.0%, peaks in slab, vac spillout 2.6%. Verdict PASS.
  Profile: `hypotheses/01_slab_validation/n_of_z.png` (+ `check_t1_interior.py`).
  Checkpoint: `shared_gs/slab_n234_L50`.
- **Phase-2 static 2 au run DONE + T3.4 PASS.** 100 steps, 51 frames; energy
  drift 2.2e-8 Ha over 2 au (machine-level) ⇒ background is static/Hermitian.
  Artifacts in `hypotheses/01_slab_validation/`: `static_xz_density.gif`,
  `energy_conservation.png`, `make_static_postproc.py`. Static run at
  `scripts/01_slab_validation/static_2au/`.
- **DELIVERY LIMITATION:** Gmail tool can only create DRAFTS (no send) and
  cannot attach files. Gifs delivered in-session via SendUserFile; a text draft
  was created. Awaiting user decision on the per-run update channel.
- **Phase-2 validation COMPLETE for T0/T1/T3.4.** Remaining in Phase 2: **T2
  (Lang–Kohn)** — surface profile/Friedel, work function, surface energy, grid
  ×½ convergence; needs r_s=4 Lang–Kohn Φ/σ pinned via docs/sources first.
- **DELIVERY = batch at phase ends** (user choice). In-session gifs at phase
  boundaries; no per-run emails.
- **PHASE 3 LAUNCHED** (background, GPU 0, ~2.5 h: ~20 min build + 880 steps ×
  ~9 s). `scripts/02_projectile_slab/wp_slab/run.cpp`: WP σ=0.5 E=100 eV +z from
  z=−23, slab background on, 20 centred plane screens, density frames + dipole +
  current + energy. Engine stock inq.
- **Phase 5 engine note:** inq-study is NOT pre-built; Phase 5 will build it via
  `inq-run` with `INQ_SOURCE=<inq-study>` per run-dir (verify inq-run honours the
  override). GS-checkpoint→inq-study load to be verified then.
- **PHASE 3 DONE + DELIVERED.** 880 steps, t=17.6, 89 frames, 20 screens. Total
  energy conserved (dE=−1.1 mHa, expected no-CAP). xz gif + response delivered
  in-session. Post-proc: `hypotheses/02_projectile_slab/make_wp_postproc.py`,
  `wp_xz_density.gif`, `wp_response.png`.
- **ADR-0008 WRITTEN** (`docs/adr/0008-...`, status accepted/validated).
- **inq-study build path SOLVED:** `INQ_SOURCE=<inq-study>` env override (per
  config.sh); pin INQ_SHARE_PATH/PSEUDOPOD_SHARE_PATH to stock inq install.
- **PHASE 5 WP+CAP LAUNCHED** (background, GPU 0, inq-study). First inq-study
  build (configuring + compiling cleanly at 83%, no errors).
  `scripts/03_cap_stopping/wp_cap/run.cpp`: WP launch z=−15.5, two-sided sin² CAP
  eta=−0.5 width0.15 mid±0.425, perturbation `sum(bg, sum(cap_lo,cap_hi))`. Logs
  total energy + num_electrons (absorbed norm) + dipole. Tests GS→inq-study load.
- **Classical twin PATTERN KNOWN** (ionic::species(UPF).mass; ions.insert;
  velocities()[0]; propagate `.ehrenfest()`). DEFERRED until WP+CAP confirms the
  inq-study path. UPF: `electron_gaussian_sigma0p35.upf`; v_z=+2.71; launch
  z=−15.5; need PROJ species symbol/mass from existing classical config.
- **PHASE 5 WP+CAP DONE + VALIDATED.** inq-study built clean; **GS loaded on
  inq-study (portability OK)**; CAP absorbs WP (235→234.169, 0.83 of WP) with
  **BATH INTACT (234.17 — no over-drain, T3.2 passes)**; energy removed 2.49 Ha.
  900 steps, 91 frames. Post-proc `hypotheses/03_cap_stopping/`:
  `make_wpcap_postproc.py`, `wpcap_traces.png`, `wpcap_xz_density.gif`.
  CAVEAT: WP only 0.83 absorbed at 18 au (0.17 residual) → clean WP bath-ΔE
  stopping needs a longer run; rely on CLASSICAL ΔKE_ion for the headline S.
- **PHASE 5 CLASSICAL+CAP LAUNCHED** (background, GPU 0, inq-study, ~40 min).
  `scripts/03_cap_stopping/classical_cap/run.cpp`: Gaussian-e ion (σ_pot=0.35,
  mass 1 a.u.), v_z=2.71, launch z=−15.5, ehrenfest, same CAP. Logs ion z/vz/KE
  → S=ΔKE_ion/25. Post-proc TODO when it lands.
- **DELIVERY:** hold WP+CAP gif; deliver both Phase-5 gifs + S comparison at the
  phase boundary (after classical).
- **PHASE 5 COMPLETE + DELIVERED.** Classical+CAP: ion decelerated KE 3.675→2.767
  Ha; across slab faces ΔKE=0.649 Ha / 25 Bohr ⇒ **S = 0.0260 Ha/Bohr = 0.706
  eV/Bohr**; bath intact (234→233.78). 91 frames. Post-proc
  `hypotheses/03_cap_stopping/make_classical_postproc.py`, `classical_stopping.png`,
  `classical_xz_density.gif`. Both Phase-5 gifs + S delivered in-session.
  CAVEAT: first-pass, SINGLE trajectory (centroid, +z through centre), COARSE
  dx=0.5 (Phase-3 showed ~1 mHa WP energy drift → grid error at high k); number is
  preliminary, not converged.

## CAMPAIGN STATUS: all 4 run-phases DONE
Phase 1 (impl, T0✓) · Phase 2 (GS, T1✓, T3.4✓) · Phase 3 (WP baseline✓) ·
Phase 5 (CAP stopping: S≈0.71 eV/Bohr classical, WP confirms mechanism✓).
Prompt relocated to `docs/campaigns/localised_jellium/` (campaigns skill).

### REMAINING (not runs)
- **T2 (Lang–Kohn)**: surface profile/Friedel (have n(z)), work function Φ,
  surface energy σ, grid ×½ convergence. Needs r_s=4 Lang–Kohn Φ/σ pinned via a
  docs/sources note + an instrumented/finer GS run. NOT done.
- **Per-phase notebooks** — ✅ DONE 2026-06-22. Executed (0 errors), outputs
  embedded, headline numbers recomputed in-notebook from the provenance CSVs:
  - `hypotheses/01_slab_validation/slab_validation_study.ipynb` (Ph1+2: impl
    sketch, T0/T1, static T3.4)
  - `hypotheses/02_projectile_slab/projectile_slab_study.ipynb` (Ph3 WP baseline)
  - `hypotheses/03_cap_stopping/cap_stopping_study.ipynb` (Ph5 WP+CAP &
    classical+CAP; S=0.7064 eV/Bohr recomputed)
  - **MASTER (all phases):**
    `hypotheses/localised_jellium_campaign_study.ipynb`
  Builders beside each: `build_*_report.py` + shared `hypotheses/_nbreport.py`.
  Rebuild: `PYTHONPATH=…/inq-stack/python venv/bin/python3 build_<name>_report.py`.
- **T0 catalogue rows** in `docs/validation/test-catalogue.md`. NOT done.
- **Convergence**: dx 0.5→0.25 + multi-trajectory/S(v) for a publishable S.

### ENGINE DECISION (refines locked Q2 "uniform inq-study")
inq-study is NOT built (no build dir) and the sin² CAP needs it. So: **GS +
Phase 2/3 (no CAP) on stock inq** (engine-identical for a GS; fast/proven);
**Phase 5 (CAP) on inq-study**, built then, with the GS checkpoint **verified to
load** on it first (portability checked, not assumed). Flagged to user for veto.

### GPU status (2026-06-21)
Both cards BUSY with the user's own jobs: GPU0 `run` PID 2069269 + orted; GPU1
`run` PID 2071679 + orted. `nvidia-smi` shows the cosmetic NVML mismatch (compute
unaffected). **No GPU job may launch until a card frees.** All remaining Phase-1
work is host-only and proceeds regardless.

### Next steps (in order; first three need NO GPU)
1. ✅ DONE — `background_perturbation.hpp` (see Done above). Real-into-complex
   handled by explicit `gpu::run` loop, not `increment` — no type clash to verify.
2. `inqkit/jellium/analytics.hpp` (currently EMPTY) — `e_self_sphere=0.6 N²/R`,
   slab per-area self-energy; plus n₀(r_s), k_F, E_F helpers for tests/reports.
3. `inq-stack/tests/include/inqkit/jellium/test_localised_background.cpp` — T0
   host test: ∫n₊=N (`operations::integral`), interior=n₀, v_bg vs analytic slab
   potential (parabolic inside, linear outside). Build + run on CPU (this is the
   FIRST compile of the two new headers — fix any API slips here). Add rows to
   `docs/validation/test-catalogue.md`.
4. `ResearchProject/systems/localised_jellium/` skeleton (ADR-0007). Slab GS+RT
   `run.cpp` on the jellium `run_template`, passing `sum(background[,CAPs])` to GS
   and RT; `extra_electrons(234)`; write n₊, v_bg, E_self.
5. [GATED on free GPU] Phase-2 GS + T1/T2/T3.4 + 2 au static run → notebook
   `hypotheses/01_slab_validation/` → email.
6. [GATED] Phase-3 WP run. [GATED on 2 GPUs] Phase-5 classical+WP.

### Key API facts (verified by source read)
- `solvers::poisson::solve(field<real_space,double>) → field<real_space,double>`
  (inq/src/solvers/poisson.hpp:217). Drops G=0 (periodic) — exact given ∫n₊=N.
- Build/fill field: `basis::field<Basis,double> f(basis); f.fill(0.0);` then
  `gpu::run(nz,ny,nx,[=]GPU_LAMBDA(iz,iy,ix){... point_op_.rvector_cartesian ...})`.
- Both `ground_state::calculate(...,pert)` and `real_time::propagate(...,pert)`
  take the perturbation → background present in SCF and RT.
- sin² CAP `perturbations::absorbing(amplitude_Ha, mid_pos_frac, width_frac)`;
  two-sided via `perturbations::sum`. mid_pos/width are FRACTIONAL: each CAP
  mid_pos=±21.25/50, width=7.5/50, amplitude=−0.5_Ha.

### Risks / to-pin
- real-into-complex potential increment on inq-study (step 1) — verify early.
- r_s=4 Lang–Kohn Φ/σ values (86.4 erg/cm² is a DIFFERENT r_s) — `docs/sources/`
  note before T2 gates.
- ADR-0008 (perturbation mechanism) still to write.
- Edge: start sharp Θ; soften (erfc) only if Gibbs ringing shows in v_bg.
