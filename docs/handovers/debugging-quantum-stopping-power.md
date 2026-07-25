# Handover: Debugging quantum stopping power — CAP-capture correction

---

## Update: 2026-07-11 (sweep extension) — correction applied to ALL aliasing-valid S(E) points

Status: done. User request: extend the corrected estimator to every point of the
quantum S(E) graph that had no aliasing problem. Valid set (per the user-confirmed
2026-06-27 aliasing finding in `docs/handovers/localised-jellium.md`): v=1.3, 2.0
(p4_wp), 3.0 clean; v=4.0 borderline (~1.1% tail, included + flagged); v=5.0 via
the h=0.35 rerun (clean grid, weak plateau, own GS anchor −71.85697 Ha); v=6.0
EXCLUDED (39% aliased, bound=lower).

New artefacts:
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/build_sweep_notebook.py`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/debugging_quantum_stopping_sweep.ipynb` (11 cells, 0 errors)

Results (E_capt = n_capt × (½v²·27.211 + 81.63) eV; S_corr = (ΔE_plateau − E_capt)/25):

| v | S_orig | E_capt (eV) | S_corr | S_Lind | f explained | verdict |
|---|---|---|---|---|---|---|
| 1.3 | 2.37 | 8.3 | 2.04 | 0.82 | 21% | NOT EXPLAINED |
| 2.0 | 2.39 | 12.4 | 1.89 | 0.44 | 26% | NOT EXPLAINED |
| 3.0 | 2.57 | 12.6 | 2.09 | 0.24 | 21% | NOT EXPLAINED |
| 4.0 | 4.50 | 8.8 | 4.25 | 0.15 | 8% | NOT EXPLAINED |
| 5.0 | 9.78 | 24.0 | 9.09 | 0.10 | 10% | NOT EXPLAINED |

- **5/5 NOT EXPLAINED**; the excess over Lindhard is systematic (1.5–87×), not a
  per-run accounting artifact — S4 books close to ≤0.3% at every v.
- All 5 S1 gates pass vs `se_state.csv` (<1e-3); code KE basis vs measured
  ⟨T_WP⟩(0) agrees to <5% (deviation grows with k0; 3.6% at v=4, aliased tail).
- Every S is an upper bound (unplateaued; v5-h035 worst, late slope −6 eV/a.u.);
  E_capt(t) also still draining at t_f for all runs (correction itself an upper bound).
- Verified en route (user question): E_input = ½v² for the ELECTRON packet —
  23 eV ⇔ v=1.3 confirmed analytically, from ⟨p_z⟩(0)=1.2975, and from the
  pre-slab centroid slope dz/dt=1.33.

Campaign task list extended to 7/7. Next: none — campaign closed; the standing
physics suspect is localisation-energy deposition + beyond-linear-response.

---

## Update: 2026-07-11 (final) — E_capt formula corrected by user to the TOTAL starting KE basis

Status: done (this supersedes the drift-only E_capt in the milestone below).
The user corrected the estimator: the captured 8% takes its share of the packet's
**total** starting KE — the code-inputted drift energy PLUS the σ-derived
localisation energy — and the retained energy is read on the late plateau:

- E_input = ½k₀² = 22.99 eV (`scripts/qsp_phase5/wp/run.cpp:64-65`, LJ_K0=1.3)
- E_loc = 3/(4σ²) = 81.63 eV (σ=0.5, `shared/configs/slab_n82_L50x50x90_E54.hpp:60`)
- basis check: E_input + E_loc = 104.6 eV == run-measured ⟨T_WP⟩(0) ✓ (asserted in-notebook)
- E_capt = 0.0798 × 104.6 = **8.35 eV** (was 1.84 drift-only)
- ΔE_plateau (last-10% mean of E_total − E_GS) = 59.4 eV
- **E_absorbed_jellium = 51.1 eV ⇒ S_corr = 2.04 eV/Bohr** (was 2.30)
- vs S_Lindhard(v=1.3) = 0.82: ratio 1.50, explained fraction f = **21.4%**
- **VERDICT: NOT EXPLAINED** (needs ratio ≤ 0.20). S2 worst case (bath CAP loss
  0.19 e⁻ → up to +20 eV on E_capt) would give S_corr ≈ 1.2, ratio ≈ 0.5 — still
  not explained, though no longer negligible.

Changed: `hypotheses/debugging_quantum_stopping/build_debugging_notebook.py` +
re-executed `debugging_quantum_stopping_v1p3.ipynb` (17 cells, 0 errors);
campaign decisions 4–5 + hypothesis updated. Next: none — campaign closed.

---

## Milestone: 2026-07-11 (later) — SCOPE CORRECTED to `p5_wp_v1p3`; re-executed; verdict NOT EXPLAINED

### Current status
COMPLETE (supersedes the p3 milestone below). The user corrected the target run:
`p5_wp_v1p3` (v=1.3, E_drift=22.99 eV, τ=153.8 a.u., 50×50×90), NOT `p3_wp`. The
campaign prompt, builder, and notebook were redone against v1p3; the superseded
`debugging_quantum_stopping_p3.ipynb` was deleted. Verdict: **NOT EXPLAINED** —
E_capt = 1.84 eV removes only 4.7% of the 38.9 eV gap; S_corr = 2.30 eV/Bohr is
2.8× the Lindhard reference (0.82 eV/Bohr at v=1.3).

### Results (executed `debugging_quantum_stopping_v1p3.ipynb`, 17 cells, 0 errors)
| quantity | value |
|---|---|
| S_WP original (S1 gate: matches `se_state.csv` exactly) | 2.37 eV/Bohr (upper bound; norm remaining 0.080) |
| n_capt = N_total(t_f) − 82 | 0.0798 e⁻ |
| E_capt = n_capt × E_drift (22.99 eV) | 1.84 eV (S3: still draining −0.28 eV over last 10% ⇒ upper bound) |
| S_corr | 2.30 eV/Bohr |
| S_Lindhard(v=1.3) | 0.82 eV/Bohr |
| explained fraction f | 4.7% |
| **verdict** | **NOT EXPLAINED** |

- **S2:** `p4_classical` (same box/CAP, no WP; ran to t=160 a.u.) loses 0.199 e⁻
  of bath to the CAP ⇒ n_capt underestimates the WP remnant by up to ~0.19 e⁻
  (≈ +4.4 eV on E_capt at this drift energy) — verdict robust to it.
- **S4:** books close: injected 109.0 eV (drift 23.0 + zero-point 81.6 + SIE 4.4)
  = CAP-removed 50.0 + retained 59.4 to +0.4 eV (0.3%). No accounting leak.
- **Physics note:** at v=1.3 the drift KE (23 eV) is *smaller* than the packet's
  zero-point energy (~82 eV); most of the retained 59 eV cannot be drift KE at
  all — the zero-point-energy fate is the dominant suspect for the excess S.

### What changed (this milestone)
- `docs/campaigns/debugging_quantum_stopping_power/debugging-quantum-stopping-power.md`
  — rewritten for v1p3 (scope-correction note kept); status done, 6/6.
- `hypotheses/debugging_quantum_stopping/build_debugging_notebook.py` — retargeted
  (direct-CSV routes; the shared ledger module hardcodes T_drift=100 eV, not reusable).
- NEW `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/debugging_quantum_stopping_v1p3.ipynb`
  — the deliverable. `debugging_quantum_stopping_p3.ipynb` DELETED (superseded).
- `docs/campaigns/INDEX.md` regenerated.

### Tests and validation
- S1 reproduce-first PASS (direct CSV = recorded `se_state.csv` 2.374 eV/Bohr to
  <1e-3; tolerance covers the Ha→eV constant precision difference);
  E_jellium(0)−E_GS = +0.36 eV. S2/S3/S4 as tabulated above. 0 execution errors.

### Assumptions still in play
- E_capt estimator (n_capt × drift KE) is user-defined (Inference), generalised
  from the user's 100-eV example to per-run E_drift.
- `p4_classical` bath-drift is an upper-bound proxy for the WP run's bath loss.

### Exact next steps
1. None — campaign closed. Follow-up suspect: zero-point-energy fate (~82 eV vs
   23 eV drift at v=1.3); a σ-sweep or the nazarov_gross comparison would probe it.
2. Optional: commit (docs commit + analysis commit, two-commit hygiene).

---
## (SUPERSEDED) Earlier same-day milestone — p3 execution (wrong run)
---

## Milestone: 2026-07-11 — Campaign planned (grilled + locked), executed, verdict NOT EXPLAINED

### Current status
COMPLETE. The campaign was authored via the campaigns skill (Stage 1–5 grill, all
decisions user-locked), then executed in the same session. The executed notebook
delivers the corrected stopping power and the binary verdict: **NOT EXPLAINED** —
the CAP-capture correction removes only E_capt = 4.4 eV of the 52.5 eV gap
(f = 8.3%), so S_corr = 2.2 eV/Bohr remains ~7.9× the Lindhard reference
(0.28 eV/Bohr at v = 2.711). Campaign frontmatter is `status: done`, 6/6 tasks;
INDEX regenerated.

### Key locked decisions (user, 2026-07-11 grill)
- **Run:** `p3_wp` (τ=100 a.u., 50×50×90, N_bath=82, r_s=5.666, σ_WP=0.5,
  E_drift=100 eV, k₀=2.7110633401, two-sided sin² CAP η=−0.7, faces ±35).
- **E_capt formula (the one free assumption, labelled Inference):**
  E_capt = n_capt × 100 eV with n_capt = N_total(t_f) − 82. Drift KE only —
  the ~81 eV zero-point KE and the measured remnant KE at t_f (2.176 Ha/norm)
  were explicitly REJECTED as alternatives by the user.
- **Reference:** point-charge Lindhard bulk ONLY (`hypotheses/qsp_phase5/lindhard_ref.npz`).
- **Verdict rule:** binary, explained ⟺ |S_corr − S_Lind|/S_Lind ≤ 0.20.
- **Output:** table-only notebook (no new S(v) figure), new folder
  `hypotheses/debugging_quantum_stopping/`.

### Results (all from the executed notebook)
| quantity | value |
|---|---|
| S_WP original (reproduced, S1 gate PASSED) | 2.38 eV/Bohr (upper bound; norm remaining 0.044, late slope −0.086 eV/a.u.) |
| n_capt = N_total(t_f) − 82 | 0.0436 e⁻ |
| E_capt | 4.4 eV (S3: still draining −0.83 eV over last 10% ⇒ upper bound) |
| S_corr = (ΔE − E_capt)/25 | 2.20 eV/Bohr |
| S_Lindhard(v=2.711) | 0.28 eV/Bohr |
| explained fraction f | 8.3% |
| **verdict** | **NOT EXPLAINED** |

- **S2 (assumption-1 check):** the classical twin `p3_classical` loses 0.136 e⁻
  of BATH to the CAP over 100 a.u. — so assumption 1 (CAP only eats the WP) is
  violated at the ~0.14 e⁻ level and n_capt = N_total − 82 systematically
  UNDERESTIMATES the WP remnant. Worst-case shift +13.6 eV on E_capt still gives
  S_corr ≈ 1.7 eV/Bohr ≈ 5× Lindhard ⇒ verdict robust. (Cross-check: the
  WP-projected `norm_check` ratio in `wp_momentum_stats.csv` gives remnant
  ≈ 0.046, consistent with 0.044.)
- **S4 (energy books):** injected 185.2 eV (drift 100 + zero-point 80.8 + SIE 4.4)
  = CAP-removed 126.1 + retained 59.4 to within +0.4 eV (0.2%) — no accounting
  leak; the excess S is physical (zero-point fate / beyond-linear-response),
  not a ledger bug.

### What changed
- `docs/campaigns/debugging_quantum_stopping_power/debugging-quantum-stopping-power.md`:
  rough draft → full autonomy-spec prompt (frontmatter tasks, locked decisions,
  guard rails, preflight); status draft → running → done, 6/6 tasks.
- NEW `ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/`
  (user-chosen location): builder + executed notebook.
- `docs/campaigns/INDEX.md`: regenerated (31 campaigns).
- Side-task in the same session: rebuilt
  `hypotheses/qsp_phase5/p5_wp_v1p3_run_notebook.ipynb` (65 cells, 0 errors) —
  old build linked figures via `../../scripts/...` which viewers refuse to serve;
  current run-notebook builder copies all 41 figures into
  `p5_wp_v1p3_run_notebook_figs/` (0 broken links).

### Files touched
- `/local/data/public/skcb2/tddft/docs/campaigns/debugging_quantum_stopping_power/debugging-quantum-stopping-power.md` — the campaign prompt (done)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/build_debugging_notebook.py` — notebook builder
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/debugging_quantum_stopping_p3.ipynb` — executed deliverable (17 cells, 0 errors)
- `/local/data/public/skcb2/tddft/docs/campaigns/INDEX.md` — regenerated
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/p5_wp_v1p3_run_notebook.ipynb` (+ `_figs/`) — rebuilt (side-task)

### Commands run
```bash
# campaign execution
venv/bin/python3 ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/build_debugging_notebook.py
# index refresh
venv/bin/python3 .claude/skills/campaigns/build_index.py docs/campaigns
# side-task: v1p3 run-notebook rebuild (original sweep args)
venv/bin/python3 .claude/skills/run-notebook/run_notebook_builder.py \
  ResearchProject/systems/localised_jellium/scripts/qsp_phase5/wp/results/p5_wp_v1p3 \
  ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/p5_wp_v1p3_run_notebook.ipynb \
  --run-cpp .../qsp_phase5/wp/run.cpp --cap-inner 35 --rs 5.666 --launch-z -23.75 \
  --v0 1.3 --e-gs-ha -70.22568216820937 --l-slab 25 --lindhard point
```

### Tests and validation
- Proposed & approved (grill): S1 reproduce-first known-case gate, S2 classical-twin
  CAP check, S3 E_capt(t) plateau, S4 energy-books closure.
- Run: all four, inside the executed notebook (0 execution errors).
- Outcomes: S1 PASS (2.377 both routes; E_jellium(0)−E_GS = +0.36 eV);
  S2 quantified 0.136 e⁻ bath loss (assumption-1 violation, verdict robust);
  S3 NOT plateaued (E_capt is itself an upper bound); S4 books close to 0.2%.
- Remaining gaps: none for this campaign's question.

### Trusted sources used
- `docs/plans/quantum-stopping-ledger-26-6-26.md` + `hypotheses/qsp_phase2/quantum_stopping_ledger_p3_26-6-26.ipynb` — retained-energy method, E_GS anchor, full-ledger assumption.
- `qsp_phase3/gs/results/run_summary.txt` — E_GS = −70.22568216820937 Ha (verified directly).
- Phase-5 `lindhard_ref.npz` / `build_se_plot.py` — Lindhard reference curve (reused verbatim).

### Attribution notes
- Notebook base ledger imported from `hypotheses/qsp_phase2/quantum_stopping_ledger.py::compute_wp_ledger('p3')` (existing validated module); only the E_capt arithmetic is new (in-notebook).

### Known issues / blockers
- None. Note only: assumption 1 is measurably violated (S2) — recorded, verdict unaffected.

### Assumptions still in play
- E_capt estimator is user-defined (Inference), not literature-grounded — by design.
- `p3_classical` bath drift is taken as an upper-bound proxy for the WP-run's bath
  loss (the classical projectile excites the bath differently).

### Exact next steps
1. Nothing required — campaign closed. If the S discrepancy is pursued further,
   the notebook's takeaway names the suspects: the σ=0.5 zero-point energy
   (~81 eV) fate and beyond-linear-response deposition (candidate follow-up
   campaign, e.g. against `nazarov_gross_comparison` or a σ-sweep).
2. Optional: commit the campaign + notebook (two-commit hygiene: docs vs sim/analysis).
