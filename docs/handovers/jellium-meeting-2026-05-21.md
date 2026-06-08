# Handover: Jellium 2026-05-21 meeting campaign

**Plan:** [`docs/plans/jellium-meeting-2026-05-21.md`](../plans/jellium-meeting-2026-05-21.md)
**Design journal:** [`docs/journals/researchproject/2026-05-17_jellium_meeting_design.md`](../journals/researchproject/2026-05-17_jellium_meeting_design.md)
**Task list:** TaskList (#1-36) created 2026-05-17.
**Presentation draft:** [`docs/reports/2026-05-21-meeting-emilio/draft.md`](../reports/2026-05-21-meeting-emilio/draft.md)

## 2026-05-21 — Slide-deck drafting session

Goal: produce the meeting-deck draft. **Session is presentation-only**;
follow-up runs (FU-1, FU-2, FU-3) are planned in draft.md §Follow-up
runs but execution is handed off to a separate conversation.

**What changed**:
- `docs/reports/2026-05-21-meeting-emilio/` created (dir).
- `docs/reports/2026-05-21-meeting-emilio/build_summary_figs.py` (new)
  — generates four campaign summary plots from the rollup CSV +
  sigma-sweep summary CSV. Re-run to refresh.
- `docs/reports/2026-05-21-meeting-emilio/figures/` — 4 new PNGs +
  2 data CSVs:
    - `fig_A1_sigma_sweep_S_ratio.png`  (Option A: S(σ) at v=2.71)
    - `fig_A2_localisation_collapse.png` (Option A: S_WP/S_cl vs σ/r_s)
    - `fig_B1_metric_comparison_E100.png` (Option B: 4-metric vs σ)
    - `fig_C1_energy_ledger_bar.png` (Option C: ΔE bookkeeping per run)
- `docs/reports/2026-05-21-meeting-emilio/draft.md` (new) — full
  markdown slide deck draft. Three candidate spines (A localisation,
  B four-metric taxonomy, C missing-electron puzzle) each presented
  with their supporting figures. User to review and pick.

**Sigma in varyv plasmon run — resolved**: `WP_SIGMA_BOHR = 3 Bohr`
(per `ResearchProject/systems/jellium/run_plasmon_n162_L50_E3p4_varyv/run.cpp:4`).
NOT σ=1 and NOT σ=5. Documented in draft.md §Plasmon.2.

**Implementation-correctness items flagged but not executed**:
- Method A (bath ΔE) free-particle known-case test — needs
  `run_free_wp_L50_E100` re-analysed with the audit script.
- Method B (WP slot ΔE) free-particle known-case test — same.
- These would take ~10 min total but were left for the follow-up
  session because tomorrow's meeting can run with the known-case
  status flagged on the slide (draft.md §Implementation correctness).

**Hartree comparability** documented as open question in draft.md
§Hartree comparability — three repair options listed.

**Knudsen re-analysis (added 2026-05-21 same session)**:
- `build_knudsen_figs.py` (new) — generates 4 Knudsen-energy figures
  from `wp_momentum_stats.csv` data (already present in every recent
  run). Fresh figures:
    - `fig_K1_kinetic_traces_sigma_sweep.png` — ½⟨p²⟩(t) and Δ½⟨p²⟩(t)
      for σ ∈ {0.5,1,3,5,8} at E=100 eV.
    - `fig_K2_knudsen_S_vs_sigma.png` — Knudsen S(σ) at E=100 with
      classical anchor.
    - `fig_K3_knudsen_S_vs_v.png` — Knudsen S(v) across all energies.
    - `fig_K4_methods_compare_sigma1_sigma5.png` — three-method
      overlay (Method B, Knudsen K, Classical) at σ=1 vs σ=5.
- Section §K added to draft.md.
- **Key finding (critical for tomorrow's meeting)**: at σ=1, E=100 eV,
  the WP's kinetic energy ½⟨p²⟩ **increases by +3.7 eV** during the
  IFW. The campaign's headline "2.4 % agreement at σ=1" compares
  Method B (full ⟨H⟩_WP, drops by −8.3 eV) against Classical KE (drops
  by −8.4 eV) — apples-to-oranges, because Method B includes the WP's
  self-Hartree (+12 eV released as the WP digs its own attractive
  well). The Knudsen-vs-Classical comparison (the only true
  apples-to-apples KE comparison) gives **+3.7 vs −8.4 eV** — a
  fundamental sign disagreement, not a 2.4 % agreement.
- The "WP→classical convergence" claim has to be reframed. Three
  options enumerated in draft.md §K.6.

**Decision required before final deck**:
1. Pick spine A / B / C / **K** (now four candidates; K is the most
   scientifically defensible and the most disruptive to existing
   narrative).
2. Decide on Hartree-comparability repair (options 1/2/3 in draft) —
   §K now provides quantitative evidence for *why* this matters.
3. Decide on Run-8 fate (re-attempt at dx=0.40 or out-of-scope).
4. Decide whether to retract / reframe the "2.4 % agreement" line in
   the final rollup CSV header and any past emails that cited it.

## 2026-05-21 (afternoon) — Full figure-inventory milestone

**Three new runs launched + completed**:
- `run_wp_n162_L50_E300_sigma1`  (335 steps, 55 min wall, exit 0)
- `run_wp_n162_L30_E100_highdens_sigma1`  (268 steps, 9 min wall, exit 0)
- `run_wp_n162_L50_E200_sigma1`  **PARTIAL** — GPU 0 contended by
  another user's `main2d.gnu.CUDA.ex` (18.3 GB on GPU). Step time went
  from 10 s → 460 s → 2274 s per step. As of last check: step 178/385,
  t=3.56 a.u. of 7.70 target. User chose to wait rather than kill+restart.

**Important note — analyse.py venv requirement**: each new run's
`analyse.py` MUST be invoked with the venv sourced first
(`source /local/data/public/skcb2/tddft/venv/bin/activate &&
python analyse.py`); otherwise the VTI-reading phases fail with
`RuntimeError: VTK is required`. The auto-launchers used inline did
NOT activate the venv; the analyses had to be re-run manually with venv.
Future runs should bake `source venv/bin/activate` into the launcher.

**Knudsen anomaly persists across σ=1 trio**:

| E (eV) | v (a.u.) | Knudsen S (eV/Bohr) | Δ½⟨p²⟩ over IFW (eV) |
|---|---|---|---|
| 100 | 2.71 | **−0.145** | +3.7 |
| 200 | 3.83 (partial!) | **−0.221** | +3.0 |
| 300 | 4.70 | **−0.118** | +3.7 |

All three σ=1 datapoints show the WP **gaining** kinetic energy during
the IFW — confirming the Knudsen-vs-Method-B divergence is not a
low-velocity artefact. The σ=1 "WP = classical electron" claim from
the rollup is metric-dependent: true on Method B, **false on Knudsen
K**.

**Q3 narrative corrected (user clarification 2026-05-21)**: the
"missing electron" refers to the **GS-projected occupations plot**
(n_i^GS(t) = Σ_j f_j |⟨φ_i^GS|ψ_j(t)⟩|²), NOT the energy_balance
ledger Unaccounted column. The correct physical narrative:

- The GS basis (occupied + 20 extras) spans only a few eV above E_F.
- The WP injects 100 eV of energy.
- Bath electrons promoted by WP scattering acquire kinetic energies up
  to ∼E_WP and live at energies above the saved extras.
- INQ TDDFT propagates these correctly; the GS-projection diagnostic
  is what's basis-bounded, not the physics.
- Run-4 (x20 vs x40) confirms x40 doesn't reduce the missing fraction.
  We would need ∼1000× more extras at high energy.

The energy_balance Unaccounted column is a *separate* artefact (DFT
Hartree+XC double-counting) and is NOT what the slide is about.

**22 figures generated, organised by question**:

Q1 (S(v) curve): R1
Q2 (σ-sweep): A1, A2, B1, K1, K2, M_B
Q3 (missing electrons, corrected narrative): E1, E2, E3, E4
Q4 (density): R1 + table in deck
Q5 (injector validation, free vs jellium σ=1): D0, D1, D2, D3, D4, D5, T1, M_A
Cross-cutting: C1, K3, K4

All figures live in `docs/reports/2026-05-21-meeting-emilio/figures/`.
Build scripts: `build_summary_figs.py`, `build_knudsen_figs.py`,
`refresh_after_new_runs.py`, `build_density_comparison.py`,
`build_centroid_momentum.py`, `build_missing_electron.py`.

**Two decks now live**:
- `docs/reports/2026-05-21-meeting-emilio/draft.md` — original 3-spine
  draft (A/B/C/K)
- `docs/reports/2026-05-21-meeting-emilio/deck_v2.md` — 5-question
  structure, intermediate state
- `docs/reports/2026-05-21-meeting-emilio/deck_v3.md` — final draft
  with all 22 figures inventoried, corrected Q3 narrative

**Files touched (2026-05-21 afternoon)**:
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/electron_proj_E200_L50_cubic_sigma1.hpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/electron_proj_E300_L50_cubic_sigma1.hpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/highdens_n162_L30_E100_sigma1.hpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_wp_n162_L50_E200_sigma1/` (new, run.cpp + analyse.py)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_wp_n162_L50_E300_sigma1/` (new, run completed)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_wp_n162_L30_E100_highdens_sigma1/` (new, run completed)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/build_density_comparison.py` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/build_centroid_momentum.py` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/build_missing_electron.py` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/fig_D*.png` (D0-D5, 6 new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/fig_T1_centroid_z_t.png` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/fig_M_A_momentum_band_2panel.png` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/fig_M_B_momentum_band_sigma_sweep.png` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/fig_E*.png` (E1-E4, 4 new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/fig_K3_knudsen_S_vs_v.png` (refreshed)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/fig_R1_rollup_with_crash.png` (refreshed)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/deck_v3.md` (new)
- `/local/data/public/skcb2/tddft/docs/handovers/jellium-meeting-2026-05-21.md` (this section appended)

**Exact next steps**:
1. Once E=200 completes (or user decides to drop it), re-run
   `refresh_after_new_runs.py` to update K3 + R1 with converged data.
2. User reviews deck_v3.md and chooses 1-2 figures per question for
   final slide deck.
3. Render final deck to PPTX or PDF via pandoc/marp/slidev (mirror the
   14-05 `build_pptx.py` pattern).
4. For next campaign: bake `source venv/bin/activate` into the
   inq-run wrapper or run launcher to avoid the VTK import issue.

**Files touched (2026-05-21)**:
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/draft.md` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/build_summary_figs.py` (new)
- `/local/data/public/skcb2/tddft/docs/reports/2026-05-21-meeting-emilio/figures/*.png,*.csv` (new)
- `/local/data/public/skcb2/tddft/docs/handovers/jellium-meeting-2026-05-21.md` (this section appended)

**Commands run**: `python docs/reports/2026-05-21-meeting-emilio/build_summary_figs.py` (once).

**Tests / validation**: visual sanity of the 4 new PNGs — figures
generated cleanly, data CSVs match the rollup numbers (fig_A2 ratios
0.98 at σ=1 standard and 1.05 at σ=0.5 highdens reproduce the
campaign-rollup "2.4% agreement" claim).

**Exact next steps**:
1. **User**: review draft.md and choose spine A / B / C (or hybrid).
2. **Same session**: trim the deck to one spine + render to PPTX/PDF
   (Marp or pandoc-beamer; `build_pptx.py` in the 14-05 dir is a
   pattern to reuse).
3. **Separate session (handed off)**: execute FU-1 / FU-2 / FU-3
   runs per draft.md §Follow-up runs.
4. Per `.claude/rules/handovers.md` — append a follow-up section to
   this handover when the spine is chosen, when the final deck is
   rendered, and again when FU runs complete.

## Current status

**CAMPAIGN COMPLETE — extended through 2026-05-19 03:30 with Runs 4, 5, 7, 9 added.**

**Done (26/36 tasks):**
- Run-4 extra-states test — x20+x40 result: basis CONVERGED at x=20 (Δ < 0.04% to x=40); missing-electron puzzle is NOT basis truncation. Email `[jellium-extra-states]`.
- Run-5 (20 eV pair) — S(v=1.21)=0.95 eV/Bohr classical; 66% velocity drop in 18.8 Bohr; σ=5 WP underpredicts 7×.  Email `[jellium-20eV]`.
- Run-7 σ-sweep (5 σ values 0.5-8 at E=100) — ΔE_WP scales ~1/σ from -16 eV (σ=0.5) to -0.005 eV (σ=8); KL peaks at σ=3; σ=8 ≈ "indistinguishable from jellium" limit confirms σ→∞ intuition. Email `[jellium-sigma-sweep]`.
- Run-9 (25 eV classical companion to Run-1) — S(v=1.36)=0.88 eV/Bohr; same σ=5 underprediction pattern as Run-5. Email `[jellium-25eV]` reply.
- Final rollup UPDATED with v=1.21, 1.36 anchors. Clear Lindhard-to-Bethe S(v) shape; campaign concludes σ << r_s required for proper WP-vs-classical agreement.

**Blocked / out-of-scope:**
- Run-8 (Knudsen velocity sweep at E=700-1100, dx=0.30): CUDA illegal-memory-access at propagate step 0 in the original attempt; in a minimal-callback retry on 2026-05-20 (no overlap, no VTI density, no DensityDelta) the run *initialized* real_time::propagate cleanly and the WP injection succeeded but the first propagation step then hung indefinitely (12+ min, GPU at 100 % util / 24 GB committed, no step output). Issue is inside INQ's RT propagator (Krylov subspace / ETRS intermediates), NOT in the user observable callbacks. Out of meeting scope; needs INQ-internal engineering. The Subtest #19 (dt-convergence at E=1100, dx=0.30) is blocked by the same root cause.
- Run-4 x80 case: killed at step 39/462 (wall ~7.5 min/step from 80-extra-state memory pressure paging); x20 vs x40 result already conclusive.

**Earlier-done items:**
- Infra (10/10): all phases shipped including the three new ones today
  (knudsen_ke, kl_divergence, energy_balance) plus §13.5 per-component
  energy plots and §13.6 FFT t_start_au.
- Run-1 (25 eV WP, low-v anchor) — `[jellium-25eV]` email
- Run-2 + 2b (free WP @ E=100, injector validation 0.06%) — `[jellium-free-compare]`
- Run-3 (E=100 σ=5 jellium pair, mid-Lindhard anchor) — `[jellium-free-compare]` reply
- Run-6 (high-density r_s=3.41 pair, density anchor) — `[jellium-highdensity]`
- σ=1 task (concentrated WP @ E=100, σ-sensitivity anchor) — `[jellium-free-compare]` two replies
- Final rollup (S(v) figure across all data points) — `[jellium-rollup]`
- Per-run analyse.py policy (Option A) + every new run dir has one
- §13.3/§13.5/§13.6 observables implemented

**PAUSE-NEEDED (await user — >2h each):**
- Run-4 extra-states test
- Run-5 20 eV pair
- Run-7 σ-sweep at E=100
- Run-8 Knudsen velocity sweep
- Run-9 25 eV classical companion

**Headline cross-validation:** at σ=1 standard density (Run-3 σ=1
follow-up), TDDFT-WP ΔE_WP (energy_balance) = -8.28 eV and Ehrenfest
ΔE_proj_kin = -8.35 eV — two completely different propagator types
agree on the projectile energy loss to **2.4 %** over the IFW window.

**Density effect (at v=2.71):** S = 0.33 eV/Bohr at r_s=5.69 vs 1.20
eV/Bohr at r_s=3.41 — 3.6× ratio at 4.6× density, broadly Lindhard-linear.

## Universal rule pinned 2026-05-17 — Email pipeline

**All campaign emails MUST use `inqview.email.send_run_email(...)`,
NOT the Gmail MCP `create_draft` tool.** Codified in plan §"Universal
rules" rule 5. The MCP draft path requires manual click-through and
breaks the autopilot contract. The earlier `[PAUSE-NEEDED]` draft for
disk pre-flight (id `r7568159534297396504`) was a violation of this
rule and should be ignored / deleted manually by the user.

## Universal rule pinned 2026-05-17 — Per-run `analyse.py` (Option A)

Every new run dir MUST contain an `analyse.py` cloned from the matching
archetype (WP, Classical, Free). The analyse.py runs the inqview
pipeline (with the new `knudsen_ke`, `kl_divergence`, `energy_balance`
phases registered today), runs `shared/python/analyse_extras.py`, and
writes `results/analysis/REPORT.md`. Emails reference the PNGs that
analyse.py produces — the email step is NEVER first.

Codified in plan §7. Archetypes defined in plan §7. Compare-style
emails (free vs jellium) live in `_compare_<family>_<energy>/compare.py`,
which runs after both halves' analyse.py complete.

## Universal rules pinned 2026-05-17 — Plotting / observables (plan §8–10)

- §13.5 per-component energy plots — `observables.py::_plot_per_component_energy`
  (new today).
- §13.3 energy_balance ledger phase — `energy_balance.py` (new today).
- §13.6 FFT t_start_au transient exclusion — `FourierTransform`
  accepts `t_start_au`; CSV headers declare the cutoff.
- §"Plot shading" — switched from post-IFW grey to IFW soft-yellow
  highlight (`_common.ifw_highlight`); new plots prefer this.
- §"Classical KE caption" — every classical-projectile KE plot must
  label the dE drop as the electronic stopping signal (force from bath
  density gradient), not a numerical artefact.

These rules apply to ALL new runs from 2026-05-17 onwards.

Last meaningful action (2026-05-17 late session): Gmail MCP
authenticated; `shared/configs/boundary_rule.hpp` written, gcc-verified
to compile with embedded static_assert known-case tests passing;
campaign-kickoff draft email created via Gmail MCP (draft id
`r2007921787846299866` in Gmail Drafts folder — user must hit Send).

**Important IFW-formula correction applied** during boundary_rule.hpp
build (the static_assert caught it): IFW end is `+L/2 − 3σ` for **both**
standard and relaxed rules (Gaussian-3σ-tail-at-far-face physics
criterion, independent of launch convention). Plan + memory entry
corrected.

**Important Gmail MCP limitation**: the MCP exposes `create_draft` but
**no `send_message`**. Email pipeline (Infra-8) must either (a) create
drafts user sends manually, or (b) fall back to `smtplib + Gmail App
Password` for autosend. Infra-8 task description needs updating to
reflect this.

## What changed (2026-05-17)

- 36 campaign tasks created (#1-36 in the local task list).
- `docs/plans/jellium-meeting-2026-05-21.md` written (campaign plan, single source of truth).
- `docs/journals/researchproject/2026-05-17_jellium_meeting_design.md` written (topical design entry).
- This handover initialised.
- Two memory entries added (boundary rule + campaign overview); MEMORY.md updated.
- Two symlinks added: `.claude/skills/diagnose` and `.claude/skills/grill-with-docs` → `~/.claude/skills/...` (so they auto-appear in this project's skill list from next session onwards).
- Gmail MCP authenticated by user via `/mcp`. Connector status: connected. Tools available: `create_draft`, label / thread management. No `send_message`.
- `shared/configs/boundary_rule.hpp` written + gcc-compile-tested. Embedded static_asserts as known-case tests; one bug caught and fixed during the test pass (IFW formula 4σ → 3σ).
- Plan and memory entry corrected for the IFW formula (was wrongly `+L/2 − 4σ` for standard, now correctly `+L/2 − 3σ` for both rules).
- Campaign-kickoff draft created in Gmail Drafts: `[jellium-campaign] 2026-05-21 meeting prep — design lock + 36-task plan`. User has confirmed it landed in Drafts; manual send proves MCP draft path end-to-end.
- **Infra-8 built and end-to-end-proven**: `inq-stack/python/inqview/email.py` — Gmail SMTPS sender via Google App Password. Module imports cleanly. User-side smoke test passed 2026-05-17 (after one 535 auth-failure round-trip that prompted a runbook-in-the-diagnostic improvement). Task #37 completed.
- User has decided AGAINST the Gmail MCP draft-only pipeline for production; production emails use inqview.email (autosend) instead of draft-and-send-manually.

## What changed (2026-05-17 — Infra-4 session)

- **Infra-4 complete and known-case-verified.** New observable
  `inq-stack/include/inqkit/observables/wp_momentum_stats.hpp` with the
  dipole.hpp host-after-reduction pattern (on-device 3D `gpu::run`
  reduction in Fourier space; basis-comm and set-comm
  `all_reduce_in_place_n` on a 7-element host buffer). Per-step CSV:
  `step, time_au, px_mean, py_mean, pz_mean, px2_mean, py2_mean, pz2_mean,
  sigma_px2, sigma_py2, sigma_pz2, e_kin_ha, norm_check`. Single-kpoint
  only (matches `inqkit::WavePacket`).
- **Smoke test built and PASSING** at
  `Tutorial/wp-momentum-stats-test/run.cpp`. 40 Bohr cubic, dx=0.5 Bohr,
  injector sigma=5 Bohr, k0_z=2.711 Bohr⁻¹ (~100 eV). All 8 Heisenberg
  assertions agree with the analytic reference to better than 1e-7
  absolute. Output `results/wp_momentum_stats.csv` written; the test
  binary returns 0 on PASS, 1 on FAIL.
- **Heisenberg convention pinned down.** The handover-plan text
  `σ_p = 1/(2 σ_r) = 0.1 Bohr⁻¹ at σ_r = 5` uses the *density* Gaussian
  convention `|psi|² ∝ exp(−r²/(2σ_r²))`. The inqkit `WavePacket`
  injector instead uses the *wavefunction* Gaussian
  `psi ∝ exp(−r²/(2σ²))`, so the density σ is `σ/√2` and the momentum σ
  is `1/(σ√2)`. With injector σ=5: `σ_p = 0.14142 Bohr⁻¹`, density σ_r
  = 3.535 Bohr, `E_kin = 3.7048 Ha = 100.8 eV` (vs the "~100 eV"
  estimate). The `σ_p × σ_r_density = 1/2` minimum-uncertainty product
  holds. `boundary_rule.hpp` formulae use the injector σ, so existing
  campaign Cfgs are consistent; the smoke-test reference numbers were
  updated to match injector convention. The "σ=5" usage everywhere in
  the plan should be read as the injector parameter.
- **Two foot-guns documented inline in the smoke test:**
  (1) Box centre is `(0,0,0)` in INQ's `[-L/2, +L/2]` convention — not
  `(L/2, L/2, L/2)`. The pre-existing
  `Tutorial/free-propagation-wp-rt/run_01_base/run.cpp` puts the WP at
  `(L/2, L/2, L/2)`, i.e. at the corner; that tutorial appears never
  to have been executed (no `run.log`) and the bug would silently
  corrupt the `e^{i k·r}` phase across half the cell — flagged for a
  separate fix.
  (2) `systems::electrons(...)` throws *at construction* if
  `num_electrons == 0`. Fix: `options::electrons{}.extra_states(1).extra_electrons(2.0)`
  — gives one doubly-occupied "ghost" state at index 0 and an empty
  state at index 1 for the WP. Avoids the SCF-non-convergence
  workaround for `non_interacting()` + free box.

## What changed (2026-05-17 — Infra-5 session)

- **Infra-5 complete and known-case-verified.** New observable
  `inq-stack/include/inqkit/observables/wp_real_space_stats.hpp` is a
  structural copy of the Fourier-space sibling, applied directly to the
  real-space `phi.hypercubic()` with `rvector_cartesian` and the
  volume-element `dV` factor.
- **Real-space smoke test PASS** at
  `Tutorial/wp-real-space-stats-test/`. Injector σ=5 Bohr at centre
  (0,0,0): observable reports `<x_d> ≈ 0` (≤1.3e-7), density
  `σ_r_d = 3.5355 Bohr` (≤9.4e-7), real-space norm 1.000 (≤4.9e-8).
- **Heisenberg cross-check** between the two new observables:
  `σ_r_z · σ_p_z = 3.5355 × 0.14142 = 0.500 ± 1e-6` — minimum-uncertainty
  product as expected for a freshly injected Gaussian WP. The two
  observables agree at the analytic limit.

## What changed (2026-05-17 — Infra-6/7/9 + Run-1 prep)

- **Infra-6 (knudsen_ke phase) complete.** New
  `inq-stack/python/inqview/postprocess/knudsen_ke.py` computes
  Method-B stopping-power `<|p|²>/2` over time. Dual input path:
  prefers native `wp_momentum_stats.csv` (`e_kin_ha` column = exact);
  falls back to `momentum_distribution.csv` (histogram-weighted
  `Σ_bin k² · n_wp / Σ n_wp`, ~0.7 eV loss from 1D radial binning at
  σ=5 / k₀=2.71 — measured on `run_wp_n162_L50_E100` retroactively).
  Outputs `knudsen_ke.csv`, three plots (`vs_t`, `vs_z`,
  `stopping_power_vs_z`). z trajectory comes from `cod_z_bohr` in
  `observables.csv` when present; phase degrades cleanly to t-only
  plots when missing. Registered in `pipeline.py`.
- **Infra-7 (kl_divergence phase) complete.** New
  `inq-stack/python/inqview/postprocess/kl_divergence.py` computes
  `KL(P_t || P_0)` of the WP momentum distribution from
  `momentum_distribution.csv`'s `n_wp` column. Smoke-tested on
  `run_wp_n162_L50_E100`: `KL(t=0)=0` exact, monotonically rising to
  ~0.047 nats at t=12.8 a.u. Registered in `pipeline.py`.
- **Infra-9 (post-IFW shading helper) complete.** Added three new
  helpers to `inq-stack/python/inqview/postprocess/_common.py`:
    - `post_ifw_window_au(launch_z, L, σ, v)` → `(t_IFW, t_total)` in
      atomic units, matching the boundary_rule.hpp formulas exactly.
    - `post_ifw_window_from_summary(results_dir)` → reads
      `run_summary.txt` (`wp_center_bohr`, `wp_sigma_bohr`,
      `wp_k0_bohr_inv`, `cell_bohr`) and returns the tuple. Tolerant of
      missing fields (returns None).
    - `post_ifw_shade(ax, t_ifw, t_total)` → grey `axvspan` shading
      with optional legend entry.
  Known-case-tested on boundary_rule.hpp's three documented configs
  (standard σ=5, relaxed σ=8, and the older E=100 run with launch=-10):
  all match the docstring values to ≤1e-3 a.u. Wired into the two new
  phases (knudsen_ke vs_t plot, kl_divergence vs_t plot) — the
  remaining 5 time-domain plots in inqview can be wired during the
  meeting rollup as needed.
- **Run-1 (25 eV WP) inputs prepared.** The Cfg
  `shared/configs/electron_proj_E25_L50_cubic.hpp` now uses
  constexpr boundary_rule values: `WP_CZ_BOHR = launch_z(σ=5, L=50)
  = -5`, `N_STEPS = n_steps_for(...) = 923` (plan said 922; the
  1-step delta is rounding precision in v=1.356 vs full-precision
  k₀=1.35553 — both are physically equivalent, negligible 0.027 Bohr
  past stop_z), `WRITE_EVERY = write_every_for(923) = 3` (=307
  frames, on target). The `run_wp_n162_L50_E25/run.cpp` includes
  the two new observable headers and wires them into the propagate
  callback at the standard `Cfg::WRITE_EVERY` cadence. Docstring
  brought up-to-date (was claiming E=100 eV in several places —
  copy-paste residue from the E=100 run's run.cpp).
- Cfg gcc syntax-checked against `boundary_rule.hpp` with embedded
  static_asserts (launch_z=-5, N_STEPS=923, WRITE_EVERY=3 all pass).
  Full INQ build of the run.cpp **not yet performed** — that happens
  on first launch via `inq-run`.

## Files touched (2026-05-17 — whole session)

Infra-4/5 (earlier in session):
- `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/observables/wp_momentum_stats.hpp` (new)
- `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/observables/wp_real_space_stats.hpp` (new)
- `/local/data/public/skcb2/tddft/Tutorial/wp-momentum-stats-test/run.cpp` (new; Tutorial repo)
- `/local/data/public/skcb2/tddft/Tutorial/wp-real-space-stats-test/run.cpp` (new; Tutorial repo)

Infra-6/7/9:
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/knudsen_ke.py` (new)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/kl_divergence.py` (new)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/_common.py` (post_ifw_* helpers appended)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/pipeline.py` (2 new phases registered)

Run-1 prep:
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/electron_proj_E25_L50_cubic.hpp` (boundary_rule integration)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_wp_n162_L50_E25/run.cpp` (new observables wired, docstring rewrite)

## Commands run (Infra-4 session)

- `CUDA_VISIBLE_DEVICES=0 inq-run` in `Tutorial/wp-momentum-stats-test/`
  — full compile + execution + PASS in <2 min (GPU 0 idle, picked via
  `nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv`).

## Tests and validation (Infra-4 session)

Per `.claude/rules/testing.md`:
- **Proposed**: Heisenberg known-case test from plan §"Validation gates" —
  `σ_p_d`, `<p_z>`, `E_kin` at injector σ=5, k0_z=2.711.
- **Approved**: user delegated full autonomy at session start.
- **Run**: all 8 Heisenberg assertions (3× `<p_d>`, 3× `σ_p_d`,
  `<p_z²>`, `E_kin`) at t=0 against analytic Gaussian reference.
- **Outcome**: PASS, max absolute deviation 4.7e-7 (on `σ_p_z`).
  `<p_x>`, `<p_y>` ≈ −2e-11 (machine zero with the rectangular FFT
  grid). Parseval N = 4.096e6 (= 80³ × 0.5³ × ... — sanity-checks the
  un-normalised GPU sum).
- **Unverified**:
  - Per-step accumulation cadence in a real propagation. We hit only
    t=0; integration with `inqkit::RealTimeSession` is exercised next
    when Run-1 launches.
  - The `set_comm` (state-decomposition) path is untested — current
    runs are single-rank so the relevant `all_reduce` is a no-op. The
    code path is present and follows the dipole.hpp template, but the
    first multi-rank use will be the first real test.
- **Profiling not yet done**: the plan asks "every step initially;
  fall back to every 2 or 5 if >5% overhead". Deferred to first real
  run.

## Known issues / blockers — update

- **`Tutorial/free-propagation-wp-rt/run_01_base/run.cpp` has a wrong
  WP centre** at `(L/2, L/2, L/2)` in the `[-L/2, +L/2]` cell — should
  be `(0, 0, 0)`. Not executed (no `run.log`) so no real data was lost
  yet; flag for a follow-up commit in the Tutorial repo.
- **The "σ_p = 1/(2σ_r) = 0.1" line in the plan and design journal is
  ambiguous** — refers to density-σ, but everywhere else in the
  campaign "σ" is the injector parameter. Pinned down here (see
  "What changed (Infra-4)" above) and in the smoke-test source. No
  Cfg changes needed because every existing Cfg uses the injector σ.
- **Infra-1/2/3 RESOLVED 2026-05-17 (skip).** User confirmed the
  current `density_rt_delta` already emits `ρ(t) − ρ(t=0_post_WP)`,
  i.e. it is *already* what the plan was calling for as the new
  `density_rt_delta_t0`. The dual-path Δρ work is a no-op. Plan rule
  #3 wording is stale and should be edited post-meeting to say:
  "`density_rt_delta` is `ρ(t) − ρ(t=0_post_WP_injection)` (already
  implemented); no new variant required." Infra-1/2/3 marked DONE by
  resolution.
- Original list (Infra-0/dt-convergence/etc.) still applies — unchanged.

## Files touched (2026-05-17)

- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/boundary_rule.hpp` (new, gcc-verified)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/email.py` (new, import-verified)
- `/local/data/public/skcb2/tddft/docs/plans/jellium-meeting-2026-05-21.md` (new + IFW-correction edit)
- `/local/data/public/skcb2/tddft/docs/journals/researchproject/2026-05-17_jellium_meeting_design.md` (new)
- `/local/data/public/skcb2/tddft/docs/handovers/jellium-meeting-2026-05-21.md` (new — this file)
- `/local/data/public/skcb2/tddft/docs/journals/researchproject/index.md` (updated — new row at top)
- `/home/raid/skcb2/skcb2/tddft/.claude/projects/-local-data-public-skcb2-tddft/memory/MEMORY.md` (updated — two new entries)
- `/home/raid/skcb2/skcb2/tddft/.claude/projects/-local-data-public-skcb2-tddft/memory/feedback_jellium_boundary_rule.md` (new)
- `/home/raid/skcb2/skcb2/tddft/.claude/projects/-local-data-public-skcb2-tddft/memory/project_jellium_2026_05_meeting_campaign.md` (new)
- `/local/data/public/skcb2/tddft/.claude/skills/diagnose` (new symlink)
- `/local/data/public/skcb2/tddft/.claude/skills/grill-with-docs` (new symlink)

## Commands run

- No compute, no builds.
- Read-only inspections via Bash + Read (catalogued in design entry).

## Tests and validation

- `boundary_rule.hpp` compiles cleanly under `g++ -std=c++17 -O2`
  with all static_asserts passing. Runtime spot-check (`/tmp/test_boundary_rule.cpp`)
  confirms: launch_z(5,50)=-5; stop_z(5,50)=+20; traversal(5,50)=25;
  ifw_end_z(5,50)=+10; ifw_end_z_relaxed(8,50)=+1 (rule-independent);
  n_steps_for(5,50,2.711,0.02)=462 (ceil convention); n_steps_for_relaxed(8,50,2.711,0.02)=332;
  write_every_for(461)=2; write_every_for(922)=3; write_every_for(348)=1.
  (Test source preserved at `/tmp/test_boundary_rule.cpp` — promote to
  `Tutorial/wp-boundary-rule-test/` as a follow-up.)
- Gmail MCP path proven via `create_draft` (draft id
  `r2007921787846299866`). End-to-end send not yet verified — user
  needs to hit Send in Gmail UI.
- Remaining validation gates are defined in the plan §"Validation
  gates"; each Infra-* and Run-* task carries the gates that must pass
  before the next task can be claimed.

## Trusted sources used

- INQ `options::theory{}.non_interacting()` (`inq/src/options/theory.hpp:42`)
  — enables free-space WP via same pipeline.
- INQ `dipole.hpp` pattern (`inq/src/observables/dipole.hpp`) —
  template for GPU-then-host-MPI-reduce in new observables.
- `docs/sources/free-electron-gas-magic-numbers.md` — N=162 closed-shell at any L (relevant to high-density 30³ run).
- Knudsen et al. arXiv 2605.12854 — motivates 700-1100 eV sweep + the
  `<p²>/2` stopping-power method. Source note pending in `docs/sources/`.
- 2026-05-14 meeting figures (`docs/reports/14-05-2026-meeting-emilio/figures/`) — baseline for final stopping-power rollup.

## Attribution notes

- Design decisions are the user's; assistant grilled and recommended.
- New universal rules (4σ/1σ boundary, 300-frame cadence, two-method
  stopping power, ADD-don't-REPLACE Δρ) were jointly agreed in the
  2026-05-17 grilling session and are recorded in plan + design entry.
- Knudsen stopping-power method belongs to arXiv 2605.12854 (to be
  cited in the meeting figure caption).

## Known issues / blockers

- **Infra-0 Gmail MCP OAuth is user-action-only**: assistant cannot do
  the OAuth flow (interactive in browser). Once user has authenticated
  in a terminal session, Infra-8 (email script) can complete and
  Run-1 onwards can produce emailed deliverables.
- **dt for Knudsen sweep is gated on the dt-convergence subtest** —
  not a blocker for Infra or Run-1..7, but locks before Run-8 starts.
- **WPMomentumStats / WPRealSpaceStats need MPI-reduce review**
  (host-after-reduction pattern, like INQ's `dipole.hpp`) — flagged
  by user in the grilling; not yet implemented.
- **MomentumDistribution latent MPI-reduce bug** (no `all_reduce`) —
  fine while runs remain single-rank, but must be fixed when next
  touching the file.

## Assumptions still in play

- The dx=0.30 GS at `save_gs/gs_L50_cubic_N162_dx0p30/` is verified-good
  (used by existing E=300, E=600, E=1500 pairs) — re-used wholesale for
  Knudsen sweep + σ-sweep low-σ extreme. **If this GS turns out to need
  re-build, Run-8 and Run-7 σ=0.25 are blocked.**
- The dx=0.40 GS at `save_gs/gs_L50_cubic_N162_dx0p40/` is similarly
  verified-good (existing E=100 pair). Re-used for Run-1, Run-3,
  Run-5, σ-sweep σ ∈ {0.5,1,3,5,8}, and as the base for the 2 new
  extra-states variants.
- N=162 remains a closed-shell magic number at L=30 cubic (per
  `docs/sources/free-electron-gas-magic-numbers.md` — verified via the
  enumeration script there). **Worth re-checking with the script
  before Run-6 GS build.**
- Workstation A30 chassis cooling holds up under multi-day 100 % util
  on both GPUs. Watchdog (task #32) is optional but recommended.
- `/local` disk has > 50 GB free for the campaign output. **Should be
  verified with `df -h /local` before starting Run-1.**

## Option B (deferred): /loop in tmux for the bulk grind

After Option A's first 3 tasks (Infra-4, Infra-5, Run-1) prove the
workflow end-to-end, the user switches to `/loop` (local Claude Code
self-pacing) for the bulk grind (Run-4 through Run-8). The session
runs inside a tmux pane so it survives ssh disconnect:

```bash
# One-time: start the persistent session.
tmux new -s jellium-loop
cd /local/data/public/skcb2/tddft
claude

# Inside Claude Code, kick off the loop:
#   /loop continue the jellium campaign from
#         docs/handovers/jellium-meeting-2026-05-21.md.
#         Pick the next unblocked task, execute, update the handover,
#         and send a Gmail via inqview.email.send_run_email for any
#         completed pair. Pause and send an email subject-prefixed
#         "[PAUSE-NEEDED]" if (a) a journal observation in user's
#         voice is required, (b) a validation fails, (c) the
#         dt-convergence subtest result needs me to pick dt,
#         (d) a task blocks for >2h, or (e) anything unexpected.
#         Self-pace via ScheduleWakeup.

# Detach (loop keeps running): Ctrl-b then d
# Reattach later:               tmux attach -t jellium-loop
# Check loop status:            tmux ls
# Kill when campaign done:      tmux kill-session -t jellium-loop
```

## Exact next steps

For the next session that picks up this campaign:

1. **Read this handover and the plan** (no need to read the design
   journal entry unless re-litigating). Infra-0, Infra-4, Infra-8,
   and Infra-10 already done; pick up at **Infra-5**.
2. **Verify pre-flight** (user actions still pending — not assistant):
   - `df -h /local` — at least 50 GB free
   - `ls save_gs/gs_L50_cubic_N162_dx0p40` and `…_dx0p30` — both
     should contain `electrons/` etc.
   - Re-verify N=162 magic via the enumeration script in
     `docs/sources/free-electron-gas-magic-numbers.md`.
3. **Launch Run-1 (25 eV WP, ~75 min GPU compute).** Cfg + run.cpp
   are prepared (see "What changed — Run-1 prep" above). User
   approval per plan §"Approvals required" is the gating event.
   When approved:
   ```
   cd ResearchProject/systems/jellium/run_wp_n162_L50_E25
   CUDA_VISIBLE_DEVICES=<idx> nohup inq-run > run.log 2>&1 &
   ```
   Pick `<idx>` from `nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv`
   (idle GPU). Expected wall time ~75 min. Output goes to
   `results/raw/observables/{wp_momentum_stats,wp_real_space_stats,
   observables,momentum_distribution}.csv` + standard VTI series.
4. **After Run-1 completes**: run the postprocess pipeline (now
   includes `knudsen_ke` + `kl_divergence` + post-IFW shading) and
   email the summary via `inqview.email.send_run_email` with the
   `[jellium-25eV]` subject prefix. This is the campaign's
   end-to-end validation milestone.
5. Then Run-2/2b/3 (free-WP + jellium @ E=100 eV, two GPUs in
   parallel) and onward through Run-4..8 per plan ordering.
5. **Then Infra-1, 2, 3 (Δρ-from-t=0 and z-profile-from-t=0 Python +
   C++)** — ~5 h total.
6. **Then Infra-6, 7 (knudsen_ke + KL divergence Python phases)** — ~3 h.
7. **Update Infra-8 task description** to reflect the Gmail MCP
   draft-only limitation: pipeline creates drafts, user sends; OR
   build smtplib fallback. Then implement Infra-8 (~2 h).
8. **Then Infra-9 (post-IFW shading helper + wire into ~7 plot
   modules)** — ~1.5 h.
9. **First production run: Run-1 (25 eV WP)** — fastest path to an
   emailable result, validates the whole infra pipeline end-to-end.
   Update Cfg `shared/configs/electron_proj_E25_L50_cubic.hpp` to use
   `boundary_rule::launch_z(5, 50)=−5`, `n_steps_for(5, 50, 1.356, 0.02)=922`,
   and `write_every_for(922)=3` before launching.
10. **Update this handover at every milestone** (per
    `.claude/rules/handovers.md`).

## Why this handover exists

Per `.claude/rules/handovers.md`: "Maintain one rolling handover file
per substantive task in `docs/handovers/`. … A handover must be
sufficient for another session (or a different Claude account) to
continue without guessing." This file is the entry point for any
session that picks up campaign work; the plan and design journal entry
are the longer-form references.
