# Handover: Stopping power from the decomposed energy ledger

Rolling handover for campaign `stopping-from-energy-decomposition` (area
`localised_jellium`). Campaign prompt:
`/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium/stopping-from-energy-decomposition.md`.
Session: `formula-sp-from-energy-decomp`. Started 2026-07-21.

---

## Goal

Derive **Definition 1** — a stopping power built from the decomposed energy ledger
(pairwise Coulomb E_PP/E_PS/E_SS/E_SB/E_PB/E_BB + kinetic split + E_xc) — for BOTH a
classical Gaussian-charge projectile and a wavepacket (WP) projectile in the
localised jellium slab, and apply it to the select plateaued+decomposed runs to make
S(E₀)/S(v₀) plots. One Jupyter notebook holds all four phases. Analysis-only — **no
GPU runs** launched by this campaign. It is the analysis dual of
`classical-highdensity-sv` (which generates the ledger data / owns Definition 2 =
E_absorbed/L_slab).

## Phase structure (6 phases, one notebook)

Expanded 4→6 phases 2026-07-22 after the user's notes
(`docs/notes/stopping-power-formula-from-energy-decomposition.md`) + the catalogued
plan (user: "add the plan for phases … everything is fine").

1. **Understand + propose formulae** — HUMAN GATE. Now also: **§C0 measure E_sp(0)
   vs E_sp(t_final)** + decide its treatment; **§D setup-simplification brainstorm**.
2. **Implement + validate the kernel** — code-test + formula-validation; E_sp term
   exposed switchably.
3. **Apply to the select runs** — USER-ASSISTED run pick. Now: C1 + **C2 full
   plotted chain** (ΔE(t), ΔKE_proj(t), fitted region, best-fit line, value) + C3
   (ZPE-subtracted) + **C1-vs-C2 classical agreement**.
4. **A1 discrepancy investigation** (NEW) — why deposit-S >> Lindhard bulk; find the
   plateaued-WP runs; **rough E_sp-subtraction recompute** (`D_corr = D_raw +
   E_sp(0)`); test hypotheses (a) localised≠bulk, (b) decomposition-sum vs E_total,
   (c) CAP anomaly; conclusion. Fable 5.
5. **Aggregate + plot** — S(E₀) headline + S(v₀) companion, classical & WP overlaid.
6. **Synthesis + setup recommendations** (NEW) — fold A1 conclusion + §D into a
   next-batch recommendation.

Model routing (user): **Fable 5** for formula brainstorm/derivation + the A1
reasoning; **Sonnet** for data-fetch + simple per-run analysis.

## E_sp rough-calc (immediate, pre-Phase-4 probe — interviewing user 2026-07-22)

User asked for a rough calc on the plateaued localised-jellium quantum-stopping
runs: since E_sp(0) is non-negligible but E_sp(t_final)≈0 (excess projectile
absorbed), subtract the initial slab–projectile interaction from `E_absorbed` and
re-derive S, checking whether it lands near Lindhard bulk. Blocking ambiguities
being interviewed BEFORE compute: (i) which baseline curve the original "too huge"
S used (E_total vs E_electronic vs decomposition-sum); (ii) whether to find the runs
now or have the user point to them (run-pick is otherwise deferred to Phase 3);
(iii) exact correction recipe + Lindhard basis (r_s of the slab). Recipe under test:
`D_corr = D_raw + E_sp(0)` with `E_sp ≡ e_ps` from `interactions.csv`.

## Status

- **DONE — grill + campaign authored.** Extended grill-with-docs settled: WP
  definition tension (brainstormed a full A–E menu; projectile-KE forbidden for WP
  headline); data scope (existing runs only; user identifies the select
  decomposition+plateau runs at Phase 3); one notebook all phases; new *separate*
  campaign bound to `classical-highdensity-sv`. Campaign file written, `status:
  draft`, 4 tasks. INDEX rebuilt (38 campaigns). CONTEXT.md glossary section added
  ("Stopping power from decomposed energy").
- **DONE — Phase 1 notebook built + executed (0 errors), NOW enriched to §A/B/C0/C/D
  (24 cells, 2026-07-22).**
  `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/stopping_from_decomposition/stopping_from_decomposition.ipynb`
  (builder `build_stopping_from_decomposition_report.py` beside it, uses the house
  `hypotheses/_nbreport.py`). §A inventory, §B closure (~1e-10 Ha), **§C0 E_sp
  measured on the r12 ledger twin** (raw E_PS≈−140 eV screened by E_PB≈+135 → net
  −5.2 eV; qsp dynamic proxy +4 eV; refutes E_sp as the ~8× cause), §C formulae
  C1/C2/C3/C4, **§D 9-idea setup brainstorm (Fable 5; top pick run-to-extinction +
  vacuum-twin)**. Ends at the HUMAN GATE.
- **DONE — E_sp rough-calc + run identification (2026-07-22).** The "first S(v) plot"
  = `qsp_phase5/figs/se_quantum_stopping.png` from `qsp_phase5/se_state.csv`; source
  runs `scripts/qsp_phase{3,4,5}/wp/results/` — **wavepacket, not psp** (classical
  gave S=0). These runs carry NO `interactions.csv`. Deposit baseline = E_GS (not
  E_total(0)); E_GS=−70.2257 Ha. E_sp correction (~4 eV net) moves S ~2.4→~2.2
  eV/Bohr, still ~8× Lindhard (0.28–0.45). Real drivers: WP zero-point KE 81.6 eV +
  drift ~180 eV + approximate CAP split + non-convergence (WP norm 0.05–0.09 at t_f).
  → Phase-4 hypotheses (b)/(c). Full numbers in `qsp_phase{3,4}/results.json`.
- **DONE — Phase 1 gate passed + §E critical analysis (2026-07-22).** User directed
  proceeding WITHOUT re-runs and answering every open question in the notebook. Added
  **§E — Critical answers** (30 cells total, 0 errors; Fable 5 stress-tested the
  verdicts). Findings, all data-grounded on the existing qsp runs:
  - **Central question (why WP S ~8× Lindhard):** decisive fact — at v=1.3/2.0 the
    deposit (59 eV) *exceeds the projectile drift KE* (23/54 eV), so the surplus is
    the WP's own internal energy, NOT stopping. **(a) localised≠bulk REFUTED** (the
    classical projectile in the SAME geometry gives S=0.25 eV/Bohr ≈0.5× Lindhard —
    geometry fine); **(b)+(c) are ONE mechanism** — `deposit=E_total(t_f)−E_GS`
    charges the slab with the WP's fixed 81.6 eV zero-point + self-interaction + drift
    KE, minus only what the CAP exports, and the CAP exports norm far better than
    energy. Deposit pinned ~0.73× the 81.6 eV zero-point at low v; retains ~96% of the
    drift KE at high v. Lindhard-expected deposit is only 7–11 eV → every row 5–50×
    too big. **No WP row is a stopping power, even as an upper bound.**
  - **Definitional:** C1-excludes-E_sp justified but the "<10%" must be re-scoped
    (few eV vs the true ~7–11 eV deposit is 30–50%); WP separability = diagnostic-grade
    not measurement-grade; C1=C2 only where CAP flux is negligible (C2 fails for a
    parked projectile — light-projectile rule); C3 vacuum-referenced can remove the
    zero-point by construction (keep as cheap cross-check).
  - **5 new hypotheses** (t=0 baseline already contaminated → measure BATH-ONLY rise;
    xc self-interaction of the extra orbital; CAP reflection of slow tails; σ_WP sweep
    as the decisive artifact test; audit the −47 eV bath entry).
  - **Bottom line:** only the classical S=0.25 eV/Bohr is trustworthy (factor-~2);
    the fix is a **rebuilt observable** (bath-only decomposed / vacuum-referenced),
    not more runs — validated by a σ_WP sweep.
- **BLOCKED (needs user re-run decision)** — Phases 2–6. §E's own conclusion: further
  progress needs a ledger-carrying *dynamic, fully-absorbing* run (§D #1/#3) — none
  exists (the gap §C0 found). Phase-2 C1 self-test (decomposition-S == Def-2 deposit)
  also needs that run; only the synthetic-fixture + closure tests are runnable now.
  User deferred re-runs → awaiting the go decision.
- **DEFERRED (by user)** — the select plateaued+decomposition runs are identified
  WITH the user at the start of Phase 3, after Phases 1–2 complete.

## Verified vs unverified

- **VERIFIED (executed in the notebook / by Sonnet audit):** closure holds to ~1e-10
  Ha for both representations — classical (`energy_hartree=E_SS`,
  `energy_external=E_SB+E_PS`) and WP (`energy_hartree=E_SS+E_PS+E_PP`,
  `energy_external=E_SB+E_PB`), plus `energy_total` sum, on the `twin_ec_rsweep`
  r12 pair and `phase5_wp/p5_null_s2_k4_wp`. Column inventory across run types
  confirmed on disk. Closure asserted in
  `inq-stack/include/inqkit/jellium/interaction_energies.hpp:17-18`.
- **UNVERIFIED / open:** the Definition-1 formula is NOT locked (Phase-1 gate
  pending); the classical S(E)/S(v) benchmark data depends on
  `classical-highdensity-sv` completing its Phase-3 sweep (currently `status:
  draft`, all tasks not done); C3 (WP projectile-partition) needs a WP-in-vacuum
  twin run + a passing gauge test before it can be reported.

## Reference runs (Phase-1 closure demo only — NOT stopping data)

- classical: `.../localised_jellium_dynamics/runs/twin_ec_rsweep/results/r12_classical`
- WP:        `.../localised_jellium_dynamics/runs/twin_ec_rsweep/results/r12_wp`
- long WP:   `.../localised_jellium_dynamics/phase5_wp/results/p5_null_s2_k4_wp`

Data-availability constraint (Sonnet audit 2026-07-21): the runs that carry
`interactions.csv` (phase5_wp, twin_ec_rsweep) are closed-boundary / short and do
NOT stop; the runs that plateau after depositing (qsp_phase3/4/5) lack the pairwise
ledger. The runs that have BOTH are the still-to-run `classical-highdensity-sv`
sweep + whatever select runs the user names at Phase 3.

## Proposed formulae (Phase 1 §C — awaiting approval)

- **C1 headline (both runs):** `S = ΔE_electronic/Δs`, split
  `ΔE_target = ΔKE_slab + ΔE_SS + ΔE_SB + ΔE_xc`; classical self-test == Def-2
  `E_absorbed/L_slab`.
- **C2 classical anchor (sanity):** `S(v0) = −d(½mv²)/ds` over v≥0.85·v0.
- **C3 WP partition (exploratory, gauge-gated):**
  `E_proj^WP = KE_proj + E_PP + E_PS + E_PB`, vacuum-baseline-subtracted;
  `S = −d[E_proj^WP − E_proj^vac]/ds`. Only if gauge test Δ(E_SS,E_SB,E_BB)≈0.
- **C4 qualifiers:** A1 late-time absorption; D1 matched-face; D2 hysteresis loop.

## Next actions (resume here)

1. **User reviews Phase 1 §A/§B/§C in the notebook and approves the formula set**
   (the gate). Record which candidates are locked.
2. **Phase 2:** implement the approved formula as a tested extraction kernel
   (skill-local or `inqview.analysis`); known-case tests (classical
   decomposition-S == Def-2 deposit within tol; closure ~1e-9; synthetic); dispatch
   a `formula-validation` agent; add a catalogue row; lock.
3. **Phase 3:** with the user, identify the select decomposition+plateau runs;
   apply the kernel per run with the gauge + conservation + plateau sanity checks.
4. **Phase 4:** aggregate by system → S(E₀)/S(v₀) figures (classical & WP overlaid)
   + Lindhard eyeball overlay + `S_decomposition_summary.csv`; execute the full
   notebook; flip frontmatter `done` flags; refresh INDEX.
