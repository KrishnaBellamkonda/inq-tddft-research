# Handover: matched WP-vs-classical induced-wake case study

## Current status
- Grill complete, plan agreed (`docs/plans/wp_vs_classical_matched_wake.md`).
- New case-study script written and being run:
  `ResearchProject/systems/jellium/scripts/wp_vs_classical_matched.py`.
- tddft-simulations SKILL §3b updated: density_wp now COMPULSORY at equal
  cadence for Jellium WP wake-analysis runs.
- Old `docs/presentations/storyline/tasks/batch2_figures/` to be DELETED by the
  user; replacement figures go in `.../tasks/wp_vs_classical_matched/`.

## CORRECTED TODO (supersedes the "wp density never saved" premise)
The premise "wavepacket density has not been saved for any run, so all
comparisons are meaningless" is **FALSE as stated**. Verified 2026-06-01 by
frame counting:
- `_wf` runs (σ0.5/3/8 @E100): density_wp at FULL cadence (== density_total).
- `_v2` runs (σ1, all energy-sweep): density_wp saved but ~10× COARSER than
  density_total (σ1_v2: 32 wp vs 317 total). → only sparse exact subtraction.
- original (non-v2/non-wf) runs: **0** density_wp frames → genuinely unusable
  for any wake / density-difference analysis.

Real action items:
1. Any FUTURE jellium WP run intended for wake analysis MUST save density_wp
   at equal cadence to density_total (now mandated in the skill). DONE (skill).
2. The batch2 WP−classical figures were built on MISMATCHED run pairs (see
   below) — not a data problem, a pairing problem. Being replaced.

## What was actually wrong with batch2 (the real bug)
Frame ORDERING was already correct (driver sorts WP times ascending). The
defect is **run-pairing**: the differenced WP and classical runs do not share
launch position / dt / end-time:
- σ-sweep WP runs launch at `boundary + 4σ` → z = −23/−21/−13/−1 Bohr for
  σ = 0.5/1/3/8, all differenced against classical_v2 at z = −21.
- σ0.5/3/8 `_wf` WP are dt=0.02; classical_v2 is dt=0.01. E20/E25 WP dt=0.01
  vs classical dt=0.02.
- Different end-times/cadences → sampling classical at WP physical times
  aliases ("classical Δn jumps around between timesteps").

The faint WP trail is mostly PHYSICAL: a σ=1–8 Gaussian electron polarizes the
bath far more weakly/diffusely than a point classical projectile.

## Chosen basis (only existing fully-matched pair)
`run_wp_n162_L50_E100_sigma1_v2` (dt=0.01, launch z=−21) ↔
`run_classical_n162_L50_E100_v2` (dt=0.01, launch z=−21). Both L=50, E=100 eV.
Sampled at the WP's 32 exact density_wp frame times, t∈[0,9.30] a.u.; classical
truncated to the WP end.

## Files touched
- NEW `docs/plans/wp_vs_classical_matched_wake.md`
- NEW `ResearchProject/systems/jellium/scripts/wp_vs_classical_matched.py`
- EDIT `.claude/skills/tddft-simulations/SKILL.md` (§3b density_wp compulsory)
- NEW `docs/presentations/storyline/tasks/wp_vs_classical_matched/` (outputs)
- NEW this handover

## Tests and validation (PASS, 2026-06-01)
Known-case checks (printed by the script):
- Δn(t0) max|.| = 0 for both runs (==0 by construction). ✓
- ∫n_system dV = 162.000 at t0/mid/end for BOTH runs → exact wp subtraction,
  charge conserved. ✓
- WP centroid −21.0 → 4.1 Bohr, monotonic (single pass). ✓

## Concrete result (the "message")
Classical projectile is a CLASSICAL ELECTRON (custom electron-ONCV-1.2.upf,
mass=m_e, charge −1) — same sign/mass/velocity/box/dt as the WP, so a fair
quantum-vs-classical comparison. Over frames with t>1 a.u.:
- peak |Δn|: WP 1.73e−2 vs classical 2.59e−3 e/Bohr → **WP wake 6.7× STRONGER
  in peak**, 11.7× larger integrated, effective width 17.2 vs 9.8 Bohr.
- The σ=1 WP trail is therefore NOT invisible — the earlier impression came
  from mismatched / large-σ batch2 figures, not missing data.
- CAVEAT: the classical electron is a pseudopotential (cusp-softened), so part
  of the gap is reduced projectile–bath coupling, NOT purely quantum
  delocalization. A bare-Coulomb / cusp-matched classical electron is needed to
  isolate the pure quantum effect. (Ties to the user's cusp-pseudopotential
  concern.)

## Known issues / blockers
- σ-sweep and energy-sweep WP−classical comparisons CANNOT be salvaged from
  existing runs (no position+dt-matched classical for σ0.5/3/8; dt mismatch for
  E20/E25). Would need new matched classical twins — user declined for now.
- The clean pair has only 32 wp frames (sparse). Adequate for a case study;
  a full-cadence twin pair would be the rigorous follow-up.

## TODO — in-depth density-feature comparison study (user-requested 2026-06-01)
**Goal.** Use the `n_system^WP − n_system^classical` infrastructure we built
(`inqview.postprocess.wake` + `scripts/wp_vs_classical_matched.py`) to gain
PHYSICAL INTUITION for how the induced bath-density features depend on:
  (A) projectile ENERGY (WP vs the matched classical electron at each E), and
  (B) WP CONCENTRATION / width σ (how delocalization reshapes the wake).
This extends the single σ1@E100 case study to a 2D (E, σ) grid of matched
WP↔classical pairs and reads off the trends.

**Method (reuse, don't reinvent).**
- For each (E, σ): build the matched difference field `D = Δn^WP − Δn^classical`
  exactly as `wp_vs_classical_matched.py` does (exact same-step n_total−n_wp,
  t0-subtracted, sampled at the WP's exact density_wp frame times, classical
  sampled at the same physical times, truncated to WP end).
- Quantify the SAME feature set already in `metrics.csv` per (E, σ, t): peak
  |Δn|, ∫|Δn|, effective width (∫/peak), trailing-oscillation wavelength,
  depletion-vs-enhancement integrals behind the centroid, and the lag of the
  wake peak behind the WP centroid. Add: leading vs trailing asymmetry, and
  the WP−classical residual amplitude as a function of E and σ.
- Produce two sweep figures: feature-vs-E (at fixed σ=1) and feature-vs-σ
  (at fixed E=100), each overlaying the WP, the classical, and the difference.

**Physical hypotheses to test (intuition we expect).**
- Energy axis: as E↑ (faster projectile), the WP wake should sharpen and the
  WP−classical difference should SHRINK toward the classical limit (ties to the
  M10 classical-limit question). At low E the difference should be largest.
- Concentration axis: as σ↓ (WP → point charge) the WP should approach the
  classical electron, so `D → 0`; as σ↑ the WP wake becomes broader/weaker in
  peak but larger in spatial extent. This isolates the delocalization effect.

**HARD PREREQUISITE (blocker — must be done first).** This study needs MATCHED
twin pairs that DO NOT YET EXIST except σ1@E100:
  - one classical-electron run + one WP run per (E, σ) cell,
  - IDENTICAL dt, launch z, N_STEPS/end-time, box, GS,
  - density_wp saved at EQUAL cadence to density_total (now mandated in the
    tddft-simulations skill §3b).
Existing runs fail this (σ-sweep WP launch at boundary+4σ ≠ classical −21; dt
mismatches; _v2 wp-density ~10× too coarse). So step 1 of this study is to
DISPATCH the matched twin grid via the tddft-simulations skill. Suggested grid:
  E ∈ {20, 50, 100, 300} eV at σ=1, and σ ∈ {0.5, 1, 3, 8} at E=100 eV
  (σ1@E100 reused as the shared corner). ~7 new twin pairs (14 runs).
Also worth pairing with the bare-Coulomb classical electron (cusp caveat above)
so the energy/σ trends are not confounded by pseudopotential softening.

## Exact next steps
1. DONE: script outputs + known-case results confirmed; metrics in this file.
2. User reviews figures (no-preview rule) in the new folder.
3. User decision: authorise the matched twin (E, σ) grid (14 runs) needed for
   the in-depth study above. Until then the study is blocked on data.
4. (Optional) bare-Coulomb classical-electron run to decompose the 6.7× wake
   gap into quantum-delocalization vs pseudopotential-softening.
