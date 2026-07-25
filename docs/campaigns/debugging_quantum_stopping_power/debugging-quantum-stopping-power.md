---
id: debugging-quantum-stopping-power
area: debugging_quantum_stopping_power
title: Debugging quantum stopping power — CAP-capture correction (p5_wp_v1p3)
status: done
hypothesis: "The excess of the p5_wp_v1p3 quantum stopping power (S_WP = 2.37 eV/Bohr, upper bound) over the point-charge Lindhard bulk reference at v=1.3 (r_s=5.666) is explained by the kinetic energy of the wavepacket fraction captured in the box: subtracting E_capt = (N_total(t_f) − 82) × (E_input + E_loc) = n_capt × 104.6 eV from the plateau retained energy brings the corrected S to within 20% of the Lindhard value."
handover: docs/handovers/debugging-quantum-stopping-power.md
tasks:
  - { name: "Reproduce original v1p3 ledger (S=2.37 vs se_state.csv)", done: true }
  - { name: "Assumption-1 CAP-on-bath check (p4_classical N_total drift)", done: true }
  - { name: "E_capt + time-resolved E_capt(t) curve", done: true }
  - { name: "CAP energy-removal ledger (where did the 109 eV go)", done: true }
  - { name: "Revised S + binary verdict vs Lindhard (within 20%)", done: true }
  - { name: "Executed notebook + frontmatter/handover update", done: true }
  - { name: "Sweep extension: correction on all aliasing-valid points (v=1.3-5.0)", done: true }
blocked_reason: ""
---

# Debugging quantum stopping power — CAP-capture correction (`p5_wp_v1p3`)

**Analysis-only campaign** (no new simulations, no GPU). Recomputes the quantum
stopping power of the `p5_wp_v1p3` run with a correction for the wavepacket
fraction captured in the box, and delivers a binary verdict against the Lindhard
bulk reference. Designed 2026-07-11 via the campaigns grill; all decisions below
are user-locked. **2026-07-11 scope correction (user):** the target run is
`p5_wp_v1p3`, NOT `p3_wp` (an earlier same-day execution against p3 was a
mis-read of the user's run choice; superseded by this version).

<identity>
Scientific-computing researcher on the localised-jellium rt-TDDFT stopping-power
project (INQ + inqkit/inqview). Adheres to repo rules: number rounding (2 s.f.,
3 max), σ_WP convention, file placement (ADR 0007), validation gates.
</identity>

## Question (hypothesis)

The v=1.3 quantum stopping power S_WP = 2.37 eV/Bohr (upper bound) sits far above
the point-charge Lindhard bulk value at the same velocity. Part of the retained
energy ΔE is not deposited energy but the drift kinetic energy still carried by
the WP fraction *captured* in the box (never absorbed by the CAP). Does removing
that captured KE bring S within 20% of Lindhard?

**Verdict rule (binary, user-locked):** explained ⟺
|S_corrected − S_Lindhard| / S_Lindhard ≤ 0.20. Report the explained energy as a
table row regardless of verdict.

## Locked decisions

1. **Run (user-corrected 2026-07-11):** `p5_wp_v1p3` — τ=153.8 a.u. (3846×0.04),
   box 50×50×90 Bohr, N_bath=82, r_s=5.666, σ_WP=0.5, E_drift=22.9936 eV
   (k₀=1.3), two-sided sin² CAP (η=−0.7, faces ±35), launch z=−23.75. Raw
   observables:
   `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/qsp_phase5/wp/results/p5_wp_v1p3/raw/observables/`
2. **Original ledger (reproduce first):** S_WP = [E_total(t_f) − E_GS]/L_z,
   L_z = 25 Bohr, E_GS = −70.22568216820937 Ha (phase-3/4/5 shared GS
   `shared_gs/slab_n82_L50x50x90`; also hardcoded in phase-5 `run_sweep.sh`).
   Known values: E_total(t_f=153.6) = −68.0444 Ha ⇒ ΔE = +59.35 eV ⇒
   S = 2.374 eV/Bohr — must match the recorded `se_state.csv` row
   (independent second route).
3. **Captured norm:** n_capt = N_total(t_f) − 82, from `electron_number.csv`
   (t=0 row is 83.000: bath 82 + WP 1). Measured: n_capt ≈ 0.0798.
4. **E_capt (user-locked formula, amended 2026-07-11):**
   E_capt = n_capt × (E_input + E_loc) — the captured density is assigned its
   share of the packet's TOTAL starting kinetic energy:
   - E_input = ½k₀² = 22.9936 eV, the energy inputted to the WP in code
     (`scripts/qsp_phase5/wp/run.cpp:64-65`, `LJ_K0=1.3` from `run_sweep.sh`);
   - E_loc = 3/(4σ²) = 81.63 eV, the localisation (zero-point) KE from σ = 0.5
     (`shared/configs/slab_n82_L50x50x90_E54.hpp:60`);
   - basis verified against the run-measured ⟨T_WP⟩(0) = 104.6 eV.
   Measured: E_capt = 0.0798 × 104.6 ≈ 8.4 eV. (The earlier drift-only variant,
   ≈1.8 eV, is superseded — history in the handover.)
5. **Corrected S:** E_absorbed_jellium = ΔE_plateau − E_capt, with ΔE_plateau =
   E_total − E_GS averaged over the late plateau (last 10% of the run);
   S_corr = E_absorbed_jellium/L_z ≈ (59.4 − 8.4)/25 ≈ 2.0 eV/Bohr.
6. **Reference (user-locked):** point-charge Lindhard bulk ONLY, r_s = 5.666.
   Source curve: `hypotheses/qsp_phase5/lindhard_ref.npz` (same curve as the
   phase-5 S(E) overlays), interpolated at v = 1.3.
7. **Output:** table-only notebook (NO new S(v) figure — user-locked) at
   `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/debugging_quantum_stopping/debugging_quantum_stopping_v1p3.ipynb`
   (new hypotheses folder, user-locked). Direct raw-CSV arithmetic (the shared
   `quantum_stopping_ledger.py` module hardcodes T_drift=100 eV and has no
   phase-5 config — not reusable here); `se_state.csv` is the second route.

## Sanity checks (all four user-locked; raw-CSV arithmetic only)

- **S1 — reproduce-first (known-case test):** recompute S_WP from
  `observables.csv` + E_GS and match the recorded `se_state.csv` value
  2.374 eV/Bohr, and E_jellium(0) − E_GS ≈ +0.4 eV
  (E_jellium(0) = E_total(0) − ⟨T_WP⟩(0) − E_SIE). Abort the campaign
  (status: blocked) if this fails — wrong number being debugged.
- **S2 — assumption 1 (CAP does not eat the jellium):** N_total(t) drift of the
  classical twin `p4_classical` (same 50×50×90 box/CAP the phase-5 sweep reused;
  no WP present; any drift = bath absorption). τ=100 vs 153.8 mismatch is
  flagged; linear extrapolation gives the worst case. Data:
  `scripts/qsp_phase4/classical/results/p4_classical/raw/observables/electron_number.csv`
- **S3 — time-resolved E_capt(t):** (N_total(t) − 82) × E_drift over the full
  run; state whether it has plateaued at t_f (else flag E_capt as upper bound).
- **S4 — CAP energy-removal ledger:** E_total(0) − E_total(t_f) vs the injected
  WP energy (drift 23.0 eV + zero-point ≈81.6 eV + E_SIE 4.4 eV): close the
  books on where the ~109 eV went (absorbed / deposited / captured).

## Guard rails

- No simulations are launched; GPU not required (pure pandas/numpy on existing
  CSVs). If any input file above is missing, set `status: blocked` with the
  path in `blocked_reason` — do not substitute another run.
- Numbers reported at 2 s.f. (3 max where a difference needs it); full precision
  stays in code cells.
- Notebook must execute end-to-end with 0 errors.
- The notebook follows the house narrative (notebook-making): title/question →
  conventions → reconstructable setup (from `run_summary.txt`) → source-file
  links → ledger tables → verdict → takeaway.

<preflight>
- Inputs exist: p5_wp_v1p3 observables.csv / electron_number.csv /
  wp_momentum_stats.csv; p4_classical electron_number.csv; lindhard_ref.npz;
  hypotheses/qsp_phase5/se_state.csv.
- E_GS anchor = −70.22568216820937 Ha (phase-5 run_sweep.sh + GS run_summary).
- S1 reproduce-first gate passes before the correction is computed.
- Verdict rule: |S_corr − S_Lind|/S_Lind ≤ 0.20 ⇒ explained; else not explained.
- On completion: flip tasks done, status → done, update handover, regenerate
  docs/campaigns/INDEX.md via the campaigns skill's build_index.py.
</preflight>

## Provenance / grounding

- Ledger method + E_GS anchor: phase-5 campaign
  `docs/campaigns/localised_jellium/qsp_phase5_velocity_sweep.md`,
  `scripts/qsp_phase5/run_sweep.sh` (EGS constant), and
  `hypotheses/qsp_phase5/se_state.csv` (recorded S per run).
- e_kin_ha in `wp_momentum_stats.csv` is per-unit-norm (verified: at t=0 it is
  3.8456 Ha = drift 0.845 Ha (23.0 eV) + zero-point ≈3.0 Ha (81.6 eV)).
- Lindhard curve: generated by phase-5 `build_se_plot.py` machinery (point-charge,
  r_s = 5.666); reused verbatim, not regenerated.
- E_SIE = 4.40 eV: same σ=0.5 packet on the same r_s≈5.67 slab as p2/p3
  (`quantum_stopping_ledger.py` p3 config).
- Inference (labelled): the captured fraction's KE share "proportional to
  density at the original drift KE" is a user-defined estimator, not a
  literature formula — recorded as the campaign's one free assumption.
