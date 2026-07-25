# Plan — Quantum stopping-power ledger + S(v) placement (26-6-26 meeting)

Status: drafted 2026-06-25, awaiting final go-ahead. Grilled via `grill-with-docs`.
Owner deliverable for the 26 Jun 2026 supervisor meeting.

## Goal

Compute the **quantum (wavepacket) electronic stopping power** at E=100 eV from the
`qsp_phase2` τ=40 a.u. twin runs, via the **retained-energy (bath) method**, carefully
tabulated; then place that point — plus the classical-slab point — on the existing
S(v) linear-response plot.

## Locked decisions (from the grill)

1. **Run data:** `qsp_phase2/{wp,classical}/results/{p2_wp,p2_classical}` — τ=40 a.u.,
   50×50×70 box, σ_WP=0.5, E=100 eV (k₀=2.711). Complete; reused (NOT the incomplete
   Phase-3 big-box run, which stalled at t≈15.7/100 a.u.).
2. **Quantum stopping formula (user's retained-energy definition):**
   - `S_WP = [E_total(t_f) − E_jellium(0)] / L_z`, `L_z = 25 Bohr` (slab thickness).
   - `E_jellium(0) ≡ E_total(0) − ⟨T_WP⟩ − E_SIE` (strip the WP kinetic + self-interaction
     from the run's own t=0 total).
   - `E_total(t_f) = E_jellium(t_f)` valid once the CAP has absorbed the WP remnants.
   - **Curiosity check (tabulate):** `E_jellium(0)` vs `E_GS` (bare-slab GS) — expect
     agreement up to a small WP–bath cross-Hartree term.
3. **Values (run-measured / locked):** `E_total(0)=−38.943 Ha`, `⟨T_WP⟩(0)=6.645 Ha`
   (measured `e_kin_ha` step 0; analytic ½k₀²+3/4σ²=6.675), `E_SIE=4.40 eV=0.162 Ha`,
   `E_total(t_f=40)=−43.253 Ha`, `E_GS=−45.759 Ha`. ⇒ `E_jellium(0)=−45.75 Ha`
   (within **0.27 eV** of E_GS ✓), `ΔE=+68 eV`, **S_WP ≈ 2.7 eV/Bohr — UPPER BOUND**
   (convergence gate NOT met at τ=40: WP norm ≈0.14, E_total slope ≈ −1 eV/a.u., not
   plateaued). Report explicitly as an upper bound with the gate diagnostics shown.
4. **Classical stopping (both methods, slab run):**
   - **Headline estimate (user choice):** lowest ion KE during the **first transit** =
     slab centre (t≈12, z≈+1.7, KE_ion≈30 eV) ⇒ `S ≈ (100−30)/23.7 ≈ 3.0 eV/Bohr`.
   - **Honest companion:** equal-potential-face loss (±12.5: 68→54 eV) ⇒ `S ≈ 0.53 eV/Bohr`
     (the conservative mean-field well makes the centre over-count). Show both.
   - The classical `E_total` energy method does NOT work (lowest E_total is t=0; the
     ion–bath interaction swamps the deposited energy; ion wraps + re-enters by t=40).
5. **Rounding:** 2 s.f. default, 3 s.f. max (new rule `.claude/rules/number-rounding.md`).
6. **Output location:** `ResearchProject/systems/localised_jellium/hypotheses/qsp_phase2/`
   (ADR 0007 hypotheses tree). Notebook + figures travel together.

## Deliverables

### D1 — Energy-ledger notebook
`hypotheses/qsp_phase2/quantum_stopping_ledger_26-6-26.ipynb` (+ `_figs/`).
- House narrative: context → formula (every term defined) → reconstructable setup
  (run paths, params) → the careful **energy ledger table** (2/3 s.f.) → S_WP (upper
  bound) + classical estimate → takeaway.
- Tables: (a) WP energy decomposition `E_total(0)→E_jellium(0)` + `E_jellium(0)` vs
  `E_GS`; (b) `ΔE`, `S_WP` with gate status; (c) classical slab-centre + face estimates.
- Diagnostics: WP norm-vs-t + E_total-slope (justify the upper-bound label); classical
  KE_ion(z) showing the well + the two evaluation points.

### D2 — S(v) Plot A
`fig_sv_quantum_point.png` — reuse `build_section1.py::fig3_stopping` machinery:
- Linear-response (point-charge Lindhard) reference curve, r_s=5.69, kF=0.337.
- Classical σ_WP=0.5 **bulk** sweep points (Method A).
- **NEW:** quantum point at v=2.711 a.u. (E=100 eV), S=2.7 eV/Bohr, drawn as an
  **upper bound** (down-arrow / open marker), clearly labelled.

### D3 — S(v) Plot B
`fig_sv_quantum_plus_classical_slab.png` — Plot A + the **classical-slab** point at
v=2.711 (value: 0.53 face / 3.0 centre — DECIDE which to mark; recommend face 0.53 as
the defensible point, with centre as a faint upper tick).

## Validation
- `stopping-power-extraction` skill for the slope/face kernels; `code-test` for any new
  numeric helper; formula grounded in the campaign + handover (cited in the notebook).
- Sanity: `E_jellium(0) ≈ E_GS` (already +0.27 eV ✓); classical face-S ≈ bulk classical
  σ0.5 S at v=2.711 (slab≈bulk cross-check).

## Micro-decisions — RESOLVED (user, 2026-06-25)
- D3 classical-slab marker = **slab-centre 3.0 eV/Bohr** (face 0.53 shown as a lower tick).
- Emit **both S(v) and S(E)** twins for each of Plot A and Plot B (4 figures total).
