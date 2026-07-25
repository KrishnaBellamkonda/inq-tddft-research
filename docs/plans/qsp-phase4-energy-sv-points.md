# Plan — Phase 4: energy-matched S(v) point (54 eV WP + redesigned classical)

Rolling plan for the next localised-jellium stopping run-set. Companion handover:
`/local/data/public/skcb2/tddft/docs/handovers/localised-jellium.md`.
Campaign doc: `/local/data/public/skcb2/tddft/docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md`.

## Goal

Add **one on-grid S(v) data point at E = 54.42 eV (v = 2.0 a.u.)** for the σ_WP = 0.5
quantum wavepacket hitting the localised jellium slab, plus a **matched classical
run**, so the localised-jellium WP S(v) can be overlaid on the bulk-classical
σ_WP = 0.5 S(v) curve at a *shared grid energy*.

The classical S(v) grid (from `ResearchProject/systems/jellium/hypotheses/06_sigma_convergence/sigma_sweep_report.py`,
figure `sv_convergence_energy.png`) is, with E = ½ m v² (m = 1 a.u.):

| v (a.u.) | 0.2 | 0.6 | 0.8 | 1.0 | 1.3 | 2.0 | 3.0 |
|---|---|---|---|---|---|---|---|
| **E (eV)** | 0.54 | 4.90 | 8.71 | 13.61 | 22.99 | **54.42** | 122.45 |

## Decisions (locked via grill, 2026-06-26)

- **ONE energy this phase: 54.42 eV (v = 2.0).** On the classical grid; zero-point-clean
  (k₀/σ_p = 2.0, where σ_p = 1/(2σ_WP) = 1.0). 122 eV (v = 3.0) deferred to a later
  phase. The existing 100 eV WP point is OFF the classical grid (v = 2.71).
- **Two runs: WP + classical**, launched concurrently on the two GPUs.
- **Same setup as Phase 3** (90-box 50×50×90, spacing 0.5, two-sided sin² CAP
  η = −0.7 Ha at ±35..±45 inner-face ±35, equidistant launch z = −23.75, τ = 100 a.u.,
  dt = 0.04, GS **reused** from `shared_gs/slab_n82_L50x50x90`). Only the projectile
  energy changes (v = 2.0).
- **WP run: mechanism UNCHANGED** from Phase 3 (it is trusted). New energy only.
- **Classical run: REDESIGNED.** Ehrenfest dynamical Gaussian-electron ion through
  the slab (real stopping physics) **until it reaches the CAP inner boundary
  (|z_ion| ≥ 35 Bohr)**; at that instant the ion is **parked** (position/velocity
  frozen) and its **Gaussian-Coulomb radial potential is replaced by zero** (it
  becomes an inert tracer). This fixes the Phase-3 classical failure where the light
  electron-mass ion was electrostatically turned at z ≈ −3.1 and reflected (never
  transited → ΔKE/x estimator broke), and prevents post-boundary re-interaction /
  spurious CAP coupling of the ion potential.

## Why the classical redesign

Phase-3 classical (100 eV): ion KE traced 100 → ~0 (turn at z = −3.1) → 18 eV; net
displacement −23.75 → −14.7 Bohr; reflected/trapped, never reached the boundary. At
54 eV (slower) it would reflect even more readily. Parking + neutralising **at the
boundary** captures the transit/turn ΔKE cleanly and stops the inert ion's long-range
Coulomb tail from perturbing the slab or coupling to the CAP afterward.

## Implementation

1. **Energy config** — new header
   `shared/configs/slab_n82_L50x50x90_E54.hpp`, struct `SlabN82_L50x50x90_E54`,
   identical to the 90-box config except `WP_EKIN_EV = 54.42` ⇒ `WP_K0 = 2.0`.
   (Convention: one Cfg header per energy — `feedback_run_pairs_layout`.) Both WP
   and classical run.cpp use `Cfg::WP_K0`, so this sets v = 2.0 consistently.

2. **WP run.cpp** — copy `qsp_phase3/wp/run.cpp`, point `#include` + `Cfg` at the
   E54 config. No other change.

3. **Classical run.cpp** — copy `qsp_phase3/classical/run.cpp`, add park+neutralise.
   **Chosen mechanism: chunked two-segment propagate** (engine-safe; no `inq/` edit):
   - **Segment 1 (Ehrenfest, projectile present):** run `real_time::propagate` in
     `WRITE_EVERY`-step chunks; after each chunk read `ions.positions()[0][2]`. When
     `|z| ≥ CAP_INNER (35)` → stop; record parked z, step, and ΔKE_ion.
   - **Segment 2 (projectile neutralised):** continue propagating the **same
     `electrons`** with the projectile removed / replaced by a zero-charge ghost so
     its potential is gone; run the remaining steps to τ = 100 so plasmons/secondaries
     keep evolving and the CAP absorbs them. Ion track continues to log the frozen z.
   - **Feasibility to confirm in the smoke test:** (a) continuing `propagate` with a
     modified `ions` + existing `electrons`; (b) ion-velocity continuity across chunks.
   - **Fallback if (a)/(b) infeasible:** callback-freeze the ion (park only) and
     **terminate at the boundary** — the ΔKE/x stopping number is already complete at
     that point. Documented as the degraded path.

4. **Stopping extraction:**
   - Classical: S = ΔKE_ion accumulated over the interacting-region path length
     (launch → boundary), and cross-checked by the energy method.
   - WP: energy method (deposited E_total(t_f) − E_GS over the slab thickness), as
     Phase 3.

## Smoke test (GPU, before production) — Tier B

dt = 0.04, short run (~200–400 steps) on **separate GPUs**, both runs:
- WP: clean injection (norm step), no NaN, smooth energy.
- Classical: ion advances, the boundary trigger fires, park + neutralise activates
  (E_total shows the ion-potential drop; ion z frozen thereafter), no NaN, two-segment
  continuity holds. If continuity breaks → fall back (park + terminate).
User approves before the τ = 100 production launch.

## Validation / tests

- t = 0 energy bookkeeping (E_system(0) − E_GS ≈ 0), as Phase 3.
- Reuse `analyse_phase3.py` → `analyse_phase4.py` (parameterised by tag + energy);
  classical S now from the clean transit ΔKE.
- Record a row in `docs/validation/test-catalogue.md` for the new classical
  park-neutralise mechanism.

## Files (under `scripts/qsp_phase4/`)

- `shared/configs/slab_n82_L50x50x90_E54.hpp`
- `scripts/qsp_phase4/wp/run.cpp`, `scripts/qsp_phase4/classical/run.cpp`
- `scripts/qsp_phase4/run_production.sh`, `post_process_phase4.sh`
- `hypotheses/qsp_phase4/analyse_phase4.py` (+ notebook builders, reused)
- GS **reused**: `shared_gs/slab_n82_L50x50x90` (no new GS).

## Guardrails

- `inq/` untouched; all projectile-neutralise logic in the wrapper / run.cpp.
- σ labelled as σ_WP = 0.5 everywhere; UPF stays `electron_gaussian_sigma0p35.upf`
  (σ_pot = σ_WP/√2). 2 s.f. reporting. Nothing committed unless asked.
