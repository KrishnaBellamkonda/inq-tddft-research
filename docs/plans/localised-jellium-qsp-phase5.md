# Plan — qsp_phase5: WP quantum stopping-power velocity sweep S(E)

Status: **locked 2026-06-26** (designed via grill-with-docs). Last phase of the
localised-jellium scattering campaign.

## Goal

Measure the **quantum (wavepacket) electronic stopping power** of a σ_WP=0.5
projectile in the localised jellium slab across a velocity grid, plot it as
**S(E)**, and overlay the existing **bulk** classical (σ_WP=0.5) and **bulk
Lindhard** references. The plot is rebuilt + emailed after *each* run, fully
autonomously (no user, no Claude in the loop).

## Decisions (all locked)

- **Grid (S vs drift energy ½k₀²·27.211 eV):** WP points at
  **{23, 54, 122, 218, 340, 490} eV** ↔ v ∈ {1.3, 2.0, 3.0, 4.0, 5.0, 6.0}.
  v=2.0 (54 eV) is **reused from phase 4** (no rerun). **5 new runs:**
  v ∈ {1.3, 3.0, 4.0, 5.0, 6.0}.
- **System:** reuse phase-4 localised GS `shared_gs/slab_n82_L50x50x90`, box
  50×50×90, slab |z|<12.5 (L_z=25 Bohr), two-sided sin² CAP (η=−0.7, faces ±35),
  launch z=−23.75, dt=0.04. **Only k₀ and τ vary** (energy is real-time-only ⇒ GS
  reuse valid). Single env-driven binary (`LJ_K0`, `LJ_N_STEPS`).
- **τ per velocity:** τ ≈ 200/v (anchored to phase-4 v=2.0→τ=100), cap 200.
  wall ≈ 0.054·τ h (anchors: phase-4 5.07 h, phase-3 5.74 h at τ=100).
- **Convergence:** WP S = [E_total(t_f) − E_GS]/L_z (energy method, **E_GS anchor**
  — the WP's drift KE lives in E_total(0); E_GS=−70.22568216820937 Ha). Gate:
  norm_f<0.02 & |late dE/dt|<0.2 eV/au. High-v converge → true values; v=1.3 (& 54)
  are **upper bounds** (late slope <0 ⇒ residual WP energy still draining). Each
  point flagged; bounds drawn as down-arrow markers.
- **Overlays (all labelled *bulk reference*):**
  - Classical σ_WP=0.5 = **σ_q=0.354 = `sigma0p35` set** (the √2 catch — NOT
    `sigma0p5`), extracted by the bulk slope method (`sigma_sweep_report.extract`).
  - Bulk Lindhard `stopping_power_point(v, kF)`, r_s=5.69 (≈ slab 5.666).
  - The one localised park-method classical point (v=2.0, **0.249 eV/Bohr**) as a
    marked geometry-matched check.
- **Guards (skill `stopping-power-extraction`):** N_total conservation
  (bath drain <2%); the WP-absorbed norm is the projectile leaving, by design.

## Files (all under `ResearchProject/systems/localised_jellium/`)

| File | Role |
|---|---|
| `scripts/qsp_phase5/wp/run.cpp` | phase-4 WP run.cpp + `LJ_K0` env override; one build for all v |
| `scripts/qsp_phase5/run_sweep.sh` | value-first 2-GPU dispatcher + smoke gate + per-run chain |
| `hypotheses/qsp_phase5/analyse_phase5.py` | per-run WP QSP (energy method, guard, bound) → `results_<tag>.json` + append `se_state.csv` |
| `hypotheses/qsp_phase5/build_se_plot.py` | reads `se_state.csv` → S(E) overlay plot + per-point table → email (threaded `[lj-wp-se-sweep]`) |
| `hypotheses/qsp_phase5/build_phase5_notebook.py` | study notebook (cumulative S(E)) |
| `.claude/skills/run-notebook/run_notebook_builder.py` | + WP energy-method QSP section (`--e-gs-ha`, `--l-slab`); soften loss-fn note |

## Dispatcher flow (autonomous, shell+Python only)

1. **Smoke gate:** v=6.0, τ≈10 a.u. (~250 steps). Confirm it propagates + writes
   observables/norm. On failure → email abort, stop. On pass → delete smoke, proceed.
2. **Production, value-first on 2 GPUs** (clean+cheap first; marginal v=1.3 last):
   - GPU0: v=6.0(490) → v=5.0(340) → v=4.0(218) → v=3.0(122)
   - GPU1: v=1.3(23)   (the 8.3 h long pole, runs in parallel)
   - ETA ≈ 10 h; partial-data-safe (high-quality converged points land first).
3. **Per-run chain** (fires the instant a run prints `done. wall`):
   `analyse_phase5.py <tag>` → append `se_state.csv` → `build_se_plot.py --email`
   → run-notebook for that run. Then the dispatcher starts the next on that GPU.
4. On all-done: build study notebook, write `POSTPROC_DONE`.

## Validation

- Smoke gate (above) is the pre-flight (Tier-A) test.
- analyse_phase5 reuses the **validated** energy-method (independently reproduced
  vs analyse_phase4 = 2.39 eV/Bohr) + the `stopping-power-extraction` guards.
- Each emailed plot states per-point convergence; bounds are not presented as
  converged values.
- Test-catalogue row added for `build_se_plot` (overlay set + Lindhard locator).

## Out of scope / caveats

- Slab-WP vs bulk-classical is a **geometry mismatch** (ADR
  `0010-localised-wp-vs-bulk-reference`) — an estimate, labelled as such.
- v=1.3 is zero-point-marginal (k₀/σ_p=1.3); v<1 excluded (k₀<σ_p, velocity
  ill-defined for σ=0.5).
- 218/340/490 eV have **no** classical point (>v=3.0); only Lindhard extends there.
