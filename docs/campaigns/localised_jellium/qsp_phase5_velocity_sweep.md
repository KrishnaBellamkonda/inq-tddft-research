---
id: locjel-qsp-phase5-se-sweep
area: localised_jellium
title: "Phase 5 — quantum WP stopping power S(E) velocity sweep"
status: ready
hypothesis: "The quantum (σ_WP=0.5 wavepacket) electronic stopping power of a localised jellium slab, measured as S=[E_total(t_f)−E_GS]/L_z across v∈{1.3,3.0,4.0,5.0,6.0} (+ reused 54 eV), exceeds the matched-width bulk classical and bulk point-charge Lindhard references across the whole velocity range, with high-v points converging to true values and slow points bounded from above."
handover: docs/handovers/localised-jellium.md
tasks:
  - { name: "Smoke gate (v=6.0 short build+run)", done: false }
  - { name: "WP run v=6.0 (490 eV)", done: false }
  - { name: "WP run v=5.0 (340 eV)", done: false }
  - { name: "WP run v=4.0 (218 eV)", done: false }
  - { name: "WP run v=3.0 (122 eV)", done: false }
  - { name: "WP run v=1.3 (23 eV)", done: false }
  - { name: "Cumulative S(E) plot + email per run", done: false }
  - { name: "Study notebook qsp_phase5_study.ipynb", done: false }
blocked_reason: ""
---

# Phase 5 — quantum WP stopping power S(E) velocity sweep

The **last phase** of the localised-jellium scattering campaign. Designed
2026-06-26 via grill-with-docs; machinery implemented and validated. This prompt
is the resumable autonomy spec — a fresh agent (or the dispatcher itself) runs it
end-to-end with **no user in the loop**.

<identity>
Scientific-computing researcher running first-principles rt-TDDFT (INQ engine +
inqkit/inqview). Adheres to the repo rules (σ_WP convention, GPU default, file
placement, validation gates, number rounding, gate-free loss functions).
</identity>

## Question (hypothesis)

Does the σ_WP=0.5 **wavepacket** deposit more energy per unit path than the
matched-width **classical** projectile and the point-charge **Lindhard** estimate,
across v? Build a single **S(E)** curve to answer it, updated + emailed after each
run.

## Run matrix (locked)

- **System (reused unchanged from phase 4):** localised jellium slab, GS
  `shared_gs/slab_n82_L50x50x90`, box 50×50×90, slab |z|<12.5 (L_z=25 Bohr),
  two-sided sin² CAP (η=−0.7, faces ±35), launch z=−23.75, dt=0.04, LDA. Only the
  drift momentum k₀ and run length τ vary (energy is real-time-only ⇒ GS reuse OK).
  Single env-driven binary (`LJ_K0`, `LJ_N_STEPS`).
- **Grid (S vs drift E=½k₀²·27.211 eV):** WP at {23, 54, 122, 218, 340, 490} eV ↔
  v∈{1.3, 2.0, 3.0, 4.0, 5.0, 6.0}. **v=2.0 (54 eV) is reused from phase 4** (no
  rerun). **5 new runs:** v∈{1.3, 3.0, 4.0, 5.0, 6.0}.
- **τ per velocity:** τ≈200/v (cap 200), anchored to phase-4 (v=2.0→τ=100). wall
  ≈0.054·τ h. Total ≈10 h on 2 GPUs, value-first.
- **Method (per run):** S=[E_total(t_f)−E_GS]/L_z, **E_GS anchor** (=−70.22568216820937
  Ha; the WP's drift KE lives in E_total(0)). Convergence gate norm_f<0.02 &
  |late dE/dt|<0.2 eV/au; else **upper bound** (down-arrow marker). N_total guard <2%.
- **Overlays (bulk references, labelled — ADR 0010):** classical σ_WP=0.5 = bulk
  σ_q=0.354 `sigma0p35` set (the √2 catch); point-charge Lindhard (r_s=5.69);
  + the lone localised park point (v=2.0, 0.25) as a geometry-matched check.

<observables_set>
The standard WP suite the phase-5 binary already emits: density VTIs
(total/system/gs/wp + wavefunction), observables (energies/current/dipole/L2),
state_energies, occupations, eigenvalues, momentum_distribution, wp_momentum_stats,
wp_real_space_stats, overlap/overlap_full (KS excitation), density_delta. The
quantum stopping power S is auto-measured from observables.csv each run.
</observables_set>

## How it runs (autonomous, shell+Python only)

One command — the dispatcher owns everything:

```bash
cd /local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/qsp_phase5
nohup bash run_sweep.sh > sweep.log 2>&1 &
```

1. **Smoke gate:** builds (`INQ_SOURCE=inq-study`) + a short v=6.0 run; verifies
   propagation + observables. On failure → emails an abort and stops (no
   production). On pass → deletes smoke, proceeds.
2. **Production, value-first on 2 GPUs:** GPU0 chains v6.0→v5.0→v4.0→v3.0 (clean,
   converged, cheap first); GPU1 runs the long marginal v1.3 alone. Partial-data-safe.
3. **Per-run chain** (the instant a run returns 0): `analyse_phase5.py` (QSP +
   convergence + guard, upserts `se_state.csv`) → `build_se_plot.py --email` (S(E)
   overlay plot, threaded `[lj-wp-se-sweep]`) → per-run notebook. State writes are
   flock-serialised.
4. **All done:** `build_phase5_notebook.py` (study notebook) + `POSTPROC_DONE`.

## Guard rails / stop conditions

- Smoke gate must pass before any production run (else abort + email).
- A failed individual run (rc≠0 or no `run_completed = true`) is **skipped**; the
  sweep continues with the remaining velocities (partial S(E) still emailed).
- Zero-point floor: v>σ_p=1.0 (k₀>σ_p) for a meaningful velocity; v=1.3 is marginal
  (flagged), v<1 excluded. Do NOT add v<1 points.
- σ convention: classical overlay is the σ_q=0.354 (`sigma0p35`) bulk set —
  **never** the `sigma0p5` set (√2 trap).
- Loss-function / spectral analysis is **never gated** by run length (see
  `feedback_fourier_loss_function_gate`).

<preflight>
- [x] GS exists: `shared_gs/slab_n82_L50x50x90`.
- [x] Binary builds (smoke gate doubles as the build; `INQ_SOURCE=inq-study`).
- [x] Both GPUs free (checked 2026-06-26; NVML mismatch is cosmetic, compute OK).
- [x] Method validated: `analyse_phase5` reproduces phase-4 energy method (2.39 eV/Bohr).
- [x] Email path: `inqview.email.send_run_email` → chiddukanna@gmail.com.
- [x] Plan: `docs/plans/localised-jellium-qsp-phase5.md`. Decision: ADR 0010.
- [x] Wall budget: ≈10 h on 2 GPUs ≪ 1-day target; ≈18.6 h on 1 GPU (still fits).
</preflight>

## Files

| File | Role |
|---|---|
| `scripts/qsp_phase5/wp/run.cpp` | env-driven WP binary (`LJ_K0`, `LJ_N_STEPS`) |
| `scripts/qsp_phase5/run_sweep.sh` | autonomous dispatcher (smoke gate + 2-GPU + per-run chain) |
| `hypotheses/qsp_phase5/analyse_phase5.py` | per-run QSP → `se_state.csv` + `results_<tag>.json` |
| `hypotheses/qsp_phase5/build_se_plot.py` | cumulative S(E) plot + threaded email |
| `hypotheses/qsp_phase5/build_phase5_notebook.py` | study notebook |
| `.claude/skills/run-notebook/run_notebook_builder.py` | + WP energy-method QSP section |
