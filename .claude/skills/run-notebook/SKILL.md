---
name: run-notebook
description: Build a deep SINGLE-RUN analysis notebook — the full standardised plot battery for one TDDFT run, for analysing it individually in depth. Use after a run completes when the user wants a per-run deep-dive (distinct from the run-SET study notebook of notebook-making). Skill-local, shippable builder that assembles over inqview.pipeline.
---

# Run Notebook

A **run-notebook** is the deep, single-run analysis artefact: every key plot for
**one** run, so a reader can analyse that run individually in depth. It is the
sibling of the **study notebook** (the `notebook-making` skill, which narrates what
a run-*SET* means) — a run-notebook narrates what **one run shows**.

It is an **assembler over `inqview.pipeline`**, never a reimplementation: it runs
the relevant pipeline phases (which compute the figures, auto-skipping what a run
can't produce), then embeds those figures into an executed `.ipynb` with the
house-narrative context. Built by this skill's **skill-local, shippable builder**
(`run_notebook_builder.py`); the generated `.ipynb` lives in the run's
`hypotheses/<sweep>/` folder (ADR 0007), never inside the skill.

## When to use

- A run completes and the user wants a **per-run deep-dive** (not the sweep study).
- The user says "make a run-notebook for <run>" or wants all key plots of one run.
- Do **not** use for a run-SET narrative (that is `notebook-making`); a run-notebook
  is per single run.

## House narrative — context first, then the battery

Sections **1–4 are the `notebook-making` context** (do not skip, do not reorder):

1. **Title + the question** — name the run (system, projectile, σ, energy, CAP/no-CAP)
   and what this single run is meant to show.
2. **Conventions & symbols** — unit system, the unified σ convention (σ = σ_wp;
   charge std = σ/√2), symbol→meaning table.
3. **Setup — fully reconstructable** — pulled from `run_summary.txt`: cell/geometry,
   grid, electrons, propagator, dt, N_STEPS, τ, CAP (η, width) or "no CAP", launch z0,
   the projectile (WP σ_wp/k0 or classical UPF), the GS reused.
4. **Source files** — linked by repo-relative path: the `run.cpp`, the dispatcher,
   this run's `analyse.py`, the inqview phases/kernels used, this builder.

Then the **standard battery** (each section AUTO-GATED on the observable existing,
so the notebook adapts to WP / classical / baseline runs):

| group | panels |
|---|---|
| **Visual intuition** | one lead **xz density GIF** (total) + per-run energetics |
| **Density-GIF battery** (standard, auto) | the matrix **{n(x,z,t), Δn=n(t)−n(0), Δn=n(t+dt)−n(t)} × {total, wp, bath}** — WP runs get **9 GIFs**, classical **3** (total; Δn = induced wake). Rendered by `inqview.visualisation.make_density_gif_battery` (bath = n_total − n_wp). **ENFORCED (2026-07-11): every battery GIF renders LINEAR \| LOG side by side** — density = linear + `LogNorm` (total/bath share ONE slab-tuned scale, low densities visible — saturate the WP blob); the two Δ kinds = symmetric diverging **linear + `SymLogNorm`** (the symlog panel exposes the low-\|Δn\| wake tail). Slab faces + CAP inner faces dashed. |
| **Density matrix** | z–t **carpets**: total n(t); Δn = n(t)−n(0); per-step Δn = n(t+dt)−n(t); + **wake** = (run − baseline). |
| **Energetics** | total-system E + components (kinetic/Hartree/xc); **WP orbital total** ⟨ψ_WP\|H\|ψ_WP⟩(t) + variance; energy balance/conservation. **ENFORCED (2026-07-11): the ΔE_total(t) plateau-test plot is a 2-panel figure with the total number density N(t)=∫n dV vs time RIGHT BESIDE it** (`delta_total_energy_fig`; N from `electron_number.csv` (classical) else ∫`density_total` VTIs (WP) — its drop is the CAP boundary absorption). **Annotate** the energy-decomposition plot with **CAP + slab-face dashed lines** and the **mean-velocity slab-exit time** (the ballistic time the projectile reaches the far slab face). |
| **Projectile & transport** | COD/centroid z(t),v(t); **stopping** (ΔKE→S, dE_kin vs z); integrated J_z(t); classical track (z,v,F) when classical; **total norm N(t) + boundary-absorption** curve. **Classical: `z(t)` + `KE(z)`** — the KE is **not monotonic** (conservative dip-and-recovery: slows to a minimum at the slab centre, recovers on exit). Mark the **equal-potential window** (symmetric slab faces) — only the net loss between equal-potential points is true stopping; window choice dominates the number. |
| **KS eigen-energy GIF** | bar-chart GIF of the KS eigen-energies E_i(t) over time (`state_energies` / `state_energy_spectra` phases): **two versions — including and excluding the wavepacket state** — plus a **Δ version E_i(t)−E_i(0)**. Example: `…/run_wp_n162_L50_E100_sigma1/results/analysis/observables/ks_energies_delta_no_wp.gif`. |
| **Collective response** | dipole_z(t) + its FFT (**plasmon ringing**); **E-field** (FFT-Poisson) carpet. **FFT-pipeline panel (standard — for EVERY FFT'd signal):** a 3×2 panel via `inqview.visualisation.fourier_panel.fft_pipeline_panel` — row 1 raw \| de-trended, row 2 windowed \| zero-padded, row 3 \|FFT\| linear \| log. The builder auto-emits it for the most dynamic dipole component; plasmon band $\hbar\omega_p$ + the $\Delta\omega=2\pi/\tau$ resolution annotated. Never show a bare spectrum without its pipeline panel. |
| **Momentum** | 1D n(k) total/bath/WP **before·after** + momentum carpet; **momentum-difference heatmaps** Δn(k_z) vs the t=0 baseline over the stored snapshots (which timestep scatters most). *2D k_z–k_⊥ scattering map requires a future 2D-momentum observable* — note, do not fake. |
| **KL metric** | wavepacket momentum **KL-divergence** KL(P_t‖P₀) vs time (`wp_integrity.kl_series`) — drift of the momentum distribution from launch. |
| **KS excitation** | `gs_projected_occupations`; `gamma_transitions` histogram; GS eigenvalue spectrum + occupation dynamics. |
| **Physical anchors & heuristics** (groups A–I) | `inqview.analysis.compute_heuristics`: HEG scales (k_F, v_F, E_F, λ_F, ω_p, **T_plasmon**, k_TF), projectile timescales (slab entry/exit/box-edge), WP zero-point KE, spreading factor, norm/absorption, Lindhard refs. Print the **T_plasmon vs τ resolution flag**. |
| **Loss function** | L(q,ω) via `spectral_weight` — **always produced**, **low-resolution NOTE** when τ short (and flag when τ < T_plasmon). |
| **Integrity / sanity** | WP norm/integrity; N(t) drainage curve. |
| **Optional (opt-in)** | KL-divergence of density; planar-density snapshots; plasmon_spectrum; bath-only energy. |

**Plot examples / spec the user demands** (`docs/notes/plots_examples.md`): energy
decomposition (with CAP + slab dashed lines), system-setup figure, GS KS-excitation
decomposition, momentum 1D + momentum-difference, E-field, KL metric, (Fourier-gated)
2D loss function, **KS eigen-energy bar GIFs (incl/excl WP + Δ)**, and **S(v)/S(E) curves
for a given density** (sweep-level — belongs in the run-SET *study* notebook, not a single
run; e.g. `systems/jellium/hypotheses/06_sigma_convergence/figures/sv_convergence{,_energy}.png`).
LEED screens are coronene-only. Keep that file as the reference for "what these plots should
look like"; every run-notebook reproduces the jellium-applicable subset.

5. **Takeaway** — 2–4 bullets: the headline number(s) for this run and what they imply.

## Carpets, GIF, format (resolved 2026-06-22)

- **z–t carpets are the default** for the density matrix (compact, whole-run-at-a-glance,
  light). **Exactly one xz-slice density GIF** (total density) as the lead visual.
  Static line plots for scalars. Do NOT GIF every panel.
- Shared fixed colour scale for any directly-compared maps/GIF frames
  (report-figures production rule 7). Figures `.png`, canonical theme
  (`inqview.visualisation.style.apply_theme()`). Do not preview generated images.
- **ENFORCED — LOG ALONGSIDE LINEAR (2026-07-11, user).** Every density map, density
  GIF, carpet, and spectrum shows **LINEAR | LOG side by side** (never log-only or
  linear-only): non-negative fields use `LogNorm`; signed/diverging fields (Δn,
  Δ-carpets) use `SymLogNorm` (symmetric log, `linthresh ≈ vmax/100`). This is the
  project "always linear AND log" shared-colorbar rule; a single-scale panel is a
  defect.
- **ENFORCED — N(t) BESIDE ΔE_total (2026-07-11, user).** The total number density
  N(t)=∫n dV vs time sits in the SAME figure as the ΔE_total(t) plateau-test plot
  (`delta_total_energy_fig`, 2-panel). Never ship the energy-change plot without N(t)
  beside it.

## Mechanics

- The builder is **skill-local**: `.claude/skills/run-notebook/run_notebook_builder.py`
  (shippable — nothing in `docs/`). Each run carries a thin `analyse.py` that imports
  and calls it; the heavy logic stays in the skill.
- It **runs `inqview.pipeline.runner.run(results_dir, phases=…)`** first (figures land
  in `results/analysis/…`, irrelevant phases auto-`[skip]`), detects the **run-type**
  from which observables exist (`momentum_distribution.csv`→WP; `electron_track.csv`→
  classical; neither→baseline), then assembles + executes the notebook (`nbformat` +
  `ExecutePreprocessor`, 0 errors), embedding each phase's figures in its battery group.
- Run with the venv + stack on path:
  `PYTHONPATH=…/inq-stack/python …/venv/bin/python3 run_notebook_builder.py <results_dir> <out.ipynb>`.
- **Reader annotations survive rebuilds — two mechanisms (the builder owns neither).**
  Rebuilds regenerate the whole notebook, so:
  1. **Harvest-before-rebuild (in-context):** every builder cell is stamped
     `metadata.gen="builder"` (+ an `anchor` slug for headed markdown). At the start
     of a rebuild the builder **harvests** any markdown a reader added (cells without
     that tag), then **re-injects** each at the same anchor (stamped `gen="user"`,
     so it round-trips every future build). If a section was renamed/removed, its
     annotation lands in a "📌 Carried-over reader annotations" section — never
     dropped. So you may **annotate the `.ipynb` directly in Jupyter** and it persists.
     Functions: `harvest_user_cells()` / `reinject()` / `tag_builder()`.
  2. **Sidecar (deliberate top notes):** `<out_stem>.notes.md` beside the `.ipynb`
     is pinned as a top "📝 Reader notes / TODOs" cell on every build.
  - **Transition:** a pre-tagging notebook (no `gen` tags) is harvested as nothing —
    rebuild each notebook **once** to tag its builder cells *before* annotating.
  - Done = a reader-added markdown cell round-trips a rebuild unchanged, at its anchor.
- **Figures are PATH-REFERENCED, not base64-embedded** (`img()` emits
  `![cap](relpath)`), so the `.ipynb` stays KB-sized and renders instantly; figures
  must travel beside the notebook. Density maps render LINEAR | LOG side by side with
  a fixed colour scale and readable (`×10ⁿ`, ≤2 s.f.) colorbars (report-figures
  rules 7–9).

## Definition of done

- Notebook executes to **0 errors**; outputs embedded; sections 1–4 precede the battery.
- Every battery section either shows real data or is **auto-skipped** (run-type adaptive)
  — never a fabricated panel; the 2D-scattering and (if absent) bath/WP-density panels
  carry the "requires future observable" note.
- Loss function present with a low-resolution note when τ is short.
- Config is fully reconstructable from `run_summary.txt`; every source file linked.
- Builder stays skill-local; generated `.ipynb` written to the run's `hypotheses/` folder.
