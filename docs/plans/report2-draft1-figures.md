# Plan: report-2 draft-1 figures

Status: **active** (opened 2026-08-03).
Style + landmines: `docs/reports/report2/drafts/draft1/CLAUDE.md` (auto-loaded).
Output: `docs/reports/report2/drafts/draft1/figures/`.
Iteration log: `docs/reports/report2/drafts/draft1/plots_draft1_log.md`.

## 0. Narrative (the user's Draft 1, as given — figures serve this)

1. **Bulk jellium** — the wavepacket's stopping power comes out *smaller* than the
   classical twin's. Two candidate causes: (a) orthogonalisation clears bath
   density near the packet, so E_PS is weaker than classical; (b) the uncancelled
   self-Hartree drives the packet's expansion.
2. **Cylindrical jellium (weak coupling)** — a geometry with no orthogonalisation
   clearing and a wide packet, to see whether the KS-orbital definitions behave.
   `S_drift` tracks the classical result; the extra channels the quantum
   projectile carries are identified. Naming: *surface friction* rather than
   stopping power.
3. **Proximity ladder** — close the rings onto the packet, walking weak coupling
   → full contact, and test whether the definition stays well behaved.
4. **Jellium slab** — abandon the KS-orbital-dependent definition; use a total-
   energy deposit definition instead. Both halves plateau once the projectile is
   far away. Sweep σ_WP and v; a concentrated packet spreads and ends wider.
5. **2D materials** — to do.
6. **Nazarov–Gross verification** — to do.

## 1. Figure inventory

Sizes are the first-cut defaults (1C = 3.5×3.0 in, 2C = 7.0×3.0 in). `S_*`
estimator symbols are defined in the draft CLAUDE.md §2.3.

### Chapter 1 — bulk jellium (figure list per user spec, 2026-08-03)

All in `figures/bulk_jellium/`. Canonical run: `bulk_ks_stopping_rs4`
(r_s = 3.99, σ_WP = 2 Bohr, 100 eV) unless stated.

| id | file | message it serves | status |
|---|---|---|---|
| B-setup | `bulk_setup_density.png` | **the setup** — t = 0 total density, xz mid-y slice, WP case. Classical launch point mentioned verbally, not drawn. | **BLOCKED** — VTIs pruned |
| B0 | `bulk_total_energy.png` | run soundness; referenced verbally ("total energy is conserved") rather than shown in the main text | done |
| B1 | `bulk_ke_vs_path.png` | **"we can do this"** — Δ`T_orb`, Δ`T_drift` (WP) + Δ½mv² (classical) vs path, three curves | done |
| B1b | `bulk_stopping_fit.png` | **"WP stopping ≪ classical"** — B1 + fit window + regression lines + S annotated | done |
| B-pauli | `bulk_bath_density_radial.png` | **orthogonalisation cause** — bath density vs radius from the projectile at t = 0, classical vs WP, plus the difference | **BLOCKED** — VTIs pruned |
| B3 | `bulk_S_vs_sigma.png` | supporting: S vs σ_WP ∈ {1,2,3}, both halves | done |
| B2 | `bulk_interaction_ledger.png` | supporting: Δ E_SS / E_PS / E_PP, classical vs WP | done |

**Self-Hartree subsection (own subsection, user decision 2026-08-03)** — figures
still to be specified; see §2 for the data blocker.

### Chapter 2 — cylindrical jellium, weak coupling (rung r10)

| id | file | content | size | data |
|---|---|---|---|---|
| C1 | `cyl_ke_vs_path.png` | `T_drift` WP vs classical KE vs path; ratio 0.80 annotated. | 1C | `hypotheses/channeling_twin` |
| C2 | `cyl_impulse_deficit.png` | cumulative impulse on WP / on classical vs time (WP receives 76 %). The mechanism. | 1C | channeling twin observables |
| C3 | `cyl_varp_channel.png` | `T_int(t)` for WP-in-medium vs WP-in-vacuum vs classical (≡ 0). **Answers the user's var(p) question visually.** | 1C | channeling twin + free-WP control |
| C4 | `cyl_sic_control.png` | var(p_z) growth: LDA vs SIC-PZ vs free. SIC removes only 21 % of the excess. | 1C | `hypotheses/channeling_sic` |

### Chapter 3 — proximity ladder

| id | file | content | size | data |
|---|---|---|---|---|
| L1 | `ladder_ratio_vs_coupling.png` | `S_wp/S_cl` vs `f_bore(t=0)` (and R_in/σ on a top axis) across r10→r00, with the bulk value 0.18 as a limit marker. **The bridge figure.** | 1C | `proximity_ladder/figures/ladder_summary.csv` |
| L2 | `ladder_S_absolute.png` | `S_wp` and `S_cl` separately vs rung — both rise, WP rises slower. | 1C | same |
| L3 | `ladder_energy_split.png` | Δ`T_drift`, Δ`T_int` (WP) and Δ KE (classical) per rung — where the energy goes as coupling grows. | 1C | same |

### Chapter 4 — jellium slab

| id | file | content | size | data |
|---|---|---|---|---|
| S1 | `slab_energy_plateau.png` | `E_tot(t) − E_GS` for a classical/WP pair, showing the plateau the definition rests on. | 1C | `sigma56_sv` σ=6 v=3.0 pair |
| S2 | `slab_eps_tail_artefact.png` | classical `E_PS(t)` vs the bare monopole `N_e/z`, and raw vs corrected S. **The methods figure that justifies the correction.** | 1C | `sigma56_sv` |
| S3 | `slab_sv_by_sigma.png` | `S_dep(v)` for the WP at σ ∈ {0.5, 2, 3, 6}. Supersedes `sv_eabsorbed_cap.png`. | 1C | `wp_highdensity_sv`, `sigma56_sv` |
| S4 | `slab_S_vs_sigma_twin.png` | `S_dep` vs σ_WP for **classical and WP together** at fixed v. The divergence figure — and it makes the missing classical σ = 2, 3 cells visible as gaps. | 1C | all slab sweeps |
| S5 | `slab_spread_vs_ratio.png` | in-transit spreading factor σ_d(exit)/σ_d(entry) vs `S_wp/S_cl`. σ=6 spreads ×1.12 and still differs by 1.9–3.4×. | 1C | `s56_S_summary.csv`, `sigma_sweep_S_deposit.csv` |

### Chapter 5 — synthesis

| id | file | content | size | data |
|---|---|---|---|---|
| X1 | `estimator_summary.png` | `S_wp/S_cl` for every system on one axis: bulk 0.18, ladder 0.48–0.80, slab 1.9–3.4. The one figure that states the thesis. | 2C | all summary CSVs |

## 2. Known gaps this plan exposes

### 2.0 RESOLVED 2026-08-03 — density fields regenerated, option C already existed

- **A + B done.** `shared/bin/run-bulk-t0.slurm`, array job **32702053_[0-5]**,
  all six COMPLETED in ~80 s each. Outputs in
  `ResearchProject/systems/jellium/scripts/bulk_t0_density/s{1,2,3}_{wp,classical}/`.
  Unmodified production binaries, 2 steps, run from a scratch CWD so the
  production `observables/` could not be clobbered.
- **C was already done** — `systems/vacuum/scripts/wp_selfinteraction/results/`
  holds 24 completed vacuum runs (σ ∈ {1,2,3,4,6,8} × 4 theories) with a
  scale-exact protocol, plus the analysis module
  `systems/vacuum/hypotheses/wp_selfinteraction/sigma_sweep.py`. No new runs.
- σ = 1 bulk `interactions.csv` is still absent, so E_PP(σ) from the *bulk* runs
  still has only two points. The vacuum sweep covers the σ-scaling instead.

The original diagnosis, kept for the record:

### 2.0a BLOCKER (historical) — the bulk density fields were pruned

`bulk_ks_stopping_rs4/{wp,classical}/run.cpp` *does* write
`results/raw/vti/density_{total,wp,delta}` at `Cfg::WRITE_EVERY = 8`, and the
GS density at step 0. **None of it is on disk** — `results/raw/` contains only
`observables/`. The same is true of every `run_free_wp_*` directory (run.cpp and
analyse.py only). So three requested items cannot be built from existing output:

- **B-setup** (t = 0 total density xz slice)
- **B-pauli** (bath density vs radius at t = 0, classical vs WP + difference)
- **self-Hartree in vacuum + its σ-scaling** (needs the free-WP controls)

Surviving field data in the whole jellium tree: two ground-state densities,
`save_gs/gs_L40x40x80_orth_N482_dx0p50/results/density_gs_system` (matches the
rs4 pair) and `save_gs/gs_L46x46x80_orth_N218_dx0p40`. These give the bath before
the projectile exists — useful as a background, useless for the t = 0 comparison,
because the whole point of B-pauli is what orthogonalising the packet does to
that bath.

**Cheapest fix — a short density-only re-run.** Both halves already resume from
`save_gs`, so a run with `BKS_N_STEPS` = a few steps and `WRITE_EVERY = 1` costs
a GS load + WP injection + 1–2 propagation steps and regenerates the t = 0
fields. Options to put to the user:

| option | what it unblocks | rough cost |
|---|---|---|
| A: rs4 wp + classical, ~2 steps, WRITE_EVERY 1 | B-setup, B-pauli | minutes on one GPU |
| B: A, plus σ = 1 and σ = 3 pairs | above, at three widths | ~3× A |
| C: free-WP vacuum runs at σ = 1, 2, 3 | self-Hartree magnitude + σ-scaling in vacuum | small (no bath) |

Note the σ = 1 bulk family also lacks `interactions.csv` entirely (predates the
retrofit), so E_PP(σ) currently has only two points, σ = 2 and 3.

## 2.1 Other gaps

- **classical slab σ = 2 and σ = 3 with a valid estimator** — the missing cells in
  S4. Existing σ=2 classical wrap twins are CAP-free multi-crossing, so
  `E_abs/L` is invalid there; needs an initial-drag extraction or a CAP re-run.
- **σ = 5 slab** — re-running; add to S3/S4/S5 when complete.
- **free-WP vacuum control matched to the cylinder run** — needed for C3/C4;
  confirm which `run_free_wp_*` matches σ and box.
- Uncertainties are quoted for bulk (`stopping_power.txt`, syst-dominated) and
  cylinder (OLS σ in `refined_stopping_summary.csv`) but **not** for the slab
  deposit values — candidate proxy is `plateau_drift_eV / L_slab`.

## 3. Order of work

B1 → B2 → L1 → C1–C3 → S1–S2 → S3–S5 → C4 → B3 → L2–L3 → X1.
Rationale: B1 and L1 are the two figures that carry the most narrative weight and
are cheapest to make (summary CSVs already exist for L1). X1 is last because it
depends on every other number being settled.
