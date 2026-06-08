# Plan: matched WP-vs-classical induced-wake case study (σ1 @ E100)

Status: agreed via grill 2026-06-01. Supersedes the batch2 wake figures
(user is deleting `docs/presentations/storyline/tasks/batch2_figures/`).

## Motivation / what was wrong with batch2

The batch2 WP−classical wake figures paired runs that are NOT physically
matched:

- σ-sweep WP runs launch at `boundary + 4σ` → z = −23/−21/−13/−1 Bohr for
  σ = 0.5/1/3/8, but were all differenced against `classical_v2` at z = −21.
- dt mismatches: σ0.5/3/8 `_wf` WP are dt=0.02; classical_v2 is dt=0.01.
  E20/E25 WP are dt=0.01; their classical companions are dt=0.02.
- Different end-times and frame cadences → sampling the classical run at WP
  physical times aliases ("classical Δn jumps around between timesteps").

Frame ORDERING was already correct (driver sorts WP times ascending); the
defect is run-pairing, not ordering, and not missing wp-density.

Wavepacket density IS saved (verified 2026-06-01):
- `_wf` runs: full cadence (density_wp == density_total frame count).
- `_v2` runs: saved but ~10× COARSER than density_total (e.g. σ1_v2: 32 wp
  vs 317 total) → only sparse exact subtraction possible.
- original (non-v2/non-wf) runs: 0 wp frames → unusable for wake.

## Chosen basis (only existing fully-matched pair)

| | run_dir | dt | launch z | tot last t | wp last t | wp frames |
|---|---|---|---|---|---|---|
| WP | `run_wp_n162_L50_E100_sigma1_v2` | 0.01 | −21 | 9.48 | 9.30 | 32 (all exact) |
| classical | `run_classical_n162_L50_E100_v2` | 0.01 | −21 | 16.56 | — | — |

Shared: dt, launch z, L=50 box, E=100 eV, σ=1 WP. Differences handled:
truncate classical to WP's last wp-frame time (~9.30 a.u.); sample both at
the 32 exact wp-frame times so `n_total − n_wp` is exact (no moving-WP dipole
residual).

## Deliverables (new folder `tasks/wp_vs_classical_matched/`)

Bath density everywhere: `n_system = n_total − n_wp` (classical: = n_total);
induced `Δn = n_system(t) − n_system(t0)`. Shared colorbar on WP & classical
panels; own colorbar on the difference. Linear AND symlog.

1. **Corrected GIF** `wake_2d.gif` / `wake_2d_log.gif`: 3-panel xz Δn
   [WP | classical | WP−classical], animated over the 32 matched frames in
   true time order, WP-centroid line, fixed scale across frames.
2. **1D GIF** `wake_1d.gif`: z-profile Δn(z,t) WP & classical overlaid +
   (WP−classical) subpanel + centroid.
3. **Report-standard static panels** (inqview.report1, usetex, 600 DPI, no
   titles): `fig_wake_sigma1E100_2d_{wp,classical,diff}.png` (+ `_log`),
   `fig_wake_sigma1E100_1d.png` (legend lower-right), plus metric panels.
4. **Quantitative metrics** (the "concrete message"):
   - wake peak |Δn| vs t (WP and classical);
   - trailing-oscillation wavelength behind the projectile (FFT of the
     z-profile tail behind the centroid);
   - depletion-vs-enhancement integral (∫ negative vs positive Δn behind
     centroid) vs t;
   - WP vs classical: ratio of integrated |wake| (quantifies the genuinely
     weaker, more diffuse WP wake vs the sharp classical wake).
   Output `metrics.csv` + `fig_wake_metrics.png` + a short `REPORT.md`.

Known-case checks (dev-feedback-loop): Δn(t0)==0 both runs; ∫n_system dV =
162 at t0/mid/late; centroid monotonic.

## Skill + record changes

- `.claude/skills/tddft-simulations/SKILL.md` §3b: add `density_wp` to the
  Jellium-WP Tier-2 compulsory list, MANDATING equal cadence
  (`wf_write_every == write_every`) for any run intended for wake /
  density-difference analysis. Note the sparse-`_v2` failure mode.
- `docs/handovers/wp_vs_classical_matched_wake.md`: corrected todo — original
  runs have 0 wp frames (unusable); `_v2` runs too coarse (sparse-only);
  future wake runs must save density_wp at full cadence.

## Not doing now (user declined)

- New classical twins for the _wf σ runs; fresh twin pair; σ-sweep /
  energy-sweep WP−classical (unmatchable from existing runs).
