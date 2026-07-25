# Handover: energy-decomposition presentation section (meeting 2026-07-16)

Task: produce **figure/image artifacts** (NOT slides — user builds slides) for the
first section of a presentation introducing the **energy-decomposition technique**
on a classical-vs-wavepacket localised-jellium case study.

Output dir: `/local/data/public/skcb2/tddft/docs/reports/meeting-16-07-2026/energy-decomposition/`
(`assets/` = PNGs, build scripts alongside).

## THE case study — one complete p2 twin pair (created 2026-07-16)

After a long search, **no pre-existing p2 run had density frames + INQ components +
pairwise interactions.csv together**. So a fresh complete pair was generated (user
approved a short re-run):

`ResearchProject/systems/localised_jellium/scripts/localised_jellium_dynamics/runs/twin_ec_r12_p2/`
- `wp/`  — phase5_wp binary, σ_WP=0.5, **k0=0 (at rest)**, periodicity 2, Lz=120,
  launch_z=−24.5 (r=12 from slab face at −12.5), N=82, dt=0.01, 4 steps,
  LJ_SAVE_EVERY=1. GS = campaign_autorun/runs/h2/gs_p2_lz120/checkpoint.
- `classical/` — proj_dyn binary (moving Gaussian charge) at **k0=0** → static;
  proj stayed at z=−24.5 (v≈0). Same everything else.
- Both emit: `raw/observables/observables.csv` (E_total/kin/hartree/xc/**external**),
  `raw/observables/interactions.csv` (pairwise e_ss/e_pp/e_ps/e_sb/e_pb/e_bb),
  `frames/total/density_t000000.vti` (t=0). Classical also `projectile.csv`
  (`energy_proj_bg_ideal`). WP also `frames/wp/`.
- **`twin_manifest.json` valid=true** (check_twin.py gate passed).
- Binaries were ALREADY built (no compile). Both GPUs free; ran WP=GPU0, cl=GPU1.
- NOTE the binaries prepend `results/` to an absolute LJ_OUT → outputs first landed
  under `results/local/data/.../wp`; they were flattened to `wp/` and `classical/`.

### Verified numbers (twin_decompose.py on the pair — reproduce the case study)
- dKin_localisation = **81.74 eV** (= 3/(4σ²), at rest) ✓
- dXC = **−16.47 eV** (WP self-XC) ✓
- residual R = d(E_H+E_ext) − U_proj_bg = **20.81 eV** (= WP self-Hartree) ✓
- SIE = R + dXC = **4.34 eV** ✓
- Pairwise gauge test = **0.0000** (E_SS,E_SB,E_BB unchanged → every Δ physical) ✓
- energy conservation 0.0000.

### Table 1 data (components, step 0, eV) — classical / WP / Δ(WP−CL)
E_total 1503.2 / 1724.0 / +220.8 · E_kin 75.0 / 156.7 / **+81.7** ·
E_H −2195.0 / −2314.0 / −119.0 · E_ext 3854.1 / 4128.7 / +274.5 ·
E_H+E_ext 1659.2 / 1814.6 / +155.5 · E_XC −230.9 / −247.3 / **−16.5** ·
E_proj_bg (cl only) **134.7** · residual **20.8**.
**CONVENTION CAVEAT (critical):** the WP cell is net −1 charged, so the Poisson G=0
term makes E_H, E_ext, E_H+E_ext AND E_total individually convention-dependent
(p2 here gives d(H+ext)=+155 vs p3 campaign −121). Only **dKin, dXC, and the
residual (= self-Hartree) are convention-invariant.** Present raw components but
highlight the residual as the physical invariant — matches the user's note
`docs/notes/energy-decomposition-skill.md`.

### Table 2 data (pairwise S/P/B, step 0, eV) — cl / WP / Δ
E_SS −2195.0/−2195.0/0.00 · E_PP 20.82/20.80/−0.01 · E_PS −139.87/−139.82/+0.04 ·
E_SB 3994.0/3994.0/0.00 · E_PB 134.69/134.65/−0.04 · E_BB −1803.0/−1803.0/0.00.
S=slab electrons, P=projectile, B=+background. Gauge test 0.

## Deliverables status
- [x] `assets/intro_density_classical.png`, `intro_density_wp.png` — t=0 total density
  xz slice from the pair; shared LOG colorbar, slab faces dashed, red dashed proj line
  at z=−24.5 annotated. Build: `build_intro_density.py`. (User previews, not me.)
- [x] Table 1 image (`table1_components.png`) + symbol PNGs (`symbol_E_{total,H,XC,ext,proj_bg}.png`)
  + residual-formula PNG (`formula_residual.png`). ‡-flagged convention-dependent rows.
- [x] Table 2 image (`table2_pairwise.png`, S/P/B + gauge-test note) + convention legend
  (`legend_spb.png`).
Build script for tables/equations: `build_tables_and_equations.py` (DONE, mathtext validated).

## r-independence proof (2026-07-16)
- Radius sweep `runs/twin_ec_rsweep/` (orchestrator `run_r_sweep.py`): at-rest twin
  pairs (Gaussian-charge classical + WP, σ=0.5, p2) at r∈{4,12,20,28,36,40}, 2 steps,
  no frames. **Gotcha:** phase5_wp does `step % SAVE_EVERY` WITHOUT a `>0` guard →
  `LJ_SAVE_EVERY=0` is div-by-zero SIGFPE (rc=−8). MUST use `LJ_SAVE_EVERY>=1`.
- Result `sweep_R_SIE.csv`: **R = 20.80–20.81 eV (flat, ±0.01), SIE = 4.33–4.34 eV
  (flat), dKin 81.7, dXC −16.5** across all r → residual/SIE are r-INDEPENDENT
  (intrinsic WP self-energy, not a projectile-slab interaction). Plot:
  `assets/residual_sie_vs_r.png` (build `build_r_independence.py`).
- [x] `assets/calc_sie.png` — SIE = R + ΔE_XC = 20.8 − 16.5 = 4.3 eV (usetex).
- [x] `assets/calc_self_hartree_analytic.png` — E_H^self = 1/(2 σ_ρ √π),
  σ_ρ=σ_WP/√2=0.354 → 0.798 Ha = 21.7 eV (usetex). Numeric self-Hartree (twin) 20.8;
  analytic 21.7; ~0.9 eV shortfall = open-z gauge (per twin_decompose note).

## Pseudopotential-vs-perturbation motivation (2026-07-16)
- [x] `assets/table_pseudopot_vs_perturbation.png` — why the Gaussian perturbation
  beats the UPF-ghost pseudopotential. VERIFIED vs `perturbation_method_study.ipynb`
  (cells 4/11/29/30) + user note. KEY: the two projectile potentials are POINTWISE
  IDENTICAL (V_loc(r)=erf(r/σ)/r, diff ~1e-16); the difference is INQ's *treatment* —
  UPF ghost (Z_val=0) puts the erf/r tail as a truncated local potential that aliases
  with r_cut → residual **7.4 eV** (grid sign-swings). Perturbation adds
  v=+poisson(n_proj) exactly → residual **20.8 eV** = WP self-Hartree (matches analytic
  21.7 within ~0.9 eV gauge AND INQ boundary-matched Poisson self-Hartree to ~0.01 eV).
  So the perturbation recovers WP−WP self-Hartree correctly; the pseudopotential does not.
- User rejected the TABLE form; wants a PLOT (writes differences in caption himself).
  → `assets/decomposition_comparison.png` (build `build_decomposition_comparison.py`):
  grouped bars of dKin/dXC/R/SIE for pseudopotential vs perturbation. dKin(81.7) &
  dXC(−16.5) identical (WP-only); R diverges 7.4 vs 20.8; SIE −9.1 (unphysical) vs
  +4.3; dashed analytic self-Hartree marker (21.7) over the R group. Table PNG deleted
  (`build_pseudopot_vs_perturbation.py` script kept for provenance).

## Density-leak / spill-out plots (2026-07-16)
- [x] Four deck plots from the semiempirical_spillout experiments (build
  `build_spillout_plots.py`, data `scripts/semiempirical_spillout/runs/`):
  `spillout_field.png` (E(z) non-zero far-field plateau), `spillout_charge.png`
  (enclosed net charge / % electrons outside), `spillout_lz_sweep.png` (Lz sweep:
  edge pile-up grows with box, interior tail Lz-flat), `spillout_w_sweep.png`
  (w sweep: w≥1 removes the pile-up).
- Reproduction VERIFIED bit-for-bit against `semiempirical_spillout.ipynb` canonical
  output. Numbers (lz160): 5.98% electrons beyond slab face (= the note's "~6%");
  Q(|z|<25)=1.43 e, E@25=0.098 eV/Bohr, Emax_vac=0.30. **The note's 0.39 e / 0.0268
  eV/Bohr was the earlier THIN-PLATE estimate (superseded) — refined full-density
  gives 1.43 e.** w=1 → edge 5e-15, Q25=0. Lz: Q25 1.29→1.43→1.55 (90→160→240),
  tail@20 Lz-flat 1.3e-6 → boundary artefact. es60 vs es20: ΔQ25 ~0.1 (not empty-states).
- p3 (PBC) runs exist (p3_lz90/160/240) — a p2-vs-p3 open-z-vs-PBC plot is available
  if the user wants the "does it disappear under PBC?" question answered.

## Next candidate steps (not yet requested)
- Workflow diagram of the decomposition method (Graphviz .dot per scientific-figures §7).
- Per-run / twin notebook for the pair (notebook-density-gif rule) if a deep-dive wanted.
- If user wants the p3 (offset-free) component numbers instead of p2 for Table 1, the
  h0_p3 runs give them (individual components still convention-dependent — see caveat).

## Rules in force
- Figures: canonical theme (`inqview.visualisation.style`), 600 dpi, PNG only,
  transparent+tight for drop-in. NEVER preview PNGs myself (user previews).
- 2 s.f. default (3 s.f. for near-equal); carry units.
- User makes slides; I only produce plots/equations/images/workflow diagrams.
- σ = σ_WP always; classical charge std σ_pot = σ_WP/√2 (methods footnote only).
