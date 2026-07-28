# Plan: qsp_phase5 momentum-KE-loss stopping notebook (orbital-free)

**Goal.** A study notebook `hypotheses/qsp_phase5/qsp_phase5_momentum_stopping.ipynb`
that extracts a *second* stopping-power estimate S_drift(v) from the qsp_phase5 WP
velocity sweep (k0 = 1.3, 3, 4, 5, 6), using only density-level (KS-label-free)
observables: the momentum-dependent kinetic-energy loss of the projectile, measured
by time-of-flight (TOF) detector planes in the vacuum corridors. Compares against
the existing deposit-based S (results_*.json, believed too high due to
localisation-energy + capture deposition).

## Method (locked in conversation 2026-07-27)

- z-profiles ρ(z,t) = ∫dxdy [n_total(t) − n_gs] from `density_total` VTIs + the
  one-off `density_gs_system` VTI. (`density_delta` is n(t)−n(0), NOT n−n_gs —
  verified numerically.)
- Longitudinal flux via 1D continuity with CAP sink:
  J(z,t) = ∫_z^{+L/2} ∂ρ/∂t dz' + A_+(t), A_± = ∫_CAP± 2W(z)ρ dz,
  W(z) = |η| sin²(π(|z|−35)/10) on 35<|z|<45, η = −0.7 Ha (inq-study
  perturbations::absorbing; sink REGION verified empirically from band-resolved
  norm decay; W shape to be validated in-notebook from late-time in-CAP decay).
- TOF detector planes at z_m = ∓15.5 Bohr (slab edge 12.5 + 3 buffer):
  crossing flux J(z_m,t), passing velocity v(t) = J/ρ, per-plane accumulations
  N = ∫J dt, P = ∫ vJ dt·? (momentum flux ∫ρv² dt), KE_z = ∫ ½v²J dt.
  Entrance plane separates incoming (J>0) from reflected (J<0) by sign/time.
- S_drift(v0) = [KE_in/N_in − KE_out/N_T] / L_slab (L=25 Bohr), transmitted
  channel; reflected + captured channels reported separately.
- Shape (localisation) term T_W = ∫|∇n|²/(8n) over vacuum lobes per frame
  (3D pass, same read as profiles). Flow-vs-shape narrative on the entrance leg.
- Cross-checks: wp_momentum_stats.csv (KS-orbital <p>, <p²>, e_kin) — measures
  "how KS-orbital-dependent" the orbital route is; momentum_distribution.csv for
  launch aliasing QC; energies ledger closure.

## Known data caveats (QC section, mandatory)

- Momentum aliasing at launch (dx=0.5 → k_Ny=6.28; σ_p=1.41): measured <p_z>(0)
  = 1.30/2.86/3.30/2.63/0.59 vs k0 = 1.3/3/4/5/6 → v4p0 corrupted, v5p0 severe,
  v6p0 catastrophic. T_drift(0) must be MEASURED (entrance TOF), never k0²/2.
  Per-point trust grades on the final S(v).
- CAP mid/width passed as fractions (40/90, 10/90) to absorbing.hpp which
  compares vs Bohr — empirically the sink IS at |z|∈[35,45]; validate W(z) shape
  from data before using the sink correction.

## Steps

1. [x] Recon: data inventory, aliasing discovery, CAP region empirically located.
2. [x] `qsp5_momentum_kinematics.py` (hypotheses/qsp_phase5/): one 3D pass per run
   → cache npz (z-profiles, lobe T_W, lobe N/∫zρ per frame). Parallel over runs.
3. [x] `build_momentum_stopping_report.py`: nbformat builder per notebook-making
   skill. METHOD UPGRADED during validation (see handover): side-adaptive flux
   integration + exceedance-matched (rank-transport) S(u) with a free-packet
   null test calibrating the trusted rank window (syst ±0.41 eV/Bohr).
4. [x] Executed to 0 errors (40 cells); S(u) curves + headline table + comparison
   plot vs deposit-based S and classical/Lindhard. Headlines in
   `momentum_stopping_summary.json`.
5. [x] Handover: docs/handovers/qsp5-momentum-stopping.md.

## Files

- Runs: ResearchProject/systems/localised_jellium/scripts/qsp_phase5/wp/results/p5_wp_v{1p3,3p0,4p0,5p0,6p0}/
- Analysis home: ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/
- Prior deposit-based S: results_p5_wp_*.json (S_eVbohr, deposited_eV)
- References: classical_sigma0p5_bulk.csv, lindhard_ref.npz
