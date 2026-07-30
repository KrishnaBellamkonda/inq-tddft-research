# Plan: norm-corrected stopping power from oscillating jellium WP runs

Date: 2026-07-28. Companion to `docs/notes/inq-energy-normalization-error.md`
and the vacuum investigation (`systems/vacuum/hypotheses/cap_norm_investigation/`).

## Goal
Post-process recent localised-jellium WP runs where "energy oscillation" (the
drain-then-rise E_total) was flagged, apply the norm correction, and get a
reasonable estimate of the electronic stopping power S — WITHOUT re-running or
editing source.

## The correction (kinetic-only, WP orbital)
INQ's `energy_kinetic` sums each orbital's kinetic divided by its own norm
(energy.hpp:83,55). Only orbitals that LOSE norm are biased. In these runs the
CAP absorbs essentially only the WP (N_total drops by ~1 = the WP; bath orbitals
keep norm≈1). So correct just the WP contribution:

  E_corr(t) = E_total(t) − occ_WP · e_kin_ha(t) · (1 − norm_WP(t))

- `E_total`      : observables.csv:energy_total (reported, Ha)
- `e_kin_ha`     : wp_momentum_stats.csv (WP per-particle kinetic ⟨T_WP⟩/norm_WP)
- `norm_WP`      : wp_real_space_stats.csv:norm_check  (∫|ψ_WP|²dV, → uses per-step
                   momentum norm_check ratio when the real-space cadence is coarse)
- `occ_WP`       : 1 (from occupations_vs_time; bath orbitals are 2)

Derivation: reported WP kinetic = occ·e_kin_ha (per-particle); extensive WP
kinetic = occ·e_kin_ha·norm_WP; subtract the difference.

## Guards (per stopping-power-extraction skill §0) — MUST pass or S is invalid
1. **N-conservation of the BATH.** The WP is *meant* to be absorbed, so N_total
   dropping by ≈1 is fine; what must hold is that the BATH (N_total − norm_WP) is
   ~constant. If the CAP also drains bath orbitals, more orbitals are biased and
   the single-orbital correction is incomplete → flag.
2. **Energy channel sanity.** Compare E_corr plateau against the bath-energy
   channel where available.

## The governing caveat (why raw ΔE_total is NOT the stopping power)
E_total_extensive = E_bath + E_WP + E_interaction. A CAP-absorbed WP removes its
remaining KE from E_total, so ΔE_total ≈ ΔE_bath − KE_carried_away. The quantity
we want for S is the energy DEPOSITED IN THE BATH, not the total (which also loses
the absorbed WP's energy). Options for the headline S (decide per run):
- **(B-slab) equal-potential slab-face window**: ΔE between the two slab faces
  z=±L_z/2 (reversible mean-field cancels; excludes the vacuum-corridor + CAP
  region) — the cleanest deposit measure for a slab with vacuum corridors.
- **(bath channel)** ΔE_bath directly if a bath-only energy observable exists.
- **raw ΔE_corr/L_z** ONLY when the WP is NOT yet absorbed at the measurement
  time (deposit measured before the WP reaches the CAP).

## Steps
1. Reusable tool `norm_corrected_stopping.py` (hypotheses/energy_oscillation_diagnosis/):
   load a run → E_corr(t), N-guards, plot E_total vs E_corr, and an S estimate
   with the window stated. Path-referenced figures, canonical theme.
2. **Validate on the CLEAN p3_wp** (E_total already decays cleanly): the
   correction must not break it; establishes a known-good S baseline.
3. Apply to the flagged OSCILLATING runs (cap_fix / pbc_open_z / the runs the
   user names) → does the correction flatten the oscillation?
4. Report S per run with the window, guards, and caveat; user owns the verdict.

## Candidate runs (confirm with user)
- clean baseline: `qsp_phase3/wp/results/p3_wp`
- meeting workhorse: `qsp_phase5/.../p5_wp_v1p3`
- flagged-oscillation: cap_fix / pbc_open_z campaign runs (identify dirs)

## Extension (2026-07-29): IN-RUN extensive kinetic observable + vacuum validation

The post-hoc correction above inherits two weaknesses: (a) it assumes bath
norms ≈ 1 (single-orbital correction), (b) it un-divides an ill-conditioned
ratio as norm→0. The in-run fix computes the BARE per-orbital kinetic before
any division exists:

- New observable `inqkit::observables::OrbitalKineticStats`
  (`inq-stack/include/inqkit/observables/orbital_kinetic_stats.hpp`): one
  `to_fourier` of the orbital set per invocation; per orbital i,
  norm_i = Σ_k |ψ̃_i|², T_i = ½ Σ_k k²|ψ̃_i|² (Parseval-exact, same math as
  `ham.kinetic_expectation_value` = `laplacian_expectation_value` at gamma,
  zero vector potential). CSV: kin_bare_total = Σ occ_i·T_i (the extensive
  kinetic), kin_normdiv_total = Σ occ_i·T_i/norm_i (must equal INQ's reported
  `energy_kinetic` — per-step identity check), norm_total, per-orbital
  norm_i/T_i columns, and wall_ms (observable cost).
- No engine edit: `inq/` and `inq-study/` untouched (`occ_sum`'s /norm at
  `energy.hpp:55` stays; SCF relies on the Rayleigh quotient).

Validation test (vacuum, DOUBLE-SIDED CAP — new geometry):
- `systems/vacuum/scripts/wp_traversal_energy/run.cpp` gains `WP_CAP2=1`
  (−z band `absorbing(η, −0.5+w/2, w)` summed with the +z band via
  `perturbations::sum`; `rvector` is contravariant/fractional, matching
  absorbing's mid_pos convention — verified in source), `WP_EXTKIN=1`
  (enables the observable), propagate wall-time in run_summary, and the
  mandatory final checkpoint (rule final-timestep-checkpoint).
- Geometry: LZ=60 box [−30,30], CAP_L=15 BOTH ends (+z [15,30], −z [−30,−15]),
  launch z=0 → 5σ₀=15 Bohr clearance to BOTH CAP inner edges (boundary rule);
  σ₀=3, E=400 eV (k₀σ₀=16.3, transit spread ~5%), h=0.4, dt=0.01, η=−3.5,
  NSTEPS≈700 (traversal 15 Bohr + full CAP transit at v=5.42).
- Runs: `dcap_extkin` (observable ON → both INQ-reported and extensive
  kinetic in one run) and `dcap_baseline` (OFF → per-step timing baseline).
- Acceptance: (1) identity kin_normdiv_total == energies.csv:kinetic to
  solver precision every step; (2) at t=0, kin_bare == kin reported (norm 1);
  (3) E_corr = total − kinetic + kin_bare decays smoothly to ~0 with the norm
  (double-sided: no wrapped remnant channel at all); (4) timing overhead
  reported (per-step % vs baseline; extrapolate to 162-orbital jellium).

## Run design (2026-07-29): extkin_plateau_E100 — first jellium run with the IN-RUN fix

Goal: a clean E_plateau (total electronic energy deposited in the slab) from ONE
CAP run, with the norm-division artifact removed in-run by
`inqkit::observables::OrbitalKineticStats` (corrected extensive total
E_corr = total − kinetic + kin_bare, per-step). Supersedes the post-hoc
`E_ext = E_reported − e_kin_ha·(1−norm)` route used by the replica campaign.

Decision log (user-interviewed 2026-07-29; each value + why):

| Parameter | Value | Justification |
|---|---|---|
| Box | 35×35×120 Bohr, z = traversal | user spec |
| Slab | N=92, r_s=4.0, thickness 20.13 Bohr (half-width 10.067, faces ±10.07), edge_width 0.5 | user: ≥20 Bohr total thickness, N = closest even (n₀·1225·20 = 91.4 → 92); r_s raised 3→4 per user (a 5-Bohr r_s=3 slab would be quasi-2D, ~1 z-subband; at 20 Bohr/r_s=4 ~2–3 subbands, εF=0.115 Ha) |
| States | 46 occ + 16 extra, T≈100 K | user: >10 extra; 35% of occupied, between n52 (38%) and replica n102 (47%) precedents |
| WP | σ=1.5, E=100 eV (k₀=2.711), launch z=−17.5, +z | user chose compact projectile over shape stability KNOWINGLY: free-dispersion ×1.58 at slab entry, ×3.0 mid-slab (τ=σ²=2.25 a.u., k₀σ=4.1, σ_E/E=35%). 5σ=7.5 Bohr to slab face; 27.5 Bohr (18σ) to −z CAP inner edge — 5σ rule met everywhere |
| CAP | two-sided, 15 Bohr each end (inner edges ±45), η=−1.0, sin², ETRS | replica η=−0.7·(20/15) rescaled to W=15; predicted survival ~8e-4 (≤0.1 eV ledger error) from vacuum calibration exp(−1.30·|η|W/v); bath 35 Bohr clear of CAP |
| Grid | h=0.5 (70×70×240 = 1.18M pts) | replica precedent at denser r_s=3.32; cutoff guard: p₀+3σ_p=4.13 ≪ k_Nyq=6.28, E_cut=537 eV > 4·E_kin |
| Time | dt=0.04, N_STEPS=1500 (t=60 a.u.) | user chose speed; H·dt=0.79 < validated 1.78 (n52 family) < 2.2 cliff. t=60 absorbs the slow (−3σ_p, v≈1.0) tail. ⚠ FIRST CAP run at dt=0.04 — absorption quality at this dt UNVERIFIED (user skipped the vacuum gate) |
| Observables | OrbitalKineticStats every step (~62 orbitals; est. ~35 ms/step ≈ 3–7%), energies + wp_momentum every step, density_total+density_wp every 10 steps (150 frames, ~1.4 GB), wavefunction_wp sparse, ckpt every 200 + final | E_plateau from the corrected series; equal-cadence rule for wp/total; /local/data at 99% → lean VTIs |
| Scope | CAP run only — NO no-CAP twin, NO dt gate | user decision (speed); recorded caveats: no in-system conservation control, single-source plateau |

Layout: scripts/extkin_plateau_E100/{gs,wp}/run.cpp + run_extkin_plateau.sh;
config shared/configs/slab_n92_L35x35x120_w0p5.hpp; GS →
shared_gs/slab_n92_L35x35x120_w0p5_h0p5; analysis →
hypotheses/extkin_plateau_E100/ (builder auto-called by dispatcher).

E_plateau extraction: E_corr(t) − E_GS, averaged over the post-absorption
window (last ~10 a.u. once norm_total(t) of the WP orbital is <1e-3 and
dE_corr/dt is flat). Cross-checks INSIDE the run: identity
kin_normdiv_total == energies.csv:kinetic per step; WP-orbital norm from
orbital_kinetic_stats vs wp_momentum_stats.
