# Plan: Screening / WP-potential test + classical-vs-WP energy deconstruction

Status: **DONE 2026-07-08** (notebook 13 figs / 0 errors; 2/2 CPU runs; see handover
`campaign-autorun-review-organisation.md` § 2026-07-08). Autonomous CPU execution.
Notebook to extend: `ResearchProject/systems/localised_jellium/hypotheses/campaign_autorun_study/theoretical_slab_model.ipynb`.

## Origin (user's notes, `docs/notes/localised-jellium-parameter-study.md` l.109–112)

- **Learning #1** — take a *comparable* classical and WP run with the full energy
  decomposition; look for the **E_xc** difference; then a holistic **total-energy
  difference** (only U_ion-ion "missing" from the classical store); check whether
  that difference agrees with the **analytical infinite-slab model**.
- **Learning #2 (special importance)** — study the **Coulomb potential produced by
  the wavepacket** vs the classical Gaussian radial potential. "Is there something
  we can do to answer these questions?"

User instruction: *deconstruct all energies correctly, then use simple arithmetic
to deduce; run autonomously on CPU; extend the notebook we've been using.*

## Physics framing (neutral — user owns verdicts)

- Two candidate reasons the WP total is flat while the classical decays: (a) genuine
  physics — screening/correlation lowering the WP's Coulomb impact; (b) bookkeeping /
  the classical UPF cutoff artifact (already shown for the classical decay).
- The **classical ghost** enters only as an external potential (no XC, no self-energy,
  E_ion=0). The **WP** is a real electron → it contributes kinetic, Hartree self-energy,
  and **exchange-correlation with the slab**. So E_xc(WP) − E_xc(classical) isolates the
  real-electron XC signature; the Hartree/external split isolates the electrostatics.
- Potential equivalence reduces to **source equivalence**: if n_WP == classical Gaussian
  charge (std σ_ρ = σ_WP/√2 = 0.354), the two projectiles are electrostatically identical
  under any solver. The only distortion is WP orthogonalisation against occupied slab states.

## Deliverables

### New CPU run — S1 (save densities)
- Modify `scripts/campaign_autorun/wp/run.cpp`: env-gated `LJ_SAVE_DENSITY=1` block that,
  at t=0 (after WP injection, before propagate), writes VTIs:
  `density_wp` (= |ψ_WP|² via `density::orbital(el, wp_idx)`),
  `density_total` (= `density::total`), `density_bath` (= `total_excluding_orbital`).
- Driver `scripts/campaign_autorun/screening_wp_test.py`: run WP insertion at **r=12**
  (clean, far from slab → intrinsic orthogonalisation distortion) and **r=4** (near the
  surface → slab-proximity distortion), p2, off the Lz=120 GS, 1 step, CPU.
  Output: `runs/screening_wp/wp_r{4,12}_p2/`.

### Analysis (Python, into the notebook)
- **B1** n_WP radial profile vs ideal Gaussian(σ_ρ=0.354); ∫n_WP; RMS/max deviation. r=4 vs r=12.
- **B2** Coulomb potential comparison: FFT-Poisson of n_WP vs the analytic classical-ghost
  Gaussian potential erf(r/(√2 σ_ρ))/r (= the UPF potential, already verified to RMS 0.000 Ha).
  Validate the Python Poisson against the analytic Gaussian first.
- **B3** Δn_bath(t=0) = n_bath − n_GS ≈ 0 (screening is dynamical, timescale ~T_plasmon ≈
  4900 steps ≫ CPU budget) — flag the dynamical-screening run as GPU/future work, honestly.
- **A1** full energy-component overlay WP vs classical at matched r (existing h0_p2 runs).
- **A2** ΔE_xc(r) = E_xc(WP) − E_xc(classical) — the XC signature.
- **A3** holistic ΔE_tot(r) with correct bookkeeping (classical omits ∫v_ghost·n₊); compare
  to the analytical infinite-slab model prediction.
- **A4** per-r arithmetic ledger: where the WP−classical difference lives (kinetic/Hartree/ext/xc).

## Validation
- C++ change: confirm the new VTIs load via `inqview.load_vti` (physical order, ∫n_WP≈1).
- Python Poisson validated against the analytic Gaussian potential before use on n_WP.
- Notebook executes to 0 errors; neutral framing; figures path-referenced; no `matplotlib.use('Agg')`.

## Rules honoured
- σ = σ_WP labels (σ_ρ = σ_WP/√2 only in the potential-generation footnote).
- Never fftshift a VTI for display; FFT-Poisson handles ordering internally and is validated.
- 2 s.f. reporting; GPU-default waived (user explicitly requested CPU).
