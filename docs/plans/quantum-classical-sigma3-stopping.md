# Plan: Quantum-vs-classical stopping at σ_wp=3, 300 eV (CAP and no-CAP twins)

Goal: isolate the **purely quantum** component of electronic stopping in interacting
LDA jellium by comparing a Gaussian **wavepacket (WP)** projectile to a **classical
Gaussian-charge** projectile at the *same* (unified) σ, with and without a complex
absorbing potential (CAP). σ-convention is unified (WP is truth; see CONTEXT.md):
σ_wp=3 ⇒ charge/density std = σ_wp/√2 = 2.121; matched classical UPF =
`electron_gaussian_wpsigma3p0.upf`.

Shared config (both studies): 50³ cubic box, N=162, r_s=5.69, LDA, reuse GS
`checkpoints/gs_L50_cubic_N162_dx0p40`, dx=0.40, **dt=0.02**, ETRS. Energy **300 eV**
(k0=v0=4.696 a.u.). Full observable suite (run.cpp `cap_baselines`): density VTI
(every 5 steps), observables.csv (E components, current, dipole, density L2),
state_energies, occupations, density_delta, electron_number; **WP runs** add
WP-resolved momentum_distribution; **classical runs** add the projectile track.
GPU: one free (GPU 1); runs are sequential; emails per run + on all-done.

PROVISIONAL until the inq-study engine regression (Task #7).

---

## Study A — WITH CAP (LAUNCHED 2026-06-21)

Two-sided sin² CAP, L=20 total (10 Bohr/side), **η=−1.0** (slabs |z|∈[15,25], free
|z|<15; reflection <10% confirmed from `twosided_cap_vs_mask_study.ipynb`). Launch
**z0=−6.5** (4σ_density inside the −z CAP edge). **N_STEPS=336** (τ=6.72 a.u. = the
rigid-particle full-exit time, launch→far CAP outer edge 31.5 Bohr).

| run | mode | purpose |
|---|---|---|
| R4 | b1 | CAP on, no projectile — bath-drainage / wake reference |
| R1 | b3 | WP σ_wp=3 — quantum stopping (coherent-peak momentum) |
| R2 | b2 | matched classical (charge std 2.121) — classical stopping (−dKE/dz) |
| R3 | (vac) | **DEFERRED** — vacuum-WP SIE control (single electron, no bath); needs new code (periodic-background subtlety). Bounds the ~7 eV SIE; mandatory before any *quantitative* quantum-component claim. |

Dispatcher: `scripts/cap_baselines/qvc_dispatch_phase1.sh`; outputs
`results/qvc_{b1_s3,wp_s3_E300,cl_s3_E300}/`. Notebook: `hypotheses/qvc_cap_sigma3/`.

---

## Study B — NO CAP (this plan; the no-CAP twin)

**Exactly Study A's config, CAP removed** (CAP_ETA=0 ⇒ no-op absorber). The projectile
now crosses the **periodic** box and wraps instead of being absorbed — so we relaunch
**closer to the boundary** for a longer clean traversal, and lengthen τ accordingly.
CAP is the *only* substantive difference, so A−B isolates the absorber's effect and B
gives a clean periodic single-pass stopping.

### Geometry / timing (RECOMMENDED — to confirm)
- **Launch z0 = −18** Bohr: the furthest launch that still keeps the WP fully inside at
  t=0 (3σ_density = 6.36 ⇒ tail at −24.4, just inside the −25 edge; no initial wrap).
  Gives a **clean traversal of ~36 Bohr** (centroid −18 → +18) vs 21.5 with the CAP —
  the "more traversal time" requested.
- **N_STEPS = 532** (τ ≈ 10.65 a.u. = one full box period, 50 Bohr at v=4.696) — covers
  the full clean crossing plus the onset of the periodic wrap. (Adjustable: 2 periods
  ≈ 1064 steps if periodic re-crossing is wanted.)
- WRITE_EVERY = 5 (same dense cadence; ~106 momentum dumps over the run).

### Runs (no CAP, same binary, CAP_ETA=0)
| run | mode | subdir | purpose |
|---|---|---|---|
| BW | b3 | qvc_nocap_wp_s3_E300 | WP σ_wp=3, no CAP — periodic quantum stopping |
| BC | b2 | qvc_nocap_cl_s3_E300 | matched classical, no CAP — periodic classical stopping |

No-projectile baseline is unnecessary without absorption (N is conserved; the static
bath is a trivial reference). Sanity check: electron_number stays flat (no drainage).

### Analysis (same as Study A)
Coherent-peak WP momentum loss (BW) and ion −dKE/dz (BC) over the clean traversal;
energetics, current, density/wake GIFs, momentum n(k) before/after, centroid track —
all with the unified-σ labels. Notebook: `hypotheses/qvc_nocap_sigma3/` with
`build_qvc_nocap_sigma3_report.py` → `qvc_nocap_sigma3_study.ipynb` (house narrative).

### Optional extension (related, to confirm)
Classical S(v) **benchmark vs analytical Lindhard** at the matched width (charge std
2.121), pure jellium, no CAP, at 150/300/450/600 eV (reusing `06_sigma_convergence`
extraction) — to see how the classical S(v) shape deviates from Lindhard at this σ.

---

## Open questions (for the user)
1. Study B launch z0 = −18 and τ = one box period (532 steps), or a different
   launch / more periods?
2. Include the multi-energy classical S(v) Lindhard benchmark (150–600 eV) in this
   batch, or keep Study B to the 300 eV twin only?
3. Should Study B auto-run after Study A Phase 1 (on GPU 1), or wait for your review
   of Study A results first?
4. R3 vacuum-WP SIE control — build it now (deferred from Study A), or after Study B?

## Test / validation status
- `gaussian_psp` unification: tests updated, 8/8 pass.
- Matched UPF `electron_gaussian_wpsigma3p0.upf`: generated + validated (charge std
  2.121, V(0)=0.376).
- Binary smoked (b1/b2/b3) at σ=3/300 eV: exit 0, no NaN, clean WP injection.
