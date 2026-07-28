# Report 2 — Runs & Campaigns tables

*Compiled 2026-07-26. Companion files (same folder): `report-2-runs.csv` (all 899 run leaves, one row
each), `report-2-run-manifest.txt` (grouped path listing). Numbers/verdicts as recorded in the
handovers/campaign files; nothing filtered.*

---

## Table 1 — Campaigns (37 tracked in docs/campaigns/INDEX.md)

*Status/x-N from the campaign INDEX. "Outcome" copies the recorded verdict/headline; where a handover and
its campaign file record different wording that is noted, not resolved. (INDEX header says "38"; its own
status counts sum to 37 — the extras are auxiliary notes, not standalone campaigns.)*

| # | Campaign | Area | Status | x/N | Objective | Recorded outcome |
|---|---|---|---|---|---|---|
| 1 | Absorbing-boundary feasibility + MFA | absorbing_boundary | done | 4/4 | sin² CAP reflectivity ε(E,L) for absorber selection | ε(E,L) surfaces mapped; production params located |
| 2 | Jellium CAP baselines B0–B3 | cap_in_jellium | done | 5/5 | CAP drainage (B1) + projectile (B2/B3) baselines in jellium | baselines built; classical B2 S≈0.24 eV/Bohr |
| 3 | Stopping-power + Fourier training | check_logic | done | 4/5 | reproduce SP + Fourier pipeline on known cases | two skills encoded; validation dossiers filled |
| 4 | Claude ecosystem modularisation | codebase_rejuvination | done | 4/4 | modularise rules/skills/hooks | implementation C1–C11 + tests complete |
| 5 | Git commits — clean history | codebase_rejuvination | done | 2/2 | clean scoped commit history | commit voice locked; committed & pushed |
| 6 | Debugging quantum SP (p5_wp_v1p3) | debugging_quantum_stopping_power | done | 7/7 | is WP-S excess over Lindhard = CAP-captured fraction? | S 2.37→2.0; excess systematic 1.5–87×, CAP-capture ruled OUT |
| 7 | Jellium σ-convergence sweep | jellium_stopping | done | 5/5 | S(v) → point-charge Lindhard as σ→0 | σ={0.15,0.25,0.35,3.0}; handover/campaign wording differs |
| 8 | Energy-oscillation diagnosis | localised_jellium | done | 8/8 | single cause of ΔE_total>0 rise (ablation) | non-Hermitian CAP artifact (Graefe 2010), conf 0.90 |
| 9 | CAP wrap-around fix | localised_jellium | done | 6/6 | CAP setup giving monotone E_total | wrap η=−2.0 w40 → rise 0.000000 eV |
| 10 | Localised-jellium impl + scattering | localised_jellium | done | 4/4 | confined slab, measurable stopping | slab validated (r_s≈5.67, N=82) |
| 11 | pbc-open-z Arm B | localised_jellium | done | 4/4 | does z-periodicity drive the oscillation? | handover/campaign wording differs (refuted vs culprit) |
| 12 | σ=0.5 classical S(v) sweep | jellium_stopping | done | 2/2 | classical −1 electron S(v): Barkas/Lindhard | Bragg peak v≈1.0; Barkas crossover confirmed |
| 13 | Annular jellium S(v) vs wall r_s | cylindrical_jellium | running | 3/7 | electron down tube bore; S(v) vs r_s | classical S(r_s=6)≈0.004 Ha/Bohr; WP rung done; r_s=2 flagged |
| 14 | Classical stopping baseline (locjel) | localised_jellium | running | 4/6 | matched classical to bracket WP | S≈0.43–0.49 eV/Bohr |
| 15 | Wide low-spread WP σ=3.5 | localised_jellium | running | 4/9 | isolate quantum stopping from spreading | high-v deposit ≈+1.7 eV (S≈0.05 eV/Bohr) |
| 16 | Energy book-keeping analysis | localised_jellium_parameter_study_2 | running | 9/9 | WP−cl gap = ZPE + E_proj_bg; ledger closes | E_SIE=4.3 eV (PBC); closure ~1e-10 Ha |
| 17 | No-CAP quantum-vs-classical twin σ=3/300 | quantum_classical_nocap | running | 0/3 | remove CAP to isolate absorber effect | design locked; not launched |
| 18 | KS-WP vs HF orbital | td-hf | blocked | 0/7 | is KS-WP a good HF-orbital approx? | gated on slab validation + TD-HF pilot |
| 19 | Graphene/coronene CAP scattering | absorbing_boundary | paused | 2/6 | perp + grazing projectile scattering | WP arm done (ε≈0.12, ~85% abs); classical blocked |
| 20 | inq-stack unit testing + restructure | codebase_rejuvination | paused | 1/5 | per-component unit tests + restructure | paused 2026-06-22 |
| 21 | Quantum-vs-classical WP stopping (locjel) | jellium_wp_stopping | paused | 4/8 | S differs WP vs classical; quantify SIE | P3.1 big-box pair built |
| 22 | Localised-jellium GS study | localised_jellium | paused | 7/9 | GS + SIE analytical mental models | E_SIE 4.3 eV locked; open-z G=0 follow-ups |
| 23 | With-CAP quantum-vs-classical twin σ=3 | quantum_classical_nocap | paused | 0/4 | production absorber reference (A−B) | after no-CAP twin |
| 24 | Loss-function feasibility for WP SP | cap_in_jellium | ready | 2/8 | L(q,ω) route to WP S(v) | 2/8; Fourier gate lifted |
| 25 | Phase-5 WP S(E) velocity sweep | localised_jellium | ready | 0/8 | quantum S(E), v∈{1.3,3,5} | ready to run |
| 26 | ML pattern-finding (WP vs classical) | ml-patterns | ready | 7/15 | quantum signatures in induced density | POD/DMD: WP rank-1 coherent vs classical rank-4 |
| 27 | Nazarov–Gross mass sweep | nazarov_gross_comparison | ready | 1/6 | does mass affect low-v stopping? | null branch only; slow branch staged |
| 28 | Li multi-k ω_peak(v) drift | quantum_kick_extension | ready | 3/6 | is ω_peak(v) drift genuine? | 3/6 single-k; multi-k ready |
| 29 | σ=0.5 production absorber baselines | absorbing_boundary | ready | 0/2 | ε(E) for σ=0.5 dispersing packet | ready to run |
| 30 | Classical-vs-WP jellium (with CAP) | cap_in_jellium | draft | 0/0 | comparison framework | sketch |
| 31 | Fixing the classical projectile (locjel) | localised_jellium | draft | 0/6 | absorbed, cutoff-corrected S_cl(v) | gated on GS study |
| 32 | High-density classical S(v) benchmark | localised_jellium | draft | 0/8 | r_s≈4.2 open-z no-CAP classical | draft design |
| 33 | Stopping from decomposed ledger | localised_jellium | draft | 0/4 | ledger-split S for WP + classical | draft scope |
| 34 | Locjel dynamics — E_proj_bg + r_cut | localised_jellium_dynamics_analysis | draft | 5/7 | complete ledger term + r_cut sweep | mostly done (r_cut 50/120) |
| 35 | Effective-mass-tuned bands | mass_tuned_bands | draft | 0/2 | m* emulates materials via E(k) | draft framework |
| 36 | Muon projectile WP-vs-classical | muon_projectile | draft | 0/5 | muon slow-spread cleaner comparison | draft design |
| 37 | Sigma effect on stopping power | sigma_effect_on_stopping_power | draft | 0/0 | (no hypothesis text) | stub |

---

## Table 2 — Runs by system → sweep (899 run leaves; full per-run rows in `report-2-runs.csv`)

| System | Sweep / group | # runs |
|---|---|---|
| localised_jellium | campaign_autorun | 216 |
| localised_jellium | localised_jellium_dynamics | 83 |
| localised_jellium | muon_mass_fork | 33 |
| localised_jellium | cap_fix | 18 |
| localised_jellium | semiempirical_spillout | 14 |
| localised_jellium | classical_highdensity_sv | 11 |
| localised_jellium | qsp_phase5 | 7 |
| localised_jellium | qsp_phase3 | 7 |
| localised_jellium | h0_base_difference | 5 |
| localised_jellium | wide_wp | 4 |
| localised_jellium | sigma1_masspair | 4 |
| localised_jellium | ke_check | 4 |
| localised_jellium | energy_oscillation_diagnosis | 4 |
| localised_jellium | wp_cap_energy_plateau | 3 |
| localised_jellium | shared_gs | 3 |
| localised_jellium | qsp_phase4 | 3 |
| localised_jellium | fullsuite_wp | 3 |
| localised_jellium | qsp_phase2 | 2 |
| localised_jellium | qsp_phase1 | 2 |
| localised_jellium | mass_pair_n162 | 2 |
| localised_jellium | classical_slab_stopping | 2 |
| localised_jellium | 03_cap_stopping | 2 |
| localised_jellium | pbc_open_z | 1 |
| localised_jellium | nazarov_gross | 1 |
| localised_jellium | fullsuite_classical | 1 |
| localised_jellium | 02_projectile_slab | 1 |
| localised_jellium | 01_slab_validation | 1 |
| **localised_jellium** | **(subtotal)** | **437** |
| vacuum | twosided_cap_vs_mask | 107 |
| vacuum | mfa_sweep | 72 |
| vacuum | cap_thin_L5 | 33 |
| vacuum | cap_real | 17 |
| vacuum | cap_monomial | 17 |
| vacuum | cap_Lopt_E10 | 15 |
| vacuum | wp_traversal_energy | 2 |
| vacuum | cap_sweep | 1 |
| **vacuum** | **(subtotal)** | **264** |
| jellium | cap_baselines | 22 |
| jellium | save_gs | 9 |
| jellium | run_classical_n162_L50_sv_sigma3p0 | 8 |
| jellium | run_sv_sigma0p5 | 7 |
| jellium | run_classical_n162_L50_sv_sigma0p15 | 7 |
| jellium | run_classical_n162_L50_sv_sigma0p35 | 6 |
| jellium | run_classical_n162_L50_sv_sigma0p25 | 6 |
| jellium | legacy_jellium | 4 |
| jellium | (individual run_* dirs, 1 each) | 73 |
| **jellium** | **(subtotal)** | **142** |
| graphene | cap_scattering | 8 |
| graphene | grazing | 4 |
| graphene | perp | 3 |
| graphene | cap_cl | 3 |
| graphene | gs | 1 |
| graphene | cap | 1 |
| **graphene** | **(subtotal)** | **20** |
| coronene | save_gs | 3 |
| coronene | (individual run_* dirs, 1 each) | 14 |
| **coronene** | **(subtotal)** | **17** |
| cylindrical_jellium | annular_sv | 13 |
| **cylindrical_jellium** | **(subtotal)** | **13** |
| QuantumKickExtension | inq-codebase (Li multi-k) | 6 |
| **QuantumKickExtension** | **(subtotal)** | **6** |
| — | **TOTAL** | **899** |
