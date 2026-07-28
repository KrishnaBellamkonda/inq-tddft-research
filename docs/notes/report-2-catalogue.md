# Report 2 — Work Catalogue (30 May – 26 Jul 2026)

*Factual catalogue of all simulation work for Report 2. No inclusion/ranking judgements — status and
each item's own recorded outcome only. Companion files (same folder):*
- *`report-2-tables.md` — campaigns table (37 rows) + runs-by-sweep rollup.*
- *`report-2-runs.csv` — all 899 run leaves, one row each (`system,sweep,run_leaf,path`).*
- *`report-2-run-manifest.txt` — grouped path listing of every run leaf.*

Status legend (state only, not a judgement): ✅ complete · 🟡 partial/provisional · 🔴 in-progress/blocked · ⚙️ infrastructure.

## Boundary — what Report 1 already covered (submitted 29 May 2026)
Report 1 delivered: the WP rt-TDDFT framework + verified **free** wavepacket propagation; the
**periodic-jellium σ-sweep** (σ ≪ r_s coupling regime); a **σ = 1 WP↔classical agreement (2.4 %)**;
and a **coronene LEED (Tsubonoya 2014) replication setup**. Report 1's outlook named Report 2's job:
*extend the σ-sweep, complete the stopping-power-definition comparison, characterise the
quantum→classical crossover.* → **Report 2 = the work from 30 May → 26 Jul 2026.**

## Catalogue by theme

### A. Localised-jellium slab system (the big methodological shift)
- New **bounded (z-open) jellium slab** replacing periodic jellium → clean total-energy book-keeping so
  ΔE_electronic can be read as stopping. ✅ implemented + scattering Phases 2–5.
- **GS parameter study**: PBC (periodicity-3) **E_SIE = 4.3 eV** locked (PZ 4.4 / far-launch 4.55 /
  PBC-r40 4.26 agree); open-z periodicity-2 gives unphysical −2.1 eV (net-charged-cell G=0 monopole) →
  **use PBC for SIE**. SIE ∝ 1/σ_WP. 🟡 (open-z G=0 fix pending)

### B. Energy-decomposition (pairwise-Coulomb) ledger + twin runs
- Full pairwise E_PP/PS/SS/SB/PB/BB + kinetic split + E_xc, **closure ~1e-10 Ha**. ✅
- Overnight **5 twin pairs** (classical Gaussian-charge vs WP) + per-pair notebooks. ✅
- Zero-point KE **81.6 / 20.4 / 5.1 eV (σ=0.5/1/2)** (1/σ²); **SIE 4.34 / 1.13 / −0.25 eV**.
  Ghost-background correction mandatory for classical subtraction. twin-run skills built. ⚙️✅

### C. Stopping power — classical baselines (localised slab)
- Classical Gaussian-charge S in slab: **≈ 0.24 eV/Bohr @100 eV (CAP baseline B2)**; initial-drag
  **S(v₀=1.3) = 0.49 ± 0.004**; const-v **0.43 ± 0.18**. Factor-2 self-consistency. ✅
- Light-projectile deceleration method (initial-drag v ≥ 0.85 v₀) validated + codified. ✅

### D. Stopping power — the wavepacket "paradox"
- WP projected S ≈ **2.37 (v=1.3) → 9.09 (v=5)** eV/Bohr, **5–87× above Lindhard**. 🟡 UPPER BOUND.
- Central finding: WP deposit (59 eV) **exceeds** projectile KE drift (23–54 eV) → surplus is **WP
  zero-point KE (81.6 eV) + self-interaction**, NOT stopping. **WP-projected S is not a stopping power.**
  CAP-capture correction ruled OUT (excess systematic 1.5–87×). ✅
- Definition-1 (ledger-split) S for WP **blocked**: needs a fully-absorbing dynamic ledger run. 🔴

### E. Absorbing boundary / CAP
- CAP built in **inq-study** (stock inq cannot compile CAP). ε(η), ε(L), monomial-vs-sin² maps; anchor
  **E=10, L=20, η=−0.5**; monomial n=1 beats sin² (8.3 % vs 20.9 %). ✅ (provisional per Task #7)
- **Energy-oscillation diagnosis**: E_total drain-then-rise = **non-Hermitian CAP + kinetic-covariance
  renormalisation (Graefe 2010)**, not electrostatics (PBC/open-z refuted as cause). ✅
- **CAP fix**: wrap-around η=−2.0 w=40 → **rise = 0.000000 eV (strictly monotone)**. ✅

### F. σ-convergence toward point-charge Lindhard (extends Report 1's σ-sweep)
- Classical σ-sweep σ={0.15,0.25,0.35,0.5} (+σ=3 launched): low/mid-v **σ-independent**; high-v weak
  σ-dependence; shortfall vs Lindhard persists → **inconclusive; σ=3 pending**. 🟡
- Overnight classical S(v) (σ=0.5, 6 pts): **Bragg peak v≈1.0, S≈0.046 Ha/Bohr; Barkas crossover
  confirmed** (sim/LR 1.32→0.66 across v=2.98→0.62). ✅

### G. Related systems / extensions
- **Cylindrical annular jellium**: classical **S(r_s=6) ≈ 0.004 Ha/Bohr** (3 pts); WP rung done;
  r_s=2 flagged. 🟡 running (3/7).
- **Muon / effective-mass fork** (per-state mass, inq-study): GPU Tier-1 + bit-for-bit PASS; σ=1
  WP-vs-classical **S_classical = 0.208 eV/Bohr**; μ=206.77 ruled out by cost; Phase-5 XC pick paused. 🟡
- **Nazarov–Gross mass-sweep** validation: designed + GS built, **not launched**. 🔴
- **mass-pair N=162 σ=1 CAP** (m=1 vs m=2 WP, full ledger): **paused mid-flight** w/ checkpoints. 🔴
- **Graphene/coronene CAP scattering**: WP arm complete (perp+grazing, ε≈0.12, ~85 % abs); classical
  arm blocked (charged-cell energy-reference bug). 🟡
- **ML pattern-finding** (Δn = n_WP − n_classical): **POD/DMD headline** — WP → **rank-1, coherent,
  low-freq (7–11 eV)**; classical point → **rank-4, incoherent, high-freq (~210 eV)** e-h. Linear-response
  residual ≈ Gaussian low-pass. PDE path deprecated (blob artifact). 🟡

### H. Infrastructure / tooling
- inqview **4-package restructure**; inqkit **Catch2 tests**; **minimum observable set** + manifest
  (ADR 0006); **canonical figure standard** (ADR 0004); **density-matrix GIFs**; Gmail notifier;
  run-catalogue; **Claude ecosystem** modularisation. ⚙️

## Plots that already have data
Classical S(v) Bragg/Barkas (`sv_stopping_extraction.ipynb`); σ-convergence family
(`figures/sv_convergence.png`); CAP ε(η)/ε(L)/monomial (absorbing-boundary notebooks); energy
oscillation + CAP-fix monotone (cap_fix / energy-oscillation notebooks); pairwise ledger + closure
(qsp_phase2 notebooks); POD/DMD rank contrast (`bathstruct_*_sweep.png`); classical-vs-WP S
(`S_comparison.png`); many density-matrix GIFs.

---
---

# DETAILED INVENTORY — every run, sweep, and campaign (compiled from disk 2026-07-26)

*Harvested by fanned-out enumeration of `ResearchProject/systems/**` (run_summary.txt / results.json /
CSVs / executed notebooks) + all 38 `docs/campaigns/**`. Numbers verbatim; "not stated" where absent.
Scale: **899 run-leaf directories** (each with a run_summary.txt) — localised_jellium 437, jellium 142,
vacuum 264, cylindrical 13, coronene 17, graphene 20, QuantumKickExtension 6. **Every leaf is enumerated
verbatim in the companion file `report-2-run-manifest.txt` (same folder) — nothing filtered.***

> **Note — points where a handover and its campaign file record differing wording (stated, not resolved):**
> (i) Arm B / z-periodicity: `pbc-open-z` handover states *"REFUTED as cause"* (open-z reproduces the
> oscillation, t_min shifts 0.2 a.u.); campaign file states *"z-periodicity culprit / G=0 convention is driver."*
> (ii) σ-convergence: `sigma-convergence` handover states *"inconclusive; low/mid-v σ-independent; σ=3
> pending"*; campaign file states *"monotonic convergence toward Lindhard."* Same σ={0.15,0.25,0.35} data.
>
> **Date note (jellium periodic):** all periodic-jellium run dirs have mtime 2026-06-08→15 (the systems-reorg
> commit), so file dates are the reorg date, not the physics date; the σ-sweep physics predates 30 May.

## SYSTEM 1 — localised_jellium (CORE of Report 2)

### qsp_phase1 — SIE / insertion energetics (t=0, σ=0.5, 90-box)
E_GS(90-box) = −70.226 Ha; E_total(0) = −38.952 Ha; launch_z = −32 Bohr; KE_WP measured 180.8 eV vs
predicted 100 (drift) + 81.6 (zero-point) = 181.6 (−0.83 eV deficit). SIE measured 4.4 eV (dual route);
a-route 85.2 eV (boundary-matched) → 80.8 eV unattributed gap = zero-point contamination.

### qsp_phase2 (100 eV, σ=0.5, τ=40 au)
S_WP(deposited) = **2.727 eV/Bohr** (NOT converged, upper bound); S_classical(face) = 0.507;
S_Lindhard(point,100eV) = 0.282. E_deposited_WP = 68.2 eV; SIE floor 4.4 eV; WP absorbed 86.4 %;
σ_z ×41.3 (0.354→14.6 Bohr); T_plasmon 48.9 au. 7 notebooks.

### qsp_phase3 (100 eV, σ=0.5, τ=100 au, 90-box)
S_WP = **2.376** (not conv.); S_classical(face) = 0.000 (ion bounced, never entered); Lindhard 0.282.
E_deposited_WP 59.4 eV; **CAP removed 126 eV** (WP-carried ~173 + bath ~−47); WP absorbed 95.4 %;
σ_z ×65. 4 notebooks, results.json on disk.

### qsp_phase4 (54 eV, σ=0.5, τ=100 au)
S_WP = **2.391**; **S_classical(deposited) = 0.249** (gauge-clean); Lindhard(54eV) 0.448.
E_dep_WP 59.8 / E_dep_classical 6.2 eV; CAP removed 81 eV; WP absorbed 90.6 %; σ_z ×67. 3 notebooks.

### qsp_phase5 — velocity ladder S_WP(v) (σ=0.5)
| v (v₀) | E (eV) | S_WP (eV/Bohr) | E_dep (eV) | norm_f | bound |
|---|---|---|---|---|---|
| 1.3 | 22.99 | **2.374** | 59.35 | 0.095 | upper |
| 2.0 | 54.42 | **2.391** | 59.77 | 0.094 | upper |
| 3.0 | 122.45 | **2.566** | 64.14 | 0.062 | upper |
| 4.0 | 217.69 | **4.505** | 112.61 | 0.030 | upper (aliasing onset) |
| 5.0 | 340.14 | **9.782** | 244.55 | 0.057 | upper (aliased) |
| 6.0 | 489.80 | **18.895** | 472.38 | 0.026 | lower (E_total rises, unphysical) |
Vetted window E≤122 eV (v≤3), aliased tail <0.1 %. All flagged NOT converged. 7 notebooks.

### classical baselines (localised slab, σ=0.5, r_s=5.667)
- **P1 Ehrenfest S(v₀=1.3) = 0.4926 ± 0.0043** (in-slab v≥1.1; initial-drag). classical_slab_stopping.
- **P2 const-v (deposit) = 0.4315 ± 0.1767**; raw E_tot 0.4042 (not recommended).
- Ehrenfest ladder: S(1.3)=0.9316, S(2.0)=0.5102, S(3.0)=0.2498 (bulk σ=0.5, in-slab drag).
- Reference: bulk_σ0.5 = 0.937; Lindhard_point(r_s=5.667) = 0.5653; Lindhard_σ=0.354 = 0.0612.

### perturbation_method — σ-ladder (r=12, L_z=120, p2) confirms scalings
| σ | E_kin_WP eV | dE_H eV | dXC eV | U_proj_bg eV | dKin eV | **SIE eV** |
|---|---|---|---|---|---|---|
| 0.35 | 252.78 | 176.04 | −25.04 | 141.16 | 177.82 | **9.84** |
| 0.5 | 156.71 | 155.49 | −16.47 | 134.69 | 81.74 | **4.34** |
| 0.7 | 116.62 | 149.09 | −12.01 | 134.65 | 41.65 | **2.44** |
| 1.0 | 95.38 | 144.44 | −8.65 | 134.65 | 20.41 | **1.14** |
ZPE ∝ 1/σ², self-Hartree ~const, SIE→0 as σ↑ (one-electron artifact). Matches book-keeping numbers.

### localised_jellium_dynamics — r_cut sweep (r_s=5.667)
r_cut=50: dE_WP 81.23, dE_CL 185.39, dKin 81.75, dXC −16.47, U_proj_bg 90.71, E_proj 141.46, self_H 20.81 eV.
r_cut=120: dE_WP 80.55, U_proj_bg 178.67, E_proj 532.81, self_H 20.81 (radius truncation ~0.68 eV on dE_WP; ghost tail pathology in E_proj).

### campaign_autorun — periodicity ledger (p2 vs p3, 216 runs)
Δ(p3−p2) in dE_WP ≈ +6.0–6.4 eV across r=4..40; d(H+E) gap grows +4.2→+7.5 eV with r; **dKin gap = 0.0**
(WP quantum self-energy periodicity-independent). Origin = G=0 Poisson convention for charged cell. E_GS
p2 = 1643.1 eV, p3 = −2953.4 eV (box-dependent background). a1/a2/a3 audit notebooks.

### muon_mass_fork — per-state mass (33 runs, 15 notebooks) + σ=1 CAP ablation
CAP study (σ=1): baseline η=−1.0 → 100 % absorbed, 0.134 % refl, N_min 52, E_end −27.67 Ha;
gap19.5 → 99.68 %/1.147 %; weak η=−0.4 → 98.75 %/4.139 %; strong η=−2.0 → 100 %/0.001 %, E_end −34.88.
Fork validated GPU Tier-1 + bit-for-bit (He GS/RT identical to 14 digits). σ=1 classical twin S=0.208 eV/Bohr.

### classical_highdensity_sv — pilot (r_s=4.18 slab, v=2.0)
S_KE_loss(in-slab) = 0.93 eV/Bohr (55.16→32.04 eV over 25 Bohr); S_E_absorbed = 1.08 (neutral-cell baseline);
−ΔKE_proj 27.8 eV; transit −30→+70.9 Bohr; plateau flat to 0.0000 (z-open, CAP-free).

### Other localised_jellium sweeps present
cap_fix (18 runs, wrap η=−2.0 w40 → rise 0.000000 eV, study notebook), energy_oscillation_diagnosis
(ablation ladder; η=−1 −138→+31 eV, weak η=−0.2 −23.4→+0.11 crosses zero), pbc_open_z (comparison.md),
mass_pair_n162 (m1 vs m2 WP, paused mid-flight step 1500/1000 ckpt), sigma1_masspair, twin_dynamics,
stopping_from_decomposition (S_deposit 59 eV > drift KE 23–54 → "WP S is not a stopping power"), wide_wp
(S_WP=0.046 eV/Bohr high-v; PBC classical self-image artefact 5.6), extend_r160, h0_base_difference,
debugging_quantum_stopping, semiempirical_spillout, ke_check, 01_slab_validation, 02_projectile_slab
(12 GIFs), 03_cap_stopping (21 GIFs), nazarov_gross (staged), wp_cap_energy_plateau, plate_model (stub).

## SYSTEM 2 — jellium (periodic; Report-1 lineage + new σ-conv/CAP/qvc)  101 runs, 84 with results
- **7-run WP-RT reference table** (jellium_run_metadata.md): run_01_base…run_07_open_shell, E∈{50,200,400} eV,
  σ∈{0.5,1.0,3.78}, N_e 38 (closed) / 40 (open), r_s 7.38, 40³ box, ±z & 45° tilt.
- **Classical family (30):** highdens L30 E∈{50,100,200,300}; L50 E∈{20,25,50,100,200,300,600}; knudsen
  E∈{700..1100} (no results); **σ-convergence sv_sigma0p15/0p25/0p35/3p0**; e1000 rect box.
- **WP family (51):** highdens L30 (σ=1); L50 E-sweep; σ-spread E=100 σ∈{0.5,1,3,8}(+_wf); dx∈{0.40,0.80};
  knudsen E700(+minimal, +mpi_inject/propagate); e1500 cubic.
- **free_wp (5), plasmon (3, E∈{3.4,15,25}), sv_sigma0p5 (1), positive_ion (1), base (4).**
- **New post-May-30 hypotheses sweeps (executed notebooks):** 06_sigma_convergence (sv_convergence.png +
  energy variant); **cap_baselines** (B1 classical/B2/B3; 27 MB nb; drainage CSV; 40+ PNG + GIFs);
  **qvc_cap_sigma3** (b1/classical/WP, E=300 σ=3; 3 nb); **qvc_nocap_sigma3** (classical E∈{150,300,450,600}
  σ=3; 4 nb). Earlier 00–05 sweeps are metadata/overlay only. 05_electron_capture: f_trap 0.443 vs 0.168.

## SYSTEM 3 — vacuum (CAP / absorber method studies)  260 runs, 5 executed study notebooks
- **cap_real** (17 runs, E=22.29 eV, L=20, η-scan): ε(η=−0.5) = **3.04e-5** (min); ε(−0.10)=0.041/86.9 %,
  ε(−0.25)=0.0024/99.2 %, ε(−1)=1.21e-5. U-shape in |η|.
- **cap_thin_L5** (33 runs, L=5, k0×η grid): η=−0.30 ε∈[0.197,0.385]; η=−0.05 ε∈[0.53,0.75]; η=−0.01 ε∈[0.77,0.90]. cap_thin_combined.csv.
- **cap_monomial** (16 runs, E=10 L=5, n∈{1..4}×η): **n=1 (linear) best** — ε(η=−0.5)=0.083/91.5 % vs
  n=2 0.197, n=3 0.315, n=4 0.412. Higher n degrades. cap_monomial_combined.csv.
- **cap_Lopt_E10** (15 runs, L∈{6..20}×η): best **L=15 η=−0.5 → ε=0.00152**; L=20 η=−0.30 → 2.6e-4. cap_Lopt_combined.csv.
- **twosided_cap_vs_mask** (105 runs, E∈{10,100,1000}×L∈{10..30}×η∈{−0.3..−1}): anchor E=100 L=20 η=−0.5
  → ε=0.0051/98.8 %; L=26→5.3e-4, L=30→1.1e-4; E=1000 harder (ε=0.0246). twosided_combined.csv.
- **mfa_sweep** (72 runs, no CAP, E 0.544→489.8 eV × L∈{5..50}): bare dephasing ε(E,L); epsilon_grid.csv.
- **wp_traversal_energy** (2 runs, E=100 k0=2.711, CAP η=−0.7 vs nocap; 100 Bohr domain; notebooks unexecuted).

## SYSTEM 4 — cylindrical_jellium (annular tube)  14 runs, annular_sv sweep
- **Classical S(r_s,v) [Ha/Bohr] (Sv_results.csv):**
  r_s=6: v0.15→0.00199, v0.30→0.00394, v0.45→0.00608 · r_s=4: 0.00296 / 0.00584 / 0.00901 ·
  r_s=2: 0.00675 / 0.01347 / 0.02168 (r_s=2 flagged: ΔE r²=0.23 vs KE-route disagreement).
- **GS battery:** tube_rs6 (24 e, r_s 6.0, Lz48, 32 states), tube_rs4 (48 e, Lz28, 44), tube_rs2 (136 e, Lz10, 88).
- **WP rung:** wp_rs6_v0p30 (σ_WP=0.5, k0=0.3) completed (S not stated). annular_sv_report.ipynb (277 cells, 180 imgs).

## SYSTEM 5 — coronene (LEED replica + WP scattering)  ~15 runs, 0 notebooks (PNG/README analysis)
- **Tsubonoya-2014 paper_replica:** run_save_gs (E_gs −150.757 Ha, 108 e, cutoff 54 Ha, 34.77³-ish cell) +
  run_propagate (WP E=200 eV k0=3.834 σ=1.0, 600 steps, 20 LEED screens). leed_total_grid.png.
- **WP σ/E/geometry sweep (mostly Report-1 era, Apr 2026):** run_base, s0p33, s3, E30, E30_s3, E800, E800_s0p33,
  35x35x40 (E∈{30,200,800}, σ∈{0.33,1,3}); larger box b6/b18/broadening (impact param z₀∈{6,18,30});
  **run_broadening 2026-05-31 (Report-2), run_cc_bond 2026-06-11 (Report-2)** (WP at C–C bond midpoint).
- GS checkpoints gs_35x35x{40,60,80}_cut40 (−150.837 Ha, 108 e). HOMO-LUMO gap not stated in summaries.

## SYSTEM 6 — graphene (CAP scattering)  ~17 runs, 1 notebook (cap_scattering_study.ipynb)
- **WP+CAP (E=100 eV, k0=2.711, σ=1.47):** centroid_cap ε=0.1186/absorbed 84.65 %; channeling_cap
  ε=0.1197/84.40 %; nocap survival 0.72. **CAP halves survival, absorbs ~85 %.** Centroid≈channeling at 100 eV.
- **Classical arm (blocked→fixed):** z_valence=0 bug → 6.7 eV vacuum drag + 103 Ha offset; fixed via
  z_valence=−1 + extra_electrons(+1) (smoke KE_loss −0.004 eV). run_cl_centroid_s1/s2/s3 = 40.91 / −76.15 /
  10.27 eV KE-loss (seed-noisy, energy-reference still imperfect). Grazing/perp smokes only.
- **GS:** graphene 3×2 (96 e, −143.942 Ha, gapless Dirac); coronene-flake perp/grazing (108 e, −150.760/−150.774 Ha).

## CAMPAIGNS — 38 total (12 done · 5 running · 6 ready · 5 paused · 1 blocked · 8 draft)

### DONE (12)
absorbing_boundary feasibility+MFA (4/4); cap_in_jellium B0–B3 baselines (5/5); check_logic
stopping+Fourier training (4/5, two skills locked); ecosystem modularisation (4/4); git-commits (2/2);
**debugging quantum SP (7/7): S 2.37→S_corr≈2.0, CAP-capture correction**; **jellium σ-convergence sweep
(5/5)** (handover/campaign wording differs — see note); **energy-oscillation diagnosis (8/8): non-Hermitian
CAP = artifact, conf 0.90**; **CAP wrap-around fix (6/6): E_total plateau clean**; localised-jellium
impl+scattering (4/4, r_s≈5.67 N=82); **pbc-open-z Arm B (4/4)** (handover/campaign wording differs — see
note); **σ=0.5 classical S(v) (2/2): Barkas crossover**.

### RUNNING (5)
cylindrical annular S(v) vs r_s (3/7); classical stopping baseline locjel (4/6, S≈0.43–0.49);
wide low-spread WP σ=3.5 (4/9); **energy book-keeping (9/9 done-but-listed-running): E_SIE=4.3 eV PBC,
ledger closes**; no-CAP q-vs-c twin σ=3/300 eV (0/3, design locked).

### READY (6)
loss-function feasibility for WP SP (2/8); **Phase-5 WP S(E) velocity sweep (0/8)**; ml pattern-finding (7/15);
Nazarov–Gross mass sweep (1/6); Li multi-k ω_peak(v) drift (3/6); σ=0.5 production absorber baselines (0/2).

### PAUSED (5)
graphene/coronene CAP scattering (2/6, paused 06-21); inq-stack unit testing (1/5); **q-vs-c WP stopping
locjel (4/8, P3.1 big-box pair)**; localised-jellium GS study (7/9, E_SIE 4.3 eV locked); with-CAP twin σ=3 (0/4).

### BLOCKED (1)
td-hf: is KS-WP orbital a good approx to HF orbital? (0/7; gated on slab validation + TD-HF pilot).

### DRAFT (8)
cap_in_jellium classical-vs-WP (0/0); classical-projectile-fix locjel (0/6); high-density classical S(v)
benchmark (0/8); stopping-from-decomposed-ledger (0/4); locjel-dynamics E_proj_bg sweep (5/7);
effective-mass-tuned bands (0/2); muon projectile WP-vs-classical (0/5); sigma-effect stub (0/0).

### Objectives vs the three project goals
- **Goal 1 (WP framework + parameter influence):** localised-slab impl + CAP feasibility + energy ledger/SIE
  DONE; σ-convergence DONE; wide-σ & classical-twin method established; Phase-5 S(E) + cylindrical extend it.
- **Goal 2 (stopping-power-definition comparison):** energy-method (deposit/L_z) locked + debugged; loss-function
  route piloted (feasibility gate pending); documented finding: WP-projected deposit exceeds projectile KE
  drift (zero-point+SIE dominated); classical and Lindhard recorded as reference anchors.
- **Goal 3 (quantum→classical crossover):** point-charge limit via σ-convergence + Barkas crossover; high-v
  no-CAP/with-CAP σ=3 twins (in flight); effective-mass/muon route to lower zero-point KE (draft).

---
---

# ADDENDUM — deep pockets (notebook numbers · meeting decks · QBall refs · validation/citations)

*Compiled 2026-07-26 from four targeted sweeps. Verbatim; "not stated" where absent.*

## A. Numbers that live only inside notebook cells

### mass_pair_n162 (m=1 vs m=2 WP, N=162 slab, full decomposition; fully-absorbed analysis)
- **m1:** v=2.711, **S=0.77 eV/Bohr** (deposit 19.2 eV), norm_f 0.003 (absorbed 1.028 e⁻), E_GS −53.839 Ha,
  t_exit 10.7 au — S_WP/Lindhard(0.28) = **2.7×**.
- **m2:** v=1.917, **S=1.02 eV/Bohr** (deposit 25.5 eV), norm_f 0.008 (UPPER BOUND, still draining), —
  S_WP/Lindhard(0.47) = **2.2×**. Gap candidates: form-factor q-cutoff / non-linear Z=−1 / finite-box+CAP.
- Heuristics: r_s 5.686, k_F 0.3375, v_F 0.3375, E_F 1.55 eV, ω_p 3.48 eV, T_plasmon 49.2 au.

### sigma1_masspair (σ=1, mass-pair)
- wp_m2_k4.5: m=2, v=2.25, E_kin 138 eV, drain 60 eV, t_min 32.6 au, **rise +180 eV** (ledger artifact).
- wp_m3_k4.5: m=3, v=1.50, E_kin 92 eV, drain 44 eV, t_min 50.2 au, **rise +125 eV**.
- p3_wp_m1_rerun (clean ref): m=1, v=2.71, E_kin 100 eV, drain 126 eV, t_min 97.9 au, rise +0.11 eV.
- **Clock law: t_min ∝ 1/v** (CAP absorption completion, NOT slow-spill); norm_wp ≈ 5 % at t_min; rise is 100 % kinetic-ledger.

### twin_dynamics — 6-pair decomposition (gauge test passes all, no inter-run gauge)
| pair | σ | k0 | regime | ZPE_dKin_loc eV | R_selfHartree eV | SIE eV | dXC eV |
|---|---|---|---|---|---|---|---|
| s0.5_k1.0 | 0.5 | 1.0 | ladder | 81.72 | 20.81 | 4.34 | −16.47 |
| s1.0_k1.1 | 1.0 | 1.1 | ladder | 20.42 | 9.79 | 1.13 | −8.66 |
| s2.0_k1.1 | 2.0 | 1.1 | capture | 5.10 | 4.37 | −0.30 | −4.66 |
| s2.0_k4.2 | 2.0 | 4.2 | null | 5.10 | 4.37 | −0.25 | −4.61 |
| s2.0_k0.4 | 2.0 | 0.4 | reflect | 5.15 | 4.34 | −0.32 | −4.66 |
| s2.0_k0.5 | 2.0 | 0.5 | tunnel | 5.92 | 3.89 | −1.39 | −5.29 |
ZPE∝1/σ², self-Hartree∝1/σ, SIE→0 as σ↑. Dispersion/reflection/tunnelling map as term-by-term divergences.

### wide_wp (σ=3.5, high-v)
- wp_per2_E300_long: E=300 eV, k0=4.696, σ_wp 3.5, N 83→82, η=−1, open-z, 1773 steps.
- Raw ΔE_total rings ±21 eV (phantom KS-state-60, ω≈0.31 Ha); **remove that orbital → plateau std 0.00,
  physical deposit ≈ +1.7 eV** (Lindhard high-v floor S≈0.05 eV/Bohr). Ledger recoverable in post-processing.

### jellium qvc σ=3 twins
- **qvc_cap_sigma3 (E=300):** WP absorbed_frac 0.3752, classical 0.3719 (~37 %, comparable; classical slips
  through slightly more). Both reach ±25 Bohr (far face).
- **qvc_nocap_sigma3 (E=150/300/450/600):** absorbed ~2.6–2.9e-8 (numerical noise) → 100 % survival, E-independent.
- **cap_baselines drainage:** η=0.50 aggressive vs η=0.10 slow over τ=140 au (drainage transient CSV).

## B. What has already been SHOWN to the supervisor (meeting decks)

### 14–15 May 2026 — first presentation (Report-1 boundary; threads A–E)
- **A Coronene LEED:** 12 runs, σ∈{0.33,1,3}×E∈{30,200,800} eV; drift ~1e-5 Ha; overlap 1e-5 separates high-E.
- **B Jellium WP scoping:** run_base_n162_L50_E1p5 (E_WP 1.5 eV, N162, r_s5.69): v 0.373→0.103; ΔKE +0.50,
  ΔE_H −0.60, ΔE_xc +0.10; WP −1.75 eV, bath +2.25 eV; persistent density hole (orthogonalisation artefact).
- **C Plasmon:** m=1 **3.533 eV** (Bohm-Gross 3.59, −1.6 %); velocity discriminator 85× below kinematic.
- **D Li-54 kick:** v=0.0626 **6.480 eV** (paper 6.5); v=0.300 **2.585 eV** (paper 2.8); collective confirmed.
- **E Classical-vs-WP stopping:** E=1500 classical v=10.5 **S=0.021 eV/Bohr**; E=100 classical v=2.71
  **S=0.217** (Bethe v⁻² within 1 %); **WP E=100 ΔE_kin 0.210 eV → classical/WP ratio 13.9× (mid-traj 17×);
  Hartree sign-flip: classical +0.112 (accumulation) vs WP −0.239 (depletion anti-wake).**
  [full detail in classical-vs-wp-case-study.md]

### 26 Jun 2026 — deck (58 figures)
Classical S(v) σ=0.15–0.5 → Lindhard; σ_WP=0.5 "best point-like compromise" locked; **CAP chosen over mask**
(L,η fixed); localised-slab GS validated, **SIE=4.40 eV (E_GS −45.759 Ha)**; energy book-keeping S=ΔE/L_z with
decomposition; WP spreading law σ(t)=σ√(1+(ħt/2mσ²)²); muon outlook (m_rel 206.77); CAP B0–B3 inside jellium.
Fourier gate lifted 2026-06-25. Figures: fig_s1/s2/s3_* (stopping, system design, CAP/mask, reflectivity, SIE, energy).

### 3 Jul 2026 — planning update
Slab E_total(0)−E_GS decomposition (periodicity 3 & 2); classical-vs-WP energy-component breakdown = primary
metric; analytical infinite-plate + semiempirical model; UPF r_c cutoff study. Figs: component_decomposition_pbc,
excess_vs_r_pbc, right_formula_pbc.

### 9 Jul 2026 — deck (30+ figs)  ← the analytical-model + energetics deck
Plate model validated vs L_z=160 DFT; **image potential −0.76 eV; Gaussian-vs-point −7.7 meV; interior dipole
barrier ~1.8 eV; net projectile-slab Coulomb <1 eV**; identity dE_WP−dE_CL = dT_zp+dE_xc+d(H+E) to machine
precision; **ZPE 3/(4σ²)=81.6 eV (measured 81.7)**; far-field (r=40) **dE_WP−dE_CL = +67.3 eV** = ZPE 81.6 +
Coulomb 2.0 − self-XC 16.5; wide WP σ=3.5 deposit ≈0; high-mass m=2.506 spreads 0.4 % vs electron 18×.
Figs: s1_* (plate/slab), s2_* (energetics/cutoff/screening), s3_* (wide WP, muon, density triptych GIFs).

### 16 Jul 2026 — working figure dir (no minutes): classical-gaussian-perturbation (9 figs), classical-vs-quantum, energy-oscillation.

### overnight-gaussian-classical-jellium/REPORT.md (11–12 Jun)
Gaussian classical σ=0.5 (V(0)=1.596 Ha repulsive), r_s5.69 N162; 15–16 s/step; **v=3.0 ⟨S⟩=0.223 eV/Bohr
(1.22× Lindhard); Bragg peak v≈1.0 S≈0.046 Ha/Bohr**; Lindhard code fixed (f-sum=1.000, 10 tests); cosine-kick
implemented (loss-function production cost-blocked ~78 h); σ=0.4 → 1.05× σ=0.5; Friedel λ≈9.3 Bohr.

### Synthesis — already shown vs still raw
- **Shown:** classical S(v)→Lindhard + Bragg/Barkas; CAP-vs-mask method; slab GS + SIE 4.40 eV; classical-vs-WP
  14–17× + Hartree sign-flip; plasmon 3.53 eV; Li kick 6.48 eV; plate model; book-keeping identity; spreading law; muon outlook.
- **Raw / not yet presented:** production quantum S(v) sweep (qsp_phase5); wide-WP quantitative; loss-function
  L(q,ω); muon fork execution; positive-ion companion; orthogonalisation-hole Gram-Schmidt rerun; cylindrical S(r_s); ml POD/DMD.

## C. QuantumKickExtension — QBall reference benchmarks (validation anchors, not INQ runs)
- **Li BCC-54 (QBall, Γ, 8 velocities v=0.0123..0.450):** energy plateau ΔE/N_uc vs Mv² reference
  (M_uc=6 m_e); **ω_peak: v=0.0626 → 6.48 eV (dipole_x), v=0.300 → 2.585 eV, v=0.450 → 2.620 eV.**
  Constants: E_F 4.7 eV, v_F 0.588, **Li plasmon DFT-RPA 6.56 eV** (arXiv:2510.07261), ecut 74 Ry, dt 0.04.
- **Li INQ multi-k (2×2×2 shifted, 400 K, 15500 steps, v=0.0123/0.0626/0.300/0.450):** run data present
  (observables.csv); **plateau/ω_peak extraction NOT yet done** — gated on fourier-analysis skill (quantum_kick_extension campaign, ready 3/6).
- **Diamond (64-atom):** raw kick .out present, no analysis. **Al (108-atom):** 2 exploratory kicks (v=0.03,0.09), stub. n2: GS only.
- Reference paper: Santervás-Arranz, Stengel, Artacho, PRResearch 7, 033292 (2025) — kick + Mv² diagnostic.

## D. Validation ledger + citation bank

### Validation notes (docs/validation/) — verdicts
| File | Validates | Verdict |
|---|---|---|
| e-proj-bg-dual-route | E_proj_bg via closed-form vs FFT-Poisson | **PASS** — max diff 0.20 eV (0.25 %); 4-term decomp reproduces d(H+E) to ±4 eV (~2 %) |
| loss-function-formula | L(q,ω)=\|n_q\|²/q² vs −Im[1/ε] | **ACCEPTED as peak-LOCATOR** (poles exact); NOT faithful line-shape/area; 3/3 tests pass |
| fft-drift-removal | baseline convention | **User verdict: uniform `mean` baseline**; peaks robust (3e-7); only ω≈0 affected |
| fft-normalization | coherent-gain / window | **Fixed**: normalize by Σwin; Hann recovers true A; Parseval ratio 1.0 |
| inqkit-errors | E01–E04 register | E01 (MPI Allreduce) FIXED; E03 (2-pass GS) FIXED; E02 (stale bath) + E04 (dx/2 COD offset) CONFIRMED/deferred |
| inqview-findings | IV-E01–E03, IV-M | IV-E01 complex-FFT fix; IV-E02 relabel power≠loss-fn; IV-E03 coherent-gain fix; baseline→mean |
| smoke-tests-launch-readiness | QKE + CAP Tier-A GPU | **PASS** (after COD-API + kick-iter fixes); campaigns launch-ready |
| test-catalogue | full test index | inqkit 26 pure+20 engine; inqview 29 pass+2 xfail; CAP 8; muon 7; book-keeping 5; cap_fix 3 |
| coronene-replication | Tsubonoya Tier-A | **NOT YET RUN** (matrix defined) |
| inqkit-tests | API coverage | locked 2026-06-10; BL-* baselines not yet captured |

### Citation bank (docs/sources/) — ~30 notes; the ones load-bearing for Report 2
correa-2018 (rt-TDDFT stopping review: Bethe/Fermi-Teller/Lindhard, energy decomp, transient-vs-steady) ·
**nazarov-gross-2025** (quantum WP stopping, mass/width-dependent friction, classical→Lindhard limit) ·
**graefe-2010** (non-Hermitian CAP: energy rate sign-indefinite — grounds the oscillation diagnosis) ·
**selsto-2010** (CAP bookkeeping gap: absorbed sector neither conserved nor monotone) ·
lindhard-1976 (Barkas Z₁³) · stopping-power-formulae (Bethe/Bohr/Bloch/Lindhard) · stopping-power-jellium-anchors
(RPA χ⁰, low-v friction Z=−1) · **tsubonoya-2014** (coronene LEED model system) · **santervas-arranz-2025**
(kick + Mv² diagnostic) · penn-1987 + matias-2025 (Penn: jellium-of-varying-r_s as material proxy) ·
quijada-2007 (finite cluster ≈ bulk at high v) · heg-mass-scaling-xc + kreibich-gross (muon mass-scaled XC / MC-DFT) ·
segui-arista-2007 + echenique-ritchie (cylindrical/tube stopping + wake/image) · kavokine-2022 + coquinot-lizee-2023
+ kral-shapiro-2001 (hydrovoltaic/quantum-friction context for tube) · ward-2024 + kononov-2025 + chiang-2025
(ML-stopping / nonlinear-stopping / density-trajectory precedents for ml-patterns) · monkhorst-pack-1976 ·
dipole_as_q0 · li_gs_xyz · free-electron-gas-magic · lee-water-dna-20ev (radiobiology energy scale) · car-parrinello.

### Citation GAPS to fill before writing
- **Primary Lindhard (1954)** — only cited via secondary sources; add a dedicated note.
- **Echenique–Nieminen–Ritchie (1981, 1986)** — nonlinear low-v jellium friction; no dedicated note.
- **POD/DMD foundations** (Tu 2014 exact-DMD; Brunton & Kutz 2019) — used in ml-patterns, no source notes.
- **Giuliani–Vignale plasmon χ(q→0,ω)** — cited inline, no note. Riss–Meyer/Muga CAP definitions (optional).
- Baseline captures (BL-coord/dens/parallel/wp) pending; E04 COD +dx/2 offset documented but not corrected.
