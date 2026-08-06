# Plan: Cylindrical (annular) jellium — KS-orbital WP stopping in the channeling limit

**Status:** IMPLEMENTED, not yet run. Branch `quantum-stopping-power`, device CSD3.
**Handover:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/handovers/cylindrical-channeling-ks-stopping.md`
**Design:** locked by the user (this file is the user's plan, expanded with the
implementation decisions taken while building it; the §2 parameter table is
unchanged from the locked version).

---

## 1. Aim (the one-sentence version)

Demonstrate that a **quantum electron wavepacket (WP)**, shot as a *channeling*
projectile down the hollow bore of an annular jellium tube, **reproduces the
classical projectile's stopping power** in the low-interaction / high-velocity
limit — thereby validating a **KS-orbital-based definition of stopping power on
the WP** against the established ΔE_total/ds classical definition.

### Physical rationale (why this should work)

- The bulk-WP contamination seen earlier was an **interaction-driven growth of
  var(p)** (the momentum-spread term of the orbital KE rose +6.8 eV while the
  drift ½⟨p⟩² stayed flat). var(p) is **conserved under free evolution**, so the
  contamination is interaction, not dispersion.
- **Channeling** (WP flies down the empty bore, couples to the wall only via its
  smooth image force) suppresses that interaction → var(p) stays ≈ frozen → the
  WP's energy change collapses to the **drift channel ½⟨p⟩²**, which *is* the
  classical projectile KE.
- **Matched charge width** (σ_pot = σ_WP/√2) makes the WP source the *same*
  potential as the classical projectile, so in the low-interaction limit the WP
  force = classical force → same stopping.
- **High velocity (50 eV)** is *required*: at low v the packet disperses
  (τ ~ σ_WP²) before it can traverse; high v crosses the region before dispersion
  sets in and keeps the projectile at ~constant velocity.

### The claim has three parts, and all three are measured

| | what it says | measured by |
|---|---|---|
| **Result** | S_WP lands on S_classical | the four S_ij fits vs the classical fit, same window |
| **Premise** | the packet stayed in the bore | `f_bore(t)` from `wp_radial_occupancy.csv` |
| **Mechanism** | channeling froze var(p) | `sigma_pz2(t)` and (T1−T2)(t) |

Result *without* premise and mechanism is a coincidence, not a validation.
`channeling_stopping.compare()` returns **AIM MET only when all three hold**, and
names the failing one otherwise.

---

## 2. Locked parameters

### Geometry (annular tube, channeling)

| Parameter | Value | Note |
|---|---|---|
| Density | **r_s = 3** | n0 = 8.84194e-3 a.u. (**exact**, see below), v_F = 0.6397, ħω_p = 9.07 eV. NEW GS required. |
| Bore inner radius | **R_in = 10 Bohr** | clears the WP's dispersed 2σ_dens (7.3 Bohr at the end of the fit window) by 2.7 Bohr |
| Wall outer radius | **R_out = 14 Bohr** | wall thickness 4 Bohr ≈ 3.6 × screening length (1/k_TF = 1.108) |
| Tube length | **L_z = 60 Bohr** | = 1.66·λ_p (λ_p = 36.1 Bohr); do NOT shorten below ~48 |
| Jellium edge smoothing | **0.5 Bohr** | erfc width at BOTH ring edges (≈ 1 grid cell) — confirmed to mean the jellium ring edge, not the projectile |
| Cell / box | **40 × 40 × 60 Bohr** | z ∈ [−30, 30]; 6 Bohr vacuum R_out→edge; fully periodic |
| Grid spacing | **dx = 0.5 Bohr** | 80×80×120; E_cut = 537 eV |
| Bath electrons | **N = 160** (104 states) | 80 occupied + 24 extra; the WP takes the last extra state |

**N = 160 is not a round number chosen for convenience** — it is the electron
count for which n0 = N/V_annulus lands on r_s = 3.000000 for this exact geometry
(V_annulus = π(R_out²−R_in²)L_z = 18095.5737 Bohr³). The binaries set
n0 = N/V so ∫n₊ = N *exactly*, which the G=0 cancellation of a periodic Poisson
solve requires.

### Projectile (BOTH runs, on-axis channeling)

| Parameter | Value | Note |
|---|---|---|
| Projectile energy | **50 eV** | v = 1.91701127 a.u.; **v/v_F = 3.00** |
| Launch | **on-axis (x=y=0), z = −28**, v_z = +1.917 | both runs identical |
| Charge / mass | −1 e / m_e | classical is Ehrenfest, WP is an occupied KS orbital |
| WP width | **σ_WP = 4 Bohr** (the label) | σ_dens = σ_WP/√2 = 2.8284; k₀σ_WP = 7.67; k₀/σ_p = 10.8 |
| Classical potential width | **σ_pot = σ_WP/√2 = 2.8284 Bohr** | MATCHED (user-confirmed). Implemented as a **moving-Gaussian CHARGE perturbation** (`v_proj = +poisson(n_proj)`), NOT UPF — no r_cut, no charge-sheet inflation. |

### Propagation

| Parameter | Value |
|---|---|
| dt | 0.02 a.u. |
| N_steps | 1500 (T = 30 a.u. → 57.5 Bohr of path, **one traversal, no wrap**) |
| Engine | **inq-study** |
| Propagator | ETRS, **no CAP** (norm- and energy-conserving) |
| Frame cadence | `CH_SAVE_EVERY = 5` → 300 density frames; wavefunctions every 75 |
| Checkpointing | interior every 500 steps + FINAL, `CH_RESUME=1` support (mandatory) |

---

## 3. Implementation decisions taken while building (not in the original plan)

These were forced by the geometry and are recorded because they change what the
analysis is allowed to claim.

### 3.1 No in-medium path correction — and that is a result, not an omission

The slab study needed `s5 = ∫f·v dt` because 25 of its 85 Bohr were vacuum, so a
centroid-path fit averaged the drag over medium *and* vacuum. **The tube is
uniform along z**: the medium fills every z the projectile visits, so the path IS
the in-medium path and `−dT/ds` is already a force. The correction the slab needed
does not exist here.

### 3.2 What the tube CAN violate is its own premise — hence `radial_occupancy`

New library observable `inqkit::observables::radial_occupancy`
(`inq-stack/include/inqkit/observables/radial_occupancy.hpp`): per step, the
fraction of |ψ_WP|² inside the bore (r⊥ < R_in), inside the wall
(R_in ≤ r⊥ < R_out) and outside, plus ⟨r⊥⟩ and σ_r⊥, minimum-image aware.

**The fit window is derived from the MEASURED f_bore(t)** — it ends where f_bore
first drops below 0.95 — not from the free-dispersion formula. The formula
(2σ_d = R_in at t = 23.3 a.u.) is plotted alongside as a cross-check. This matters
because the packet stops being Gaussian as soon as it scatters, and because a
*measured* breach earlier than the formula would itself be a physical finding
(the packet being pushed into the wall, not merely spreading into it).

### 3.3 Minimum image is mandatory in BOTH halves

The launch point is 2 Bohr from the −z face = **0.71 σ_pot**. A plain Cartesian
Gaussian keeps only Φ(2/2.83) = 76 % of its charge there, and the loss is
**asymmetric between the +δ and −δ finite-difference evaluations**, so it does not
cancel out of the force: the classical projectile would feel a fake force at
launch. The wavepacket has no such problem — a KS orbital lives on a plain 3-D FFT
basis and wraps exactly.

So `moving_gaussian_projectile_perturbation(..., minimum_image=true)` **and** a
matching minimum-image option added to `projectile_force_axis`/`projectile_force_z`
(new, defaulted false so every published run keeps its behaviour).

**And the wavepacket needs it too — this was got wrong first time round.** The
claim above that "the wavepacket has no such problem" is true of the PROPAGATION
and false of the INJECTION: `inqkit::WavePacket` built its Gaussian from a plain
Cartesian displacement, so the packet was TRUNCATED at the face and normalisation
hid it in the norm. Six of nine t=0 gates failed, all in z (var(p_z) fifteen times
too large from the truncation's sharp edge), while x and y were perfect.
`WavePacket::minimum_image(true)` (new, defaulted false) fixes it, and the twin
REQUIRES it: the classical half already wraps its charge, so a clipped packet is
not its twin at the very boundary this study introduces on purpose.

Equally, the WP's ~24 % of density on the far side of the cell at t = 0 makes the
**circular (Resta phase) centroid mandatory**; the naive ⟨z⟩ slides smoothly to a
wrong answer rather than jumping, so it cannot be repaired afterwards. ⟨p_z⟩, T1
and T2 are momentum-space expectation values and are unaffected, so the *primary*
measurement is clean regardless.

### 3.4 Full 3-D Ehrenfest for the classical half

The force is computed in all three directions and the projectile is free to leave
the axis. By the tube's symmetry the transverse force vanishes at r⊥ = 0, so
x(t) ≈ 0 becomes a **measured** statement about channeling stability instead of a
constraint imposed by integrating only z. F_x is written every step.

### 3.5 Correctness gates that this geometry makes available

- **WP:** no CAP ⇒ H is Hermitian and time-independent ⇒ `energy_total` is
  CONSERVED. Reported at the end of the run; a real drift invalidates every S.
- **Classical:** `E_electronic + KE_proj + U_proj_bg` must be flat; a drift means
  the Hellmann–Feynman force and the perturbation potential disagree.
- **Transverse images never overlap.** 6σ_d = L_xy only at t = 34.1 a.u., *after*
  the run ends — unlike the slab study, this run never drags a periodic array.

---

## 4. Required observables (BOTH runs)

**KS-stopping channels (the primary deliverable), two definitions each:**
- **Position:** (a) `s4 = ∫⟨p_z⟩ dt`, (b) `s3 =` circular centroid ⟨z⟩(t).
  Sources: `wp_momentum_stats.csv:pz_mean`, `wp_real_space_stats.csv:z_mean_circ`.
- **ΔKE:** (a) `T2 = ½⟨p_z⟩²`, (b) `T1 = ⟨p²⟩/2m`.
  `S_24` (drift-vs-drift) is the **headline**: built from ⟨p_z⟩ on both sides, so
  it is a stopping power whether or not var(p) is frozen. The other three are
  cross-checks, and their agreement with S_24 is itself evidence.
- `var(p) = sigma_pz²(t)` — the cleanliness/mechanism diagnostic.

**Energy decomposition / interactions (full pairwise ledger, both runs):**
`interactions.csv` with E_SS, E_PP, E_PS, E_SB, E_PB, E_BB + the closure checks,
every step, in an **identical 12-column schema in both halves** so one loader
reads both (`.claude/rules/decomposed-interaction-energies.md`).
**E_PP is the WP self-Hartree** — the only term with no classical counterpart, and
the leading suspect for any residual discrepancy in S.

**Density decomposition frames:** `density_total`, `density_wp`, `density_delta`
(+ coarse), `density_gs_system` baseline.

---

## 5. Analysis / deliverables

| Artefact | Path |
|---|---|
| Analysis engine (deterministic; all arithmetic) | `ResearchProject/systems/cylindrical_jellium/hypotheses/channeling_twin/channeling_stopping.py` |
| Engine tests (13 cases, exact) | `.../channeling_twin/tests/` |
| Per-run deep dives (2) | `.../channeling_twin/{classical,wp}_*.ipynb` via `build_run_notebooks.py` |
| **Phase notebook (the deliverable)** | `.../channeling_twin/channeling_twin_comparison.ipynb` via `build_comparison_notebook.py` |
| Summary table | `.../channeling_twin/stopping_summary.csv` |

The **phase notebook** carries, in this order: the twin density-matrix GIF
(mandatory, `.claude/rules/notebook-density-gif.md`), kinematics, **the S bar
chart with uncertainties**, the fits and their residuals, **f_bore/⟨r⊥⟩ vs the
bore radius**, **var(p_z) and T1−T2 vs their free-evolution values**, the
s3-vs-s4 consistency check, the pairwise ledger with E_PP and the closure
residuals, the correctness gates, and a **computed verdict**.

---

## 6. Prerequisites & gates (in order — this IS the dispatch chain)

`./shared/bin/submit-channeling-twin.sh` submits all of it:

1. **`chan-tests`** — library gate: `radial_occupancy`, minimum-image charge AND
   force, projectile wrap, slab occupancy. Nothing runs if these fail.
2. **`chan-gs`** — the shared r_s = 3 tube GS. Gates on exact neutrality, state
   count, and **bore depletion < 0.5** (if the bore is not electron-poor there is
   no channel to channel down). Idempotent.
3. **`chan-twin {wp,classical} smoke`** — builds each binary, 20 steps, t=0
   analytic gates.
4. **`chan-twin {wp,classical}`** — the two 1500-step production halves,
   concurrent.
5. **`chan-nb`** — `check_twin.py --dynamic` parity gate → 2 run notebooks →
   the comparison notebook.

**`cutoff_guard.py`: PASS** for both halves, run 2026-08-01 —
WP aliased tail 0.00 % (σ_p = 0.177 vs k_Nyq = 6.28), classical E_cut = 537 eV ≥
1.10 × 50 eV.

**No separate classical pilot.** The `chan-twin classical smoke` stage plus the
full-run conservation gate cover what the plan's pilot was for; the classical half
is not the expensive one and holds no GPU that the WP half needs.

---

## 7. Scope boundary (honest)

This validates the classical↔quantum bridge in the **fast (v/v_F ≈ 3) channeling**
regime, at **one velocity point**. The near-peak (v ≈ v_F) physics is inaccessible
to a non-dispersing WP — the WP method is intrinsically a high-velocity probe. The
Lindhard curves drawn in the per-run notebooks are **bulk** response functions
while the projectile is in a vacuum bore, so they are an upper reference, not a
prediction for this geometry.

A follow-on σ_WP or velocity sweep is a straightforward re-parameterisation of the
same three binaries (every physical parameter is an env var), but is not in scope
here.

---

## 8. Refined-analysis notebook (added 2026-08-02, user-directed)

The 222 MB comparison notebook answered "was the aim met?" with a window the
ANALYSIS chose (f_bore >= 0.95). The user's judgement is that the window should
instead be chosen BY EYE from the diagnostics, and separately for each half. This
section adds a **lightweight** notebook whose job is to show those diagnostics.

**Deliverable:** `hypotheses/channeling_twin/refined_analysis.ipynb`, built by
`build_refined_notebook.py` on top of a thin, tested data layer `refined.py`.
Lightweight = no embedded density GIFs (they are already in
`channeling_twin_comparison.ipynb`; re-embedding costs 222 MB and adds nothing).

**Section order is the user's, not mine:**

1. **Position.** Classical z(t) (trivial) against the WP's TWO position
   definitions — the circular (Resta-phase) centroid and `z0 + \int <p_z> dt`.
   Agreement is the Ehrenfest consistency check; divergence localises where the
   packet stops behaving like a particle.
2. **Kinetic energy.**
   - classical: `Delta E_total`(bath) and `1/2 m v^2`(projectile) on one axis,
     plus their sum as the closure residual. Verified 2026-08-02: closes to
     2.2e-5 eV over 1501 steps.
   - WP: `T_1 = <p>^2/2m` and `T_2 = <p>^2/2m + var(p)/2m` vs t.
     **NOTE THE LABEL SWAP vs `ks_stopping.py`** (T1<->T2 there); the notebook
     prints the mapping so no reader can conflate them.
   - a classical-vs-WP overlay.
3. **Momentum distribution at three times** — the direct picture of what the
   `var(p)` term is doing (broadening vs shifting vs growing a scattered tail).
4. **Interaction-energy deltas** — `Delta E_SS`, `Delta E_PS`, `Delta E_PP`,
   `Delta E_SB`, `Delta E_PB` for BOTH halves
   (`.claude/rules/decomposed-interaction-energies.md`).
5. **Window selection — DEFERRED to the user.** A single clearly-marked
   parameter cell (`T_WIN_CL`, `T_WIN_WP`, both `None` on first build). While
   unset the notebook prints guidance instead of fitting. Once set, the SAME
   fitting code in `channeling_stopping.py`/`ks_stopping.py` runs on the
   user-chosen windows, per-half.

**Validation:** `tests/test_refined.py` pins the algebraic identities the
notebook's claims rest on (T_2 - T_1 == var(p)/2m exactly; classical closure;
cumulative-trapezoid path recovers a known constant-velocity track; interaction
deltas are zero at t=0).
