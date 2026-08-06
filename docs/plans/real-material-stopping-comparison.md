# Plan: Two stopping-power definitions on a real material

**Created:** 2026-08-05. **Status:** material-selection investigation complete;
simulation campaign not yet designed in detail, nothing launched.

## Goal

Apply the project's two stopping-power definitions to a real material and
compare against experimental / analytical references:

1. **KS-orbital-dependent definition** — per-orbital energy bookkeeping
   (bulk-jellium KS stopping line of work; see
   docs/handovers/slab-ks-orbital-stopping-wrap.md).
2. **Slab definition** — ΔE_total / L_z across a finite slab
   (localised-jellium line of work).

Then discuss where the two definitions agree/disagree when the electron gas
is replaced by real bands.

## Material selection (the primary decision)

Candidates given by the user: bilayer graphene and "BNO3".
**"BNO3" is not a standard material** — no stopping-power literature exists
under that formula (verified by search, 2026-08-05). Interpreted as
**h-BN (hexagonal boron nitride)**, graphene's standard 2D companion.
(If B2O3 or a perovskite oxide was meant, note the IAEA data concentrate on
oxides — Al2O3/SiO2/Ta2O5/TiO2 — but those are 3D bulk targets, a poor fit
to the slab definition; say so and re-open only if the user confirms.)

### Decision: **bilayer graphene** (recommended)

Selection criterion #1 (user): availability of experimental evidence.
Graphene wins on every experimental channel; h-BN loses on all of them.

| Evidence channel | Bilayer graphene | h-BN |
|---|---|---|
| Electron projectile, OUR energy range (≈15–300 eV) | eV-TEM transmission + LEEM reflection through free-standing 1–4 layers, 0–25 eV: measured IMFP, layer-resolved (Geelen PRL 123, 086802 (2019)); LEEM layer-counting oscillations (n−1 minima, Hibino 2008); ~70% transparency at 50–250 eV (LEED holography); TOF IMFP (Nanomaterials 11, 2435 (2021)) | none found (only EELS on crystals) |
| Electron projectile, high energy | NIST ESTAR graphite, 1 keV–10 GeV (ICRU 37, re-evaluated) — high-energy anchor only | not in ESTAR default list; Bragg-rule estimate only |
| Ion projectile | carbon among the best-measured targets in the IAEA database (H/He, foils); HCI/heavy-ion transmission through free-standing 1–2 layer graphene with measured energy loss per layer (TU Wien: Nat. Commun. 7, 13948 (2016); PRA 93, 052708 (2016); Commun. Phys. 2019, 2021) | absent from IAEA compound highlights; irradiation studies measure defects, not energy loss |
| TDDFT theory baseline | Ojanperä PRB 89, 035120 (2014) graphitic targets; graphite anisotropy (arXiv:1905.07200); trajectory-sampling methodology (npj Comput. Mater. 9, 205 (2023)) | none found for stopping |
| Analytical reference | 2D dielectric / plasmon stopping models for graphene (Mišković group; RPA π-band dielectric function literature) | sparse |

Source notes: docs/sources/geelen-2019-evtem-graphene.md,
gruber-2016-hci-graphene.md, ojanpera-2014-graphite-stopping.md,
montanari-2024-iaea-database.md.

### Motivation (secondary criterion — writes itself)

- Bilayer graphene is the **thinnest possible slab**: the slab definition
  ΔE/L_z maps onto measured energy-loss-per-layer, and the mono→bi-layer
  step is an experimental L_z-knob matching our slab-thickness scaling runs.
- Graphene's π/σ valence electrons are quasi-free → the jellium picture is
  a controlled approximation; the KS-orbital definition gains real,
  physically distinct orbitals (π vs σ) to decompose over.
- Semi-metal (no gap) → screening is jellium-like. h-BN's ~6 eV gap breaks
  the jellium analogy — which makes h-BN the natural *second* material
  later ("what does a gap do to each definition?"), not the first.
- Novelty: no published layer-dependent TDDFT stopping study of bilayer
  graphene was found (2026-08-05 search; "not found", not "does not exist").

## Honest caveats (must appear in any write-up)

1. **ESTAR does not reach our regime.** ESTAR is semi-empirical Bethe-based
   *calculation*, valid ≥ 1 keV; our electron runs are 15–300 eV. Below
   1 keV the right observables are transmission/IMFP (eV-TEM) and
   ELF-derived stopping. Use ESTAR only as the ≥1 keV anchor.
2. **IMFP ≠ S.** Experiments give attenuation/transmission per layer;
   S is mean energy loss per path. Our runs can emit both (WP norm loss
   per layer ↔ transmission; ΔE per layer ↔ S), so compare each channel
   to its own experimental counterpart, never IMFP directly to S.
3. **HCI energy loss includes charge exchange** — compare our fixed-charge
   projectiles to low-charge-state / heavy-ion or proton-foil data.
4. **Trajectory dependence.** Real lattice → S depends on impact point
   (hollow/bridge/top) and incidence angle; runs must state and sample
   trajectories (unlike jellium where the trajectory is irrelevant).

## Campaign design v1 (2026-08-05 brainstorm; user confirmed bilayer graphene)

Decisions locked: slab ΔE/L_z definition is the PRIMARY aim; KS-orbital
definition extracted for free from the same runs; CAPs (absorbing
boundaries) in the box; electron projectile (WP + classical twin),
continuing the corpus. Working numbers: a = 2.46 Å = 4.6 Bohr; Bernal
interlayer d = 3.35 Å = 6.3 Bohr; v = √(2E): 15→300 eV ≙ v = 1.0→4.7 a.u.,
λ = 2π/v = 6.0→1.3 Bohr.

### Stage A — Ground state build + verification

- LDA (PZ81) + norm-conserving pseudopod carbon (consistent with jellium
  corpus). Extra justification to pin with a source note BEFORE runs
  (literature-check TODO): LDA fortuitously binds graphite layers near the
  experimental spacing; plain PBE does not.
- AB (Bernal) bilayer. **Supercell must be 3n×3n so the K point folds to
  Γ** (Γ-only RT-TDDFT would otherwise miss the Dirac-derived metallic
  states). Pilot tier: 3×3 (36 atoms, 144 valence e⁻, 72 doubly-occupied
  orbitals ≈ n162-jellium cost); convergence tier: 6×6 (144 atoms).
- Verification battery (Tier A, cheap): cutoff sweep 30→50 Ha to
  < 1 mHa/atom; vacuum sweep; a and d vs literature LDA (≈2.45/3.3 Å);
  Γ eigenvalue sanity — bilayer Dirac-state split γ₁ ≈ 0.4 eV; occupation/
  smearing choice recorded. GS saved to shared_gs/.

### Periodicity (user-locked, 2026-08-05)

The slab definition requires a CAP ⇒ the box is **periodic in-plane (xy)
only, non-periodic along z**. Implementation routes: (a) true 2D/slab
periodicity if INQ's cell supports it cleanly, or (b) the corpus route —
formally periodic z with a large vacuum gap in which the CAP absorbs all
flux before wraparound (z-periodicity physically moot). VERIFY against
docs/inq_tutorial.md which route INQ supports; if (b), add a Hartree
image-slab (dipole/monopole correction) check to the Stage A battery.

### Why an in-plane supercell is unavoidable (user question, answered)

NOT for bulk properties. Two reasons:
1. The in-plane grid is periodic, so the projectile is repeated in xy.
   A 1×1 primitive cell (4.6 Bohr period) would simulate an infinite
   array of simultaneous projectiles with overlapping wakes. The
   supercell separates the projectile from its images (14 Bohr @ 3×3,
   28 Bohr @ 6×6).
2. A 3n×3n supercell at Γ ≡ a 3n×3n k-mesh on the primitive cell — this
   IS the k-sampling, and folds the Dirac K states to Γ.
Cost check: 3×3 bilayer = 72 occupied orbitals < n162 jellium (81) —
standard corpus cost on one A100. Only the single 6×6 finite-size check
(~288 orbitals) is expensive.

### Stage B — Box geometry (z-budget, outside-in)

| Zone | Extent (Bohr) | Rationale |
|---|---|---|
| CAP top | 15–20 | ≥ ~3λ of slowest flux (λ = 6.0 at 15 eV) |
| launch + relax | 12–15 | WP starts ~12 Bohr above slab; wake decays pre-CAP |
| slab | 6.3 (+tails) | the bilayer |
| exit relax | 12–15 | transmitted WP + secondaries clear slab |
| CAP bottom | 15–20 | transmitted flux |

Total L_z ≈ 70–90 Bohr (familiar: σ=5/6 jellium slabs used L_z=105).
In-plane: 3×3 → 14 Bohr between periodic track images (tight for wake);
6×6 → 28 Bohr — the finite-size check (cf. Kononov npj Comput. Mater. 9,
205 (2023)).

### Stage C — Projectile / RT parameters

- Energy grid E = 15, 25, 50, 100, 200, 300 eV: 15–25 eV overlaps Geelen
  eV-TEM (layer-resolved transmission, direct comparison); 50–300 eV
  overlaps LEED transparency; expected stopping maximum near v matching
  the π+σ plasmon (~27 eV) response, around ~100 eV — grid resolves a
  peak, not a slope.
- **σ_WP = 2.0 Bohr — PRODUCTION (revised 2026-08-05, Phase 1 measurement).**
  The σ=0.5 CAP scan measured η-INDEPENDENT low-k residues ≈ 0.2: at k0=1.05,
  σ_p = 1/(√2·0.5) = 1.41 EXCEEDS k0 — ~20% of the packet has near-zero/backward
  p_z (never beam-like below ~250 eV), and the packet carries ½·3σ_p² ≈ 82 eV
  internal spread energy, swamping a 15 eV drift. At σ_WP=2.0: σ_p=0.354
  (k0=1.05 ≈ 3σ_p — clean), spread energy 5.1 eV (caveat at E=15–25 eV,
  recorded), and transverse spreading toward plane-wave character MATCHES the
  transmission-experiment geometry (Geelen uses plane-wave-like electrons).
  σ_WP=0.5 remains a high-E-only (≥250 eV) Phase 3 stretch variant. Classical
  twin at σ_pot = σ_WP/√2 per convention.
- Trajectories: NO through-hollow channel exists in AB stacking (hexagon
  centre of layer 1 sits atop an atom of layer 2). Pilot: hexagon-centre
  of top layer; then top-site and bridge; report site-resolved S before
  averaging.
- dt: cutoff-limited; take corpus-standard, verify energy conservation
  CAP-off before production (do NOT quote a dt until tested). Run length
  ~40–60 Bohr of flight → few thousand steps. Final-timestep checkpoint +
  resume mandatory (rule).

### Stage D — CAP design + energy bookkeeping

- Two-sided sin² CAP (`perturbations::absorbing`, inq-study engine).
  **VERDICT (Phase 1 scan, 2026-08-05, σ_WP=2.0): η = 1 Ha, W = 20 Bohr/face,
  L_z = 80 — ONE CAP for the whole campaign** (per-energy CAPs would break
  cross-energy CAP-ledger comparability; lower η leaks 12% at 300 eV).
  Measured vacuum residues: k=1.36: 6.3e-3 (documented norm-accounting
  uncertainty for E=25 eV; partly in-transit slow tail); k=1.92: 1.9e-4;
  k=2.71: 7.5e-7; k=3.83: 3.0e-5; k=4.70: 2.1e-4. k=1.05 (15 eV): ~3%
  reflection floor — E=15 CAP treatment is a Phase 3 decision.
  σ=0.5 scan (kept in capscan/results, no _s2 suffix) documents WHY σ=0.5
  is unusable below ~250 eV (η-independent ~20% lingering-weight floor).
- **CAP scan is over the ELECTRON k-range only (user question, 2026-08-05
  resolved):** the CAP acts on electronic orbitals; classical Ehrenfest
  ions do not interact with a CAP at all. Ion containment is kinematic:
  max binary-collision transfer electron→carbon is 4mM/(m+M)²·E ≈
  1.8e-4·E ≈ 0.05 eV at 300 eV, far below graphene's ~20 eV displacement
  threshold — no carbon reaches the CAP. Flux to absorb = transmitted WP
  (k ≤ 4.7), backscatter, slow secondaries (low-k = hard end of scan).
  A heavier-mass scan would be needed ONLY if the muon-mass-fork WP trick
  is imported (k = m·v shifts) — not in this campaign.
- **Bookkeeping subtlety:** with a CAP, ΔE_total(box) ≠ deposited energy;
  the slab observable is ΔE_slab = ΔE_total(box) + E_removed_by_CAP. Port
  the jellium m8/m9 loss-function CAP accounting — do not reinvent.
- **Ehrenfest ions ON (user decision, 2026-08-05):** consequences —
  (i) deposited energy splits: ΔE_slab = ΔE_electronic + ΔKE_ions +
  E_CAP_removed → ion kinetic energy is an explicit observables column
  (measured, expected small per kinematics above);
  (ii) Phase 0 gate: residual GS forces ≈ 0 (relaxed lattice) or ions
  drift spuriously during RT and pollute ΔKE_ions;
  (iii) rt_state.txt persists ALL ion positions/velocities (final-
  timestep-checkpoint rule), not just projectile state.
- decomposed-interaction-energies rule adaptation: real ions replace the
  jellium background group B; closure gates become lattice-external terms
  (same CSV schema; document adapted gates in run.cpp header).

### Definition mapping (the science deliverable)

- Slab definition: per-layer ΔE via bilayer−monolayer subtraction (the
  L_z-derivative) → monolayer runs are REQUIRED, not optional; transmission
  per layer at 15–25 eV compares directly to Geelen.
- KS-orbital definition: free from same runs; orbitals now have physical
  identity (π vs σ) → "which band absorbs the energy" — an axis jellium
  could not offer.

### Phase order (v2, after 2026-08-05 user feedback)

- **Phase R — reference freeze (DONE 2026-08-05, BLINDED):** experimental
  reference data verified extractable and frozen BEFORE run design:
  ESTAR graphite table (1 keV–1 GeV, parsed CSV), Geelen 2019 quoted
  anchors (λ_inel ≈ 3→1 layers over 0→25 eV; mono-vs-bi 5–15 eV
  signature), CXRO f1/f2 (ELF bridge route, caveated), proof figure.
  Pack: `docs/validation/reference_data_bilayer_graphene/`; protocol +
  inventory: `docs/validation/bilayer-graphene-stopping-reference.md`.
  **Do not consult during Phases 0–3; unbox in Phase 4 only.**
- **Phase 0 — GS build + verify:** mono + bilayer GS (Stage A battery).
- **Phase 1 — CAP reflection scan:** vacuum-only WP shots at each v into
  the CAP; fixes L_z, CAP width/strength (residue < 1e-3, v = 1.0–4.7).
- **Phase 2 — classical-twin sanity (EARLY, user-requested):** 3 twin
  pairs on the bilayer at E = 25, 100, 300 eV (grid ends + expected peak).
  User hypothesis: WP ≈ classical here ("expect no difference").
  CAVEAT (recorded): bulk jellium measured WP/classical ≈ 2.2
  (docs/handovers/bulk-jellium-ks-stopping.md); plausible it shrinks for
  a 6.3 Bohr slab (transit time 2–6 a.u. → little time for
  self-interaction wake feedback — Inference, unverified). Gate: ratio
  MEASURED and documented per energy, not asserted; ratio ≈ 1 would
  itself be a result (thin-slab suppression of the WP SIE excess) —
  either outcome informs Phase 3, neither aborts it.
- **Phase 3 — WP production sweep:** mono + bilayer × 6 energies ×
  sites (hexagon-centre first, then top/bridge).
- **Phase 4 — definitions + experiment comparison:** slab ΔE/layer and
  transmission vs Geelen (15–25 eV) / LEED transparency (50–250 eV);
  KS-orbital (π vs σ) decomposition; write-up.

### Run-output contract (user-locked, 2026-08-05 — applies to EVERY run)

| Output | Cadence | Notes |
|---|---|---|
| observables.csv + full energy decomposition | every step | includes explicit ΔKE_ions column (Ehrenfest) |
| interactions.csv (pairwise decomposition) | every step | rule-mandated; B group = ionic lattice (adapted closure gates) |
| CAP-removed energy + norm ledger | every step | needed for ΔE_slab = ΔE_box + E_CAP |
| WP KS-orbital momentum distribution |ψ̃_WP(k)|² | LOW cadence (~every 100 steps) | user: keep runs fast; inqview.analysis momentum machinery |
| VTI density frames | LOW cadence (~every 50–100 steps → ~30–60 frames) | user: manageable disk; still enough frames for the mandatory density GIF |
| Checkpoints | 2–3 per run: interior at ~1/3 and ~2/3 + FINAL at last step | user decision overrides the ~200-step interior default; final-timestep rule holds; rt_state persists all ion R/V + wp_idx |
| Run notebook | per run, auto-built post-run | run-notebook skill; density GIF mandatory |
| Phase notebook | per phase | notebook-making skill; synthesises the phase's runs against its gate |

### Phase 4 plot battery (designed to answer THE question:
how do the two definitions compare — to each other, and to the real world)

1. S(E) overlay: S_slab (ΔE/layer) vs S_KS (orbital definition), same
   runs, mono + bilayer — the definition-vs-definition plot.
2. WP/classical ratio vs E (Phase 2 twins) — representation dependence.
3. Transmission per layer vs Geelen λ_inel anchors (via T+R = e^(−d/λ_inel))
   — UNBOX reference pack here.
4. Mono-vs-bilayer 5–15 eV feature vs Geelen interlayer-bandgap signature.
5. S(E) high-E trend against the ESTAR graphite curve (≥ 1 keV anchor).
6. π vs σ orbital-resolved energy absorption (KS-definition-only axis).

### Open questions / user decision points (remaining)

1. Supercell tier for production: 3×3 throughout with one 6×6 check
   (recommended), or 6×6 production if cost allows?
2. INQ periodicity route (a) vs (b) — resolve from docs/inq_tutorial.md
   at Phase 0 implementation.

Resolved 2026-08-05: Ehrenfest ions ON (free C lattice) — user decision;
periodicity xy-only + absorbing CAP on z — confirmed; CAP scan is
electron-k only, no heavier-mass scan (see Stage D).

## Cross-references

- Handover: docs/handovers/real-material-stopping.md
- Prior definitions: docs/handovers/slab-ks-orbital-stopping-wrap.md,
  docs/handovers/wp-localised-jellium-solving-cap.md
