# Muon through electronic jellium — inq-study engineering notes

> Working notes built during a grill-with-docs session (2026-07-06).
> Purpose: capture *exactly* what `inq-study` would need so a muon can traverse
> electronic jellium, with open questions annotated inline (**Q:**). Companion to
> the campaign spec `docs/campaigns/muon_projectile/muon_projectile.md`.
> Nothing here is locked until moved to the campaign spec via `/campaigns`.

## 0. Grounded facts (verified in code this session)

- The WP projectile is **injected as the last extra-state orbital of the INQ
  `electrons` object** — `inqkit/wavepacket/wavepacket.hpp`
  (`inject_into_last_extra_state`). Projectile + jellium KS electrons share **one
  `orbital_set`, one Hamiltonian, one density**.
- The kinetic operator applies a **single scalar prefactor** to the whole set:
  `inq-study/src/hamiltonian/ks_hamiltonian.hpp`
  - `:202` `operations::laplacian(phi_fs, -0.5, …)` (real-space apply)
  - `:235` `operations::laplacian_add(phi, hphi, -0.5, …)` (fourier apply)
  - `:245` `operations::laplacian_expectation_value(phi, -0.5, …)` (kinetic energy)
  - `-0.5` = −ℏ²/2m with m = 1 (electron mass, atomic units).
- `operations::laplacian(...)` (`inq-study/src/operations/laplacian.hpp`) takes
  `FactorType factor` as a **scalar**; its GPU kernel multiplies every state
  index `ist` by the same `factor`. **No per-orbital hook exists today.**
- Classical heavy projectile needs **zero engine work** — ionic mass is a
  `run.cpp` parameter.

## 1. The root branch — is the muon quantum or classical?  (RESOLVED — Q1)

**Decision (2026-07-06): BOTH, for a comparison.**
- **Classical muon** = heavy point charge (mass 206.77 m_e, charge ±1) under
  Ehrenfest. **Zero engine work** — ionic mass is a `run.cpp` parameter. Serves
  as the reference baseline.
- **Quantum muon** = Gaussian WP orbital propagating with −ℏ²/2m_μ. **Requires the
  inq-study per-orbital-mass fork** (§2). This is the research contribution.
- The comparison (classical vs quantum muon) is the deliverable, mirroring the
  existing electron WP-vs-classical pairs.

## 2. inq-study per-orbital-mass fork — the specifics (Strategy A, recommended)

**Design in one line:** carry a per-state *inverse mass* (default 1.0) on
`electrons` — mirroring `occupations_` — read it in `ks_hamiltonian`, and apply
it as a per-state prefactor in the three `operations::laplacian*` kernels. The
ground state stays untouched because the constructor defaults to all-ones; only
the RT propagator opts in. Each site below carries a `// PROPOSAL:` marker in the
source.

| # | File (all under `inq-study/src/` unless noted) | Layer | Change (marked `PROPOSAL:` in file) |
|---|---|---|---|
| 1 | `systems/electrons.hpp` | Data / interface | Add `gpu::array<double,2> inverse_mass_` (shape `[kpin][local_state]`, mirrors `occupations_`); default-fill 1.0 at both `reextent` sites; add `inverse_mass()` accessors. This is the knob `run.cpp` writes. |
| 2 | `operations/laplacian.hpp` | Compute | Add **per-state-factor overloads** of `laplacian`, `laplacian_add`, `laplacian_expectation_value`; kernel `factor` → `factor_arr[ist]`. Keep the scalar versions for all other callers. (`laplacian_in_place` NOT on the muon path.) |
| 3 | `hamiltonian/ks_hamiltonian.hpp` | Router | Add member `inverse_mass_`; constructor gains a defaulted (all-ones) inverse-mass arg; call sites **`:202`** (real apply), **`:235`** (fourier apply), **`:245`** (energy ledger) pass the per-state factor `-0.5*inverse_mass_[ist]`. |
| 4 | `real_time/propagate.hpp` `:79` | Opt-in (RT) | Pass `electrons.inverse_mass()[0]` into the RT `ks_hamiltonian` ctor. Needed whenever mass affects **dynamics** (muon projectile, or all-muon RT). |
| 5 | `ground_state/calculator.hpp` `:97` | Opt-in (GS) | Pass `electrons.inverse_mass()[0]` into the SCF `ks_hamiltonian`. **Needed for ALL-MUON jellium + band-structure** (bath mass is in the GS). NOT needed for electron-jellium + RT-injected projectile. |
| — | `ground_state/initial_guess.hpp` | **No change** | Builds a Hamiltonian only for the mass-independent overlap operator (orthogonalisation); SCF re-converges with the correct mass regardless. Annotated "NO CHANGE" in-file. |
| 6 | `ResearchProject/.../<muon_wp>/run.cpp` *(consumer, NOT inq-study)* | Consumer | Set `electrons.inverse_mass()`: for a muon **projectile** → only `[0][wp_idx]=1/M_MU`; for an all-**muon jellium** → all entries `=1/M_MU` (bath + projectile). Classical-muon runs never touch this. |

**GS-affected rule:** mass enters the GS *iff* the particle whose mass changed is
present in the GS. Projectile injected at RT → GS untouched (file 5 stays default).
All-muon jellium / band structure → GS carries the mass (file 5 opts in). This is
why the design threads `inverse_mass` through BOTH the GS (`calculator.hpp`) and RT
(`propagate.hpp`) construction sites, each defaulting to all-ones.

### Correctness notes / caveats baked into the annotations
- **Energy ledger is consistent by construction:** the same per-state factor used
  in the apply kernels is reused at `:245` → `energy.hpp:82` sums it (occupation-
  weighted) into `E_kinetic`. This is the muon route's *clean-ledger* payoff.
- **Robustness across the propagator:** keying the mass off `phi.set_part()`
  global indices (on the Hamiltonian) is immune to the ephemeral-`orbital_set`
  problem — internal ETRS/CN temporaries preserve the partition, so the muon slot
  never silently reverts to `m_e`. (This is why we did NOT store mass on
  `orbital_set` — Strategy A1, rejected.)
- **Vector-potential caveat:** the `+A²/2m` term added via `scalar_potential_add`
  is NOT mass-scaled by this fork. It is **zero at gamma / zero vector potential**,
  which is the WP-muon regime (momentum lives in the orbital phase `exp(ik₀·r)`,
  not in a vector potential). A vector-potential *kick* would need extra work.
- **Checkpoint (optional):** `inverse_mass_` need not be serialised for a fresh RT
  run (run.cpp sets it). Only RT-restart-from-mid-run would need save/load.

### Validation gates the fork must pass (before any quantum-muon run)
- **Free-muon dispersion:** ⟨T⟩ = ħ²k₀²/2m_μ for a plane-wave/WP in vacuum at the
  muon slot (formula-validation).
- **Spreading rate:** τ_spread ∝ m — a muon WP spreads ~207× slower than an
  electron WP at the same σ.
- **Electrons unchanged:** with the muon slot present, the electron orbitals'
  kinetic energies and the GS are bit-for-bit identical to the unforked engine.
- Add a catalogue row (`docs/validation/test-catalogue.md`).

## 3. XC / Hartree / functional considerations

### The core problem (muon-through-ELECTRON-jellium)
The muon WP is one orbital in the shared `orbital_set`, so |ψ_μ|² enters the ONE
total density that feeds `V_H[n]` and `V_xc[n]`:
- **Hartree:** muon↔electron Coulomb is correct physics; muon self-Hartree is
  spurious (SIE) — same SIE the electron-WP campaign already brackets.
- **XC:** the muon is swept into the electrons' LDA. A muon is a *distinguishable
  heavy lepton* → it has **no exchange/correlation-as-an-electron** with the bath
  and no XC self-energy. For the muon this is *purely* spurious (worse than the
  electron-WP case, where the projectile at least IS an electron).
- **"Different functional?"** No — you don't switch XC *family* (LDA stays right
  for the jellium electrons). The issue is that the muon shouldn't be *inside* the
  electron XC / self-Hartree term at all.

### Option (b) — RECORDED FOR DELIBERATION (user, 2026-07-06): muon in its own spin channel
Give the muon WP its **own spin channel** so it is a distinguishable species: it
then has no same-channel exchange with the jellium electrons. This mimics
distinguishability at the exchange level. Caveats to work through later:
- Removes cross-species *exchange*, but **not** self-Hartree SIE.
- INQ spin channels are electron spin, not species labels — using a spin channel
  as a species tag is a modelling choice to validate (does the LDA correlation
  between channels then misrepresent muon↔electron correlation?).
- Engine cost: spin-polarised setup + ensuring the muon channel carries muon mass
  (still the §2 per-orbital fork). **Deliberation open — not decided.**

### Option (c) — user's simpler starting system: muon projectile through MUON jellium
**Idea (user, 2026-07-06):** make *every* particle a muon — bath + projectile.
Then projectile and bath are the **same species → indistinguishable → exchange is
legitimate**, and the distinguishability artifact is escaped entirely. Start here
as the simplest *consistent* system, before tackling the harder electron-jellium
(distinguishable) case.

Consequences to weigh (these are the live grill points):
1. **Engine work: SAME `inverse_mass_` fork (§2), set GLOBALLY.** Muon jellium =
   set `inverse_mass_` to `1/m_μ` for *all* orbitals instead of only the WP slot.
   No new machinery beyond §2. BUT —
2. **The ground state must ALSO be muon-mass.** Unlike electron-jellium+muon-
   projectile (GS = pure electron jellium, muon injected only at RT, GS pristine),
   a muon jellium GS is itself a muon HEG. So `ground_state/calculator.hpp` and
   `ground_state/initial_guess.hpp` must ALSO receive the inverse mass — the
   "opt-in only in propagate" simplification of §2 no longer holds. More files,
   but same array.
3. **XC parameterisation (the subtle one).** *Inference (verify via
   literature-review):* a homogeneous gas of identical mass-m charge-e fermions is,
   in **effective atomic units** (a₀\* = ħ²/m e², Ha\* = m e⁴/ħ²), identical to the
   electron gas at the same r_s — this is why r_s is THE parameter. But INQ runs in
   *real* atomic units (m_e ≡ 1), and libxc's LDA is electron-parameterised. So a
   muon jellium at a given *physical* density, run with electron-LDA as-is, is a
   **consistent fictitious "heavy-fermion" gas**, not a quantitatively exact muon
   gas, unless XC is mass-rescaled. For a *controlled quantum-vs-classical
   comparison* the XC error likely **cancels in the difference** (both share the
   same bath), so option (c) may be perfectly adequate as a first system — TBC.

**Staging that emerges:** (c) muon-jellium [simplest, escapes distinguishability,
global mass] → later (a)/(b) electron-jellium [per-orbital fork + distinguishability
handling]. The §2 fork serves ALL of these — it's a superset.

### CRITICAL rescaling fact (2026-07-06) — reframes (c) vs (a)/(b)
A one-component Coulomb system rescales EXACTLY with mass. In effective units
(a* = ħ²/m e², E* = m e⁴/ħ²) the muon-jellium Hamiltonian is the electron-jellium
Hamiltonian; the bracket depends only on r_s (textbook HEG scaling, Giuliani &
Vignale). Therefore, at fixed dimensionless r_s:
- length ∝ 1/m, energy ∝ m, time ∝ 1/m, **velocity is mass-invariant** (a*/t* = a₀/t₀).
- **Stopping power: S ∝ m² in real units, IDENTICAL in effective units.**

**Implication:** muon-through-muon-jellium at fixed r_s = the electron problem in
disguise → NO new stopping physics (S_muon = m²·S_electron). And the muon's whole
*motivation* (slower spreading / cleaner ledger *than an electron WP*) only exists
when projectile mass ≠ bath mass. So:
- **All-muon jellium (c) = VALIDATION phase** — its value is precisely that it MUST
  reproduce S_electron·m², an end-to-end correctness test of the global-mass fork.
  At fixed r_s, electron-LDA is EXACT (same dimensionless r_s), so the Q4 XC worry
  dissolves — no rescaling needed if r_s is held fixed.
- **Muon-in-electron-jellium (a)/(b) = PHYSICS phase** — the only setting with novel
  physics; unavoidably distinguishable → needs SIE-bound (a) or spin-channel (b).

## 5. Call chain from run.cpp → kinetic operator (answers "which file calls which")

Verified in code 2026-07-06. Two chains share the same Hamiltonian + kinetic path;
they differ only in who constructs the Hamiltonian and how it is applied.

### Data flow of the mass (the short version)
```
run.cpp:  electrons.inverse_mass()[0][ist] = 1/m        (systems/electrons.hpp — stored)
   │
   ▼  passed into the ks_hamiltonian constructor
ks_hamiltonian.inverse_mass_   (hamiltonian/ks_hamiltonian.hpp — member)
   │
   ▼  read inside operator()(phi), built into factor_arr = -0.5 * inverse_mass_
operations::laplacian*(phi_fs, factor_arr, …)   (operations/laplacian.hpp — applied)
   │
   ▼  GPU kernel:  lapl = factor_arr[ist] * (-g² + …)
```

### A. Ground-state chain (band-structure demo, muon jellium — mass affects the GS)
```
run.cpp
  → ground_state::calculate(...)                         [ground_state/calculator.hpp]
      → constructs hamiltonian::ks_hamiltonian(basis, states, pot, ions, inverse_mass)
      → SCF eigensolver applies ham(phi) many times:
          ks_hamiltonian::operator()(phi)                [hamiltonian/ks_hamiltonian.hpp:194/218]
            → KINETIC:  operations::laplacian(phi_fs, factor_arr, …)  [operations/laplacian.hpp]   (in FOURIER space)
            → POTENTIAL: hamiltonian::scalar_potential_add(…)                                       (in REAL space)
            → NON-LOCAL: projectors_all_.apply(…)
```
(For the muon-projectile case the GS is pure electron jellium and is instead just
`electrons.load(GS_DIR)` — the muon appears only in the RT chain below.)

### B. Real-time chain (muon projectile — mass affects dynamics)
```
run.cpp:199  real_time::propagate(ions, electrons, step_fn, theory.lda(), rt_opts, pert)
  → real_time::propagate(...)                            [real_time/propagate.hpp:39]
      → constructs ham(..., electrons.inverse_mass()[0]) [propagate.hpp:79]
      → per step, calls the propagator:
          etrs(...)            [real_time/etrs.hpp:23]      OR
          crank_nicolson(...)  [real_time/crank_nicolson.hpp]
            → operations::exponential_2_for_1(ham, …) / exponential_in_place(ham, …)  [operations/exponential.hpp]
                → applies ham(phi) repeatedly (Taylor series of exp(-iH dt)):
                    ks_hamiltonian::operator()(phi)       [ks_hamiltonian.hpp:194/218]
                      → KINETIC:  operations::laplacian / laplacian_add (factor_arr)   [laplacian.hpp]
                      → POTENTIAL: scalar_potential_add(...)
                      → NON-LOCAL: projectors
      → energy ledger each step: energy.calculate(ham, electrons)   [hamiltonian/energy.hpp:82]
            → ham.kinetic_expectation_value(phi)          [ks_hamiltonian.hpp:245]
                → operations::laplacian_expectation_value(phi, factor_arr, …)  [laplacian.hpp]
```

## 6. Corrections to the `MY THOUGHTS:` comments (2026-07-06)

**(ks_hamiltonian.hpp) "kinetic = laplacian in real space; Coulomb = laplacian_add
in fourier space" — INCORRECT, on two counts:**
1. **`laplacian` and `laplacian_add` are BOTH the kinetic operator** (-½∇²).
   Neither is the Coulomb/potential term. They differ only by call style:
   `laplacian` RETURNS a new field (used when the input ψ is in real space —
   `operator()` at :194); `laplacian_add` ADDS in place to `hphi` (used when the
   input ψ is already in fourier space — `operator()` at :218).
2. **The kinetic term is applied in FOURIER space in BOTH cases** (∇² is diagonal
   there: multiply by -g² — note the `static_assert(... fourier_space ...)` in both
   functions). The **Coulomb / local potential** (Hartree + XC + local pseudo, i.e.
   the scalar potential V(r)) is applied in **REAL space** by
   `hamiltonian::scalar_potential_add` (:206/:226), because V(r)·ψ(r) is diagonal in
   real space. This dual-space split (kinetic in G-space, potential in R-space, FFT
   between) is the standard plane-wave pseudopotential method. → your guess had the
   two spaces swapped.

**(ks_hamiltonian.hpp) "adding an additional parameter to the functions that
specifies their mass" — CORRECT.** Yes: a per-state factor (= -0.5·inverse_mass)
replaces the scalar -0.5.

**(laplacian.hpp) "the laplacian is applied here using emplace" — INCORRECT.** It is
applied via `gpu::run(...)` with a `GPU_LAMBDA` (a GPU parallel kernel launch), not
`emplace`. (`emplace_back` appears elsewhere — building projector lists in
ks_hamiltonian — unrelated to applying the laplacian.)

**(laplacian.hpp) "multiply each of the states ist by the inverse mass before the
Laplacian is applied" — PARTLY RIGHT, one important correction.** Do NOT multiply
the STATE (the wavefunction ψ) by 1/m. Multiply the **prefactor** by 1/m, per state,
inside the kernel: `factor_arr[ist]*(-g² + …)`. Scaling the prefactor is
mathematically equivalent (the operator is linear) but scaling ψ itself would
corrupt its normalisation and every other term that reads ψ (potential, projectors,
density). So: per-state factor, not per-state wavefunction.

**(laplacian.hpp) "what other changes are made here / what purpose" — none beyond
the prefactor.** The ONLY change in this file is scalar `factor` → per-state
`factor_arr[ist]`, in three functions (as overloads). The single purpose is: give
each orbital its own -ħ²/2m. The "pre-factor" IS the (inverse) mass. Nothing else
in this file changes.

**(propagate.hpp) "pass inverse mass through the pipeline; ensure the mass change
persists" — CORRECT, with the mechanism spelled out.** Persistence is guaranteed
because (1) the Hamiltonian stores `inverse_mass_` once at construction and (2) it is
keyed by the **state-partition index**, so it survives every temporary `orbital_set`
the exponential/CN spawns internally. **BUT** — now that you want mass to be
*generalisable* AND to drive band structure (a ground-state property), the earlier
"opt in only at propagate, leave GS untouched" is too narrow: the GS calculator must
ALSO receive `inverse_mass` whenever the *bath* mass changes (muon jellium,
band-structure demo). Thread `electrons.inverse_mass()` through ALL ks_hamiltonian
construction sites, defaulting to all-ones. "GS untouched" now holds ONLY for the
special case of an electron bath with a heavy projectile injected at RT.

## 4. Open questions log

- **Q1 (RESOLVED):** Both classical + quantum muon, for comparison.
- **Q2 (RESOLVED):** Strategy A — per-state inverse mass on `electrons`; 4-file fork.
- **Design update:** mass must be GENERALISABLE (tunable per orbital, any value) —
  drives both the muon fork and the new mass-tuned-bands campaign.
- **Q3 (DELIBERATION OPEN):** muon-in-electron-jellium XC/distinguishability —
  options (a) accept+bound SIE, (b) own spin channel, (c) all-muon jellium escape.
- **Q4 (OPEN):** muon-jellium XC — electron-LDA as-is vs mass-rescaled; same r_s vs
  same physical density.
- **New:** mass-tuned-bands campaign drafted (`docs/campaigns/mass_tuned_bands/`).
