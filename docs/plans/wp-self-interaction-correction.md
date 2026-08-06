# Plan: removing the projectile self-interaction from KS-orbital stopping

**Status:** REVIEWED + AMENDED 2026-08-02 — implementation proceeding
(user directive: review against literature, amend, implement, validate,
run notebook, autonomous CSD3 submission).
**Original design:** 2026-08-02 (kept below; amendments marked ⟦AMENDED⟧;
the review is §0).
**Motivated by:** `docs/handovers/cylindrical-channeling-ks-stopping.md`,
2026-08-02 (5th) — the T1 deficit is dominated by LDA self-interaction error.
**Sources:** `docs/sources/messud-2008-tdsic.md`,
`docs/sources/perdew-zunger-1981-sic.md`,
`docs/sources/mundt-2007-zero-force-tdkli.md`,
`docs/sources/nazarov-2025-quantum-projectile-stopping.md`.

---

## 0. REVIEW (2026-08-02) — what the literature check found

The plan was reviewed against the time-dependent-SIC literature (Perdew–Zunger
1981; Messud, Dinh, Reinhard, Suraud PRL 101, 096404 (2008); Mundt, Kümmel,
van Leeuwen, Reinhard PRA 75, 050501(R) (2007); Legrand, Suraud, Reinhard
J. Phys. B 35, 1115 (2002) [ADSIC]; Ullrich, Gossmann, Gross PRL 74, 872
(1995) [TDOEP-KLI]) and against the actual INQ code paths. Verdict: **the
core design is sound and is a recognised member of the family of simplified
TD-SIC schemes, but the plan as written had two defects that would have made
its own validation tier fail, and one gate that is theoretically impossible to
satisfy as stated.** All are fixed below.

### D1 (MAJOR, would have failed Tier V): the xc subtraction was spin-inconsistent

The plan specified `v_xc[n_wp, 0]` "fully spin-polarised" (canonical PZ). But
INQ's production and vacuum runs are spin-RESTRICTED (`options::theory{}.lda()`
= libxc `XC_LDA_X` + `XC_LDA_C_PZ`, evaluated UNPOLARISED on the total
density). In the one-electron vacuum tier the Hamiltonian therefore contains
`v_H[n_wp] + v_xc^unpol[n_wp]`. Subtracting the POLARISED `v_xc[n_wp,0]`
leaves a residual `v_xc^pol − v_xc^unpol` — for LDA exchange alone a factor
`2^{1/3}`, i.e. a ~26 % remnant of the exchange potential — and the
free-dispersion gates (σ_z to 0.5 %) would FAIL for a *correct* implementation
of the *specified* scheme. **Fix: the xc self-term is evaluated with the SAME
spin treatment as the run** (unpolarised libxc `LDA_X` + `LDA_C_PZ` on
`n_wp`), through INQ's own `hamiltonian::xc_term::evaluate_functional`, so the
vacuum cancellation is exact by construction. Variants renamed:

| variant | subtracts from the WP's Hamiltonian | vacuum expectation |
|---|---|---|
| **none** | — | over-fast spreading (the artefact, isolated) |
| **SIC-H** | `Q v_H[n_wp] Q` | *under*-spreads: the attractive `v_xc^unpol[n_wp]` remains and self-binds |
| **SIC-PZrun** | `Q (v_H[n_wp] + v_xc^unpol[n_wp]) Q` | exact free dispersion — the closed-form gates apply to THIS variant only |

Canonical polarised PZ (SIC-PZpol) is noted for completeness but NOT
implemented: in a spin-restricted run it corresponds to no term actually
present in the Hamiltonian. (See `docs/sources/perdew-zunger-1981-sic.md`.)

### D2 (MAJOR, wrong gate): "E_corrected conservation is the primary correctness gate" is not achievable in the jellium, and the literature says so

Messud et al. prove that the only schemes conserving energy + zero-force +
orthonormality simultaneously are the full variational ones (Lagrange
multipliers on BOTH sides + the symmetry condition
`⟨ψ_β|U_β − U_α|ψ_α⟩ = 0`, propagated with a double set). Our projected kick
`V = QvQ` is exactly the ONE-SIDED Lagrange-multiplier scheme: it preserves
orthonormality EXACTLY (that is what Q is for) but is not the variational
scheme, so `E_corrected` has a residual drift channel

```
dE_corr/dt = 2 Im Σ_j ⟨ψ_wp| ĥ_LDA |ψ_j⟩ ⟨ψ_j| v_SIC |ψ_wp⟩   (our derivation)
```

which vanishes in vacuum (no occupied ψ_j) but not in the jellium. The
TD-KLI precedent (Mundt et al. 2007: zero-force violation → unphysical
self-excitation) shows such residuals can grow secularly. **Fix:**
`E_corrected` is exactly conserved as a GATE only in Tier V (vacuum). In
Tier B and production it is a MEASURED diagnostic with two conditions:
(i) |drift| over the run ≪ E_PP(0) = 1.94 eV (soft gate: < 5 % of it, i.e.
< 0.1 eV, calibrated in Tier B); (ii) NON-SECULAR — no steady linear growth.
A violation is reported, not silently absorbed.

### D3 (defect in the decision rule): SIC-H and SIC-PZrun are *expected* to differ in vacuum

The original rule "if indistinguishable use SIC-H" was written as if the xc
self-term might be negligible; with the vacuum Hamiltonian above they MUST
differ (`v_xc^unpol[n_wp]` is not small for a σ = 4 packet). The vacuum tier
separates the Hartree and xc parts of the self-interaction rather than
providing a degenerate check. Amended decision rule in §4.

### D4 (limitation now stated): one-orbital SIC does not remove every representation asymmetry

Even a perfect correction leaves the bath feeling `v_xc[n_S + n_wp]` while the
classical twin's bath feels `v_xc[n_S]` (its projectile is not in `n`). This
is arguably physics (a real electron xc-couples to the medium), but it means
the corrected WP run is still not term-by-term identical to the classical
ledger. Recorded as an interpretation caveat for the notebook, not corrected.

### D5 (limitation now stated): tagged-orbital identity and unitary non-invariance

PZ-SIC is not invariant under unitary mixing of occupied orbitals (Messud
Eq. 3e discussion); "correct only the projectile" presupposes the projectile
stays one KS column. INQ's RT propagators evolve columns independently (no
occupied-subspace rotation during propagation), so the tagged identity is
preserved by construction; this assumption is stated, and it is one more
reason the full double-set scheme (which deliberately mixes orbitals) is the
WRONG tool here — it would smear the correction across projectile and bath.
ADSIC and TDOEP-KLI are likewise rejected because they modify EVERY orbital's
Hamiltonian, breaking parity with the classical twin's plain-LDA bath.

### D6 (validation gap closed): Tier V never exercises the projection

In vacuum there are no occupied bath states, so Q = 1 and Tier V validates
ONLY the kick + bookkeeping. The projection machinery is validated in Tier B
by its own gates (max overlap, norm removed, bath orthonormality spot-check)
plus the engine test with a synthetic occupied manifold.

### D7 (stale): the §7 disk blocker is resolved

Measured 2026-08-02: /rds at 88 %, **127 GB free**. The production run
(~11 GB) fits with an order of magnitude of headroom. No user decision needed.

### D8 (expectation calibrated by new literature): a residual deficit is EXPECTED

Nazarov & Gross (arXiv:2510.26222, 2025) show from exact-factorization theory
that a finite-mass quantum projectile GENUINELY stops differently from a
classical point charge of the same charge and velocity. So the corrected run
does not test "does the WP become classical" but "how does the 20 % deficit
split between SIE artefact and genuine quantum kinematics". A nonzero
post-SIC deficit is a result, not a failure.

### Also fixed in passing

- `E_corrected` is variant-dependent: `E_KS − U[n_wp]` (SIC-H),
  `E_KS − U[n_wp] − E_xc^unpol[n_wp]` (SIC-PZrun).
- Strang bookkeeping across RESUME: `rt_state.txt` records whether the saved
  state is at a full-kick or closing-half-kick boundary (`sic_boundary=`);
  the resume branch applies the compensating half-kick when needed. Interior
  checkpoints sit at full-kick boundaries and resume with no action.
- Observable timing convention: all per-step observables are written at the
  post-propagator, PRE-kick state (uniformly), so WP curves stay directly
  comparable with the uncorrected run and the classical twin. The single
  endpoint half-kick after the final observable write is O(dt·v_SIC) and
  documented.

---

## 1. The problem, in one paragraph

The wavepacket projectile is an occupied Kohn–Sham orbital, so its own charge
density `n_wp` enters the Hartree potential `v_H[n]` that the same orbital then
feels. In exact DFT a one-electron density has **zero** self-interaction — the
Hartree self-repulsion is cancelled identically by exchange. LDA does not
cancel it. Measured on the production channeling run: `E_PP = 1.936 eV` at
`t = 0`, decaying to `0.292 eV` as the packet is blown apart by its own field.
The packet expands **1.47×** faster than free dispersion, and its measured
stopping power sits ~20 % below its classical twin. Whether the first causes the
second is the hypothesis set out in §5, not something established here.

**The classical twin needs no correction and is already the right reference.**
Its Gaussian charge is an *external potential*: it acts on the bath electrons
but never on itself, and its Ehrenfest force comes from `phi_drag`, which
excludes its own field. That asymmetry is the reason the twin comparison isolates
this term at all: at `t = 0` the two halves agree to 8e-12 Ha on every pairwise
energy, so self-interaction is the one thing that differs between them.

---

## 2. The correction  ⟦AMENDED per D1/D2⟧

Perdew–Zunger-type SIC **restricted to the projectile orbital**, with the xc
self-term evaluated run-consistently (unpolarised). The WP evolves under

```
H_wp  =  H_KS  -  Q ( v_H[n_wp]  [+ v_xc^unpol[n_wp]] ) Q
```

while every bath orbital keeps the unmodified `H_KS`. Variants **SIC-H** and
**SIC-PZrun** as defined in §0/D1; the vacuum tier discriminates.

`n_wp` is the occupation-weighted orbital density `occ_wp · |ψ_wp|²` with
`occ_wp = 1` asserted at start-up (the PZ one-electron logic presupposes one
electron in the orbital).

**`n_wp` stays in the total density.** The bath *must* keep feeling the
projectile's field — that coupling IS the stopping interaction. Only the
projectile stops feeling its own.

### Why this needs no engine modification

Both subtracted terms are local multiplicative potentials, so the correction
is applied as a per-step kick on column `ist_wp` of `phi.matrix()` inside the
`real_time::propagate` callback (the established wrapper pattern — cf. the
mask absorber). INQ calls the callback after every full ETRS step with mutable
`electrons`; a real multiplicative kick does not change `|ψ_wp|²`, so the
Hamiltonian INQ rebuilds from the density at the next step is unaffected by
the kick itself (the projection changes the density only at O((dt·v)²),
logged).

| need | existing API |
|---|---|
| mutable orbital access | `electrons.kpin()[0].matrix()` (as `wavepacket.hpp` writes) |
| `v_H[n_wp]` | `inq::solvers::poisson::solve(n_wp)` (as `interaction_energies.hpp` calls) |
| `n_wp` | `inqkit::jellium::orbital_density_field(electrons, wp_idx)` |
| `v_xc^unpol[n_wp]`, `E_xc^unpol[n_wp]` | `inq::hamiltonian::xc_functional{XC_LDA_X,1}` + `{XC_LDA_C_PZ,1}` via the static `hamiltonian::xc_term::evaluate_functional` (INQ's own libxc path — cancellation exact by construction) |

**`inq/` is not touched, and `inq-study/` is not needed.** New code lives in
`inq-stack/include/inqkit/wavepacket/self_interaction_correction.hpp`.

**Cost:** one extra Poisson solve + one LDA evaluation + ~80 overlap/axpy pairs
per step — small against the propagator's 104 per-orbital FFT sets. The
existing `interactions.csv` already does two Poisson solves every step.

### Orthogonality — the constraint that decides the design  ⟦VERIFIED against Messud 2008⟧

In ordinary TDDFT every orbital obeys the SAME `H_KS`, so the overlap matrix is
conserved exactly. An orbital-dependent Hamiltonian breaks that:

```
d/dt <psi_wp | psi_i>  =  -i <psi_wp| V |psi_i>
```

With a plain `V = v_H[n_wp]` this element is non-zero wherever packet and bath
overlap — and in this system it grows (`f_bore` 0.998 → 0.457). Losing
orthogonality double-counts density and violates Pauli exclusion. The
uncorrected runs have NO orthogonality error; a plain SIC would trade a known
artefact for a new one.

**Fix: project the correction into the unoccupied subspace.**

```
Q  =  1  -  sum_{j occupied} |psi_j><psi_j|          (Q psi_j = 0)
V  =  Q  v_SIC  Q                                     (Hermitian)
```

Then `V|ψ_i> = 0` for every occupied i, the leak rate is identically zero, and
because the injector enforces `ψ_wp ⊥ {ψ_j}` at t = 0, `Qψ_wp = ψ_wp` and
`<ψ_wp|V|ψ_wp> = <ψ_wp|v_SIC|ψ_wp>` — the projected operator removes the same
self-interaction energy. Nothing is given up.

**Literature status of this construction (review outcome):** it is exactly the
one-sided Lagrange-multiplier scheme — the WP equation keeps its multiplier
terms `Σ_j |ψ_j⟩⟨ψ_j|v|ψ_wp⟩`, the bath equations drop theirs. Messud et al.'s
variational TDSIC would add the bath-side terms and the symmetry condition; we
deliberately do not, because (i) orthonormality is already exact one-sided,
(ii) the bath Hamiltonian must stay IDENTICAL to the classical twin's, and
(iii) the double-set unitary mixing would destroy the tagged projectile
column. The price, per D2: `E_corrected` is conserved only approximately in
the jellium — measured, gated softly, and reported.

There is a physical reading as well: a projectile electron added to a filled
Fermi sea must scatter into states ABOVE it; Q enforces Pauli blocking.

**Cost of the projection:** ~80 inner products + 80 axpys per step (the same
GPU reduction pattern as the injector's Gram–Schmidt), and the overlaps it
computes are themselves the leak diagnostic — logged for free.

**Practical scheme** (kick, then project, then renormalise):

```
psi_wp  <-  N . Q . exp(+i dt_eff v_SIC) psi_wp
```

which agrees with `exp(+i dt QvQ)` to O(dt²) per step given `ψ_wp = Qψ_wp`.
The norm removed by Q each step is logged (`norm_removed`), and its cumulative
sum is a gate: if it is not ≪ 1 the correction is doing violence to the state
and the run says so, rather than absorbing it into the renormalisation.

**Kick scheduling (Strang):** opening half-kick at the t = 0 callback, full
kick after every interior step, closing half-kick at the final step. Interior
checkpoints therefore sit at full-kick boundaries (resume needs no action);
the final checkpoint sits after the closing half-kick and `rt_state.txt`
records `sic_boundary=closed` so an extension run applies the compensating
half-kick on load. Splitting error O(dt²) globally, matching ETRS at dt = 0.02.

### Energy bookkeeping  ⟦AMENDED per D2⟧

The corrected conserved functional is variant-dependent:

```
SIC-H:      E_corrected = E_KS[n] - U[n_wp]
SIC-PZrun:  E_corrected = E_KS[n] - U[n_wp] - E_xc^unpol[n_wp]
```

Written every step to `sic.csv` together with `u_self = U[n_wp]`, `exc_self`,
`max_overlap_pre`, `norm_removed`, `cum_norm_removed`. Gates: exact
conservation (< 1e-5 eV drift) in VACUUM only; in the jellium, soft gate
|drift| < 0.1 eV over the run AND non-secular (D2). `E_PP` in
`interactions.csv` remains, now a diagnostic of packet size, not a term the
WP's Hamiltonian contains.

---

## 3. Alternatives considered and rejected  ⟦extended by the review⟧

| option | why not |
|---|---|
| Exact exchange for the WP orbital | identical to SIC-H for a single orbital, far costlier in RT |
| Larger `sigma_WP` | reduces (~1/σ) but never removes the artefact; degrades wake resolution. Kept as the independent cross-check (§6) |
| Heavier projectile | changes the physics being studied |
| Separate spin channel | the Hartree term is built from the total density either way |
| **ADSIC** (Legrand 2002) | corrects EVERY orbital → bath Hamiltonian no longer matches the classical twin's plain LDA; wrong tool for a single-orbital artefact |
| **TDOEP-KLI** (Ullrich 1995) | common-potential approximation; violates zero-force + energy conservation (Mundt 2007) while ALSO not being exact for the one orbital we care about |
| **Full variational TDSIC / double-set** (Messud 2008) | conserves everything but mixes the WP with bath orbitals via the symmetrizing unitary — destroys the tagged projectile identity the measurement rests on |

---

## 4. Validation — vacuum first, where the answer is closed-form  ⟦AMENDED per D1/D3/D6⟧

**Tier V — one electron, empty box, no jellium, no background.**
Implemented by EXTENDING the existing `vacuum/scripts/wp_selfinteraction`
difference measurement (built earlier the same day; its smokes passed and its
`noninteracting` reference production run is complete): `sigma_WP = 4`,
**k0 = 0 (stationary — spreading is frame-independent and the box then never
wraps)**, spacing 0.5 (production), box 72³, `theory{}.lda()`, dt = 0.02, 1500
steps. Five runs in one job: noninteracting / hartree / lda (the uncorrected
difference triplet) + **(b) lda+SIC-H, (c) lda+SIC-PZrun.**

| quantity | exact free value | gated for | why |
|---|---|---|---|
| `var(p_z)` | constant `1/(2σ²) = 0.03125` | (c) only | conserved under free evolution — sharpest discriminator |
| `sigma_z(t)` | `sqrt(σ²/2 + t²/(2σ²))` | (c) only | direct over/under-spreading measure |
| `<p_z>` | constant `k0` (= 0 here) | (b) and (c) | net self-force must vanish — the zero-force check (Mundt) |
| `E_corrected` | constant | (b) and (c) | split-operator bookkeeping; exact in vacuum (Q = 1) |

Expected qualitative ordering (reported, not gated): (a) spreads FASTER than
free (Hartree self-repulsion), (b) spreads SLOWER than free (xc self-binding
remains), (c) free. The ordering is itself a three-point consistency check
that each subtracted term does what its sign says.

**Decision rule for the production variant:** SIC-PZrun **iff** it passes its
closed-form gates (var(p_z) drift < 0.1 %; σ_z within 0.5 % of analytic;
<p_z> drift < 0.01 %; E_corrected drift < 1e-5 eV). SIC-H is not a candidate
for the production headline (vacuum run (b) measures the xc self-term it
leaves behind) but stays available as an env-switch sensitivity check. If (c)
fails its gates the implementation is wrong: stop, diagnose, no production.

Cost: one electron, one state, 512k grid points — minutes per run on one A100.

**Tier B — 200 steps of the production jellium system, SIC-PZrun.**
Gates: `max_i |<psi_wp|psi_i>|` logged every step (the projection's own
overlaps); `cum_norm_removed < 1e-3`; bath orthonormality spot-check
(`max_{i≠j} |<ψ_i|ψ_j>| < 1e-8` at steps 0/100/200); `E_corrected` drift
measured, soft-gated at 200/1500 of the 0.1 eV budget and non-secular;
`E_PP(t)` decaying more slowly than in the uncorrected run (qualitative,
reported).

**Engine test** `test_wp_sic_engine.cpp` (runs in the tests SLURM stage):
(i) the kick applies exactly `exp(i dt v)` pointwise on a known Gaussian;
(ii) with a synthetic occupied manifold, the projected kick leaves every
`<ψ_j|ψ_wp>` at ≤ 1e-12 and removes the predicted O(dt²) norm;
(iii) `E_xc^unpol[n_wp]` from the header matches INQ's own `energy_xc` for a
1-electron system whose total density IS `n_wp` — pinning D1 run-consistency.

---

## 5. Production re-run  (intent unchanged; variant = Tier V winner)

Re-run only the **WP half** of `channeling_twin` with SIC-PZrun (pending Tier
V), identical in every physical parameter (`channeling_tube_rs3.hpp`, same GS,
same dt / N_STEPS / cadences). The classical half is unchanged and not re-run.
The new run lives in `scripts/channeling_sic/wp/` so the uncorrected run stays
untouched for the three-way comparison (classical / WP / WP+SIC).

### The hypothesis, as a chain of arguments

Links (1)–(7) unchanged from the original design (see git history for the full
chain): established theory (LDA self-interaction, no classical counterpart,
twin-verified at t = 0 to 8e-12 Ha) + measurements from the completed run
(E_PP releases 1.64 eV while var(p_perp)/2m gains 1.95 eV; 1.47× excess
expansion vs the exact Rayleigh baseline; impulse ratio 0.920 → 0.764 with
r = +0.98 against f_bore) + two physical arguments (a repulsive self-field
does expansion work; a more diffuse charge couples more weakly).

**Conditional:** *if* the self-interaction drives the excess expansion, and
*if* expansion weakens the coupling, *then* removing it changes the expansion
rate and, through that route, the measured stopping ratio.

⟦AMENDED per D8⟧ The refined expectation is a SPLIT, not a collapse: the
corrected run's residual deficit measures the genuine quantum-kinematic part
(Nazarov–Gross 2025 predicts one exists); the removed part measures the SIE
artefact. Three outcomes, all informative: (i) ratio → ~1: the deficit was
almost pure SIE; (ii) ratio rises but stays < 1: the split is quantified — the
headline result; (iii) ratio unchanged: the causal chain is refuted at its
first link.

### Observables recorded for the comparison

- transverse expansion vs the free-dispersion baseline; `f_bore(t)`,
  `<r_perp>(t)`, `sigma_r(t)`; `var(p_z)` and `var(p_perp)` separately;
  instantaneous + cumulative impulse ratio vs the classical twin; `S(T1)`,
  `S(T2)` over the user-chosen windows (9–25; 21–30; 5–20);
- NEW: `sic.csv` — `E_corrected`, `u_self`, `exc_self`, `max_overlap_pre`,
  `norm_removed`, `cum_norm_removed`.

---

## 6. Independent cross-check — DEFERRED

`sigma_WP` sweep on the UNCORRECTED binaries (`E_PP ~ 1/σ` separates SIE from
kinematics by a different route). Not launched in this campaign (GPU budget
goes to tests + Tier V + Tier B + production); remains the recommended hedge
if the production outcome is ambiguous.

---

## 7. Disk  ⟦RESOLVED per D7⟧

127 GB free measured 2026-08-02. The production WP run needs ~11 GB including
checkpoints and notebooks. The 5.1 GB of channeling_twin smoke checkpoints can
still be reclaimed but no longer block anything; left untouched.

---

## 8. Work breakdown  ⟦AMENDED: CSD3 chain⟧

| # | item | output |
|---|---|---|
| 1 | `inqkit/wavepacket/self_interaction_correction.hpp` | kick + projection + diagnostics + E_corrected |
| 2 | engine test `test_wp_sic_engine.cpp` (+ registration) | the three known-case checks of §4 |
| 3 | `scripts/wp_sic_vacuum/run.cpp` + SLURM stage | Tier V, 3 runs (env `SIC_MODE=none/h/pzrun`) |
| 4 | `scripts/channeling_sic/wp/run.cpp` (clone of channeling_twin/wp + SIC + sic.csv) | Tier B (200 steps) and production (1500) from ONE binary |
| 5 | `shared/bin/run-chan-sic.slurm` + `submit-channeling-sic.sh` | chain: tests → vacuum → tierB → prod |
| 6 | `hypotheses/channeling_sic/` | wp run notebook + refined-analysis notebook (structure identical to channeling_twin's), three-way comparison |
| 7 | handover + `docs/validation/test-catalogue.md` rows | resumability + validation record |

The GS is reused (`shared_gs/tube_rs3_R10_14_L60_dx0p5`); no SCF job needed.

---

## 9. Resolved questions (were "open questions for the user")

1. Disk: resolved (D7), nothing deleted.
2. Variant: Tier V decides; SIC-PZrun is the only candidate that can pass the
   closed-form gates (D1/D3).
3. σ_WP sweep hedge: deferred (§6), not launched.
4. Reviewer check of the projected-SIC derivation: DONE — the derivation is
   confirmed and identified with the one-sided Lagrange-multiplier form of
   variational TDSIC (Messud 2008); its energy-conservation consequence (D2)
   is now part of the design. The kick-project-renormalise concern resolves
   as O(dt²) per step within the Q subspace (`ψ_wp = Qψ_wp` holds identically
   under the scheme).
