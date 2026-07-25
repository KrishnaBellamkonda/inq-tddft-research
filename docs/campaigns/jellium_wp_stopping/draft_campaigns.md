# DRAFT — Jellium wavepacket-vs-classical stopping: three campaigns

> **Rough draft (2026-06-24).** To be refined / added to / subtracted from. Built
> from the `03_cap_stopping` baselines + this brainstorming session. Campaign 1
> will be formalised next via the `/campaigns` skill. Supporting analysis:
> `brainstorming-jellium-campaigns.ipynb`; restrictions:
> `notes_campaign1_sigma05_restrictions.md`. Literature anchor:
> `docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`.

## Shared premise (validated)
A projectile of identical charge & velocity but treated **quantum vs classically**
gives **different** electronic stopping — a real, publishable effect (Nazarov &
Gross 2025, Exact Factorization). The classical point-charge limit recovers
Lindhard; the quantum projectile's **width** matters. System: localised jellium
**slab** (r_s≈4, |z|<12.5, 50³ Bohr cell, two-sided sin² CAP). For every run the
WP and the classical Gaussian potential use the **same σ** (charge std =
σ_WP/√2), so the only intended difference is quantum-ness.

---

## Campaign 1 — `quantum-stopping-power` (small σ, point-like)  [formalise next]
**Aim.** Extract the **quantum stopping power** of a near-point-like wavepacket
(σ_WP = 0.5 Bohr) and compare WP vs matched classical vs point-Lindhard.

**Why σ small.** σ_WP=0.5 is the most point-like ⇒ most relevant to Lindhard. The
classical σ_WP=0.5 run is **validated ≈ point-Lindhard** (0.706 vs 0.716 eV/Bohr).

**Known restrictions (recorded, NOT blockers — see notes file).** 72× spreading;
stalled centroid ⇒ Δz ambiguous; no-wrap vs full-absorption incompatible in a
50-Bohr box; E_total ledger contaminated by zero-point (82 eV) + SIE (4.5 eV).

**Method (to settle in the formal campaign).**
- Primary estimator: **smoke-test the force/work-integral** `S=(1/Δz)∫⟨ψ|−∂_zV_ind|ψ⟩v dt`
  on existing `p5_wp` data first; if stable, use it (survives spreading). Fallback:
  E_total ledger with a **bigger box (≥80 Bohr) + run to full absorption**.
- Subtract zero-point KE; bound SIE with a per-σ **vacuum-WP control**.
- Matched classical at σ_pot=0.354 (`electron_gaussian_sigma0p35.upf`, exists).

**σ decision (2026-06-24): σ_WP = 0.5 LOCKED for the first test run.** Revisit
raising σ **only if the first run's results motivate it**. The σ_WP=1.0 question is
**undecidable from existing data** (the analytical finite-σ Lindhard over-suppresses;
no σ_WP=1.0 classical sweep exists) — so the deciding **σ_WP=1.0 classical S(v) vs
point-Lindhard** check is **deferred** until after the first σ=0.5 run, not a
pre-requisite. (See notebook.)

---

## Campaign 2 — large rigid σ (isolate quantum-vs-classical at fixed σ)
**Aim.** With a **rigid** (non-spreading) wavepacket, attribute any WP−classical
difference **purely to quantum-ness** (Pauli + interference), since spreading is
removed and σ is matched.

**Operating point (from the sanity-check agent).** σ_WP ≈ **4 Bohr at E ≥ 300 eV**:
full-run spread ≤12%, zero-point KE ~1.3 eV, SIE ~0.6 eV, transit ~8 a.u., fits
the 50-Bohr box (clears before wrap). σ_WP=3 has a ready UPF; σ_WP=4 needs one
generated (trivial, √2 convention). **Pick an energy for which a low-σ classical
reference run already exists** for comparability (to confirm in the Campaign-2
brainstorm).

**Caveats to bake in.** "Matched σ" does **not** fully isolate quantum-ness — the
WP carries **Pauli + SIE**, the classical ghost neither; quantify/bound the SIE.
The CAP was tuned for v≈2.7; a faster rigid packet absorbs less per length ⇒
**pilot-check CAP completeness** (residual norm→0), maybe raise η.

**Status.** Specifics deferred to a dedicated Campaign-2 brainstorm.

---

## Campaign 3 — muon classical vs wavepacket  [FUTURE — engine work required]
**Aim.** A heavier (muon-mass) **quantum** projectile spreads far more slowly
(τ_spread ∝ mass), giving a near-rigid quantum packet at *small* σ — so WP vs
classical isolates quantum-ness with minimal spreading.

**Initial observations (record for when this is revisited).**
- **Not supported in stock INQ.** A KS orbital's kinetic mass is hardwired to
  m_e (`inq/src/hamiltonian/ks_hamiltonian.hpp:202`, `operations::laplacian(...,−0.5,...)`,
  no per-state mass anywhere). A *quantum* muon needs an **`inq-study` engine
  fork**: a per-orbital ħ²/2m prefactor at the laplacian call sites (well-scoped,
  but a real engine modification — never edit `inq/`).
- A **classical** heavier-mass point charge is **trivial today** (ionic mass at
  `run.cpp` is tunable) — but that is not a quantum wavepacket.
- Tradeoff to remember (Nazarov-Gross): heavier mass also shrinks the de Broglie
  wavelength (more classical diffraction) — mass *maps* the quantum↔classical
  tradeoff, it doesn't escape it.
- **Decision deferred.** Revisit after Campaigns 1–2; the fork is the gate.
