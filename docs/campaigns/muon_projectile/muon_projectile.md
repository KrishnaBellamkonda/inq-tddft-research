---
# ROUGH DRAFT — consolidated 2026-06-27 from the existing muon spec
# (draft_campaigns.md "Campaign 3" + localised-jellium.md handover) so a later
# `/campaigns` session can build the executable plan. NOT autonomy-ready.
# Fork ownership, mass grid, and run specifics are deferred to that pass.
id: muon-projectile
area: muon_projectile
title: "Muon projectile — heavier-mass quantum WP vs classical (slow-spreading route)"
status: draft
hypothesis: "A muon-mass quantum wavepacket spreads ~m times more slowly and carries ~m times less zero-point KE than an electron WP, so WP-vs-classical stopping on the jellium slab can be compared with minimal dispersion and a clean energy ledger — but because heavier mass also shrinks the de Broglie wavelength (more classical), the campaign actually maps the quantum→classical crossover vs projectile mass rather than escaping it."
handover: docs/handovers/muon-projectile.md
tasks:
  - { name: "Confirm fork ownership + sequencing (after Campaigns 1–2); reframe hypothesis around the de Broglie tradeoff if needed", done: false }
  - { name: "Classical muon scoping runs (tunable ionic mass; map S vs mass / crossover, no engine work)", done: false }
  - { name: "inq-study per-orbital-mass fork (per-orbital ħ²/2m at the laplacian call sites; ONLY the WP state) — code-test + formula-validation", done: false }
  - { name: "Quantum muon WP vs classical muon — isolate quantum-ness; relate to Campaign 2 wide-σ results", done: false }
  - { name: "Analysis: quantum→classical crossover vs mass; WP−classical quantum component with SIE bounded", done: false }
blocked_reason: ""
---

# Muon projectile — heavier-mass quantum WP vs classical (slow-spreading route)

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository. σ always means the wavepacket width σ_WP.
</identity>

<rough_draft_banner>
This is a ROUGH DRAFT consolidating the muon spec that was already written
(scattered across `draft_campaigns.md` and the `localised-jellium.md` handover).
Its job is to gather it in one place, record the 2026-06-27 verification + spec
check, and tee it up for the careful `/campaigns` pass. Do NOT execute it as-is.
"(deferred to the /campaigns pass)" marks intent, not locked decisions.
</rough_draft_banner>

<prior_art>
This campaign already existed in sketch form — this file consolidates it:
- `docs/campaigns/jellium_wp_stopping/draft_campaigns.md` → **"Campaign 3 — muon
  classical vs wavepacket [FUTURE — engine work required]"**.
- `docs/handovers/localised-jellium.md`:686–688 (muon quantum WP infeasible in
  stock INQ; classical-muon trivial; "Muon = Campaign 3 (FUTURE)") and :305
  ("clean routes: large-σ / muon (zp→0)").
- `docs/campaigns/jellium_wp_stopping/brainstorming-jellium-campaigns.ipynb`.
- Recast in the supervisor presentation as an energy×σ broadening outlook
  (`docs/handovers/thesis-supervisor-meeting.md`:209).
</prior_art>

<description>
**Aim (user's words + existing spec).** A heavier (muon-mass) **quantum**
projectile spreads far more slowly (τ_spread ∝ mass), giving a near-rigid quantum
packet at *small* σ — so WP vs classical isolates quantum-ness with minimal
spreading. This is a **second, independent route** to the same goal as Campaign 2
(wide low-spread WP): Campaign 2 reaches near-rigidity by raising σ; the muon
reaches it by raising mass.

**Two extra advantages of the muon route (record these):**
1. **Clean energy ledger.** The zero-point KE that contaminates the electron WP's
   E_total stopping ledger (≈ 3/4σ² ≈ 82 eV at σ=0.5) scales as **1/m**, so a muon
   (m_μ ≈ 206.77 m_e) carries ~207× less (~0.4 eV) — the ledger method becomes far
   cleaner (handover :305).
2. **Small-σ rigidity** keeps the interaction range close to the classical point
   charge without needing the large box that Campaign 2's wide σ forces.

**Decision this informs.** Whether projectile *mass* is a usable knob to access the
quantum-vs-classical stopping difference cleanly — and where, along the mass axis,
the quantum packet effectively becomes classical.
</description>

<spec_check>
Resolved / verified 2026-06-27:

1. **Engine gate — VERIFIED.** A *quantum* muon is not supported in stock INQ. The
   KS-orbital kinetic prefactor is hardwired to m_e at
   `inq/src/hamiltonian/ks_hamiltonian.hpp:202`
   (`operations::laplacian(phi_fs, -0.5, …)`), and likewise at `:235`
   (`laplacian_add(…, -0.5, …)`) and `:245`
   (`laplacian_expectation_value(…, -0.5, …)`). The `-0.5` = −ℏ²/2m with m=1. It is
   applied to the **whole orbital_set at once**, so a global change would make the
   *electrons* muon-mass too. The fix must be a **per-orbital mass** so only the WP
   state gets −1/(2m_μ): a well-scoped but genuine `inq-study` modification (NEVER
   edit `inq/`).
2. **Classical muon — trivial.** The ionic mass at `run.cpp` is tunable, so a
   classical heavy-mass point charge needs no engine work — but it is not a quantum
   wavepacket.
3. **SPEC CAVEAT to confirm/reframe (the central one).** Heavier mass cuts both
   ways: it slows spreading (near-rigid packet) **but also shrinks the de Broglie
   wavelength → less quantum diffraction/interference → a *more classical*
   projectile** (Nazarov–Gross). Mass therefore *maps* the quantum↔classical
   tradeoff; it does not escape it. The hypothesis is framed accordingly (map the
   crossover vs mass). Confirm or reframe in the careful pass.
</spec_check>

<observables_set>
(deferred to the /campaigns pass) Reuse the full observable suite + current cadence
(ADR-0006, `minimum_observable_set.hpp`). Especially: WP momentum distribution
n(k,t) (coherent-peak centroid → S), WP real-space + momentum spread stats
(`wp_real_space_stats`, `wp_momentum_stats`) to confirm slow spreading, E_total
ledger (now clean), classical ion track (z,v,F).
</observables_set>

<resolved_decisions>
Locked in this rough-draft session (2026-06-27):
- **Route:** projectile *mass* as the rigidity knob (sibling to Campaign 2's σ knob).
- **Engine gate verified** (per-orbital ħ²/2m fork in inq-study; line-refs above).
- **System:** localised jellium slab (same as Campaigns 1/2), σ-convention σ_WP.
- **Hypothesis framing:** crossover-mapping, accounting for the de Broglie tradeoff.
- **Placement:** `docs/campaigns/muon_projectile/` (own folder; sibling of the
  jellium_wp_stopping programme it grew out of as "Campaign 3").

Deferred to the /campaigns pass:
- **Fork ownership** — candidates: (a) staged (classical scoping first, then fork +
  quantum muon); (b) fork in-scope/central; (c) fork as external prerequisite
  (campaign blocked until it lands).
- **Mass grid** — single muon mass (206.77 m_e) vs a mass sweep to map the
  crossover (the de Broglie tradeoff argues for a sweep).
- σ/energy/box/dt/duration; sequencing relative to Campaigns 1–2.
</resolved_decisions>

<guard_rails>
(deferred to the /campaigns pass)
- **Fork is the gate:** no quantum-muon run before the per-orbital-mass fork passes
  code-test + formula-validation (free-muon KE = ℏ²k²/2m_μ; spreading rate τ∝m;
  electron orbitals' mass unchanged).
- Sequenced after Campaigns 1–2 (the σ=0.5 + wide-σ campaigns).
- SIE bounded via vacuum-WP control before any "quantum component" is reported.
- Boundary 4σ/1σ + 300-frame VTI; abort on NaN / complex energy / GPU occupied.
- PROVISIONAL until the fork is validated and the crossover framing is confirmed.
</guard_rails>

<tasks>
(rough — done-criteria sharpened via /campaigns)
1. **Confirm fork ownership + sequencing**; reframe the hypothesis around the de
   Broglie tradeoff if the user wants. *Done when:* ownership + framing locked.
2. **Classical muon scoping** — tunable-ionic-mass runs mapping S vs mass /
   crossover (no engine work). *Done when:* a classical S(mass) trend exists.
3. **inq-study per-orbital-mass fork** — per-orbital ħ²/2m at `:202/:235/:245`,
   applied only to the WP state; pre-gate. *Done when:* code-test +
   formula-validation pass; catalogue row added.
4. **Quantum muon WP vs classical muon** — isolate quantum-ness; relate to Campaign
   2. *Done when:* matched WP/classical pair run + compared.
5. **Analysis** — quantum→classical crossover vs mass; WP−classical quantum
   component with SIE bounded; executed notebook. *Done when:* notebook + verdict.
</tasks>

<rules>
- NEVER edit `inq/`; the per-orbital-mass change goes in `inq-study` only.
- The mass change is **per-orbital** — the electron orbitals must keep m_e.
- σ always means σ_WP; report S at 2 s.f. (3 s.f. only for near-equalities).
- NEVER report a WP−classical "quantum component" without the SIE bounded.
- Ground the mass/de-Broglie tradeoff per Nazarov–Gross (literature-review).
</rules>

<preflight>
(rough draft — NOT yet autonomy-ready; reminder of what /campaigns must satisfy)
- [ ] Fork ownership + sequencing + hypothesis framing locked.
- [ ] inq-study per-orbital-mass fork implemented + pre-gated (free-muon KE +
      spreading + electrons-unchanged tests); catalogue row.
- [ ] Mass grid, σ/energy/box/dt/duration locked with values.
- [ ] SIE-bounding control defined; coherent-peak S estimator fixed.
- [ ] Notebook contract + handover pointer; engine claims carry source line-refs.
</preflight>
