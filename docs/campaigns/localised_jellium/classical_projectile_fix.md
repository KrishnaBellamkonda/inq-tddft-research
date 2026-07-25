---
# ROUGH DRAFT — authored interactively on 2026-06-27 to capture intent only.
# NOT autonomy-ready. To be built out carefully later via the `campaigns` skill
# (the five gated stages). Numbers, run matrix, and guard rails below are
# placeholders/intent, not locked decisions.
id: classical-projectile-fix
area: localised_jellium
title: Fixing the classical projectile (localised jellium stopping power)
status: draft
hypothesis: "A correctly-absorbed, cutoff-corrected classical Gaussian projectile through the localised jellium slab yields an equilibrated stopping power S_cl(v) — referenced to E_jellium(t=0)=E_total(t=0)-<T_WP>-SIE — that agrees with known results (Lindhard/RPA, our own bulk-jellium classical S(v), published DFT/TDDFT ion-in-jellium, and the slow-ion friction limit) within a tolerance to be fixed."
handover: docs/handovers/classical-projectile-fix.md
tasks:
  - { name: "Consume campaign #5 outputs (cutoff prescription, SIE, energy reference, dE-oscillation explanation)", done: false }
  - { name: "Gather benchmark S(v) targets (Lindhard/RPA + bulk-jellium classical + published DFT/TDDFT + slow-ion friction)", done: false }
  - { name: "Artefact-free projectile removal/absorption (no box re-entry, no abrupt-stop artefact)", done: false }
  - { name: "Equilibration: E_total(t) plateau before E_final is read", done: false }
  - { name: "Compute S_cl(v) by both methods (dE-slope and energy/L_z) and cross-check agreement", done: false }
  - { name: "Benchmark S_cl(v) against the four known-result targets; verdict + figure", done: false }
blocked_reason: ""
---

# Fixing the classical projectile (localised jellium stopping power)

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository.
</identity>

<rough_draft_banner>
This is a ROUGH DRAFT. Its job is to capture the user's intent and the simple
questions already resolved (2026-06-27), so that a later `/campaigns` session can
turn it into an autonomy-ready prompt. Do NOT execute it as-is. Sections marked
"(rough — to be locked via /campaigns)" are intent, not locked decisions.
</rough_draft_banner>

<dependency>
**HARD GATE: this campaign is sequenced AFTER campaign #5** ("Thorough exploration
of the ground state of the jellium slab system"). The user's decision (2026-06-27):
*"this needs understanding of the system carefully before fixing the classical
projectile."* Do not start the fixing/benchmark runs until #5 has delivered:

1. The **energy reference** for the localised jellium system, which is NOT the bare
   jellium ground state. It must be taken as
   `E_jellium(t=0) = E_total(t=0) - <T_WP> - SIE`
   (and, for the classical run, the analogous reference accounting for the
   classical projectile–jellium repulsion — #5 examines whether the current SIE
   estimate wrongly omits this classical repulsion term).
2. The **long-range cutoff prescription** for the classical projectile's radial
   Coulomb potential — owned entirely by #5 (literature search + electrostatics).
   #1 inherits whatever #5 concludes, to avoid loop-around / PBC self-interaction.
3. An **explanation of the ΔE-vs-time oscillation** seen in the classical runs
   (see <description>), so the gradient method can be applied with understanding.
</dependency>

<description>
**Why this campaign exists (user's words, lightly edited).** In the previous
"quantum stopping power" campaign the classical projectile runs had a few
problems, and the aim here is to come up with a classical version of the
localised jellium run such that **stopping power can be calculated effectively**,
then **benchmark it against known results** that must be sought out and compared.

There are two ways to extract the stopping power in the classical (Ehrenfest)
projectile scenario:
- **Gradient method** — take ΔE over time and find the slope of ΔE(t). (This is
  the method already written down in the `stopping-power-extraction` skill.)
- **Quantum-stopping-like method** — find the total energy gained by the localised
  jellium system and divide by the length of traversal L_z. **Care is required**
  with the reference: `E_jellium(t=0) = E_total(t=0) - <T_WP> - SIE`, not the bare
  jellium GS (see <dependency>).

Both definitions currently have problems:

1. **The CAP does not absorb the classical projectile.** Unlike the wavepacket
   (which is electron density and is absorbed cleanly), the classical projectile
   lives in `ions` and the density CAP cannot see it. In several runs the
   projectile **re-entered the simulation box**, changing the total-energy
   dynamics. This was worked around by **parking the classical electron at one of
   the edges** (`ions.remove(0)` at the CAP inner face), but the abrupt stop
   produced **artefacts that are not real**.
2. **The ΔE-vs-time graph looks like an oscillation**, and the cause is not yet
   understood. Diagnosing this is deferred to campaign #5; combining that
   understanding with screening / wake / long-range-Coulomb effects in the box
   (and the long-range cutoff of the projectile's radial potential, owned by #5)
   should explain it.
3. **The energy method needs the system to reach equilibrium before E_final is
   read.** In the classical runs seen so far this equilibration has not happened,
   chiefly because of problems (1) and (2).

**What "fixed" means here.** A classical localised-slab run that produces an
equilibrated, artefact-free stopping power S_cl(v), extractable by *both* methods
in mutual agreement, and consistent with the benchmark ladder below.

**Decision this informs.** Whether the classical projectile is a trustworthy
baseline for the quantum-vs-classical stopping-power comparison (the purpose of
the wider jellium stopping-power programme).
</description>

<benchmark_targets>
The classical localised-slab S(v) must be compared against ALL of the following
(user-selected, 2026-06-27), forming a layered validation ladder from
cheapest/in-repo to literature:

1. **Lindhard / RPA S(v)** — linear-response electron-gas stopping at matched
   r_s. Already in-repo via `inqview.analysis.lindhard_elf` and the existing
   point-Lindhard cross-reference. Cheapest, immediate.
2. **Our own bulk-jellium classical S(v)** — the project's periodic (bulk)
   jellium classical Gaussian-projectile S(v) at matched r_s/σ. Tests
   slab-vs-bulk consistency using identical machinery.
3. **Published DFT/TDDFT ion-in-jellium stopping** — literature
   (e.g. Echenique–Nagy–Nieminen–Puska nonlinear; Pruneda / Correa / Schleife
   TDDFT). Requires a `literature-review` pass to extract comparable numbers and
   record `docs/sources/` entries.
4. **Slow-ion friction limit** — the low-v analytic limit S ≈ Q·v (friction
   coefficient / transport cross-section), as a shape/asymptotics check.
</benchmark_targets>

<observables_set>
(rough — to be locked via /campaigns) Reuse the existing minimal + derived
observable set (ADR-0006, `minimum_observable_set.hpp`) at the current cadence.
Most relevant here:
- `E_total(t)` (for both extraction methods, with the corrected reference);
- Ehrenfest ion kinetic energy / ΔKE_ion and ion trajectory `z_ion(t)`
  (kinetic-channel cross-check);
- per-orbital energies; density VTI at the 300-frame target cadence;
- `N(t)` (bath-conservation guard — if the CAP drains the bath, E_total is
  garbage).
No new observable/kernel is anticipated; if one appears it is pre-gated
(`code-test` + `formula-validation` + catalogue row before any expensive run).
</observables_set>

<resolved_decisions>
Locked in this rough-draft session (2026-06-27):
- **Scope/sequencing:** standalone campaign, HARD-GATED behind campaign #5
  (see <dependency>). User: "I would block it until campaign 5 is not done."
- **Energy reference:** `E_jellium(t=0) = E_total(t=0) - <T_WP> - SIE` — not the
  bare jellium GS. The gradient method is the one written in the
  `stopping-power-extraction` skill.
- **Benchmark targets:** all four in <benchmark_targets>.
- **Long-range cutoff:** owned by campaign #5; inherited here.
- **Area/placement:** `docs/campaigns/localised_jellium/`.

To be locked later via /campaigns: geometry/N/r_s/box, the σ_WP/energy grid for
the classical sweep (likely matched to the existing WP runs), propagator + dt +
duration, the specific artefact-free removal/absorption scheme, and the numeric
agreement tolerance for "benchmark passes".

**Brainstorm progress (2026-06-30 /campaigns session — Stage 2, partial; PAUSED
to develop the wide-WP campaign first).** Locked-in-discussion (not yet written
into the run matrix; carry forward when this resumes):
1. **Extraction regime = BOTH methods, cross-checked.** Initial-drag −dKE/ds over
   the early v≥0.85·v₀ window AND traversal energy/L_z, required to *agree*.
   Insight: the two only agree when ΔE_deposited ≪ KE_proj (mild deceleration), so
   **"do the methods agree?" IS the diagnostic** for the linear-response regime —
   expected to hold at high v (490/340 eV), break at low v (23 eV, light electron
   self-stops). That breakdown is physics and sets the velocity floor of validity.
2. **Exit fix = NO artificial removal.** Drop the abrupt `ions.remove(0)` park
   entirely. Verified engine fact: INQ folds ion coords unconditionally in all 3
   dirs (`cell.hpp:219-227` `position_in_cell`; used by `ionic/interaction.hpp:41`
   + `ions.hpp:202`); there is **no "ion left the box" mechanism** — periodicity
   only switches the Poisson kernel (3D vs `poisson_solve_2d`). So re-entry is
   avoided by **truncating the energy read before the ion reaches the CAP face**,
   not by removal.
3. **Pilot matrix = A vs Z (two no-removal transit cells only; ramp-down dropped).**
   - **A** = periodicity-3 (PBC) — reuses #5's validated reference E_SIE ≈ 4.3 eV
     directly.
   - **Z** = periodicity-2 (open-z) — kills the projectile's z self-image (this IS
     #5's "thread-D cutoff"), but **requires re-deriving the energy reference under
     open-z** (#5 found −2.1 eV — unphysical — from a net-charge G=0 compensation
     term in the 2D kernel, and explicitly DEFERRED the fix). That re-derivation is
     a pre-gated code/research task for the Z cell.
   Keep the winner by cleanest E_total(t) + best two-method agreement.
4. **Concern (2) ΔE-oscillation = diagnosed BY the A-vs-Z comparison.** Leading
   hypothesis (backed by #5's static image-error finding + arXiv:2307.03213
   "periodic-image re-crossing into excited density"): a PBC z self-image artefact.
   If the oscillation collapses under open-z (Z) → confirmed artefact, open-z is
   the fix; if it survives → real physics (plasmon/wake ringing or finite-size
   standing mode), characterise vs ω_p / box modes. #5 did NOT dynamically diagnose
   it (GS-static work only), so it stays an open prediction tested here.
5. **CAP + box = MATCH THE WP (phase-5) RUN EXACTLY** (user, 2026-06-30): same box
   (50×50×90, slab |z|<12.5), same two-sided sin² CAP (η/faces ±35 as phase-5), so
   the slab GS and the overflow density are identical in WP and classical runs →
   apples-to-apples WP-vs-classical S(v). This SUPERSEDES the "elongate the box"
   idea: re-entry is handled by truncate-before-CAP-face within the matched box,
   not by a longer box. The neutral-slab GS is BC-independent (#5), so A-vs-Z is a
   projectile-electrostatics knob *within* the matched geometry.
   - **Open flag:** phase-5 WP is PBC; if Z (open-z) wins for the classical, decide
     whether the final WP-vs-classical overlay needs a matched-BC (open-z) WP point.
   - **Equilibration (concern 3):** with the CAP present E_total is NOT conserved
     (non-Hermitian drain) — must bookkeep CAP-drained energy and keep the N(t)≈const
     bath-drain guard; read E_final only after the post-slab plateau.

Still to lock when resumed: velocity/energy grid (match phase-5 {23,54,122,218,
340,490} eV), truncation criterion (z-threshold vs CAP face), the open-z reference
re-derivation (code-test + formula-validation + catalogue row), two-method
agreement tolerance, benchmark-ladder targets (Stage 3 research), notebook +
dispatch + handover mechanics.
</resolved_decisions>

<candidate_fix_approaches>
(rough — brainstorm to be sharpened via /campaigns) Candidate ways to "fix" the
projectile so the energy bookkeeping is clean:
- **Transit-only / elongated box** so the projectile never needs absorbing within
  the measurement window (avoids re-entry entirely; pairs with the long-range
  cutoff from #5).
- **Artefact-free ion removal** — replace the abrupt `ions.remove(0)` park with a
  smooth ramp-down / damped removal so no spurious energy is injected at the stop.
- **Ion-side absorbing scheme** — an absorber that acts on the ion channel (not
  just the density CAP), if feasible in the wrapper without editing `inq/`.
The choice depends on #5's cutoff prescription and on which keeps E_total(t) clean
through the equilibration window.
</candidate_fix_approaches>

<guard_rails>
(rough — to be locked via /campaigns)
- **Do not execute before campaign #5 is done** (hard gate, <dependency>).
- **Energy reference** must use `E_jellium(t=0)=E_total(t=0)-<T_WP>-SIE`, with the
  classical projectile–jellium repulsion handled per #5's conclusion.
- **Equilibration gate:** E_total(t) must plateau before E_final is read; if it
  does not, the result is an upper bound only ("not_converged").
- **N(t) ≈ const** — if the bath drains, both extraction methods are invalid.
- **Periodic-wrap truncation:** truncate before the projectile re-enters the box.
- **Two methods must agree** (gradient vs energy/L_z) within tolerance, else flag.
- Abort on NaN / complex energy / GPU occupied by another user; boundary 4σ/1σ
  and 300-frame VTI cadence per the always-on rules.
- Results are PROVISIONAL until they clear the benchmark ladder.
</guard_rails>

<tasks>
(rough — done-criteria to be sharpened via /campaigns)
1. **Consume campaign #5 outputs** — transcribe #5's cutoff prescription, SIE
   magnitude, the energy reference, and the ΔE-oscillation explanation into this
   campaign. *Done when:* #5 is `done` and its conclusions are recorded here.
2. **Gather benchmark targets** — assemble the four S(v) references
   (<benchmark_targets>); Lindhard/RPA from repo, published DFT/TDDFT via
   `literature-review` with `docs/sources/` entries, slow-ion friction limit,
   and our own bulk-jellium classical S(v). *Done when:* a benchmark table exists.
3. **Artefact-free projectile removal/absorption** — implement and verify a scheme
   (see <candidate_fix_approaches>) where the ion exits cleanly with no box
   re-entry and no abrupt-stop artefact in E_total(t). *Done when:* a run shows
   clean exit, N(t)≈const, no removal spike.
4. **Equilibration** — ensure E_total(t) plateaus before E_final. *Done when:* the
   equilibration gate passes.
5. **Compute S_cl(v) by both methods** — gradient ΔE(t) and energy
   [E_total(t_f)-E_jellium ref]/L_z; cross-check. *Done when:* both agree within
   tolerance.
6. **Benchmark** — compare S_cl(v) to the four targets; produce verdict + figure.
   *Done when:* comparison table/figure + pass/fail verdict recorded.
(Optional follow-on: σ-matched classical S(E) sweep once the single-point fix is
trusted.)
</tasks>

<rules>
- ALWAYS keep `inq/` immutable; engine-level experimentation goes in `inq-study/`,
  wrapper logic in `inqkit`.
- ALWAYS report S(v) numbers at 2 s.f. (3 s.f. only for genuine near-equalities).
- σ always means σ_WP; surface σ_pot=σ_WP/√2 only in a methods footnote.
- NEVER claim "fixed"/"benchmarked" on compile alone — only after the benchmark
  ladder is cleared.
</rules>

<preflight>
(rough draft — NOT yet autonomy-ready; this block is a reminder of what /campaigns
must still satisfy before flipping status away from draft)
- [ ] Dependency on campaign #5 resolved (cutoff + energy reference + ΔE
      explanation transcribed in).
- [ ] Geometry/N/r_s/box, σ/energy grid, propagator/dt/duration locked with values.
- [ ] Artefact-free removal scheme chosen and validated by a pilot.
- [ ] Equilibration + N(t) + wrap guards have numeric criteria.
- [ ] Benchmark tolerance defined; the four targets reproducible.
- [ ] Notebook output contract + handover pointer + Gmail/dispatch mechanics set.
</preflight>
