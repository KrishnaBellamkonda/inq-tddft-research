# 2026-05-17 — Jellium 2026-05-21 meeting campaign: design lock

**Type:** Topical (no single run; exempt from `run_summary.txt` verbatim rule per `.claude/rules/journal-entries.md` §2).
**Status:** Design phase complete. No compute started.
**Plan:** [`docs/plans/jellium-meeting-2026-05-21.md`](../../plans/jellium-meeting-2026-05-21.md)
**Handover:** [`docs/handovers/jellium-meeting-2026-05-21.md`](../../handovers/jellium-meeting-2026-05-21.md)
**Linked entries:** none yet — per-run entries will be created post-run
and linked back here.

## Why this entry exists

This is the design-lock journal entry preceding a multi-day jellium
campaign that targets the 2026-05-21 meeting with Emilio. It captures:

- The motivations for each run family.
- The new universal rules adopted this campaign (boundary, cadence, Δρ
  baselines, post-IFW shading, two-method stopping power).
- The open TODOs deferred past the meeting.
- The grilling dialogue that produced these decisions (key resolved
  ambiguities, so future readers don't re-relitigate them).

Per `.claude/skills/journal-writing`, observation text in run-based
entries is the user's voice and must not be invented. This is a
*design* entry — its content is the locked design, not run observations
— so it can be written before any compute.

## Motivations

### Why a campaign now

- Prior meeting (2026-05-14) presented a stopping-power figure with a
  classical-WP comparison at a limited energy spread. Emilio raised
  questions about (i) the σ-dependence (professor's σ→∞
  "indistinguishable from jellium" intuition), (ii) the unaccounted-
  overlap puzzle (missing electrons in the GS-basis projection), and
  (iii) whether the classical-WP stopping-power agreement claimed in
  Knudsen et al., arXiv 2605.12854, holds at our setup.
- 4-day window to the 2026-05-21 meeting is tight but feasible if
  pairs run concurrently on the two A30s.

### Why this energy menu

- **20-25 eV**: DNA radiolysis scale; low-v anchor.
- **100 eV**: existing canonical pair on disk; comparison partner for
  σ-sweep, extra-states test, free-space WP, high-density.
- **300, 600, 1500 eV**: existing pairs; retroactive Knudsen analysis
  via momentum_distribution.csv.
- **700-1100 eV** (Knudsen sweep): brackets the regime where the paper
  claims classical-WP convergence (v ≈ 5-7 Bohr/atu).

### Why a separate free-space WP

To distinguish *what jellium does to the WP* from *what the WP does
anyway*. Free-particle Gaussian dynamics has an exact analytical
solution (σ_r(t)² = σ_r(0)² + (σ_p(0)·t)²). Three independent traces
of σ_r²(t) — INQ non-interacting, Python split-step, analytical —
also validate the WP injector + propagator under controlled conditions.

## Key resolved ambiguities (from the design grilling)

These are decisions locked in conversation that aren't otherwise
discoverable from the codebase or plan:

1. **σ in the "5σ from boundary" rule** ⇒ σ ≡ `WP_SIGMA_BOHR` (real-
   space Gaussian width); not FWHM, not σ_z. Standard launch is now
   `−L/2 + 4σ`, stop is `+L/2 − σ` ⇒ traversal `L − 5σ`.

2. **Existing E=100 eV pair** (legacy convention) is *not* re-run as
   the canonical reference — it's preserved as the meeting-deck
   reference, and a *new* `_v2` pair is run under the new 4σ/1σ rule
   as the comparison partner for the free-space run and σ-sweep.

3. **σ=8 large-σ outlier** uses a relaxed `−L/2 + 3σ` launch (and
   `+L/2 − σ` stop) because the standard 4σ rule is geometrically
   infeasible in L=50 above σ=5. Documented in plan with a TODO to
   test the σ→∞ claim in a bigger box post-meeting.

4. **Δρ baselines** are ADD, not REPLACE: keep `density_rt_delta`
   (GS-subtracted) for comparability with prior figures; add
   `density_rt_delta_t0` (t=0-subtracted) to remove the static WP
   density feature from the colormap.

5. **Knudsen stopping-power method**: implemented as a new C++
   `WPMomentumStats` observable for new runs (decomposed
   `<p_d>`/`σ_p_d`), plus a degraded Python retroactive path
   (`<|k|²>(t)` only) from the existing `momentum_distribution.csv`.
   Both methods plotted alongside the user's existing INQ orbital
   eigenvalue method on the final stopping-power curve.

6. **Extra-states test**: current default is 20, not 8 (as initially
   stated). The useful experiment is *going higher* (test 40, 80) to
   see whether the unaccounted-overlap drops monotonically toward
   zero, which would confirm the basis-completeness hypothesis. Two
   new GSes required.

7. **Email infrastructure**: Gmail MCP route (Option A from the
   grilling) — needs one-time interactive OAuth from the user.
   Self-to-self at `chiddukanna@gmail.com`. One email per pair,
   threaded per family via In-Reply-To, plus a final rollup per
   family.

8. **GPU time budget**: no hardware wall-time limit on this
   workstation; the "9 h GPU-0 budget" mentioned in
   `electron_proj_E25_L50_cubic.hpp` was self-imposed daily-turnaround
   prudence, lifted explicitly for this campaign so long classical
   companions can run.

9. **WP orbital VTI** is *not* saved by current runs (decision in
   `run_wp_n162_L50_E100/run.cpp` due to slow `density::orbital`
   extraction). For Knudsen-style estimates, `WPMomentumStats` works
   in k-space directly and does not require the orbital VTI.

10. **VTI cadence rule**: target ≈ 300 frames per run via
    `WRITE_EVERY = max(1, round(N_STEPS / 300))`. Comparable gif
    durations and disk usage across the campaign.

## Open TODOs (not in critical path)

- **σ→∞ bigger-box test** (task #31): L=50 caps σ at 8 Bohr; a
  proper σ→∞ "indistinguishable from jellium" test needs L=100+
  (expensive GS, deferred post-meeting).
- **Read KL divergence + 2-Wasserstein theory** (task #29): user
  learning goal.
- **Fix `MomentumDistribution` MPI-reduce** (task #30): currently
  per-rank-writes a shared CSV without `all_reduce`. Fine while
  single-rank, latent bug at scale. Fix in Infra-4 when nearby.
- **25 eV classical companion** (task #25): explicitly post-Knudsen
  per user direction.
- **Workstation watchdog** (task #32): optional thermal/Xid monitor.

## Files associated with this campaign (created or to be created)

| File | Status | Purpose |
|---|---|---|
| `docs/plans/jellium-meeting-2026-05-21.md` | created today | full design lock |
| `docs/handovers/jellium-meeting-2026-05-21.md` | created today | rolling status |
| this file | created today | design journal entry |
| `docs/reports/2026-05-21-meeting-emilio/figures/stopping_power_v2.png` | pending Final task | meeting deliverable |
| `docs/sources/knudsen-2025-electron-graphene.md` | pending | source note for arXiv 2605.12854 |
| `shared/configs/boundary_rule.hpp` | pending Infra-10 | universal launch/stop helpers |
| 6+ new Cfg headers | pending per-run tasks | one per new energy |
| `inq-stack/include/inqkit/observables/wp_momentum_stats.hpp` | pending Infra-4 | Knudsen-method observable |
| `inq-stack/include/inqkit/observables/wp_real_space_stats.hpp` | pending Infra-5 | sibling |
| `Tutorial/wp-momentum-stats-test/` | pending Infra-4 | known-case smoke test |

## What was *not* decided (will need future input)

- **Exact threading mechanism for Gmail MCP**: the MCP's support for
  `In-Reply-To`/`References` headers vs. subject-prefix-only threading
  will be verified once Infra-0 OAuth lands and Infra-8 is being built.
- **dt for Knudsen sweep**: gated on the dt-convergence subtest result
  (10-min test). Either dt=0.01 or dt=0.005.
- **Whether `WPRealSpaceStats` cadence can be every step in
  production**: gated on profiling first Knudsen run.

## Attribution

This design entry summarises decisions reached in conversation between
the user and the assistant on 2026-05-17 (grilling session via
`grill-with-docs` skill). Sources cited inline:

- Knudsen et al., *Ultrafast electron dynamics of electron-irradiated
  graphene*, arXiv 2605.12854 (motivates §700-1100 eV sweep + the
  `<p²>/2` stopping-power method).
- 2026-05-14 meeting figures and scripts in
  `docs/reports/14-05-2026-meeting-emilio/figures/` (baseline for the
  final stopping-power rollup).
- `docs/sources/free-electron-gas-magic-numbers.md` (N=162 still
  closed-shell at L=30).
- INQ `options::theory{}.non_interacting()` (`inq/src/options/theory.hpp:42`)
  — enables the free-space WP run on the same INQ pipeline.
