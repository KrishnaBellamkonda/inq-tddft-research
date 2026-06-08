# Plan: Jellium WP-jellium scattering campaign for 2026-05-21 Emilio meeting

**Status:** Design phase complete 2026-05-17. Locked.
**Owner:** chiddukanna (skcb2)
**Deadline:** 2026-05-21 (Thursday, +4 days).
**Related:** [Design journal entry](../journals/researchproject/2026-05-17_jellium_meeting_design.md), [Rolling handover](../handovers/jellium-meeting-2026-05-21.md), prior meeting `docs/reports/14-05-2026-meeting-emilio/`.

## Goal

Deliver a multi-energy stopping-power curve for an electron projectile in
jellium covering **wave-packet (WP)** and **classical-point-electron**
projectiles, compared via two independent stopping-power estimators
(INQ orbital eigenvalue trajectory vs Knudsen-style `<p²>/2` from the
momentum-space WP), and contextualised against the 14-05 meeting figures
and Knudsen et al., *Ultrafast electron dynamics of electron-irradiated
graphene* (arXiv 2605.12854) which claims classical-WP stopping-power
agreement near v ≈ 5-7 Bohr/atu (E ≈ 800 eV).

Secondary objectives:

- Probe σ-dependence at fixed E=100 eV (6 σ values) to address the
  professor's σ→∞ "indistinguishable from jellium" intuition.
- Test the basis-completeness hypothesis for the unaccounted-overlap
  ("missing electrons") issue at E=100 eV via `extra_states ∈ {20, 40, 80}`.
- Add anchors at E=20 eV (DNA energy scale, WP+classical pair) and
  high-density (L=30³ at fixed N=162, r_s ≈ 3.41) point.
- Validate the WP injector + propagator against free-particle dynamics
  via an INQ-non-interacting free-space run + Python Schrödinger toy +
  analytical Gaussian-spread reference.

## Universal rules (new this campaign)

### 1. Boundary rule (`shared/configs/boundary_rule.hpp`)

| Rule | Standard | Large-σ relaxed |
|---|---|---|
| Launch z (centroid) | `−L/2 + 4σ` | `−L/2 + 3σ` |
| Stop z (centroid) | `+L/2 − σ` | `+L/2 − σ` |
| Traversal | `L − 5σ` | `L − 4σ` |
| t_IFW end (centroid at) | `+L/2 − 3σ` | `+L/2 − 3σ` (same — Gaussian-tail-at-far-face criterion) |
| Applies to | every new run | σ=8 outlier only |

**Sanity check at σ=5, L=50, v=2.711, dt=0.02:** launch z=−5, stop z=+20,
traversal=25 Bohr, **N_STEPS=462** at dt=0.02 (boundary_rule.hpp ceil
convention), IFW end at centroid +10, **t_IFW=5.53 a.u.**, t_total=9.22
a.u. ⇒ post-IFW shaded region `[5.53, 9.22]` a.u. (40% of trajectory).

**Sanity check at σ=8 relaxed, L=50, v=2.711:** launch=−1, stop=+17,
traversal=18, N_STEPS_relaxed=332, IFW end at centroid +1 (only 2 Bohr
past launch — ~10% of trajectory). Most of σ=8 is post-IFW; informative
for σ→∞ TODO #31 anyway.

Codified as constexpr `launch_z(σ, L)`, `stop_z(σ, L)`, `n_steps_for(L, σ, v, dt)`.

### 2. Cadence rule

`write_every_for(n_steps, target=300)` — target ≈ 300 VTI frames per run
across the campaign for comparable gif lengths and disk usage.

### 3. Δρ baselines (ADD don't REPLACE)

- **Keep** `density_rt_delta` (= ρ(t) − ρ_GS) for legacy comparability.
- **Add** `density_rt_delta_t0` (= ρ(t) − ρ(t=0_post_WP_injection)) to
  remove the static WP-density feature from the colour map. Dual path:
  Python phase for retroactive analysis on existing runs; C++ writer in
  inqkit for new runs.
- Same dual treatment for `delta_density_z_profile_t0` (1D
  along-projectile-axis variant).

### 4. Post-IFW shading

Universal helper in `inqview/postprocess/_common.py` shades the
post-interference-free-window region on every time-domain plot. Derived
per Cfg from launch_z, stop_z, v.

### 5. Email pipeline (MUST USE inqview.email)

**All campaign emails — including `[PAUSE-NEEDED]` alerts, completed-run
reports, and the final rollup — MUST be sent via
`inqview.email.send_run_email(...)` (the Python library autosend path).
The Gmail MCP `create_draft` tool MUST NOT be used.**

Reason: `inqview.email` autosends via Gmail SMTP + app-password and
returns a `Message-ID` that subsequent emails in the same family
(e.g. `[jellium-25eV]`) thread under via `in_reply_to`. The MCP
draft path requires manual user click-through, breaks the auto-loop
contract, and produces threads that the rollup cannot stitch back
together programmatically. Per the 2026-05-17 handover decision,
production = autosend, not draft-and-send-manually.

Conventions when calling `send_run_email`:

- `subject` prefix `[PAUSE-NEEDED]` for any blocker that needs user
  judgement (validation fail, >2 h compute about to launch, design
  ambiguity surfaced mid-run).
- `subject` prefix `[jellium-<family>]` for completed-run reports
  (e.g. `[jellium-25eV]`, `[jellium-knudsen-sweep]`, `[jellium-rollup]`).
- Attach the relevant figures (`knudsen_ke_vs_t.png`,
  `stopping_power_vs_v.png`, …) via the `attachments` list.
- Use `in_reply_to=<prev_msg_id>` to keep the family threaded.

### 6. Stopping-power estimators (both used in the final rollup)

- **Method A (existing)**: WP orbital eigenvalue ε_WP(t) — the user's
  established INQ-native method.
- **Method B (new, Knudsen)**: `<|k|²>(t)/2` from the new
  `WPMomentumStats` observable (or `momentum_distribution.csv`
  histogram in degraded retroactive mode).
- **Classical**: `½m_e v_ion(t)²` from `ehrenfest()` trajectory.

All three plotted per energy on the final stopping-power vs v curve.

### 7. Per-run `analyse.py` (mandatory)

**Every new run directory MUST ship a self-contained `analyse.py` that
is invoked once after the run completes.** This file imports the
campaign's inqview pipeline, runs the archetype's phase list, runs
custom analyses via `shared/python/analyse_extras.py`, and writes
`results/analysis/REPORT.md`. The email step happens AFTER analyse.py
finishes — emails reference `analysis/observables/*.png` artefacts that
only exist if the pipeline ran.

Three archetypes are recognised; each new run dir clones the matching
one:

- **WP archetype** (jellium WP runs): 17 phases incl. `wp_trajectory`,
  `gamma_transitions`, `eigenvalues_gs`, `state_energies`,
  `state_energy_spectra`, `occupations`, `momentum`, `knudsen_ke`,
  `kl_divergence`, `energy_balance`, `density`, `overlap`, `orbitals`.
- **Classical archetype** (jellium classical runs): 18 phases — drops
  `wp_trajectory`/`gamma_transitions`, adds `bath_energy`, `stopping`,
  `gs_projected_occupations`. Includes `knudsen_ke` + `kl_divergence` +
  `energy_balance` which silently skip on classical inputs that lack
  the WP data.
- **Free archetype** (non-interacting box + WP): 8 phases —
  `summary, layout, observables, occupations, momentum, knudsen_ke,
  kl_divergence, density`. No GS, no eigenvalues, no overlap.

Compare-style emails (free vs jellium WP at the same energy) live in a
separate `_compare_<family>_<energy>/compare.py` next to the run dirs,
NOT inside any one run dir. compare.py runs AFTER both halves'
analyse.py have completed.

### 8. Observables-reference §13 implementation (pinned 2026-05-17)

The following observables-reference §13 rules MUST be honoured by every
new run + every pipeline phase from 2026-05-17 onwards:

- **§13.1.1** ScalarFormatter(useOffset=False) — already enforced by
  `inqview/plots.py::plot_spectrum`; carried into new phases.
- **§13.1.3** HOMO dashed line — enforced by `state_energies.py`.
- **§13.1.4** Fixed colour scale across animation frames — enforced
  by `density.py::_global_vmin_vmax`.
- **§13.1.5** Smooth FFT spectra — `FourierTransform(zero_pad=4)` default.
- **§13.3** Energy-balance ledger (ΔE_WP, ΔE_bath, ΔE_total,
  Unaccounted) — implemented as the new `energy_balance` phase.
- **§13.5** Per-component energy time series — implemented in
  `observables.py::_plot_per_component_energy` (separate PNG per
  energy_total/energy_kinetic/energy_hartree/energy_xc).
- **§13.6** FFT transient-region exclusion — `FourierTransform` accepts
  `t_start_au`; spectra produced from `observables.csv` declare the
  cutoff used in their CSV header.

The legacy E={50, 100, 300, 600} eV WP runs lack the new
`WPMomentumStats` / `WPRealSpaceStats` CSVs — they are kept on the
post-meeting backlog (task #25) for re-running with the new C++
observables.

### 9. Plot shading convention (pinned 2026-05-17)

Time-domain plots HIGHLIGHT the interference-free window (IFW) with a
soft-yellow `axvspan(0, t_IFW)` (helper `_common.ifw_highlight`). The
older post-IFW grey shading helper (`_common.post_ifw_shade`) is kept
for backwards compatibility but NEW plots should prefer the highlight
form — the IFW is the region of interest, not the noise.

### 10. Classical kinetic-energy caption (pinned 2026-05-17)

Every plot showing the classical projectile's kinetic energy MUST
carry a caption explaining that the energy DROP under Ehrenfest
dynamics IS the electronic stopping signal (force from bath density
gradient transferring energy to electronic excitations), not a
numerical artefact. Hand-rolled comparison scripts and inqview phases
both observe this rule.

### 11. New C++ observables (`inq-stack/include/inqkit/observables/`)

- `wp_momentum_stats.hpp` — per-step `<p_d>`, `<p_d²>`, `σ_p_d²`,
  `E_kin`. Host-after-reduction pattern (mirror INQ `dipole.hpp`).
  Heisenberg known-case test in `Tutorial/wp-momentum-stats-test/`.
- `wp_real_space_stats.hpp` — sibling: `<x_d>`, `<x_d²>`, `σ_r_d²`.
- Cadence: every step initially; profile after first Knudsen run,
  fall back to every 2 or 5 if > 5 % overhead.

## Run inventory

| # | Run | dx | σ (Bohr) | N_STEPS | dt | WE | GS | Compute | Email family |
|---|---|---|---|---|---|---|---|---|---|
| Infra-0..10 | Universal infra | — | — | — | — | — | — | ~16 h dev | — |
| 1 | 25 eV WP | 0.40 | 5 | 922 | 0.02 | 3 | dx0p40 | ~75 min | `[jellium-25eV]` |
| 2 | Free-space WP @ 100 eV (INQ non-int) | 0.40 | 5 | 461 | 0.02 | 2 | fresh trivial | ~10 min | `[jellium-free-compare]` |
| 2b | Python Schrödinger toy @ 100 eV | 0.40 | 5 | 461 | 0.02 | 2 | n/a | ~5 min CPU | (in family) |
| 3 | E=100 v2 jellium (WP + Classical) | 0.40 | 5 | 461 | 0.02 | 2 | dx0p40 | ~30 min | `[jellium-free-compare]` |
| 4 | Extra-states test (WP, extra ∈ {20,40,80}) | 0.40 | 5 | 461 | 0.02 | 2 | new x40, x80 | 60+90 min + 2×30 min GS | `[jellium-extra-states]` |
| 5 | 20 eV WP+Classical | 0.40 | 5 | 1030 | 0.02 | 3 | dx0p40 | ~6 h classical | `[jellium-20eV]` |
| 6 | High-density 30³ WP+Classical | 0.40 | 0.5 | 507 | 0.02 | 2 | new gs_L30_N162 | ~60 min + 30 min GS | `[jellium-highdensity]` |
| 7 | σ-sweep @ E=100 (WP, σ ∈ {0.25,0.5,1,3,5,8}) | 0.30/0.40 | varies | 332-899 | 0.02 | 1-3 | dx0p30 + dx0p40 | ~4-5 h | `[jellium-sigma-sweep]` |
| Subtest | dt-convergence @ E=1100 (dt=0.005 vs 0.01) | 0.30 | 5 | 100 each | varies | 1 | dx0p30 | ~10 min | (subtest) |
| 8 | Knudsen sweep (WP+Classical pairs, E ∈ {700-1100}) | 0.30 | 5 | 278-348 | 0.01 or 0.005 | 1-2 | dx0p30 | ~3-5 h | `[jellium-knudsen-sweep]` |
| 9 | 25 eV classical companion (post-Knudsen) | 0.40 | n/a | 922 | 0.02 | 3 | dx0p40 | ~6 h | `[jellium-25eV-classical]` |
| Final | Retroactive Knudsen on E={50,300,600,1500} | (existing) | 5 | (existing) | — | — | — | ~30 min analysis | (in rollup) |
| Final | Stopping-power rollup + meeting figure | — | — | — | — | — | — | ~1 h | `[jellium-rollup]` |

**Total compute estimate**: ~25-30 h on 2 GPUs running pairs concurrently.
**Total disk estimate**: ~30 GB (12 new runs × 300 frames × 6 MB + dx=0.30 runs heavier).

## Per-family motivations

### Run-1 — 25 eV WP

Existing planned WP run not yet executed (`run_wp_n162_L50_E25/` has
only `run.cpp` + `analyse.py`). Energy κ = 2/v = 1.47 (Bohr classical
regime), v/v_F = 4.02. Anchor for low-v end of the stopping-power curve.
Cfg already exists at `shared/configs/electron_proj_E25_L50_cubic.hpp`;
update to new 4σ/1σ rule (launch z=−5, N_STEPS=922) before launch.

### Run-2 + 2b + 3 — Free vs jellium @ E=100 eV

Validate the WP injector against free-particle physics; isolate the
*jellium effect* on the WP. Three independent traces of σ_r²(t):
INQ-non-interacting (Run-2), Python split-step (Run-2b), analytical
Gaussian (closed form). Compare against jellium (Run-3). Comparison
plots: side-by-side density gif, diff gif, σ_r²(t) growth, σ_p²(t)
per direction, KL divergence over time, ⟨z⟩(t).

### Run-4 — Extra-states test @ E=100 eV

Tests the basis-completeness hypothesis for the
unaccounted-overlap/missing-electrons puzzle. Current default
EXTRA_STATES=20; test {20, 40, 80}. Going *higher* fixes the issue if
my-projection-basis-is-incomplete hypothesis is right. Going lower
(<20) is uninteresting (known to make it worse). New GSes required at
40 and 80 extra states.

### Run-5 — 20 eV DNA-scale pair

Water in DNA radiolysis sits near this energy. WP + classical pair
maximally informative; classical alone is ~6 h but adds the
matched-energy comparison point and a low-v anchor for the stopping
curve below 25 eV.

### Run-6 — High-density 30³ pair

L=30 cubic at fixed N=162 raises density: r_s drops from 5.69 → 3.41
Bohr (Na → between Li and Al). σ=0.5 Bohr WP, E=100 eV. Direction-finder
for the density-dependence of stopping power. N=162 remains magic at
any L (filled shell at `|G|² ≤ 6`). Classical companion adds context
(denser bath ⇒ stronger classical retardation).

### Run-7 — σ-sweep at E=100 eV

Six WP-only runs at σ ∈ {0.25, 0.5, 1, 3, 5, 8} Bohr. Tests how the
WP's spatial extent governs scattering vs elastic propagation. σ=0.25
needs dx=0.30 for Nyquist (σ_p=2.0 Bohr⁻¹); σ=0.5/1/3/5 use dx=0.40;
σ=8 uses the relaxed 3σ/1σ boundary rule (geometric necessity). The
σ=0.5 case is the cross-validation point for the high-density Run-6.

### Subtest — dt-convergence @ E=1100

10-min check that dt=0.01 gives < 1 % drift vs dt=0.005 at the highest
sweep energy. Locks dt for the Knudsen sweep.

### Run-8 — Knudsen sweep E ∈ {700, 800, 900, 1000, 1100} eV

5 pairs (10 runs) at dx=0.30 with the new observables enabled. Bracket
the regime where Knudsen et al. claim classical-WP stopping-power
convergence. Both estimators evaluated.

### Run-9 — 25 eV classical companion (post-Knudsen)

Closes the low-v end of the stopping-power curve with a matched-energy
classical pair to the existing 25 eV WP. Deferred to last because it
costs ~6 h alone.

### Final — Retroactive + rollup

Apply the new retroactive Python Knudsen pipeline to existing
E={50, 300, 600, 1500} eV pairs to harvest their `<|k|²>(t)` without
re-running. Combine with new-observable runs to produce the final
stopping-power vs v curve (12 anchors, 3 traces per WP run, Lindhard
overlay, shaded paper-claim band) in
`docs/reports/2026-05-21-meeting-emilio/figures/stopping_power_v2.png`.

## Sequencing graph

```
Infra-0  Gmail MCP OAuth          ┐
Infra-1..10  observables + helpers ├── parallel dev (~16 h human-driven)
                                   │
       ▼
Run-1 (25 eV WP)  ──── after Infra-1..10 ───┐
                                             │
Run-2 (free) + Run-2b (toy) + Run-3 (v2)  ───┼─ parallel on GPU 0&1
                                             │
       ▼ comparison plots + email
                                             │
Run-4 prep (2 GSes) ──┐                      │
Run-4 (extra-states) ─┘                      │
                                             │
Run-5 (20 eV pair) ──── GPU 0&1 concurrent ──┤
                                             │
Run-6 prep (L=30 GS) ─┐                      │
Run-6 (highdens pair) ─┘                     │
                                             │
Run-7 (σ-sweep pairs of 2) ──────────────────┤
                                             │
Subtest (dt-convergence) ──┐                 │
Run-8 (Knudsen sweep) ─────┘                 │
                                             │
Run-9 (25 eV classical, post-Knudsen) ───────┤
                                             │
Final retroactive analyses ──────────────────┤
                                             ▼
                                Final stopping-power rollup + email
```

User-priority ordering (which families to start first):
1 → 2/2b/3 → 4 → 5 → 6 → 7 → 8 → 9 → final.

## Compute & disk budget

- ~25-30 h GPU wall on 2 A30s running pairs concurrently
- ~30 GB disk for new VTI + observables (300 frames × 6-12 MB × ~20 runs)
- `df -h /local` check before launch
- Long-running classicals (Run-5, Run-9, parts of Run-8) ~6 h each
- No SLURM/wall-time enforcement; workstation A30s with `nohup`-based
  launches survive ssh disconnect

## Validation gates

Per `.claude/rules/testing.md` and `.claude/skills/physics-correctness`:

- **Infra-4, Infra-5**: known-case Heisenberg test required before any
  production run uses the new observables. `σ_r=5 ⇒ σ_p=0.1 Bohr⁻¹`;
  `<p_z>=k₀=2.711 Bohr⁻¹`; `E_kin ≈ 100 eV`.
- **Infra-1, 2, 3**: smoke check on a tiny pair of frames before
  full-run registration.
- **Run-2 vs Run-2b vs analytic**: σ_r²(t) agreement to ≤ 1 % across
  all three traces — validates injector + propagator at the same time.
- **Subtest**: dt-convergence < 1 % drift before locking Knudsen dt.
- **Every run**: standard verification (energy drift, cod_z slope,
  density_l2(0)=0, momentum peak at k₀) per
  `scripts/verify_smoke_outputs.py`.

## Open TODOs (deferred / off-critical-path)

- **σ→∞ bigger-box test** (task #31) — analytical claim needs L=100+
  box to test directly. Deferred post-meeting; meanwhile the σ=8
  data point + analytical argument in this plan addresses it
  qualitatively.
- **Read KL divergence + 2-Wasserstein theory** (task #29) —
  understanding the distance metrics chosen for jellium-vs-free
  comparison. Helps interpret KL(t) plot output.
- **Fix latent MomentumDistribution MPI-reduce** (task #30) — currently
  writes per-rank to shared CSV without `all_reduce`. Fine if all
  runs remain single-rank, latent bug at scale. Fix in Infra-4.
- **25 eV classical companion** (task #25) — explicitly scheduled
  post-Knudsen-sweep per user direction.
- **Workstation watchdog** (task #32) — optional 5-line script for
  thermal/Xid monitoring during the multi-day campaign.

## Attribution

- Project design: chiddukanna (skcb2) + meeting context with Emilio
  Artacho (14-05-2026 meeting).
- Stopping-power-via-momentum-spread method: Knudsen et al., arXiv
  2605.12854 (*Ultrafast electron dynamics of electron-irradiated
  graphene*) — to be added to `docs/sources/` as a source note.
- INQ TDDFT: alphataubio.com/inq, header-only C++17 (gitignored at
  `inq/`).
- Pre-existing INQ MomentumDistribution observable (jellium runs since
  2026-05-08) provides retroactive `<|k|²>(t)` extraction without
  re-running.
- Magic-number reference: `docs/sources/free-electron-gas-magic-numbers.md`.

## Approvals required

- [x] Design tree locked (user, 2026-05-17 conversation)
- [ ] Gmail MCP OAuth (user, in terminal)
- [ ] Approval to start Infra-1..10 development
- [ ] Approval to start Run-1
- [ ] dt-convergence subtest result (gates Run-8)

## Why this plan exists

Per `CLAUDE.md`: "Before any substantive implementation, create or
update a plan in `docs/plans/`." This plan is the single source of
truth for the campaign; the rolling handover tracks status, the
design journal entry captures motivations + open questions, and
per-run journal entries (created post-run) record observations.
