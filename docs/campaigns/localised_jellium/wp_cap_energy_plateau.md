---
title: "WP + localised-jellium CAP energy-plateau diagnostic"
status: complete
launched: 2026-07-22
completed: 2026-07-23
engine: inq-study
gpu: 0 (autonomous, setsid-detached)
family: "[wp-cap-energy-plateau]"
headline: "no-CAP plateau sits 86 eV ABOVE the CAP plateau — the trapped, should-have-radiated energy"
---

# WP-CAP energy-plateau diagnostic

## Question

High-density runs show a plateauing `energy_total`, but the retained energy at
the plateau seems too large for the system. **Does the CAP change the energy
scale / plateau of a WP-in-localised-jellium run?** Compare one run WITHOUT and
one WITH a two-sided CAP; the gap between their `energy_total` plateaus is the
energy radiated to the boundaries (which the CAP drains and the closed box keeps).

## System (locked)

- Localised jellium SLAB: 25×25 periodic face, 25 Bohr thick (half-width 12.5,
  faces ±12.5), centred in a 25×25×140 box. erfc edge softening w = 0.5 Bohr.
- N = 102 electrons → r_s ≈ 3.32, E_F ≈ 4.5 eV, ħω_p ≈ 7.8 eV.
- Grid h = 0.5 Bohr; LDA; T ≈ 100 K smearing; inq-study engine.
- WP projectile: σ_WP = 1 Bohr, E = 100 eV (k₀ ≈ 2.711), mass 1, launched z=−20.5
  (8 Bohr from the −12.5 face) moving +z.
- CAP (run 2 only): two-sided sin² absorber, η = −0.7 Ha, 10 Bohr/side at the far
  ends z∈[±60,±70] (fractional mid ±0.4643, width 0.07143). Functional only on
  inq-study (stock inq drops the imaginary term — see
  memory `reference_stock_inq_cannot_compile_cap`).

## Runs

1. **vacuum warm-up** (`systems/vacuum/scripts/wp_traversal_energy`,
   non-interacting): single WP full traversal, no-CAP vs CAP. no-CAP conserves
   E_total (4.42 Ha); CAP drains it → 0. Validates the machinery + the CAP energy
   removal. DONE.
2. **jellium no-CAP** (`scripts/wp_cap_energy_plateau/wp`, WP_CAP_ETA=0), 100 a.u.
3. **jellium CAP** (WP_CAP_ETA=−0.7), 100 a.u.

Each records ALL KS energy components every step, the momentum distribution every
step, the projectile wavefunction every 10 steps, and density frames; checkpointed
every 200 steps + resumable. Per-run notebook (energy decomposition, momentum,
density GIFs) via `analyse.py`; headline no-CAP-vs-CAP plateau overlay via
`compare.py`.

## What to look for

1. Difference in energy scales between the two runs.
2. Whether `energy_total` plateaus (and at what value) in each.
3. Plateau gap = energy radiated to the boundaries.

## Orchestration

Autonomous, `setsid`-detached on GPU 0; emails at each stage. See
`docs/handovers/wp-localised-jellium-solving-cap.md` for launch/kill/resume.

## Results — examine in this order (why → evidence)

All four propagations completed autonomously (`run_completed = true`); the CAP run
finished 2026-07-23 19:55 (wall ≈ 13.6 h). Open the notebooks/figures below **top to
bottom** — the order retraces the reasoning from "is the system built right?" to the
final answer.

**0 · Why we did this.** The plateauing `energy_total` in high-density runs looked
too energetic — a periodic box has nowhere to shed the energy a wavepacket dumps
into / reflects off the slab. Test: identical WP-through-slab run **without** vs
**with** a CAP; if no-CAP sits above CAP, the gap is the trapped energy.

| # | Step (why) | Open this | Headline |
|---|---|---|---|
| 1 | **Is the slab + launch built right?** (shared start of both RT runs) | `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/gs/report/gs_report.ipynb` | E_GS = −830 Ha (−22586 eV), r_s 3.3; WP launches in vacuum (0.013 % of centre density), 8.0 Bohr from the face, force-free |
| 2 | **Does the machinery + CAP work at all?** (vacuum, no slab) — no-CAP | `/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/scripts/wp_traversal_energy/results/nocap/report/run_report.ipynb` | E_total conserved = 120.4 eV |
| 3 | vacuum — CAP drains it | `/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/scripts/wp_traversal_energy/results/cap/report/run_report.ipynb` | E_total → 83.8 eV; **gap 36.6 eV** removed by the CAP |
| 4 | **The real experiment** — jellium, no-CAP | `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/wp/results/nocap/report/run_report.ipynb` | closed box, energy plateaus at −22462 eV |
| 5 | jellium — CAP | `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/wp/results/cap/report/run_report.ipynb` | ΔE_total −93.5 eV; plateau at −22548 eV |
| 6 | **The answer** — plateau overlay | `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/wp/results/jellium_energy_compare.png` | **no-CAP plateau is 86 eV ABOVE CAP** → 86 eV was trapped in the box |

**Finding.** The hypothesis holds: without an absorber the localised-jellium box
retains ≈ 86 eV that should have radiated away. The gap grows from the vacuum case
(36.6 eV) because the interacting slab scatters far more of the wavepacket than empty
space does. Both jellium runs share the identical GS of step 1, so the 86 eV is
attributable to the CAP alone — nothing in the setup differs.

*Provisional:* the 86 eV needs a physical decomposition (reflected-WP KE vs plasmon
excitation vs trapped bound state) before over-claiming — a `scientific-panel` read
of step 6 is the natural next move.

## Source & machinery paths

- GS run: `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/gs/run.cpp`
- RT run (both): `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/wp/run.cpp`
- Config header: `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/shared/configs/slab_n102_L25x25x140_w0p5.hpp`
- Orchestrator (whole chain): `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/orchestrate.sh`
- GS notebook builder: `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/gs/make_gs_report.py`
- Run log: `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/orchestrate.log`
- Handover: `/local/data/public/skcb2/tddft/docs/handovers/wp-localised-jellium-solving-cap.md`
