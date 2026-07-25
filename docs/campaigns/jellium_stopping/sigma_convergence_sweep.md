---
id: jstop-sigma-convergence
area: jellium_stopping
title: "Jellium S(v) sigma-convergence sweep (sigma=0.15 -> 3.0)"
status: done
hypothesis: "As the classical-projectile width sigma -> point-charge, jellium S(v) converges toward the point-charge Lindhard stopping."
handover: docs/handovers/sigma-convergence-stopping.md
tasks:
  - { name: "sigma=0.15 run", done: true }
  - { name: "sigma=0.25 run", done: true }
  - { name: "sigma=0.35 run", done: true }
  - { name: "sigma=3.0 run", done: true }
  - { name: "convergence study (sigma_sweep_report.py + sv_convergence figures)", done: true }
blocked_reason: ""
---

# Jellium S(v) sigma-convergence sweep (sigma=0.15 -> 3.0)

<!-- Retroactive campaign record (backfilled 2026-06-22). This work predates the
campaigns skill; authoritative detail lives in the linked plan/handover below.
Expand into the full template (docs/campaigns/template.md) via the campaigns skill
if this campaign is revived/extended. -->

**Plan / handover:** docs/plans/sigma-convergence-stopping.md ; docs/handovers/sigma-convergence-stopping.md

> NOTE: sigma=0.5 is covered by the separate `jstop-sv-sigma0p5-classical`
> campaign. sigma=3.0 launched 2026-06-15 (detached). Hypotheses folder:
> `ResearchProject/systems/jellium/hypotheses/06_sigma_convergence/`.
