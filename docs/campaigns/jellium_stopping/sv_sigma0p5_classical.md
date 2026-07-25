---
id: jstop-sv-sigma0p5-classical
area: jellium_stopping
title: "Jellium S(v) - sigma=0.5 classical projectile sweep"
status: done
hypothesis: "A classical -1 / sigma=0.5 Bohr erf-Gaussian electron's stopping S(v) in r_s=5.69 jellium traces the Barkas-sign crossover and trends toward Lindhard at high v."
handover: docs/handovers/overnight-gaussian-classical-jellium.md
tasks:
  - { name: "6-velocity S(v) classical runs", done: true }
  - { name: "S(v) extraction + executed notebook", done: true }
blocked_reason: ""
---

# Jellium S(v) - sigma=0.5 classical projectile sweep

<!-- Retroactive campaign record (backfilled 2026-06-22). This work predates the
campaigns skill; authoritative detail lives in the linked plan/handover below.
Expand into the full template (docs/campaigns/template.md) via the campaigns skill
if this campaign is revived/extended. -->

**Plan / handover:** docs/plans/overnight-gaussian-classical-jellium.md ; docs/handovers/overnight-gaussian-classical-jellium.md (branch `overnight-gaussian-classical`)

> NOTE: classical S(v) is complete (6 velocities, Barkas crossover). The
> quantum WP / loss-function *production* extension was COST-BLOCKED (~20k steps,
> GPU budget) and k-points deferred - that extension feeds `cap-jellium-loss-function`.
> Low-velocity points are a noisy friction tail (graceful degradation).
