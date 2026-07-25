# ADR 0011 — Held-out split protocol for the ml-patterns agentic retry loop

- **Status:** accepted
- **Date:** 2026-07-01
- **Context scope:** `docs/campaigns/ml-patterns/pattern-finding-in-wp-classical-runs.md`
  (the autonomous induced-density ML-discovery campaign and any successor that
  uses an automated pipeline-config search)

## Context

The ml-patterns campaign is **fully autonomous** and, after T1, wraps its Rung-1
analysis in a **bounded agentic validation+retry loop (≤4 tries, keep best, then
proceed)** — the executor inspects the result, adjusts the analysis pipeline
(q-window, DMD rank, subtraction order, mode count…), and re-runs.

The campaign's headline claim rests on **parameter-free** predictions: the
form-factor q-ratio should follow `exp(−q²σ_pot²/2)` (σ_pot read from the UPF) and
the DMD wake frequency should equal `ω_p=√(4πn)` (fixed by r_s), each within a
**±20%** band. This parameter-free quality is the campaign's entire defence against
the "ML found what I wanted" failure that the deep research named as the dominant
risk (spurious discovery from un-normalised fields + tunable analysis choices).

The naive loop objective — *"retry pipeline configs until agreement is within ±20%,
keep the best agreement"* — silently re-introduces a fitted parameter: **the
pipeline configuration itself**, selected to maximise agreement on the *same* cells
then cited as confirmation. That is p-hacking; the result would not survive review.
Reverting to "no retry / method-validity only" was considered but rejected — it
discards the legitimate need to fix an analysis config that *masks a real signal*
(wrong q-window, unconverged POD, miswindowed DMD).

## Decision

The retry loop **may** optimise for ±20% agreement, but under a **pre-registered
train/test (calibration/held-out) split**:

1. **Pin** a non-overlapping split of the analysis cells **in the campaign prompt**
   (not chosen at runtime): *form-factor cut* (E=100 eV) — calibration σ_WP∈{1,5},
   held-out σ_WP∈{0.5,3,8}; *wake cut* (σ_WP=5) — calibration = even-velocity-index
   energies, held-out = odd. Echoed in the notebook.
2. The ≤4-try loop tunes the pipeline config to maximise agreement **on the
   calibration cells only**.
3. The winning config is **frozen** and the **verdict is read exclusively from the
   held-out cells**; that held-out number is what the campaign reports.
4. **All ≤4 attempts are logged** (configs + scores), not just the winner.
5. CONFIRM / REFUTE / INCONCLUSIVE are all valid reported outcomes; a refute is
   never retried into a confirm. INCONCLUSIVE = method validity unreachable in ≤4.

## Consequences

- **Defensible discovery.** Config selection happens on data independent of the
  data that confirms/refutes, so "keep best of 4" is honest model selection, not
  cherry-picking — the parameter-free guarantee survives.
- **Keeps the user's intent.** The autonomous "retry, keep best, move on" behaviour
  is preserved; only the *scoring data* changes.
- **Binding on the executing agent.** A future agent must NOT "simplify" the loop
  into agreement-on-all-cells, and must NOT report agreement computed on the
  held-out split's tuning. This is encoded in the campaign `<guard_rails>` and
  `<rules>`; this ADR is the rationale a reader will otherwise wonder about.
- **Cost.** Splitting reduces the cells available for each role; with only 5
  σ-cells in the form-factor cut and ~7 in the wake cut, the split is **within-cut
  but coarse** (2 calibration / 3 held-out σ; even/odd-velocity-index energies), and
  the **shared** pipeline config is tuned on the *union* of both cuts' calibration
  cells. The held-out form-factor test thus rests on 3 σ-points — thin, so the
  verdict is reported with that caveat, and a borderline result is INCONCLUSIVE
  rather than a forced confirm/refute.
