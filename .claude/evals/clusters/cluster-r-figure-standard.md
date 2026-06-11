# Eval: Cluster R — global figure standard (LOCKED 2026-06-11)

Single global standard = `inqview.visualisation.style` (ADR 0004). Every
figure-producing code path imports it; `report-figures` skill owns the standard
+ the 5 project annotation rules; `tufte` = timeless principles only.

## Programmatic evals (.claude/evals/programmatic/)

1. **Theme-import enforcement.** Grep every figure-producing Python path
   (`inq-stack/python/inqview/**`, run `analyse.py` templates in
   `ResearchProject/.../shared/python/`, comparison `scripts/`) for plotting:
   any module that calls `plt.subplots`/`savefig` MUST import
   `inqview.visualisation.style` (or go through a factory that does). FAIL on a
   rogue `plt.rcParams[...] = …` / `plt.style.use(...)` / hardcoded `figsize`
   outside the theme. (Allowlist the theme module itself.)
2. **Units drift.** The units stated in the `report-figures` skill (annotation
   rule 5) and any prose == the theme module's canonical unit map
   (eV/Bohr, eV, Bohr, fs, Bohr⁻¹). FAIL on disagreement.
3. **No-restate structural check.** `tufte` SKILL.md must NOT contain the column-
   widths table or the 5 project annotation rules (they moved); grep sentinel.

## Behavioural eval (.claude/evals — LLM-judge)

- Trigger: "make a stopping-power comparison plot" / "render the density wake"
  → the model uses the canonical theme (semantic cmap role, figure factory),
  not ad-hoc styling. (Auto-plot path: theme yes, interactive grill no.)
- Trigger: "produce the report figure panel for Fig 3" → full `report-figures`
  5-phase workflow fires (grill → plot → minipage → compile).
- Distinct-core check: a pure critique request ("is this chart honest?") →
  `tufte`; a prose request ("draft the methods section") → `report-writing`.

## What must NOT break

- The 3 skills stay separate (distinct cores). This eval does not merge them.
- Only the *standard's facts* are centralised; each skill's workflow is untouched.
