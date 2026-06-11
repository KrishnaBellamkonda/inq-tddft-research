# Rule: Scientific Grounding

Apply to: `inq-stack/`, `ResearchProject/`, `Tutorial/`, `QuantumKickExtension/`, `docs/reports/`, `docs/plans/`, `docs/handovers/`, `docs/sources/`

PROCEDURE (writing a `docs/sources/<author-year>.md` note + its template + the
search strategy) lives in the `literature-review` skill. This rule is the
always-on policy.

## Invariants

1. Ground scientific, numerical, and algorithmic claims in trustworthy sources:
   peer-reviewed papers > authoritative textbooks (Parr & Yang; Engel &
   Dreizler; Marx & Hutter) > official INQ/libxc/pseudopod/spglib docs >
   reputable, clearly-authored lecture notes. Treat forums, blogs, unreviewed
   snippets, and vendor marketing as low-trust — verify independently.
2. Distinguish direct source statements from your own inferences; label
   inferences explicitly ("Inference: ..." / "This suggests...").
3. When uncertain about a physical claim, say so and propose how to verify.
4. Credit adapted formulas/algorithms near the code (comment) and in the
   plan/handover.
5. DFT/TDDFT: justify functional (e.g. PBE for metals, LDA for jellium),
   pseudopotential (norm-conserving vs PAW; pseudopod database), and
   time-step/energy-cutoff choices by literature or convergence tests.
