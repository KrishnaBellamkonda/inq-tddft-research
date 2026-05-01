# Rule: Scientific Grounding

Apply to: `inq/src/`, `ResearchProject/`, `docs/reports/`, `docs/plans/`, `docs/handovers/`, `docs/sources/`

## Rules

1. All scientific explanations, implementation choices, and modelling recommendations must be grounded in trustworthy sources.

   Acceptable sources:
   - Peer-reviewed journal papers
   - Authoritative textbooks (e.g. Parr & Yang, Engel & Dreizler, Marx & Hutter)
   - Official INQ documentation (`docs/inq-docs/` mirror or alphataubio.com/inq)
   - Official documentation for libxc, pseudopod, spglib
   - Lecture notes from clearly-authored, reputable university courses

   Low-trust sources (verify independently):
   - Stack Overflow, forums, unreviewed snippets
   - Blog posts, vendor marketing

2. Distinguish direct source statements from your own inferences. Label inferences explicitly: "Inference: ..." or "This suggests...".

3. When uncertain about a physical claim, say so and propose how to verify.

4. Record important sources in `docs/sources/` and cross-reference them from plans, handovers, reports, and code comments.

5. If code structure, formulas, or algorithms are adapted from a paper, package, or prior note, credit the source near the adapted logic (code comment) and in the plan/handover.

6. For DFT/TDDFT claims specifically:
   - Functional choices should be justified by literature (e.g. PBE for metals, LDA for jellium)
   - Pseudopotential choices (norm-conserving vs PAW, pseudopod database) should be documented
   - Time-step and energy cutoff choices should reference convergence tests or literature values
