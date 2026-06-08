---
name: literature-review
description: Use when grounding scientific claims in trustworthy sources, choosing functionals/pseudopotentials/timesteps, citing benchmarks, or writing a source note in `docs/sources/`. Required whenever a physical, numerical, or algorithmic choice in `inq-stack/`, `ResearchProject/`, or `Tutorial/` needs justification.
---

# Literature Review

## 1. Identify what needs grounding

Before writing or coding, list the specific claims that require a citation:
- Physical approximations (why this functional, why this pseudopotential)
- Numerical parameters (cutoffs, time steps, convergence criteria)
- Validation benchmarks (what energy / spectrum does the literature report?)
- Algorithm choices (why ETRS vs Crank-Nicolson, why Broyden mixing?)

## 2. Gather sources

Preferred source types (in order):
1. Peer-reviewed papers (PRL, PRB, J. Chem. Phys., J. Phys.: Condens. Matter, …)
2. Authoritative textbooks (Parr & Yang; Engel & Dreizler; Marx & Hutter)
3. Official INQ docs: `docs/inq-docs/` (local mirror) or alphataubio.com/inq
4. libxc documentation for functional choices
5. Reputable university lecture notes (clearly authored)

Primary reference for this project:
- Santervás-Arranz, Stengel, Artacho, *Phys. Rev. Research* **7**, 033292 (2025) — quantum kick / TDDFT energy absorption.

## 3. Distinguish source vs. inference

- **Direct statement**: "According to [source], the PBE functional gives X."
- **Inference**: "Inference: this suggests Y would also apply, but this has not been verified here."

Never present an inference as a direct source statement.

## 4. Write a source note

Create `/local/data/public/skcb2/tddft/docs/sources/<author-year-keyword>.md`:

```md
# Source: <Author Year — Title fragment>

## Full citation
<author(s)>, <title>, <journal>, <volume>, <page/DOI>, <year>

## Relevance to this project
<one paragraph>

## Key claims used
- Claim 1 (page/section X)
- Claim 2 (page/section Y)

## Limitations / uncertainties
<anything the source does not cover or explicitly caveats>

## Cross-references
- Plan: docs/plans/<task>.md
- Code: inq-stack/include/inqkit/<module>/<file>:<line> (if code adapted)
       or ResearchProject/systems/<material>/<task>/<file>:<line>
```

## 5. Propagate source awareness

After writing the source note:
- Reference it in the relevant plan (`docs/plans/`)
- Reference it in the handover (`docs/handovers/`)
- Add a brief comment in code where logic is adapted: `// Adapted from <Author Year>, <equation/section>`
- Cite it in any report or manuscript

## Common sources for INQ / TDDFT projects

| Topic | Key references |
|---|---|
| TDDFT theory | Runge & Gross (1984); Casida & Huix-Rotllant Rev. (2012) |
| ETRS propagator | Castro et al., J. Chem. Phys. **121**, 3425 (2004) |
| Crank-Nicolson | Standard numerical methods texts |
| PBE functional | Perdew, Burke, Ernzerhof, PRL **77**, 3865 (1996) |
| LDA functional | Perdew & Zunger, PRB **23**, 5048 (1981) |
| Quantum kick / electronic stopping | Artacho, J. Phys.: Condens. Matter **19**, 275211 (2007) |
| INQ code | See alphataubio.com/inq |
