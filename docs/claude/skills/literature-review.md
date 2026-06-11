# Skill: Literature Review

Use this skill when gathering sources, extracting evidence, tracking uncertainty, or creating source notes in `docs/sources/`.

---

## Protocol

### 1. Identify what needs grounding

Before writing, list the specific claims that require a citation or source:
- Physical approximations (why this functional, why this pseudopotential)
- Numerical parameters (cutoffs, time steps, convergence criteria)
- Validation benchmarks (what energy / spectrum does the literature report?)
- Algorithm choices (why ETRS vs Crank-Nicolson, why Broyden mixing?)

### 2. Gather sources

Preferred source types (in order):
1. Peer-reviewed papers (PRL, PRB, J. Chem. Phys., J. Phys.: Condens. Matter, etc.)
2. Authoritative textbooks
3. Official INQ docs: `docs/inq-docs/` (local mirror) or alphataubio.com/inq
4. libxc documentation for functional choices
5. Reputable university lecture notes (clearly authored)

# TODO: There are no primary references now. However, there is a library of related documents being gathered in literature/ folder. Also, new documents uploaded to the drive folder can be downloaded to this folder. This becomes the primary repository of search. If not found, then the internet. Sometimes references from these texts can be helpful. 
For this project, the primary reference is:
- Santervás-Arranz, Stengel, Artacho, *Phys. Rev. Research* 7, 033292 (2025) — quantum kick / TDDFT energy absorption

### 3. Distinguish source vs. inference

- **Direct statement**: "According to [source], the PBE functional gives X."
- **Inference**: "Inference: this suggests Y would also apply, but this has not been verified here."

Never present an inference as a direct source statement.

### 4. Write a source note

Create `docs/sources/<author-year-keyword>.md`:

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
- Code: inq/src/<file>:<line> (if code was adapted)
```

### 5. Propagate source awareness

After writing the source note:
- Reference it in the relevant plan (`docs/plans/`)
- Reference it in the handover (`docs/handovers/`)
- Add a brief comment in code where logic is adapted: `// Adapted from <Author Year>, <equation/section>`
- Cite it in any report or manuscript

---

## Common sources for INQ / TDDFT projects

| Topic | Key references |
|---|---|
| TDDFT theory | Runge & Gross (1984); Casida & Huix-Rotllant Rev. (2012) |
| ETRS propagator | Castro et al., J. Chem. Phys. 121, 3425 (2004) |
| Crank-Nicolson | Standard numerical methods texts |
| PBE functional | Perdew, Burke, Ernzerhof, PRL 77, 3865 (1996) |
| LDA functional | Perdew & Zunger, PRB 23, 5048 (1981) |
| Quantum kick / electronic stopping | Artacho, J. Phys.: Condens. Matter 19, 275211 (2007) |
| INQ code | See alphataubio.com/inq for INQ-specific references |
