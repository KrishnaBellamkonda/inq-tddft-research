# Skill: Report Writing

Use this skill when converting plans, handovers, and source notes into report-ready prose, or when drafting manuscripts, slides, or figure captions.

---

## Protocol

### 1. Assemble source material first

Before drafting, collect:
- The relevant handovers from `docs/handovers/`
- The relevant source notes from `docs/sources/`
- The relevant validation notes from `docs/validation/`
- Any existing report fragments in `docs/reports/`

### 2. Structure

Follow the standard scientific paper structure unless the target format dictates otherwise:

1. **Introduction** — motivation, gap in literature, what this work does
2. **Methods** — computational setup (cell, functional, cutoff, k-points, time step, pseudopotentials, code/version)
3. **Results** — key findings with figures; cite source notes for physical interpretation
4. **Discussion** — comparison to literature, limitations, caveats
5. **Conclusions**
6. **References**

### 3. Attribution rules

- Every adapted method, formula, or algorithm must cite the original source.
- Distinguish "method from X, applied here" from "method developed in this work".
- Never imply original authorship for borrowed ideas.
- Use the target journal's citation style. If none specified, use author-year (APA) in notes, and numeric in manuscripts.

### 4. Computational methods section template

```
All calculations were performed using INQ [INQ-ref], a GPU-accelerated
DFT/TDDFT code. Ground-state calculations used the [LDA/PBE/...] exchange-
correlation functional [Perdew-ref]. Norm-conserving pseudopotentials from
the pseudopod library were employed. The plane-wave kinetic energy cutoff was
[X] Ry, and the simulation cell was [describe]. The Brillouin zone was sampled
using a [n×n×n] Monkhorst-Pack grid [MP-ref]. Real-time propagation used the
ETRS algorithm [Castro-2004] with a time step of [X] atomic units over [N]
steps ([total time] atomic units). Perturbations: [describe kick/laser].
```

### 5. Figure captions

- State what is plotted (quantity, units, axes).
- State the system and calculation parameters if not already in Methods.
- Mention if GPU/CPU agreement was verified.
- If comparing to a reference, cite it.

### 6. Uncertainty language

- "The calculation converged to X Ha." — only if convergence was verified.
- "The result is consistent with..." — qualified comparison.
- "We infer that..." — explicit inference label.
- Avoid "shows", "proves", "demonstrates" for unverified claims.

---

## Report files

Place drafts in `docs/reports/<topic>.md` or `docs/reports/<topic>.tex`.
Store figure-caption working notes in `docs/reports/figure-captions.md`.
