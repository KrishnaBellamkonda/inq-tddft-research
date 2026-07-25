# Rule: Number rounding / significant figures

Apply to: every reported number in `docs/reports/`, `docs/presentations/`,
`docs/handovers/`, `docs/validation/`, run `analyse.py`, notebook tables,
figure annotations/captions, and any tabulated result a human reads. Always on.

## The one rule

**Round reported numbers to 2 significant figures by default; 3 s.f. at most.**

- Headline / spoken / slide numbers: **2 s.f.** (e.g. `S = 2.7 eV/Bohr`,
  `E_SIE = 4.4 eV`, `ΔE = 68 eV`).
- Use **3 s.f.** only when 2 would hide a meaningful distinction — e.g. two
  quantities being differenced that agree to within their leading digits
  (`E_jellium(0) = −45.75 Ha` vs `E_GS = −45.76 Ha`), or a value whose third
  digit carries the physics.
- Never report a chain of 6–12 raw solver digits in human-facing text or tables.
  Keep full precision only inside code / intermediate variables; round at the
  point of *presentation*.

## Why

A presentation or report number with 8 digits reads as false precision and
buries the result a reader actually needs. 2 s.f. is the legible default; 3 s.f.
is the escape hatch for genuine near-equalities. (User decision, 2026-06-25, for
the supervisor presentations and all downstream reporting.)

## How to apply

- When tabulating an energy ledger, round each row consistently (same s.f. across
  a column) and show the differenced rows at the s.f. needed to see the gap.
- Carry units explicitly (Ha / eV) and convert with 1 Ha = 27.211 eV.
- If a number must keep more digits for reproducibility, put the full value in a
  code cell / provenance line, not in the headline table.
