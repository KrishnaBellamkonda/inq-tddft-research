---
name: scientific-tables
description: Use when adding a data/parameter/results TABLE to a presentation slide or report — ground-state ledgers (E_GS, E_total(0), E_H(0)), run/parameter tables (jellium + projectile), stopping-power result tables. Produces a native editable PPTX table by default (clear symbol-or-quantity headers with units, coloured header row, black text, no zebra) and a LaTeX/mathtext-rendered PNG fallback when symbols render poorly natively. Self-contained builder: make_table.py.
---

# Scientific tables — clear, editable, symbol-headed

The single source of truth for **table making** in decks and reports. It LAYERS
ON `scientific-figures` §6 (header row is the only coloured row; body rows white,
no zebra striping) and the AE `assertive-evidence-presentation` colour rule (all
text black). Everything here is mandatory unless a rule says otherwise.

Always: venv python (`/local/data/public/skcb2/tddft/venv/bin/python3`); the
skill-local builder is `make_table.py` (ships with this skill — import it).

## The rule set

1. **Headers name the quantity precisely.** Use the **symbol** when it is
   standard and widely recognised (`S`, `E_total(0)`, `E_H(0)`, `r_s`, `σ`,
   `ω_p`, `Z`, `v`); otherwise spell out the **quantity name** ("Interior
   density", "Box length"). Never a bare ambiguous label ("value", "x").
2. **Units are mandatory in the header**, in parentheses: `S (eV/Bohr)`,
   `E_total(0) (eV)`, `L_z (Bohr)`. Round the cell numbers to **2 s.f.** by
   default (project number-rounding rule); keep a column's s.f. consistent.
3. **Header row coloured, body white.** One pale accent fill on the header row +
   bold; body rows plain white, thin light-grey rules, **no zebra striping, no
   per-cell colour** (`scientific-figures` §6).
4. **All text black** (AE no-grey rule) — header and body alike.
5. **Prefer a NATIVE editable table in the slide.** A native PPTX table lets the
   user fix a number in PowerPoint without a rebuild — always the default.
6. **Fall back to a rendered PNG only when native symbol rendering is poor** —
   Greek letters, subscripts/superscripts, fractions, or unit exponents that a
   native table box mangles. Then render the table to a high-DPI PNG and insert
   it as a picture. State in the slide/report that this table is an image.

## Which path — decision

| Situation | Path | Function |
|---|---|---|
| Plain text + simple units (`eV`, `Bohr`, integers) | **native** (default) | `add_native_table` |
| Headers/cells need real maths (`$\sigma_\mathrm{WP}$`, `$e/\mathrm{Bohr}^3$`, fractions) | **PNG** | `table_to_png(engine="mpl")` |
| Publication-grade booktabs typography wanted | **PNG (LaTeX)** | `table_to_png(engine="latex")` |

`engine="mpl"` (matplotlib mathtext) is the robust default PNG path — no external
LaTeX needed. `engine="latex"` uses `pdflatex` + `booktabs`/`siunitx` and a PDF
rasteriser (`pdftoppm` or ImageMagick `convert`); higher quality, more deps.

## Usage

```python
import sys; sys.path.insert(0, ".claude/skills/scientific-tables")
import make_table as mt

# (a) native, editable table onto a python-pptx slide
mt.add_native_table(
    slide,
    header=["Run", "E_total(0) (eV)", "S (eV/Bohr)"],
    rows=[["23 eV", "-45.76", "0.021"],
          ["50 eV", "-45.71", "0.038"]],
    x=1.0, y=1.5, w=8.0)          # Inches; h auto if omitted

# (b) maths-heavy header -> PNG, then add_picture it yourself
mt.table_to_png(
    header=["Quantity", r"$E_\mathrm{total}(0)$ (eV)", r"$S$ (eV/Bohr)"],
    rows=[[r"$\sigma_\mathrm{WP}=0.5$", "-45.76", "0.021"]],
    out_png="assets/gs_ledger.png", engine="mpl")   # or engine="latex"
```

Cells are inserted **as given** — format numbers to 2 s.f. before calling. The
header row and thin borders are styled for you; do not re-colour body cells.

## Verify (before shipping any table)

Run `make_table.py --demo` (writes `/tmp/demo_mpl.png`, `/tmp/demo_latex.png`,
and a native-table smoke via the snippet in the module). Then check: symbol or
quantity in every header, units present, 2 s.f. numbers, header the only coloured
row, black text. A table failing any of these is not ready.

## Files

| File | Role |
|---|---|
| `make_table.py` | builder — `add_native_table`, `table_to_png` (mpl / latex) |
