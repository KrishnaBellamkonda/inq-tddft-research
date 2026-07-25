#!/usr/bin/env python3
"""Scientific-table helpers for presentations and reports.

Two rendering paths, matching the house rule (native editable table by default;
LaTeX-rendered PNG only when symbol rendering in a native table is poor):

  1. add_native_table(slide, header, rows, ...)
        Insert a *native, editable* python-pptx table onto a slide. Header row is
        accent-filled + bold black; body rows are plain white with thin light
        rules (scientific-figures SKILL §6). Text stays black (AE no-grey rule).
        Use this by DEFAULT — the user can edit numbers in PowerPoint directly.

  2. table_to_png(header, rows, out_png, engine="mpl", ...)
        Render the table to a high-DPI PNG when the header/cells carry maths that
        a native table renders poorly (Greek, subscripts, fractions, units with
        exponents). engine="mpl" (default) uses matplotlib mathtext — no external
        LaTeX install needed, robust everywhere. engine="latex" uses a real
        booktabs tabular via pdflatex (higher typographic quality; needs a LaTeX
        install + a PDF->PNG rasteriser).

Header convention (both paths): use the SYMBOL when it is standard and widely
recognised (S, E_total, r_s, sigma, omega_p), otherwise the spelled-out QUANTITY
name. Always carry units in the header cell, e.g. "S (eV/Bohr)".

Run:  /local/data/public/skcb2/tddft/venv/bin/python3 make_table.py --demo
"""
from __future__ import annotations
import os
import subprocess
import tempfile

# ---- house style ----------------------------------------------------------
FONT = "Calibri"
BLACK = (0x00, 0x00, 0x00)
WHITE = (0xFF, 0xFF, 0xFF)
HEADER_FILL = (0xD9, 0xE1, 0xF2)   # muted pale blue accent (header row only)
RULE = (0xBF, 0xBF, 0xBF)          # thin light-grey cell rules


# ==========================================================================
# 1. Native python-pptx table (default, editable)
# ==========================================================================
def add_native_table(slide, header, rows, x, y, w, h=None,
                     font=FONT, body_pt=14, header_pt=14,
                     header_fill=HEADER_FILL, col_widths=None):
    """Add a styled, editable native table to a python-pptx *slide*.

    header : list[str]           column headers (symbol-or-quantity + units)
    rows   : list[list[str]]     body cells (already formatted to 2 s.f. by you)
    x,y,w,h: Inches (floats). h defaults to a sensible per-row height.
    col_widths : optional list[float] Inches per column (else even split).

    Returns the inserted GraphicFrame.
    """
    from pptx.util import Inches, Pt
    from pptx.dml.color import RGBColor
    from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

    nrows, ncols = len(rows) + 1, len(header)
    if h is None:
        h = 0.34 * nrows
    gf = slide.shapes.add_table(nrows, ncols,
                                Inches(x), Inches(y), Inches(w), Inches(h))
    tbl = gf.table
    tbl.first_row = True          # header styling on
    tbl.horz_banding = False      # NO zebra striping (house rule)

    if col_widths:
        for j, cw in enumerate(col_widths):
            tbl.columns[j].width = Inches(cw)

    def _style(cell, text, *, bold, fill):
        cell.fill.solid()
        cell.fill.fore_color.rgb = RGBColor(*fill)
        cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        cell.margin_left = Inches(0.06); cell.margin_right = Inches(0.06)
        cell.margin_top = Inches(0.02); cell.margin_bottom = Inches(0.02)
        tf = cell.text_frame
        tf.word_wrap = True
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        run = p.add_run(); run.text = str(text)
        run.font.name = font
        run.font.size = Pt(header_pt if bold else body_pt)
        run.font.bold = bold
        run.font.color.rgb = RGBColor(*BLACK)          # always black text

    for j, htext in enumerate(header):
        _style(tbl.cell(0, j), htext, bold=True, fill=header_fill)
    for i, row in enumerate(rows, start=1):
        for j, ctext in enumerate(row):
            _style(tbl.cell(i, j), ctext, bold=False, fill=WHITE)

    _thin_borders(tbl, RULE)
    return gf


def _thin_borders(tbl, rgb):
    """Apply thin uniform light-grey borders to every cell (python-pptx has no
    high-level border API, so we edit the cell XML directly)."""
    from pptx.oxml.ns import qn
    hexc = "%02X%02X%02X" % rgb
    for row in tbl.rows:
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            for tag in ("a:lnL", "a:lnR", "a:lnT", "a:lnB"):
                for old in tcPr.findall(qn(tag)):
                    tcPr.remove(old)
            for tag in ("a:lnL", "a:lnR", "a:lnT", "a:lnB"):
                ln = tcPr.makeelement(qn(tag),
                                      {"w": "6350", "cap": "flat"})  # 0.5 pt
                fill = ln.makeelement(qn("a:solidFill"), {})
                clr = ln.makeelement(qn("a:srgbClr"), {"val": hexc})
                fill.append(clr); ln.append(fill); tcPr.append(ln)


# ==========================================================================
# 2. Table -> PNG (fallback for heavy maths in headers/cells)
# ==========================================================================
def table_to_png(header, rows, out_png, engine="mpl", dpi=300,
                 col_align=None, fontsize=16):
    """Render a table to PNG. Cells may contain matplotlib/LaTeX maths in $...$.

    engine="mpl"   : matplotlib mathtext (default, no external deps).
    engine="latex" : real booktabs tabular via pdflatex (needs LaTeX + a PDF
                     rasteriser: pdftoppm or ImageMagick `convert`).
    """
    if engine == "latex":
        return _table_to_png_latex(header, rows, out_png, dpi=dpi,
                                   col_align=col_align)
    return _table_to_png_mpl(header, rows, out_png, dpi=dpi, fontsize=fontsize)


def _table_to_png_mpl(header, rows, out_png, dpi=300, fontsize=16):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ncols = len(header)
    nrows = len(rows) + 1
    fig_w = max(1.6 * ncols, 4.0)
    fig_h = 0.5 * nrows + 0.3
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.axis("off")
    cells = [list(header)] + [list(r) for r in rows]
    tbl = ax.table(cellText=cells, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.scale(1.0, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#BFBFBF")
        cell.set_linewidth(0.6)
        cell.get_text().set_color("black")
        if r == 0:
            cell.set_facecolor("#D9E1F2")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor("white")
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    plt.close(fig)
    return out_png


def _table_to_png_latex(header, rows, out_png, dpi=300, col_align=None):
    align = col_align or ("l" + "c" * (len(header) - 1))
    def esc(s):  # keep $...$ maths verbatim; escape stray % and &
        return str(s).replace("%", r"\%")
    lines = [r"\documentclass[border=4pt]{standalone}",
             r"\usepackage{booktabs}\usepackage{amsmath}\usepackage{siunitx}",
             r"\begin{document}",
             r"\begin{tabular}{%s}" % align, r"\toprule",
             " & ".join(r"\textbf{%s}" % esc(h) for h in header) + r" \\",
             r"\midrule"]
    for row in rows:
        lines.append(" & ".join(esc(c) for c in row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{document}"]
    tex = "\n".join(lines)

    with tempfile.TemporaryDirectory() as d:
        texf = os.path.join(d, "t.tex")
        with open(texf, "w") as fh:
            fh.write(tex)
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "t.tex"],
                       cwd=d, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        pdf = os.path.join(d, "t.pdf")
        # rasterise: prefer pdftoppm, else ImageMagick convert
        if _have("pdftoppm"):
            base = os.path.join(d, "out")
            subprocess.run(["pdftoppm", "-png", "-r", str(dpi), pdf, base],
                           check=True)
            png = base + "-1.png"
        elif _have("convert"):
            png = os.path.join(d, "out.png")
            subprocess.run(["convert", "-density", str(dpi), pdf,
                            "-quality", "100", png], check=True)
        else:
            raise RuntimeError("engine='latex' needs pdftoppm or ImageMagick "
                               "convert to rasterise; use engine='mpl' instead.")
        os.replace(png, out_png)
    return out_png


def _have(prog):
    from shutil import which
    return which(prog) is not None


if __name__ == "__main__":
    import sys
    if "--demo" in sys.argv:
        hdr = ["Quantity", r"$E_\mathrm{total}(0)$ (eV)", r"$S$ (eV/Bohr)"]
        rows = [["23 eV run", "-45.76", "0.021"],
                ["50 eV run", "-45.71", "0.038"]]
        table_to_png(hdr, rows, "/tmp/demo_mpl.png", engine="mpl")
        print("wrote /tmp/demo_mpl.png")
        try:
            table_to_png(hdr, rows, "/tmp/demo_latex.png", engine="latex")
            print("wrote /tmp/demo_latex.png")
        except Exception as e:
            print("latex engine unavailable:", e)
