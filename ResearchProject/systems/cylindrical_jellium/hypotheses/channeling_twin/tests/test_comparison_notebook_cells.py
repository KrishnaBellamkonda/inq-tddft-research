#!/usr/bin/env python3
"""Static guards on the comparison-notebook builder's generated code cells.

WHY THIS EXISTS
---------------
build_comparison_notebook.py writes Python SOURCE as strings and executes it, so
its failure modes are not caught by importing it. Two of them are silent:

1. A non-raw string literal containing a LaTeX command whose escape is VALID
   Python. ``"$\\approx 0$"`` inside a generated cell becomes ``"$\x07pprox 0$"``
   — the bell character — and matplotlib then dies with a mathtext parse error
   pointing at "pprox", 200 lines into an executed notebook. ``\\alpha``,
   ``\\beta``, ``\\nu``, ``\\tau``, ``\\rho``, ``\\theta``, ``\\varphi`` are all in
   this class (\\a \\b \\n \\t \\r \\v \\f). This bit the builder once already
   (2026-08-01) and cost a full notebook execution to find.
2. A mathtext command matplotlib does not implement (``\\le`` rather than
   ``\\leq``). Same cost, same lateness.

Both are cheap to catch statically, so they are.

Run:  <repo>/venv/bin/python3 -m pytest <this file> -q
"""
from __future__ import annotations

import io
import re
import sys
import tokenize
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import build_comparison_notebook as B  # noqa: E402

# LaTeX commands matplotlib's mathtext does NOT provide but that a LaTeX habit
# reaches for. Extend as they turn up.
UNSUPPORTED_MATHTEXT = (r"\le ", r"\ge ", r"\le$", r"\ge$", r"\text{", r"\mbox{")


def code_cells() -> list[tuple[int, str]]:
    return [(i, src) for i, (kind, src) in enumerate(B.cells("wp", "classical"))
            if kind == "code"]


def test_every_code_cell_compiles():
    for i, src in code_cells():
        compile(src, f"<cell{i}>", "exec")     # raises SyntaxError on failure


def test_no_non_raw_backslash_literals():
    """Every generated string literal containing a backslash must be raw.

    The whitelist is for literals whose backslash IS meant as a Python escape
    (a real newline in a title), not as LaTeX.
    """
    WHITELIST = (r"\n",)
    offenders = []
    for i, src in code_cells():
        for tok in tokenize.generate_tokens(io.StringIO(src).readline):
            if tok.type != tokenize.STRING:
                continue
            text = tok.string
            body_start = text.index(text.lstrip("rRbBfFuU")[0]) if text[0] not in "\"'" else 0
            prefix = text[:body_start]
            if "r" in prefix.lower() or "\\" not in text:
                continue
            stripped = text
            for w in WHITELIST:
                stripped = stripped.replace(w, "")
            if "\\" in stripped:
                offenders.append(f"cell{i} L{tok.start[0]}: {text[:100]}")
    assert not offenders, (
        "non-raw string literal(s) containing a backslash — a valid Python escape "
        "would silently eat the LaTeX command:\n  " + "\n  ".join(offenders))


def test_no_unsupported_mathtext_commands():
    offenders = []
    for i, src in code_cells():
        for cmd in UNSUPPORTED_MATHTEXT:
            # The builder writes \\le into its own source, which reaches the cell
            # as \le; search the CELL text.
            if cmd in src:
                for ln, line in enumerate(src.splitlines(), 1):
                    if cmd in line:
                        offenders.append(f"cell{i} L{ln}: {cmd.strip()!r} in {line.strip()[:90]}")
    assert not offenders, (
        "mathtext command(s) matplotlib does not implement:\n  " + "\n  ".join(offenders))


def test_math_spans_are_balanced():
    """An odd number of unescaped '$' in a label means an unterminated math span."""
    offenders = []
    for i, src in code_cells():
        for ln, line in enumerate(src.splitlines(), 1):
            if "$" not in line:
                continue
            # count $ inside string literals only, crudely but adequately
            n = len(re.findall(r"(?<!\\)\$", line))
            if n % 2:
                offenders.append(f"cell{i} L{ln}: {line.strip()[:90]}")
    assert not offenders, ("odd number of '$' — unterminated math span:\n  "
                           + "\n  ".join(offenders))


def test_the_three_verdict_blocks_are_present():
    """The notebook must still contain a figure block for each part of the aim.

    Result / premise / mechanism: dropping one silently would leave a notebook
    that looks complete and cannot support its own conclusion.
    """
    md = "\n".join(src for kind, src in B.cells("wp", "classical") if kind == "md")
    for marker in ("THE RESULT", "THE PREMISE", "THE MECHANISM",
                   "Correctness gates", "Verdict"):
        assert marker in md, f"the '{marker}' section is missing from the notebook"

    code = "\n".join(src for _, src in code_cells())
    for fig in ("03_stopping_bars", "05_channeling", "06_var_p_freeze", "09_gates"):
        assert fig in code, f"figure {fig} is no longer produced"
    # The density GIF is mandatory (.claude/rules/notebook-density-gif.md).
    assert "make_twin_density_matrix" in code
    assert "display(Image(filename=" in code, "the GIF must be DISPLAYED, not just written"


if __name__ == "__main__":
    import traceback
    fails = 0
    tests = [(n, o) for n, o in sorted(globals().items())
             if n.startswith("test_") and callable(o)]
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception:
            fails += 1
            print(f"  FAIL  {name}")
            traceback.print_exc()
    print(f"\n{len(tests) - fails}/{len(tests)} passed")
    raise SystemExit(1 if fails else 0)
