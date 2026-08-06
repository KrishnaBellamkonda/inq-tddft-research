#!/usr/bin/env python3
"""Static guards on EVERY bulk-jellium notebook builder's generated code cells.

WHY THIS EXISTS
---------------
These builders write Python SOURCE as strings and execute it, so their failure
modes survive an ordinary import and only surface deep inside a slow notebook
execution. Three are silent or very late:

1. A **non-raw** builder string containing a backslash. The escape is applied when
   the BUILDER is parsed, so the emitted cell is already wrong:
     * ``\\r`` in ``\\rangle`` / ``\\rm``  -> carriage return -> SyntaxError
     * ``\\n`` in an f-string annotation -> REAL newline -> unterminated literal
     * ``\\a`` in ``\\approx``            -> BEL -> compiles, renders WRONG
   ``\\a \\b \\f \\n \\r \\t \\v`` are all in this class, so ``\\alpha``, ``\\beta``,
   ``\\rangle``, ``\\rm``, ``\\tau``, ``\\nu``, ``\\theta``, ``\\varphi`` all bite.
2. A mathtext command matplotlib does not implement. ``\\frac12`` is valid Python
   AND valid TeX and still fails at render: mathtext wants ``\\frac{num}{den}``.
   No blocklist catches this, so the spans are parsed by matplotlib itself.
3. An unbalanced ``$`` span.

COST OF NOT HAVING THIS (2026-08-01, one day, three separate hits):
  * phase builder, ``\\rangle``/``\\rm``  -> 5 broken cells, one GIF-stage execution
  * phase builder, ``\\frac12``          -> caught here, minutes before it ran
  * run builders, ``\\n`` in annotate    -> **all 8 run-notebook array tasks failed**
The third landed in the ONE builder the first version of this file did not cover,
which is why it now covers every builder rather than the one that broke first.

Run:  venv/bin/python -m pytest \
        ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/tests/ -q
"""
from __future__ import annotations

import importlib.util
import io
import re
import sys
import tokenize
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
HYP = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))

FAMILIES = ["bulk_ks_stopping", "bulk_ks_stopping_sigma3",
            "bulk_ks_stopping_rs4", "bulk_ks_stopping_rs4_sigma3"]

# LaTeX matplotlib's mathtext does NOT provide but a LaTeX habit reaches for.
UNSUPPORTED_MATHTEXT = (r"\le ", r"\ge ", r"\le$", r"\ge$", r"\text{", r"\mbox{")


def _load(family: str, module: str):
    """Import a builder under a unique name (several are same-named clones)."""
    path = HYP / family / f"{module}.py"
    spec = importlib.util.spec_from_file_location(f"_{module}_{family}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_nb(family: str, half: str):
    return _load(family, "build_run_notebook").build(half)


def _phase_nb(family: str):
    return _load("bulk_ks_stopping", "build_phase_notebook").build(family)


# (label, thunk) over every notebook this hypotheses tree can emit.
BUILDS: list[tuple[str, object]] = (
    [(f"run:{f}/{h}", (lambda f=f, h=h: _run_nb(f, h)))
     for f in FAMILIES for h in ("wp", "classical")]
    + [(f"phase:{f}", (lambda f=f: _phase_nb(f))) for f in FAMILIES]
)
IDS = [lbl for lbl, _ in BUILDS]


def code_cells(thunk) -> list[tuple[int, str]]:
    nb = thunk()
    return [(i, c.source) for i, c in enumerate(nb.cells) if c.cell_type == "code"]


@pytest.mark.parametrize("label,thunk", BUILDS, ids=IDS)
def test_every_code_cell_compiles(label, thunk):
    """The bug that shipped twice: an emitted cell that is not valid Python."""
    for i, src in code_cells(thunk):
        compile(src, f"<{label} cell{i}>", "exec")   # raises SyntaxError on failure


def literal_strings(src: str) -> list[str]:
    """RUNTIME values of every string literal in a cell, via the AST.

    Reading the raw source text instead would be wrong twice over: a correct
    non-raw literal ``"$\\\\langle x$"`` appears as two backslashes in the source
    but is ONE at runtime, and an f-string's placeholders are not literal text.
    The AST gives what the cell will actually hold. For f-strings only the
    literal segments are returned (a math span split across a placeholder is not
    checked — acceptable, none are written that way).
    """
    import ast
    out = []
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return out                                     # the compile test owns this
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            out.append(node.value)
        elif isinstance(node, ast.JoinedStr):
            out.append("".join(v.value for v in node.values
                               if isinstance(v, ast.Constant)
                               and isinstance(v.value, str)))
    return out


# Control characters that are NEVER intentional in these cells and are exactly
# what a prematurely-applied escape leaves behind:
#   \a BEL (\approx)  \b BS (\beta)   \f FF (\frac!)  \v VT (\varphi)  \r CR (\rangle, \rm)
# \t and \n ARE legitimately written by these builders, so they are excluded --
# an eaten \n leaves a real newline inside a literal, which is a SyntaxError and
# is caught by test_every_code_cell_compiles instead.
BAD_CONTROL = {chr(c) for c in range(32)} - {"\t", "\n"}
CTRL_NAMES = {"\a": r"\a (from \approx …)", "\b": r"\b (from \beta …)",
              "\f": r"\f (from \frac …)", "\v": r"\v (from \varphi …)",
              "\r": r"\r (from \rangle, \rm …)"}


@pytest.mark.parametrize("label,thunk", BUILDS, ids=IDS)
def test_no_control_characters_in_string_literals(label, thunk):
    """A control char inside a cell literal means an escape already fired.

    This is the CELL-LEVEL signature of the builder bug: a non-raw builder string
    containing LaTeX has its escape applied when the BUILDER is parsed, so the
    emitted cell holds the control character rather than the LaTeX command.
    Precise (no false positives on correctly-escaped ``\\\\langle`` or a genuine
    ``\\n``) and it catches the SILENT variants that still compile.
    """
    offenders = []
    for i, src in code_cells(thunk):
        for s in literal_strings(src):
            bad = sorted(BAD_CONTROL & set(s))
            if bad:
                names = ", ".join(CTRL_NAMES.get(c, repr(c)) for c in bad)
                offenders.append(f"{label} cell{i}: {names} in {s[:70]!r}")
    assert not offenders, (
        "control character(s) inside cell string literals — a LaTeX escape was "
        "applied when the BUILDER was parsed; use code(r\"\"\"...\"\"\"):\n  "
        + "\n  ".join(offenders))


@pytest.mark.parametrize("label,thunk", BUILDS, ids=IDS)
def test_every_math_span_actually_renders(label, thunk):
    """Parse every $...$ span with matplotlib's own parser.

    Stronger than any blocklist: it catches unsupported commands nobody thought
    to enumerate. This is what caught \\frac12 before it reached a GPU node.
    Spans are taken from RUNTIME literal values, not source text.
    """
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib import mathtext

    parser = mathtext.MathTextParser("agg")
    offenders = []
    spans = {(i, s) for i, src in code_cells(thunk)
             for lit in literal_strings(src)
             for s in re.findall(r"\$[^$\n]+\$", lit)}
    for i, span in sorted(spans):
        try:
            parser.parse(span, dpi=72, prop=None)
        except Exception as e:                        # noqa: BLE001
            offenders.append(f"{label} cell{i}: {span[:60]!r} -> {type(e).__name__}: {e}")
    assert not offenders, ("mathtext span(s) matplotlib cannot render:\n  "
                           + "\n  ".join(offenders))


@pytest.mark.parametrize("label,thunk", BUILDS, ids=IDS)
def test_no_unsupported_mathtext_commands(label, thunk):
    offenders = []
    for i, src in code_cells(thunk):
        for cmd in UNSUPPORTED_MATHTEXT:
            for ln, line in enumerate(src.splitlines(), 1):
                if cmd in line:
                    offenders.append(f"{label} cell{i} L{ln}: {cmd.strip()!r}")
    assert not offenders, ("mathtext command(s) matplotlib does not implement:\n  "
                           + "\n  ".join(offenders))


@pytest.mark.parametrize("label,thunk", BUILDS, ids=IDS)
def test_math_spans_are_balanced(label, thunk):
    offenders = []
    for i, src in code_cells(thunk):
        for ln, line in enumerate(src.splitlines(), 1):
            if "$" in line and len(re.findall(r"(?<!\\)\$", line)) % 2:
                offenders.append(f"{label} cell{i} L{ln}: {line.strip()[:80]}")
    assert not offenders, ("odd number of '$' — unterminated math span:\n  "
                           + "\n  ".join(offenders))


# ---------------------------------------------------------------------------
# Content contracts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("half", ["wp", "classical"])
def test_run_notebook_has_the_interaction_energy_section(family, half):
    """.claude/rules/decomposed-interaction-energies.md — every run notebook must
    show the pairwise decomposition AND gate its closure against INQ's scalars."""
    nb = _run_nb(family, half)
    joined = "\n".join(c.source for c in nb.cells)
    assert "K.load_interactions(RUN, HALF)" in joined, "no interactions.csv load"
    assert "Pairwise interaction energies" in joined, "no section-7 heading"
    assert "closure vs INQ E_H" in joined, "closure gate not printed"
    for term in ("e_ss", "e_pp", "e_ps"):
        assert term in joined, f"{term} never plotted"


@pytest.mark.parametrize("family", FAMILIES)
def test_phase_notebook_displays_the_density_gif(family):
    """.claude/rules/notebook-density-gif.md — produced AND displayed inline.

    Writing the file and printing "wrote ...gif" is explicitly NOT sufficient.
    """
    joined = "\n".join(c.source for c in _phase_nb(family).cells)
    assert "make_twin_density_matrix" in joined, "no density-matrix GIF cell"
    assert "display(Image(filename=path))" in joined, \
        "density GIF is produced but never DISPLAYED inline"


def test_builder_compile_guard_is_wired():
    """check_cells_compile() must stay on the phase build path, not just exist."""
    src = (HYP / "bulk_ks_stopping" / "build_phase_notebook.py").read_text()
    assert "check_cells_compile(nb)" in src, \
        "build_one() no longer calls check_cells_compile()"


# ---------------------------------------------------------------------------
# Negative self-tests — the guards must FIRE on the three bugs that shipped.
# A guard that passes on everything is worthless; these pin its sensitivity by
# reconstructing each real 2026-08-01 defect exactly as it appeared in a cell.
# ---------------------------------------------------------------------------

def _fake(src: str):
    """A one-code-cell notebook thunk, for driving the guards directly."""
    import nbformat as nbf
    return lambda: nbf.v4.new_notebook(cells=[nbf.v4.new_code_cell(src)])


def _check_ctrl(label, thunk):
    """Undecorated body of the control-character guard, for the self-tests."""
    offenders = []
    for i, src in code_cells(thunk):
        for s in literal_strings(src):
            if BAD_CONTROL & set(s):
                offenders.append(f"{label} cell{i}: {s[:70]!r}")
    assert not offenders, "control character(s) inside cell string literals"


def _check_math(label, thunk):
    """Undecorated body of the mathtext-render guard, for the self-tests."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib import mathtext
    parser = mathtext.MathTextParser("agg")
    offenders = []
    for i, src in code_cells(thunk):
        for lit in literal_strings(src):
            for span in re.findall(r"\$[^$\n]+\$", lit):
                try:
                    parser.parse(span, dpi=72, prop=None)
                except Exception:                     # noqa: BLE001
                    offenders.append(f"{label} cell{i}: {span!r}")
    assert not offenders, "mathtext span(s) matplotlib cannot render"


def test_guard_fires_on_eaten_carriage_return():
    r"""BUG 1: `\rangle` in a non-raw builder string -> CR in the cell.

    The literal below holds a REAL \r, which is what the builder emitted. Python
    treats a bare CR as a LINE TERMINATOR, so this is an unterminated string
    literal -- caught by the COMPILE guard, like the eaten-\n case, not by the
    control-character guard (which owns the escapes that still compile: \a \b
    \f \v). That is precisely the SyntaxError the 5 broken phase cells showed.
    """
    thunk = _fake('ax.set_ylabel("$\\langle p_z\rangle$")')
    with pytest.raises(SyntaxError):
        for i, src in code_cells(thunk):
            compile(src, f"<cell{i}>", "exec")


def test_guard_fires_on_eaten_bel_from_approx():
    r"""BUG 1, SILENT variant: `\approx` -> BEL. Compiles fine, renders wrong."""
    thunk = _fake('ax.set_title("$x \x07pprox 0$")')
    with pytest.raises(AssertionError):
        _check_ctrl("fake", thunk)


def test_guard_fires_on_eaten_form_feed_from_frac():
    r"""BUG 1, SILENT variant: `\frac` in a non-raw string -> FF + 'rac'."""
    thunk = _fake('ax.set_ylabel("$\x0crac{1}{2}mv^2$")')
    with pytest.raises(AssertionError):
        _check_ctrl("fake", thunk)


def test_guard_fires_on_frac12_mathtext():
    r"""BUG 3: `\frac12` is valid Python AND valid TeX, but not valid mathtext."""
    thunk = _fake(r'ax.set_ylabel(r"$\frac12mv^2$")')
    with pytest.raises(AssertionError):
        _check_math("fake", thunk)


def test_guard_fires_on_unterminated_literal_from_eaten_newline():
    r"""BUG 2: `\n` in a non-raw builder string -> real newline -> SyntaxError.

    This is what failed all 8 run-notebook array tasks.
    """
    thunk = _fake('ax.annotate(f"cloud clips the +z face\nt = 1.0 a.u.")')
    with pytest.raises(SyntaxError):
        for i, src in code_cells(thunk):
            compile(src, f"<cell{i}>", "exec")


def test_guards_pass_correct_cells():
    """The counterpart: correctly-written cells must NOT trip any guard.

    Covers both legitimate patterns the first version of this test wrongly
    flagged — an escaped backslash in a non-raw literal, and a genuine \\n.
    """
    thunk = _fake(
        'print("\\n--- header ---")\n'
        'ax.set_ylabel("$T_1=\\\\langle p^2\\\\rangle/2m$  [eV]")\n'
        'ax.set_title(r"$\\frac{1}{2}mv^2$ and $\\sigma_{WP}$")\n'
    )
    _check_ctrl("fake", thunk)
    _check_math("fake", thunk)
    for i, src in code_cells(thunk):
        compile(src, f"<cell{i}>", "exec")
