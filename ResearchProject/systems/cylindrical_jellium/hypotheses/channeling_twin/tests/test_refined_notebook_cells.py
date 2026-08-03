"""Static guards on the cells build_refined_notebook.py emits.

These run in milliseconds and catch the failure modes that otherwise only show
up after a multi-minute notebook execution — or, worse, do not show up at all
and silently corrupt a label.

THE BUG CLASS THIS EXISTS FOR. Notebook cells are written as Python strings
inside the builder, so a LaTeX label passes through TWO levels of string
escaping. ``\\frac`` inside a NON-raw builder string becomes a FORM FEED plus
"rac", and ``\\rangle`` becomes a CARRIAGE RETURN plus "angle" — silently, with
no error, producing a mangled axis label in a figure nobody re-reads. Both were
present on first write of build_refined_notebook.py and were caught by the
scanner below, not by eye.

Matplotlib mathtext is also a SUBSET of LaTeX: ``\\le`` is not a command (``\\leq``
is), and an unknown command raises only when the figure is rendered.
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

BUILDER = PKG / "build_refined_notebook.py"


@pytest.fixture(scope="module")
def emitted():
    import build_refined_notebook as B
    return B.cells("wp", "classical", [0.0, 15.0, 30.0])


def test_every_code_cell_parses(emitted):
    for i, (kind, src) in enumerate(emitted):
        if kind == "code":
            try:
                ast.parse(src)
            except SyntaxError as exc:
                pytest.fail(f"cell {i} is not valid Python: {exc}\n{src[:400]}")


def test_no_control_characters_leaked_into_any_cell(emitted):
    """No FORM FEED / CARRIAGE RETURN / VERTICAL TAB anywhere.

    Their presence means a LaTeX command was written into a non-raw builder
    string and Python ate the backslash: \\frac -> \\x0c + 'rac',
    \\rangle -> \\r + 'angle', \\vec -> \\x0b + 'ec'.
    """
    for i, (kind, src) in enumerate(emitted):
        for ch, name in (("\f", "FORM FEED (\\frac?)"),
                         ("\r", "CARRIAGE RETURN (\\rangle?)"),
                         ("\v", "VERTICAL TAB (\\vec?)"),
                         ("\a", "BEL (\\approx?)"),
                         ("\b", "BACKSPACE (\\beta?)")):
            assert ch not in src, (
                f"cell {i} ({kind}) contains a raw {name}. A LaTeX command was "
                f"written into a NON-raw builder string — make that code(...) or "
                f"md(...) argument a raw string, or double the backslash.")


def test_builder_source_has_no_swallowable_escapes():
    """Scan the BUILDER too, so the bug is caught at its source.

    A single backslash followed by one of abfnrtv starts a valid Python escape.
    Inside a raw string that is fine; the emitted-cell test above is what proves
    it. This test additionally flags any such sequence so a reviewer sees it.
    """
    src = BUILDER.read_text()
    hits = sorted({m.group(0) for m in
                   re.finditer(r'(?<!\\)\\([abfnrtv])[a-zA-Z]{2,}', src)})
    # Whatever survives here must be inside raw strings — proven by the
    # control-character test. Recorded so the list cannot grow unnoticed.
    assert set(hits) <= {"\\frac", "\\rangle"}, (
        f"new swallowable escape sequences in the builder: {hits}. "
        f"Confirm they sit inside raw strings, then add them here.")


def test_no_unsupported_mathtext_commands(emitted):
    """Matplotlib mathtext is a LaTeX subset; these are the usual mistakes."""
    forbidden = (r"\le ", r"\le}", r"\le$", r"\ge ", r"\ge}", r"\ge$",
                 r"\text{", r"\begin{", r"\eqnarray", r"\substack")
    for i, (kind, src) in enumerate(emitted):
        if kind != "code":
            continue          # markdown is rendered by MathJax, which allows more
        for f in forbidden:
            assert f not in src, (
                f"cell {i} uses {f!r}, which matplotlib mathtext does not "
                f"support (use \\leq / \\geq / \\mathrm{{}} instead)")


def test_dollar_signs_are_balanced_in_code_cells(emitted):
    for i, (kind, src) in enumerate(emitted):
        if kind != "code":
            continue
        n = len(re.findall(r"(?<!\\)\$", src))
        assert n % 2 == 0, f"cell {i} has an odd number ({n}) of '$' delimiters"


def test_required_sections_are_present(emitted):
    """The user's requested section list, in order. Pins scope."""
    md_text = "\n".join(src for kind, src in emitted if kind == "md")
    for needle in ("## 1. Where is the projectile?",
                   "## 2. Kinetic energy",
                   "### 2a. Classical",
                   "### 2b. Wavepacket",
                   "### 2c. Classical against wavepacket",
                   "## 3. The momentum distribution, at three times",
                   "## 3b. The 2-D momentum map",
                   "## 4. Pairwise interaction energies",
                   "## 4b. $\\Delta(E_{PS} + E_{PB})$",
                   "## 5. Can the $\\Delta T_1$ vs classical difference be explained?",
                   "## 6. Choose the fit window"):
        assert needle in md_text, f"missing section: {needle}"


def test_the_t1_t2_convention_table_is_stated(emitted):
    """The label swap MUST be documented in the notebook itself.

    refined.py uses the user's convention (T1 = drift), which is the reverse of
    ks_stopping.py's. A reader who assumes the engine convention reads the
    study's conclusion backwards, so the table is not optional decoration.
    """
    md_text = "\n".join(src for kind, src in emitted if kind == "md")
    assert "swapped" in md_text.lower()
    assert "ks_stopping.py" in md_text


def test_windows_are_the_user_chosen_ones(emitted):
    """The three USER-CHOSEN windows (2026-08-02), pinned verbatim.

    These are not defaults and must not be silently 'improved': the user picked
    9-25 for T1 off the local-slope plateau, and asked for BOTH 21-30 and 5-20
    for T2 explicitly to compare them. Changing them changes the reported
    stopping powers, so the change has to be deliberate enough to update a test.
    """
    code_text = "\n".join(src for kind, src in emitted if kind == "code")
    assert '("T1  9-25",        "T1", ( 9.0, 25.0))' in code_text
    assert '("T2  21-30",       "T2", (21.0, 30.0))' in code_text
    assert '("T2  5-20",        "T2", ( 5.0, 20.0))' in code_text


def test_classical_is_fitted_over_the_same_window_as_each_wp_row(emitted):
    """A WP number is only comparable to a classical one from the SAME window.

    Both projectiles decelerate, so cross-window comparison compares different
    velocities. Guards against the classical fit drifting back to a single
    global window.
    """
    code_text = "\n".join(src for kind, src in emitted if kind == "code")
    assert "cl.t.to_numpy(), t0, t1" in code_text
    assert "S_cl" in code_text and "ratio" in code_text
