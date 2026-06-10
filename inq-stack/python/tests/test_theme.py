"""Numeric theme-config test (ADR 0004 / IV-M10).

Guards the DESIGNED constants without rendering any figure: cmap-role mapping,
the fixed-dimension geometry, and the rcParams installed by apply_theme(). No
pixel comparison — these are data assertions on the visual standard.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")          # headless; no display needed

import matplotlib.pyplot as plt
import pytest

from inqview.visualisation import style

pytestmark = pytest.mark.theme


def test_cmap_roles_are_the_designed_values():
    assert style.cmap_for("sequential") == "inferno"
    assert style.cmap_for("diverging") == "RdBu_r"      # zero-centred
    assert style.cmap_for("phase") == "twilight_shifted"


def test_unknown_cmap_role_raises():
    with pytest.raises(ValueError):
        style.cmap_for("rainbow")


def test_fixed_dimension_constants():
    assert style.ONE_COL_IN == (3.5, 3.0)
    assert style.TWO_COL_W_IN == 7.0


def test_figure_one_col_has_exact_size_and_fixed_axes():
    fig, ax = style.figure_one_col()
    try:
        assert tuple(fig.get_size_inches()) == pytest.approx((3.5, 3.0))
        # the axes box is the fixed rectangle (panels align) — not auto-laid-out
        assert ax.get_position().bounds == pytest.approx((0.180, 0.160, 0.785, 0.805))
    finally:
        plt.close(fig)


def test_figure_two_col_is_seven_inches_wide():
    fig, ax = style.figure_two_col()
    try:
        assert fig.get_size_inches()[0] == pytest.approx(7.0)
    finally:
        plt.close(fig)


def test_apply_theme_installs_designed_rcparams():
    style.apply_theme()
    rc = matplotlib.rcParams
    assert rc["font.size"] == 10
    assert rc["axes.linewidth"] == pytest.approx(0.8)
    assert rc["xtick.direction"] == "in"
    assert rc["image.cmap"] == "inferno"


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
