"""Data-contract tests for the energy-component renderers (IV-M07 viz half).

Renderers are not pixel-tested (ADR 0005), but rule #6 still demands a test: we
assert the rendered artists carry EXACTLY the numbers off the
``EnergyComponents`` dataclass — bar heights == ``breakdown(...)``, line y-data
== ``dE_* * HA_TO_EV``. That catches the real failure mode (a renderer that
recomputes or mis-wires the data) without comparing images, and proves the
compute→render contract (render consumes, never recomputes).

Expected values are fixed up front from a tiny synthetic observables frame with
an analytically known decomposition — NOT captured from the renderer's output.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless, portable

import numpy as np
import pandas as pd
import pytest

from inqview.analysis.energy_components import HA_TO_EV, compute
from inqview.visualisation import energy_components as ecr

pytestmark = pytest.mark.theme

# Tiny known system. E_ext = E_total - (kin + H + xc), so:
#   E_ext = [13-12, 15-13, 16-14] = [1, 2, 2]
_DF = pd.DataFrame(
    {
        "time_au": [0.0, 1.0, 2.0],
        "energy_total": [13.0, 15.0, 16.0],
        "energy_kinetic": [10.0, 11.0, 13.0],
        "energy_hartree": [5.0, 5.0, 4.0],
        "energy_xc": [-3.0, -3.0, -3.0],
    }
)
_COMPONENTS = ("kinetic", "hartree", "xc", "external")


@pytest.fixture()
def ec():
    return compute(_DF)


def test_kernel_decomposition_known(ec):
    # Guard the inputs the renderer tests rely on (sum invariant + residual).
    assert np.allclose(ec.E_ext, [1.0, 2.0, 2.0])
    assert np.allclose(ec.component_sum(), ec.E_total)


def test_bars_heights_equal_breakdown(ec):
    fig, ax = ecr.render_initial_vs_final_bars(ec)
    init, final = ec.breakdown("initial"), ec.breakdown("final")
    # Two BarContainers: [0] = initial group, [1] = final group, in _COMPONENTS order.
    assert len(ax.containers) == 2
    init_heights = [p.get_height() for p in ax.containers[0]]
    final_heights = [p.get_height() for p in ax.containers[1]]
    assert np.allclose(init_heights, [init[c] for c in _COMPONENTS])
    assert np.allclose(final_heights, [final[c] for c in _COMPONENTS])
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_flow_lines_ydata_equal_dE_in_eV(ec):
    fig, ax = ecr.render_flow_lines(ec)
    by_label = {ln.get_label(): ln for ln in ax.get_lines()}
    expected = {
        "kinetic": ec.dE_kin,
        "Hartree": ec.dE_hartree,
        "xc": ec.dE_xc,
        "external": ec.dE_ext,
        "total": ec.dE_total,
    }
    for label, dE in expected.items():
        assert label in by_label, f"missing line {label!r}"
        assert np.allclose(by_label[label].get_ydata(), dE * HA_TO_EV)
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_breakdown_gif_writes_file(ec, tmp_path):
    out = tmp_path / "energy_flow.gif"
    try:
        result = ecr.render_breakdown_gif(ec, out, fps=4)
    except RuntimeError as exc:  # no Pillow writer in this environment
        pytest.skip(f"GIF writer unavailable: {exc}")
    assert result == out
    assert out.exists() and out.stat().st_size > 0
