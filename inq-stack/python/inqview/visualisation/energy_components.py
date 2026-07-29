"""Renderers for the functional energy-component flow (IV-M07 viz half).

Consumes an :class:`inqview.analysis.energy_components.EnergyComponents` result
(the compute half) and answers the question "where did the energy go?" — into
kinetic, Hartree, exchange-correlation, or the external/electron-ion term. Three
views, per the design ask:

- :func:`render_initial_vs_final_bars` — grouped bars, each component at t0 vs
  the final step (absolute energies, Ha).
- :func:`render_flow_lines` — ΔE(t) time series per component (eV), the
  "energy flow": positive = that store gained energy.
- :func:`render_breakdown_gif` — animated ΔE bars over time, saved as a GIF.

Renderers NEVER recompute: every plotted number comes straight off the dataclass
(``breakdown``/``dE_*``/``redistribution_ev``). This is the visualisation layer,
so matplotlib is imported here (never in ``analysis``). Theme: ADR-0004
(``inqview.visualisation.style``).
"""
from __future__ import annotations

from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt

from ..analysis.energy_components import HA_TO_EV, EnergyComponents
from . import style

# Stable categorical colours, shared across all three renderers so a component
# keeps its colour between the bar chart, the line plot and the GIF.
COMPONENT_COLORS: dict[str, str] = {
    "kinetic": "#1b9e77",    # teal
    "hartree": "#d95f02",    # orange
    "xc": "#7570b3",         # violet
    "external": "#666666",   # grey
    "total": "#000000",      # black (sum / reference)
}
# Component order used on every axis (external last; total handled separately).
_COMPONENTS = ("kinetic", "hartree", "xc", "external")
_LABEL = {"kinetic": "kinetic", "hartree": "Hartree", "xc": "xc",
          "external": "external", "total": "total"}


def render_initial_vs_final_bars(ec: EnergyComponents):
    """Grouped bar chart: each energy component at t0 vs the final step (Ha).

    Returns ``(fig, ax)``. Bar heights are exactly ``ec.breakdown('initial'|'final')``.
    """
    style.apply_theme()
    init = ec.breakdown("initial")
    final = ec.breakdown("final")

    fig, ax = style.figure_one_col()
    x = range(len(_COMPONENTS))
    w = 0.38
    ax.bar([i - w / 2 for i in x], [init[c] for c in _COMPONENTS], width=w,
           label="initial", color=[COMPONENT_COLORS[c] for c in _COMPONENTS],
           alpha=0.55, edgecolor="black", linewidth=0.5)
    ax.bar([i + w / 2 for i in x], [final[c] for c in _COMPONENTS], width=w,
           label="final", color=[COMPONENT_COLORS[c] for c in _COMPONENTS],
           alpha=1.0, edgecolor="black", linewidth=0.5)
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xticks(list(x))
    ax.set_xticklabels([_LABEL[c] for c in _COMPONENTS])
    ax.set_ylabel("energy (Ha)")
    ax.set_title("Energy components: initial vs final")
    ax.legend(frameon=False, loc="best")
    return fig, ax


def render_flow_lines(ec: EnergyComponents):
    """ΔE(t) per component in eV (change from t0 — the energy flow).

    Returns ``(fig, ax)``. Each line's y-data is exactly ``ec.dE_<c> * HA_TO_EV``;
    a dashed reference line is the total drift ``ec.dE_total``.
    """
    style.apply_theme()
    t = ec.time_au
    series = {
        "kinetic": ec.dE_kin,
        "hartree": ec.dE_hartree,
        "xc": ec.dE_xc,
        "external": ec.dE_ext,
    }
    fig, ax = style.figure_one_col()
    for c in _COMPONENTS:
        ax.plot(t, series[c] * HA_TO_EV, color=COMPONENT_COLORS[c], label=_LABEL[c])
    ax.plot(t, ec.dE_total * HA_TO_EV, color=COMPONENT_COLORS["total"],
            linestyle="--", linewidth=1.0, label="total")
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\Delta E$ (eV)")
    ax.set_title("Energy flow (change from $t_0$)")
    ax.legend(frameon=False, ncol=2, loc="best")
    return fig, ax


def render_breakdown_gif(
    ec: EnergyComponents,
    path: Union[str, Path],
    *,
    fps: int = 12,
    stride: int = 1,
) -> Path:
    """Animate the per-component ΔE bars over time and save a GIF.

    Each frame shows ``[dE_kin, dE_hartree, dE_xc, dE_ext](t)`` in eV with a fixed
    y-range, so bars growing/shrinking show where energy accumulates. ``stride``
    subsamples frames. Returns the written ``Path``.

    Raises ``RuntimeError`` if no GIF writer (Pillow) is available.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    style.apply_theme()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    dE = {
        "kinetic": ec.dE_kin * HA_TO_EV,
        "hartree": ec.dE_hartree * HA_TO_EV,
        "xc": ec.dE_xc * HA_TO_EV,
        "external": ec.dE_ext * HA_TO_EV,
    }
    frames = list(range(0, len(ec.time_au), max(1, stride)))
    # Fixed y-range across the whole trajectory so frames are comparable.
    lo = min(float(arr.min()) for arr in dE.values())
    hi = max(float(arr.max()) for arr in dE.values())
    pad = 0.05 * max(hi - lo, 1e-12)
    colors = [COMPONENT_COLORS[c] for c in _COMPONENTS]

    fig, ax = style.figure_one_col()
    x = range(len(_COMPONENTS))
    bars = ax.bar(list(x), [0.0] * len(_COMPONENTS), color=colors,
                  edgecolor="black", linewidth=0.5)
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xticks(list(x))
    ax.set_xticklabels([_LABEL[c] for c in _COMPONENTS])
    ax.set_ylabel(r"$\Delta E$ (eV)")
    ax.set_ylim(lo - pad, hi + pad)
    title = ax.set_title("")

    def _update(fi: int):
        for bar, c in zip(bars, _COMPONENTS):
            bar.set_height(dE[c][fi])
        title.set_text(f"Energy flow  t = {ec.time_au[fi]:.2f} a.u.")
        return (*bars, title)

    anim = FuncAnimation(fig, _update, frames=frames, blit=False)
    try:
        anim.save(str(path), writer=PillowWriter(fps=fps))
    except Exception as exc:  # pragma: no cover - environment dependent
        plt.close(fig)
        raise RuntimeError(f"could not write GIF {path}: {exc}") from exc
    plt.close(fig)
    return path
