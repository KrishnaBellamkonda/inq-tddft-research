"""
plots.py — reusable scientific plotting functions for inqview.

Covers time-series observables (energy, current, dipole) produced by
inqkit::io::ObservablesWriter, and frequency-domain spectra from fourier.py.
Each function accepts the path to an observables.csv (or a FourierResult) and
returns a matplotlib Figure.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

import matplotlib.pyplot as plt
import pandas as pd

from .config import DEFAULT_THEME

if TYPE_CHECKING:
    from .fourier import FourierResult
    from .fields import RealField3D
    from .screens import LeedPattern

PathLike = Union[str, Path]


def load_observables(csv_path: PathLike) -> pd.DataFrame:
    """Load an observables CSV into a DataFrame, stripping whitespace from column names."""
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()
    return df


def plot_total_energy_vs_time(
    csv_path: PathLike,
    ax: "plt.Axes | None" = None,
    **kwargs,
) -> plt.Figure:
    """Plot ONLY the total energy vs time (TODO 1e).

    Use ``plot_all_energy_components_vs_time`` for the multi-component
    version (kinetic, Hartree, XC, total all on the same axes).
    """
    df = load_observables(csv_path)
    defaults = DEFAULT_THEME.plot
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=defaults.figsize, dpi=defaults.dpi)
    else:
        fig = ax.figure
    if "energy_total" not in df.columns:
        raise ValueError(f"{csv_path}: energy_total column missing")
    ax.plot(df["time_au"], df["energy_total"],
            color=defaults.line_colors[0], label="E_total", **kwargs)
    ax.set_xlabel("Time (a.u.)")
    ax.set_ylabel("Total energy (Ha)")
    ax.legend()
    ax.set_title("Total energy vs Time")
    fig.tight_layout()
    return fig


def plot_all_energy_components_vs_time(
    csv_path: PathLike,
    ax: "plt.Axes | None" = None,
    **kwargs,
) -> plt.Figure:
    """Plot every available energy component on a single axes.

    Components rendered when present in the CSV:
    energy_total, energy_kinetic, energy_hartree, energy_xc.
    """
    df = load_observables(csv_path)
    defaults = DEFAULT_THEME.plot

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=defaults.figsize, dpi=defaults.dpi)
    else:
        fig = ax.figure

    time = df["time_au"]
    colors = defaults.line_colors

    if "energy_total" in df.columns:
        ax.plot(time, df["energy_total"], color=colors[0], label="E_total", **kwargs)
    if "energy_kinetic" in df.columns:
        ax.plot(time, df["energy_kinetic"], color=colors[1], label="E_kinetic", **kwargs)
    if "energy_hartree" in df.columns:
        ax.plot(time, df["energy_hartree"], color=colors[2], label="E_hartree", **kwargs)
    if "energy_xc" in df.columns:
        ax.plot(time, df["energy_xc"], color=colors[3], label="E_xc", **kwargs)

    ax.set_xlabel("Time (a.u.)")
    ax.set_ylabel("Energy (Ha)")
    ax.legend()
    ax.set_title("All energy components vs Time")
    fig.tight_layout()
    return fig


# Backwards-compatible alias: existing callers passing plot_energy_vs_time
# get the multi-component figure (the historical behaviour).
plot_energy_vs_time = plot_all_energy_components_vs_time


def plot_current_vs_time(
    csv_path: PathLike,
    ax: "plt.Axes | None" = None,
    **kwargs,
) -> plt.Figure:
    """Plot all available current components vs time."""
    df = load_observables(csv_path)
    defaults = DEFAULT_THEME.plot

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=defaults.figsize, dpi=defaults.dpi)
    else:
        fig = ax.figure

    time = df["time_au"]
    colors = defaults.line_colors
    components = [("current_x", "J_x", colors[0]),
                  ("current_y", "J_y", colors[1]),
                  ("current_z", "J_z", colors[2])]

    for col, label, color in components:
        if col in df.columns:
            ax.plot(time, df[col], color=color, label=label, **kwargs)

    ax.set_xlabel("Time (a.u.)")
    ax.set_ylabel("Current (a.u.)")
    ax.legend()
    ax.set_title("Electronic Current vs Time")
    fig.tight_layout()
    return fig


def plot_dipole_vs_time(
    csv_path: PathLike,
    ax: "plt.Axes | None" = None,
    **kwargs,
) -> plt.Figure:
    """Plot all available dipole components vs time."""
    df = load_observables(csv_path)
    defaults = DEFAULT_THEME.plot

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=defaults.figsize, dpi=defaults.dpi)
    else:
        fig = ax.figure

    time = df["time_au"]
    colors = defaults.line_colors
    components = [("dipole_x", "d_x", colors[0]),
                  ("dipole_y", "d_y", colors[1]),
                  ("dipole_z", "d_z", colors[2])]

    for col, label, color in components:
        if col in df.columns:
            ax.plot(time, df[col], color=color, label=label, **kwargs)

    ax.set_xlabel("Time (a.u.)")
    ax.set_ylabel("Dipole (a.u.)")
    ax.legend()
    ax.set_title("Dipole Moment vs Time")
    fig.tight_layout()
    return fig


def plot_observables_summary(
    csv_path: PathLike,
    output_path: PathLike | None = None,
    **kwargs,
) -> plt.Figure:
    """3-panel summary: energy | current | dipole vs time.

    Panels are omitted if their columns are absent from the CSV.
    """
    df = load_observables(csv_path)
    defaults = DEFAULT_THEME.plot

    has_energy  = any(c in df.columns for c in ("energy_total", "energy_kinetic"))
    has_current = any(c in df.columns for c in ("current_x", "current_y", "current_z"))
    has_dipole  = any(c in df.columns for c in ("dipole_x", "dipole_y", "dipole_z"))

    n_panels = sum([has_energy, has_current, has_dipole])
    if n_panels == 0:
        raise ValueError("No plottable columns found in the CSV.")

    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(defaults.figsize[0], defaults.figsize[1] * n_panels),
        dpi=defaults.dpi,
        sharex=True,
    )
    if n_panels == 1:
        axes = [axes]

    time = df["time_au"]
    colors = defaults.line_colors
    panel = 0

    if has_energy:
        ax = axes[panel]
        if "energy_total" in df.columns:
            ax.plot(time, df["energy_total"], color=colors[0], label="E_total", **kwargs)
        if "energy_kinetic" in df.columns:
            ax.plot(time, df["energy_kinetic"], color=colors[1], label="E_kinetic", **kwargs)
        ax.set_ylabel("Energy (Ha)")
        ax.legend(fontsize="small")
        ax.set_title("Observables Summary")
        panel += 1

    if has_current:
        ax = axes[panel]
        for col, label, color in [("current_x", "J_x", colors[0]),
                                   ("current_y", "J_y", colors[1]),
                                   ("current_z", "J_z", colors[2])]:
            if col in df.columns:
                ax.plot(time, df[col], color=color, label=label, **kwargs)
        ax.set_ylabel("Current (a.u.)")
        ax.legend(fontsize="small")
        panel += 1

    if has_dipole:
        ax = axes[panel]
        for col, label, color in [("dipole_x", "d_x", colors[0]),
                                   ("dipole_y", "d_y", colors[1]),
                                   ("dipole_z", "d_z", colors[2])]:
            if col in df.columns:
                ax.plot(time, df[col], color=color, label=label, **kwargs)
        ax.set_ylabel("Dipole (a.u.)")
        ax.legend(fontsize="small")

    axes[-1].set_xlabel("Time (a.u.)")
    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=defaults.dpi)
        plt.close(fig)
    return fig


# ---------------------------------------------------------------------------
# Frequency-domain plots
# ---------------------------------------------------------------------------

def plot_spectrum(
    result: "FourierResult",
    output_path: PathLike,
    x_max_au: float | None = None,
    log_scale: bool = False,
) -> plt.Figure:
    """Plot the amplitude spectrum from a FourierResult.

    Parameters
    ----------
    result      : FourierResult returned by FourierTransform.transform*().
    output_path : path where the PNG is saved.
    x_max_au    : if set, truncate the x-axis to this frequency (atomic units).
    log_scale   : if True, use a log y-axis.
    """
    defaults = DEFAULT_THEME.plot
    fig, ax = plt.subplots(figsize=defaults.figsize, dpi=defaults.dpi)

    freq = result.frequency_au
    amp = result.amplitude
    if x_max_au is not None:
        mask = freq <= x_max_au
        freq = freq[mask]
        amp = amp[mask]

    ax.plot(freq, amp, color=defaults.line_colors[0], linewidth=1.2)
    ax.set_xlabel("Frequency (Ha/\u0127)")
    ax.set_ylabel("|FFT| (normalised)")
    ax.set_title(f"Spectrum: {result.column}")
    if log_scale:
        ax.set_yscale("log")

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    return fig


def plot_spectrum_summary(
    results: list["FourierResult"],
    output_path: PathLike,
    x_max_au: float | None = None,
    log_scale: bool = False,
) -> plt.Figure:
    """Multi-panel spectrum figure, one subplot per FourierResult, shared x-axis."""
    if not results:
        raise ValueError("results list is empty.")

    defaults = DEFAULT_THEME.plot
    n = len(results)
    fig, axes = plt.subplots(
        n, 1,
        figsize=(defaults.figsize[0], defaults.figsize[1] * n),
        dpi=defaults.dpi,
        sharex=True,
    )
    if n == 1:
        axes = [axes]

    for ax, result in zip(axes, results):
        freq = result.frequency_au
        amp = result.amplitude
        if x_max_au is not None:
            mask = freq <= x_max_au
            freq = freq[mask]
            amp = amp[mask]

        ax.plot(freq, amp, color=defaults.line_colors[0], linewidth=1.2)
        ax.set_ylabel("|FFT| (normalised)")
        ax.set_title(f"Spectrum: {result.column}")
        if log_scale:
            ax.set_yscale("log")

    axes[-1].set_xlabel("Frequency (Ha/\u0127)")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    return fig


# ---------------------------------------------------------------------------
# Field slice plots
# ---------------------------------------------------------------------------

def plot_density_slice(
    field: "RealField3D",
    output_path: PathLike,
    axis: int = 2,
    slice_index: int | None = None,
) -> plt.Figure:
    """Plot a 2D colour map of a density slice through a RealField3D.

    Parameters
    ----------
    field       : RealField3D loaded via inqview.load_real_field().
    output_path : path where the PNG is saved.
    axis        : axis normal to the slice plane: 0=x, 1=y, 2=z.
    slice_index : index along that axis. None uses the midplane.
    """
    import numpy as np

    defaults = DEFAULT_THEME.plot
    data = field.array  # shape (nx, ny, nz), x_slowest_z_fastest

    shape = data.shape
    if slice_index is None:
        slice_index = shape[axis] // 2

    ox, oy, oz = field.meta.origin_bohr
    dx, dy, dz = field.meta.spacing_bohr
    nx, ny, nz = shape

    axis_labels = ["x", "y", "z"]

    if axis == 0:
        slc = data[slice_index, :, :]   # shape (ny, nz)
        xlabel = "z (bohr)"
        ylabel = "y (bohr)"
        extent = [oz, oz + nz * dz, oy, oy + ny * dy]
        pos_bohr = ox + slice_index * dx
    elif axis == 1:
        slc = data[:, slice_index, :]   # shape (nx, nz)
        xlabel = "z (bohr)"
        ylabel = "x (bohr)"
        extent = [oz, oz + nz * dz, ox, ox + nx * dx]
        pos_bohr = oy + slice_index * dy
    else:  # axis == 2
        slc = data[:, :, slice_index]   # shape (nx, ny)
        xlabel = "y (bohr)"
        ylabel = "x (bohr)"
        extent = [oy, oy + ny * dy, ox, ox + nx * dx]
        pos_bohr = oz + slice_index * dz

    fig, ax = plt.subplots(figsize=defaults.figsize, dpi=defaults.dpi)
    im = ax.imshow(
        slc,
        origin="lower",
        extent=extent,
        cmap=defaults.scalar_cmap,
        aspect="equal",
        interpolation="nearest",
    )
    plt.colorbar(im, ax=ax, label="Density (bohr\u207b\u00b3)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    field_label = field.meta.field_name or "density"
    ax.set_title(
        f"{field_label}  |  {axis_labels[axis]} = {pos_bohr:.2f} bohr"
        f"  (slice {slice_index})"
    )
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    return fig


# ---------------------------------------------------------------------------
# LEED screen plots
# ---------------------------------------------------------------------------

def plot_leed_pattern(
    pattern: "LeedPattern",
    output_path: PathLike,
    log_scale: bool = False,
    vmax: float | None = None,
    show_colorbar: bool = True,
) -> plt.Figure:
    """2D colourmap of a LeedPattern with physical axes in bohr.

    Parameters
    ----------
    pattern     : LeedPattern loaded via inqview.load_leed_pattern().
    output_path : path where the PNG is saved.
    log_scale   : use logarithmic colour scale (adds 1 to avoid log(0)).
    vmax        : upper colour limit; None uses data max.
    show_colorbar : whether to show the colour bar.
    """
    import numpy as np

    defaults = DEFAULT_THEME.plot
    data = pattern.data.copy()

    if log_scale:
        data = np.log1p(data)
        clabel = "log(1 + \u03c1\u00b7dt)  [bohr\u207b\u00b3\u00b7a.u.]"
    else:
        clabel = "\u03c1\u00b7dt  [bohr\u207b\u00b3\u00b7a.u.]"

    fig, ax = plt.subplots(figsize=defaults.figsize, dpi=defaults.dpi)
    im = ax.imshow(
        data,
        origin="lower",
        extent=pattern.extent_bohr,
        cmap=defaults.scalar_cmap,
        aspect="equal",
        interpolation="nearest",
        vmin=0.0,
        vmax=vmax,
    )
    if show_colorbar:
        plt.colorbar(im, ax=ax, label=clabel)
    ax.set_xlabel("x (bohr)")
    ax.set_ylabel("y (bohr)")
    ax.set_title(
        f"{pattern.label}  |  z = {pattern.z_bohr:.2f} bohr"
        f"  |  T = {pattern.total_time_au:.2f} a.u."
    )
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    return fig
