from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class PlotDefaults:
    """
    Global plotting defaults for the Python side of inqview.

    These are intentionally semantic rather than ParaView-specific.
    Later, the same semantic roles can be mapped to Matplotlib and ParaView.
    """

    dpi: int = 160
    figsize: tuple[float, float] = (6.4, 4.2)

    # Semantic colormap roles
    scalar_cmap: str = "cividis"  # densities, amplitudes, positive scalar fields
    signed_cmap: str = "coolwarm"  # real/imag orbital parts, signed differences
    phase_cmap: str = "twilight_shifted"  # phase / cyclic data

    # General plot appearance
    facecolor: str = "white"
    axes_facecolor: str = "white"
    savefig_facecolor: str = "white"
    grid: bool = False
    transparent_background: bool = False

    # Sensible default cycle for line plots
    line_colors: tuple[str, ...] = (
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
    )

    # Missing / invalid values
    nan_color: str = "#bdbdbd"


@dataclass(frozen=True)
class RenderDefaults:
    """
    Rendering defaults for future VTI / ParaView pipelines.

    These remain semantic for now. We do not hardwire ParaView transfer-function
    names at this stage.
    """

    scalar_cmap: str = "cividis"
    signed_cmap: str = "coolwarm"
    phase_cmap: str = "twilight_shifted"

    background_rgb: tuple[float, float, float] = (1.0, 1.0, 1.0)
    image_size: tuple[int, int] = (1600, 1200)

    # Future rendering policy knobs
    use_log_scale_for_density: bool = False
    rescale_to_all_timesteps: bool = True
    invert_opacity_for_density: bool = False


@dataclass(frozen=True)
class Theme:
    """
    Top-level theme object.
    """

    name: str = "scientific_cividis"
    plot: PlotDefaults = field(default_factory=PlotDefaults)
    render: RenderDefaults = field(default_factory=RenderDefaults)


DEFAULT_THEME = Theme()
