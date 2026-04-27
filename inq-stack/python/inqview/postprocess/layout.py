"""Phase: ``layout`` — per-run xz layout diagram.

Produces ``results/analysis/layout/layout_xz.png`` summarising the
geometric setup of the run: cell extent, target plane, WP starting
position, and every screen z-position with its index annotated.

The diagram is purely a 2D xz projection. Inputs are read from
``results/run_summary.txt`` (cell, WP centre, σ) and either
``results/raw/screens/screen_config.csv`` (preferred) or — when the
config file is missing because of the early-Phase-1 ofstream-parent
bug — ``run_summary.txt``'s ``screen_z[k]`` lines.

Generalisable: takes ``(Lx, Lz, b, sigma, screen_z_list,
target_extent_x_bohr)``; the coronene wrapper supplies the molecule's
in-plane extent (from the canonical ``shared/geometry/coronene.xyz``).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

from . import _common
from . import pipeline as _pipeline


# Approximate in-plane radius of coronene from
# ResearchProject/systems/coronene/shared/geometry/coronene.xyz
# (max |x| = 4.578674 Å = 8.65 Bohr; max |y| = 4.758297 Å = 8.99 Bohr).
_CORONENE_HALF_EXTENT_BOHR = 8.65


def _parse_run_summary(path: Path) -> dict:
    """Parse the simple ``key = value`` layout of run_summary.txt.

    Returns the keys we need:
      cell_bohr, wp_center_bohr, wp_sigma_bohr, screen_z (list).
    """
    out: dict = {"screen_z": []}
    if not path.exists():
        return out
    z_re = re.compile(r"^screen_z\[(\d+)\]\s*=\s*(\S+)")
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = z_re.match(line)
            if m:
                out["screen_z"].append(float(m.group(2)))
                continue
            if "=" in line:
                k, _, v = line.partition("=")
                out[k.strip()] = v.strip()
    return out


def _parse_screen_config_csv(path: Path) -> list[float]:
    """Return [z_bohr_screen_0, ...] from screen_config.csv if it exists."""
    if not path.exists():
        return []
    z: list[float] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n")
        cols = [c.strip() for c in header.split(",")]
        try:
            iz = cols.index("z_bohr")
        except ValueError:
            return []
        for line in fh:
            parts = [c.strip() for c in line.rstrip("\n").split(",")]
            if len(parts) > iz:
                try:
                    z.append(float(parts[iz]))
                except ValueError:
                    pass
    return z


def render_layout_xz(
    out_path: Path,
    *,
    run_name: str,
    Lx_bohr: float,
    Lz_bohr: float,
    wp_cz_bohr: float,
    wp_sigma_bohr: float,
    screen_z_bohr: Iterable[float],
    target_half_extent_x_bohr: float = _CORONENE_HALF_EXTENT_BOHR,
    target_z_bohr: float = 0.0,
) -> None:
    """Render a single xz-layout PNG to ``out_path``.

    The diagram contains:

    * Cell rectangle: x ∈ [-Lx/2, +Lx/2], z ∈ [-Lz/2, +Lz/2].
    * Target line: thick horizontal at z = ``target_z_bohr`` spanning
      x ∈ [-target_half_extent_x_bohr, +target_half_extent_x_bohr].
    * WP marker: a horizontal Gaussian-like band at z = ``wp_cz_bohr``,
      vertical half-width = ``wp_sigma_bohr``.
    * 20 horizontal lines, one per screen, with the screen index
      labelled to the right of the cell.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    out_path.parent.mkdir(parents=True, exist_ok=True)

    half_x = 0.5 * Lx_bohr
    half_z = 0.5 * Lz_bohr
    fig, ax = plt.subplots(figsize=(7, 7), dpi=120)

    # Cell rectangle.
    ax.add_patch(Rectangle(
        (-half_x, -half_z), Lx_bohr, Lz_bohr,
        edgecolor="black", facecolor="none", linewidth=1.5))

    # Target plane (the molecule's in-plane footprint at z=target_z).
    ax.plot([-target_half_extent_x_bohr, +target_half_extent_x_bohr],
            [target_z_bohr, target_z_bohr],
            color="forestgreen", linewidth=4, solid_capstyle="round",
            label=f"target (z = {target_z_bohr:g})")

    # WP starting position: centroid line + ±σ band.
    ax.fill_between(
        [-half_x, +half_x],
        wp_cz_bohr - wp_sigma_bohr, wp_cz_bohr + wp_sigma_bohr,
        color="orange", alpha=0.20)
    ax.plot([-half_x, +half_x], [wp_cz_bohr, wp_cz_bohr],
            color="darkorange", linewidth=2,
            label=f"WP centroid (z = {wp_cz_bohr:g}, σ = {wp_sigma_bohr:g})")

    # Screens — one horizontal line per screen, with the index annotated.
    screen_z_list = list(screen_z_bohr)
    for k, z in enumerate(screen_z_list):
        ax.plot([-half_x, +half_x], [z, z],
                color="steelblue", linewidth=0.6, alpha=0.85)
        ax.annotate(
            f"screen_{k:02d}", xy=(half_x, z), xytext=(4, 0),
            textcoords="offset points",
            va="center", ha="left", fontsize=6, color="steelblue")

    # Bring the WP and target above the screens visually.
    ax.set_xlim(-half_x * 1.05, half_x * 1.40)  # right pad for screen labels
    ax.set_ylim(-half_z * 1.05, half_z * 1.05)
    ax.set_xlabel("x (bohr)")
    ax.set_ylabel("z (bohr)")
    ax.set_aspect("equal")
    ax.set_title(_common.title(
        run_name,
        f"layout (xz):  {Lx_bohr:g} × {Lz_bohr:g} bohr cell, "
        f"{len(screen_z_list)} screens",
        multiline=False))
    ax.legend(loc="upper left", fontsize="small", framealpha=0.85)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    summary_path = results_dir / "run_summary.txt"
    if not summary_path.exists():
        _pipeline.skip(f"run_summary.txt missing at {summary_path}")

    summary = _parse_run_summary(summary_path)
    cell = summary.get("cell_bohr", "")
    parts = cell.split()
    if len(parts) < 3:
        _pipeline.skip(f"could not parse cell_bohr from {summary_path!s}")
    Lx_bohr = float(parts[0])
    Lz_bohr = float(parts[2])

    wp_cx, wp_cy, wp_cz = (
        summary.get("wp_center_bohr", "0 0 0").split()[:3])
    wp_cz_bohr = float(wp_cz)
    wp_sigma_bohr = float(summary.get("wp_sigma_bohr", "1.0"))

    # Prefer screen_config.csv (richer, includes labels). Fall back to
    # run_summary.txt's parsed screen_z[k] entries.
    z_list = _parse_screen_config_csv(
        results_dir / "raw" / "screens" / "screen_config.csv")
    if not z_list:
        z_list = list(summary.get("screen_z", []))

    out_dir = _common.ensure_dir(results_dir / "analysis" / "layout")
    out_path = out_dir / "layout_xz.png"
    if not _common.need_rebuild(out_path, rebuild):
        return {"layout_xz": str(out_path), "cached": True,
                "n_screens": len(z_list)}

    render_layout_xz(
        out_path, run_name=run_name,
        Lx_bohr=Lx_bohr, Lz_bohr=Lz_bohr,
        wp_cz_bohr=wp_cz_bohr, wp_sigma_bohr=wp_sigma_bohr,
        screen_z_bohr=z_list,
    )
    return {"layout_xz": str(out_path), "n_screens": len(z_list)}
