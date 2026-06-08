"""Fixed-palette GIF writer — guarantees an IDENTICAL value->colour mapping in
every frame.

Why this exists
---------------
matplotlib's PillowWriter lets pillow pick a SEPARATE optimised 256-colour
palette per frame. With a fixed colorbar (clim frozen) the underlying
value->RGBA mapping is constant, but the GIF quantiser then re-bins each frame's
RGBA to a DIFFERENT 256-colour table, so the same Δn value can render as a
slightly different colour from frame to frame ("the gradation changes").

`save_gif_fixed_palette` instead:
  1. renders every frame to an RGB array off the Agg canvas,
  2. builds ONE master 256-colour palette from ALL frames at once,
  3. quantises every frame to that SAME palette (no dithering),
so the colour gradation is bit-for-bit identical across the whole animation.
"""
from __future__ import annotations
import numpy as np
from PIL import Image


def _frame_rgb(fig):
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()


def save_gif_fixed_palette(fig, update, n_frames, path, duration_ms=140, dpi=None):
    """Render n_frames via update(k) and write a GIF with a single global
    palette shared by every frame (constant colour gradation).

    update(k) must mutate the figure's artists for frame k (return value
    ignored). dpi, if given, sets the render resolution.
    """
    if dpi is not None:
        fig.set_dpi(dpi)
    rgb = []
    for k in range(n_frames):
        update(k)
        rgb.append(_frame_rgb(fig))
    # one master palette from the WHOLE animation (stack all frames, quantize once)
    stack = np.concatenate(rgb, axis=0)
    master = Image.fromarray(stack).quantize(colors=256, method=Image.MEDIANCUT)
    pil = [Image.fromarray(f).quantize(palette=master, dither=Image.NONE) for f in rgb]
    pil[0].save(path, save_all=True, append_images=pil[1:],
                duration=duration_ms, loop=0, disposal=2, optimize=False)
    return str(path)
