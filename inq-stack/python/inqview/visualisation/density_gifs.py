"""xz-slice density GIF battery for jellium-slab projectile runs.

Produces the matrix  {density, delta0, dstep} × {total, wp, bath}:
  - density : n(x,z,t)
  - delta0  : Δn = n(t) − n(0)
  - dstep   : Δn = n(t+dt) − n(t)
over the categories
  - total : the full electron density (density_total VTIs)
  - wp    : the wavepacket orbital |ψ_WP|² (density_wp VTIs)
  - bath  : n_total − n_wp, the slab response with the projectile removed
            (the canonical bath density, see CONTEXT.md → "Density decomposition")

A WP run (density_wp present) yields up to 9 GIFs; a classical run has no WP
orbital, so it yields 3 GIFs of the **total** density, with delta0 read as the
projectile-**induced wake** (n(t) − n(0)).

VTIs are loaded in PHYSICAL order via ``inqview.load_vti`` — never fftshift'd
(see the vti-coordinate-mapping rule). Each GIF marks the slab faces and CAP
inner faces.

Colour scheme (shared-colorbar rule — every GIF shows LINEAR | LOG side by side):
  - density : linear | log (LogNorm), viridis. ``total`` and ``bath`` share ONE
              range tuned to the slab density (so low densities are visible — the
              WP blob is allowed to saturate); ``wp`` gets its own range. Pass
              ``density_vmax`` to share a single scale across runs (e.g. apply the
              classical-total scale to the WP-total GIF, as requested).
  - delta0 / dstep : symmetric diverging (RdBu_r), linear | symlog (SymLogNorm),
              each owns a robust |Δn| range. The symlog panel exposes the
              low-|Δn| wake tail the linear panel saturates near the projectile.
"""
from __future__ import annotations

import glob
import os
import re
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm, SymLogNorm

from .field_io import load_vti

_TPAT = re.compile(r"_t(\d+)\.vti$")


def _frame_time(path: str, dt: float) -> float:
    m = _TPAT.search(os.path.basename(path))
    return (int(m.group(1)) * dt) if m else float("nan")


def _slice_stack(vti_dir: str, idx, dt: float):
    """Return (times[T], slices[T, nz, nx]) of the xz mid-y plane for frames idx."""
    files = sorted(glob.glob(os.path.join(vti_dir, "*.vti")))
    if not files:
        return None, None, None
    files = [files[k] for k in idx]
    first = load_vti(files[0])
    iy = first.data.shape[1] // 2
    times = np.array([_frame_time(f, dt) for f in files])
    sl = np.empty((len(files), first.data.shape[2], first.data.shape[0]))  # (T,nz,nx)
    for t, f in enumerate(files):
        sl[t] = load_vti(f).data[:, iy, :].T
    return times, sl, (first.x, first.z)


def _save_gif(slices, times, axes, out_path, *, title, cap_lines, kind,
              vmax=None, vmin=None, fps=10, per_frame_norm=False):
    """Render one density-battery GIF as LINEAR | LOG side-by-side panels.

    Every kind carries BOTH scales (shared-colorbar rule "always linear AND log"):
      - ``density`` : linear (left) | log ``LogNorm`` (right), viridis. A shared
        ``vmax`` locks the scale across frames and (for total/bath) across runs.
      - diverging ``delta0``/``dstep`` : linear (left) | ``SymLogNorm`` (right),
        symmetric RdBu_r. The log panel exposes low-|Δn| structure (the wake tail)
        that the linear panel saturates away near the projectile.
    Scales are LOCKED from the WHOLE stack (not the middle frame) so the colorbars
    are fixed and representative of every frame (scientific-figures GIF rule).

    ``per_frame_norm`` (density kind only): normalise the LINEAR panel to each
    frame's own max — n/nₘₐₓ(t) on a fixed 0..1 scale — so a strongly-DISPERSING
    feature (a free wavepacket whose peak collapses ~1/σ³) stays visible as it
    moves, instead of vanishing under a single stack-wide vmax. The LOG panel keeps
    an ABSOLUTE, widened (~4-decade) scale so the physical fade (e.g. CAP
    absorption) remains legible. Off by default → slab-run comparability unchanged."""
    x, z = axes
    ext = [x[0], x[-1], z[0], z[-1]]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.4))
    if kind == "density":
        gmax = float(np.nanmax(slices)) or 1e-12
        vmax = vmax if vmax is not None else float(np.percentile(slices, 99.5))
        if per_frame_norm:
            # linear: per-frame normalised (0..1); log: absolute, widened to span
            # the whole collapse (global max → 4 decades below).
            vmax_lin, vmin_log = 1.0, gmax * 1e-4
            def _lin(fr):
                m = float(np.nanmax(fr)) or 1e-12
                return np.clip(fr, 0.0, None) / m
        else:
            vmax_lin, vmin_log = vmax, (vmin if vmin is not None else vmax * 1e-3)
            def _lin(fr):
                return np.clip(fr, 0.0, None)
        imL = axL.imshow(_lin(slices[0]), origin="lower", aspect="auto",
                         extent=ext, cmap="viridis", vmin=0.0, vmax=vmax_lin)
        imR = axR.imshow(np.clip(slices[0], vmin_log, None), origin="lower", aspect="auto",
                         extent=ext, cmap="viridis", norm=LogNorm(vmin=vmin_log, vmax=gmax))
        cbl = "n / nₘₐₓ(t)  (per-frame)" if per_frame_norm else "n (a₀⁻³, linear)"
        cbr, clip_lo = "n (a₀⁻³, log)", vmin_log
    else:  # diverging difference
        vmax = vmax if vmax is not None else float(np.percentile(np.abs(slices), 99.0)) or 1e-12
        lthr = vmax / 100.0
        imL = axL.imshow(slices[0], origin="lower", aspect="auto",
                         extent=ext, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        imR = axR.imshow(slices[0], origin="lower", aspect="auto", extent=ext,
                         cmap="RdBu_r", norm=SymLogNorm(linthresh=lthr, vmin=-vmax, vmax=vmax, base=10))
        cbl, cbr, clip_lo = "Δn (a₀⁻³, linear)", "Δn (a₀⁻³, symlog)", None
    for ax in (axL, axR):
        for zz in cap_lines:
            ax.axhline(zz, ls="--", lw=0.7, color="0.4")
        ax.set_xlabel("x (Bohr)")
    axL.set_ylabel("z (Bohr)")
    fig.colorbar(imL, ax=axL, fraction=0.046, pad=0.02).set_label(cbl, fontsize=8)
    fig.colorbar(imR, ax=axR, fraction=0.046, pad=0.02).set_label(cbr, fontsize=8)
    sup = fig.suptitle("", fontsize=9)

    def upd(k):
        if kind == "density":
            imL.set_data(_lin(slices[k]))
            imR.set_data(np.clip(slices[k], clip_lo, None))
        else:
            imL.set_data(slices[k])
            imR.set_data(slices[k])
        sup.set_text(f"{title} — t = {times[k]:.1f} a.u.")
        return [imL, imR, sup]

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    an = animation.FuncAnimation(fig, upd, frames=len(slices), blit=False)
    an.save(out_path, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
    return vmax


def make_density_gif_battery(
    results_dir: str,
    out_dir: str,
    *,
    run_label: str,
    dt: float,
    slab_face: float,
    cap_inner: float,
    frames_max: int = 30,        # keep GIFs small enough for notebook viewers to render
    density_vmax: Optional[float] = None,
    fps: int = 10,
    run_title: Optional[str] = None,
    per_frame_norm_wp: bool = False,
    cap_lines: Optional[tuple] = None,
    cats_filter: Optional[list] = None,   # restrict to these categories (e.g. ["total"])
    kinds_filter: Optional[list] = None,  # restrict to these kinds (e.g. ["density"])
):
    """Build the density-GIF battery for one run.

    ``run_label`` keys the output **filenames**; ``run_title`` is the
    human-readable run name baked into the figure title (defaults to a
    prettified ``run_label`` — ``"wp"→"Wavepacket run"``, ``"classical"→
    "Classical run"``).

    Returns ``(gifs, density_vmax)`` where ``gifs`` is an ordered list of
    ``(category, kind, path, caption)`` and ``density_vmax`` is the shared
    total/bath density scale actually used (reuse it across runs to lock scales).
    """
    if run_title is None:
        run_title = {"wp": "Wavepacket run",
                     "classical": "Classical run"}.get(
                         run_label, run_label.replace("_", " ").title())
    raw = os.path.join(results_dir, "raw", "vti")
    os.makedirs(out_dir, exist_ok=True)
    # default: symmetric two-sided slab+CAP lines. A one-sided CAP (e.g. the
    # vacuum +z-only absorber) passes an explicit cap_lines so no spurious -z line
    # is drawn.
    if cap_lines is None:
        cap_lines = (slab_face, -slab_face, cap_inner, -cap_inner)

    tot_dir = os.path.join(raw, "density_total")
    wp_dir = os.path.join(raw, "density_wp")
    has_wp = os.path.isdir(wp_dir) and bool(glob.glob(os.path.join(wp_dir, "*.vti")))

    nfiles = len(glob.glob(os.path.join(tot_dir, "*.vti")))
    if nfiles == 0:
        return [], density_vmax
    idx = list(range(0, nfiles, max(1, nfiles // frames_max)))

    times, tot, axes = _slice_stack(tot_dir, idx, dt)
    cats = {"total": tot}
    if has_wp:
        _, wp, _ = _slice_stack(wp_dir, idx, dt)
        cats["wp"] = wp
        cats["bath"] = tot - wp

    # shared total/bath density scale (slab-tuned). If not supplied, derive from
    # the bath (slab only, no WP spike) when available, else from total.
    base = cats.get("bath", tot)
    if density_vmax is None:
        # global stack (not middle frame) → clim fixed + representative of every frame
        density_vmax = float(np.percentile(base, 99.7))

    KIND_TTL = {"density": "density  n(x,z,t)",
                "delta0": "Δn = n(t) − n(0)",
                "dstep": "Δn = n(t+dt) − n(t)"}
    CAT_TTL = {"total": "Total system", "wp": "Wavepacket |ψ|²", "bath": "Bath (slab only)"}
    gifs = []
    for cat, stack in cats.items():
        if cats_filter is not None and cat not in cats_filter:
            continue
        for kind in ("density", "delta0", "dstep"):
            if kinds_filter is not None and kind not in kinds_filter:
                continue
            if kind == "density":
                series = stack
                vmax = density_vmax if cat in ("total", "bath") else None  # wp: own scale
            elif kind == "delta0":
                series = stack - stack[0][None]
                vmax = None
            else:  # dstep
                series = np.diff(stack, axis=0)
                vmax = None
            label_kind = KIND_TTL[kind]
            if cat == "total" and kind == "delta0" and not has_wp:
                label_kind = "Δn = n(t) − n(0)  (induced wake)"
            title = f"{run_title} · {CAT_TTL[cat]} · {label_kind}"
            fname = f"{run_label}_{cat}_{kind}.gif"
            out = os.path.join(out_dir, fname)
            tt = times[1:] if kind == "dstep" else times
            pfn = per_frame_norm_wp and cat == "wp" and kind == "density"
            _save_gif(series, tt, axes, out, title=title, cap_lines=cap_lines,
                      kind=("density" if kind == "density" else "diff"), vmax=vmax,
                      fps=fps, per_frame_norm=pfn)
            gifs.append((cat, kind, out, title))
    return gifs, density_vmax


def make_twin_density_matrix(
    classical_dir: str,
    wp_dir: str,
    out_dir: str,
    *,
    dt: float,
    slab_face: float,
    cap_inner: Optional[float] = None,
    frames_max: int = 30,
    fps: int = 10,
    total_subpath: str = "frames/total",
):
    """Classical-vs-WP density MATRIX on the mid-y xz **total-density** slice.

    Builds the 3×3 matrix requested for a twin pair — rows {classical, wavepacket,
    WP−classical} × columns {density, induced, instantaneous} — as separate GIFs,
    each a LINEAR | LOG (or symlog) panel pair with the slab faces marked. Reuses
    the canonical slice/render primitives (:func:`_slice_stack`, :func:`_save_gif`)
    so the conventions match the per-run battery exactly (physical-order VTIs, no
    fftshift; global-stack colour scales).

    Column definitions (locked with the user, 2026-07-15):
      - ``density``       : n(x,z,t)
      - ``induced``       : Δn = n(t) − n(0)                     (delta0)
      - ``instantaneous`` : Δn = n(t) − n(t−Δt)                  (dstep)

    Shared scales (shared-colorbar rule): the classical and WP rows share ONE scale
    per column (they are directly compared); the WP−classical difference row owns
    its own symmetric scale. ``total_subpath`` locates the total-density VTIs under
    each run dir (these runs use ``frames/total``, not the battery's ``raw/vti``).

    Returns an ordered list of ``(row, column, path, title)``.
    """
    os.makedirs(out_dir, exist_ok=True)
    cl_tot = os.path.join(classical_dir, total_subpath)
    wp_tot = os.path.join(wp_dir, total_subpath)
    ncl = len(glob.glob(os.path.join(cl_tot, "*.vti")))
    nwp = len(glob.glob(os.path.join(wp_tot, "*.vti")))
    if ncl == 0 or nwp == 0:
        return []
    n = min(ncl, nwp)
    idx = list(range(0, n, max(1, n // frames_max)))

    times, NC, axes = _slice_stack(cl_tot, idx, dt)
    _, NW, _ = _slice_stack(wp_tot, idx, dt)
    cap_lines = (slab_face, -slab_face) if cap_inner is None \
        else (slab_face, -slab_face, cap_inner, -cap_inner)

    # column series for each row
    def series(stack, col):
        if col == "density":
            return stack, times
        if col == "induced":
            return stack - stack[0][None], times
        return np.diff(stack, axis=0), times[1:]      # instantaneous

    COL_TTL = {"density": "density  n(x,z,t)",
               "induced": "induced  Δn = n(t) − n(0)",
               "instantaneous": "instantaneous  Δn = n(t) − n(t−Δt)"}
    ROW_TTL = {"classical": "Classical", "wp": "Wavepacket", "wp_minus_cl": "WP − classical"}
    gifs = []
    for col in ("density", "induced", "instantaneous"):
        Scl, tt = series(NC, col)
        Swp, _ = series(NW, col)
        Sdif = Swp - Scl
        # shared classical/WP scale for this column (difference row owns its own)
        if col == "density":
            shared = float(np.percentile(np.concatenate([NC, NW]), 99.7))
            cl_kind = wp_kind = "density"
        else:
            shared = float(np.percentile(np.abs(np.concatenate([Scl, Swp])), 99.0)) or 1e-12
            cl_kind = wp_kind = "diff"
        rows = [("classical", Scl, cl_kind, shared),
                ("wp", Swp, wp_kind, shared),
                ("wp_minus_cl", Sdif, "diff", None)]   # diff row: own robust scale
        for row, S, kind, vmax in rows:
            title = f"{ROW_TTL[row]} · {COL_TTL[col]}"
            out = os.path.join(out_dir, f"matrix_{row}_{col}.gif")
            _save_gif(S, tt, axes, out, title=title, cap_lines=cap_lines,
                      kind=kind, vmax=vmax, fps=fps)
            gifs.append((row, col, out, title))
    return gifs
