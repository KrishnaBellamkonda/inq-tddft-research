#!/usr/bin/env python3
"""build_rs6_velocity_notebook.py — r_s=6 three-velocity deep-look notebook.

User-specified (2026-07-08): ONE notebook for the r_s=6 tube runs with a section per
velocity v∈{0.15,0.30,0.45}. Per velocity:
  * the FULL density-GIF matrix {n, Δn=n(t)−n(0), Δn=n(t+dt)−n(t)} rendered in ALL
    THREE orthogonal planes {xy (z=mid, face-on annulus+bore), xz (y=0), yz (x=0)};
  * total current density current_{x,y,z}(t) on one axis;
  * total energy energy_total(t).
NO stopping-power section (deferred, per the user).

TUBE-aware (this is a PERIODIC annular tube, NOT a slab / no CAP): the wall is the
annulus R_in=5, R_out=13 Bohr — drawn as CIRCLES in the xy plane and as vertical
bands at ±5,±13 in the axial (xz, yz) planes. Densities load via the canonical
`inqview.load_vti` (PHYSICAL order — NEVER fftshift; vti-coordinate-mapping rule).

Run (venv + stack on path):
    PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_rs6_velocity_notebook.py
Writes ``rs6_velocity_sweep.ipynb`` + figures under ``rs6_velocity_figs/`` beside it
(file-placement: run-tied analysis lives in hypotheses/annular_sv/).
"""
from __future__ import annotations
import glob
import os
import sys
from pathlib import Path

import numpy as np

STACK = "/local/data/public/skcb2/tddft/inq-stack/python"
if STACK not in sys.path:
    sys.path.insert(0, STACK)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402
import nbformat as nbf  # noqa: E402
import pandas as pd  # noqa: E402

from inqview import load_vti  # noqa: E402  canonical physical-order loader
from inqview.visualisation import style  # noqa: E402
try:
    style.apply()
except Exception:
    pass

SYS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium")
SWEEP = SYS / "annular_sv"
HERE = SYS / "hypotheses" / "annular_sv"
RS = int(sys.argv[1]) if len(sys.argv) > 1 else 6      # wall-density index (6, 4, 2)
GEO = {6: dict(N=24, L_z=48, launch=-23),
       4: dict(N=48, L_z=28, launch=-13),
       2: dict(N=136, L_z=10, launch=-4)}
FIGROOT = HERE / f"rs{RS}_velocity_figs"
R_IN, R_OUT = 5.0, 13.0
DT = 0.02
FRAMES = 36            # GIF frames
FPS = 8
VELS = [0.15, 0.30, 0.45]
RUN = {v: f"rs{RS}_v{v:.2f}".replace(".", "p") for v in VELS}


def find(root: Path, name: str):
    return next(Path(root).glob(f"**/{name}"), None)


# ------------------------------------------------------------- slice extraction
def load_planes(dens_dir: Path, n_frames: int):
    """Load n_frames evenly-spaced density VTIs ONCE each and carve the three
    orthogonal mid-plane slices. Returns (times, {plane: stack[T,a,b]}, axes)."""
    files = sorted(glob.glob(os.path.join(str(dens_dir), "*.vti")))
    if not files:
        return None, None, None
    idx = list(range(0, len(files), max(1, len(files) // n_frames)))[:n_frames]
    first = load_vti(files[idx[0]])
    x, y, z = np.asarray(first.x), np.asarray(first.y), np.asarray(first.z)
    ix0, iy0, izm = int(np.argmin(np.abs(x))), int(np.argmin(np.abs(y))), len(z) // 2
    times = np.array([int(_step(files[k])) * DT for k in idx])
    xy, xz, yz = [], [], []
    for k in idx:
        d = np.asarray(load_vti(files[k]).data)   # (nx, ny, nz), physical order
        xy.append(d[:, :, izm])                   # [x, y]
        xz.append(d[:, iy0, :])                   # [x, z]
        yz.append(d[ix0, :, :])                   # [y, z]
    planes = {"xy": (np.array(xy), x, y), "xz": (np.array(xz), x, z),
              "yz": (np.array(yz), y, z)}
    return times, planes, (x, y, z, z[izm])


def _step(path: str) -> int:
    import re
    m = re.search(r"_t?(\d+)\.vti$", os.path.basename(path))
    return int(m.group(1)) if m else 0


# ------------------------------------------------------------------- wall + ion
PLANE_LABEL = {
    "xy": ("x (Bohr)", "y (Bohr)", "z = {zmid:.0f} — face-on annulus + hollow bore"),
    "xz": ("x (Bohr)", "z — tube axis (Bohr)", "y = 0 — axial cut"),
    "yz": ("y (Bohr)", "z — tube axis (Bohr)", "x = 0 — axial cut"),
}


def _draw_wall(ax, plane):
    if plane == "xy":
        for r in (R_IN, R_OUT):
            ax.add_patch(Circle((0, 0), r, fill=False, ls="--", lw=0.8, ec="0.35"))
    else:  # axial planes: wall bands where the cut crosses the annulus
        for a in (R_IN, -R_IN, R_OUT, -R_OUT):
            ax.axvline(a, ls="--", lw=0.7, color="0.4")


def _ion_ab(plane, times, trk):
    """Projectile (a,b) in the plane's coords, interpolated onto frame times."""
    t = trk["time_au"].to_numpy()
    xi = np.interp(times, t, trk["x"].to_numpy())
    yi = np.interp(times, t, trk["y"].to_numpy())
    zi = np.interp(times, t, trk["z"].to_numpy())
    if plane == "xy":
        return np.column_stack([xi, yi]), zi           # + z for in-plane test
    if plane == "xz":
        return np.column_stack([xi, zi]), None
    return np.column_stack([yi, zi]), None


def save_gif(stack, times, aax, bax, plane, kind, out, title, ion_ab, ion_z, zmid):
    """Animate one plane×kind GIF with tube walls + the on-axis projectile."""
    fig, ax = plt.subplots(figsize=(4.8, 4.4) if plane == "xy" else (4.4, 4.8))
    ext = [aax.min(), aax.max(), bax.min(), bax.max()]
    if kind == "density":
        vmax = float(np.percentile(stack[len(stack) // 2], 99.7)) or 1e-12
        vmin = vmax * 1e-3
        im = ax.imshow(np.clip(stack[0].T, vmin, None), origin="lower", extent=ext,
                       aspect="equal" if plane == "xy" else "auto",
                       cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax))
        cbl = "n (a₀⁻³, log)"
    else:
        vmax = float(np.percentile(np.abs(stack), 99.0)) or 1e-12
        im = ax.imshow(stack[0].T, origin="lower", extent=ext,
                       aspect="equal" if plane == "xy" else "auto",
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        cbl = "Δn (a₀⁻³)"
    _draw_wall(ax, plane)
    mk, = ax.plot([], [], "o", mfc="cyan", mec="k", mew=0.8, ms=9, label="projectile")
    trail = None
    if plane != "xy":
        trail, = ax.plot([], [], "-", color="cyan", lw=1.3, alpha=0.8)
    ax.legend(loc="upper right", fontsize=7, framealpha=0.5)
    xl, yl, sub = PLANE_LABEL[plane]
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label(cbl, fontsize=8)
    ttl = ax.set_title("", fontsize=8)

    def upd(k):
        frame = np.clip(stack[k].T, vmin, None) if kind == "density" else stack[k].T
        im.set_data(frame)
        ttl.set_text(f"{title} — t = {times[k]:.1f} a.u.")
        arts = [im, ttl]
        if plane == "xy":
            near = abs(ion_z[k] - zmid) < 1.5      # projectile piercing this plane?
            mk.set_data([ion_ab[k, 0]], [ion_ab[k, 1]])
            mk.set_alpha(1.0 if near else 0.25)
            arts.append(mk)
        else:
            mk.set_data([ion_ab[k, 0]], [ion_ab[k, 1]])
            trail.set_data(ion_ab[:k + 1, 0], ion_ab[:k + 1, 1])
            arts += [mk, trail]
        return arts

    an = animation.FuncAnimation(fig, upd, frames=len(stack), blit=False)
    an.save(out, writer=animation.PillowWriter(fps=FPS))
    plt.close(fig)


KIND = {"density": ("n(x·,t)", "density"),
        "delta0": ("Δn = n(t) − n(0)  (induced)", "diff"),
        "dstep": ("Δn = n(t+dt) − n(t)", "diff")}


def velocity_figs(v):
    label = RUN[v]
    rd = SWEEP / label
    dens = find(rd, "density_system")
    obs = find(rd, "observables.csv")
    trk = find(rd, "electron_track.csv")
    if not (dens and obs and trk):
        print(f"[{label}] MISSING inputs (dens={bool(dens)} obs={bool(obs)} trk={bool(trk)})")
        return None
    outdir = FIGROOT / label
    outdir.mkdir(parents=True, exist_ok=True)
    times, planes, axes = load_planes(dens, FRAMES)
    zmid = axes[3]
    T = pd.read_csv(trk).drop_duplicates("step")

    gifs = {}   # plane -> [(kind, relpath, caption)]
    for plane, (raw, aax, bax) in planes.items():
        ion_ab, ion_z = _ion_ab(plane, times, T)
        series = {"density": raw, "delta0": raw - raw[0][None],
                  "dstep": np.diff(raw, axis=0)}
        gifs[plane] = []
        for kind, arr in series.items():
            tt = times[1:] if kind == "dstep" else times
            iab = ion_ab[1:] if kind == "dstep" else ion_ab
            iz = (ion_z[1:] if (ion_z is not None and kind == "dstep") else ion_z)
            klab, _ = KIND[kind]
            f = outdir / f"{plane}_{kind}.gif"
            save_gif(arr, tt, aax, bax, plane, kind, str(f),
                     f"{label} · {plane} · {klab}", iab, iz, zmid)
            rel = os.path.relpath(f, HERE)
            gifs[plane].append((kind, rel, f"{plane} plane — {klab}"))
        print(f"[{label}] {plane}: 3 GIFs")

    # current density current_{x,y,z}(t)
    O = pd.read_csv(obs).drop_duplicates("step").sort_values("time_au")
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for c, col in [("C0", "current_x"), ("C1", "current_y"), ("C3", "current_z")]:
        if col in O:
            ax.plot(O["time_au"], O[col], c, lw=1.3, label=col)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("total current density (a.u.)")
    ax.set_title(f"{label}: induced total current density vs time")
    ax.legend(); ax.grid(alpha=.25)
    fcur = outdir / "current_xyz.png"
    fig.tight_layout(); fig.savefig(fcur, dpi=150); plt.close(fig)

    # total energy energy_total(t)
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.plot(O["time_au"], O["energy_total"], "C2-", lw=1.4)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$E_\mathrm{total}$ (Ha)")
    ax.set_title(f"{label}: total electronic energy vs time")
    ax.grid(alpha=.25)
    fen = outdir / "energy_total.png"
    fig.tight_layout(); fig.savefig(fen, dpi=150); plt.close(fig)
    print(f"[{label}] current + energy plots")

    return dict(label=label, gifs=gifs, current=os.path.relpath(fcur, HERE),
                energy=os.path.relpath(fen, HERE), zmid=zmid)


# --------------------------------------------------------------------- assemble
def build():
    import math
    g = GEO[RS]
    n0 = 3.0 / (4 * math.pi * RS ** 3)
    dens_word = {6: "dilute", 4: "intermediate", 2: "dense"}[RS]
    FIGROOT.mkdir(parents=True, exist_ok=True)
    results = {v: velocity_figs(v) for v in VELS}

    cells = [nbf.v4.new_markdown_cell(
        f"# Cylindrical (annular) jellium — r_s = {RS}, three-velocity deep look\n\n"
        f"*The {dens_word}-wall (r_s = {RS}) tube run at all three projectile velocities "
        "v ∈ {0.15, 0.30, 0.45} a.u., one section each. For every velocity: the full "
        "density-GIF matrix in the three orthogonal planes, the induced total current "
        "density, and the total electronic energy. Stopping power is deferred (a later "
        "conversation).*")]

    cells.append(nbf.v4.new_markdown_cell(
        "## Setup & conventions\n\n"
        f"**System.** Periodic annular jellium tube, axis ∥ z. Wall between R_in = 5 and "
        f"R_out = 13 Bohr (8 Bohr wall); hollow bore along z. r_s = {RS} ⇒ N = {g['N']} "
        f"electrons, L_z = {g['L_z']} Bohr, L_xy ≈ 40, dx = 0.5, n₀ = {n0:.5f} a₀⁻³. "
        f"**Projectile:** classical electron (Gaussian UPF σ_pot = 0.354 = σ_WP/√2, mass "
        f"mₑ, **free Ehrenfest**), launched on-axis (x=y=0) at z ≈ {g['launch']}, drifting "
        f"+z. LDA, dt = 0.02.\n\n"
        "**Planes.** `xy` (z = mid): the annulus face-on — wall = two dashed **circles** "
        "r = 5, 13; the projectile pierces this plane on-axis (marker brightens as it "
        "crosses). `xz` (y = 0) and `yz` (x = 0): axial cuts — the wall becomes vertical "
        "dashed bands at ±5, ±13, and the projectile glides up z (cyan marker + trail).\n\n"
        "**Density kinds.** `n` = total density (log scale); `Δn = n(t) − n(0)` = the "
        "**induced** density; `Δn = n(t+dt) − n(t)` = per-step change (instantaneous flux). "
        "VTIs load via `inqview.load_vti` (physical order, never fftshift).\n\n"
        "**Source files.** run: "
        "[`scripts/annular_sv/classical/run.cpp`](../../scripts/annular_sv/classical/run.cpp) · "
        "geometry: [`shared/configs/annular_tube.hpp`](../../shared/configs/annular_tube.hpp) · "
        f"builder: [`build_velocity_notebook.py`](build_velocity_notebook.py) (RS={RS}) · "
        "run-SET report: [`annular_sv_report.ipynb`](annular_sv_report.ipynb) · "
        "guided index: [`annular_sv_index.ipynb`](annular_sv_index.ipynb)."))

    for v in VELS:
        r = results[v]
        if r is None:
            cells.append(nbf.v4.new_markdown_cell(
                f"## v = {v:.2f} a.u.\n\n*(inputs missing — run not built)*"))
            continue
        cells.append(nbf.v4.new_markdown_cell(
            f"# Section — v₀ = {v:.2f} a.u.  (`{r['label']}`)\n\n"
            f"Projectile launched at v₀ = {v:.2f}; the light electron decelerates under "
            f"the wall drag as it glides down the bore."))

        cells.append(nbf.v4.new_markdown_cell(
            "## Density-GIF matrix — three planes × three kinds\n\n"
            "Each row is one plane; each column one density kind. Read the wall by the "
            "dashed circles (xy) / bands (xz, yz)."))
        for plane in ("xy", "xz", "yz"):
            _, _, subt = PLANE_LABEL[plane]
            head = f"### {plane} plane — {subt.format(zmid=r['zmid'])}\n"
            imgs = "\n\n".join(
                f"*{cap}*\n\n![{cap}]({rel})" for kind, rel, cap in r["gifs"][plane])
            cells.append(nbf.v4.new_markdown_cell(head + "\n" + imgs))

        cells.append(nbf.v4.new_markdown_cell(
            "## Total current density — current_x, current_y, current_z\n\n"
            f"![current xyz]({r['current']})"))
        cells.append(nbf.v4.new_markdown_cell(
            "## Total electronic energy vs time\n\n"
            f"![energy total]({r['energy']})"))

    cells.append(nbf.v4.new_markdown_cell(
        "## Note\n\nStopping-power extraction is intentionally omitted here (deferred to a "
        "later conversation). See the run-SET report / guided index for S(v) and β(r_s)."))

    nb = nbf.v4.new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})
    out = HERE / f"rs{RS}_velocity_sweep.ipynb"
    nbf.write(nb, str(out))
    print(f"wrote {out}  ({len(cells)} cells)")


if __name__ == "__main__":
    build()
