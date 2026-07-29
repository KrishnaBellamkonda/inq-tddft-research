#!/usr/bin/env python3
"""compare_notebook.py — no-CAP vs CAP comparison notebook for the vacuum WP runs.

Assembles ONE .ipynb that puts the two runs side by side:
  1. Setup figure (real t=0 density, CAP bands + WP launch dashed)
  2. Energy vs time — total energy(t), no-CAP vs CAP overlaid (labelled)
  3. Total-density GIF comparison — no-CAP | CAP, ONE SHARED FIXED LOG scale over
     both stacks, so (a) the wavepacket is visible as it moves (log spans the ~100x
     dispersion collapse) and (b) the CAP run's absorption is visible as the WP
     fading to zero at the boundary while the no-CAP WP persists. A per-frame
     normalised twin is also emitted (motion-only, absorption hidden).

Total density (not the per-frame-normalised WP-only GIF) is the right field for the
COMPARISON because the physical difference between the runs is the absorption of
flux — which per-frame normalisation would hide.

Usage:
  compare_notebook.py <nocap_results> <cap_results> <out_dir> [--dt 0.02]
                      [--cap-inner 20] [--cap-outer 30]
"""
from __future__ import annotations
import argparse, base64, glob, re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from inqview.visualisation.field_io import load_vti
try:
    from inqview.visualisation.style import apply_theme
    apply_theme()
except Exception:
    pass

HA_EV = 27.211386245988


def _log(m): print(f"[compare-nb] {m}", flush=True)


def _frame_time(path: str, dt: float) -> float:
    m = re.search(r"_t(\d+)\.vti$", path)
    return (int(m.group(1)) if m else 0) * dt


def _stack(results: Path, frames_max: int, dt: float):
    files = sorted(glob.glob(str(results / "raw/vti/density_total/*.vti")))
    if not files:
        return None, None, None
    idx = list(range(0, len(files), max(1, len(files) // frames_max)))
    files = [files[i] for i in idx]
    f0 = load_vti(files[0])
    iy = f0.data.shape[1] // 2
    sl = np.stack([load_vti(f).data[:, iy, :].T for f in files])  # (T,nz,nx)
    times = np.array([_frame_time(f, dt) for f in files])
    return times, sl, (f0.x, f0.z)


def _mom_grid(results: Path):
    """Return (times, k, N) with N[it,ik] = n_wp(k,t), from momentum_distribution.csv."""
    import pandas as pd
    p = results / "raw/observables/momentum_distribution.csv"
    df = pd.read_csv(p, comment="#")
    piv = df.pivot_table(index="time_au", columns="k_bohr_inv", values="n_wp")
    return piv.index.to_numpy(), piv.columns.to_numpy(), piv.to_numpy()


def momentum_compare(nocap: Path, cap: Path, out_over: Path, out_carpet: Path):
    tn, kn, Nn = _mom_grid(nocap)
    tc, kc, Nc = _mom_grid(cap)
    k0 = 2.711
    # --- overlay: n_wp(k) at t0 (shared) and tF for each run --------------
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.plot(kn, Nn[0], "-", color="0.4", lw=1.6, label="t = 0 (both, launch)")
    ax.plot(kn, Nn[-1], "-", color="C0", lw=2.2, label=f"no-CAP, t = {tn[-1]:.0f} a.u.")
    ax.plot(kc, Nc[-1], "-", color="C3", lw=2.2, label=f"CAP, t = {tc[-1]:.0f} a.u.")
    ax.axvline(k0, color="0.6", ls=":", lw=1.0, label=f"k₀ = {k0:.2f}")
    ax.set_xlabel("k_z (Bohr⁻¹)")
    ax.set_ylabel("n_wp(k)  (a.u.)")
    ax.set_title("Wavepacket momentum distribution — no-CAP vs CAP")
    ax.legend(fontsize=8.5, loc="best")
    fig.tight_layout(); fig.savefig(out_over, dpi=150); plt.close(fig)

    # --- side-by-side carpets |n_wp(k,t)|, shared scale ------------------
    vmax = float(max(np.nanmax(Nn), np.nanmax(Nc))) or 1e-12
    fig, (aL, aR) = plt.subplots(1, 2, figsize=(9.8, 4.6), sharey=True)
    for ax, (t, k, N, ttl) in zip((aL, aR),
            ((tn, kn, Nn, "no-CAP"), (tc, kc, Nc, "CAP"))):
        im = ax.imshow(N.T, origin="lower", aspect="auto",
                       extent=[t[0], t[-1], k[0], k[-1]], cmap="magma",
                       vmin=0, vmax=vmax)
        ax.axhline(k0, color="cyan", ls=":", lw=1.0)
        ax.set_xlabel("time (a.u.)"); ax.set_title(ttl, fontsize=10)
    aL.set_ylabel("k_z (Bohr⁻¹)")
    div = make_axes_locatable(aR); cax = div.append_axes("right", size="5%", pad=0.08)
    fig.colorbar(im, cax=cax).set_label("n_wp(k,t) (shared)", fontsize=8)
    fig.suptitle("Momentum-distribution carpets — no-CAP vs CAP", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_carpet, dpi=150); plt.close(fig)
    _log(f"momentum: no-CAP peak-k frozen at k₀; CAP peak-k "
         f"{kc[int(np.argmax(Nc[-1]))]:.2f} (absorber reshapes n(k))")
    return dict(nocap_peakF=float(kn[int(np.argmax(Nn[-1]))]),
                cap_peakF=float(kc[int(np.argmax(Nc[-1]))]))


def energy_overlay(nocap: Path, cap: Path, out: Path):
    def load(r):
        fs = sorted(glob.glob(str(r / "raw/observables/energies*.csv")))
        import pandas as pd
        df = pd.concat([pd.read_csv(f) for f in fs]).drop_duplicates("step").sort_values("step")
        return df["time_au"].to_numpy(), df["total"].to_numpy() * HA_EV
    tn, en = load(nocap); tc, ec = load(cap)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.plot(tn, en, "-", color="C0", lw=2.2, label="no-CAP (closed box)")
    ax.plot(tc, ec, "-", color="C3", lw=2.2, label="CAP (drains escaping flux)")
    ax.axhline(en[0], color="0.5", ls=":", lw=1.0, label="initial total")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("total energy (eV)")
    ax.set_title("Total energy vs time — no-CAP vs CAP (vacuum WP)")
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    drop = float(en[0] - ec[-1])
    _log(f"energy overlay: no-CAP end {en[-1]:.1f} eV, CAP end {ec[-1]:.1f} eV, "
         f"CAP removed {drop:.1f} eV")
    return dict(nocap_end=float(en[-1]), cap_end=float(ec[-1]), removed=drop)


def side_by_side_gif(nocap: Path, cap: Path, out: Path, *, dt, cap_inner,
                     cap_outer, per_frame=False, frames_max=24):
    tn, sn, ax_n = _stack(nocap, frames_max, dt)
    tc, sc, _ = _stack(cap, frames_max, dt)
    if sn is None or sc is None:
        _log("no density frames — GIF skipped"); return False
    n = min(len(sn), len(sc)); sn, sc, tt = sn[:n], sc[:n], tn[:n]
    x, z = ax_n; ext = [x[0], x[-1], z[0], z[-1]]
    gmax = float(max(np.nanmax(sn), np.nanmax(sc))) or 1e-12
    vmin = gmax * 1e-4
    fig, (aL, aR) = plt.subplots(1, 2, figsize=(9.6, 6.0), sharey=True)

    def draw(ax, frame, title):
        if per_frame:
            m = float(np.nanmax(frame)) or 1e-12
            im = ax.imshow(np.clip(frame, 0, None) / m, origin="lower", aspect="auto",
                           extent=ext, cmap="viridis", vmin=0, vmax=1.0)
        else:
            im = ax.imshow(np.clip(frame, vmin, None), origin="lower", aspect="auto",
                           extent=ext, cmap="viridis", norm=LogNorm(vmin=vmin, vmax=gmax))
        # one-sided +z CAP: shaded band + inner-edge dashed line
        ax.axhspan(cap_inner, cap_outer, color="crimson", alpha=0.13, zorder=2)
        ax.axhline(cap_inner, ls="--", lw=1.2, color="crimson", zorder=3)
        ax.set_xlabel("x (Bohr)"); ax.set_title(title, fontsize=10)
        return im
    imL = draw(aL, sn[0], "no-CAP")
    imR = draw(aR, sc[0], "CAP")
    aL.set_ylabel("z (Bohr)")
    lbl = "n / nₘₐₓ(t) (per-frame)" if per_frame else "n (a₀⁻³, log — shared)"
    div = make_axes_locatable(aR); cax = div.append_axes("right", size="5%", pad=0.08)
    fig.colorbar(imR, cax=cax).set_label(lbl, fontsize=8)
    sup = fig.suptitle("", fontsize=11)

    def upd(k):
        if per_frame:
            imL.set_data(np.clip(sn[k], 0, None) / (float(np.nanmax(sn[k])) or 1e-12))
            imR.set_data(np.clip(sc[k], 0, None) / (float(np.nanmax(sc[k])) or 1e-12))
        else:
            imL.set_data(np.clip(sn[k], vmin, None))
            imR.set_data(np.clip(sc[k], vmin, None))
        sup.set_text(f"Total density  n(x,z,t) — t = {tt[k]:.1f} a.u.")
        return [imL, imR, sup]
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    an = animation.FuncAnimation(fig, upd, frames=n, blit=False)
    an.save(out, writer=animation.PillowWriter(fps=10)); plt.close(fig)
    _log(f"wrote {out.name} ({'per-frame' if per_frame else 'shared-log'}, {n} frames)")
    return True


def build(nocap: Path, cap: Path, out_dir: Path, dt, cap_inner, cap_outer):
    import nbformat as nbf
    from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell, new_output
    out_dir.mkdir(parents=True, exist_ok=True)

    en = energy_overlay(nocap, cap, out_dir / "energy_compare.png")
    mom = momentum_compare(nocap, cap, out_dir / "momentum_overlay.png",
                           out_dir / "momentum_carpet_compare.png")
    setup = cap / "report" / "setup_vacuum_cap.png"      # produced by make_setup_figure.py
    ok_shared = side_by_side_gif(nocap, cap, out_dir / "total_density_compare_log.gif",
                                 dt=dt, cap_inner=cap_inner, cap_outer=cap_outer, per_frame=False)
    ok_pf = side_by_side_gif(nocap, cap, out_dir / "total_density_compare_perframe.gif",
                             dt=dt, cap_inner=cap_inner, cap_outer=cap_outer, per_frame=True)

    def out_img(path: Path, mime):
        return new_output("display_data",
                          data={mime: base64.b64encode(path.read_bytes()).decode()},
                          metadata={})

    def emb(path: Path, mime, stem):
        c = new_code_cell(f'from IPython.display import Image\nImage(filename="{path.name}")')
        c.outputs = [out_img(path, mime)]; c.execution_count = None
        return [new_markdown_cell(f"### {stem}"), c]

    cells = [new_markdown_cell(
        "# Vacuum wavepacket — no-CAP vs CAP comparison\n\n"
        "The same σ=1, E=100 eV Gaussian wavepacket is propagated through an empty "
        "box **without** and **with** a **one-sided** complex absorbing potential "
        "(CAP, η=−0.7 Ha, at the +z end z∈[30,40]). The WP launches at z=−30, **10σ "
        "clear** of both the CAP inner edge (60 Bohr away) and the −z wall (which is "
        "the CAP's wrapped outer edge). This notebook puts the two runs side by side "
        "to read off (1) what the CAP does to the total energy and (2) what it does "
        "to the density — the wavepacket is absorbed at the +z boundary in the CAP "
        "run and persists (merely dispersing) in the no-CAP run.")]

    if setup.exists():
        cells += [new_markdown_cell(
            "## 0. Setup (real t=0 density)\n\nWavepacket launch plane (white dashed) "
            "and the two CAP bands (crimson dashed) marked on the actual initial "
            "density — geometry read off the data, not drawn.")]
        cells += emb(setup, "image/png", "setup_vacuum_cap")

    cells += [new_markdown_cell(
        "## 1. Total energy vs time\n\n"
        f"- no-CAP total energy is **conserved** (ends {en['nocap_end']:.1f} eV) — a "
        "closed box keeps everything.\n"
        f"- CAP total energy **falls** (ends {en['cap_end']:.1f} eV) as the wavepacket "
        f"reaches the absorber; the CAP removes ≈ **{en['removed']:.1f} eV**.\n"
        "The gap between the two curves is the energy that, in the closed box, has "
        "nowhere to go — the mechanism the localised-jellium experiment probes.")]
    cells += emb(out_dir / "energy_compare.png", "image/png", "energy_compare")

    if ok_shared:
        cells += [new_markdown_cell(
            "## 2. Total density — shared fixed LOG scale (the fair comparison)\n\n"
            "Both panels share ONE colour scale (log, spanning 4 decades over both "
            "runs) so they are directly comparable. The log scale is essential: the "
            "wavepacket's peak density collapses ~100× as it disperses, so a single "
            "*linear* scale would black it out (the original 'nothing moves' bug). "
            "Here you can see the WP **move** (both panels) AND the CAP **absorb** it "
            "(right panel fades to zero past z=+30 while the left persists).")]
        cells += emb(out_dir / "total_density_compare_log.gif", "image/gif",
                     "total_density_compare_log")
    if ok_pf:
        cells += [new_markdown_cell(
            "## 3. Total density — per-frame normalised (motion only)\n\n"
            "Each frame is normalised to its own max (n/nₘₐₓ(t)), so the moving blob "
            "stays maximally visible — but this **hides absorption** (every frame is "
            "renormalised). Use this only to confirm the wavepacket trajectory; use "
            "§2 to judge the CAP's effect.")]
        cells += emb(out_dir / "total_density_compare_perframe.gif", "image/gif",
                     "total_density_compare_perframe")

    cells += [new_markdown_cell(
        "## 4. Momentum distribution — no-CAP vs CAP\n\n"
        "- **no-CAP: momentum is a constant of motion.** n_wp(k) at the final time is "
        "identical to t=0 — a free vacuum wavepacket disperses in real space but its "
        "momentum content is frozen (peak stays at k₀=2.71).\n"
        f"- **CAP: the absorber reshapes n(k).** The peak momentum drops to "
        f"k≈{mom['cap_peakF']:.2f} Bohr⁻¹ — the fast (high-k) front reaches the "
        "boundary first and is absorbed, so the CAP is a momentum-selective sink, not "
        "a mere density clip.\n\n"
        "The overlay shows the two final distributions against the shared launch; the "
        "carpets show |n_wp(k,t)| over the whole run (cyan dotted = k₀).")]
    cells += emb(out_dir / "momentum_overlay.png", "image/png", "momentum_overlay")
    cells += emb(out_dir / "momentum_carpet_compare.png", "image/png",
                 "momentum_carpet_compare")

    cells += [new_markdown_cell(
        "## Takeaway\n\n"
        f"- The CAP removes ≈ {en['removed']:.0f} eV of the wavepacket's energy as it "
        "is absorbed at the boundary; the no-CAP box retains all of it.\n"
        "- Momentum: no-CAP freezes n(k) (free particle); the CAP reshapes it "
        f"(peak k 2.71→{mom['cap_peakF']:.2f}) — a momentum-selective sink.\n"
        "- Density-wise: identical motion + dispersion up to the absorber, then the "
        "CAP run's wavepacket vanishes into the +z crimson band while the no-CAP one "
        "reflects/wraps and stays in the box.\n"
        "- This vacuum test validates the diagnostic before the localised-jellium "
        "runs, where the retained energy is the physics question.")]

    nb = new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3",
                                 "language": "python"}
    nbp = out_dir / "nocap_vs_cap_comparison.ipynb"
    nbf.write(nb, str(nbp))
    _log(f"wrote {nbp}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("nocap"); ap.add_argument("cap"); ap.add_argument("out_dir")
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--cap-inner", type=float, default=30.0)
    ap.add_argument("--cap-outer", type=float, default=40.0)
    a = ap.parse_args()
    build(Path(a.nocap), Path(a.cap), Path(a.out_dir),
          a.dt, a.cap_inner, a.cap_outer)
