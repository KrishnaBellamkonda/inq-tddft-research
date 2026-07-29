#!/usr/bin/env python3
"""build_wp_notebook.py — wavepacket run (wp_rs6_v0p30) density + energy/current,
with its matched CLASSICAL analogue (rs6_v0p30) for side-by-side comparison.

User (2026-07-08): density plots + energy/current curves for the WP run, plus the
classical analogue to compare against.

WP-specific density decomposition (the point of a WP run):
  total = density_system (25 e = 24 wall + 1 WP) · wp = density_wp (|ψ_WP|², 1 e) ·
  bath  = total − wp (the 24-electron wall response; canonical bath density).
Each channel gets the full tube-aware matrix {n, Δn=n(t)−n(0), Δn=n(t+dt)−n(t)} in the
three planes {xy (z=mid), xz (y=0), yz (x=0)} → 27 GIFs, with the WP CENTROID
(wp_real_space_stats) overlaid. Plus current_{x,y,z}(t) and energy_total(t).

The classical analogue rs6_v0p30 (same r_s, v) is referenced from its already-built
3-plane figures in rs6_velocity_figs/ (total density only — a classical run has no WP
channel), plus its own current/energy.

Reuses the tube-plane machinery from build_velocity_notebook.py (DRY). Run:
    PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_wp_notebook.py
Writes wp_rs6_v0p30_planes.ipynb + figures under wp_velocity_figs/wp_rs6_v0p30/.
"""
from __future__ import annotations
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import nbformat as nbf

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import build_velocity_notebook as bvn   # reuse load_planes/save_gif/_ion_ab/PLANE_LABEL/KIND

SYS = bvn.SYS
WP_RUN = SYS / "annular_sv" / "wp_rs6_v0p30"
CL_RUN = SYS / "annular_sv" / "rs6_v0p30"          # matched classical analogue (same r_s, v)
FIGROOT = HERE / "wp_velocity_figs" / "wp_rs6_v0p30"
CMPDIR = FIGROOT / "compare"                        # WP-vs-classical comparison figures
CH_LABEL = {"total": "Total system (wall + WP)", "wp": "Wavepacket |ψ|²",
            "bath": "Bath = total − WP (wall response)"}


def _find(root, name):
    return next(Path(root).glob(f"**/{name}"), None)


def wp_centroid_track():
    """A track-like df (time_au, x, y, z) from the WP centroid stats."""
    p = _find(WP_RUN, "wp_real_space_stats.csv")
    R = pd.read_csv(p, comment="#")
    return pd.DataFrame({"time_au": R["time_au"], "x": R["x_mean"],
                         "y": R["y_mean"], "z": R["z_mean"]})


def build_wp_figs():
    FIGROOT.mkdir(parents=True, exist_ok=True)
    dsys = _find(WP_RUN, "density_system")
    dwp = _find(WP_RUN, "density_wp")
    times, tot, axes = bvn.load_planes(dsys, bvn.FRAMES)
    _, wp, _ = bvn.load_planes(dwp, bvn.FRAMES)
    zmid = axes[3]
    trk = wp_centroid_track()
    # bath = total − wp per plane (stacks share grid + frame times)
    bath = {p: (tot[p][0] - wp[p][0], tot[p][1], tot[p][2]) for p in tot}
    channels = {"total": tot, "wp": wp, "bath": bath}

    gifs = {}   # (channel, plane) -> [(kind, relpath, caption)]
    for ch, planes in channels.items():
        for plane, (raw, aax, bax) in planes.items():
            ion_ab, ion_z = bvn._ion_ab(plane, times, trk)
            series = {"density": raw, "delta0": raw - raw[0][None],
                      "dstep": np.diff(raw, axis=0)}
            lst = []
            for kind, arr in series.items():
                tt = times[1:] if kind == "dstep" else times
                iab = ion_ab[1:] if kind == "dstep" else ion_ab
                iz = (ion_z[1:] if (ion_z is not None and kind == "dstep") else ion_z)
                klab = bvn.KIND[kind][0]
                f = FIGROOT / f"{ch}_{plane}_{kind}.gif"
                if not f.exists():                     # reuse existing GIFs (expensive)
                    bvn.save_gif(arr, tt, aax, bax, plane, kind, str(f),
                                 f"WP · {ch} · {plane} · {klab}", iab, iz, zmid)
                lst.append((kind, os.path.relpath(f, HERE), f"{plane} — {klab}"))
            gifs[(ch, plane)] = lst
            print(f"[wp] {ch}/{plane}: 3 GIFs {'(reused)' if all((FIGROOT / f'{ch}_{plane}_{k}.gif').exists() for k in ('density','delta0','dstep')) else ''}")

    # current + energy
    obs = _find(WP_RUN, "observables.csv")
    O = pd.read_csv(obs).drop_duplicates("step").sort_values("time_au")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for c, col in [("C0", "current_x"), ("C1", "current_y"), ("C3", "current_z")]:
        if col in O:
            ax.plot(O["time_au"], O[col], c, lw=1.3, label=col)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("total current density (a.u.)")
    ax.set_title("wp_rs6_v0p30: induced total current density vs time")
    ax.legend(); ax.grid(alpha=.25)
    fcur = FIGROOT / "current_xyz.png"
    fig.tight_layout(); fig.savefig(fcur, dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.plot(O["time_au"], O["energy_total"], "C2-", lw=1.4)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$E_\mathrm{total}$ (Ha)")
    ax.set_title("wp_rs6_v0p30: total electronic energy vs time"); ax.grid(alpha=.25)
    fen = FIGROOT / "energy_total.png"
    fig.tight_layout(); fig.savefig(fen, dpi=150); plt.close(fig)
    print("[wp] current + energy plots")
    return dict(gifs=gifs, current=os.path.relpath(fcur, HERE),
                energy=os.path.relpath(fen, HERE), zmid=zmid)


def build_comparison_figs():
    """Direct WP-vs-classical comparison figures (same r_s=6, v=0.30).

    Physics of the density comparison: the WP `density_system` = wall(24 e) + WP |ψ|²
    (25 e); the classical projectile is a pseudo-ion contributing NO electron density,
    so its `density_system` = wall(24 e). The matched wall-response pair is therefore
    **WP bath = total − wp** vs **classical total** — both the 24-e wall responding to
    the projectile, both starting from the identical wall GS (canonical-bath-density).
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    CMPDIR.mkdir(parents=True, exist_ok=True)
    HA_EV = 27.21138625
    out = {}

    Owp = pd.read_csv(_find(WP_RUN, "observables.csv")).drop_duplicates("step").sort_values("time_au")
    Ocl = pd.read_csv(_find(CL_RUN, "observables.csv")).drop_duplicates("step").sort_values("time_au")

    # ---- Fig A: induced current density, WP vs classical ------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.2))
    a1.plot(Owp["time_au"], Owp["current_z"], "C0-", lw=1.6, label="WP")
    a1.plot(Ocl["time_au"], Ocl["current_z"], "C3--", lw=1.6, label="classical")
    a1.set_xlabel("time (a.u.)"); a1.set_ylabel(r"axial current $J_z$ (a.u.)")
    a1.set_title("axial induced current $J_z(t)$"); a1.legend(frameon=False); a1.grid(alpha=.25)
    for O, c, ls, tag in [(Owp, "C0", "-", "WP"), (Ocl, "C3", "--", "cl")]:
        a2.plot(O["time_au"], O["current_x"], c, ls=ls, lw=1.1, alpha=.9, label=f"{tag} $J_x$")
        a2.plot(O["time_au"], O["current_y"], c, ls=":", lw=1.1, alpha=.7, label=f"{tag} $J_y$")
    a2.set_xlabel("time (a.u.)"); a2.set_ylabel(r"transverse current (a.u.)")
    a2.set_title("transverse $J_x, J_y(t)$ (≈0 by on-axis symmetry)")
    a2.legend(frameon=False, fontsize=7, ncol=2); a2.grid(alpha=.25)
    fA = CMPDIR / "compare_current.png"
    fig.tight_layout(); fig.savefig(fA, dpi=160); plt.close(fig)
    out["current"] = os.path.relpath(fA, HERE)

    # ---- Fig B: total electronic energy ΔE_total(t) — WP conserved vs classical rise
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    for O, c, ls, tag in [(Owp, "C0", "-", "WP"), (Ocl, "C3", "--", "classical")]:
        dE = (O["energy_total"] - O["energy_total"].iloc[0]).to_numpy()
        rng = f"range [{dE.min()*HA_EV:+.2f}, {dE.max()*HA_EV:+.2f}] eV"
        ax.plot(O["time_au"], dE * HA_EV, c, ls=ls, lw=1.7, label=f"{tag}  ({rng})")
    ax.axhline(0, color="0.6", lw=0.7)
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\Delta E_\mathrm{total} = E(t)-E(0)$  (eV)")
    ax.set_title("Total electronic energy: WP conserved vs classical deposit "
                 "(r_s=6, v=0.30)", fontsize=10)
    ax.annotate("WP total energy is CONSERVED\n(projectile is inside the quantum\n"
                "system → its KE is in E_total)", xy=(0.03, 0.05),
                xycoords="axes fraction", va="bottom", fontsize=8, color="C0",
                bbox=dict(boxstyle="round", fc="white", ec="C0", alpha=.8))
    ax.legend(title="projectile", loc="upper left"); ax.grid(alpha=.25)
    fB = CMPDIR / "compare_energy_deposit.png"
    fig.tight_layout(); fig.savefig(fB, dpi=160); plt.close(fig)
    out["energy"] = os.path.relpath(fB, HERE)

    # ---- density stacks (xz propagation plane) ----------------------------------
    times, tot, _ = bvn.load_planes(_find(WP_RUN, "density_system"), bvn.FRAMES)
    _, wp, _ = bvn.load_planes(_find(WP_RUN, "density_wp"), bvn.FRAMES)
    _, cl, _ = bvn.load_planes(_find(CL_RUN, "density_system"), bvn.FRAMES)
    tot_xz, x, z = tot["xz"]; wp_xz, _, _ = wp["xz"]; cl_xz, _, _ = cl["xz"]
    bath_xz = tot_xz - wp_xz                              # WP wall response (24 e)
    ext = [x.min(), x.max(), z.min(), z.max()]

    def _walls(ax):
        for a in (5, -5, 13, -13):
            ax.axvline(a, ls="--", lw=0.6, color="0.4")

    # ---- Fig C: induced Δn grid — WP wall-response vs classical, matched times ---
    bath_ind = bath_xz - bath_xz[0]
    cl_ind = cl_xz - cl_xz[0]
    idx = [8, 17, 26, 34]                                 # 4 representative frames (skip t≈0)
    vmax = float(np.percentile(np.abs(np.concatenate([bath_ind[idx], cl_ind[idx]])), 99.5))
    dvmax = float(np.percentile(np.abs((bath_ind - cl_ind)[idx]), 99.5)) or vmax
    fig, axs = plt.subplots(len(idx), 3, figsize=(11.0, 2.7 * len(idx)),
                            constrained_layout=True)
    col_t = ["WP wall response  Δn(bath)", "classical  Δn(total)", "difference  WP − cl"]
    for r, k in enumerate(idx):
        panels = [(bath_ind[k], vmax, "RdBu_r"), (cl_ind[k], vmax, "RdBu_r"),
                  ((bath_ind - cl_ind)[k], dvmax, "PuOr_r")]
        for cc, (arr, vm, cmap) in enumerate(panels):
            ax = axs[r, cc]
            im = ax.imshow(arr.T, origin="lower", extent=ext, aspect="auto",
                           cmap=cmap, vmin=-vm, vmax=vm)
            _walls(ax)
            if r == 0:
                ax.set_title(col_t[cc], fontsize=9)
            if cc == 0:
                ax.set_ylabel(f"t = {times[k]:.1f}\nz (Bohr)", fontsize=8)
            ax.set_xlabel("x (Bohr)" if r == len(idx) - 1 else "")
            if cc == 2:
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label("Δn diff", fontsize=7)
            elif cc == 1:
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label("Δn (shared)", fontsize=7)
    fig.suptitle("Induced density Δn = n(t)−n(0), xz plane — WP wall-response vs classical "
                 "(cols 1–2 share one linear colour scale)", fontsize=10)
    fC = CMPDIR / "compare_induced_density_xz.png"
    fig.savefig(fC, dpi=150); plt.close(fig)
    out["induced_density"] = os.path.relpath(fC, HERE)

    # ---- Fig D: total density (log), decomposition at one late time -------------
    k = idx[-1]
    fig, axs = plt.subplots(1, 3, figsize=(12.0, 4.2), constrained_layout=True)
    fields = [(tot_xz[k], "WP total (wall + WP, 25 e)"),
              (bath_xz[k], "WP bath = total − WP (24 e)"),
              (cl_xz[k], "classical total (24 e)")]
    vmx = float(np.percentile(np.concatenate([tot_xz[k], bath_xz[k], cl_xz[k]]), 99.7))
    vmn = vmx * 1e-3
    for ax, (arr, ttl) in zip(axs, fields):
        im = ax.imshow(np.clip(arr.T, vmn, None), origin="lower", extent=ext, aspect="auto",
                       cmap="viridis", norm=LogNorm(vmin=vmn, vmax=vmx))
        _walls(ax); ax.set_title(ttl, fontsize=9); ax.set_xlabel("x (Bohr)")
    axs[0].set_ylabel("z (Bohr)")
    fig.colorbar(im, ax=axs, fraction=0.025, pad=0.02).set_label("n (a₀⁻³, log, shared)", fontsize=8)
    fig.suptitle(f"Total density (log), xz plane, t = {times[k]:.1f} a.u. — the WP blob "
                 "on-axis is the only difference between WP-total and classical-total", fontsize=10)
    fD = CMPDIR / "compare_total_density_log_xz.png"
    fig.savefig(fD, dpi=150); plt.close(fig)
    out["total_density"] = os.path.relpath(fD, HERE)

    print(f"[cmp] 4 comparison figures → {CMPDIR}")
    return out


def build_notebook(wp, cmp):
    cells = [nbf.v4.new_markdown_cell(
        "# Wavepacket run wp_rs6_v0p30 — density, energy & current, vs the classical analogue\n\n"
        "*The quantum wavepacket projectile (σ_WP=0.5, k0=0.30) gliding down the r_s=6 "
        "tube bore, compared against its matched CLASSICAL electron ghost (rs6_v0p30, "
        "same r_s and v=0.30). Both are on-axis, free Ehrenfest, LDA, dt=0.02.*\n\n"
        "**WP density decomposition:** total = wall + WP (25 e); wp = |ψ_WP|² (1 e); "
        "bath = total − WP (the 24-electron wall response). Planes: xy (face-on annulus, "
        "walls = dashed circles r=5,13), xz/yz (axial, wall bands at ±5,±13). The WP "
        "centroid is overlaid (cyan). Stopping power deferred.")]

    # WP density matrix, grouped by channel then plane
    cells.append(nbf.v4.new_markdown_cell("# Part 1 — Wavepacket run (wp_rs6_v0p30)"))
    for ch in ("total", "wp", "bath"):
        cells.append(nbf.v4.new_markdown_cell(f"## Density matrix — {CH_LABEL[ch]}"))
        for plane in ("xy", "xz", "yz"):
            _, _, subt = bvn.PLANE_LABEL[plane]
            imgs = "\n\n".join(f"*{cap}*\n\n![{cap}]({rel})"
                               for _, rel, cap in wp["gifs"][(ch, plane)])
            cells.append(nbf.v4.new_markdown_cell(
                f"### {plane} plane — {subt.format(zmid=wp['zmid'])}\n\n" + imgs))
    cells.append(nbf.v4.new_markdown_cell(
        "## WP — total current density current_{x,y,z}(t)\n\n"
        f"![wp current]({wp['current']})"))
    cells.append(nbf.v4.new_markdown_cell(
        "## WP — total electronic energy vs time\n\n"
        f"![wp energy]({wp['energy']})"))

    # Classical analogue — reference the already-built rs6_v0p30 3-plane figures
    cl = "rs6_velocity_figs/rs6_v0p30"
    cells.append(nbf.v4.new_markdown_cell(
        "# Part 2 — Classical analogue (rs6_v0p30) to compare against\n\n"
        "Same r_s=6, v=0.30, but a classical Gaussian electron (no WP channel → total "
        "density only). Figures reused from the r_s=6 velocity notebook."))
    cells.append(nbf.v4.new_markdown_cell(
        "## Classical density matrix — total density, three planes\n\n" +
        "\n\n".join(
            f"### {plane} plane\n\n" + "\n\n".join(
                f"*{plane} — {k}*\n\n![{plane} {k}]({cl}/{plane}_{kf}.gif)"
                for k, kf in [("n(x·,t)", "density"), ("Δn = n(t) − n(0)", "delta0"),
                              ("Δn = n(t+dt) − n(t)", "dstep")])
            for plane in ("xy", "xz", "yz"))))
    cells.append(nbf.v4.new_markdown_cell(
        "## Classical — total current density\n\n"
        f"![cl current]({cl}/current_xyz.png)"))
    cells.append(nbf.v4.new_markdown_cell(
        "## Classical — total electronic energy vs time\n\n"
        f"![cl energy]({cl}/energy_total.png)"))

    # -------- Part 3 — DIRECT WP vs classical comparison --------
    cells.append(nbf.v4.new_markdown_cell(
        "# Part 3 — Direct WP vs classical comparison (r_s = 6, v = 0.30)\n\n"
        "Parts 1–2 show the two runs separately; here they are overlaid on the same "
        "axes / panels so the quantum wavepacket and its matched classical ghost can be "
        "read against each other directly. Same tube, same launch velocity v = 0.30, "
        "same σ (σ_WP = 0.5; the classical UPF uses the derived σ_pot = σ_WP/√2).\n\n"
        "**Density comparison — the matched pair.** The WP `density_system` = wall (24 e) "
        "**+** the WP's own |ψ|² (1 real electron) = 25 e; the classical projectile is a "
        "*pseudo-ion* that adds **no** electron density, so its `density_system` = wall "
        "(24 e). The physically comparable wall response is therefore **WP bath "
        "(= total − WP)** vs **classical total** — both the 24-electron wall reacting to "
        "the projectile, both starting from the identical wall ground state."))
    cells.append(nbf.v4.new_markdown_cell(
        "## Induced current density — WP vs classical\n\n"
        "Axial $J_z(t)$ (left) is the dominant induced current; the transverse $J_x, J_y$ "
        "(right) stay ~0 in both runs by on-axis symmetry.\n\n"
        f"![compare current]({cmp['current']})"))
    cells.append(nbf.v4.new_markdown_cell(
        "## Total electronic energy — why it compares the two runs differently\n\n"
        "$\\Delta E_\\mathrm{total}(t) = E(t) - E(0)$. **Read this carefully — the "
        "total-energy channel is NOT a like-for-like deposit here:**\n\n"
        "- **Classical:** the projectile is an *external* classical ion, so its kinetic "
        "energy is **not** part of `energy_total`. The KE it loses to the electrons "
        "appears as a **rise** in `energy_total` (~+1.2 eV) — that rise *is* the "
        "deposited energy (this is the stopping-power signal used in the classical S "
        "extraction).\n"
        "- **WP:** the wavepacket electron is *inside* the quantum system, so "
        "`energy_total` **includes the WP's own kinetic energy** and is **conserved** "
        "by the unitary (no-absorber) TDDFT evolution — ΔE ≈ 0. The energy merely "
        "redistributes internally (WP KE → wall excitation); it does not show up as a "
        "net `energy_total` change.\n\n"
        "So for the WP the meaningful energy-transfer signature is **not** the total "
        "energy but the **wake it drives in the wall density** and the **induced "
        "current** — precisely the two comparisons below, which *are* directly "
        "like-for-like. (The WP centroid also spreads across the periodic tube — "
        "σ_z² grows ~300× — so a WP 'KE loss' read off the centroid is unreliable; the "
        "density/current channels are the robust comparison.)\n\n"
        f"![compare energy]({cmp['energy']})"))
    cells.append(nbf.v4.new_markdown_cell(
        "## Induced density Δn — WP wall-response vs classical (xz plane)\n\n"
        "Rows are four matched snapshot times; columns are the WP wall response "
        "$\\Delta n(\\mathrm{bath})$, the classical $\\Delta n(\\mathrm{total})$, and their "
        "difference. Columns 1–2 **share one linear colour scale** (directly comparable, "
        "shared-colorbar rule); the difference column owns its own scale. Dashed lines "
        "mark the wall at |x| = 5, 13 Bohr. (Log is not shown here — Δn is signed.)\n\n"
        f"![compare induced density]({cmp['induced_density']})"))
    cells.append(nbf.v4.new_markdown_cell(
        "## Total density (log) — the decomposition made visible\n\n"
        "WP total (25 e) · WP bath (24 e) · classical total (24 e) at one late time, log "
        "scale, one shared colour bar. The on-axis WP blob is the **only** difference "
        "between WP-total and classical-total; once it is removed (WP bath), the WP wall "
        "density and the classical wall density are the like-for-like comparison.\n\n"
        f"![compare total density]({cmp['total_density']})"))

    nb = nbf.v4.new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})
    out = HERE / "wp_rs6_v0p30_planes.ipynb"
    nbf.write(nb, str(out))
    # verify every image reference resolves
    missing = [ref.split(")")[0] for c in cells for ref in c.source.split("](")[1:]
               if not (HERE / ref.split(")")[0]).exists()]
    print(f"wrote {out}  ({len(cells)} cells)")
    print(f"image refs missing: {missing if missing else 'none — all resolve'}")


if __name__ == "__main__":
    wp = build_wp_figs()
    cmp = build_comparison_figs()
    build_notebook(wp, cmp)
