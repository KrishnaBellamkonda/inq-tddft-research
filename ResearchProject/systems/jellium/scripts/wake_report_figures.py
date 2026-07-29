#!/usr/bin/env python3
"""Report-standard STATIC wake-difference panels (applications.report1 style).

Remakes the induced-wake WP-vs-classical comparison as publication figures per
the report-figures skill: apply_style() (usetex/Computer Modern, 600 DPI), NO
titles/suptitles (captions go in LaTeX), INDIVIDUAL plots (one PNG per panel,
composed in LaTeX minipages), RdBu_r diverging maps, shared colour scale for the
directly-compared WP & classical panels, own scale for the difference.

Per case (a WP run + the matched classical run) at a representative time t*:
  fig_wake_<tag>_2d_wp.png        Δn_system (WP), RdBu_r, shared clim
  fig_wake_<tag>_2d_classical.png Δn_system (classical), SAME clim
  fig_wake_<tag>_2d_diff.png      Δn_WP − Δn_classical, own clim
  fig_wake_<tag>_1d.png           z-profile: WP, classical, difference + centroid

Three steps (see inqview.pipeline.wake): n_system=n_total−n_wp (exact step) →
Δn=n_system(t)−n_system(t0) → cross-system Δn_WP−Δn_classical at matched t.

Output dir: docs/presentations/storyline/tasks/batch2_figures/report_standard/
Usage: wake_report_figures.py <case>   (case = sigma1 | sigma0p5 | ... | E20 ...)
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, SymLogNorm
from inqview.pipeline import wake
from applications.report1 import apply_style, panel_label, palette_regime3, references
from applications.report1._shared_style import ONE_COL_IN, STYLE_CONFIG

apply_style()
JB = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
OUT = Path("/local/data/public/skcb2/tddft/docs/presentations/storyline/tasks/batch2_figures/report_standard")
OUT.mkdir(parents=True, exist_ok=True)
CL_E100 = f"{JB}/run_classical_n162_L50_E100_v2"
SIGMA = {"0p5": f"{JB}/run_wp_n162_L50_E100_sigma0p5_wf",
         "1":   f"{JB}/run_wp_n162_L50_E100_sigma1_v2",
         "3":   f"{JB}/run_wp_n162_L50_E100_sigma3_wf",
         "8":   f"{JB}/run_wp_n162_L50_E100_sigma8_wf"}
ENERGY = {"20": (f"{JB}/run_wp_n162_L50_E20_sigma1_v2", f"{JB}/run_classical_n162_L50_E20"),
          "25": (f"{JB}/run_wp_n162_L50_E25_sigma1_v2", f"{JB}/run_classical_n162_L50_E25"),
          "50": (f"{JB}/run_wp_n162_L50_E50_sigma1_v2", f"{JB}/run_classical_n162_L50_E50_v2"),
          "100": (f"{JB}/run_wp_n162_L50_E100_sigma1_v2", CL_E100),
          "300": (f"{JB}/run_wp_n162_L50_E300_sigma1_v2", f"{JB}/run_classical_n162_L50_E300_v2")}
CLIP = 99.5


def induced(run, t, t0frames):
    """Δn 2D slab + 1D line at the WP-exact frame nearest t (classical: any)."""
    if wake.has_wp(run):
        ts = np.array(wake.wp_frame_times(run)); t = ts[np.argmin(np.abs(ts - t))]
    nsys, o, sp, ta, _ = wake.bath_volume(run, t)
    ny = nsys.shape[1]
    slab = nsys[:, ny // 2, :] - t0frames[run]["slab"]
    line = nsys.sum(axis=(1, 2)) * sp[0] * sp[1] - t0frames[run]["line"]
    x = o[0] + sp[0] * np.arange(nsys.shape[2]); z = o[2] + sp[2] * np.arange(nsys.shape[0])
    return x, z, slab, line, ta


def t0_cache(run):
    nsys, o, sp, ta, _ = wake.bath_volume(run, 0.0)
    ny = nsys.shape[1]
    return {"slab": nsys[:, ny // 2, :], "line": nsys.sum(axis=(1, 2)) * sp[0] * sp[1]}


def imshow_panel(slab, x, z, vmax, tag, label, norm=None):
    fig, ax = plt.subplots(figsize=ONE_COL_IN)
    ext = [z[0], z[-1], x[0], x[-1]]
    im = ax.imshow(slab.T, origin="lower", extent=ext, aspect="auto", cmap="RdBu_r",
                   norm=norm, vmin=None if norm else -vmax, vmax=None if norm else vmax)
    ax.set_xlabel(r"$z$ (Bohr)"); ax.set_ylabel(r"$x$ (Bohr)")
    panel_label(ax, label)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(r"$\Delta n_{\mathrm{system}}$ (e/Bohr$^3$)")
    fig.subplots_adjust(left=0.16, right=0.86, bottom=0.16, top=0.97)
    fig.savefig(OUT / f"{tag}.png", dpi=STYLE_CONFIG["save_dpi"])
    plt.close(fig)
    print(f"  wrote {tag}.png")


def make_case(wp_run, cl_run, tag, tfrac=0.6):
    t0c = {wp_run: t0_cache(wp_run), cl_run: t0_cache(cl_run)}
    ts = np.array(wake.wp_frame_times(wp_run)); tstar = ts[int(tfrac * (len(ts) - 1))]
    xw, zw, sw, lw, taw = induced(wp_run, tstar, t0c)
    xc, zc, sc, lc, tac = induced(cl_run, tstar, t0c)
    cent = wake.wp_centroid_z(wp_run, tstar)
    print(f"[{tag}] t*={taw:.2f} a.u.; centroid={cent:.1f} Bohr")
    # shared clim for WP & classical (directly compared); own for diff
    vmax = max(np.percentile(np.abs(sw), CLIP), np.percentile(np.abs(sc), CLIP))
    dvmax = np.percentile(np.abs(sw - sc), CLIP)
    imshow_panel(sw, xw, zw, vmax, f"fig_wake_{tag}_2d_wp", "(a)")
    imshow_panel(sc, xc, zc, vmax, f"fig_wake_{tag}_2d_classical", "(b)")
    imshow_panel(sw - sc, xw, zw, dvmax, f"fig_wake_{tag}_2d_diff", "(c)")
    # 1D z-profile (individual panel)
    fig, ax = plt.subplots(figsize=ONE_COL_IN)
    ax.plot(zw, lw, color=palette_regime3[0], label="WP")
    ax.plot(zc, lc, color=palette_regime3[1], label="classical")
    ax.plot(zw, lw - lc, color=palette_regime3[2], lw=1.0, label="WP $-$ classical")
    ax.axhline(0, **references["asymptote"])
    if cent is not None and np.isfinite(cent):
        ax.axvline(cent, color=references["annotation"]["color"], ls=":", lw=0.9)
    ax.set_xlabel(r"$z$ (Bohr)"); ax.set_ylabel(r"induced $\Delta n(z)$ (e/Bohr)")
    ax.legend(loc="lower right")
    panel_label(ax, "(d)")
    fix = dict(left=0.185, right=0.97, bottom=0.16, top=0.97); fig.subplots_adjust(**fix)
    fig.savefig(OUT / f"fig_wake_{tag}_1d.png", dpi=STYLE_CONFIG["save_dpi"])
    plt.close(fig); print(f"  wrote fig_wake_{tag}_1d.png")


if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "sigma1"
    if case == "all":
        for s in SIGMA:
            make_case(SIGMA[s], CL_E100, f"sigma{s}_E100")
        for e, (wp, cl) in ENERGY.items():
            make_case(wp, cl, f"E{e}_sigma1")
    elif case.startswith("sigma"):
        s = case[len("sigma"):]
        make_case(SIGMA[s], CL_E100, f"sigma{s}_E100")
    elif case.startswith("E"):
        e = case[1:]; wp, cl = ENERGY[e]
        make_case(wp, cl, f"E{e}_sigma1")
