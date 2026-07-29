#!/usr/bin/env python3
"""SUPERSEDED prototype — do not use for report figures.

The authoritative S(v) extraction and the sv_preexisting_extraction.png figure
are produced by the energy-gain regression notebook
``docs/reports/overnight-gaussian-classical-jellium/build_sv_extraction_notebook.py``
(Method A, one point per run, single point-charge Lindhard reference).

This prototype used a local binned-by-v slope (rejected: it produced a velocity
"spread" rather than one point per run) and a per-sigma Lindhard family. Both
are obsolete. The reference here is now the SINGLE point-charge Lindhard curve
(``stopping_power_point``); the earlier "Barkas" interpretation is withdrawn
(a Z=-1 electron is identical to Z=+1 in Z^2 linear response — no Barkas from a
single sign).

Usage (venv): python3 make_sv_comparison.py
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.analysis.stopping_extract import load_track, stopping_vs_v
from inqview.analysis import lindhard_elf as E
from inqview.visualisation import style as S

RUN_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(RUN_DIR, "results")
FIGDIR = ("/local/data/public/skcb2/tddft/docs/reports/"
          "overnight-gaussian-classical-jellium/figures")

# --- jellium density + projectile constants (see REPORT.md sec 1) ---
RS = 5.69
KF = E.kF_from_rs(RS)            # = v_F (a.u.)
WP = E.omega_p(KF)              # plasmon frequency
MASS = 1.0                      # electron projectile, m_e (no fictitious mass)

# label / sigma per run subdir
RUN_LABEL = {
    "v3p0": ("v0=3.0", 0.5),
    "v2p0": ("v0=2.0", 0.5),
    "v1p3": ("v0=1.3", 0.5),
    "v0p8": ("v0=0.8", 0.5),
    "v0p6": ("v0=0.6", 0.5),
    "sig0p4_v1p0": ("v0=1.0, sigma=0.4", 0.4),
}


def lindhard_curve(sigma: float, vgrid: np.ndarray) -> np.ndarray:
    """Single point-charge Lindhard reference on a v-grid (sigma arg ignored)."""
    return np.array([E.stopping_power_point(float(v), KF) for v in vgrid])


def main():
    S.apply_theme()

    # ---- gather completed simulation tracks ----
    tracks = sorted(glob.glob(os.path.join(RESULTS, "*", "electron_track.csv")))
    runs = []  # (subdir, label, sigma, v, Sv)
    for tpath in tracks:
        sub = os.path.basename(os.path.dirname(tpath))
        if sub not in RUN_LABEL:
            continue
        summ = os.path.join(os.path.dirname(tpath), "run_summary.txt")
        if not (os.path.exists(summ) and "run_completed  = true" in open(summ).read()):
            print(f"  {sub}: not completed, skipping"); continue
        label, sigma = RUN_LABEL[sub]
        tr = load_track(tpath, mass=MASS, axis="z")
        v, Sv = stopping_vs_v(tr, transient_bohr=3.0, window=21)
        if v.size == 0:
            print(f"  {sub}: too short, skipping"); continue
        runs.append((sub, label, sigma, v, Sv))
        print(f"  {sub}: sigma={sigma} v in [{v.min():.2f},{v.max():.2f}] "
              f"<S>={Sv.mean():.4f} n={v.size}")

    # order runs by launch velocity (descending) for a stable colour assignment
    order = sorted(range(len(runs)), key=lambda i: -runs[i][3].max())
    runs = [runs[i] for i in order]

    # ---- analytic Lindhard reference curves (this density + projectile) ----
    vgrid05 = np.linspace(0.12, 3.25, 64)
    slr05 = lindhard_curve(0.5, vgrid05)
    vgrid04 = np.linspace(0.12, 1.35, 30)
    slr04 = lindhard_curve(0.4, vgrid04)
    # Bragg peak of the sigma=0.5 linear-response curve
    ipk = int(np.argmax(slr05))
    v_peak, s_peak = vgrid05[ipk], slr05[ipk]

    # ---- figure: 2 stacked panels sharing the v axis ----
    fig, (ax, axr) = S.plt.subplots(
        2, 1, sharex=True, figsize=(6.6, 6.4),
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.08},
    )

    cmap = S.plt.get_cmap(S.cmap_for("sequential"))
    n = len(runs)
    for i, (sub, label, sigma, v, Sv) in enumerate(runs):
        col = cmap(0.12 + 0.74 * i / max(n - 1, 1))
        marker = "s" if sigma == 0.4 else "o"
        ax.scatter(v, Sv, s=11, color=col, marker=marker, alpha=0.75,
                   edgecolors="none", label=f"{label}")
        # ratio vs the analytic curve at THIS run's sigma
        if sigma == 0.4:
            slr_at = np.interp(v, vgrid04, slr04)
        else:
            slr_at = np.interp(v, vgrid05, slr05)
        ratio = Sv / slr_at
        axr.scatter(v, ratio, s=11, color=col, marker=marker, alpha=0.75,
                    edgecolors="none")

    # single point-charge Lindhard reference (one curve, all runs)
    ax.plot(vgrid05, slr05, "-", color="k", lw=1.6, zorder=3,
            label="Lindhard S$_{LR}$ (point charge)")
    ax.scatter([v_peak], [s_peak], marker="*", s=90, color="k", zorder=4)
    ax.annotate(f"LR Bragg peak\nv={v_peak:.2f}, S={s_peak:.3f}",
                xy=(v_peak, s_peak), xytext=(v_peak + 0.35, s_peak),
                fontsize=6, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.4"))

    # reference markers
    for a in (ax, axr):
        a.axvline(KF, ls=":", color="gray", lw=0.8)
    ax.text(KF, ax.get_ylim()[1] * 0.02, " v$_F$=k$_F$", color="gray",
            fontsize=6, rotation=90, va="bottom", ha="left")
    axr.axhline(1.0, ls="-", color="0.3", lw=0.9)
    axr.text(2.9, 0.62, "below point-charge LR\n(finite-width suppression)",
             fontsize=5.5, color="0.4", ha="center", va="center")

    ax.set_ylabel("S(v)  (Ha/Bohr)")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=6, loc="upper right", framealpha=0.9, ncol=1)
    ax.set_title(
        f"rt-TDDFT electron vs Lindhard LR  (jellium r$_s$={RS}, "
        f"k$_F$={KF:.3f}, $\\omega_p$={WP:.3f})", fontsize=7.5)

    axr.set_ylabel("S$_{sim}$ / S$_{LR}$")
    axr.set_xlabel("v  (a.u.)")
    axr.set_xlim(0, 3.35)
    axr.set_ylim(0.4, 1.7)

    out = os.path.join(FIGDIR, "sv_preexisting_extraction.png")
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches=None)
    print(f"\nwrote {out}")

    # ---- print a compact comparison table (median ratio per run) ----
    print(f"\n{'run':16} {'sigma':>5} {'v_mid':>6} {'S_sim':>8} {'S_LR':>8} {'ratio':>6}")
    for sub, label, sigma, v, Sv in runs:
        grid = (vgrid04, slr04) if sigma == 0.4 else (vgrid05, slr05)
        slr_at = np.interp(v, *grid)
        ratio = Sv / slr_at
        vm = np.median(v)
        print(f"{sub:16} {sigma:>5} {vm:>6.2f} {np.median(Sv):>8.4f} "
              f"{np.median(slr_at):>8.4f} {np.median(ratio):>6.2f}")


if __name__ == "__main__":
    main()
