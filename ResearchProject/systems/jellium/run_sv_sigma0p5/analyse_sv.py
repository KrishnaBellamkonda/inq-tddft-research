#!/usr/bin/env python3
"""Assemble the S(v) curve from completed run_sv_sigma0p5 velocity runs.

Scans results/<subdir>/electron_track.csv, extracts local S(v) (binned by
instantaneous v) per run, overlays the corrected Lindhard S_LR(v; sigma=0.5)
reference, and writes the money plot + a summary table. Tolerant of partial
runs: processes whatever has produced a track so far.

Usage (venv): python3 analyse_sv.py
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
REPORT = ("/local/data/public/skcb2/tddft/docs/reports/"
          "overnight-gaussian-classical-jellium/figures")
RS = 5.69
KF = E.kF_from_rs(RS)
SIGMA = 0.5


def main():
    S.apply_theme()
    tracks = sorted(glob.glob(os.path.join(RESULTS, "*", "electron_track.csv")))
    tracks = [t for t in tracks if "smoke" not in t]
    print(f"found {len(tracks)} track(s)")

    fig, ax = S.figure_one_col()
    rows = []
    cmap = S.cmap_for("sequential")
    for i, tpath in enumerate(tracks):
        sub = os.path.basename(os.path.dirname(tpath))
        sig = 0.4 if "sig0p4" in sub else 0.5
        try:
            tr = load_track(tpath, mass=1.0, axis="z")
            v, Sv = stopping_vs_v(tr, transient_bohr=3.0, window=21)
        except Exception as e:
            print(f"  {sub}: SKIP ({e})"); continue
        if v.size == 0:
            print(f"  {sub}: too short yet ({tr.s.max():.1f} Bohr)"); continue
        col = cmap(0.15 + 0.7 * i / max(len(tracks) - 1, 1))
        marker = "s" if sig == 0.4 else "o"
        ax.scatter(v, Sv, s=10, color=col, marker=marker,
                   label=f"{sub} (σ={sig})", alpha=0.7)
        rows.append((sub, sig, v.max(), v.min(), Sv.mean(), v.size))
        print(f"  {sub}: σ={sig} v∈[{v.min():.2f},{v.max():.2f}] "
              f"S∈[{Sv.min():.4f},{Sv.max():.4f}] n={v.size}")

    # Lindhard reference S_LR(v; sigma=0.5)
    vgrid = np.linspace(0.1, 3.2, 40)
    slr = np.array([E.stopping_power_sigma(v, KF, SIGMA) for v in vgrid])
    ax.plot(vgrid, slr, "-", color="k", lw=1.2,
            label="Lindhard S_LR(v;σ=0.5)", zorder=1)

    ax.axvline(KF, ls=":", color="gray", lw=0.8)
    ax.set_xlabel("v  (a.u.)")
    ax.set_ylabel("S(v)  (Ha/Bohr)")
    ax.legend(fontsize=5, loc="best")
    os.makedirs(REPORT, exist_ok=True)
    out = os.path.join(REPORT, "sv_money_plot.png")
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")

    print("\n=== summary ===")
    print(f"{'run':18} {'σ':>4} {'vmax':>6} {'vmin':>6} {'<S>':>8} {'n':>4}")
    for r in rows:
        print(f"{r[0]:18} {r[1]:>4} {r[2]:>6.2f} {r[3]:>6.2f} {r[4]:>8.4f} {r[5]:>4}")


if __name__ == "__main__":
    main()
