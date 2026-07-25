#!/usr/bin/env python3
"""S(v) / S(E) placement plots for the 26-6-26 meeting.

Plot A: linear-response (point-charge Lindhard) curve + classical sigma_WP=0.5 BULK
        sweep points + the QUANTUM (WP) 100 eV point.
Plot B: Plot A + the classical-through-slab point (p2_classical).

The classical bulk S(v) is extracted with the SAME Method-A kernel as the
presentation's Section-1 figure (fixed 20%-of-time transient cut + free-intercept
slope fit). The Lindhard reference uses kF=0.337 (r_s=5.69), matching that figure.

WP point: from quantum_stopping_ledger (energy method). It is an UPPER BOUND at
tau=40 (the packet is not fully absorbed and is zero-point-contaminated), so it is
drawn with a downward arrow and a guide-band toward the comparable-to-LR zone where
the *converged* value is physically expected. Classical-slab point: slab-centre
KE-min estimate (user's choice; not trusted to be quantitatively right).

Numbers rounded to 2 s.f. in annotations per .claude/rules/number-rounding.md.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/local/data/public/skcb2/tddft")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "inq-stack/python"))
sys.path.insert(0, str(HERE))

from inqview.visualisation import style          # noqa: E402
from inqview.analysis import lindhard_elf as LR   # noqa: E402
import quantum_stopping_ledger as L               # noqa: E402

_spk = ROOT / ".claude/skills/stopping-power-extraction/stopping_power.py"
_spec = importlib.util.spec_from_file_location("spk", _spk)
spk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(spk)

style.apply_theme()
HA_TO_EV = L.HA_TO_EV
FIGDIR = HERE / "sv_plots_26-6-26_figs"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Linear-response reference at the SLAB density (the medium the WP traverses):
# n0 = 82 e / (50*50*25) = 1.312e-3, r_s=5.667, kF=0.3387 (slab_n82_L50x50x70.hpp).
# This is 0.41% from the bulk S(v) density (162/50^3, kF=0.3373) — matched by design,
# so one curve serves both the bulk-classical points and the quantum-slab point.
# Lindhard 1954 / Lindhard & Winther 1964 (docs/sources/stopping-power-formulae.md);
# finite-slab applicability is velocity-gated (Quijada 2007).
N_DENS = 82 / (50.0 * 50.0 * 25.0)            # slab density 1.312e-3 a0^-3
KF = (3 * np.pi**2 * N_DENS) ** (1.0 / 3.0)   # 0.3387 Bohr^-1 (slab)
JB = ROOT / "ResearchProject/systems/jellium"

# classical BULK sigma_WP=0.5 series (run dir labels sigma0p35 = charge std)
RUN_DIR = "run_classical_n162_L50_sv_sigma0p35"
PSP_TAG = "sigma0p35"
VTAGS = [(0.2, "v0p2"), (0.6, "v0p6"), (0.8, "v0p8"), (1.0, "v1p0"),
         (1.3, "v1p3"), (2.0, "v2p0"), (3.0, "v3p0")]


def _psp_ok(vtag):
    rs = JB / RUN_DIR / "results" / vtag / "run_summary.txt"
    if not rs.exists():
        return False
    for line in rs.read_text().splitlines():
        if line.strip().startswith("psp"):
            return f"electron_gaussian_{PSP_TAG}.upf" in line
    return False


def _extract_S(track_csv):
    df = pd.read_csv(track_csv).drop_duplicates(subset="step", keep="last")
    df = df.dropna(subset=["vx", "vy", "vz", "z", "time_au"]).reset_index(drop=True)
    if len(df) < 8:
        return None
    t = df["time_au"].to_numpy(); z = df["z"].to_numpy()
    KE = 0.5 * (df["vx"]**2 + df["vy"]**2 + df["vz"]**2).to_numpy() * HA_TO_EV
    fit = spk.fixed_time_fraction(t, z, KE[0] - KE, frac=0.20)
    if fit.get("status") != "ok":
        return None
    return dict(v=abs(df["vz"].to_numpy()[0]), S=fit["S"], se=fit["se"])


def classical_bulk_sv():
    rows = []
    for v0, vtag in VTAGS:
        trk = JB / RUN_DIR / "results" / vtag / "raw" / "observables" / "electron_track.csv"
        if not trk.exists() or not _psp_ok(vtag):
            continue
        r = _extract_S(trk)
        if r:
            rows.append(r)
            print(f"  classical bulk sigma_WP=0.5  {vtag}: v={r['v']:.3f}  S={r['S']:.3f} eV/Bohr")
    return pd.DataFrame(rows).sort_values("v").reset_index(drop=True)


def lindhard_curve(v_lo, v_hi):
    cache = FIGDIR / f"_lr_kF{KF:.4f}_{v_lo:.3f}_{v_hi:.3f}.npz"
    v = np.logspace(np.log10(v_lo), np.log10(v_hi), 50)
    if cache.exists():
        z = np.load(cache); return z["v"], z["S"]
    S = np.array([LR.stopping_power_point(float(vv), KF) for vv in v]) * HA_TO_EV
    np.savez(cache, v=v, S=S)
    return v, S


def build():
    cls = classical_bulk_sv()
    # both r_s=5.67 WP runs (same density as this S(v) plot): p2 (tau=40) and the
    # big-box converged p3 (tau=100). Both full-ledger UPPER BOUNDS.
    wp_p2 = L.compute_wp_ledger("p2")
    wp_p3 = L.compute_wp_ledger("p3")
    cl_slab = L.compute_classical_slab("p2")  # p3 classical is anomalous (trapped)
    v_pt = L.V_POINT_AU
    S_wp_p2 = wp_p2["S_wp_ev_per_bohr"]       # 2.7
    S_wp_p3 = wp_p3["S_wp_ev_per_bohr"]       # 2.4
    S_slab = cl_slab["S_center_ev_per_bohr"]  # p2 slab-centre 3.0 (user choice)

    v_all = np.concatenate([cls["v"].to_numpy(), [v_pt]])
    v_lo, v_hi = min(v_all.min() * 0.8, 0.18), max(v_all.max() * 1.15, 3.5)
    v_grid, S_lr = lindhard_curve(v_lo, v_hi)
    # LR value at the WP velocity (for the guide-band bottom = "comparable" zone)
    S_lr_at_pt = float(np.interp(v_pt, v_grid, S_lr))

    def _panel(ax, x_of_v, xlabel, with_slab):
        ax.plot(x_of_v(v_grid), S_lr, "-", color="k", lw=1.8, zorder=2,
                label=r"linear response (slab $r_s$=5.67)")
        ax.errorbar(x_of_v(cls["v"].to_numpy()), cls["S"].to_numpy(),
                    yerr=cls["se"].to_numpy(), fmt="o", color="C0", ms=5.5,
                    capsize=2.5, mec="k", mew=0.4, zorder=4,
                    label=r"classical $\sigma_{\rm WP}=0.5$ (bulk)")
        # quantum points — UPPER BOUNDS (both r_s=5.67): p2 (tau=40) and p3 (tau=100).
        # They barely differ -> the offset above LR is the zero-point, not convergence.
        xpt = x_of_v(np.array([v_pt]))[0]
        ax.scatter([xpt], [S_wp_p2], marker="v", s=80, color="C3", ec="k", lw=0.6,
                   zorder=6, label=rf"quantum WP $\tau$=40 (upper bound, {S_wp_p2:.1f})")
        ax.scatter([xpt], [S_wp_p3], marker="^", s=80, color="C1", ec="k", lw=0.6,
                   zorder=6, label=rf"quantum WP $\tau$=100 (upper bound, {S_wp_p3:.1f})")
        if with_slab:
            ax.scatter([xpt], [S_slab], marker="s", s=70, color="C4", ec="k",
                       lw=0.5, zorder=6,
                       label="classical through slab (not trusted)")
        ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("stopping power  S  (eV/Bohr)")
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=6.8, frameon=False, loc="upper left")

    def make(with_slab, tag, title):
        # S(v)
        fig, ax = style.figure_one_col()
        _panel(ax, lambda v: v, "velocity  v  (a.u.)", with_slab)
        ax.set_xlim(v_lo, v_hi)
        ax.set_title(title + r" — $S(v)$")
        style.save_presentation(fig, FIGDIR / f"fig_sv_{tag}.png")
        plt.close(fig)
        # S(E)
        fig, ax = style.figure_one_col()
        _panel(ax, lambda v: 0.5 * v**2 * HA_TO_EV, "projectile energy  E  (eV)", with_slab)
        ax.set_xlim(0.5 * v_lo**2 * HA_TO_EV, 0.5 * v_hi**2 * HA_TO_EV)
        ax.set_title(title + r" — $S(E)$")
        style.save_presentation(fig, FIGDIR / f"fig_se_{tag}.png")
        plt.close(fig)
        print(f"  wrote fig_sv_{tag}.png / fig_se_{tag}.png")

    make(False, "A_quantum", "Quantum WP points (tau=40 & 100)")
    make(True, "B_quantum_plus_slab", "Quantum WP + classical slab")
    print(f"\n  WP points: v={v_pt:.3f}  S_p2={S_wp_p2:.2f}  S_p3={S_wp_p3:.2f} "
          f"(both upper bounds); LR@v={S_lr_at_pt:.3f}; p2 classical-slab={S_slab:.2f}")
    print(f"  figures in {FIGDIR}")


if __name__ == "__main__":
    build()
