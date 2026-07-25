#!/usr/bin/env python3
"""qsp_phase5 — CROSS-RUN comparison figures for the study notebook.

Lays the six WP runs (5 new + reused 54 eV) side by side so the velocity trend —
and the high-v grid-aliasing artifact — is visible at a glance:
  • cmp_energy.png        — E_total(t)−E_GS and E_kinetic(t), all runs overlaid
  • cmp_norm.png          — WP orbital norm(t) (absorption), all runs
  • cmp_momentum_kz.png    — n_wp(k_z) at t=0 vs the grid Nyquist k=π/h=6.28 (THE aliasing evidence)
  • cmp_centroid_sigma.png — ⟨z⟩(t) transit + σ_z(t) spreading, all runs
  • figs/cmp_density_*eV.gif — the xz total-density GIFs gathered + downsized

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_phase5_comparisons.py
"""
from __future__ import annotations
import os, glob, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

HA   = 27.211386
E_GS = -70.22568216820937
HBAR_H = 0.5                       # spacing (Bohr)
K_NYQ = np.pi / HBAR_H             # = 6.2832 a.u.
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)
LJ = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
try:
    from inqview.visualisation import style as STYLE; STYLE.apply_theme()
except Exception as exc:  # noqa: BLE001
    print(f"[cmp] theme unavailable ({exc})")

# (tag, v, E_eV, results_dir, total_density_gif, verdict)
RUNS = [
    ("v1.3", 1.3,  23, f"{LJ}/scripts/qsp_phase5/wp/results/p5_wp_v1p3",
     f"{HERE}/p5_wp_v1p3_run_notebook_figs/p5_wp_v1p3_total_density.gif", "clean"),
    ("v2.0", 2.0,  54, f"{LJ}/scripts/qsp_phase4/wp/results/p4_wp",
     f"{HERE}/../qsp_phase4/p4wp_run_notebook_figs/p4_wp_total_density.gif", "clean"),
    ("v3.0", 3.0, 122, f"{LJ}/scripts/qsp_phase5/wp/results/p5_wp_v3p0",
     f"{HERE}/p5_wp_v3p0_run_notebook_figs/p5_wp_v3p0_total_density.gif", "clean"),
    ("v4.0", 4.0, 218, f"{LJ}/scripts/qsp_phase5/wp/results/p5_wp_v4p0",
     f"{HERE}/p5_wp_v4p0_run_notebook_figs/p5_wp_v4p0_total_density.gif", "borderline"),
    ("v5.0", 5.0, 340, f"{LJ}/scripts/qsp_phase5/wp/results/p5_wp_v5p0",
     f"{HERE}/p5_wp_v5p0_run_notebook_figs/p5_wp_v5p0_total_density.gif", "aliased_rerun"),
    ("v6.0", 6.0, 490, f"{LJ}/scripts/qsp_phase5/wp/results/p5_wp_v6p0",
     f"{HERE}/p5_wp_v6p0_run_notebook_figs/p5_wp_v6p0_total_density.gif", "aliased"),
]
COLORS = cm.viridis(np.linspace(0.05, 0.92, len(RUNS)))
LS = {"clean": "-", "borderline": "--", "aliased": ":", "aliased_rerun": ":"}


def _label(r):
    # NB: the SHOWN v5 run is the original h=0.5 (aliased) one — the coarse-grid
    # diagnostic; it was re-run on a finer h=0.35 grid for the S(E) point.
    suf = {"clean": "", "borderline": " [borderline]", "aliased": " [ALIASED]",
           "aliased_rerun": " [h=0.5 aliased → re-run @0.35]"}[r[5]]
    return f"{r[2]} eV ({r[0]}){suf}"


def _obs(rdir):
    return pd.read_csv(os.path.join(rdir, "raw", "observables", "observables.csv"))


def _stats(rdir):
    return pd.read_csv(os.path.join(rdir, "raw", "observables", "wp_real_space_stats.csv"), comment="#")


def energy_fig():
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for r, c in zip(RUNS, COLORS):
        try:
            o = _obs(r[3]); t = o["time_au"].values
            a.plot(t, (o["energy_total"].values - E_GS) * HA, LS[r[5]], color=c, lw=1.7, label=_label(r))
            b.plot(t, o["energy_kinetic"].values * HA, LS[r[5]], color=c, lw=1.7, label=_label(r))
        except Exception as exc:  # noqa: BLE001
            print(f"[cmp energy] {r[0]}: {exc}")
    a.set_xlabel("time (a.u.)"); a.set_ylabel("E_total(t) − E_GS  (eV)")
    a.set_title("Retained energy  E_total(t) − E_GS  (= S·L_z at convergence)", fontsize=9)
    a.grid(alpha=.25); a.legend(fontsize=7, frameon=False, ncol=2)
    b.set_xlabel("time (a.u.)"); b.set_ylabel("E_kinetic(t)  (eV)")
    b.set_title("Total kinetic energy — WP KE entering/leaving the box", fontsize=9)
    b.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/cmp_energy.png", dpi=170); plt.close(fig)
    print("[cmp] cmp_energy.png")


def norm_fig():
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for r, c in zip(RUNS, COLORS):
        try:
            s = _stats(r[3])
            ax.plot(s["time_au"].values, s["norm_check"].values, LS[r[5]], color=c, lw=1.8, label=_label(r))
        except Exception as exc:  # noqa: BLE001
            print(f"[cmp norm] {r[0]}: {exc}")
    ax.axhline(0.02, ls="--", color="0.6", lw=.8, label="convergence gate 0.02")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("WP orbital norm")
    ax.set_title("WP absorption: norm(t)→0 as the packet exits into the CAP", fontsize=9.5)
    ax.grid(alpha=.25); ax.legend(fontsize=7.5, frameon=False, ncol=2)
    fig.tight_layout(); fig.savefig(f"{FIGS}/cmp_norm.png", dpi=170); plt.close(fig)
    print("[cmp] cmp_norm.png")


def momentum_fig():
    """n_wp(k_z) at t=0 vs the grid Nyquist — the aliasing evidence."""
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    for r, c in zip(RUNS, COLORS):
        try:
            m = pd.read_csv(os.path.join(r[3], "raw", "observables", "momentum_distribution.csv"), comment="#")
            m0 = m[m["step"] == m["step"].min()]
            k, nwp = m0["k_bohr_inv"].values, m0["n_wp"].values
            ax.plot(k, nwp / max(nwp.max(), 1e-30), LS[r[5]], color=c, lw=1.8, label=_label(r))
        except Exception as exc:  # noqa: BLE001
            print(f"[cmp momentum] {r[0]}: {exc}")
    # Nyquist wavevector shown as a labelled reference line (legend, not on-plot text).
    ax.axvline(K_NYQ, color="crimson", lw=1.6, ls="-", label=f"grid Nyquist k=π/h={K_NYQ:.2f}")
    ax.axvspan(K_NYQ, K_NYQ + 3, color="crimson", alpha=0.10)
    ax.set_xlabel("k_z  (a.u.)"); ax.set_ylabel("n_wp(k_z) at t=0  (peak-normalised)")
    ax.set_title("WP momentum  n_wp(k_z) at t=0  vs grid Nyquist", fontsize=9.5)
    ax.set_xlim(0, 6.5); ax.set_ylim(0, 1.15)
    ax.grid(alpha=.25); ax.legend(fontsize=7.5, frameon=False, loc="center left")
    fig.tight_layout(); fig.savefig(f"{FIGS}/cmp_momentum_kz.png", dpi=175); plt.close(fig)
    print("[cmp] cmp_momentum_kz.png")


def centroid_fig():
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for r, c in zip(RUNS, COLORS):
        try:
            s = _stats(r[3]); t = s["time_au"].values
            zc = s["z_mean"].values if "z_mean" in s else None
            sz = np.sqrt(s["sigma_z2"].values) if "sigma_z2" in s else None
            if zc is not None: a.plot(t, zc, LS[r[5]], color=c, lw=1.7, label=_label(r))
            if sz is not None: b.plot(t, sz, LS[r[5]], color=c, lw=1.7, label=_label(r))
        except Exception as exc:  # noqa: BLE001
            print(f"[cmp centroid] {r[0]}: {exc}")
    for zz, lab in [(12.5, "slab face"), (-12.5, None), (35, "CAP inner"), (-35, None)]:
        a.axhline(zz, ls="--", lw=.7, color="0.6")
    a.set_xlabel("time (a.u.)"); a.set_ylabel("⟨z⟩ (Bohr)")
    a.set_title("WP centroid ⟨z⟩(t) — slope = velocity; faces ±12.5, CAP ±35", fontsize=8.8)
    a.grid(alpha=.25); a.legend(fontsize=7, frameon=False, ncol=2)
    b.set_xlabel("time (a.u.)"); b.set_ylabel("σ_z (Bohr)")
    b.set_title("WP spreading σ_z(t)", fontsize=9); b.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/cmp_centroid_sigma.png", dpi=170); plt.close(fig)
    print("[cmp] cmp_centroid_sigma.png")


def gather_gifs(target_w=300):
    from PIL import Image, ImageSequence
    out = []
    for r in RUNS:
        src = r[4]
        if not os.path.exists(src):
            print(f"[cmp gif] MISSING {r[0]}: {src}"); continue
        dst = f"{FIGS}/cmp_density_{r[2]}eV.gif"
        try:
            im = Image.open(src); fr = [f.copy() for f in ImageSequence.Iterator(im)]
            du = [f.info.get("duration", 100) for f in fr]; w, h = fr[0].size
            if w > target_w:
                sc = target_w / w
                fr = [f.convert("RGB").resize((int(w*sc), int(h*sc)), Image.LANCZOS) for f in fr]
            fr = [f.convert("RGB").quantize(colors=64, method=Image.FASTOCTREE) for f in fr]
            fr[0].save(dst, save_all=True, append_images=fr[1:], loop=0, duration=du, optimize=True, disposal=2)
            out.append((r, dst)); print(f"[cmp gif] {r[0]} -> {os.path.basename(dst)}")
        except Exception as exc:  # noqa: BLE001
            print(f"[cmp gif] {r[0]} failed ({exc}); copying raw"); shutil.copy(src, dst); out.append((r, dst))
    return out


if __name__ == "__main__":
    energy_fig(); norm_fig(); momentum_fig(); centroid_fig(); gather_gifs()
    print("[cmp] done")
