#!/usr/bin/env python3
"""build_energy_and_stopping_notebook.py — energy overlays + stopping-power walkthrough.

User-specified (2026-07-09). ONE study notebook for the annular-tube runs with:

  PART A  three total-energy-vs-time plots, one per r_s ∈ {6,4,2}; each overlays the
          three projectile velocities v ∈ {0.15,0.30,0.45} as traces.
  PART B  three total-energy-vs-time plots, one per velocity; each overlays the three
          r_s values as traces (same v, different r_s).
  PART C  the ΔE = E(t) − E(0) version of every plot in A and B (6 more).
  PART D  for the r_s = 6, v0 = 0.45 run: a STEP-BY-STEP walkthrough of how the
          stopping power was defined/measured — every step rendered as its own plot,
          one after the other, faithfully reproducing the pipeline in `per_run.py`
          (`stopping_analysis`) which calls the stopping-power-extraction skill kernels
          (`stopping_power.py`; Correa 2018 Method A, light-projectile initial-drag).

Data (per run):  observables.csv  (time_au, energy_total, current_*, ...)
                 electron_track.csv (z, vz, ke_ion_ha)
                 electron_number.csv (N_total)  — for the N-conservation guard.

Writes `energy_and_stopping.ipynb` + figures under `energy_stopping_figs/` beside it
(file-placement: run-tied analysis lives in hypotheses/annular_sv/).

Run (venv + stack + skill on path):
    PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_energy_and_stopping_notebook.py
"""
from __future__ import annotations
import glob
import os
import sys
from pathlib import Path

import numpy as np

STACK = "/local/data/public/skcb2/tddft/inq-stack/python"
SKILL_SP = "/home/raid/skcb2/skcb2/tddft/.claude/skills/stopping-power-extraction"
for p in (STACK, SKILL_SP):
    if p not in sys.path:
        sys.path.insert(0, p)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nbformat as nbf  # noqa: E402
import pandas as pd  # noqa: E402

import stopping_power as sp  # noqa: E402  the skill kernels (free_fit, kinetic_channel, …)
from inqview.visualisation import style  # noqa: E402
try:
    style.apply()
except Exception:
    pass

HA_EV = 27.21138625

SYS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium")
RUNROOT = SYS / "annular_sv"
HERE = SYS / "hypotheses" / "annular_sv"
FIGROOT = HERE / "energy_stopping_figs"
STEPDIR = FIGROOT / "stopping_steps"

RS_LIST = [6, 4, 2]
VELS = [0.15, 0.30, 0.45]
# density word + N per r_s (from the run configs)
RS_META = {6: dict(word="dilute", N=24, L_z=48, launch=-23),
           4: dict(word="intermediate", N=48, L_z=28, launch=-13),
           2: dict(word="dense", N=136, L_z=10, launch=-4)}
# consistent colour per member of a family
VCOL = {0.15: "C0", 0.30: "C1", 0.45: "C3"}
RSCOL = {6: "C0", 4: "C2", 2: "C3"}


def run_name(rs, v):
    return f"rs{rs}_v{v:.2f}".replace(".", "p")


def _find(rs, v, name):
    hits = glob.glob(str(RUNROOT / run_name(rs, v) / "**" / name), recursive=True)
    return hits[0] if hits else None


def load_obs(rs, v):
    f = _find(rs, v, "observables.csv")
    return pd.read_csv(f).drop_duplicates("step").sort_values("time_au") if f else None


def load_trk(rs, v):
    f = _find(rs, v, "electron_track.csv")
    return pd.read_csv(f).drop_duplicates("step").sort_values("time_au") if f else None


# =========================================================== PART A / B / C plots
def energy_by_rs(rs, delta: bool):
    """One figure: E_total(t) (or ΔE) for a fixed r_s, traces = velocity."""
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for v in VELS:
        O = load_obs(rs, v)
        if O is None:
            continue
        e = O["energy_total"].to_numpy()
        y = e - e[0] if delta else e
        lab = f"v = {v:.2f}"
        if delta:
            lab += f"  (peak ΔE = {y.max():.3f} Ha)"
        ax.plot(O["time_au"], y, VCOL[v], lw=1.6, label=lab)
    ax.set_xlabel("time (a.u.)")
    if delta:
        ax.set_ylabel(r"$\Delta E_\mathrm{total} = E(t)-E(0)$  (Ha)")
        ax.set_title(f"r_s = {rs} ({RS_META[rs]['word']}) tube: energy deposited vs time")
        ax.axhline(0, color="0.6", lw=0.7)
    else:
        ax.set_ylabel(r"$E_\mathrm{total}$  (Ha)")
        ax.set_title(f"r_s = {rs} ({RS_META[rs]['word']}) tube: total electronic "
                     "energy vs time")
    ax.legend(title="projectile velocity"); ax.grid(alpha=.25)
    tag = "delta" if delta else "raw"
    f = FIGROOT / f"energy_{tag}_rs{rs}.png"
    fig.tight_layout(); fig.savefig(f, dpi=160); plt.close(fig)
    return os.path.relpath(f, HERE)


def energy_by_velocity(v, delta: bool):
    """One figure: E_total(t) (or ΔE) for a fixed velocity, traces = r_s."""
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for rs in RS_LIST:
        O = load_obs(rs, v)
        if O is None:
            continue
        e = O["energy_total"].to_numpy()
        y = e - e[0] if delta else e
        lab = f"r_s = {rs} ({RS_META[rs]['word']})"
        if delta:
            lab += f"  (peak ΔE = {y.max():.3f} Ha)"
        else:
            lab += f"  (E(0) = {e[0]:.1f} Ha)"
        ax.plot(O["time_au"], y, RSCOL[rs], lw=1.6, label=lab)
    ax.set_xlabel("time (a.u.)")
    if delta:
        ax.set_ylabel(r"$\Delta E_\mathrm{total} = E(t)-E(0)$  (Ha)")
        ax.set_title(f"v = {v:.2f} a.u.: energy deposited vs time, by wall density r_s")
        ax.axhline(0, color="0.6", lw=0.7)
    else:
        ax.set_ylabel(r"$E_\mathrm{total}$  (Ha)")
        ax.set_title(f"v = {v:.2f} a.u.: total electronic energy vs time, by r_s\n"
                     "(absolute baselines differ by ~200 Ha — see the ΔE version)",
                     fontsize=10)
    ax.legend(title="wall density"); ax.grid(alpha=.25)
    tag = "delta" if delta else "raw"
    f = FIGROOT / (f"energy_{tag}_v{v:.2f}".replace(".", "p") + ".png")
    fig.tight_layout(); fig.savefig(f, dpi=160); plt.close(fig)
    return os.path.relpath(f, HERE)


# ===================================================== PART D — stopping-power steps
def stopping_steps(rs=6, v0=0.45, vfrac=0.85):
    """Reproduce `per_run.stopping_analysis` step-by-step, one figure per step, for the
    r_s=6 v0=0.45 run. Returns an ordered list of (rel_path, title, method_md)."""
    STEPDIR.mkdir(parents=True, exist_ok=True)
    O = load_obs(rs, v0)
    T = load_trk(rs, v0)
    numf = _find(rs, v0, "electron_number.csv")
    N = pd.read_csv(numf)["N_total"].to_numpy() if numf else None

    tE = O["time_au"].to_numpy()
    e_abs = O["energy_total"].to_numpy()
    E = e_abs - e_abs[0]                                  # electronic deposit ΔE_total(t)
    tT = T["time_au"].to_numpy()
    z = T["z"].to_numpy(); z0 = z[0]
    s_tr = np.abs(z - z0)                                 # path length s = |z − z0|
    vz_tr = T["vz"].to_numpy()
    ke_tr = T["ke_ion_ha"].to_numpy()

    # align track → energy sample times
    x = np.interp(tE, tT, s_tr)                           # path at each energy sample
    vz = np.interp(tE, tT, vz_tr)

    # early near-constant-velocity window (widen if sparse) — light-projectile rule
    used_vf = vfrac
    for vf in (vfrac, 0.70, 0.50):
        m = vz >= vf * v0
        m[:max(2, int(0.03 * len(x)))] = False
        used_vf = vf
        if m.sum() >= 20:
            break
    x0, xT = float(x[m].min()), float(x[m].max())

    prim = sp.free_fit(x, E, x0, xT)                      # PRIMARY: dE_total slope
    kin = sp.kinetic_channel(s_tr, ke_tr, x, x0, xT)      # sanity: −dKE_ion/dx
    nguard = (sp.conservation_guard(N) if N is not None
              else {"ok": None, "drained_frac": float("nan")})
    S, Sk = prim["S"], kin["S"]
    ratio = S / Sk if Sk else float("nan")
    dE_win = float(np.interp(xT, x, E) - np.interp(x0, x, E))
    dKE_win = float(np.interp(x0, s_tr, ke_tr) - np.interp(xT, s_tr, ke_tr))
    econs = dE_win / dKE_win if dKE_win else float("nan")

    steps = []

    def emit(name, fig, title, md):
        f = STEPDIR / name
        fig.tight_layout(); fig.savefig(f, dpi=160); plt.close(fig)
        steps.append((os.path.relpath(f, HERE), title, md))

    # ---- Step 1: the raw signals -------------------------------------------------
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.0))
    axs[0].plot(tE, e_abs, "C2-", lw=1.4)
    axs[0].set_xlabel("time (a.u.)"); axs[0].set_ylabel(r"$E_\mathrm{total}$ (Ha)")
    axs[0].set_title("electronic total energy (raw)")
    axs[1].plot(tT, z, "C0-", lw=1.4)
    axs[1].axhline(z0, ls="--", color="0.6", lw=0.8, label=f"launch z₀ = {z0:.0f}")
    axs[1].set_xlabel("time (a.u.)"); axs[1].set_ylabel("projectile z (Bohr)")
    axs[1].set_title("projectile position"); axs[1].legend(fontsize=8, frameon=False)
    axs[2].plot(tT, vz_tr, "C3-", lw=1.4, label=r"$v_z$")
    ax2b = axs[2].twinx()
    ax2b.plot(tT, ke_tr, "C1-", lw=1.2, label="KE_ion")
    axs[2].set_xlabel("time (a.u.)"); axs[2].set_ylabel(r"$v_z$ (a.u.)", color="C3")
    ax2b.set_ylabel("KE_ion (Ha)", color="C1")
    axs[2].set_title("projectile velocity + KE (it decelerates)")
    for a in axs:
        a.grid(alpha=.25)
    emit("step1_raw_signals.png", fig, "Step 1 — the raw measured signals",
         "**What we start from.** Two independent data streams for this run:\n\n"
         "- the **electronic** channel `observables.csv` → total electronic energy "
         "`E_total(t)` (left);\n"
         "- the **projectile** channel `electron_track.csv` → position `z(t)` (middle) "
         "and velocity `v_z(t)` + kinetic energy `KE_ion(t)` (right).\n\n"
         "The projectile is a classical mass-mₑ electron under **free Ehrenfest** "
         f"dynamics, launched on-axis at z₀ = {z0:.0f} Bohr with v₀ = {v0:.2f} a.u. "
         "Being light, it **decelerates strongly** — v_z and KE_ion collapse toward "
         "zero within the tube (the light-projectile rule): so S must be read from the "
         "*early* drag, not a steady-state plateau.")

    # ---- Step 2: N-conservation guard -------------------------------------------
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    if N is not None:
        tN = np.linspace(tE[0], tE[-1], len(N))
        ax.plot(tN, N, "C4-", lw=1.5)
        ax.axhline(N[0], ls="--", color="0.6", lw=0.8, label=f"N(0) = {N[0]:.2f}")
        drained = nguard["drained_frac"] * 100
        ax.set_title(f"Guard 1 — electron number N(t) "
                     f"(drained {drained:.2f}%, ok={nguard['ok']})")
    else:
        ax.text(0.5, 0.5, "electron_number.csv not found", ha="center")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$N_\mathrm{total}$ (electrons)")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=.25)
    emit("step2_N_guard.png", fig, "Step 2 — Guard: electron number N(t) ≈ const",
         "**Why.** The whole method assumes `E_total` changes only because the "
         "projectile deposits energy. If a CAP (or any sink) drained electrons, "
         "`E_total` would be dominated by that, not by the deposit, and **both** the "
         "energy and kinetic methods would be invalid (stopping-power-extraction §0). "
         "This is a **periodic tube with no absorber**, so N should be flat. We check "
         "`conservation_guard(N)`: drain ≤ 2% ⇒ pass.")

    # ---- Step 3: energy-conservation guard (both channels vs time) --------------
    dKE_loss_t = ke_tr[0] - ke_tr                          # KE lost by projectile
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.plot(tE, E * HA_EV, "C0-", lw=1.6, label=r"electronic gain $\Delta E_\mathrm{total}$")
    ax.plot(tT, dKE_loss_t * HA_EV, "C3--", lw=1.6,
            label=r"projectile loss $-\Delta KE_\mathrm{ion}$")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("energy (eV)")
    ax.set_title("Guard 2 — energy conservation: electronic gain ≈ projectile KE loss")
    ax.legend(fontsize=9, frameon=False); ax.grid(alpha=.25)
    emit("step3_energy_conservation.png", fig,
         "Step 3 — Guard: ΔE_total(t) ≈ −ΔKE_ion(t)",
         "**Why.** The two channels are independent (electronic `observables.csv` vs "
         "the classical track). If the energy the *gas* gains tracks the energy the "
         "*projectile* loses, the signal is trustworthy — this agreement **is** energy "
         "conservation. They need not be identical (some KE goes to the mean field), "
         "but they should rise together and stay within a few %.")

    # ---- Step 4: convert to path, deposit vs s ----------------------------------
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.plot(x, E * HA_EV, ".", color="0.5", ms=4)
    ax.set_xlabel(r"path  $s = |z - z_0|$  (Bohr)")
    ax.set_ylabel(r"$\Delta E_\mathrm{total}$ electronic deposit (eV)")
    ax.set_title("Step 4 — re-express the deposit against path travelled (all steps)")
    ax.grid(alpha=.25)
    emit("step4_deposit_vs_path.png", fig,
         "Step 4 — deposit vs path  ΔE_total(s)",
         "**The quantity S acts on.** Stopping power is energy per unit **path**, "
         "S = dE/dx (Correa 2018 Eq. 10). So we map each energy sample onto the "
         "distance the projectile has travelled, `s = |z − z₀|`, by interpolating the "
         "track onto the energy-sample times. The **slope** of this curve is S — but "
         "only over the right window (next step).")

    # ---- Step 5: window selection (light-projectile rule) -----------------------
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.plot(tE, vz, "C0-", lw=1.5)
    thr = used_vf * v0
    ax.axhline(thr, ls="--", color="0.5", lw=0.9,
               label=f"threshold {used_vf:.2f}·v₀ = {thr:.3f}")
    ax.axhspan(thr, vz.max() * 1.02, color="C0", alpha=0.10, label="early window (kept)")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$v_z$ (a.u.)")
    ax.set_title("Step 5 — select the EARLY near-constant-velocity window "
                 f"(v_z ≥ {used_vf:.2f}·v₀)")
    ax.legend(fontsize=9, frameon=False); ax.grid(alpha=.25)
    emit("step5_window_selection.png", fig,
         "Step 5 — the early v_z ≥ 0.85·v₀ window",
         "**The light-projectile rule.** A mass-mₑ electron loses most of its KE in a "
         "few Bohr, so its velocity sweeps from v₀ down to ~0 across the run. A slope "
         "fit over the *whole* run would average S over that entire velocity range — "
         "the wrong quantity. We want S **at v₀**, so we keep only the early segment "
         f"where v_z ≥ {vfrac:.2f}·v₀ (widened to 0.70/0.50·v₀ only if fewer than ~20 "
         f"points survive; here {used_vf:.2f}·v₀ was used, "
         f"{int(m.sum())} points). This defines the fit window "
         f"s ∈ [{x0:.2f}, {xT:.2f}] Bohr.")

    # ---- Step 6: the free-intercept fit (PRIMARY) -------------------------------
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.plot(x, E * HA_EV, ".", color="0.78", ms=3, label="all steps")
    ax.plot(x[m], E[m] * HA_EV, "C0.", ms=5, label="early window (fitted)")
    xs = np.linspace(x0, xT, 50)
    ax.plot(xs, (prim["S"] * xs + prim["E0"]) * HA_EV, "k-", lw=1.6,
            label=(f"free-intercept fit\nS = {S:.4f} Ha/Bohr = {S*HA_EV:.3f} eV/Bohr"
                   f"\nr² = {prim['r2']:.2f}"))
    ax.set_xlabel(r"path  $s = |z - z_0|$  (Bohr)")
    ax.set_ylabel(r"$\Delta E_\mathrm{total}$ electronic deposit (eV)")
    ax.set_title("Step 6 — PRIMARY method: S = slope of ΔE_total(s), free intercept",
                 weight="bold", fontsize=10)
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=.25)
    emit("step6_free_fit.png", fig,
         "Step 6 — fit the slope  →  the stopping power S",
         "**The defined stopping power.** A least-squares line "
         "`ΔE_total = S·s + E₀` over the early window, with the intercept **free** "
         "(never forced through the origin — the transient deposits a fixed E₀; forcing "
         "the origin biases S high, stopping-power-extraction §2). The **slope is S**; "
         "the error bar is the regression standard error. This slope — the electronic "
         "energy deposited per unit path — **is** the stopping power (Method A, "
         "continuous-traversal geometry).")

    # ---- Step 7: sanity channel (kinetic) ---------------------------------------
    dKE_loss = ke_tr[0] - np.interp(x, s_tr, ke_tr)
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.plot(x[m], E[m] * HA_EV, "C0.", ms=5,
            label=f"PRIMARY ΔE_total  → S = {S:.4f} Ha/Bohr")
    ax.plot(x[m], dKE_loss[m] * HA_EV, "C3.", ms=5,
            label=f"SANITY −ΔKE_ion  → S = {Sk:.4f} Ha/Bohr")
    ax.set_xlabel(r"path  $s$  (Bohr)"); ax.set_ylabel("deposited energy (eV)")
    ax.set_title(f"Step 7 — cross-check: kinetic channel vs primary "
                 f"(ratio {ratio:.2f})")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=.25)
    emit("step7_sanity_channel.png", fig,
         "Step 7 — the independent sanity channel  −dKE_ion/ds",
         "**Cross-check, never the headline.** The same slope is measured from the "
         "*projectile's* KE loss `−dKE_ion/ds` over the same window. This is a fully "
         "independent channel (the track, not `E_total`); their agreement is energy "
         "conservation. The energy-deposit slope from Step 6 is **THE** stopping power "
         "— the kinetic number is only a check (stopping-power-extraction §5). A "
         "divergence > 10% is flagged for the user, not silently averaged.")

    # ---- Step 8: the result card -------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    ax.plot(x[m], E[m] * HA_EV, "C0.", ms=5, label="early-window deposit")
    ax.plot(xs, (prim["S"] * xs + prim["E0"]) * HA_EV, "k-", lw=1.8)
    ax.set_xlabel(r"path  $s = |z - z_0|$  (Bohr)")
    ax.set_ylabel(r"$\Delta E_\mathrm{total}$ deposit (eV)")
    ax.set_title(f"Step 8 — result:  r_s = {rs}, v₀ = {v0:.2f} a.u.", fontsize=11,
                 weight="bold")
    flags = []
    if abs(ratio - 1) > 0.10:
        flags.append(f"channel divergence {abs(ratio-1)*100:.0f}%")
    if prim["r2"] < 0.80:
        flags.append(f"poor fit r²={prim['r2']:.2f}")
    if nguard.get("ok") is False:
        flags.append(f"N drained {nguard['drained_frac']*100:.1f}%")
    card = (f"S (headline)   = {S:.4f} Ha/Bohr\n"
            f"               = {S*HA_EV:.3f} eV/Bohr  (±{prim['se']*HA_EV:.3f})\n"
            f"r²             = {prim['r2']:.2f}\n"
            f"window s       = [{x0:.2f}, {xT:.2f}] Bohr\n"
            f"mean v in win  = {vz[m].mean():.3f} a.u.\n"
            f"KE sanity S    = {Sk*HA_EV:.3f} eV/Bohr (ratio {ratio:.2f})\n"
            f"energy cons.   = {econs:.2f}\n"
            f"N drained      = {nguard['drained_frac']*100:.2f}%\n"
            f"flags          = {', '.join(flags) if flags else 'none'}")
    ax.text(0.03, 0.97, card, transform=ax.transAxes, va="top", ha="left",
            family="monospace", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.6"))
    ax.grid(alpha=.25)
    emit("step8_result.png", fig,
         "Step 8 — the headline number, with its guards and flags",
         "**Reporting.** The headline S is the energy-deposit slope with its regression "
         "error; alongside it we report the window, the guards (N-conservation, energy "
         "conservation), the KE sanity ratio, and any flags. If a flag fires, both "
         "numbers are surfaced and the accept/reject verdict is the user's.")

    result = dict(S_ha=S, S_ev=S * HA_EV, se_ev=prim["se"] * HA_EV, r2=prim["r2"],
                  Sk_ev=Sk * HA_EV, ratio=ratio, econs=econs,
                  N_drained=nguard["drained_frac"], window=[x0, xT],
                  v_mean=float(vz[m].mean()), npts=int(m.sum()),
                  used_vf=used_vf, flags=flags)
    return steps, result


# ===================================================================== assemble
def build():
    FIGROOT.mkdir(parents=True, exist_ok=True)
    # generate figures
    raw_rs = {rs: energy_by_rs(rs, delta=False) for rs in RS_LIST}
    raw_v = {v: energy_by_velocity(v, delta=False) for v in VELS}
    del_rs = {rs: energy_by_rs(rs, delta=True) for rs in RS_LIST}
    del_v = {v: energy_by_velocity(v, delta=True) for v in VELS}
    steps, res = stopping_steps(rs=6, v0=0.45)

    cells = []
    md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))

    # 1 — title + question
    md("# Annular jellium tube — energy overlays & stopping-power walkthrough\n\n"
       "*How does the electronic total energy of the tube respond to the projectile, "
       "across the three wall densities r_s ∈ {6, 4, 2} and the three launch "
       "velocities v ∈ {0.15, 0.30, 0.45} a.u.? And — for one run (r_s = 6, "
       "v₀ = 0.45) — exactly how was the **stopping power** defined and measured, "
       "step by step?*\n\n"
       "This notebook has four parts:\n\n"
       "- **A** — total energy `E_total(t)`, one plot per r_s (traces = velocity);\n"
       "- **B** — total energy `E_total(t)`, one plot per velocity (traces = r_s);\n"
       "- **C** — the `ΔE = E(t) − E(0)` version of every plot in A and B;\n"
       "- **D** — the stopping-power measurement for r_s = 6, v₀ = 0.45, broken into "
       "its individual steps, each as its own plot.")

    # 2 — conventions
    md("## Conventions & symbols\n\n"
       "| symbol | meaning | units |\n|---|---|---|\n"
       "| `E_total(t)` | total electronic energy (from `observables.csv`) | Ha |\n"
       "| `ΔE_total` | deposit `E(t) − E(0)` | Ha (plotted in eV where noted) |\n"
       "| `s = \\|z − z₀\\|` | projectile path length | Bohr |\n"
       "| `v_z`, `KE_ion` | projectile axial velocity, kinetic energy | a.u., Ha |\n"
       "| `S` | stopping power `dE_total/ds` | Ha/Bohr (≡ 27.211 eV/Bohr) |\n"
       "| `r_s` | wall Wigner–Seitz radius (density index) | Bohr |\n\n"
       "**Baselines matter.** The three r_s runs contain different electron numbers "
       "(N = 24 / 48 / 136 for r_s = 6 / 4 / 2), so their **absolute** E(0) differ by "
       "~200 Ha (≈ +1.4, −1.5, −213 Ha). Raw `E_total` overlays across r_s are "
       "therefore not directly comparable — that is exactly why Part C re-references "
       "every curve to its own E(0). Unit: 1 Ha = 27.211 eV. Method source: "
       "**Correa 2018**, *Comput. Mater. Sci.* 150, 291 (Eq. 10, S = dE/dx).")

    # 3 — setup / 4 — sources (kept short; full geometry lives in the sweep notebooks)
    md("## Setup & source files\n\n"
       "Periodic annular jellium tube (axis ∥ z, wall between R_in = 5 and R_out = 13 "
       "Bohr); classical Gaussian-electron projectile (mass mₑ, free Ehrenfest, LDA, "
       "dt = 0.02) launched on-axis. Full geometry per r_s: N and L_z = "
       "24/48 · 48/28/10 Bohr for r_s = 6/4/2. See the per-r_s deep-look notebooks for "
       "the density evolution.\n\n"
       "| file | role |\n|---|---|\n"
       "| [`build_energy_and_stopping_notebook.py`](build_energy_and_stopping_notebook.py) | this builder |\n"
       "| [`per_run.py`](per_run.py) `stopping_analysis` | the extraction reproduced in Part D |\n"
       "| `stopping-power-extraction/stopping_power.py` | skill kernels (`free_fit`, `kinetic_channel`, `conservation_guard`) |\n"
       "| [`rs6_velocity_sweep.ipynb`](rs6_velocity_sweep.ipynb) · [`rs2_velocity_sweep.ipynb`](rs2_velocity_sweep.ipynb) | per-r_s density deep-looks |\n"
       "| [`annular_sv_report.ipynb`](annular_sv_report.ipynb) · [`annular_sv_index.ipynb`](annular_sv_index.ipynb) | run-SET report + guided index |\n\n"
       "Data per run: `observables.csv` (energy, current) and `electron_track.csv` "
       "(z, v_z, KE_ion); `electron_number.csv` for the N-conservation guard.")

    # PART A
    md("# Part A — total energy vs time, one plot per r_s\n\n"
       "Each plot fixes the wall density r_s and overlays the three projectile "
       "velocities. Raw `E_total(t)` in Ha (absolute).")
    for rs in RS_LIST:
        md(f"### r_s = {rs} ({RS_META[rs]['word']} wall)\n\n![E_total rs{rs}]({raw_rs[rs]})")

    # PART B
    md("# Part B — total energy vs time, one plot per velocity\n\n"
       "Each plot fixes the projectile velocity and overlays the three wall densities. "
       "Raw `E_total(t)` in Ha. **Note:** the three r_s traces sit at very different "
       "absolute levels (baselines ~200 Ha apart), so on a shared axis the deposit is "
       "invisible — the physically comparable view is the ΔE version in Part C.")
    for v in VELS:
        md(f"### v = {v:.2f} a.u.\n\n![E_total v{v:.2f}]({raw_v[v]})")

    # PART C
    md("# Part C — ΔE = E(t) − E(0) versions\n\n"
       "Subtracting each run's own baseline exposes the **deposited** energy on a "
       "common scale. First the per-r_s family (traces = velocity), then the "
       "per-velocity family (traces = r_s).")
    md("## C.1 — ΔE, one plot per r_s (traces = velocity)")
    for rs in RS_LIST:
        md(f"### r_s = {rs} ({RS_META[rs]['word']})\n\n![dE rs{rs}]({del_rs[rs]})")
    md("## C.2 — ΔE, one plot per velocity (traces = r_s)")
    for v in VELS:
        md(f"### v = {v:.2f} a.u.\n\n![dE v{v:.2f}]({del_v[v]})")

    # PART D
    md("# Part D — how the stopping power was defined  (r_s = 6, v₀ = 0.45)\n\n"
       "The stopping power is the electronic energy the projectile deposits per unit "
       "path, **S = dE_total/ds** (Correa 2018, Eq. 10). For this light, decelerating "
       "projectile it is measured by **Method A on the early drag window** — the exact "
       "procedure coded in `per_run.stopping_analysis`, which calls the "
       "stopping-power-extraction skill kernels. Below, each step of that procedure is "
       "rendered as its own plot, in order.")
    for rel, title, method in steps:
        md(f"## {title}\n\n{method}\n\n![{title}]({rel})")

    # Takeaway
    r = res
    md("# Takeaway\n\n"
       f"- **Deposit grows with both v and wall density.** The ΔE overlays (Part C) "
       "show a larger, faster energy deposit at higher velocity and at denser walls "
       "(smaller r_s) — consistent with stronger screening/drag in the dense tube.\n"
       "- **Raw vs ΔE.** Absolute `E_total` is dominated by the electron-count baseline "
       "(~200 Ha spread across r_s); only ΔE = E(t) − E(0) makes the deposit "
       "comparable — that is why both are shown.\n"
       f"- **Stopping power (r_s = 6, v₀ = 0.45):** the energy-deposit slope gives "
       f"**S = {r['S_ev']:.3f} eV/Bohr** (= {r['S_ha']:.4f} Ha/Bohr, ±{r['se_ev']:.3f}, "
       f"r² = {r['r2']:.2f}) over s ∈ [{r['window'][0]:.1f}, {r['window'][1]:.1f}] Bohr, "
       f"mean v = {r['v_mean']:.3f}. Independent KE sanity channel: "
       f"{r['Sk_ev']:.3f} eV/Bohr (ratio {r['ratio']:.2f}); N drained "
       f"{r['N_drained']*100:.2f}%. "
       f"Flags: {', '.join(r['flags']) if r['flags'] else 'none'}.\n"
       "- The **energy-deposit slope is the defined stopping power**; the kinetic "
       "channel is only the conservation cross-check (stopping-power-extraction §5).")

    nb = nbf.v4.new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})
    out = HERE / "energy_and_stopping.ipynb"
    nbf.write(nb, str(out))
    # verify every image reference resolves
    missing = []
    for c in cells:
        for tok in c.source.split("]("):
            if tok.startswith(("energy_stopping_figs", "./energy_stopping_figs")):
                rel = tok.split(")")[0]
                if not (HERE / rel).exists():
                    missing.append(rel)
    print(f"wrote {out}  ({len(cells)} cells)")
    print(f"image refs missing: {missing if missing else 'none — all resolve'}")
    print(f"stopping result: S = {res['S_ev']:.3f} eV/Bohr (r²={res['r2']:.2f}, "
          f"ratio {res['ratio']:.2f})")


if __name__ == "__main__":
    build()
