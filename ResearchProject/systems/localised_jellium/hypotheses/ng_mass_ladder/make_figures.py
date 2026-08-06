#!/usr/bin/env python3
"""Nazarov-Gross validation figures.

Plan: docs/plans/nazarov-gross-slab-mass-ladder.md

EVERY FIGURE HERE ANSWERS A SPECIFIC PART OF THE CLAIM. The campaign is not a
survey; it is a test of two statements, so the figure set is organised around
them rather than around what happens to be plottable:

  ng_01_S_vs_mass          THE CLAIM. S at fixed charge and velocity, versus M,
                           with the classical M->inf run as the anchor.
  ng_02_S_vs_width         THE MECHANISM. S versus the MEASURED mid-transit
                           width, mass ladder and sigma sweep overlaid. If mass
                           acts only through width they lie on one curve.
  ng_03_width_vs_time      mass -> width, the first link of the chain, measured.
  ng_04_deposit_vs_path    the raw observable S is fitted from, fit window shaded.
  ng_05_interactions       the P/S/B pairwise decomposition (E_PP is the quantum
                           residual with no classical counterpart).
  ng_06_kinetic_channels   where the projectile's energy went: bath vs its own
                           spreading. T2 - T1 = var(p)/2M has no classical twin.
  ng_07_ledger_closure     the correctness panel: energy conservation and the
                           interactions-vs-INQ closure gates.

House standard: canonical theme, no on-canvas titles, PNG only, one panel per
file, annotated with its own headline number (.claude/skills/report-figures).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
REPO = HERE.parents[4]
sys.path.insert(0, str(REPO / "inq-stack/python"))

import ng_analysis as NG                  # noqa: E402

try:
    from inqview.visualisation import style
    style.apply_theme()
except Exception:                          # noqa: BLE001
    pass                                   # figures still render with defaults

CL_COLOR, WP_COLOR = "#1f4e9c", "#c0392b"   # classical = blue, wavepacket = red
SWEEP_COLOR = "#d68910"


def _fig(w=5.0, h=3.8):
    f, a = plt.subplots(figsize=(w, h))
    return f, a


def _save(fig, out: Path, name: str) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{name}.png"
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p.name}")
    return p


# --------------------------------------------------------------- 01 THE CLAIM
def fig_S_vs_mass(tbl, out: Path):
    q = tbl[(tbl.half == "wp") & np.isfinite(tbl.S_ev_per_bohr)].sort_values("mass")
    cl = tbl[(tbl.half == "classical") & np.isfinite(tbl.S_ev_per_bohr)]
    if q.empty:
        return None
    fig, ax = _fig()
    ax.errorbar(q.mass, q.S_ev_per_bohr, yerr=q.S_stderr, marker="o", ms=6,
                lw=1.6, capsize=3, color=WP_COLOR, label="quantum wavepacket", zorder=3)
    for _, row in cl.iterrows():
        inf = row.mass > 1e5
        ax.axhline(row.S_ev_per_bohr, ls="--" if inf else ":", lw=1.4, color=CL_COLOR,
                   label=("classical, $M\\to\\infty$" if inf else "classical, $M=1$"), zorder=2)
    ax.set_xscale("log")
    ax.set_xlabel(r"projectile mass $M$  ($m_e$)")
    ax.set_ylabel(r"$S = \mathrm{d}E_{\mathrm{bath}}/\mathrm{d}s$  (eV/Bohr)")
    v = NG.does_S_depend_on_mass(tbl)
    if "n_sigma" in v:
        ax.annotate(f"spread {v['spread']:.4f} eV/Bohr = {v['n_sigma']:.1f}$\\sigma$\n{v['verdict']}",
                    xy=(0.03, 0.97), xycoords="axes fraction", va="top", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.9))
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    return _save(fig, out, "ng_01_S_vs_mass")


# ----------------------------------------------------------- 02 THE MECHANISM
def fig_S_vs_width(tbl, out: Path):
    d = tbl[(tbl.half == "wp") & np.isfinite(tbl.S_ev_per_bohr)
            & np.isfinite(tbl.sigma_iso_mid)]
    if len(d) < 2:
        return None
    is_ladder = np.isclose(d.sigma_WP_nominal, NG.SIGMA_WP_DEFAULT)
    fig, ax = _fig()
    ax.errorbar(d[is_ladder].sigma_iso_mid, d[is_ladder].S_ev_per_bohr,
                yerr=d[is_ladder].S_stderr, ls="none", marker="o", ms=7, capsize=3,
                color=WP_COLOR, label=r"mass ladder (fixed $\sigma_0$)")
    ax.errorbar(d[~is_ladder].sigma_iso_mid, d[~is_ladder].S_ev_per_bohr,
                yerr=d[~is_ladder].S_stderr, ls="none", marker="s", ms=7, capsize=3,
                mfc="none", color=SWEEP_COLOR, label=r"$\sigma$ sweep (fixed $M=1$)")
    c = NG.does_S_collapse_on_width(tbl)
    if "power_law_exponent" in c and len(d) >= 3:
        xs = np.linspace(d.sigma_iso_mid.min(), d.sigma_iso_mid.max(), 50)
        lnA = np.log(np.maximum(d.S_ev_per_bohr, 1e-12)).mean() - c["power_law_exponent"] * np.log(d.sigma_iso_mid).mean()
        ax.plot(xs, np.exp(lnA) * xs ** c["power_law_exponent"], lw=1.2, color="0.4", ls="-",
                label=fr"$S\propto\sigma^{{{c['power_law_exponent']:.2f}}}$, $r^2={c['r2']:.3f}$")
        ax.annotate(c["verdict"], xy=(0.03, 0.05), xycoords="axes fraction", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.9))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"measured mid-transit width  $\sigma_{\mathrm{iso}}$  (Bohr)")
    ax.set_ylabel(r"$S$  (eV/Bohr)")
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, out, "ng_02_S_vs_width")


# -------------------------------------------------------- 03 mass -> width
def fig_width_vs_time(runs, out: Path):
    fig, ax = _fig()
    any_curve = False
    for r in runs:
        w = NG.wp_width(r)
        if w.empty:
            continue
        any_curve = True
        ax.plot(w.time_au, w.sigma_iso, lw=1.5, label=f"$M={r.mass:g}$")
    if not any_curve:
        plt.close(fig)
        return None
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\sigma_{\mathrm{iso}}(t)$  (Bohr)")
    ax.annotate("lighter packets disperse faster — this is the\nfirst link of the NG chain, measured directly",
                xy=(0.03, 0.97), xycoords="axes fraction", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.9))
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, out, "ng_03_width_vs_time")


# ------------------------------------------------- 04 the raw fitted quantity
def fig_deposit_vs_path(runs, out: Path):
    fig, ax = _fig()
    drew = False
    for r in runs:
        dep = NG.bath_energy_ev(r)
        trk = NG.projectile_track(r)
        if dep is None or trk.empty:
            continue
        z = np.interp(r.obs["time_au"], trk["time_au"], trk["z"])
        col = CL_COLOR if r.half == "classical" else WP_COLOR
        ls = "--" if r.half == "classical" else "-"
        lbl = "classical" if r.half == "classical" else f"WP $M={r.mass:g}$"
        ax.plot(z, dep, lw=1.4, color=col, ls=ls, label=lbl, alpha=0.9)
        drew = True
    if not drew:
        plt.close(fig)
        return None
    ax.axvspan(-NG.SLAB_HALF, NG.SLAB_HALF, color="0.85", alpha=0.5, zorder=0)
    ax.annotate("slab", xy=(0, ax.get_ylim()[1]), ha="center", va="top", fontsize=8, color="0.35")
    ax.set_xlabel(r"projectile position $z$  (Bohr)")
    ax.set_ylabel(r"$\Delta E_{\mathrm{bath}}$  (eV)")
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, out, "ng_04_deposit_vs_path")


# ------------------------------------------- 05 pairwise energy decomposition
def fig_interactions(r, out: Path, suffix=""):
    if r.inter.empty:
        return None
    fig, ax = _fig(5.6, 3.8)
    t = r.inter["time_au"]
    for col, lab in [("e_ss", r"$E_{SS}$ bath-bath"), ("e_pp", r"$E_{PP}$ projectile self"),
                     ("e_ps", r"$E_{PS}$ projectile-bath"), ("e_sb", r"$E_{SB}$ bath-background"),
                     ("e_pb", r"$E_{PB}$ projectile-background")]:
        if col in r.inter:
            ax.plot(t, (r.inter[col] - r.inter[col].iloc[0]) * NG.HA_EV, lw=1.3, label=lab)
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\Delta E$ from $t=0$  (eV)")
    ax.annotate(r"$E_{PP}$ is the projectile self-Hartree — it has"
                "\nNO classical counterpart (a rigid cloud's is constant)",
                xy=(0.03, 0.03), xycoords="axes fraction", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))
    ax.legend(frameon=False, fontsize=7.5, ncol=2)
    return _save(fig, out, f"ng_05_interactions{suffix}")


# ------------------------------------------------- 06 where the energy went
def fig_kinetic_channels(r, out: Path, suffix=""):
    k = NG.kinetic_channels(r)
    if k.empty:
        return None
    fig, ax = _fig()
    ax.plot(k.time_au, k.T1_drift_ev - k.T1_drift_ev.iloc[0], lw=1.5, color=WP_COLOR,
            label=r"$\Delta T_1=\Delta\langle p\rangle^2/2M$ (drift)")
    ax.plot(k.time_au, k.T2_total_ev - k.T2_total_ev.iloc[0], lw=1.5, color="0.35",
            label=r"$\Delta T_2=\Delta\langle p^2\rangle/2M$ (total)")
    ax.plot(k.time_au, k.var_p_over_2m_ev - k.var_p_over_2m_ev.iloc[0], lw=1.5,
            color=SWEEP_COLOR, ls="--", label=r"$\Delta\,\mathrm{var}(p)/2M$ (own spreading)")
    dep = NG.bath_energy_ev(r)
    if dep is not None:
        ax.plot(r.obs["time_au"], dep, lw=1.5, color=CL_COLOR, ls=":",
                label=r"$\Delta E_{\mathrm{bath}}$ (what actually stopped it)")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\Delta E$  (eV)")
    ax.annotate("the gap between $-\\Delta T_1$ and $\\Delta E_{bath}$ is energy the\n"
                "packet absorbed into ITSELF and never delivered",
                xy=(0.03, 0.03), xycoords="axes fraction", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))
    ax.legend(frameon=False, fontsize=7.5)
    return _save(fig, out, f"ng_06_kinetic_channels{suffix}")


# --------------------------------------------------------- 07 correctness
def fig_ledger_closure(runs, out: Path):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    drew = False
    for r in runs:
        if r.obs.empty or "energy_total" not in r.obs:
            continue
        e = r.obs["energy_total"].to_numpy()
        a1.plot(r.obs["time_au"], (e - e[0]) * NG.HA_EV, lw=1.2,
                label=f"{r.tag}", alpha=0.9)
        drew = True
        if r.half == "wp" and not r.inter.empty and "e_hartree_check" in r.inter:
            t = r.inter["time_au"]
            hart = np.interp(t, r.obs["time_au"], r.obs.get("energy_hartree", r.obs["energy_total"]))
            a2.plot(t, (r.inter["e_hartree_check"] - hart) * NG.HA_EV, lw=1.2, label=r.tag)
    if not drew:
        plt.close(fig)
        return None
    a1.set_xlabel("time (a.u.)"); a1.set_ylabel(r"$\Delta E_{\mathrm{total}}$ (eV)")
    a1.legend(frameon=False, fontsize=7)
    a2.set_xlabel("time (a.u.)")
    a2.set_ylabel(r"$E_{SS}\!+\!E_{PS}\!+\!E_{PP}-E_{\mathrm{Hartree}}^{\mathrm{INQ}}$ (eV)")
    a2.annotate("closure gate: must be ~0", xy=(0.04, 0.92), xycoords="axes fraction", fontsize=7.5)
    a2.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    return _save(fig, out, "ng_07_ledger_closure")


# --------------------------------------------------------------------- driver
def build_all(scripts_dir: Path, run_specs: list[tuple[str, str]], out: Path) -> list[Path]:
    """run_specs: [(half, tag), ...]. Missing/incomplete runs are skipped."""
    runs = [NG.load_run(scripts_dir, h, t) for h, t in run_specs]
    runs = [r for r in runs if r.complete]
    if not runs:
        print("  no completed runs — nothing to plot")
        return []
    tbl = NG.ladder_table(scripts_dir, run_specs)
    tbl.to_csv(out.parent / "ng_ladder_table.csv", index=False)

    made = []
    for fn in (lambda: fig_S_vs_mass(tbl, out),
               lambda: fig_S_vs_width(tbl, out),
               lambda: fig_width_vs_time([r for r in runs if r.half == "wp"], out),
               lambda: fig_deposit_vs_path(runs, out),
               lambda: fig_ledger_closure(runs, out)):
        try:
            p = fn()
            if p:
                made.append(p)
        except Exception as exc:                                # noqa: BLE001
            print(f"  [warn] a figure failed: {exc!r}")
    for r in runs:
        for fn in (lambda r=r: fig_interactions(r, out, f"_{r.tag}"),
                   lambda r=r: fig_kinetic_channels(r, out, f"_{r.tag}")):
            try:
                p = fn()
                if p:
                    made.append(p)
            except Exception as exc:                            # noqa: BLE001
                print(f"  [warn] per-run figure failed for {r.tag}: {exc!r}")
    return made


if __name__ == "__main__":
    scripts = REPO / "ResearchProject/systems/localised_jellium/scripts/ng_mass_ladder"
    specs = [("classical", "cl_inf"), ("classical", "cl_m1"),
             ("wp", "wp_m3"), ("wp", "wp_m1p2"), ("wp", "wp_m1"), ("wp", "wp_m0p5"),
             ("wp", "wp_m1_s2p0"), ("wp", "wp_m1_s3p0"), ("wp", "wp_m1_s6p0")]
    made = build_all(scripts, specs, HERE / "figures")
    print(f"{len(made)} figures written to {HERE/'figures'}")
