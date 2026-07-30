#!/usr/bin/env python3
"""Build the CAP energy-normalization investigation PHASE notebooks.

Companion to docs/plans/cap-energy-normalization-validation.md and
docs/notes/inq-energy-normalization-error.md. Consumes the completed vacuum WP
runs on GPU 1 (results/ under scripts/wp_traversal_energy) and the cross-run
aggregation investigation_summary.csv.

Produces, under this hypotheses folder:
  figures/*.png                      per-phase comparison figures (canonical theme)
  phase0_baseline_and_proof.ipynb    Phase 0 + 6 : the phenomenon + the /norm proof
  phase1_geometry_independence.ipynb Phase 1a    : eta-sweep, residual not reflection
  phase2_partial_absorption.ipynb    Phase 2     : E_ext/E0 vs norm identity line
  phase3_decisive_mask.ipynb         Phase 3     : mask ETRS vs CN, E_ext resolves it

Each phase notebook is a study (run-SET) notebook: context -> what was run ->
the phase figure -> the cross-run table -> the verdict against the plan falsifiers.
Study notebooks tabulate a run-SET; per-run deep dives are the run_report.ipynb.
"""
from __future__ import annotations
import sys, glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/local/data/public/skcb2/tddft")
sys.path.insert(0, str(ROOT / "inq-stack" / "python"))
from inqview.visualisation.style import apply_theme  # noqa: E402

apply_theme()
HA_EV = 27.211386
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)
RESULTS = ROOT / "ResearchProject/systems/vacuum/scripts/wp_traversal_energy/results"


def obs(tag: str, name: str) -> pd.DataFrame | None:
    p = RESULTS / tag / "raw" / "observables" / name
    return pd.read_csv(p, comment="#") if p.exists() else None


def norm_series(tag: str):
    """physical norm(t) from momentum-stats norm_check (write_every=1), = <psi|psi>(t)/<psi|psi>(0)."""
    m = obs(tag, "wp_momentum_stats.csv")
    if m is None:
        return None, None
    t = m["time_au"].to_numpy(float)
    nc = m["norm_check"].to_numpy(float)
    return t, nc / nc[0]


def ekin_series(tag: str):
    e = obs(tag, "energies.csv")
    if e is None:
        return None, None
    return e["time_au"].to_numpy(float), e["kinetic"].to_numpy(float)


SUMMARY = pd.read_csv(RESULTS / "investigation_summary.csv")


# ---------------------------------------------------------------- Phase 0 + 6
def fig_phase0():
    """The phenomenon (E_reported rises under CAP while E_ext decays) + the /norm proof."""
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for tag, lbl, c in [("nocap", "no CAP", "C0"), ("cap", "1-sided CAP", "C3")]:
        t, ek = ekin_series(tag)
        tn, nrm = norm_series(tag)
        if t is None:
            continue
        ax[0].plot(t, ek * HA_EV, c, label=f"E_reported ({lbl})")
        # E_ext = E_reported*norm on the (coarser) norm grid
        eki = np.interp(tn, t, ek)
        ax[0].plot(tn, eki * nrm * HA_EV, c, ls="--",
                   label=f"E_ext=E_rep·norm ({lbl})")
        ax[1].plot(tn, nrm, c, label=lbl)
    ax[0].set(xlabel="t (au)", ylabel="energy (eV)",
              title="Phase 0: reported vs extensive energy")
    ax[0].legend(fontsize=8)
    ax[1].set(xlabel="t (au)", ylabel="orbital norm  ⟨ψ|ψ⟩",
              title="Phase 0: norm decay under the CAP")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "phase0_phenomenon.png", dpi=140)
    plt.close(fig)

    # Phase 6 proof: energies.csv:kinetic == e_kin_ha (per-particle) != e_kin_ha*norm
    e = obs("cap", "energies.csv")
    m = obs("cap", "wp_momentum_stats.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    t = e["time_au"].to_numpy(float)
    k_en = e["kinetic"].to_numpy(float)
    ek_mom = np.interp(t, m["time_au"].to_numpy(float), m["e_kin_ha"].to_numpy(float))
    nrm = np.interp(t, *reversed(norm_series("cap")))  # placeholder; recompute below
    tn, ns = norm_series("cap")
    nrm = np.interp(t, tn, ns)
    ax.plot(t, k_en * HA_EV, "C0", lw=3, label="energies.csv : kinetic")
    ax.plot(t, ek_mom * HA_EV, "C1", ls="--", lw=1.6, label="wp_momentum_stats : e_kin_ha")
    ax.plot(t, ek_mom * nrm * HA_EV, "C3", ls=":", lw=2, label="e_kin_ha · norm  (extensive)")
    ax.set(xlabel="t (au)", ylabel="kinetic energy (eV)",
           title="Phase 6: INQ prints the /norm (per-particle) kinetic")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "phase6_norm_proof.png", dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------- Phase 1
def fig_phase1():
    df = SUMMARY[SUMMARY.group == "eta_sweep"].copy()
    df["eta"] = df["param"].str.extract(r"eta=(-?[0-9.]+)").astype(float)
    df = df.sort_values("eta")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(-df.eta, df.dE_rep_eV, "o-", color="C3")
    ax[0].set(xlabel="|η|  (CAP strength, Ha)", ylabel="ΔE_reported (eV)",
              title="Phase 1: reported rise grows with |η|…")
    ax[1].plot(df.frac_absorbed, df.dE_rep_eV, "o-", color="C0")
    ax[1].set(xlabel="fraction absorbed  (1 − norm_T)", ylabel="ΔE_reported (eV)",
              title="…because it tracks the absorbed fraction, not reflection")
    for a in ax:
        a.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(FIG / "phase1_eta_sweep.png", dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------- Phase 2
def fig_phase2():
    df = SUMMARY[SUMMARY.group.isin(["partial_abs", "eta_sweep"])].copy()
    fig, ax = plt.subplots(figsize=(6.0, 5.6))
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="identity  E_ext/E0 = norm")
    pa = df[df.group == "partial_abs"]
    es = df[df.group == "eta_sweep"]
    ax.scatter(pa.norm_T, pa.E_ext_frac, s=70, color="C0", zorder=3,
               label="Phase 2 partial-absorption")
    ax.scatter(es.norm_T, es.E_ext_frac, s=45, color="C3", marker="s", zorder=3,
               label="Phase 1 η-sweep")
    ax.set(xlabel="orbital norm  norm_T", ylabel="E_ext(T) / E0",
           title="Phase 2: extensive energy IS the absorbed norm")
    ax.legend(fontsize=8)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(FIG / "phase2_identity_line.png", dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------- Phase 3
def fig_phase3():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for tag, lbl, c in [("exp3a_mask_etrs", "mask + ETRS (norm-losing)", "C3"),
                        ("exp3b_mask_cn", "mask + CN (norm-preserving)", "C0")]:
        t, ek = ekin_series(tag)
        tn, nrm = norm_series(tag)
        if t is None:
            continue
        ax[0].plot(t, ek * HA_EV, c, label=lbl)
        eki = np.interp(tn, t, ek)
        ax[1].plot(tn, eki * nrm * HA_EV, c, label=f"E_ext  {lbl}")
    ax[0].set(xlabel="t (au)", ylabel="E_reported (eV)",
              title="Phase 3: both rise ~+19 eV in the reported energy")
    ax[0].legend(fontsize=8)
    ax[1].set(xlabel="t (au)", ylabel="E_ext = E_reported·norm (eV)",
              title="…but E_ext SEPARATES them: ETRS→0, CN pumps +4.8%")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "phase3_decisive.png", dpi=140)
    plt.close(fig)


def main():
    fig_phase0()
    fig_phase1()
    fig_phase2()
    fig_phase3()
    print("[phasefig] wrote", *(p.name for p in sorted(FIG.glob("*.png"))))


if __name__ == "__main__":
    main()
