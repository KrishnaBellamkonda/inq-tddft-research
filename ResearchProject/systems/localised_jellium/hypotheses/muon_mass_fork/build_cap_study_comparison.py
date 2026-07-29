#!/usr/bin/env python3
"""Build the CAP-parameter comparison STUDY notebook for the σ=1 WP run.

Overlays the four CAP settings (baseline η=−1.0, R1 wider-gap, R2 weak, R3 strong)
on the decision metrics — total N(t)=∫n dV (the primary metric, per the user:
absorption is judged on the whole-cell electron count, not the WP orbital norm),
reflection, WP absorption, wake-under-CAP, and E_total(t). Extracts the time series
once (VTI + CSV), writes a combined CSV, renders LINEAR|LOG comparison figures
(enforced house rule), and assembles+executes the notebook.

Run: PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_cap_study_comparison.py
"""
from __future__ import annotations
import glob, re, os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
try:
    from inqview.visualisation.style import apply_theme; apply_theme()
except Exception:
    pass
from inqview.visualisation.field_io import load_vti

ROOT = Path("/local/data/public/skcb2/tddft")
BASE = ROOT / "ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/effmass_sigma1/wp/results"
OUT  = ROOT / "ResearchProject/systems/localised_jellium/hypotheses/muon_mass_fork"
FIG  = OUT / "cap_study_comparison_figs"; FIG.mkdir(exist_ok=True)
DT   = 0.04
LAUNCH_Z = -16.5

# key, label, results-dir, CAP inner |z|, η, region, gap
RUNS = [
    ("baseline", "baseline η=−1.0 [25,40]",  "sigma1",      25.0, -1.0, "[25,40]", 12.5),
    ("gap19p5",  "R1 gap19.5 η=−1.0 [32,40]","cap_gap19p5", 32.0, -1.0, "[32,40]", 19.5),
    ("eta0p4",   "R2 weak η=−0.4 [25,40]",   "cap_eta0p4",  25.0, -0.4, "[25,40]", 12.5),
    ("eta2p0",   "R3 strong η=−2.0 [25,40]", "cap_eta2p0",  25.0, -2.0, "[25,40]", 12.5),
]
COL = {"baseline": "C0", "gap19p5": "C1", "eta0p4": "C3", "eta2p0": "C2"}
# per-run density GIFs already rendered by the per-run notebooks (linear|log) — reuse them.
GIFSRC = {
    "baseline": ("effmass_sigma1_wp_run_figs",          "sigma1"),
    "gap19p5":  ("effmass_sigma1_cap_gap19p5_wp_run_figs", "cap_gap19p5"),
    "eta0p4":   ("effmass_sigma1_cap_eta0p4_wp_run_figs",  "cap_eta0p4"),
    "eta2p0":   ("effmass_sigma1_cap_eta2p0_wp_run_figs",  "cap_eta2p0"),
}
_zt = lambda f: int(re.search(r"_t(\d+)", f).group(1))


def df_to_md(df):
    """Minimal DataFrame → GitHub markdown table (no 'tabulate' dependency)."""
    cols = list(df.columns)
    head = "| " + " | ".join(str(c) for c in cols) + " |"
    sep  = "| " + " | ".join("---" for _ in cols) + " |"
    body = "\n".join("| " + " | ".join(str(v) for v in row) + " |"
                     for row in df.itertuples(index=False, name=None))
    return "\n".join([head, sep, body])


def extract():
    """Return {key: dict of time series} and the combined tidy DataFrame."""
    series, rows = {}, []
    for key, label, d, inner, eta, region, gap in RUNS:
        rd = BASE / d
        # ---- VTI series: N(t)=∫n dV, wake under CAP (|z|≥inner), reflection ----
        tot = sorted(glob.glob(str(rd / "raw/vti/density_total/*.vti")), key=_zt)
        wp  = sorted(glob.glob(str(rd / "raw/vti/density_wp/*.vti")),    key=_zt)
        tV, N, wake = [], [], []
        for f in tot:
            v = load_vti(f); s = v.spacing; z = v.z; dV = s[0]*s[1]*s[2]
            tV.append(_zt(f)*DT); N.append(v.data.sum()*dV)
            wake.append(v.data[:, :, np.abs(z) >= inner].sum()*dV)
        tW, refl = [], []
        n0 = None
        for i, f in enumerate(wp):
            v = load_vti(f); s = v.spacing; z = v.z; dV = s[0]*s[1]*s[2]
            if i == 0:
                n0 = v.data.sum()*dV
            tW.append(_zt(f)*DT)
            refl.append(v.data[:, :, z < LAUNCH_Z].sum()*dV / max(n0, 1e-30))
        # ---- CSV series: WP orbital norm, E_total ----
        try:
            rs = pd.read_csv(rd / "raw/observables/wp_real_space_stats.csv", comment="#")
            t_norm, wpnorm = rs.time_au.to_numpy(), rs.norm_check.to_numpy()
        except Exception:
            t_norm, wpnorm = np.array([]), np.array([])
        try:
            o = pd.read_csv(rd / "raw/observables/observables.csv")
            t_E, E = o.time_au.to_numpy(), o.energy_total.to_numpy()
        except Exception:
            t_E, E = np.array([]), np.array([])
        series[key] = dict(label=label, inner=inner, eta=eta, region=region, gap=gap,
                           tV=np.array(tV), N=np.array(N), wake=np.array(wake),
                           tW=np.array(tW), refl=np.array(refl),
                           t_norm=t_norm, wpnorm=wpnorm, t_E=t_E, E=E)
        # tidy rows (VTI cadence) for the combined CSV
        for k in range(len(tV)):
            rows.append(dict(run=key, label=label, eta=eta, region=region, gap=gap,
                             t_au=tV[k], N_total=N[k], wake_under_cap=wake[k],
                             refl_frac=refl[k] if k < len(refl) else np.nan))
    return series, pd.DataFrame(rows)


def _lin_log(ax_pair, plot_fn, title_l, title_r, ylabel, log_ylabel=None):
    """Apply plot_fn(ax, logscale) to a (linear, log) axis pair — enforced house rule."""
    axL, axR = ax_pair
    plot_fn(axL, False); axL.set_title(title_l); axL.set_ylabel(ylabel); axL.grid(alpha=.25)
    plot_fn(axR, True);  axR.set_title(title_r); axR.set_yscale("log")
    axR.set_ylabel(log_ylabel or ylabel); axR.grid(alpha=.25, which="both")
    for a in (axL, axR):
        a.set_xlabel("time (a.u.)"); a.legend(fontsize=7, frameon=False)


def figures(S):
    figs = {}
    # --- 1. total N(t): the DECISION metric (linear full | linear zoom to bath) ---
    f, (aL, aR) = plt.subplots(1, 2, figsize=(12.4, 4.2))
    for k, s in S.items():
        aL.plot(s["tV"], s["N"], color=COL[k], lw=1.6, marker="o", ms=2.5, label=s["label"])
        aR.plot(s["tV"], s["N"], color=COL[k], lw=1.6, marker="o", ms=2.5, label=s["label"])
    aL.axhline(53, ls=":", c="0.6"); aL.axhline(52, ls=":", c="0.6")
    aL.set_title("Total electron number N(t)=∫n dV"); aL.set_ylabel("N (electrons)")
    aR.set_ylim(51.9, 53.05); aR.axhline(52, ls="--", c="0.5", lw=.8)
    aR.set_title("Zoom on the bath floor (N=52)  →  no CAP eats slab charge")
    for a in (aL, aR):
        a.set_xlabel("time (a.u.)"); a.grid(alpha=.25); a.legend(fontsize=7, frameon=False)
    f.tight_layout(); p = FIG / "cmp_N_of_t.png"; f.savefig(p, dpi=150); plt.close(f); figs["N"] = p

    # --- 2. WP orbital norm(t): absorption completeness (linear | log) ---
    f, axp = plt.subplots(1, 2, figsize=(12.4, 4.2))
    def _norm(ax, log):
        for k, s in S.items():
            if len(s["t_norm"]):
                y = np.clip(s["wpnorm"], 1e-6, None) if log else s["wpnorm"]
                ax.plot(s["t_norm"], y, color=COL[k], lw=1.5, label=s["label"])
    _lin_log(axp, _norm, "WP orbital norm (linear)", "WP orbital norm (log)",
             "‖ψ_WP‖²", "‖ψ_WP‖² (log)")
    f.tight_layout(); p = FIG / "cmp_wp_norm.png"; f.savefig(p, dpi=150); plt.close(f); figs["norm"] = p

    # --- 3. reflection: backward WP density past launch (linear | log) ---
    f, axp = plt.subplots(1, 2, figsize=(12.4, 4.2))
    def _refl(ax, log):
        for k, s in S.items():
            if len(s["tW"]):
                y = 100*np.clip(s["refl"], 1e-6, None) if log else 100*s["refl"]
                ax.plot(s["tW"], y, color=COL[k], lw=1.5, label=s["label"])
    _lin_log(axp, _refl, "Reflected fraction (linear)", "Reflected fraction (log)",
             "back of launch  z<−16.5  (% of N_wp0)", "reflected % (log)")
    f.tight_layout(); p = FIG / "cmp_reflection.png"; f.savefig(p, dpi=150); plt.close(f); figs["refl"] = p

    # --- 4. wake charge under the CAP region (linear) + E_total(t) ---
    f, (aL, aR) = plt.subplots(1, 2, figsize=(12.4, 4.2))
    for k, s in S.items():
        aL.plot(s["tV"], s["wake"], color=COL[k], lw=1.5, label=f"{s['label']} (|z|≥{s['inner']:.0f})")
        if len(s["t_E"]):
            aR.plot(s["t_E"], s["E"], color=COL[k], lw=1.5, label=s["label"])
    aL.set_title("Charge under the CAP region (WP+wake)"); aL.set_ylabel("∫_{|z|≥inner} n dV (e⁻)")
    aR.set_title("E_total(t) — CAP energy accounting"); aR.set_ylabel("E_total (Ha)")
    for a in (aL, aR):
        a.set_xlabel("time (a.u.)"); a.grid(alpha=.25); a.legend(fontsize=7, frameon=False)
    f.tight_layout(); p = FIG / "cmp_wake_energy.png"; f.savefig(p, dpi=150); plt.close(f); figs["wake"] = p
    return figs


def per_run_en_fig(key, s):
    """Per-run 2-panel [ E_total(t) | N(t)=∫n dV ] figure."""
    f, (aE, aN) = plt.subplots(1, 2, figsize=(11.6, 4.0))
    if len(s["t_E"]):
        aE.plot(s["t_E"], s["E"], color=COL[key], lw=1.5)
    aE.set_xlabel("time (a.u.)"); aE.set_ylabel(r"$E_\mathrm{total}$ (Ha)")
    aE.set_title("Total energy vs time"); aE.grid(alpha=.25)
    aN.plot(s["tV"], s["N"], color=COL[key], lw=1.6, marker="o", ms=2.5)
    aN.axhline(52, ls="--", c="0.5", lw=.8)
    aN.set_xlabel("time (a.u.)"); aN.set_ylabel(r"$N(t)=\int n\,dV$  (electrons)")
    aN.set_title("Total electron number vs time"); aN.set_ylim(51.9, 53.05); aN.grid(alpha=.25)
    f.suptitle(s["label"], fontsize=10)
    f.tight_layout(rect=(0, 0, 1, 0.96))
    p = FIG / f"run_{key}_EN.png"; f.savefig(p, dpi=150); plt.close(f)
    return p


def summary_table(S):
    def late_refl(s):
        if not len(s["tW"]):
            return np.nan
        m = s["tW"] > 0.5*s["tW"].max()
        return 100*np.nanmax(s["refl"][m]) if m.any() else np.nan
    rows = []
    for k, s in S.items():
        rows.append(dict(Run=s["label"], eta=s["eta"], region=s["region"], gap=s["gap"],
                         N_min=round(float(np.min(s["N"])), 3),
                         N_end=round(float(s["N"][-1]), 2),
                         WP_absorbed_pct=round(100*(1-float(s["wpnorm"][-1])) if len(s["wpnorm"]) else np.nan, 2),
                         refl_late_pct=round(late_refl(s), 3),
                         E_end_Ha=round(float(s["E"][-1]) if len(s["E"]) else np.nan, 2)))
    return pd.DataFrame(rows)


def build():
    S, tidy = extract()
    csv = OUT / "cap_study_timeseries.csv"; tidy.to_csv(csv, index=False)
    figs = figures(S)
    tab = summary_table(S); tab.to_csv(OUT / "cap_study_summary.csv", index=False)
    rel = lambda p: os.path.relpath(p, OUT)

    nb = new_notebook()
    C = nb.cells
    C.append(new_markdown_cell(
        "# CAP-parameter comparison — σ=1 wavepacket through the r_s≈5.68 slab\n\n"
        "**Question.** Is the complex absorbing potential (CAP) of the σ=1 WP run contaminating the "
        "physics — reflecting the packet, or *eating slab/wake charge*? Four CAP settings of the **same** "
        "run are overlaid. **Decision metric = total electron number N(t)=∫n dV** (absorption judged on the "
        "whole-cell count, NOT the WP orbital norm alone): if the CAP only removes the 1-e⁻ projectile, "
        "N(t) settles at 52 (the slab count); if it eats the bath/wake, N dips below 52."))
    C.append(new_markdown_cell(
        "## Conventions & the four runs\n"
        "Box 40×40×80 Bohr, dx=0.333; slab faces at |z|=12.5; WP σ=1, k₀=5.693, m=2.10, launched z=−16.5. "
        "CAP is a two-sided sin² absorber `perturbations::absorbing(η·Ha, ±center, width)`. "
        "**σ = σ_WP** throughout. N(t) from ∫`density_total` VTIs; reflection = backward WP density past the "
        "launch point (z<−16.5) in the late half of the run; WP absorption from the orbital norm.\n\n"
        + df_to_md(tab)))
    C.append(new_markdown_cell(
        "## Setup / provenance (reconstructable)\n"
        "- Runs: `scripts/muon_mass_fork/effmass_sigma1/wp/results/{sigma1,cap_gap19p5,cap_eta0p4,cap_eta2p0}`\n"
        "- Parametrised binary `effmass_sigma1/wp/run.cpp` (env `EM_CAP_ETA/CENTER_BOHR/WIDTH_BOHR`); shared GS "
        "`shared_gs/slab_n52_L40x40x80_dx0p333` (E=−36.9405 Ha).\n"
        "- Orchestrator `effmass_sigma1/cap_study/orchestrate.sh`; this builder "
        "`hypotheses/muon_mass_fork/build_cap_study_comparison.py`; data `cap_study_timeseries.csv`.\n"
        "- Per-run deep notebooks: `effmass_sigma1_cap_{gap19p5,eta0p4,eta2p0}_wp_run.ipynb` + baseline "
        "`effmass_sigma1_wp_run.ipynb`."))
    C.append(new_markdown_cell(
        "## 1. Total electron number N(t) — the decision metric\n"
        "Left: full range (53 at launch → 52 after the 1-e⁻ WP is absorbed). Right: zoom on the bath floor. "
        "**All four settings hold N=52 (N_min ∈ [51.999, 52.012]) — no CAP eats slab/wake charge.**"))
    C.append(new_markdown_cell(f"![N]({rel(figs['N'])})"))
    C.append(new_markdown_cell(
        "## 2. WP absorption completeness (linear | log)\n"
        "The orbital norm →0 as the packet drains into the CAP. Baseline and strong η reach 0 (complete); "
        "**weak η=−0.4 plateaus at ~1.25% un-absorbed** (optical depth too low)."))
    C.append(new_markdown_cell(f"![norm]({rel(figs['norm'])})"))
    C.append(new_markdown_cell(
        "## 3. Reflection — backward WP density past the launch point (linear | log)\n"
        "Genuine reflection (z<−16.5, behind launch). **Weak η=−0.4 reflects 4.1%** (leaky branch: the packet "
        "reaches the box edge and bounces); baseline 0.13%, strong η=−2.0 ≈0.001%. Reflection *falls* with |η| "
        "here — we are on the under-damped side; η≥1 is clean."))
    C.append(new_markdown_cell(f"![refl]({rel(figs['refl'])})"))
    C.append(new_markdown_cell(
        "## 4. Charge under the CAP region, and E_total(t)\n"
        "Left: transient charge (WP + induced wake) under |z|≥inner. Wider gap (R1, inner=32) sees least. "
        "Right: **E_total is CAP-strength-dependent** — strong η=−2.0 ends at −34.9 Ha vs −27.7 elsewhere, "
        "because the CAP removes strength-dependent energy. ⇒ the WP **energy-method stopping is a CAP-dependent "
        "bound, not a clean observable**; take S from the classical twin / near-field WP momentum."))
    C.append(new_markdown_cell(f"![wake]({rel(figs['wake'])})"))

    # ---- per-run detail: density GIFs + E(t) + N(t), one section each ----
    C.append(new_markdown_cell(
        "# Per-run detail\n"
        "Each CAP setting below: the **total-density xz evolution** (linear | log) and the "
        "**induced Δn = n(t)−n(0)** (linear | symlog, the wake), then **E_total(t)** and "
        "**N(t)=∫n dV** side by side. Dashed lines = slab faces (±12.5) and CAP inner faces."))
    trow = {r.Run: r for r in tab.itertuples(index=False)}
    for key, s in S.items():
        en = per_run_en_fig(key, s)
        gdir, gpref = GIFSRC[key]
        dens = f"{gdir}/{gpref}_total_density.gif"
        delt = f"{gdir}/{gpref}_total_delta0.gif"
        r = trow.get(s["label"])
        stat = (f" — N_min **{r.N_min}**, WP absorbed **{r.WP_absorbed_pct}%**, "
                f"reflection **{r.refl_late_pct}%**, E_end **{r.E_end_Ha} Ha**" if r is not None else "")
        C.append(new_markdown_cell(
            f"## {s['label']}\n"
            f"CAP η={s['eta']}, region {s['region']}, slab–CAP gap {s['gap']} Bohr{stat}."))
        C.append(new_markdown_cell(f"**Total density n(x,z,t)** (linear | log)\n\n![dens]({dens})"))
        C.append(new_markdown_cell(f"**Induced Δn = n(t)−n(0)** (linear | symlog — the wake)\n\n![delta]({delt})"))
        C.append(new_markdown_cell(f"**Energy & electron number vs time**\n\n![en]({rel(en)})"))

    C.append(new_markdown_cell(
        "## Takeaway\n"
        "- **No bath-eating problem:** total N(t) holds at 52 for every CAP → the CAP absorbs exactly the "
        "1-e⁻ projectile, nothing else. The suspected problem is absent.\n"
        "- **Keep the current CAP (η=−1.0, [25,40]):** complete absorption, 0.13% reflection, no bath loss.\n"
        "- **Do not weaken:** η=−0.4 is the only bad setting (4.1% edge reflection + 1.25% residue).\n"
        "- **Stopping is not the WP energy method:** E_total depends on η → use the classical twin / WP momentum.\n"
        "- Wider gap (R1) protects the wake marginally (0.33→0.22 e⁻ under CAP) at a small adiabaticity cost "
        "(1.15% reflection); adopt only if a wake analysis shows |z|=25 clipping matters."))
    C.append(new_code_cell(
        "import pandas as pd\n"
        "df = pd.read_csv('cap_study_summary.csv'); df"))

    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": str(OUT)}})
    nbf = OUT / "cap_study_comparison.ipynb"
    nbformat.write(nb, nbf)
    print(f"wrote {nbf}  ({len(nb.cells)} cells)")
    print(f"data  {csv}")
    print(tab.to_string(index=False))


if __name__ == "__main__":
    build()
