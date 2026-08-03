#!/usr/bin/env python3
"""THE PHASE NOTEBOOK for the cylindrical proximity ladder.

Plan: docs/plans/cylindrical-proximity-ladder.md
Run AFTER build_ladder_figures.py (needs figures/ladder_summary.csv).

WHAT THIS IS, AND WHY IT IS SEPARATE FROM THE RUN NOTEBOOKS
-----------------------------------------------------------
The per-rung notebooks (`rung_*.ipynb`) show ONE run in depth. This notebook is
the SYNTHESIS: it answers the question the whole campaign was built to answer —

    how does a KS-orbital definition of electronic stopping power compare with
    the classical dE/ds definition, as the system moves from WEAK to STRONG
    interaction, where weak->strong is parameterised by the rung radius R_in?

Everything here is computed live from `figures/ladder_summary.csv` and the raw
per-step frames, so re-running after a rung is added or re-run updates the
conclusions rather than leaving stale prose. The narrative cells state what the
numbers mean; the code cells produce the numbers. If the two ever disagree, the
code is right and the prose needs editing — the notebook is built to make that
disagreement visible rather than to hide it.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent


def md(t): return nbf.v4.new_markdown_cell(t)
def code(s): return nbf.v4.new_code_cell(s)


SETUP = '''\
%matplotlib inline
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

HERE = Path.cwd()
SYS  = HERE.parents[1]                       # systems/cylindrical_jellium
sys.path.insert(0, str(SYS / "hypotheses/channeling_twin"))
from inqview.visualisation import style
style.apply_theme()

T = pd.read_csv("figures/ladder_summary.csv")
ORDER = ["r10", "r08", "r06", "r04", "r00"]
T["k"] = T.rung.apply(lambda r: ORDER.index(r) if r in ORDER else 99)
T = T.sort_values("k").drop(columns="k").reset_index(drop=True)

# the estimator columns the fit stage produced
SW = [c for c in T.columns if c.startswith("S_wp[")]
PRIMARY = next((c for c in SW if "T1" in c), SW[0])       # the drift channel
PRIM_CL = PRIMARY.replace("S_wp[", "S_cl[")
T["ratio_T1"] = T[PRIMARY] / T[PRIM_CL]
print(f"{len(T)} rungs: {list(T.rung)}")
print(f"primary estimator column: {PRIMARY}")
'''

LADDER_TABLE = '''\
view = T[["rung", "R_in_over_sigma", "fw_fit_mean", "fw_fit_drift",
          "frac_loss_cl", PRIMARY, PRIM_CL, "ratio_T1"]].copy()
view.columns = ["rung", "R_in/sigma", "<f_wall> over fit", "coupling drift",
                "frac. energy lost", "S_wp (T1)", "S_classical", "S_wp/S_cl"]
view.round({"R_in/sigma": 2, "<f_wall> over fit": 3, "coupling drift": 2,
            "frac. energy lost": 3, "S_wp (T1)": 4, "S_classical": 4, "S_wp/S_cl": 3})
'''

WEAK_STRONG = '''\
fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.3), constrained_layout=True)
x = np.arange(len(T)); lab = T.rung

ax[0].plot(x, 100*T.frac_loss_cl, "o-", color="tab:blue", label="classical")
ax[0].plot(x, 100*T.fw_fit_mean, "s--", color="tab:green", label=r"$\\langle f_\\mathrm{wall}\\rangle$")
ax[0].set_xticks(x); ax[0].set_xticklabels(lab)
ax[0].set_ylabel("per cent"); ax[0].set_xlabel("rung:  weak  -->  strong")
ax[0].set_title("the coupling axis is real", fontsize=9)
ax[0].legend(frameon=False, fontsize=8); ax[0].grid(alpha=0.3)

ax[1].plot(x, T[PRIM_CL], "o-", color="tab:blue", label=r"$S_\\mathrm{classical}$")
ax[1].plot(x, T[PRIMARY], "s-", color="tab:red", label=r"$S_\\mathrm{WP}$ ($T_1$ drift)")
ax[1].set_xticks(x); ax[1].set_xticklabels(lab)
ax[1].set_ylabel("S  (eV/Bohr)"); ax[1].set_xlabel("rung:  weak  -->  strong")
ax[1].set_title("the two definitions diverge", fontsize=9)
ax[1].legend(frameon=False, fontsize=8); ax[1].grid(alpha=0.3)
plt.show()

g_cl = T[PRIM_CL].iloc[-1] / T[PRIM_CL].iloc[0]
g_wp = T[PRIMARY].iloc[-1] / T[PRIMARY].iloc[0]
print(f"across the ladder:  S_classical x{g_cl:.2f}   S_wp(T1) x{g_wp:.2f}")
print(f"the classical definition responds {g_cl/g_wp:.2f}x more strongly to the coupling")
'''

RATIO = '''\
fig, ax = plt.subplots(figsize=(4.6, 3.4), constrained_layout=True)
ax.errorbar(T.fw_fit_mean, T.ratio_T1,
            xerr=[T.fw_fit_mean - T.fw_fit_lo, T.fw_fit_hi - T.fw_fit_mean],
            fmt="o-", color="tab:purple", lw=1.4, ms=6, capsize=3)
for _, r in T.iterrows():
    ax.annotate(r.rung, (r.fw_fit_mean, r.ratio_T1),
                textcoords="offset points", xytext=(7, 5), fontsize=8)
ax.axhline(1.0, color="k", lw=0.8, ls=":")
ax.set_xlabel(r"$\\langle f_\\mathrm{wall}\\rangle$ over the fit window  (measured coupling)")
ax.set_ylabel(r"$S_\\mathrm{WP}/S_\\mathrm{classical}$")
ax.set_title("bars = coupling drift within the window", fontsize=9)
ax.grid(alpha=0.3)
plt.show()

d = T.ratio_T1.to_numpy()
print("ratio:", "  ".join(f"{v:.3f}" for v in d))
print("step :", "  ".join(f"{v:+.3f}" for v in np.diff(d)))
if len(d) >= 4 and abs(np.diff(d)[-1]) < 0.5*abs(np.diff(d)[-2]):
    print(f"\\n-> SATURATION: the last step ({np.diff(d)[-1]:+.3f}) is less than half the")
    print(f"   previous ({np.diff(d)[-2]:+.3f}). The under-reporting bottoms out near {d[-1]:.3f},")
    print(f"   i.e. the KS drift estimator's error is BOUNDED at ~{100*(1-d[-1]):.0f} %, not unbounded.")
'''

MECHANISM = '''\
# Where the missing stopping power goes. var(p) is EXACTLY conserved under free
# evolution, so any growth is interaction and nothing else.
import re
vp = {}
for tag in T.rung:
    p = (SYS/"scripts/channeling_twin/wp/results/wp/run_summary.txt" if tag == "r10"
         else SYS/f"scripts/proximity_ladder/wp/results/{tag}/run_summary.txt")
    if p.is_file():
        m = re.search(r"var_pz_growth_pct\\s*=\\s*([-\\d.eE+]+)", p.read_text())
        if m: vp[tag] = float(m.group(1))
T["var_p_growth_pct"] = T.rung.map(vp)

fig, ax = plt.subplots(figsize=(4.6, 3.4), constrained_layout=True)
sh = 100*(1 - T.ratio_T1)
ax.plot(T.var_p_growth_pct, sh, "o-", color="tab:orange", lw=1.4, ms=6)
for _, r in T.iterrows():
    ax.annotate(r.rung, (r.var_p_growth_pct, 100*(1-r.ratio_T1)),
                textcoords="offset points", xytext=(7, -9), fontsize=8)
ax.set_xlabel(r"var$(p_z)$ growth over the run  (%)")
ax.set_ylabel(r"$T_1$ shortfall  $100(1 - S_\\mathrm{WP}/S_\\mathrm{cl})$  (%)")
ax.set_title("the shortfall tracks the momentum-spread growth", fontsize=9)
ax.grid(alpha=0.3)
plt.show()
print(T[["rung", "var_p_growth_pct", "ratio_T1"]].assign(
      shortfall_pct=sh.round(1)).to_string(index=False))
'''

ESTIMATORS = '''\
# All estimators/windows the fit stage produced, as ratios. T2 = drift + var(p)/2m.
rows = []
for c in SW:
    cc = c.replace("S_wp[", "S_cl[")
    if cc in T.columns:
        rows.append(pd.Series((T[c]/T[cc]).to_numpy(), index=T.rung, name=c[6:-1].strip()))
R = pd.DataFrame(rows)
display(R.round(3))
print()
# A saturating trend has ONE tiny reversal; an erratic one swings. Judge by the
# largest reversal RELATIVE TO the total range, not by a strict monotonicity test
# — otherwise a plateau (the physically interesting case) is mislabelled erratic.
for name, r in R.iterrows():
    v = r.to_numpy(); d = np.diff(v); rng = v.max() - v.min()
    up = d[d > 0].sum() if (d > 0).any() else 0.0
    frac = up / rng if rng > 0 else 0.0
    if frac < 0.05:   verdict = "monotonic"
    elif frac < 0.20: verdict = f"monotonic then SATURATING (reversal {100*frac:.1f}% of range)"
    else:             verdict = f"ERRATIC (reversals {100*frac:.0f}% of range)"
    print(f"{name:16s} range {v.min():.3f}-{v.max():.3f}   {verdict}")
'''


def build(out: Path) -> None:
    nb = nbf.v4.new_notebook()
    c = [
        md("# Proximity ladder — phase analysis\n\n"
           "## The question\n\n"
           "**How does a Kohn–Sham orbital definition of electronic stopping power "
           "compare with the classical ΔE/Δs definition, as the system moves from "
           "weak to strong interaction?**\n\n"
           "Weak → strong is parameterised by the tube bore radius `R_in`: the "
           "wavepacket is fired on-axis down a periodic r_s = 3 jellium tube, and "
           "the wall is brought inward from 2.5 σ_WP to a filled cylinder. "
           "Everything else is held fixed — r_s = 3.000000, projectile 50 eV "
           "(v/v_F = 3.00), σ_WP = 4 Bohr, 40×40×60 Bohr cell, dx = 0.5, "
           "1500 steps × dt 0.02.\n\n"
           "| rung | R_in | R_in/σ | N_e | shape |\n|---|---|---|---|---|\n"
           "| r10 | 10 | 2.5 | 160 | annulus |\n| r08 | 8 | 2.0 | 220 | annulus |\n"
           "| r06 | 6 | 1.5 | 266 | annulus |\n| r04 | 4 | 1.0 | 300 | annulus |\n"
           "| r00 | 0 | — | 326 | **filled cylinder** |\n\n"
           "Each rung is a matched **twin**: a wavepacket (occupied KS orbital) and a "
           "classical Gaussian charge at σ_pot = σ_WP/√2, identical in every physical "
           "parameter. A difference in S between the halves is a quantum effect and "
           "nothing else.\n\n"
           "### The two estimators\n\n"
           "$$T_1 = \\frac{|\\langle \\mathbf{p}\\rangle|^2}{2m} \\quad\\text{(drift)}"
           "\\qquad T_2 = T_1 + \\frac{\\mathrm{var}(\\mathbf{p})}{2m} = "
           "\\frac{\\langle p^2\\rangle}{2m} \\quad\\text{(total)}$$\n\n"
           "$S = -\\mathrm{d}T/\\mathrm{d}s$ for the wavepacket, "
           "$S = +\\mathrm{d}E_\\mathrm{total}/\\mathrm{d}s$ for the classical twin. "
           "$\\mathrm{var}(\\mathbf{p})$ is **exactly conserved under free evolution**, "
           "so any growth in it is interaction and nothing else — that is what makes "
           "$T_2 - T_1$ a direct readout of the contamination."),
        code(SETUP),

        md("## 1. The ladder\n\n"
           "`<f_wall>` is the **measured** fraction of wavepacket charge inside the "
           "jellium over the fit window — the coupling coordinate. It is used in "
           "preference to the nominal `R_in` because the bore is not empty: rung r10's "
           "ground state already has 16 of its 160 electrons inside it, so *distance to "
           "the background edge* is not *distance to the electrons*."),
        code(LADDER_TABLE),

        md("## 2. Is the weak → strong axis real?\n\n"
           "Before comparing estimators, confirm the ladder actually spans a coupling "
           "range. If classical energy loss and measured occupancy both rise strongly "
           "and monotonically, the axis is doing its job."),
        code(WEAK_STRONG),

        md("## 3. The headline: how the two definitions compare\n\n"
           "Plotted against the **measured** coupling, with horizontal bars showing how "
           "much the coupling drifts *within* each fit window. That drift is intrinsic: "
           "the packet spreads into the wall as it flies, so each S is an average over a "
           "coupling range rather than a value at a point. Constraining the range to "
           "1.5× would need a ~120-step window, too short for a stable fit. The bars are "
           "the honest representation of that, and they narrow toward strong coupling — "
           "so the strong end is the *better*-conditioned part of the ladder."),
        code(RATIO),

        md("## 4. Mechanism — where does the missing stopping power go?\n\n"
           "If the drift channel under-reports, the momentum transfer must be going "
           "somewhere else. The candidate is the packet's own **momentum spread**: "
           "energy that broadens $|\\psi(p)|^2$ instead of decelerating its centroid. "
           "Because var$(p)$ is exactly conserved for a free packet, its growth is an "
           "unambiguous interaction signal."),
        code(MECHANISM),

        md("## 5. All estimators and windows\n\n"
           "$T_1$ is the drift channel; $T_2$ adds var$(p)/2m$. A useful estimator "
           "should vary smoothly and monotonically with coupling — an erratic one is "
           "unusable regardless of how close any single value sits to 1."),
        code(ESTIMATORS),

        md("## 6. What this establishes\n\n"
           "Read the printed output above rather than trusting this prose — the cells "
           "recompute from the CSV, so if a rung is re-run the numbers move and this "
           "text may need updating.\n\n"
           "**The expected findings, as of the 5-rung campaign:**\n\n"
           "1. **The axis works.** Fractional energy loss ~10 % → ~42 %, measured "
           "coupling ~0.12 → ~0.99. Genuinely perturbative to strongly non-linear, at "
           "fixed density, energy, packet width, cell and grid.\n"
           "2. **The definitions diverge systematically.** Classical S responds ~2× more "
           "strongly to the coupling than the KS drift channel does. The ratio falls "
           "monotonically from ~0.80 to ~0.41.\n"
           "3. **The error is BOUNDED.** The ratio saturates: the last step is less than "
           "half the previous one, so the drift estimator under-reports by ~60 % once "
           "immersed and stops degrading. That is a more useful statement than an "
           "unbounded decline — it means a KS drift measurement can be *corrected* "
           "rather than merely distrusted.\n"
           "4. **The mechanism is identified.** The shortfall tracks var$(p)$ growth "
           "(~45 % → ~115 %). The missing stopping power is momentum going into packet "
           "width rather than centroid deceleration.\n"
           "5. **$T_2$ is unusable** — erratic rather than smoothly wrong.\n\n"
           "### The limit this campaign does NOT cross\n\n"
           "The projectile's Gaussian form factor $\\exp(-q^2\\sigma_\\mathrm{pot}^2/2)$ "
           "is 0.018 at $q = 1$ and $3\\times10^{-26}$ at $q = 2v_0 = 3.83$. The plasmon "
           "pole sits at $q_\\mathrm{min} = \\omega_p/v = 0.174$; the electron–hole "
           "continuum runs to $q = 2v$. **This projectile couples to the collective "
           "response and essentially nothing else, at every rung.** Shrinking $R_\\mathrm{in}$ "
           "scales how much medium responds; it does not harden the projectile.\n\n"
           "So this ladder spans **weak-collective → strong-collective**. Reaching the "
           "pair channel is a **σ_WP axis**, not an $R_\\mathrm{in}$ axis. Every "
           "conclusion above is a statement about collective stopping.\n\n"
           "### Other caveats\n\n"
           "* Each S is an average over a coupling range (bars in §3), widest at the "
           "weak end (~4.9×) and narrowest at the strong end (~1.04×).\n"
           "* Three things move together along the ladder — proximity, electron count "
           "(160 → 326) and the target's mode spectrum (thin annulus → solid nanowire). "
           "They are inseparable in this geometry because the electrons added *are* the "
           "close ones. The `r04n160` same-N control (not yet run) is what would "
           "disentangle them.\n"
           "* Pauli blocking is negligible throughout (`max_overlap` ≲ 5e-7 even fully "
           "immersed) because $k_0 = 1.92$ sits 7.2 σ_p above $k_F = 0.64$. The "
           "wavepacket/classical difference is therefore *not* an exclusion effect."),
    ]
    nb["cells"] = c
    nb.metadata.kernelspec = {"name": "python3", "display_name": "Python 3",
                              "language": "python"}
    out.write_text(nbf.writes(nb))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "phase_analysis.ipynb"))
    a = ap.parse_args()
    build(Path(a.out))
    print(f"wrote {a.out}")
