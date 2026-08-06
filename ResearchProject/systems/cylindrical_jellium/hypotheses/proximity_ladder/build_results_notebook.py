#!/usr/bin/env python3
"""Emit `results.ipynb` — the plain-language answer to "what did the ladder show?".

Plan:     docs/plans/cylindrical-proximity-ladder.md
Handover: docs/handovers/cylindrical-proximity-ladder.md

SCOPE
-----
There are now three notebooks over this campaign and they are deliberately not
interchangeable:

  rung_<tag>.ipynb   one run pair in full detail (24 panels + density GIFs)
  phase_analysis.ipynb   the campaign's internal synthesis: gates, costs,
                         validation records, the estimator classifier
  results.ipynb (this)   the READING: S(T1) and S(T2) against interaction
                         strength, five figures, one table, stated in words

This one is written to be handed to someone who did not run the campaign. It
therefore leads with the pictures, states the window and why, and puts the
caveats where they cannot be skipped rather than in a footnote.

DENSITY GIF (.claude/rules/notebook-density-gif.md)
--------------------------------------------------
The rule requires an animated n(r,t) displayed inline near the top. The per-rung
GIFs already exist under figures/<rung>/ (the raw VTIs they were rendered from
have since been deleted to reclaim quota), so this notebook EMBEDS the two that
carry the comparison — weakest and strongest coupling, WP minus classical — as
the visual intuition section, rather than re-rendering.
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
import json
from pathlib import Path
import numpy as np, pandas as pd
from IPython.display import Image, display

HERE = Path.cwd()
FIG  = HERE / "figures/results"

T = pd.read_csv(FIG / "results_summary.csv")
W = json.loads((FIG / "window.json").read_text())

# one human label per rung, used by every cell below
NICE = {"r10": "R_in = 2.5 sigma", "r08": "R_in = 2.0 sigma",
        "r06": "R_in = 1.5 sigma", "r04": "R_in = 1.0 sigma",
        "r00": "filled tube"}

print(f"{len(T)} rungs: {list(T.rung)}")
print(f"fit window: t = {W['t0_au']:.0f} - {W['t1_au']:.0f} a.u.")
print(f"rationale: {W['rationale']}")
'''

TABLE = '''\
tab = pd.DataFrame({
    "bore R_in":      T.rung.map(NICE),
    "coupling f_wall": T.f_wall_mean.round(2),
    "S(T1)":          T.S_T1.round(3),
    "S(T2)":          T.S_T2.round(3),
    "S classical":    T.S_CL.round(3),
    "S(T1)/S_cl":     T.ratio_T1.round(2),
    "S(T2)/S_cl":     T.ratio_T2.round(2),
    "r2 of T2 fit":   T.r2_T2.round(2),
})
print("stopping power in eV/Bohr, fitted over t = 11-20 a.u.\\n")
display(tab)
print("\\nvelocity criterion over the window (rule: stay above 0.85 v0):")
print(f"  classical   min v/v0 = {T.v_over_v0_min_cl.min():.3f}")
print(f"  wavepacket  min v/v0 = {T.v_over_v0_min_wp.min():.3f}")
'''

NUMBERS = '''\
w = T[T.f_wall_mean < 0.2].iloc[0]          # weakest coupling
s = T[T.f_wall_mean > 0.9]                  # strongest two
print(f"weak coupling ({NICE[w.rung]}):  S(T1)/S_cl = {w.ratio_T1:.2f}"
      f"   S(T2)/S_cl = {w.ratio_T2:.2f}   -> they differ by {w.ratio_T1/w.ratio_T2:.1f}x")
print(f"strong coupling:                S(T1)/S_cl = {s.ratio_T1.mean():.2f}"
      f"   S(T2)/S_cl = {s.ratio_T2.mean():.2f}   -> they agree to "
      f"{100*abs(s.ratio_T1.mean()-s.ratio_T2.mean())/s.ratio_T1.mean():.0f}%")
print()
print(f"classical S rises {T.S_CL.iloc[-1]/T.S_CL.iloc[0]:.2f}x across the ladder")
print(f"S(T1)     rises {T.S_T1.iloc[-1]/T.S_T1.iloc[0]:.2f}x")
print(f"S(T2)     rises {T.S_T2.iloc[-1]/T.S_T2.iloc[0]:.1f}x  (from a near-zero base)")
print()
frac = 100 * T.var_growth_ev / T.dE_T1_ev
for r, f in zip(T.rung, frac):
    print(f"  {NICE[r]:>10s}:  {f:5.1f}% of the drift loss went into momentum spread")
'''


def build(out: Path) -> None:
    nb = nbf.v4.new_notebook()
    c = [
        md("# Proximity ladder — stopping power vs interaction strength\n\n"
           "**The question.** A 50 eV electron flies down the axis of a cylindrical "
           "jellium tube ($r_s = 3$). We bring the tube wall inwards in five steps — "
           "bore radius $R_{\\rm in} = 2.5\\sigma_{\\rm WP} \\to 2.0 \\to 1.5 \\to 1.0 \\to$ "
           "filled — and ask how the **stopping power** behaves. Everything else is held "
           "fixed: $\\sigma_{\\rm WP} = 4$, $v_0 = 1.917$, jellium density, grid, timestep.\n\n"
           "**Two definitions of the projectile's kinetic energy** are in play, and they "
           "are not the same thing for a wavepacket:\n\n"
           "| | definition | what it counts |\n|---|---|---|\n"
           "| $T_1$ | $\\lvert\\langle p\\rangle\\rvert^2/2m$ | energy of the **drift** — the "
           "centroid's motion. The direct analogue of $\\tfrac12 mv^2$. |\n"
           "| $T_2$ | $\\langle p^2\\rangle/2m$ | the **total** kinetic energy, drift plus "
           "spread. This is the projectile's actual KE. |\n\n"
           "They differ by exactly $T_2 - T_1 = \\mathrm{var}(p)/2m$. Under free evolution "
           "$\\mathrm{var}(p)$ is *exactly conserved*, so any growth in it is unambiguously "
           "caused by the interaction — that is what makes this decomposition clean.\n\n"
           "Each wavepacket run has a **classical twin**: a point projectile with the same "
           "energy and a Gaussian charge cloud of matched width "
           "($\\sigma_{\\rm pot} = \\sigma_{\\rm WP}/\\sqrt2$), in the same tube."),
        code(SETUP),

        md("## Visual intuition — what the quantum projectile does differently\n\n"
           "Induced density, wavepacket **minus** classical, mid-$y$ slice of the $xz$ "
           "propagation plane. Where the two representations agree this is black; the "
           "structure that survives is the quantum-specific response.\n\n"
           "Left-hand case is the weakest coupling in the ladder, right-hand the strongest."),
        code('for tag in ("r10", "r00"):\n'
             '    p = HERE / f"figures/{tag}/matrix_wp_minus_cl_induced.gif"\n'
             '    print(f"===== {tag}:  {NICE[tag]} =====")\n'
             '    if p.is_file():\n'
             '        display(Image(filename=str(p)))   # embeds the bytes in the .ipynb\n'
             '    else:\n'
             '        print("  (no GIF on disk)")'),

        md("## 1. The raw evidence: energy actually lost\n\n"
           "Same five runs, three different answers, because \"the projectile's kinetic "
           "energy\" is three different quantities. Note the **middle panel dips below "
           "zero** at intermediate coupling: measured by $T_2$, the wavepacket *gains* "
           "kinetic energy over the first ~30 Bohr."),
        code('display(Image(filename=str(FIG / "F1_energy_loss_vs_path.png")))'),

        md("## 2. Where a stopping power may legitimately be read off\n\n"
           "$S = -\\mathrm{d}E/\\mathrm{d}s$ at each instant. Two shaded regions:\n\n"
           "* **pink** — the wake is still building. At $r_s = 3$, "
           "$\\omega_p = \\sqrt{4\\pi n} = 0.333$ a.u., so the bath cannot polarise faster "
           "than a quarter period, $\\pi/2\\omega_p \\approx 4.7$ a.u. Every curve rises from "
           "zero over exactly this scale. **No stopping power exists here.**\n"
           "* **grey** — the window used for every number quoted below, "
           "$t \\in [11, 20]$ a.u. Its upper edge is set by the light-projectile rule "
           "(`.claude/rules/light-projectile-stopping.md`): the classical projectile in the "
           "filled tube drops below $0.85\\,v_0$ at $t = 20.6$.\n\n"
           "The classical panel shows a genuine **plateau** across this window — a "
           "well-defined steady-state stopping power. The $T_1$ panel peaks and then "
           "decays. The $T_2$ panel starts **negative**."),
        code('display(Image(filename=str(FIG / "F2_local_stopping_vs_time.png")))'),

        md("> **Why the earlier numbers looked erratic.** The windows inherited from the "
           "channeling twin were `T2 5-20` and `T2 21-30`. The first straddles the $T_2$ "
           "zero-crossing at every rung, the second sits in the late decay. Both averaged "
           "across a sign change, which is why $T_2$ appeared to jump around with no "
           "trend. That was a property of the window, not of $T_2$."),

        md("## 3. The numbers"),
        code(TABLE),

        md("## 4. Headline — stopping power against interaction strength\n\n"
           "$x$ is the **measured** coupling $f_{\\rm wall}$: the fraction of $|\\psi|^2$ "
           "actually sitting inside the jellium, averaged over the fit window, with its "
           "min–max drawn as horizontal bars. This is used rather than nominal $R_{\\rm in}$ "
           "because the bore is not empty — the weakest rung's ground state already has 16 "
           "of its 160 electrons inside it."),
        code('display(Image(filename=str(FIG / "F3_stopping_vs_coupling.png")))'),

        md("Both quantum definitions lie **below** the classical curve everywhere, and both "
           "rise with coupling — but they approach it from opposite directions. Dividing "
           "out the classical value makes that explicit:"),
        code('display(Image(filename=str(FIG / "F4_ratio_convergence.png")))'),
        code(NUMBERS),

        md("## 5. Why the two definitions differ — and why they stop differing\n\n"
           "The gap is exactly $\\mathrm{var}(p)/2m$: energy taken out of the directed "
           "motion but retained by the packet as momentum spread. Panel (b) is the whole "
           "explanation of the convergence in §4."),
        code('display(Image(filename=str(FIG / "F5_variance_mechanism.png")))'),

        md("## 6. Reading\n\n"
           "**a. Both definitions understate the classical stopping power, at every "
           "coupling.** The wavepacket is never stopped as hard as its classical twin.\n\n"
           "**b. At weak coupling the two definitions disagree violently** — 0.84 vs 0.13 "
           "of the classical value, a factor of ~6.7. $T_1$ nearly reproduces the classical "
           "answer; $T_2$ says the projectile has barely lost any kinetic energy at all. "
           "Both are correct about what they measure. The drift *is* decelerating, and the "
           "total kinetic energy *is* nearly unchanged, because the loss from the drift is "
           "being reinvested in momentum spread rather than given to the bath.\n\n"
           "**c. At strong coupling they converge** — 0.44/0.48 at $1.0\\sigma$ and "
           "0.46/0.45 filled, agreeing to within a few percent. Once the packet is fully "
           "immersed, the spread channel accounts for only ~2% of the drift loss, so \"drift "
           "energy\" and \"total kinetic energy\" become the same measurement. **This is "
           "where a wavepacket stopping power is a well-posed quantity.**\n\n"
           "**d. The converged value is ~0.46, not 1.** Fully immersed, with both "
           "definitions agreeing, the wavepacket still loses energy at less than half the "
           "classical rate. That residual is a genuine quantum deficit, not an artefact of "
           "which kinetic energy was chosen — the two independent definitions agreeing is "
           "what rules that out.\n\n"
           "**e. The mechanism is the momentum-spread channel, and it dies out.** The "
           "fraction of the drift loss diverted into $\\mathrm{var}(p)$ falls 54% → 43% → "
           "21% → 7% → 2% across the ladder, essentially linearly in the coupling. At weak "
           "coupling the projectile interacts with a thin sliver of its own tail, which "
           "shears the packet in momentum space more than it decelerates it. Fully immersed, "
           "the force acts on the whole packet at once and simply decelerates it."),

        md("## 7. What this does *not* establish\n\n"
           "**The ladder is a collective-coupling axis, not a route to the full stopping "
           "power.** The projectile's Gaussian form factor is $\\exp(-q^2\\sigma_{\\rm pot}^2/2)$: "
           "0.37 at $q = 0.5$, 0.018 at $q = 1$, and $3\\times10^{-26}$ at $q = 2v_0 = 3.83$. "
           "The electron–hole pair channel needs momentum transfers of order $2v_0$, so it is "
           "suppressed by 26 orders of magnitude **at every rung**. Bringing the wall in "
           "changes *how much* collective response the projectile sees; it does not open a new "
           "channel. Reaching the pair channel requires narrowing $\\sigma_{\\rm WP}$ — a "
           "different axis entirely.\n\n"
           "**Coupling is not constant within the fit window.** $f_{\\rm wall}$ drifts by "
           "2.9× at the weakest rung and 1.02× at the filled one, because the packet spreads "
           "transversally as it flies. Each $S$ is therefore an average over a *range* of "
           "couplings, and that range is widest exactly where the ladder is weakest. This is "
           "intrinsic to a spreading wavepacket: holding the coupling to within 1.5× would "
           "require $t < 2.4$ a.u., far too short for a stable fit. It is drawn honestly as "
           "the horizontal error bars, not hidden.\n\n"
           "**The $r08$ $T_2$ fit is the weakest number in the table** ($r^2 = 0.66$): its "
           "zero crossing at $t = 14.1$ falls inside the window. Its ratio (0.15) should be "
           "read as \"small\", not as a precise value.\n\n"
           "**Two variables move together along the ladder.** Bringing the wall in also adds "
           "electrons (160 → 326) to hold the density fixed, so \"closer wall\" and \"more "
           "wall\" are not separated. The `r04n160` same-N control (~2 GPU-h, not yet run) is "
           "the experiment that would disentangle them.\n\n"
           "**Pauli blocking is not the explanation.** $k_0 = 1.92$ sits 7.2 σ above "
           "$k_F = 0.64$ and the measured `max_overlap` is $\\le 5\\times10^{-7}$ even fully "
           "immersed."),

        md("---\n*Figures: `figures/results/` (`build_results_figures.py`). Per-run detail: "
           "`rung_<tag>.ipynb`. Campaign gates and costs: `phase_analysis.ipynb`.*"),
    ]
    nb["cells"] = c
    nb.metadata.kernelspec = {"name": "python3", "display_name": "Python 3",
                              "language": "python"}
    out.write_text(nbf.writes(nb))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "results.ipynb"))
    a = ap.parse_args()
    build(Path(a.out))
    print(f"wrote {a.out}")
