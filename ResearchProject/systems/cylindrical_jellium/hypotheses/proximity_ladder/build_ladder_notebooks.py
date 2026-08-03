#!/usr/bin/env python3
"""Per-rung and cross-rung notebooks for the cylindrical proximity ladder.

Plan: docs/plans/cylindrical-proximity-ladder.md
Run AFTER `build_ladder_figures.py`, which produces `figures/<rung>/*.png`,
`figures/comparison/*.png` and `figures/ladder_summary.csv`.

WHAT EACH RUNG NOTEBOOK CONTAINS (the user's list, 2026-08-02)
  1. density-matrix GIF, inline and at the TOP (.claude/rules/notebook-density-gif.md)
  2. interaction energies, classical and wavepacket
  3. projectile position / trajectory
  4. T1 (drift), T2 (total), and the var(p) term T2 - T1
  5. classical 1/2 m v^2
  6. classical vs wavepacket energy-loss definitions, overlaid
  7. wavepacket momentum loss at several times
  8. stopping power from T1, T2 and the classical definition, with uncertainties

The panels are already drawn by `build_ladder_figures.py` under shared axis limits
across every rung — that sharing is the point, so a reader flipping between rungs
sees a real change rather than an autoscale artefact. This module ASSEMBLES them,
and executes only light cells (the summary table), so a notebook build cannot fail
on a heavy recomputation the figure stage already did.

WHY THE FIGURES ARE EMBEDDED AS FILES BUT THE GIF IS EMBEDDED AS BYTES
The density GIF must ANIMATE when the notebook is reopened without its sidecar,
so it is base64'd into the stored output via IPython.display.Image
(.claude/rules/notebook-density-gif.md). Static PNGs are referenced by relative
path, which keeps the .ipynb small and lets a rebuilt figure appear without
regenerating the notebook.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
FIGS = HERE / "figures"

RUNG_TITLE = {
    "r10": "R_in = 10 Bohr = 2.5 sigma_WP — the weak (channeling) limit",
    "r08": "R_in = 8 Bohr = 2.0 sigma_WP",
    "r06": "R_in = 6 Bohr = 1.5 sigma_WP",
    "r04": "R_in = 4 Bohr = 1.0 sigma_WP — 37 % of the packet starts inside the wall",
    "r00": "Filled cylinder — the projectile is immersed from step 0",
    "r04n160": "R_in = 4 Bohr at N = 160 — the same-N control",
}

# Section order and the figure-stem prefixes that belong to each. Prefixes are
# matched against the START of the file name, so new panels with the same prefix
# are picked up without editing this table.
SECTIONS = [
    ("Trajectory and projectile position",
     "How far the projectile actually travelled, and whether the wavepacket's two "
     "position definitions (density centroid vs integral of <p_z> dt) still agree. "
     "Where they part company is where 'the wavepacket has a trajectory' stops "
     "being true — a diagnostic, not an error.",
     ("02", "03")),
    ("Kinetic-energy channels: T1, T2 and var(p)",
     "T1 = |<p>|^2/2m is the DRIFT channel; T2 = <p^2>/2m adds var(p)/2m. "
     "T2 - T1 is exactly 1.2755 eV for a free packet, so its drift is a direct "
     "readout of interaction with the bath rather than of dispersion.",
     ("04",)),
    ("Energy loss: classical and wavepacket definitions together",
     "The classical 1/2 m v^2 loss alongside the wavepacket's T1 and T2 losses. "
     "These are different definitions of the same physical question, and the gap "
     "between them is the result, not a discrepancy to be averaged away.",
     ("05",)),
    ("Pairwise interaction energies",
     "The P/S/B decomposition (.claude/rules/decomposed-interaction-energies.md). "
     "E_PP is the wavepacket's self-Hartree and has NO classical counterpart; "
     "E_PB is likewise zero for a rigid charge in a z-uniform background. Those "
     "absences are physics, not missing data.",
     ("09",)),
    ("Momentum-space: what the projectile lost, and to where",
     "The wavepacket's momentum distribution and its change from t = 0. k_z is on "
     "the exact FFT grid (dk_z = 2 pi / 60 = 0.105) — that is a hard resolution "
     "limit, not a binning choice.",
     ("10", "11", "12")),
    ("Stopping power from all three estimators",
     "S from T1, T2 and the classical dE/ds, with fit uncertainties. Per "
     ".claude/rules/light-projectile-stopping.md this is the INITIAL DRAG over the "
     "early near-constant-velocity window; a full-run regression would average S "
     "over every velocity from v0 down to whatever the projectile ends at.",
     ("13",)),
]


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(text)


def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(src)


def figures_for(rung_dir: Path, prefixes: tuple[str, ...]) -> list[Path]:
    if not rung_dir.is_dir():
        return []
    return sorted(p for p in rung_dir.glob("*.png")
                  if p.name.startswith(prefixes))


def gif_cell(rung: str) -> list[nbf.NotebookNode]:
    """The mandatory density-matrix GIF, inline and near the top.

    Emits the cell even when no GIF exists — it then prints why, which is the
    signal that the run needs re-running with density frames enabled, rather
    than a silently missing section.
    """
    return [
        md("## 1. Density evolution — read this first\n\n"
           "The real-space density n(r,t) in the propagation x-z plane (mid-y "
           "slice), for the classical projectile, the wavepacket, and their "
           "difference. This is the most direct picture of what the quantum "
           "projectile does differently: dispersion, reflection, capture. A "
           "static carpet compresses time onto one axis and hides it."),
        code(
            "from pathlib import Path\n"
            "from IPython.display import Image, display, Markdown\n"
            "\n"
            "# The density MATRIX: three representations x three kinds\n"
            "#   representation: classical | wp | wp_minus_cl\n"
            "#   kind:           density n(x,z,t) | induced n(t)-n(0) | instantaneous n(t)-n(t-dt)\n"
            "# Mid-y xz slice, PHYSICAL order (never fftshifted — vti-coordinate-mapping rule).\n"
            "# Image(filename=...) base64-embeds the bytes into the stored output, so these\n"
            "# animate on reopen even though the source VTIs have been reclaimed.\n"
            f"d = Path('figures/{rung}')\n"
            "order = ['classical', 'wp', 'wp_minus_cl']\n"
            "kinds = ['density', 'induced', 'instantaneous']\n"
            "shown = 0\n"
            "for kind in kinds:\n"
            "    for rep in order:\n"
            "        g = d / f'matrix_{rep}_{kind}.gif'\n"
            "        if g.is_file():\n"
            "            display(Markdown(f'**{rep} — {kind}**'))\n"
            "            display(Image(filename=str(g)))\n"
            "            shown += 1\n"
            "if shown == 0:\n"
            "    print(f'no density GIFs under {d}')\n"
            "else:\n"
            "    print(f'{shown} density-matrix animations')\n"
        ),
    ]


def build_rung_notebook(rung: str, out: Path) -> int:
    rung_dir = FIGS / rung
    nb = nbf.v4.new_notebook()
    cells = [
        md(f"# Proximity ladder — rung `{rung}`\n\n"
           f"**{RUNG_TITLE.get(rung, rung)}**\n\n"
           "Plan: `docs/plans/cylindrical-proximity-ladder.md`\n\n"
           "A matched classical / wavepacket pair fired on-axis down a periodic "
           "r_s = 3 jellium tube at 50 eV (v/v_F = 3.00), sigma_WP = 4 Bohr. Every "
           "rung of the ladder shares the cell, grid, projectile and density; only "
           "the bore radius changes.\n\n"
           "> **Axis limits are shared across every rung.** A difference you see "
           "between two rungs' figures is a real difference, not an autoscale.\n\n"
           "> **What this rung cannot show.** The projectile's Gaussian form factor "
           "is 0.018 at q = 1 and 3e-26 at q = 2v_0, so it couples to the collective "
           "response and essentially nothing else. The ladder spans weak-collective "
           "to strong-collective coupling; the electron-hole pair channel is a "
           "sigma_WP axis, not an R_in axis."),
    ]
    cells += gif_cell(rung)

    n_fig = 0
    for i, (title, blurb, prefixes) in enumerate(SECTIONS, start=2):
        figs = figures_for(rung_dir, prefixes)
        cells.append(md(f"## {i}. {title}\n\n{blurb}"))
        if not figs:
            cells.append(md(f"_No figures matching {prefixes} under `{rung_dir}`._"))
            continue
        n_fig += len(figs)
        body = "\n".join(
            f"![{p.stem}](figures/{rung}/{p.name})\n\n*{p.stem}*\n" for p in figs)
        cells.append(md(body))

    # the rung's row of the ladder table, executed so it reflects current data
    cells.append(md("## 8. This rung in the ladder\n\n"
                    "The row this run contributes to the cross-rung comparison."))
    cells.append(code(
        "import pandas as pd\n"
        "t = pd.read_csv('figures/ladder_summary.csv')\n"
        f"row = t[t.rung == '{rung}']\n"
        "display(row.T if len(row) else 'this rung is not in the summary table')\n"
    ))

    nb["cells"] = cells
    nb.metadata.kernelspec = {"name": "python3", "display_name": "Python 3",
                              "language": "python"}
    out.write_text(nbf.writes(nb))
    return n_fig


def build_comparison_notebook(rungs: list[str], out: Path) -> int:
    figs = sorted((FIGS / "comparison").glob("*.png")) if (FIGS / "comparison").is_dir() else []
    nb = nbf.v4.new_notebook()
    cells = [
        md("# Proximity ladder — cross-rung comparison\n\n"
           "Plan: `docs/plans/cylindrical-proximity-ladder.md`\n\n"
           "How the classical/quantum agreement in the stopping power behaves as "
           "the jellium wall is brought in from 2.5 sigma_WP to a filled cylinder, "
           "at fixed r_s = 3, fixed 50 eV projectile and fixed sigma_WP = 4.\n\n"
           "## How to read the x-axis\n\n"
           "Rungs are labelled by R_in/sigma_WP because that is legible, but the "
           "coupling is **exponential** in R_in, not linear: the wavepacket charge "
           "inside the wall runs 0.19 %, 1.8 %, 10.5 %, 37 %, 100 % across the "
           "ladder at t = 0. It is also **time-dependent** — the packet spreads "
           "from sigma_d = 2.83 to 6.01 Bohr over the run, so by t = 30 those same "
           "rungs span only 25 % to 100 %. **The rungs are distinct only early**, "
           "which is exactly why S is fitted in an early window keyed to a common "
           "measured coupling rather than a common time.\n\n"
           "## Three things move together, by construction\n\n"
           "Shrinking R_in at fixed r_s changes proximity, electron count "
           "(160 -> 326) and the target's mode spectrum (thin annulus with two "
           "coupled surfaces -> solid nanowire) simultaneously. In this geometry "
           "they are inseparable, because the electrons added *are* the close ones. "
           "The `r04n160` same-N control is what separates 'the wall is closer' "
           "from 'there is more wall'."),
    ]
    if figs:
        cells.append(md("\n".join(
            f"![{p.stem}](figures/comparison/{p.name})\n\n*{p.stem}*\n" for p in figs)))
    else:
        cells.append(md("_No comparison figures yet — run `build_ladder_figures.py`._"))

    cells.append(md("## The ladder table"))
    cells.append(code(
        "import pandas as pd\n"
        "t = pd.read_csv('figures/ladder_summary.csv')\n"
        "display(t)\n"
    ))
    cells.append(md(
        "## Coverage\n\n"
        "Which rungs made it into this comparison, and why any are missing. "
        "A silently short ladder would read as a completed one."))
    cells.append(code(
        "import json, pathlib\n"
        "m = json.loads(pathlib.Path('figures/manifest.json').read_text())\n"
        "print('missing    :', m.get('missing', {}))\n"
        "print('load failed:', m.get('load_failed', {}))\n"
    ))

    nb["cells"] = cells
    nb.metadata.kernelspec = {"name": "python3", "display_name": "Python 3",
                              "language": "python"}
    out.write_text(nbf.writes(nb))
    return len(figs)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rungs", default="r10,r08,r06,r04,r00")
    a = ap.parse_args()

    if not FIGS.is_dir():
        print(f"No figures directory at {FIGS} — run build_ladder_figures.py first.")
        return 1

    manifest = {}
    try:
        manifest = json.loads((FIGS / "manifest.json").read_text())
    except Exception:                                    # noqa: BLE001
        pass

    built = []
    for rung in [r.strip() for r in a.rungs.split(",") if r.strip()]:
        if not (FIGS / rung).is_dir():
            print(f"[skip] {rung}: no figures/{rung}/ — the rung is missing or failed")
            continue
        out = HERE / f"rung_{rung}.ipynb"
        n = build_rung_notebook(rung, out)
        built.append(rung)
        print(f"wrote {out.name}  ({n} figures)")

    out = HERE / "ladder_comparison.ipynb"
    n = build_comparison_notebook(built, out)
    print(f"wrote {out.name}  ({n} comparison figures)")

    if manifest.get("missing") or manifest.get("load_failed"):
        print("\nNOT a complete ladder:")
        for k in ("missing", "load_failed"):
            for tag, why in (manifest.get(k) or {}).items():
                print(f"  {tag}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
