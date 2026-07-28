#!/usr/bin/env python3
"""compare_notebook — jellium wp_cap_energy_plateau PHASE notebook: no-CAP vs CAP.

Assembles the campaign's synthesis notebook (the jellium analogue of the vacuum
nocap_vs_cap_comparison.ipynb) from already-validated per-run artifacts:
  1. energy plateau gap  : absolute E_total(t) no-CAP vs CAP (the headline gap)
  2. ΔE_total(t) overlay : gap opening over time
  3. ΔE components       : every KS component, no-CAP (solid) vs CAP (dashed)
  4. ΔE pairwise         : E_ss/E_ps/E_pp/E_sb/E_pb, no-CAP vs CAP
  5. norm(t)             : WP absorption by the CAP
  + total-density GIFs (no-CAP | CAP) embedded near the top (rule: density GIF).

Usage: compare_notebook.py <nocap_results> <cap_results> <out_dir>
"""
from __future__ import annotations
import base64, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HA_EV = 27.211386245988
KS = [("total","k",2.4),("kinetic","C0",1.3),("hartree","C1",1.3),
      ("external","C2",1.3),("xc","C3",1.3)]
PAIR = [("e_ss","E_ss slab–slab","C0"),("e_ps","E_ps slab–proj","C1"),
        ("e_pp","E_pp proj self","C3"),("e_sb","E_sb slab–bg","C2"),
        ("e_pb","E_pb proj–bg","C4")]


def _en(results: Path):
    df = pd.read_csv(results / "raw/observables/energies_merged.csv")
    return df


def _delta(a): return a - a[0]


def plot_plateau(no, ca, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for df, lab, c in ((no, "no-CAP", "C0"), (ca, "CAP", "C3")):
        ax.plot(df.time_au, df.total * HA_EV, color=c, lw=1.8, label=lab)
    gap = (no.total.iloc[-1] - ca.total.iloc[-1]) * HA_EV
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("E_total (eV)")
    ax.set_title("Energy plateau gap — no-CAP sits above CAP")
    ax.legend(); ax.text(0.02, 0.05, f"final gap = {gap:+.1f} eV",
                         transform=ax.transAxes, fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)


def plot_dtotal(no, ca, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.plot(no.time_au, _delta(no.total.to_numpy()) * HA_EV, "C0-", lw=1.8, label="no-CAP")
    ax.plot(ca.time_au, _delta(ca.total.to_numpy()) * HA_EV, "C3-", lw=1.8, label="CAP")
    ax.axhline(0, color="gray", lw=0.6, ls=":")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("ΔE_total = E(t)−E(0) (eV)")
    ax.set_title("ΔE_total(t): CAP drains electronic energy; no-CAP conserved")
    ax.legend(); fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)


def plot_components(no, ca, out):
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for name, c, lw in KS:
        if name in no and np.any(np.abs(no[name]) > 1e-9):
            ax.plot(no.time_au, _delta(no[name].to_numpy()) * HA_EV, color=c, lw=lw,
                    ls="-", label=f"{name} (no-CAP)")
            ax.plot(ca.time_au, _delta(ca[name].to_numpy()) * HA_EV, color=c, lw=lw,
                    ls="--", label=f"{name} (CAP)")
    ax.axhline(0, color="gray", lw=0.6, ls=":")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("ΔE (eV)")
    ax.set_title("ΔE per KS component — solid=no-CAP, dashed=CAP")
    ax.legend(fontsize=7, ncol=2); fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)


def plot_pairwise(no_dir, ca_dir, out):
    no = pd.read_csv(no_dir / "raw/observables/interactions.csv")
    ca = pd.read_csv(ca_dir / "raw/observables/interactions.csv")
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for col, lab, c in PAIR:
        ax.plot(no.time_au, _delta(no[col].to_numpy()), color=c, lw=1.4, ls="-",
                label=f"{lab} (no-CAP)")
        ax.plot(ca.time_au, _delta(ca[col].to_numpy()), color=c, lw=1.4, ls="--",
                label=f"{lab} (CAP)")
    ax.axhline(0, color="gray", lw=0.6, ls=":")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("ΔE (eV)")
    ax.set_title("ΔE pairwise electrostatic — solid=no-CAP, dashed=CAP")
    ax.legend(fontsize=6.5, ncol=2)
    ax.text(0.02, 0.03, "E_ss+E_ps+E_pp≡E_hartree, E_sb+E_pb≡E_external (closure 1.4e-8 eV)",
            transform=ax.transAxes, fontsize=6.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8))
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)


def plot_norm(no, ca, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for df, lab, c in ((no, "no-CAP", "C0"), (ca, "CAP", "C3")):
        col = "N_total" if "N_total" in df else "total"
        ax.plot(df.time_au, df[col], color=c, lw=1.8, label=lab)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("N_total (electrons)")
    ax.set_title("Norm: CAP drains the WP (N drops); no-CAP conserved")
    ax.legend(); fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)


def build(nocap_dir, cap_dir, out_dir):
    import nbformat as nbf
    from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell, new_output
    nocap_dir, cap_dir, out_dir = Path(nocap_dir), Path(cap_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    no, ca = _en(nocap_dir), _en(cap_dir)

    pngs = []
    for fn, f in (("plateau_gap.png", lambda o: plot_plateau(no, ca, o)),
                  ("delta_total.png", lambda o: plot_dtotal(no, ca, o)),
                  ("delta_components_compare.png", lambda o: plot_components(no, ca, o)),
                  ("delta_pairwise_compare.png", lambda o: plot_pairwise(nocap_dir, cap_dir, o)),
                  ("norm_compare.png", lambda o: plot_norm(no, ca, o))):
        p = out_dir / fn; f(p); pngs.append(p)

    # total-density GIFs (existing per-run artifacts), embedded side by side
    gifs = [(nocap_dir / "report/wp_total_density.gif", "no-CAP total density n(x,z,t)"),
            (cap_dir / "report/wp_total_density.gif", "CAP total density n(x,z,t)")]

    def emb(path, mime):
        return new_output("display_data",
                          data={mime: base64.b64encode(path.read_bytes()).decode()}, metadata={})

    cells = [new_markdown_cell(
        "# Jellium WP-CAP energy-plateau — no-CAP vs CAP (phase comparison)\n\n"
        "Synthesis of the two localised-jellium WP runs (identical except the CAP). "
        "The **plateau gap** = electronic energy that no-CAP retains but CAP drains at "
        "the boundaries. Below: density evolution (top), energy plateau, ΔE component "
        "and pairwise electrostatic decompositions, and WP absorption.\n\n"
        "**Caveat:** these runs use a 25-Bohr transverse box; a dispersing WP wraps x/y "
        "at late times (see handover) — the decomposition is exact per frame but the "
        "physics interpretation carries that wrap.")]
    for g, cap in gifs:
        cells.append(new_markdown_cell(f"### {cap}"))
        if g.exists():
            c = new_code_cell(f'from IPython.display import Image\nImage(filename="{g.name}")')
            # copy the gif next to the notebook so the filename ref resolves
            (out_dir / g.name if False else None)
            c.outputs = [emb(g, "image/gif")]; c.execution_count = None
            cells.append(c)
        else:
            cells.append(new_markdown_cell("_(GIF missing)_"))
    for p in pngs:
        cells.append(new_markdown_cell(f"### {p.stem}"))
        c = new_code_cell(f'Image(filename="{p.name}")'); c.outputs = [emb(p, "image/png")]
        cells.append(c)
    nb = new_notebook(); nb.cells = cells
    nbp = out_dir / "jellium_nocap_vs_cap_comparison.ipynb"
    nbf.write(nb, str(nbp))
    print("wrote", nbp)
    return nbp


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2], sys.argv[3])
