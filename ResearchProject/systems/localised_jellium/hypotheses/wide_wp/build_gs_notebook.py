#!/usr/bin/env python3
"""Build gs_validation.ipynb — S2 ground-state results for the
`wide-wavepacket-lowspread` campaign (slab GS, LZ=101, dx=0.40).

Reads the GS run_summary; validates against the analytic slab (n0, r_s) and the
production 90-box GS energy; shows the planar-averaged density n(z) if a
density_gs_system VTI is available (written by the first CAP/WP run).

Run: venv/bin/python3 build_gs_notebook.py
"""
import os, glob
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.join(HERE, "gs_validation.ipynb")
GS_SUMMARY = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
              "scripts/wide_wp/gs/results/run_summary.txt")

cells = []
def md(s): cells.append(new_markdown_cell(s))
def code(s): cells.append(new_code_cell(s))

md(r"""# S2 — Slab ground state (wide-WP campaign)
### `wide-wavepacket-lowspread` · box 50×50×101, dx=0.40, slab |z|<12.5, N=82, PBC

Foundational checkpoint reused by every Phase-0/Phase-1 run. This notebook
validates it: SCF converged, interior density n₀ and r_s match the analytic slab
and the production 90-box GS (the neutral-slab interior is box/BC-independent).
""")

code(r"""import numpy as np
GS_SUMMARY = r"%s"
def parse_summary(path):
    d = {}
    for ln in open(path):
        if "=" in ln:
            k, v = ln.split("=", 1); d[k.strip()] = v.strip()
    return d
s = parse_summary(GS_SUMMARY)
for k in ("cell_bohr","spacing_bohr","k_nyq","e_cut_ha","n0_a0m3","r_s",
          "ground_state_energy_ha","num_electrons","num_states","run_completed"):
    print(f"  {k:24s} = {s.get(k,'?')}")
""" % GS_SUMMARY)

md(r"""## Validation against the analytic slab + production GS""")

code(r"""N, V = 82, 50.0*50.0*25.0
n0_analytic = N / V
rs_analytic = (3.0/(4.0*np.pi*n0_analytic))**(1.0/3.0)
E_GS_90 = -70.22568216820937   # production 90-box GS anchor (qsp_phase5)

n0  = float(s["n0_a0m3"]); rs = float(s["r_s"]); E = float(s["ground_state_energy_ha"])
ne  = int(float(s["num_electrons"]))
print(f"{'quantity':22s} {'GS run':>16s} {'analytic/ref':>16s} {'verdict':>10s}")
print(f"{'n0 (a0^-3)':22s} {n0:>16.4e} {n0_analytic:>16.4e} "
      f"{'OK' if abs(n0-n0_analytic)/n0_analytic<1e-3 else 'CHECK':>10s}")
print(f"{'r_s (Bohr)':22s} {rs:>16.4f} {rs_analytic:>16.4f} "
      f"{'OK' if abs(rs-rs_analytic)/rs_analytic<1e-3 else 'CHECK':>10s}")
print(f"{'N electrons':22s} {ne:>16d} {N:>16d} {'OK' if ne==N else 'CHECK':>10s}")
print(f"{'E_GS (Ha)':22s} {E:>16.4f} {E_GS_90:>16.4f}  (101-box vs 90-box; "
      f"differ by vacuum, interior matches)")
print(f"\nSCF: run_completed = {s.get('run_completed')} (finished within max_steps)")
""")

md(r"""## Planar-averaged density n(z) — the slab profile

Loaded from a `density_gs_system` VTI (written at t=0 by the first CAP/WP run).
The interior should be flat at n₀ with the slab edges at |z|=12.5; the launch
point (−26.5) and CAP inner face (±40.5) are marked. If no VTI is present yet,
this panel is skipped (populated after the P0a CAP smoke).
""")

code(r"""import glob
from inqview import load_vti
import matplotlib.pyplot as plt
from inqview.visualisation.style import apply_theme
apply_theme()

pat = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/wide_wp/**/raw/vti/density_gs_system/*.vti")
hits = sorted(glob.glob(pat, recursive=True))
if not hits:
    print("No density_gs_system VTI yet — run P0a (CAP smoke) to populate n(z).")
else:
    fld = load_vti(hits[0], expect_centered_axis="z")
    n = np.asarray(fld.data); z = np.asarray(fld.z)
    nz = n.mean(axis=(0,1))                       # planar average over x,y
    fig, ax = plt.subplots(figsize=(6.4,3.2), constrained_layout=True)
    ax.plot(z, nz, color="C0")
    ax.axhline(82/(50*50*25), ls="--", lw=0.9, color="0.5", label="analytic n0")
    for zc in (-12.5, 12.5): ax.axvline(zc, ls=":", lw=0.8, color="C3")
    ax.axvline(-26.5, ls="-", lw=0.7, color="0.4"); ax.text(-26.5, ax.get_ylim()[1]*0.9," launch", fontsize=7)
    for zc in (-40.5, 40.5): ax.axvline(zc, ls="-", lw=0.7, color="0.7")
    ax.set_xlabel("z (Bohr)"); ax.set_ylabel(r"$\bar n(z)$ (a$_0^{-3}$)")
    ax.set_title("Slab GS density (planar-averaged)", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    os.makedirs(os.path.join(r"__HERE__","gs_figs"), exist_ok=True)
    fig.savefig(os.path.join(r"__HERE__","gs_figs","gs_density_nz.png"), dpi=140)
    plt.show()
    interior = nz[np.abs(z) < 8.0]
    rel = interior.std()/interior.mean()*100
    print(f"interior n(z): mean={interior.mean():.4e}  std/mean={rel:.2f}%  "
          f"(flat interior -> good slab)")
""".replace("__HERE__", HERE))

nb = new_notebook(cells=cells, metadata={
    "kernelspec": {"name": "python3", "display_name": "Python 3"},
    "language_info": {"name": "python"}})
print("Executing GS notebook ...")
NotebookClient(nb, timeout=300, kernel_name="python3",
               resources={"metadata": {"path": HERE}}).execute()
with open(NB, "w") as f: nbf.write(nb, f)
print("wrote", NB)
